// BEST H-matrix solver: the C interface (best_hmatrix.h) that BEST calls.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "best_hmatrix.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "cluster.h"
#include "gmres.h"
#include "hmatrix.h"
#include "kernel.h"
#include "mat.h"
#include "parallel.h"

#ifdef __GLIBC__
#include <malloc.h>
#endif

namespace bh {

namespace {

// Settings that are not BEST options.  Each can be changed with the environment variable named.
const int kLeaf = 32;         // triangles per leaf cluster
const double kEta = 2.0;      // admissibility
const int kRestart = 30;      // GMRES restart length
const int kMaxSteps = 500;    // GMRES steps (batched products) per solve
const double kTolFactor = 0.1;  // GMRES stops at a relative residual of kTolFactor * eps
const int kForkLeaves = 2;      // the H-LU forks for blocks of at least this many leaf clusters a side

std::deque<std::string> g_log;

void say( const char* fmt, ... ) __attribute__( ( format( printf, 1, 2 ) ) );
void say( const char* fmt, ... ) {
   char buf[ 512 ];
   va_list ap;
   va_start( ap, fmt );
   std::vsnprintf( buf, sizeof buf, fmt, ap );
   va_end( ap );
   g_log.emplace_back( buf );
}

double now() {
   using namespace std::chrono;
   return duration<double>( steady_clock::now().time_since_epoch() ).count();
}

double gb( double bytes ) { return bytes / 1e9; }

double env_double( const char* name, double def ) {
   const char* s = std::getenv( name );
   if ( !s || !*s ) {
      return def;
   }
   char* end = nullptr;
   const double v = std::strtod( s, &end );
   return end != s ? v : def;
}

int env_threads() {
   const char* s = std::getenv( "OMP_NUM_THREADS" );
   if ( !s ) {
      return 0;
   }
   const int n = std::atoi( s );
   return n > 0 ? n : 0;
}

struct Solver {
   Mesh mesh;
   std::unique_ptr<Kernel> K;
   ClusterTree ct;
   std::unique_ptr<HMatrix> H;  // GR + the deflation term; with direct, its H-LU factors
   std::unique_ptr<HMatrix> P;  // GMRES: the preconditioner's H-LU factors
   bool direct = false;
   double eps = 0, tol = 0;
   int restart = kRestart, maxsteps = kMaxSteps;
};

std::unique_ptr<Solver> g;

// Memory per triangle in KB, fitted to the prototype's runs on sesmesh meshes of lysozyme
// (30,000 to 1,000,000 triangles): the H-matrix at eps 1e-4 takes -33 + 15.8 ln N, at 1e-6
// -124 + 32.5 ln N (interpolated in log eps between); the H-LU factors about 1.15 times that;
// the GMRES preconditioner at 1e-1 with the near field compressed -13 + 3.5 ln N.
double h_kb_per_tri( int n, double eps ) {
   const double l = -std::log10( eps ), ln = std::log( double( n ) );
   const double a = -33.0 - 45.5 * ( l - 4.0 ), b = 15.8 + 8.35 * ( l - 4.0 );
   return std::max( 45.0, a + b * ln );  // at least the dense near field
}

double estimate_bytes( int n, double eps, bool direct, double eps_pre ) {
   double kb = h_kb_per_tri( n, eps ) * ( direct ? 1.15 : 1.0 );
   if ( !direct ) {
      kb += eps_pre >= 0.1 ? std::max( 3.0, -13.0 + 3.5 * std::log( double( n ) ) ) : h_kb_per_tri( n, eps_pre );
   }
   // Heap overhead, and GMRES's vectors for six right-hand sides.
   kb = 1.3 * kb + ( direct ? 0.5 : 4.5 );
   return kb * 1e3 * n;
}

template <class F>
int guarded( F&& f ) {
   try {
      return f();
   } catch ( const std::bad_alloc& ) {
      say( "error: out of memory" );
   } catch ( const std::exception& e ) {
      say( "error: %s", e.what() );
   }
   return 1;
}

void to_tree( const ClusterTree& ct, const double* B, int q, double* Bt ) {
   const int n = int( ct.perm.size() );
   const size_t n3 = 3 * size_t( n );
   parallel_for( n, [&]( long p, int ) {
      const size_t s = 3 * size_t( ct.perm[ p ] ), d = 3 * size_t( p );
      for ( int c = 0; c < q; ++c ) {
         for ( int a = 0; a < 3; ++a ) {
            Bt[ c * n3 + d + a ] = B[ c * n3 + s + a ];
         }
      }
   }, 1024 );
}

void from_tree( const ClusterTree& ct, const double* Bt, int q, double* B ) {
   const int n = int( ct.perm.size() );
   const size_t n3 = 3 * size_t( n );
   parallel_for( n, [&]( long p, int ) {
      const size_t s = 3 * size_t( p ), d = 3 * size_t( ct.perm[ p ] );
      for ( int c = 0; c < q; ++c ) {
         for ( int a = 0; a < 3; ++a ) {
            B[ c * n3 + d + a ] = Bt[ c * n3 + s + a ];
         }
      }
   }, 1024 );
}

void say_stats( const char* what, const HMatrix& h, int n, double seconds ) {
   const HStats st = h.stats();
   const double dense = 8.0 * 9.0 * double( n ) * n;
   say( "%s in %.1f s: %.2f GB (%.0f KB per triangle, %.2f%% of dense), %ld dense and %ld low-rank blocks, "
        "rank max %d mean %.1f",
        what, seconds, gb( st.bytes() ), st.bytes() / n / 1e3, 100.0 * st.bytes() / dense, st.n_dense, st.n_lowrank,
        st.kmax, st.kmean );
}

// glibc keeps the memory freed by the truncations in its heaps; hand what it can back to the system.
void trim() {
#ifdef __GLIBC__
   malloc_trim( 0 );
#endif
}

void say_memory() {
   say( "memory: resident %.2f GB, peak %.2f GB, available %.2f GB", gb( double( current_rss() ) ),
        gb( double( peak_rss() ) ), gb( double( mem_available() ) ) );
}

}  // namespace

}  // namespace bh

using namespace bh;

extern "C" {

int bhm_available( void ) { return 1; }

int bhm_create( int n, const double* v1, const double* v2, const double* v3, const double* area, double afac,
                double cg, int nthreads ) {
   return guarded( [&] {
      if ( n < 2 ) {
         say( "error: the H-matrix solver needs at least 2 triangles" );
         return 1;
      }
      g = std::make_unique<Solver>();
      set_threads( nthreads > 0 ? nthreads : env_threads() );
      g->mesh = mesh_from_best( n, v1, v2, v3, area, afac );
      g->K = std::make_unique<Kernel>( g->mesh, cg );
      say( "H-matrix solver: %d triangles, %d threads", n, threads() );
      return 0;
   } );
}

void bhm_incenters( double* x ) {
   const int n = g->mesh.n;
   const std::vector<double>& xc = g->K->x();
   for ( int i = 0; i < n; ++i ) {
      for ( int c = 0; c < 3; ++c ) {
         x[ size_t( c ) * n + i ] = xc[ 3 * size_t( i ) + c ];
      }
   }
}

void bhm_entries( int i0, int m, int j0, int nc, double* a, int lda ) {
   if ( m <= 0 || nc <= 0 ) {
      return;
   }
   // One kernel block per pair of triangles covering the range; columns of triangles in parallel.
   const int t0 = i0 / 3, t1 = ( i0 + m - 1 ) / 3, s0 = j0 / 3, s1 = ( j0 + nc - 1 ) / 3;
   parallel_for( s1 - s0 + 1, [&]( long ds, int ) {
      const int n2 = s0 + int( ds );
      double blk[ 9 ];
      for ( int n1 = t0; n1 <= t1; ++n1 ) {
         g->K->block( n1, n2, blk );
         for ( int bb = 0; bb < 3; ++bb ) {
            const int j = 3 * n2 + bb - j0;
            if ( j < 0 || j >= nc ) {
               continue;
            }
            for ( int aa = 0; aa < 3; ++aa ) {
               const int i = 3 * n1 + aa - i0;
               if ( i >= 0 && i < m ) {
                  a[ size_t( j ) * lda + i ] = blk[ 3 * aa + bb ];
               }
            }
         }
      }
   } );
}

double bhm_diag_mean( void ) { return g->K->diag_mean(); }

int bhm_build( double eps, int direct, double eps_pre, double gdef, const double* xmu ) {
   return guarded( [&] {
      Solver& s = *g;
      const int n = s.mesh.n, n3 = 3 * n;
      s.direct = direct != 0;
      s.eps = eps;
      s.tol = env_double( "BEST_HMATRIX_GMRES_TOL", kTolFactor * eps );
      // The H-LU forks for blocks of at least this many leaf clusters a side.
      HMatrix::set_fork( int( env_double( "BEST_HMATRIX_FORK", kForkLeaves ) ) );
      double t = now();
      s.ct = build_cluster_tree( s.mesh, s.K->x(), kLeaf );
      say( "cluster tree: depth %d, %zu leaf clusters of %d to %d triangles, %.2f s", s.ct.depth, s.ct.leaves.size(),
           n >> s.ct.depth, ( n + ( 1 << s.ct.depth ) - 1 ) >> s.ct.depth, now() - t );

      const double need = estimate_bytes( n, eps, s.direct, eps_pre ), avail = double( mem_available() );
      say( "memory: about %.1f GB needed (estimate), %.1f GB available", gb( need ), gb( avail ) );
      if ( avail > 0 && need > avail && env_double( "BEST_HMATRIX_MEMCHECK", 1 ) != 0 ) {
         say( "error: not enough memory for the H-matrix solver at eps %.0e; a larger -he eps needs less "
              "(set BEST_HMATRIX_MEMCHECK=0 to try anyway)",
              eps );
         return 1;
      }

      HOptions opt;
      opt.eps = eps;
      opt.eta = kEta;
      long calls = 0;
      t = now();
      s.H = std::make_unique<HMatrix>( HMatrix::build( s.ct, *s.K, opt, &calls ) );
      std::vector<double> mt( n3 );
      to_tree( s.ct, xmu, 1, mt.data() );
      s.H->add_rank1( gdef, mt.data(), mt.data() );
      say_stats( "H-matrix", *s.H, n, now() - t );
      say( "H-matrix: eps %.0e, eta %.1f, %.3g kernel blocks per N^2", eps, kEta, double( calls ) / ( double( n ) * n ) );
      trim();
      say_memory();

      long lrc = 0;
      double lrs = 0;
      lowrank_stats( lrc, lrs, true );
      HMatrix::profile( true );
      const double cpu0 = cpu_seconds();
      t = now();
      if ( s.direct ) {
         const double piv = s.H->lu( eps );
         say_stats( "H-LU factors", *s.H, n, now() - t );
         say( "H-LU: eps %.0e, smallest pivot ratio %.2e", eps, piv );
      } else {
         s.P = std::make_unique<HMatrix>( s.H->copy( eps_pre, true ) );
         const double piv = s.P->lu( eps_pre );
         say_stats( "GMRES preconditioner (H-LU factors)", *s.P, n, now() - t );
         say( "GMRES preconditioner: eps %.0e, near field compressed, smallest pivot ratio %.2e; GMRES(%d) to a "
              "relative residual of %.0e",
              eps_pre, piv, s.restart, s.tol );
      }
      lowrank_stats( lrc, lrs, true );
      say( "factorization: %.0f CPU s, %ld truncations taking %.0f CPU s; by site: %s", cpu_seconds() - cpu0, lrc, lrs,
           HMatrix::profile( true ).c_str() );
      trim();
      say_memory();
      return 0;
   } );
}

int bhm_solve( int nrhs, double* b ) {
   return guarded( [&] {
      Solver& s = *g;
      const int n3 = 3 * s.mesh.n;
      const double t = now();
      std::vector<double> Bt( size_t( n3 ) * nrhs );
      to_tree( s.ct, b, nrhs, Bt.data() );
      if ( s.direct ) {
         s.H->lu_solve( Bt.data(), n3, nrhs );
         from_tree( s.ct, Bt.data(), nrhs, b );
         say( "solve: %d right-hand sides with the H-LU factors, %.2f s", nrhs, now() - t );
         return 0;
      }
      std::vector<double> Xt( size_t( n3 ) * nrhs );
      const HMatrix& H = *s.H;
      const HMatrix& P = *s.P;
      double ta = 0, tm = 0;
      int na = 0, nm = 0;
      Apply A = [&]( const double* X, int ldx, int q, double* Y, int ldy ) {
         const double t0 = now();
         H.matvec( X, ldx, q, Y, ldy );
         ta += now() - t0;
         na += q;
      };
      Precond M = [&]( double* X, int ldx, int q ) {
         const double t0 = now();
         P.lu_solve( X, ldx, q );
         tm += now() - t0;
         nm += q;
      };
      const GmresResult r = gmres( n3, nrhs, A, M, Bt.data(), n3, Xt.data(), n3, s.tol, s.restart, s.maxsteps );
      from_tree( s.ct, Xt.data(), nrhs, b );
      int imin = 1 << 30, imax = 0;
      double worst = 0;
      for ( int c = 0; c < nrhs; ++c ) {
         imin = std::min( imin, r.iters[ c ] );
         imax = std::max( imax, r.iters[ c ] );
         worst = std::max( worst, r.relres[ c ] );
      }
      say( "solve: %d right-hand sides by GMRES, %d to %d iterations, largest relative residual %.1e, %.2f s "
           "(products %.2f s for %d vectors, preconditioner %.2f s for %d)",
           nrhs, imin, imax, worst, now() - t, ta, na, tm, nm );
      if ( !( worst <= s.tol ) ) {
         say( "error: GMRES did not reach a relative residual of %.0e in %d steps", s.tol, s.maxsteps );
         return 1;
      }
      return 0;
   } );
}

int bhm_log_line( char* buf, int len ) {
   if ( g_log.empty() ) {
      return -1;
   }
   const std::string& s = g_log.front();
   const int k = std::min( len, int( s.size() ) );
   std::memcpy( buf, s.data(), size_t( k ) );
   g_log.pop_front();
   return k;
}

void bhm_destroy( void ) { g.reset(); }

}  // extern "C"
