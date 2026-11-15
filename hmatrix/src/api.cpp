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
#include "meshcheck.h"
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
// The GMRES preconditioner is rebuilt 10 times finer, down to kPreMin, when its H-LU breaks down: a
// pivot that is 0 or not finite, or a smallest pivot ratio below kPivotCollapse (smooth meshes give
// 1e-5 to 0.2; the probe-0 and tunnel meshes that broke it, 1e-29 to 1e-41), or when GMRES with it
// breaks down or does not converge.
const double kPivotCollapse = 1e-10;
const double kPreMin = 1e-3;
// Collocation points closer than this (A), on triangles with opposite normals, are reported with
// such a breakdown (meshcheck.h).
const double kClose = 0.05;

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

template <class Elt>
struct Mats {
   std::unique_ptr<HMatrixT<Elt>> H;  // GR + the deflation term; with direct, its H-LU factors
   std::unique_ptr<HMatrixT<Elt>> P;  // GMRES: the preconditioner's H-LU factors
};

struct Solver {
   Mesh mesh;
   std::unique_ptr<Kernel> K;
   ClusterTree ct;
   bool single = false;  // the blocks in single precision: f, otherwise d
   Mats<double> d;
   Mats<float> f;
   bool direct = false;
   double eps = 0, tol = 0;
   double eps_pre = 0;     // GMRES: the preconditioner's eps now
   bool reported = false;  // the mesh's defects have been logged
   int restart = kRestart, maxsteps = kMaxSteps;
};

template <class Elt>
Mats<Elt>& mats( Solver& s );
template <>
Mats<double>& mats<double>( Solver& s ) {
   return s.d;
}
template <>
Mats<float>& mats<float>( Solver& s ) {
   return s.f;
}

std::unique_ptr<Solver> g;

// The normal-field deflation: gamma[k] for surface k; comp, the surface of each triangle in tree
// order, is empty for one surface.
struct Deflation {
   std::vector<double> gamma;
   std::vector<int> comp;
};

// Memory per triangle in KB, fitted to the prototype's runs on sesmesh meshes of lysozyme
// (30,000 to 1,000,000 triangles): the H-matrix at eps 1e-4 takes -33 + 15.8 ln N, at 1e-6
// -124 + 32.5 ln N (interpolated in log eps between); the H-LU factors about 1.15 times that;
// the GMRES preconditioner at 1e-1 with the near field compressed -13 + 3.5 ln N.
double h_kb_per_tri( int n, double eps ) {
   const double l = -std::log10( eps ), ln = std::log( double( n ) );
   const double a = -33.0 - 45.5 * ( l - 4.0 ), b = 15.8 + 8.35 * ( l - 4.0 );
   return std::max( 45.0, a + b * ln );  // at least the dense near field
}

double estimate_bytes( int n, double eps, bool direct, double eps_pre, bool single ) {
   double kb = h_kb_per_tri( n, eps ) * ( direct ? 1.15 : 1.0 );
   if ( !direct ) {
      kb += eps_pre >= 0.1 ? std::max( 3.0, -13.0 + 3.5 * std::log( double( n ) ) ) : h_kb_per_tri( n, eps_pre );
   }
   // Single precision halves the blocks.
   if ( single ) {
      kb *= 0.5;
   }
   // Heap overhead and GMRES's vectors; the direct H-LU's peak is about 1.7 times its factors (the
   // products of large hierarchical blocks build full-height temporaries).
   kb = ( direct ? 1.75 : 1.3 ) * kb + ( direct ? 0.5 : 4.5 );
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

template <class Elt>
void say_stats( const char* what, const HMatrixT<Elt>& h, int n, double seconds ) {
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

// The triangle (BEST's numbering, from 1) of an unknown in tree order.
int tri( const Solver& s, int row ) { return row >= 0 ? s.ct.perm[ row / 3 ] + 1 : 0; }

// Logs, once, the features of the mesh that make GR ill-conditioned (meshcheck.h).
void report_mesh( Solver& s ) {
   if ( s.reported ) {
      return;
   }
   s.reported = true;
   const MeshDefects d = mesh_defects( s.mesh, s.K->x(), kClose );
   say( "mesh: smallest angle %.2g deg (triangle %d); %ld triangles with an angle below 1 deg", d.min_angle,
        d.min_angle_tri + 1, d.n_sliver );
   if ( d.n_close > 0 ) {
      say( "mesh: %ld pairs of triangles with opposite normals have collocation points closer than %.2f A (closest "
           "%.2g A: triangles %d and %d)",
           d.n_close, kClose, d.close_min, d.close_a + 1, d.close_b + 1 );
   } else {
      say( "mesh: no two triangles with opposite normals have collocation points closer than %.2f A", kClose );
   }
   say( "mesh: slivers folded onto their neighbors and sheets facing each other across narrow gaps make GR "
        "ill-conditioned; probe-0 surfaces and narrow tunnels have them, smoother meshes (a larger probe radius, "
        "tunnels closed) do not" );
}

// The factorization's profile line.
void say_profile( double cpu0 ) {
   long lrc = 0;
   double lrs = 0;
   lowrank_stats( lrc, lrs, true );
   say( "factorization: %.0f CPU s, %ld truncations taking %.0f CPU s; by site: %s", cpu_seconds() - cpu0, lrc, lrs,
        HMatrix::profile( true ).c_str() );
}

void reset_profile() {
   long lrc = 0;
   double lrs = 0;
   lowrank_stats( lrc, lrs, true );
   HMatrix::profile( true );
}

// -hlu: the H-LU factorization of the H-matrix at eps, in place.
template <class Elt>
int factor_direct( Solver& s ) {
   const int n = s.mesh.n;
   Mats<Elt>& m = mats<Elt>( s );
   reset_profile();
   const double cpu0 = cpu_seconds(), t = now();
   double piv = 0;
   try {
      piv = m.H->lu( s.eps );
   } catch ( const ZeroPivot& z ) {
      report_mesh( s );
      say( "error: H-LU: %s pivot at triangle %d; without -hlu BEST solves by GMRES, refining its preconditioner as "
           "needed",
           z.finite ? "a zero" : "a NaN or Inf", tri( s, z.row ) );
      return 1;
   }
   say_stats( "H-LU factors", *m.H, n, now() - t );
   say( "H-LU: eps %.0e, smallest pivot ratio %.2e (triangle %d)", s.eps, piv, tri( s, HMatrixT<Elt>::pivot_row() ) );
   say_profile( cpu0 );
   trim();
   say_memory();
   if ( !( piv >= kPivotCollapse ) ) {
      report_mesh( s );
      say( "error: H-LU: the factorization broke down (smallest pivot ratio %.2e); without -hlu BEST solves by GMRES, "
           "refining its preconditioner as needed",
           piv );
      return 1;
   }
   return 0;
}

// GMRES: builds the preconditioner, the H-LU of a copy of the H-matrix at s.eps_pre with the near
// field compressed too.  Returns 0; 1 if its smallest pivot ratio is below kPivotCollapse; 2 if a
// pivot was 0 or not finite (the factors are then unusable).
template <class Elt>
int build_pre( Solver& s ) {
   const int n = s.mesh.n;
   Mats<Elt>& m = mats<Elt>( s );
   m.P.reset();
   trim();
   reset_profile();
   const double cpu0 = cpu_seconds(), t = now();
   m.P = std::make_unique<HMatrixT<Elt>>( m.H->copy( s.eps_pre, true ) );
   int st = 0;
   try {
      const double piv = m.P->lu( s.eps_pre );
      say_stats( "GMRES preconditioner (H-LU factors)", *m.P, n, now() - t );
      say( "GMRES preconditioner: eps %.0e, near field compressed, smallest pivot ratio %.2e (triangle %d); GMRES(%d) "
           "to a relative residual of %.0e",
           s.eps_pre, piv, tri( s, HMatrixT<Elt>::pivot_row() ), s.restart, s.tol );
      if ( !( piv >= kPivotCollapse ) ) {
         say( "GMRES preconditioner: broke down (smallest pivot ratio below %.0e)", kPivotCollapse );
         st = 1;
      }
   } catch ( const ZeroPivot& z ) {
      say( "GMRES preconditioner: eps %.0e, near field compressed: broke down, %s pivot at triangle %d", s.eps_pre,
           z.finite ? "a zero" : "a NaN or Inf", tri( s, z.row ) );
      st = 2;
   }
   say_profile( cpu0 );
   trim();
   say_memory();
   return st;
}

// GMRES: moves the preconditioner to an eps 10 times smaller, unless it is at kPreMin already or the
// memory for it is not available (then it says so and returns false, and the old one stays).  Logs
// the mesh's defects the first time.
template <class Elt>
bool finer( Solver& s ) {
   report_mesh( s );
   if ( s.eps_pre <= kPreMin * 1.001 ) {
      return false;
   }
   const double next = s.eps_pre / 10;
   // At 1e-3 the preconditioner took 0.73 times the H-matrix's memory on a probe-0 lysozyme mesh; at
   // 1e-2 less.  The H-matrix is there already, and the old preconditioner goes first.
   const HMatrixT<Elt>& H = *mats<Elt>( s ).H;
   const double need = ( next < 5e-3 ? 1.0 : 0.5 ) * 1.3 * H.stats().bytes();
   const HMatrixT<Elt>* P = mats<Elt>( s ).P.get();
   const double avail = double( mem_available() ) + ( P ? P->stats().bytes() : 0.0 );
   if ( avail > 0 && need > avail && env_double( "BEST_HMATRIX_MEMCHECK", 1 ) != 0 ) {
      say( "GMRES preconditioner: one at eps %.0e would need about %.1f GB, %.1f GB available (set "
           "BEST_HMATRIX_MEMCHECK=0 to try anyway)",
           next, gb( need ), gb( avail ) );
      return false;
   }
   say( "GMRES preconditioner: rebuilding it at eps %.0e", next );
   s.eps_pre = next;
   return true;
}

// bhm_build after the cluster tree and the memory check, with the blocks in Elt.
template <class Elt>
int build_t( Solver& s, double eps, double eps_pre, const Deflation& def, const double* xmu ) {
   const int n = s.mesh.n, n3 = 3 * n;
   Mats<Elt>& m = mats<Elt>( s );
   HOptions opt;
   opt.eps = eps;
   opt.eta = kEta;
   long calls = 0;
   double t = now();
   m.H = std::make_unique<HMatrixT<Elt>>( HMatrixT<Elt>::build( s.ct, *s.K, opt, &calls ) );
   std::vector<double> mt( n3 );
   to_tree( s.ct, xmu, 1, mt.data() );
   if ( def.comp.empty() ) {
      m.H->add_rank1( def.gamma[ 0 ], mt.data(), mt.data() );
   } else {
      m.H->add_deflation( def.gamma, def.comp, mt.data() );
   }
   say_stats( "H-matrix", *m.H, n, now() - t );
   say( "H-matrix: eps %.0e, eta %.1f, %.3g kernel blocks per N^2, %s precision", eps, kEta,
        double( calls ) / ( double( n ) * n ), s.single ? "single" : "double" );
   trim();
   say_memory();

   if ( s.direct ) {
      return factor_direct<Elt>( s );
   }
   // GMRES: the preconditioner, finer while it breaks down.
   s.eps_pre = eps_pre;
   for ( ;; ) {
      const int st = build_pre<Elt>( s );
      if ( st == 0 ) {
         return 0;
      }
      if ( !finer<Elt>( s ) ) {
         // At the finest level a small pivot ratio may still do; a zero pivot may not.
         if ( st == 2 ) {
            say( "error: the GMRES preconditioner broke down at eps %.0e", s.eps_pre );
            return 1;
         }
         return 0;
      }
   }
}

// True if every solution is finite; otherwise says so.
bool finite( const std::vector<double>& X ) {
   for ( double x : X ) {
      if ( !std::isfinite( x ) ) {
         say( "error: the solution is not finite (NaN or Inf)" );
         return false;
      }
   }
   return true;
}

// bhm_solve with the blocks in Elt: Bt holds the right-hand sides in tree order, and the solutions
// on return.
template <class Elt>
int solve_t( Solver& s, int nrhs, std::vector<double>& Bt ) {
   const int n3 = 3 * s.mesh.n;
   const double t = now();
   Mats<Elt>& m = mats<Elt>( s );
   if ( s.direct ) {
      m.H->lu_solve( Bt.data(), n3, nrhs );
      say( "solve: %d right-hand sides with the H-LU factors, %.2f s", nrhs, now() - t );
      return finite( Bt ) ? 0 : 1;
   }
   std::vector<double> Xt( size_t( n3 ) * nrhs );
   const HMatrixT<Elt>& H = *m.H;
   for ( ;; ) {
      const HMatrixT<Elt>& P = *m.P;
      const double t1 = now();
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
      int imin = 1 << 30, imax = 0;
      double worst = 0;
      for ( int c = 0; c < nrhs; ++c ) {
         imin = std::min( imin, r.iters[ c ] );
         imax = std::max( imax, r.iters[ c ] );
         // Not std::max, which would pass over a NaN; a NaN stays.
         if ( std::isnan( r.relres[ c ] ) || r.relres[ c ] > worst ) {
            worst = r.relres[ c ];
         }
      }
      say( "solve: %d right-hand sides by GMRES, %d to %d iterations, largest relative residual %.1e, %.2f s "
           "(products %.2f s for %d vectors, preconditioner %.2f s for %d)",
           nrhs, imin, imax, worst, now() - t1, ta, na, tm, nm );
      if ( worst <= s.tol ) {
         Bt.swap( Xt );
         return finite( Bt ) ? 0 : 1;
      }
      if ( r.broke ) {
         say( "GMRES broke down: the residual grew instead of shrinking (the preconditioner is ill-conditioned)" );
      }
      // Again with a finer preconditioner, as long as there is one.
      bool again = false;
      while ( finer<Elt>( s ) ) {
         const int st = build_pre<Elt>( s );
         if ( st == 0 || ( st == 1 && s.eps_pre <= kPreMin * 1.001 ) ) {
            again = true;
            break;
         }
      }
      if ( !again ) {
         say( "error: GMRES did not reach a relative residual of %.0e in %d steps%s", s.tol, s.maxsteps,
              r.broke ? " (it broke down)" : "" );
         return 1;
      }
   }
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

int bhm_build( double eps, int direct, double eps_pre, int ndef, const double* gdef, const int* comp, const double* xmu,
               int single ) {
   return guarded( [&] {
      Solver& s = *g;
      const int n = s.mesh.n;
      s.direct = direct != 0;
      s.single = single != 0;
      s.eps = eps;
      s.tol = env_double( "BEST_HMATRIX_GMRES_TOL", kTolFactor * eps );
      eps_pre = env_double( "BEST_HMATRIX_PRE_EPS", eps_pre );
      // The H-LU forks for blocks of at least this many leaf clusters a side.
      HMatrix::set_fork( int( env_double( "BEST_HMATRIX_FORK", kForkLeaves ) ) );
      double t = now();
      s.ct = build_cluster_tree( s.mesh, s.K->x(), kLeaf );
      say( "cluster tree: depth %d, %zu leaf clusters of %d to %d triangles, %.2f s", s.ct.depth, s.ct.leaves.size(),
           n >> s.ct.depth, ( n + ( 1 << s.ct.depth ) - 1 ) >> s.ct.depth, now() - t );

      const double need = estimate_bytes( n, eps, s.direct, eps_pre, s.single ), avail = double( mem_available() );
      say( "memory: about %.1f GB needed (estimate), %.1f GB available", gb( need ), gb( avail ) );
      if ( avail > 0 && need > avail && env_double( "BEST_HMATRIX_MEMCHECK", 1 ) != 0 ) {
         say( "error: not enough memory for the H-matrix solver at eps %.0e; %s (set BEST_HMATRIX_MEMCHECK=0 to try "
              "anyway)",
              eps, s.single ? "a larger -he eps needs less, but 1e-3 saves only about 20% over 1e-4, at about 50 "
                              "times the error"
                            : "-hs (single precision) needs about half" );
         return 1;
      }
      Deflation def;
      def.gamma.assign( gdef, gdef + std::max( 1, ndef ) );
      if ( ndef > 1 ) {
         def.comp.resize( n );
         for ( int p = 0; p < n; ++p ) {
            def.comp[ p ] = comp[ s.ct.perm[ p ] ] - 1;  // tree order, from 0
         }
      }
      return s.single ? build_t<float>( s, eps, eps_pre, def, xmu ) : build_t<double>( s, eps, eps_pre, def, xmu );
   } );
}

int bhm_solve( int nrhs, double* b ) {
   return guarded( [&] {
      Solver& s = *g;
      const int n3 = 3 * s.mesh.n;
      std::vector<double> Bt( size_t( n3 ) * nrhs );
      to_tree( s.ct, b, nrhs, Bt.data() );
      const int info = s.single ? solve_t<float>( s, nrhs, Bt ) : solve_t<double>( s, nrhs, Bt );
      from_tree( s.ct, Bt.data(), nrhs, b );
      return info;
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
