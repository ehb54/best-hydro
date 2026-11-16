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
#include "memory.h"
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
const double kSampleFrac = 0.02;  // the memory estimate's sample: this fraction of the leaves of each kind
const long kSampleMin = 8;        // ... but at least this many (or all)
const double kFloorFrac = 0.05;   // the guard stops the solver below this fraction of the node's memory
const double kFloorMin = 1e9;     // ... or below this many bytes, whichever is more
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
   double floor = 0;       // the memory guard's floor (bytes)
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

// The memory estimate.  The block tree is built before any entry is computed, so the counts and
// sizes of its dense and low-rank leaves are known exactly; only the ranks of the low-rank leaves
// are not.  They are measured on a sample of the leaves, filled first and kept for the H-matrix
// (so the sample costs no extra work): per block level, the ranks weighted by nr + nc, which is what
// the memory of a level's leaves goes with.  The preconditioner's sizes come from the same sample,
// truncated again at its eps.  Rough surfaces (unions of atoms, tunnels) have several times as
// many blocks per triangle as smooth ones, and so need several times the memory per triangle; the
// ranks change much less.
//
// Sizes below are in matrix entries; a block's width is nr + nc.
struct MemPlan {
   long n_dense = 0, n_lowrank = 0, n_hblocks = 0;
   double dense_entries = 0;  // the dense leaves: sum of nr nc
   double diag_entries = 0;   // ... those on the diagonal
   double near_entries = 0;   // ... those off it, the near field
   double lr_width = 0;       // the low-rank leaves: sum of widths
   double lr_entries = 0;     // ... sum of k times the width (estimated)
   double lr_ranks = 0;       // ... sum of k (estimated)
   double pre_entries = 0;    // the preconditioner before its H-LU (estimated)
   long n_sample = 0;         // leaves in the sample
   double sample_entries = 0; // their entries (already in memory when the estimate is made)
   double sample_s = 0;       // seconds taken to fill it
};

template <class Elt>
long count_hblocks( const BlockT<Elt>* b ) {
   if ( b->type != BType::H ) {
      return 0;
   }
   long c = 1;
   for ( const auto& ch : b->ch ) {
      c += count_hblocks( ch.get() );
   }
   return c;
}

// Takes kSampleFrac of v (at least kSampleMin, at most all), evenly spaced, into out.
void pick_sample( const std::vector<long>& v, std::vector<long>& out ) {
   const long n = long( v.size() );
   const long k = std::min( n, std::max( kSampleMin, long( std::ceil( kSampleFrac * double( n ) ) ) ) );
   for ( long j = 0; j < k; ++j ) {
      out.push_back( v[ size_t( ( ( 2 * j + 1 ) * n ) / ( 2 * k ) ) ] );
   }
}

// Fills a sample of H's leaves (marked in filled) and measures the block tree.  Without direct,
// it also measures the preconditioner's ranks at eps_pre: the sampled low-rank leaves truncated
// again, the sampled near-field leaves compressed (as HMatrix::copy does).
template <class Elt>
MemPlan plan_memory( HMatrixT<Elt>& H, const Kernel& K, const HOptions& opt, bool direct, double eps_pre,
                     std::vector<char>& filled, long* calls ) {
   MemPlan p;
   const std::vector<BlockT<Elt>*>& L = H.leaves();
   const int depth = H.tree().depth;
   std::vector<std::vector<long>> lr( depth + 1 );
   std::vector<double> lr_w( depth + 1, 0.0 );
   std::vector<long> near;
   for ( size_t i = 0; i < L.size(); ++i ) {
      const BlockT<Elt>* b = L[ i ];
      const double e = double( b->nr ) * b->nc;
      if ( b->type == BType::R ) {
         ++p.n_lowrank;
         p.lr_width += b->nr + b->nc;
         lr_w[ b->level ] += b->nr + b->nc;
         lr[ b->level ].push_back( long( i ) );
      } else {
         ++p.n_dense;
         p.dense_entries += e;
         if ( b->t == b->s ) {
            p.diag_entries += e;
         } else {
            p.near_entries += e;
            near.push_back( long( i ) );
         }
      }
   }
   p.n_hblocks = count_hblocks( H.root() );
   std::vector<long> sample;
   std::vector<size_t> lr_beg( depth + 2, 0 );
   for ( int l = 0; l <= depth; ++l ) {
      lr_beg[ l ] = sample.size();
      pick_sample( lr[ l ], sample );
   }
   lr_beg[ depth + 1 ] = sample.size();
   const size_t near_beg = sample.size();
   if ( !direct ) {
      pick_sample( near, sample );
   }
   const double t = now();
   H.fill( K, opt, sample, calls );
   p.sample_s = now() - t;
   p.n_sample = long( sample.size() );
   for ( long i : sample ) {
      filled[ size_t( i ) ] = 1;
      const BlockT<Elt>* b = L[ size_t( i ) ];
      p.sample_entries += b->type == BType::R ? double( b->rank() ) * ( b->nr + b->nc ) : double( b->nr ) * b->nc;
   }
   if ( memory_low() ) {
      return p;
   }
   // The preconditioner's entries for each sampled leaf.
   std::vector<double> pre( sample.size(), 0.0 );
   if ( !direct ) {
      parallel_for( long( sample.size() ), [&]( long j, int ) {
         const BlockT<Elt>* b = L[ size_t( sample[ size_t( j ) ] ) ];
         Mat U, V;
         if ( b->type == BType::R ) {
            U = to_double( b->U );
            V = to_double( b->V );
            truncate( U, V, eps_pre );
            pre[ size_t( j ) ] = double( U.n ) * ( b->nr + b->nc );
         } else {
            Mat D = to_double( b->D );
            dense_to_lowrank( D, eps_pre, U, V );
            pre[ size_t( j ) ] = std::min( double( U.n ) * ( b->nr + b->nc ), double( b->nr ) * b->nc );
         }
      } );
   }
   for ( int l = 0; l <= depth; ++l ) {
      double w = 0, kw = 0, k = 0, pw = 0;
      for ( size_t j = lr_beg[ l ]; j < lr_beg[ l + 1 ]; ++j ) {
         const BlockT<Elt>* b = L[ size_t( sample[ j ] ) ];
         w += b->nr + b->nc;
         kw += double( b->rank() ) * ( b->nr + b->nc );
         k += b->rank();
         pw += pre[ j ];
      }
      if ( w > 0 ) {
         p.lr_entries += lr_w[ l ] * kw / w;
         p.pre_entries += lr_w[ l ] * pw / w;
         p.lr_ranks += double( lr[ l ].size() ) * k / double( lr_beg[ l + 1 ] - lr_beg[ l ] );
      }
   }
   if ( !direct ) {
      double e = 0, pe = 0;
      for ( size_t j = near_beg; j < sample.size(); ++j ) {
         const BlockT<Elt>* b = L[ size_t( sample[ j ] ) ];
         e += double( b->nr ) * b->nc;
         pe += pre[ j ];
      }
      p.pre_entries += p.diag_entries + ( e > 0 ? p.near_entries * pe / e : 0.0 );
   }
   return p;
}

// The peak memory of the process, in bytes: what it uses when the estimate is made (BEST's arrays,
// the mesh, the kernel's tables, the cluster tree, the sample) and what is still to come.
struct MemEstimate {
   double used = 0;     // resident now
   double h = 0;        // the H-matrix (with -hlu, its factors)
   double pre = 0;      // the GMRES preconditioner at the peak (in its H-LU, or after it with the vectors)
   double vectors = 0;  // GMRES: the Krylov bases and the other vectors, and the products' temporaries
   double other = 0;    // the blocks' bookkeeping, the threads' heaps, the factorization's temporaries
   double peak = 0;
};

// The model's constants, from runs in 11-2026 on meshes of 26,000 to 667,000 triangles: SOMO's
// SES meshes of spheres and proteins (probe 0 to 3 A), and sesmesh's rough probe-0 meshes of
// proteins (unions of atoms), in single and double precision, with 8 to 64 threads.
//
// With GMRES the peak comes in one of two places.  During the preconditioner's H-LU its blocks from
// before the H-LU and its factors are in memory together: up to 2.35 times its size before the
// H-LU, as its factors came out 1.1 to 1.25 times that size.  (1ADO's union of atoms, whose H-LU
// broke down to factors of 0.3 times, took 1.2 times.)  After the H-LU the heap keeps much of what
// it freed, and GMRES adds its vectors: 1.1 to 1.7 times the preconditioner's size stays resident.
// The estimate takes the larger of the two.  With 64 threads the peak was 2% to 4% higher than
// with 32: each thread's heap arena keeps some of what it frees.  With these constants the
// estimate came out 1% to 10% above the peak with GMRES, and 8% to 14% with -hlu.
//
// The estimate covers the first preconditioner, at eps_pre.  If its H-LU breaks down, finer()
// checks the memory again before each finer one.
const double kBlockBytes = 300;     // per block of each matrix: the block and its allocations
const double kPreLU = 2.4;          // the peak in the preconditioner's H-LU, over its size before it
const double kPreResident = 1.5;    // the preconditioner after its H-LU, over its size before it
const double kThreadBytes = 4e6;    // per thread
const double kDirectGrowth = 1.15;  // -hlu: the H-LU factors over the H-matrix
const double kDirectPeak = 1.7;     // -hlu: the peak over the factors (the products' temporaries)

MemEstimate estimate_memory( const MemPlan& p, int n, int nrhs, bool direct, double elt, double used ) {
   MemEstimate e;
   e.used = used;
   const double blocks = double( p.n_dense + p.n_lowrank + p.n_hblocks );
   const double per_thread = kThreadBytes * threads();
   // The deflation's rank-one term adds one to every low-rank leaf's rank.
   const double h = elt * ( p.dense_entries + p.lr_entries + p.lr_width );
   if ( direct ) {
      e.h = kDirectGrowth * h;
      // The right-hand sides: BEST's, and the tree-ordered copy solved in place.
      e.vectors = 8.0 * nrhs * 3.0 * n * 3;
      e.other = kBlockBytes * blocks + ( kDirectPeak - 1 ) * e.h + per_thread;
   } else {
      e.h = h;
      // Per right-hand side, GMRES keeps kRestart + 1 basis vectors and about 7 more (its own, the
      // tree-ordered copies, BEST's); each product also keeps V^T X for every low-rank leaf.
      e.vectors = 8.0 * nrhs * ( 3.0 * n * ( kRestart + 8 ) + p.lr_ranks + p.n_lowrank );
      // The preconditioner at the peak: in its H-LU, before the vectors exist, or with them in GMRES.
      const double pre = elt * p.pre_entries;
      e.pre = std::max( kPreLU * pre - e.vectors, kPreResident * pre );
      e.other = 2 * kBlockBytes * blocks + per_thread;
   }
   // The sample's entries are in use already, and part of the H-matrix.
   e.peak = e.used + e.h + e.pre + e.vectors + e.other - elt * p.sample_entries;
   return e;
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

// Frees the matrices after the memory guard stopped the solver, and says so.  Returns -1.
template <class Elt>
int stop_low( Solver& s, const char* where ) {
   Mats<Elt>& m = mats<Elt>( s );
   m.H.reset();
   m.P.reset();
   trim();
   say( "error: available memory was down to %.1f GB %s, below the floor of %.1f GB; the solver stopped there and "
        "freed its memory, so that the node does not run out (BEST_HMATRIX_MEMFLOOR sets the floor in GB, 0 turns "
        "this off)",
        gb( memory_low_seen() ), where, gb( s.floor ) );
   return -1;
}

// For calibrating the estimate (BEST_HMATRIX_MEMLOG=1): its parts next to the sizes built.
bool memlog() { return env_double( "BEST_HMATRIX_MEMLOG", 0 ) != 0; }

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

// -hlu: the H-LU factorization of the H-matrix at eps, in place.  Returns 0, 1 if it broke down,
// -1 if the memory guard stopped it (stop_low).
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
      if ( memory_low() ) {
         return stop_low<Elt>( s, "during the H-LU factorization" );
      }
      report_mesh( s );
      say( "error: H-LU: %s pivot at triangle %d; without -hlu BEST solves by GMRES, refining its preconditioner as "
           "needed",
           z.finite ? "a zero" : "a NaN or Inf", tri( s, z.row ) );
      return 1;
   }
   if ( memory_low() ) {
      return stop_low<Elt>( s, "during the H-LU factorization" );
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
// pivot was 0 or not finite (the factors are then unusable); -1 if the memory guard stopped it
// (stop_low).  With plan (the first one), BEST_HMATRIX_MEMLOG=1 logs its size against the estimate's.
template <class Elt>
int build_pre( Solver& s, const MemPlan* plan = nullptr ) {
   const int n = s.mesh.n;
   Mats<Elt>& m = mats<Elt>( s );
   m.P.reset();
   trim();
   reset_profile();
   const double cpu0 = cpu_seconds(), t = now();
   m.P = std::make_unique<HMatrixT<Elt>>( m.H->copy( s.eps_pre, true ) );
   if ( memory_low() ) {
      return stop_low<Elt>( s, "while the GMRES preconditioner was set up" );
   }
   if ( plan && memlog() ) {
      const HStats st = m.P->stats(), sh = m.H->stats();
      say( "memory check: preconditioner before its H-LU %.4g entries (estimate %.4g); H-matrix low-rank %.4g "
           "entries, ranks %.4g (estimates %.4g, %.4g)",
           st.bytes() / double( sizeof( Elt ) ), plan->pre_entries, sh.bytes_lowrank / double( sizeof( Elt ) ),
           sh.kmean * double( sh.n_lowrank ), plan->lr_entries + plan->lr_width,
           plan->lr_ranks + double( plan->n_lowrank ) );
      say_memory();  // no trim() here: the log must not change the run's memory
   }
   int st = 0;
   try {
      const double piv = m.P->lu( s.eps_pre );
      if ( memory_low() ) {
         return stop_low<Elt>( s, "during the GMRES preconditioner's H-LU factorization" );
      }
      say_stats( "GMRES preconditioner (H-LU factors)", *m.P, n, now() - t );
      say( "GMRES preconditioner: eps %.0e, near field compressed, smallest pivot ratio %.2e (triangle %d); GMRES(%d) "
           "to a relative residual of %.0e",
           s.eps_pre, piv, tri( s, HMatrixT<Elt>::pivot_row() ), s.restart, s.tol );
      if ( !( piv >= kPivotCollapse ) ) {
         say( "GMRES preconditioner: broke down (smallest pivot ratio below %.0e)", kPivotCollapse );
         st = 1;
      }
   } catch ( const ZeroPivot& z ) {
      if ( memory_low() ) {
         return stop_low<Elt>( s, "during the GMRES preconditioner's H-LU factorization" );
      }
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
   // Above the memory guard's floor, as the guard would stop it below.  Without /proc/meminfo there
   // is nothing to check against.
   const HMatrixT<Elt>* P = mats<Elt>( s ).P.get();
   const double free_now = double( mem_available() );
   const double avail = free_now + ( P ? P->stats().bytes() : 0.0 ) - s.floor;
   if ( free_now > 0 && need > avail && env_double( "BEST_HMATRIX_MEMCHECK", 1 ) != 0 ) {
      say( "GMRES preconditioner: one at eps %.0e would need about %.1f GB, %.1f GB available above the floor of "
           "%.1f GB (set BEST_HMATRIX_MEMCHECK=0 to try anyway)",
           next, gb( need ), gb( std::max( 0.0, avail ) ), gb( s.floor ) );
      return false;
   }
   say( "GMRES preconditioner: rebuilding it at eps %.0e", next );
   s.eps_pre = next;
   return true;
}

// bhm_build after the cluster tree, with the blocks in Elt.
template <class Elt>
int build_t( Solver& s, double eps, double eps_pre, const Deflation& def, const double* xmu, int nrhs ) {
   const int n = s.mesh.n, n3 = 3 * n;
   Mats<Elt>& m = mats<Elt>( s );
   const MemoryGuard guard( s.floor );
   const bool check = env_double( "BEST_HMATRIX_MEMCHECK", 1 ) != 0;
   if ( check && memory_low() ) {
      say( "error: not enough memory for the H-matrix solver: %.1f GB available, below the floor of %.1f GB "
           "(BEST_HMATRIX_MEMFLOOR sets the floor in GB)",
           gb( memory_low_seen() ), gb( s.floor ) );
      return 1;
   }
   HOptions opt;
   opt.eps = eps;
   opt.eta = kEta;
   long calls = 0;
   double t = now();
   m.H = std::make_unique<HMatrixT<Elt>>( HMatrixT<Elt>::structure( s.ct, kEta ) );
   std::vector<char> filled( m.H->leaves().size(), 0 );
   const MemPlan plan = plan_memory( *m.H, *s.K, opt, s.direct, eps_pre, filled, &calls );
   if ( memory_low() ) {
      stop_low<Elt>( s, "while the memory estimate's sample was computed" );
      return 1;
   }
   const double avail = double( mem_available() ), used = double( current_rss() );
   const MemEstimate est = estimate_memory( plan, n, nrhs, s.direct, double( sizeof( Elt ) ), used );
   say( "memory: about %.1f GB needed (estimate: %.1f in use, H-matrix%s %.1f, preconditioner %.1f, vectors %.1f, "
        "other %.1f), %.1f GB available",
        gb( est.peak ), gb( used ), s.direct ? " factors" : "", gb( est.h ), gb( est.pre ), gb( est.vectors ),
        gb( est.other ), gb( avail ) );
   say( "memory estimate: %ld dense and %ld low-rank blocks, their ranks from a sample of %ld blocks (%.1f s)",
        plan.n_dense, plan.n_lowrank, plan.n_sample, plan.sample_s );
   if ( memlog() ) {
      say( "memory estimate parts: dense %.4g, diagonal %.4g, near %.4g entries; low-rank width %.4g, entries %.4g, "
           "ranks %.4g; preconditioner entries %.4g; %ld hierarchical blocks; sample %.4g entries; in use %.4g GB; "
           "floor %.1f GB",
           plan.dense_entries, plan.diag_entries, plan.near_entries, plan.lr_width, plan.lr_entries, plan.lr_ranks,
           plan.pre_entries, plan.n_hblocks, plan.sample_entries, gb( used ), gb( s.floor ) );
   }
   if ( check && avail > 0 && est.peak - used > avail - s.floor ) {
      m.H.reset();
      trim();
      say( "error: not enough memory for the H-matrix solver at eps %.0e: about %.1f GB more needed, %.1f GB available "
           "above the floor of %.1f GB; %s (set BEST_HMATRIX_MEMCHECK=0 to try anyway: the solver then stops if "
           "available memory falls below the floor)",
           eps, gb( est.peak - used ), gb( std::max( 0.0, avail - s.floor ) ), gb( s.floor ),
           s.single ? "a larger -he eps needs less, but 1e-3 saves only about 20% over 1e-4, at about 50 times the error"
                    : "-hs (single precision) needs about half" );
      return 1;
   }
   std::vector<long> rest;
   rest.reserve( filled.size() );
   for ( size_t i = 0; i < filled.size(); ++i ) {
      if ( !filled[ i ] ) {
         rest.push_back( long( i ) );
      }
   }
   m.H->fill( *s.K, opt, rest, &calls );
   if ( memory_low() ) {
      stop_low<Elt>( s, "while the H-matrix was built" );
      return 1;
   }
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
      return factor_direct<Elt>( s ) != 0 ? 1 : 0;
   }
   // GMRES: the preconditioner, finer while it breaks down.
   s.eps_pre = eps_pre;
   for ( bool first = true;; first = false ) {
      const int st = build_pre<Elt>( s, first ? &plan : nullptr );
      if ( st < 0 ) {
         return 1;
      }
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
      if ( r.stopped ) {
         stop_low<Elt>( s, "during GMRES" );
         return 1;
      }
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
         if ( st < 0 ) {
            return 1;
         }
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
               int single, int nrhs ) {
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

      // The memory guard's floor: kFloorFrac of the node's memory, at least kFloorMin.
      s.floor = std::max( kFloorMin, kFloorFrac * double( mem_total() ) );
      if ( std::getenv( "BEST_HMATRIX_MEMFLOOR" ) ) {
         s.floor = 1e9 * env_double( "BEST_HMATRIX_MEMFLOOR", gb( s.floor ) );
      }
      nrhs = std::max( 1, nrhs );
      Deflation def;
      def.gamma.assign( gdef, gdef + std::max( 1, ndef ) );
      if ( ndef > 1 ) {
         def.comp.resize( n );
         for ( int p = 0; p < n; ++p ) {
            def.comp[ p ] = comp[ s.ct.perm[ p ] ] - 1;  // tree order, from 0
         }
      }
      return s.single ? build_t<float>( s, eps, eps_pre, def, xmu, nrhs )
                      : build_t<double>( s, eps, eps_pre, def, xmu, nrhs );
   } );
}

int bhm_solve( int nrhs, double* b ) {
   return guarded( [&] {
      Solver& s = *g;
      const int n3 = 3 * s.mesh.n;
      const MemoryGuard guard( s.floor );
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
