// BEST H-matrix solver: H-matrix build (ACA), matvec and H-LU.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "hmatrix.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <string>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <stdexcept>

#include "parallel.h"

namespace bh {

// ---------------------------------------------------------------------------------------------
// Structure and fill
// ---------------------------------------------------------------------------------------------

namespace {

std::unique_ptr<Block> make_blocks( const ClusterTree& ct, int t, int s, double eta ) {
   auto b = std::make_unique<Block>();
   const Cluster& ctl = ct.c[ t ];
   const Cluster& csl = ct.c[ s ];
   b->t = t;
   b->s = s;
   b->level = ctl.level;
   b->r0 = 3 * ctl.beg;
   b->nr = 3 * ctl.n;
   b->c0 = 3 * csl.beg;
   b->nc = 3 * csl.n;
   if ( t != s && admissible( ctl, csl, eta ) ) {
      b->type = BType::R;
   } else if ( ctl.is_leaf() || csl.is_leaf() ) {
      if ( !( ctl.is_leaf() && csl.is_leaf() ) ) {
         throw std::runtime_error( "cluster tree is not level-consistent" );
      }
      b->type = BType::D;
   } else {
      b->type = BType::H;
      for ( int i = 0; i < 2; ++i ) {
         for ( int j = 0; j < 2; ++j ) {
            b->ch[ 2 * i + j ] = make_blocks( ct, ctl.son[ i ], csl.son[ j ], eta );
         }
      }
   }
   return b;
}

void collect_leaves( Block* b, std::vector<Block*>& out ) {
   if ( b->type == BType::H ) {
      for ( auto& c : b->ch ) {
         collect_leaves( c.get(), out );
      }
   } else {
      out.push_back( b );
   }
}

// Adaptive cross approximation with partial pivoting.  Each step takes the three rows of one row
// triangle (one 3x3 kernel block per column triangle) and makes up to three rank-one terms from
// them; columns are computed three at a time too (one column triangle).  Stops after a row step
// in which every new term had |u| |v| <= eps |S|_F.  Returns false if the rank reached kmax.
bool aca( const Kernel& K, const int* rt, int m, const int* ct, int n, double eps, int kmax, Mat& U, Mat& V,
          long& calls ) {
   const int M = 3 * m, N = 3 * n;
   std::vector<double> us, vs;
   int k = 0;
   std::vector<char> row_used( m, 0 );
   std::vector<double> R( 3 * size_t( N ) ), u( M ), v( N ), cache;
   std::vector<int> cache_tri;
   double norm2 = 0, piv_max = 0;
   int i = m / 2, zero_rows = 0;
   double g[ 9 ];
   for ( ;; ) {
      for ( int jt = 0; jt < n; ++jt ) {
         K.block( rt[ i ], ct[ jt ], g );
         for ( int a = 0; a < 3; ++a ) {
            for ( int b = 0; b < 3; ++b ) {
               R[ size_t( a ) * N + 3 * jt + b ] = g[ 3 * a + b ];
            }
         }
      }
      calls += n;
      row_used[ i ] = 1;
      for ( int l = 0; l < k; ++l ) {
         const double* ul = &us[ size_t( l ) * M ];
         const double* vl = &vs[ size_t( l ) * N ];
         for ( int a = 0; a < 3; ++a ) {
            const double ua = ul[ 3 * i + a ];
            if ( ua != 0 ) {
               double* Ra = &R[ size_t( a ) * N ];
               for ( int j = 0; j < N; ++j ) {
                  Ra[ j ] -= ua * vl[ j ];
               }
            }
         }
      }
      cache_tri.clear();
      cache.clear();
      int added = 0;
      bool small = true;
      for ( int step = 0; step < 3 && k < kmax; ++step ) {
         int ab = -1, jb = 0;
         double best = 0;
         for ( int a = 0; a < 3; ++a ) {
            const double* Ra = &R[ size_t( a ) * N ];
            for ( int j = 0; j < N; ++j ) {
               if ( std::fabs( Ra[ j ] ) > best ) {
                  best = std::fabs( Ra[ j ] );
                  ab = a;
                  jb = j;
               }
            }
         }
         if ( ab < 0 || best <= 1e-15 * piv_max ) {
            break;
         }
         piv_max = std::max( piv_max, best );
         const int jt = jb / 3, bc = jb % 3;
         int slot = -1;
         for ( size_t q = 0; q < cache_tri.size(); ++q ) {
            if ( cache_tri[ q ] == jt ) {
               slot = int( q );
            }
         }
         if ( slot < 0 ) {
            slot = int( cache_tri.size() );
            cache_tri.push_back( jt );
            cache.resize( size_t( slot + 1 ) * 3 * M );
            for ( int it = 0; it < m; ++it ) {
               K.block( rt[ it ], ct[ jt ], g );
               for ( int c = 0; c < 3; ++c ) {
                  for ( int a = 0; a < 3; ++a ) {
                     cache[ ( size_t( slot ) * 3 + c ) * M + 3 * it + a ] = g[ 3 * a + c ];
                  }
               }
            }
            calls += m;
         }
         const double* col = &cache[ ( size_t( slot ) * 3 + bc ) * M ];
         for ( int r = 0; r < M; ++r ) {
            u[ r ] = col[ r ];
         }
         for ( int l = 0; l < k; ++l ) {
            const double vjl = vs[ size_t( l ) * N + jb ];
            if ( vjl != 0 ) {
               const double* ul = &us[ size_t( l ) * M ];
               for ( int r = 0; r < M; ++r ) {
                  u[ r ] -= vjl * ul[ r ];
               }
            }
         }
         const double piv = R[ size_t( ab ) * N + jb ];
         for ( int r = 0; r < M; ++r ) {
            u[ r ] /= piv;
         }
         for ( int j = 0; j < N; ++j ) {
            v[ j ] = R[ size_t( ab ) * N + j ];
         }
         for ( int a = 0; a < 3; ++a ) {
            const double ua = u[ 3 * i + a ];
            double* Ra = &R[ size_t( a ) * N ];
            for ( int j = 0; j < N; ++j ) {
               Ra[ j ] -= ua * v[ j ];
            }
         }
         double nu = 0, nv = 0, cross = 0;
         for ( int r = 0; r < M; ++r ) {
            nu += u[ r ] * u[ r ];
         }
         for ( int j = 0; j < N; ++j ) {
            nv += v[ j ] * v[ j ];
         }
         for ( int l = 0; l < k; ++l ) {
            const double* ul = &us[ size_t( l ) * M ];
            const double* vl = &vs[ size_t( l ) * N ];
            double du = 0, dv = 0;
            for ( int r = 0; r < M; ++r ) {
               du += ul[ r ] * u[ r ];
            }
            for ( int j = 0; j < N; ++j ) {
               dv += vl[ j ] * v[ j ];
            }
            cross += du * dv;
         }
         norm2 += 2 * cross + nu * nv;
         us.insert( us.end(), u.begin(), u.end() );
         vs.insert( vs.end(), v.begin(), v.end() );
         ++k;
         ++added;
         if ( std::sqrt( nu * nv ) > eps * std::sqrt( std::max( norm2, 0.0 ) ) ) {
            small = false;
         }
      }
      zero_rows = added == 0 ? zero_rows + 1 : 0;
      if ( added > 0 && small ) {
         break;
      }
      if ( k >= kmax ) {
         return false;
      }
      if ( zero_rows >= 3 ) {
         break;
      }
      int next = -1;
      double bu = 0;
      if ( k > 0 ) {
         const double* ul = &us[ size_t( k - 1 ) * M ];
         for ( int it = 0; it < m; ++it ) {
            if ( !row_used[ it ] ) {
               const double val =
                   std::max( std::fabs( ul[ 3 * it ] ), std::max( std::fabs( ul[ 3 * it + 1 ] ), std::fabs( ul[ 3 * it + 2 ] ) ) );
               if ( val > bu ) {
                  bu = val;
                  next = it;
               }
            }
         }
      }
      if ( next < 0 ) {
         for ( int it = 0; it < m; ++it ) {
            if ( !row_used[ it ] ) {
               next = it;
               break;
            }
         }
      }
      if ( next < 0 ) {
         break;
      }
      i = next;
   }
   U.resize( M, k );
   V.resize( N, k );
   std::copy( us.begin(), us.end(), U.a.begin() );
   std::copy( vs.begin(), vs.end(), V.a.begin() );
   return true;
}

void fill_leaf( Block* b, const ClusterTree& ct, const Kernel& K, const HOptions& opt, long& calls ) {
   const Cluster& t = ct.c[ b->t ];
   const Cluster& s = ct.c[ b->s ];
   const int* rt = ct.perm.data() + t.beg;
   const int* st = ct.perm.data() + s.beg;
   if ( b->type == BType::D ) {
      b->D.resize( b->nr, b->nc );
      K.fill( rt, t.n, st, s.n, b->D.data(), b->nr );
      calls += long( t.n ) * s.n;
      if ( opt.near && b->t != b->s ) {
         Mat P = b->D, U, V;
         dense_to_lowrank( P, opt.eps, U, V );
         if ( double( U.n ) * ( b->nr + b->nc ) < double( b->nr ) * b->nc ) {
            b->type = BType::R;
            b->U = std::move( U );
            b->V = std::move( V );
            b->D.release();
         }
      }
      return;
   }
   const int kmax = std::max( 3, std::min( b->nr, b->nc ) / 2 );
   bool done = false;
   if ( long( t.n ) * s.n > opt.dense_small ) {
      done = aca( K, rt, t.n, st, s.n, opt.aca_factor * opt.eps, kmax, b->U, b->V, calls );
      if ( done ) {
         truncate( b->U, b->V, opt.eps );
      }
   }
   if ( !done ) {
      Mat P( b->nr, b->nc );
      K.fill( rt, t.n, st, s.n, P.data(), b->nr );
      calls += long( t.n ) * s.n;
      dense_to_lowrank( P, opt.eps, b->U, b->V );
   }
}

}  // namespace

HMatrix::HMatrix( const ClusterTree& ct, std::unique_ptr<Block> root ) : ct_( &ct ), root_( std::move( root ) ) {
   index();
}

void HMatrix::index() {
   leaves_.clear();
   collect_leaves( root_.get(), leaves_ );
   rows_.assign( ct_->leaves.size(), {} );
   for ( size_t i = 0; i < leaves_.size(); ++i ) {
      const Cluster& t = ct_->c[ leaves_[ i ]->t ];
      for ( int l = t.leaf_first; l < t.leaf_first + t.leaf_count; ++l ) {
         rows_[ l ].push_back( int( i ) );
      }
   }
}

HMatrix HMatrix::build( const ClusterTree& ct, const Kernel& K, const HOptions& opt, long* kernel_calls ) {
   HMatrix h( ct, make_blocks( ct, 0, 0, opt.eta ) );
   // Largest blocks first, so the long ACA runs start early.
   std::vector<Block*> order = h.leaves_;
   std::stable_sort( order.begin(), order.end(),
                     []( const Block* a, const Block* b ) { return double( a->nr ) * a->nc > double( b->nr ) * b->nc; } );
   std::vector<long> calls( threads(), 0 );
   parallel_for( long( order.size() ), [&]( long i, int thread ) { fill_leaf( order[ i ], ct, K, opt, calls[ thread ] ); } );
   if ( kernel_calls ) {
      *kernel_calls = 0;
      for ( long c : calls ) {
         *kernel_calls += c;
      }
   }
   return h;
}

namespace {

std::unique_ptr<Block> copy_block( const Block* b ) {
   auto c = std::make_unique<Block>();
   c->r0 = b->r0;
   c->nr = b->nr;
   c->c0 = b->c0;
   c->nc = b->nc;
   c->t = b->t;
   c->s = b->s;
   c->level = b->level;
   c->type = b->type;
   if ( b->type == BType::H ) {
      for ( int i = 0; i < 4; ++i ) {
         c->ch[ i ] = copy_block( b->ch[ i ].get() );
      }
   }
   return c;  // the leaves' data is copied by HMatrix::copy, in parallel
}

}  // namespace

HMatrix HMatrix::copy( double eps, bool near ) const {
   HMatrix h( *ct_, copy_block( root_.get() ) );
   // Same structure, so the leaves are in the same order.
   parallel_for( long( leaves_.size() ), [&]( long i, int ) {
      const Block* s = leaves_[ i ];
      Block* d = h.leaves_[ i ];
      if ( s->type == BType::D ) {
         if ( near && eps > 0 && s->t != s->s ) {
            Mat P = s->D, U, V;
            dense_to_lowrank( P, eps, U, V );
            if ( double( U.n ) * ( s->nr + s->nc ) < double( s->nr ) * s->nc ) {
               d->type = BType::R;
               d->U = std::move( U );
               d->V = std::move( V );
               return;
            }
         }
         d->D = s->D;
      } else {
         d->U = s->U;
         d->V = s->V;
         if ( eps > 0 ) {
            truncate( d->U, d->V, eps );
         }
      }
   } );
   return h;
}

HStats HMatrix::stats() const {
   HStats st;
   st.kmax_level.assign( ct_->depth + 1, 0 );
   st.nlr_level.assign( ct_->depth + 1, 0 );
   double ksum = 0;
   for ( const Block* b : leaves_ ) {
      if ( b->type == BType::D ) {
         ++st.n_dense;
         st.bytes_dense += double( b->nr ) * b->nc * sizeof( double );
      } else {
         ++st.n_lowrank;
         const int k = b->rank();
         st.bytes_lowrank += double( k ) * ( b->nr + b->nc ) * sizeof( double );
         st.kmax = std::max( st.kmax, k );
         ksum += k;
         st.kmax_level[ b->level ] = std::max( st.kmax_level[ b->level ], k );
         ++st.nlr_level[ b->level ];
      }
   }
   st.kmean = st.n_lowrank ? ksum / st.n_lowrank : 0;
   return st;
}

// ---------------------------------------------------------------------------------------------
// Matvec
// ---------------------------------------------------------------------------------------------

void HMatrix::matvec( const double* X, int ldx, int p, double* Y, int ldy ) const {
   // V^T X for every low-rank leaf first, then each leaf cluster's rows, adding the blocks that
   // cover them in a fixed order.
   std::vector<Mat> Z( leaves_.size() );
   parallel_for( long( leaves_.size() ), [&]( long i, int ) {
      const Block* b = leaves_[ i ];
      if ( b->type == BType::R && b->rank() > 0 ) {
         Z[ i ].resize( b->rank(), p );
         gemm( 'T', 'N', b->rank(), p, b->nc, 1.0, b->V.data(), b->nc, X + b->c0, ldx, 0.0, Z[ i ].data(), b->rank() );
      }
   } );
   const ClusterTree& ct = *ct_;
   parallel_for( long( ct.leaves.size() ), [&]( long l, int ) {
      const Cluster& cl = ct.c[ ct.leaves[ l ] ];
      const int r0 = 3 * cl.beg, nr = 3 * cl.n;
      for ( int j = 0; j < p; ++j ) {
         for ( int r = 0; r < nr; ++r ) {
            Y[ size_t( j ) * ldy + r0 + r ] = 0.0;
         }
      }
      for ( int i : rows_[ l ] ) {
         const Block* b = leaves_[ i ];
         const int off = r0 - b->r0;
         if ( b->type == BType::D ) {
            gemm( 'N', 'N', nr, p, b->nc, 1.0, b->D.data() + off, b->nr, X + b->c0, ldx, 1.0, Y + r0, ldy );
         } else if ( b->rank() > 0 ) {
            gemm( 'N', 'N', nr, p, b->rank(), 1.0, b->U.data() + off, b->nr, Z[ i ].data(), b->rank(), 1.0, Y + r0, ldy );
         }
      }
   }, 4 );
}

// ---------------------------------------------------------------------------------------------
// H-LU
// ---------------------------------------------------------------------------------------------

namespace {

// Truncations by call site: count and nanoseconds (all threads).
const int kNSites = 5;
const char* const kSites[ kNSites ] = { "upd", "quarter", "merge", "dd", "pend" };
std::atomic<long> g_site_n[ kNSites ], g_site_ns[ kNSites ];

struct SiteTimer {
   int s;
   std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
   explicit SiteTimer( int site ) : s( site ) {}
   ~SiteTimer() {
      ++g_site_n[ s ];
      g_site_ns[ s ] += long( std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now() - t0 ).count() );
   }
};

double g_fork_min = 0;          // fork only for blocks with at least this many entries
int g_fork_leaves = 8;          // ... (g_fork_leaves leaf clusters)^2
std::mutex g_piv_m;
double g_piv_min = 1e300, g_piv_max = 0;

template <class A, class B>
void par2( bool doit, A&& a, B&& b ) {
   if ( doit ) {
      fork_join( a, b );
   } else {
      a();
      b();
   }
}

bool big( const Block* b ) { return double( b->nr ) * double( b->nc ) >= g_fork_min; }

// Y (nr x p) += alpha A X (nc x p).
void mul_block_dense( const Block* A, double alpha, const double* X, int ldx, int p, double* Y, int ldy ) {
   switch ( A->type ) {
      case BType::D:
         gemm( 'N', 'N', A->nr, p, A->nc, alpha, A->D.data(), A->nr, X, ldx, 1.0, Y, ldy );
         break;
      case BType::R: {
         const int k = A->rank();
         if ( k == 0 ) {
            return;
         }
         Mat T( k, p );
         gemm( 'T', 'N', k, p, A->nc, 1.0, A->V.data(), A->nc, X, ldx, 0.0, T.data(), k );
         gemm( 'N', 'N', A->nr, p, k, alpha, A->U.data(), A->nr, T.data(), k, 1.0, Y, ldy );
         break;
      }
      case BType::H:
         for ( int i = 0; i < 2; ++i ) {
            for ( int j = 0; j < 2; ++j ) {
               const Block* c = A->c( i, j );
               mul_block_dense( c, alpha, X + ( c->c0 - A->c0 ), ldx, p, Y + ( c->r0 - A->r0 ), ldy );
            }
         }
         break;
   }
}

// Y (p x nc) += alpha X (p x nr) A.
void mul_dense_block( double alpha, const double* X, int ldx, int p, const Block* A, double* Y, int ldy ) {
   switch ( A->type ) {
      case BType::D:
         gemm( 'N', 'N', p, A->nc, A->nr, alpha, X, ldx, A->D.data(), A->nr, 1.0, Y, ldy );
         break;
      case BType::R: {
         const int k = A->rank();
         if ( k == 0 ) {
            return;
         }
         Mat T( p, k );
         gemm( 'N', 'N', p, k, A->nr, 1.0, X, ldx, A->U.data(), A->nr, 0.0, T.data(), p );
         gemm( 'N', 'T', p, A->nc, k, alpha, T.data(), p, A->V.data(), A->nc, 1.0, Y, ldy );
         break;
      }
      case BType::H:
         for ( int i = 0; i < 2; ++i ) {
            for ( int j = 0; j < 2; ++j ) {
               const Block* c = A->c( i, j );
               mul_dense_block( alpha, X + size_t( c->r0 - A->r0 ) * ldx, ldx, p, c,
                                Y + size_t( c->c0 - A->c0 ) * ldy, ldy );
            }
         }
         break;
   }
}

// X (nr x p) <- L^-1 X, L unit lower triangular (a diagonal block after lu).
void trsm_ll_dense( const Block* L, double* X, int ldx, int p ) {
   if ( L->type == BType::D ) {
      trsm( 'L', 'L', 'N', 'U', L->nr, p, 1.0, L->D.data(), L->nr, X, ldx );
      return;
   }
   const int n0 = L->c( 0, 0 )->nr;
   trsm_ll_dense( L->c( 0, 0 ), X, ldx, p );
   mul_block_dense( L->c( 1, 0 ), -1.0, X, ldx, p, X + n0, ldx );
   trsm_ll_dense( L->c( 1, 1 ), X + n0, ldx, p );
}

// X (nr x p) <- U^-1 X, U upper triangular.
void trsm_ul_dense( const Block* U, double* X, int ldx, int p ) {
   if ( U->type == BType::D ) {
      trsm( 'L', 'U', 'N', 'N', U->nr, p, 1.0, U->D.data(), U->nr, X, ldx );
      return;
   }
   const int n0 = U->c( 0, 0 )->nr;
   trsm_ul_dense( U->c( 1, 1 ), X + n0, ldx, p );
   mul_block_dense( U->c( 0, 1 ), -1.0, X + n0, ldx, p, X, ldx );
   trsm_ul_dense( U->c( 0, 0 ), X, ldx, p );
}

// X (p x nc) <- X U^-1, U upper triangular.
void trsm_ur_dense( const Block* U, double* X, int ldx, int p ) {
   if ( U->type == BType::D ) {
      trsm( 'R', 'U', 'N', 'N', p, U->nc, 1.0, U->D.data(), U->nc, X, ldx );
      return;
   }
   const int n0 = U->c( 0, 0 )->nc;
   trsm_ur_dense( U->c( 0, 0 ), X, ldx, p );
   mul_dense_block( -1.0, X, ldx, p, U->c( 0, 1 ), X + size_t( n0 ) * ldx, ldx );
   trsm_ur_dense( U->c( 1, 1 ), X + size_t( n0 ) * ldx, ldx, p );
}

Mat transpose( const Mat& A ) {
   Mat T( A.n, A.m );
   for ( int j = 0; j < A.n; ++j ) {
      for ( int i = 0; i < A.m; ++i ) {
         T( j, i ) = A( i, j );
      }
   }
   return T;
}

// C (R) <- C + alpha U V^T, truncated.  U: C->nr rows (ldu), V: C->nc rows (ldv), k columns.
void add_lowrank_r( Block* C, double alpha, const double* U, int ldu, const double* V, int ldv, int k, double eps ) {
   if ( k == 0 ) {
      return;
   }
   const int kc = C->rank(), m = C->nr, n = C->nc;
   Mat Un( m, kc + k ), Vn( n, kc + k );
   std::copy( C->U.a.begin(), C->U.a.end(), Un.a.begin() );
   std::copy( C->V.a.begin(), C->V.a.end(), Vn.a.begin() );
   for ( int j = 0; j < k; ++j ) {
      for ( int i = 0; i < m; ++i ) {
         Un( i, kc + j ) = alpha * U[ size_t( j ) * ldu + i ];
      }
      for ( int i = 0; i < n; ++i ) {
         Vn( i, kc + j ) = V[ size_t( j ) * ldv + i ];
      }
   }
   {
      SiteTimer st( 0 );
      truncate( Un, Vn, eps );
   }
   C->U = std::move( Un );
   C->V = std::move( Vn );
}

// C (any type) += alpha U V^T.
void add_lowrank( Block* C, double alpha, const double* U, int ldu, const double* V, int ldv, int k, double eps ) {
   if ( k == 0 ) {
      return;
   }
   switch ( C->type ) {
      case BType::D:
         gemm( 'N', 'T', C->nr, C->nc, k, alpha, U, ldu, V, ldv, 1.0, C->D.data(), C->nr );
         break;
      case BType::R:
         add_lowrank_r( C, alpha, U, ldu, V, ldv, k, eps );
         break;
      case BType::H: {
         // Combine with the pending update at this level; flush() pushes it down later.
         const int kp = C->Pu.n, m = C->nr, n = C->nc;
         Mat Un( m, kp + k ), Vn( n, kp + k );
         std::copy( C->Pu.a.begin(), C->Pu.a.end(), Un.a.begin() );
         std::copy( C->Pv.a.begin(), C->Pv.a.end(), Vn.a.begin() );
         for ( int j = 0; j < k; ++j ) {
            for ( int i = 0; i < m; ++i ) {
               Un( i, kp + j ) = alpha * U[ size_t( j ) * ldu + i ];
            }
            for ( int i = 0; i < n; ++i ) {
               Vn( i, kp + j ) = V[ size_t( j ) * ldv + i ];
            }
         }
         if ( kp > 0 ) {
            SiteTimer st( 4 );
            truncate( Un, Vn, eps );
         }
         C->Pu = std::move( Un );
         C->Pv = std::move( Vn );
         break;
      }
   }
}

// Pushes an H block's pending update down to its children (children that are H keep it pending).
void flush( Block* C, double eps ) {
   if ( C->type != BType::H || C->Pu.n == 0 ) {
      return;
   }
   Mat U = std::move( C->Pu ), V = std::move( C->Pv );
   C->Pu.release();
   C->Pv.release();
   auto child = [&]( int i, int j ) {
      Block* c = C->c( i, j );
      add_lowrank( c, 1.0, U.data() + ( c->r0 - C->r0 ), U.m, V.data() + ( c->c0 - C->c0 ), V.m, U.n, eps );
   };
   const bool f = big( C );
   par2( f, [&] { par2( f, [&] { child( 0, 0 ); }, [&] { child( 0, 1 ); } ); },
         [&] { par2( f, [&] { child( 1, 0 ); }, [&] { child( 1, 1 ); } ); } );
}

// A B ~ U V^T, A (r,k) and B (k,c) at the same level.
void product_lowrank( const Block* A, const Block* B, double eps, Mat& U, Mat& V ) {
   const int nr = A->nr, nc = B->nc;
   if ( A->type == BType::R ) {
      const int ka = A->rank();
      U = A->U;
      V.resize( nc, ka );
      if ( ka == 0 ) {
         return;
      }
      Mat VaT = transpose( A->V );  // ka x |k|
      Mat T( ka, nc );
      mul_dense_block( 1.0, VaT.data(), ka, ka, B, T.data(), ka );
      V = transpose( T );
      return;
   }
   if ( B->type == BType::R ) {
      const int kb = B->rank();
      V = B->V;
      U.resize( nr, kb );
      if ( kb == 0 ) {
         return;
      }
      mul_block_dense( A, 1.0, B->U.data(), B->nr, kb, U.data(), nr );
      return;
   }
   if ( A->type == BType::D && B->type == BType::D ) {
      Mat P( nr, nc );
      gemm( 'N', 'N', nr, nc, A->nc, 1.0, A->D.data(), nr, B->D.data(), B->nr, 0.0, P.data(), nr );
      SiteTimer st( 3 );
      dense_to_lowrank( P, eps, U, V );
      return;
   }
   if ( A->type != BType::H || B->type != BType::H ) {
      throw std::runtime_error( "product_lowrank: mixed dense and hierarchical blocks" );
   }
   // Each quarter of the product, then the four together.
   Mat Uq[ 4 ], Vq[ 4 ];
   auto quarter = [&]( int i, int j ) {
      Mat& Ui = Uq[ 2 * i + j ];
      Mat& Vi = Vq[ 2 * i + j ];
      Mat U1, V1;
      product_lowrank( A->c( i, 0 ), B->c( 0, j ), eps, Ui, Vi );
      product_lowrank( A->c( i, 1 ), B->c( 1, j ), eps, U1, V1 );
      if ( U1.n > 0 ) {
         const int k0 = Ui.n;
         Mat Un( Ui.m, k0 + U1.n ), Vn( Vi.m, k0 + V1.n );
         std::copy( Ui.a.begin(), Ui.a.end(), Un.a.begin() );
         std::copy( U1.a.begin(), U1.a.end(), Un.a.begin() + Ui.a.size() );
         std::copy( Vi.a.begin(), Vi.a.end(), Vn.a.begin() );
         std::copy( V1.a.begin(), V1.a.end(), Vn.a.begin() + Vi.a.size() );
         {
            SiteTimer st( 1 );
            truncate( Un, Vn, eps );
         }
         Ui = std::move( Un );
         Vi = std::move( Vn );
      }
   };
   const bool f = double( nr ) * nc >= g_fork_min;
   par2( f, [&] { par2( f, [&] { quarter( 0, 0 ); }, [&] { quarter( 0, 1 ); } ); },
         [&] { par2( f, [&] { quarter( 1, 0 ); }, [&] { quarter( 1, 1 ); } ); } );
   int ktot = 0;
   for ( int q = 0; q < 4; ++q ) {
      ktot += Uq[ q ].n;
   }
   U.resize( nr, ktot );
   V.resize( nc, ktot );
   int col = 0;
   for ( int i = 0; i < 2; ++i ) {
      for ( int j = 0; j < 2; ++j ) {
         const Mat& Ui = Uq[ 2 * i + j ];
         const Mat& Vi = Vq[ 2 * i + j ];
         const int roff = A->c( i, 0 )->r0 - A->r0;
         const int coff = B->c( 0, j )->c0 - B->c0;
         for ( int l = 0; l < Ui.n; ++l, ++col ) {
            for ( int r = 0; r < Ui.m; ++r ) {
               U( roff + r, col ) = Ui( r, l );
            }
            for ( int r = 0; r < Vi.m; ++r ) {
               V( coff + r, col ) = Vi( r, l );
            }
         }
      }
   }
   SiteTimer st( 2 );
   truncate( U, V, eps );
}

// C += alpha A B, A (r,k), B (k,c), C (r,c), all at one level.
void addmul( Block* C, double alpha, const Block* A, const Block* B, double eps ) {
   if ( C->type == BType::D ) {
      if ( A->type == BType::R ) {
         const int ka = A->rank();
         if ( ka == 0 ) {
            return;
         }
         Mat VaT = transpose( A->V );
         Mat T( ka, C->nc );
         mul_dense_block( 1.0, VaT.data(), ka, ka, B, T.data(), ka );
         gemm( 'N', 'N', C->nr, C->nc, ka, alpha, A->U.data(), A->nr, T.data(), ka, 1.0, C->D.data(), C->nr );
      } else if ( B->type == BType::R ) {
         const int kb = B->rank();
         if ( kb == 0 ) {
            return;
         }
         Mat T( C->nr, kb );
         mul_block_dense( A, 1.0, B->U.data(), B->nr, kb, T.data(), C->nr );
         gemm( 'N', 'T', C->nr, C->nc, kb, alpha, T.data(), C->nr, B->V.data(), B->nc, 1.0, C->D.data(), C->nr );
      } else if ( A->type == BType::D && B->type == BType::D ) {
         gemm( 'N', 'N', C->nr, C->nc, A->nc, alpha, A->D.data(), A->nr, B->D.data(), B->nr, 1.0, C->D.data(), C->nr );
      } else {
         throw std::runtime_error( "addmul: hierarchical factor for a dense block" );
      }
      return;
   }
   if ( C->type == BType::R || A->type == BType::R || B->type == BType::R ) {
      Mat U, V;
      product_lowrank( A, B, eps, U, V );
      add_lowrank( C, alpha, U.data(), U.m, V.data(), V.m, U.n, eps );
      return;
   }
   if ( A->type != BType::H || B->type != BType::H ) {
      throw std::runtime_error( "addmul: dense factor for a hierarchical block" );
   }
   auto job = [&]( int i, int j ) {
      addmul( C->c( i, j ), alpha, A->c( i, 0 ), B->c( 0, j ), eps );
      addmul( C->c( i, j ), alpha, A->c( i, 1 ), B->c( 1, j ), eps );
   };
   const bool f = big( C );
   par2( f, [&] { par2( f, [&] { job( 0, 0 ); }, [&] { job( 0, 1 ); } ); },
         [&] { par2( f, [&] { job( 1, 0 ); }, [&] { job( 1, 1 ); } ); } );
}

// B <- L^-1 B, L (t,t) unit lower, B (t,s).
void trsm_ll( const Block* L, Block* B, double eps ) {
   flush( B, eps );
   switch ( B->type ) {
      case BType::R:
         if ( B->rank() > 0 ) {
            trsm_ll_dense( L, B->U.data(), B->nr, B->rank() );
         }
         break;
      case BType::D:
         trsm_ll_dense( L, B->D.data(), B->nr, B->nc );
         break;
      case BType::H: {
         auto col = [&]( int j ) {
            trsm_ll( L->c( 0, 0 ), B->c( 0, j ), eps );
            addmul( B->c( 1, j ), -1.0, L->c( 1, 0 ), B->c( 0, j ), eps );
            trsm_ll( L->c( 1, 1 ), B->c( 1, j ), eps );
         };
         par2( big( B ), [&] { col( 0 ); }, [&] { col( 1 ); } );
         break;
      }
   }
}

// B <- B U^-1, U (s,s) upper, B (t,s).
void trsm_ur( const Block* U, Block* B, double eps ) {
   flush( B, eps );
   switch ( B->type ) {
      case BType::R:
         if ( B->rank() > 0 ) {
            Mat Y = transpose( B->V );  // k x nc
            trsm_ur_dense( U, Y.data(), Y.m, Y.m );
            B->V = transpose( Y );
         }
         break;
      case BType::D:
         trsm_ur_dense( U, B->D.data(), B->nr, B->nr );
         break;
      case BType::H: {
         auto row = [&]( int i ) {
            trsm_ur( U->c( 0, 0 ), B->c( i, 0 ), eps );
            addmul( B->c( i, 1 ), -1.0, B->c( i, 0 ), U->c( 0, 1 ), eps );
            trsm_ur( U->c( 1, 1 ), B->c( i, 1 ), eps );
         };
         par2( big( B ), [&] { row( 0 ); }, [&] { row( 1 ); } );
         break;
      }
   }
}

void getrf_nopiv( Mat& A ) {
   const int n = A.m;
   double pmin = 1e300, pmax = 0;
   for ( int k = 0; k < n; ++k ) {
      const double p = A( k, k );
      pmin = std::min( pmin, std::fabs( p ) );
      pmax = std::max( pmax, std::fabs( p ) );
      if ( p == 0 ) {
         throw std::runtime_error( "H-LU: zero pivot" );
      }
      const double inv = 1.0 / p;
      for ( int i = k + 1; i < n; ++i ) {
         A( i, k ) *= inv;
      }
      for ( int j = k + 1; j < n; ++j ) {
         const double akj = A( k, j );
         if ( akj != 0 ) {
            double* cj = A.col( j );
            const double* ck = A.col( k );
            for ( int i = k + 1; i < n; ++i ) {
               cj[ i ] -= ck[ i ] * akj;
            }
         }
      }
   }
   std::lock_guard<std::mutex> lk( g_piv_m );
   g_piv_min = std::min( g_piv_min, pmin );
   g_piv_max = std::max( g_piv_max, pmax );
}

void hlu_block( Block* A, double eps ) {
   flush( A, eps );
   if ( A->type == BType::D ) {
      getrf_nopiv( A->D );
      return;
   }
   if ( A->type != BType::H ) {
      throw std::runtime_error( "H-LU: low-rank diagonal block" );
   }
   hlu_block( A->c( 0, 0 ), eps );
   par2( big( A->c( 0, 1 ) ), [&] { trsm_ll( A->c( 0, 0 ), A->c( 0, 1 ), eps ); },
         [&] { trsm_ur( A->c( 0, 0 ), A->c( 1, 0 ), eps ); } );
   addmul( A->c( 1, 1 ), -1.0, A->c( 1, 0 ), A->c( 0, 1 ), eps );
   hlu_block( A->c( 1, 1 ), eps );
}

}  // namespace

void HMatrix::set_fork( int leaves ) { g_fork_leaves = std::max( 1, leaves ); }

std::string HMatrix::profile( bool reset ) {
   std::string s;
   for ( int i = 0; i < kNSites; ++i ) {
      char buf[ 96 ];
      std::snprintf( buf, sizeof buf, "%sn_%s=%ld cpu_%s=%.3f", i ? " " : "", kSites[ i ], g_site_n[ i ].load(), kSites[ i ],
                     1e-9 * double( g_site_ns[ i ].load() ) );
      s += buf;
      if ( reset ) {
         g_site_n[ i ] = 0;
         g_site_ns[ i ] = 0;
      }
   }
   return s;
}

double HMatrix::lu( double eps ) {
   // Fork for blocks of at least g_fork_leaves x g_fork_leaves leaf-cluster sizes.
   const double leafdofs = 3.0 * ct_->c[ ct_->leaves[ 0 ] ].n;
   g_fork_min = ( g_fork_leaves * leafdofs ) * ( g_fork_leaves * leafdofs );
   g_piv_min = 1e300;
   g_piv_max = 0;
   hlu_block( root_.get(), eps );
   // Every H block was flushed when the factorization reached it; check, and flush any left.
   long left = 0;
   std::vector<Block*> stack = { root_.get() };
   while ( !stack.empty() ) {
      Block* b = stack.back();
      stack.pop_back();
      if ( b->type == BType::H ) {
         if ( b->Pu.n > 0 ) {
            ++left;
            flush( b, eps );
         }
         for ( auto& c : b->ch ) {
            stack.push_back( c.get() );
         }
      }
   }
   if ( left > 0 ) {
      std::fprintf( stderr, "H-LU: %ld blocks still had pending updates\n", left );
   }
   return g_piv_max > 0 ? g_piv_min / g_piv_max : 0;
}

void HMatrix::add_rank1( double gamma, const double* u, const double* v ) {
   parallel_for( long( leaves_.size() ), [&]( long i, int ) {
      Block* b = leaves_[ i ];
      const double* ur = u + b->r0;
      const double* vc = v + b->c0;
      if ( b->type == BType::D ) {
         for ( int j = 0; j < b->nc; ++j ) {
            for ( int r = 0; r < b->nr; ++r ) {
               b->D( r, j ) += gamma * ur[ r ] * vc[ j ];
            }
         }
      } else {
         const int k = b->rank();
         Mat Un( b->nr, k + 1 ), Vn( b->nc, k + 1 );
         std::copy( b->U.a.begin(), b->U.a.end(), Un.a.begin() );
         std::copy( b->V.a.begin(), b->V.a.end(), Vn.a.begin() );
         for ( int r = 0; r < b->nr; ++r ) {
            Un( r, k ) = gamma * ur[ r ];
         }
         for ( int r = 0; r < b->nc; ++r ) {
            Vn( r, k ) = vc[ r ];
         }
         b->U = std::move( Un );
         b->V = std::move( Vn );
      }
   } );
}

void HMatrix::lu_solve( double* X, int ldx, int p ) const {
   parallel_for( p, [&]( long j, int ) {
      trsm_ll_dense( root_.get(), X + size_t( j ) * ldx, ldx, 1 );
      trsm_ul_dense( root_.get(), X + size_t( j ) * ldx, ldx, 1 );
   } );
}

}  // namespace bh
