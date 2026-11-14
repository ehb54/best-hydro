// BEST H-matrix solver: small dense kernels (see dense.h).
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// gemm packs op(A) and op(B) into thread-local panels and runs an 8 x 4 register-blocked kernel.
// On x86-64 the kernel is compiled twice, for AVX2 with FMA and for the baseline, and the first is
// used when the processor has both.  Products with few columns (the matrix-vector products) skip
// the packing.  Compile with -O3 so that the kernels are vectorized.
#include "dense.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <numeric>

namespace bh {

namespace {

constexpr int kMR = 8, kNR = 4;              // register block
constexpr int kKC = 256, kMC = 128, kNC = 1024;  // cache blocks

thread_local std::vector<double> t_pack_a, t_pack_b;

// op(A)[ i0 .. i0 + mc, l0 .. l0 + kc ) into strips of kMR rows: strip s, then l, then i.
void pack_a( char ta, const double* A, int lda, int i0, int mc, int l0, int kc, double* p ) {
   for ( int s = 0; s < mc; s += kMR ) {
      const int mr = std::min( kMR, mc - s );
      if ( ta == 'N' ) {
         for ( int l = 0; l < kc; ++l, p += kMR ) {
            const double* a = A + size_t( l0 + l ) * lda + i0 + s;
            int i = 0;
            for ( ; i < mr; ++i ) {
               p[ i ] = a[ i ];
            }
            for ( ; i < kMR; ++i ) {
               p[ i ] = 0.0;
            }
         }
      } else {
         for ( int l = 0; l < kc; ++l, p += kMR ) {
            int i = 0;
            for ( ; i < mr; ++i ) {
               p[ i ] = A[ size_t( i0 + s + i ) * lda + l0 + l ];
            }
            for ( ; i < kMR; ++i ) {
               p[ i ] = 0.0;
            }
         }
      }
   }
}

// op(B)[ l0 .. l0 + kc, j0 .. j0 + nc ) into strips of kNR columns: strip s, then l, then j.
void pack_b( char tb, const double* B, int ldb, int l0, int kc, int j0, int nc, double* p ) {
   for ( int s = 0; s < nc; s += kNR ) {
      const int nr = std::min( kNR, nc - s );
      for ( int l = 0; l < kc; ++l, p += kNR ) {
         int j = 0;
         if ( tb == 'N' ) {
            for ( ; j < nr; ++j ) {
               p[ j ] = B[ size_t( j0 + s + j ) * ldb + l0 + l ];
            }
         } else {
            const double* b = B + size_t( l0 + l ) * ldb + j0 + s;
            for ( ; j < nr; ++j ) {
               p[ j ] = b[ j ];
            }
         }
         for ( ; j < kNR; ++j ) {
            p[ j ] = 0.0;
         }
      }
   }
}

// C[ mc x nc ] += alpha (packed A) (packed B).
inline __attribute__( ( always_inline ) ) void macro_body( int mc, int nc, int kc, double alpha, const double* pa,
                                                            const double* pb, double* C, int ldc ) {
   for ( int jb = 0; jb < nc; jb += kNR ) {
      const int nr = std::min( kNR, nc - jb );
      const double* bp = pb + size_t( jb / kNR ) * kc * kNR;
      for ( int ib = 0; ib < mc; ib += kMR ) {
         const int mr = std::min( kMR, mc - ib );
         const double* ap = pa + size_t( ib / kMR ) * kc * kMR;
         double acc[ kNR ][ kMR ] = {};
         for ( int l = 0; l < kc; ++l ) {
            const double* a = ap + size_t( l ) * kMR;
            const double* b = bp + size_t( l ) * kNR;
            for ( int j = 0; j < kNR; ++j ) {
               const double bj = b[ j ];
               for ( int i = 0; i < kMR; ++i ) {
                  acc[ j ][ i ] += a[ i ] * bj;
               }
            }
         }
         for ( int j = 0; j < nr; ++j ) {
            double* c = C + size_t( jb + j ) * ldc + ib;
            for ( int i = 0; i < mr; ++i ) {
               c[ i ] += alpha * acc[ j ][ i ];
            }
         }
      }
   }
}

// C[ m x nj ] += alpha op(A) op(B)[ :, 0 .. nj ) for nj <= kNR columns, A not packed.  bp: op(B)
// packed as one kNR-column strip (pack_b).  'N': kMR-row strips of A against all nj columns at
// once; 'T': for each column of A, nj dot products at once.
inline __attribute__( ( always_inline ) ) void thin_body( char ta, int m, int nj, int k, double alpha, const double* A,
                                                           int lda, const double* bp, double* C, int ldc ) {
   if ( ta == 'N' ) {
      int ib = 0;
      for ( ; ib + kMR <= m; ib += kMR ) {
         double acc[ kNR ][ kMR ] = {};
         for ( int l = 0; l < k; ++l ) {
            const double* a = A + size_t( l ) * lda + ib;
            const double* b = bp + size_t( l ) * kNR;
            for ( int j = 0; j < kNR; ++j ) {
               const double bj = b[ j ];
               for ( int i = 0; i < kMR; ++i ) {
                  acc[ j ][ i ] += a[ i ] * bj;
               }
            }
         }
         for ( int j = 0; j < nj; ++j ) {
            double* c = C + size_t( j ) * ldc + ib;
            for ( int i = 0; i < kMR; ++i ) {
               c[ i ] += alpha * acc[ j ][ i ];
            }
         }
      }
      for ( ; ib < m; ++ib ) {
         double acc[ kNR ] = {};
         for ( int l = 0; l < k; ++l ) {
            const double a = A[ size_t( l ) * lda + ib ];
            for ( int j = 0; j < kNR; ++j ) {
               acc[ j ] += a * bp[ size_t( l ) * kNR + j ];
            }
         }
         for ( int j = 0; j < nj; ++j ) {
            C[ size_t( j ) * ldc + ib ] += alpha * acc[ j ];
         }
      }
   } else {
      for ( int i = 0; i < m; ++i ) {
         const double* a = A + size_t( i ) * lda;
         double acc[ kNR ][ 4 ] = {};
         int l = 0;
         for ( ; l + 4 <= k; l += 4 ) {
            const double* b = bp + size_t( l ) * kNR;
            for ( int j = 0; j < kNR; ++j ) {
               for ( int u = 0; u < 4; ++u ) {
                  acc[ j ][ u ] += a[ l + u ] * b[ u * kNR + j ];
               }
            }
         }
         for ( ; l < k; ++l ) {
            for ( int j = 0; j < kNR; ++j ) {
               acc[ j ][ 0 ] += a[ l ] * bp[ size_t( l ) * kNR + j ];
            }
         }
         for ( int j = 0; j < nj; ++j ) {
            C[ size_t( j ) * ldc + i ] += alpha * ( ( acc[ j ][ 0 ] + acc[ j ][ 1 ] ) + ( acc[ j ][ 2 ] + acc[ j ][ 3 ] ) );
         }
      }
   }
}

using MacroFn = void ( * )( int, int, int, double, const double*, const double*, double*, int );
using ThinFn = void ( * )( char, int, int, int, double, const double*, int, const double*, double*, int );

void macro_base( int mc, int nc, int kc, double alpha, const double* pa, const double* pb, double* C, int ldc ) {
   macro_body( mc, nc, kc, alpha, pa, pb, C, ldc );
}
void thin_base( char ta, int m, int nj, int k, double alpha, const double* A, int lda, const double* bp, double* C,
                int ldc ) {
   thin_body( ta, m, nj, k, alpha, A, lda, bp, C, ldc );
}

#if defined( __x86_64__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
__attribute__( ( target( "avx2,fma" ) ) ) void macro_avx2( int mc, int nc, int kc, double alpha, const double* pa,
                                                           const double* pb, double* C, int ldc ) {
   macro_body( mc, nc, kc, alpha, pa, pb, C, ldc );
}
__attribute__( ( target( "avx2,fma" ) ) ) void thin_avx2( char ta, int m, int nj, int k, double alpha, const double* A,
                                                          int lda, const double* bp, double* C, int ldc ) {
   thin_body( ta, m, nj, k, alpha, A, lda, bp, C, ldc );
}
bool has_avx2() {
   __builtin_cpu_init();
   return __builtin_cpu_supports( "avx2" ) && __builtin_cpu_supports( "fma" );
}
const bool g_avx2 = has_avx2();
const MacroFn g_macro = g_avx2 ? macro_avx2 : macro_base;
const ThinFn g_thin = g_avx2 ? thin_avx2 : thin_base;
#else
const MacroFn g_macro = macro_base;
const ThinFn g_thin = thin_base;
#endif

}  // namespace

void gemm( char ta, char tb, int m, int n, int k, double alpha, const double* A, int lda, const double* B, int ldb,
           double beta, double* C, int ldc ) {
   if ( m <= 0 || n <= 0 ) {
      return;
   }
   if ( beta != 1.0 ) {
      for ( int j = 0; j < n; ++j ) {
         double* c = C + size_t( j ) * ldc;
         if ( beta == 0.0 ) {
            std::fill( c, c + m, 0.0 );
         } else {
            for ( int i = 0; i < m; ++i ) {
               c[ i ] *= beta;
            }
         }
      }
   }
   if ( k <= 0 || alpha == 0.0 ) {
      return;
   }
   std::vector<double>& pa = t_pack_a;
   std::vector<double>& pb = t_pack_b;
   if ( n <= 4 * kNR ) {
      // Few columns (matrix-vector products): A is read once per kNR columns, not packed.
      if ( pb.size() < size_t( k ) * kNR ) {
         pb.resize( size_t( k ) * kNR );
      }
      for ( int j0 = 0; j0 < n; j0 += kNR ) {
         const int nj = std::min( kNR, n - j0 );
         pack_b( tb, B, ldb, 0, k, j0, nj, pb.data() );
         g_thin( ta, m, nj, k, alpha, A, lda, pb.data(), C + size_t( j0 ) * ldc, ldc );
      }
      return;
   }
   for ( int jc = 0; jc < n; jc += kNC ) {
      const int nc = std::min( kNC, n - jc );
      for ( int pc = 0; pc < k; pc += kKC ) {
         const int kc = std::min( kKC, k - pc );
         const size_t nb = size_t( ( nc + kNR - 1 ) / kNR ) * kNR * kc;
         if ( pb.size() < nb ) {
            pb.resize( nb );
         }
         pack_b( tb, B, ldb, pc, kc, jc, nc, pb.data() );
         for ( int ic = 0; ic < m; ic += kMC ) {
            const int mc = std::min( kMC, m - ic );
            const size_t na = size_t( ( mc + kMR - 1 ) / kMR ) * kMR * kc;
            if ( pa.size() < na ) {
               pa.resize( na );
            }
            pack_a( ta, A, lda, ic, mc, pc, kc, pa.data() );
            g_macro( mc, nc, kc, alpha, pa.data(), pb.data(), C + size_t( jc ) * ldc + ic, ldc );
         }
      }
   }
}

// ---------------------------------------------------------------------------------------------
// Triangular solves, blocked: small solves on diagonal blocks, gemm for the rest.
// ---------------------------------------------------------------------------------------------

namespace {
constexpr int kTB = 64;
}

void trsm_llnu( int m, int n, const double* L, int ldl, double* B, int ldb ) {
   for ( int kb = 0; kb < m; kb += kTB ) {
      const int nb = std::min( kTB, m - kb );
      for ( int j = 0; j < n; ++j ) {
         double* b = B + size_t( j ) * ldb;
         for ( int kk = kb; kk < kb + nb; ++kk ) {
            const double x = b[ kk ];
            if ( x == 0.0 ) {
               continue;
            }
            const double* l = L + size_t( kk ) * ldl;
            for ( int i = kk + 1; i < kb + nb; ++i ) {
               b[ i ] -= l[ i ] * x;
            }
         }
      }
      if ( kb + nb < m ) {
         gemm( 'N', 'N', m - kb - nb, n, nb, -1.0, L + size_t( kb ) * ldl + kb + nb, ldl, B + kb, ldb, 1.0, B + kb + nb,
               ldb );
      }
   }
}

void trsm_lunn( int m, int n, const double* U, int ldu, double* B, int ldb ) {
   const int nblk = ( m + kTB - 1 ) / kTB;
   for ( int bk = nblk - 1; bk >= 0; --bk ) {
      const int kb = bk * kTB, nb = std::min( kTB, m - kb );
      for ( int j = 0; j < n; ++j ) {
         double* b = B + size_t( j ) * ldb;
         for ( int kk = kb + nb - 1; kk >= kb; --kk ) {
            const double* u = U + size_t( kk ) * ldu;
            b[ kk ] /= u[ kk ];
            const double x = b[ kk ];
            if ( x == 0.0 ) {
               continue;
            }
            for ( int i = kb; i < kk; ++i ) {
               b[ i ] -= u[ i ] * x;
            }
         }
      }
      if ( kb > 0 ) {
         gemm( 'N', 'N', kb, n, nb, -1.0, U + size_t( kb ) * ldu, ldu, B + kb, ldb, 1.0, B, ldb );
      }
   }
}

void trsm_runn( int m, int n, const double* U, int ldu, double* B, int ldb ) {
   for ( int jb = 0; jb < n; jb += kTB ) {
      const int nb = std::min( kTB, n - jb );
      if ( jb > 0 ) {
         gemm( 'N', 'N', m, nb, jb, -1.0, B, ldb, U + size_t( jb ) * ldu, ldu, 1.0, B + size_t( jb ) * ldb, ldb );
      }
      for ( int j = jb; j < jb + nb; ++j ) {
         double* bj = B + size_t( j ) * ldb;
         const double* u = U + size_t( j ) * ldu;
         for ( int l = jb; l < j; ++l ) {
            const double ulj = u[ l ];
            if ( ulj == 0.0 ) {
               continue;
            }
            const double* bl = B + size_t( l ) * ldb;
            for ( int i = 0; i < m; ++i ) {
               bj[ i ] -= bl[ i ] * ulj;
            }
         }
         const double inv = 1.0 / u[ j ];
         for ( int i = 0; i < m; ++i ) {
            bj[ i ] *= inv;
         }
      }
   }
}

// ---------------------------------------------------------------------------------------------
// Householder QR and one-sided Jacobi SVD.
// ---------------------------------------------------------------------------------------------

int householder_qr( int m, int k, double* A, std::vector<double>& tau, bool pivot, double stop,
                    std::vector<int>* perm ) {
   const int q = std::min( m, k );
   tau.assign( std::max( q, 0 ), 0.0 );
   std::vector<double> cn, cn0;
   if ( pivot ) {
      cn.resize( k );
      for ( int l = 0; l < k; ++l ) {
         const double* a = A + size_t( l ) * m;
         double s = 0;
         for ( int i = 0; i < m; ++i ) {
            s += a[ i ] * a[ i ];
         }
         cn[ l ] = s;
      }
      cn0 = cn;
      if ( perm ) {
         perm->resize( k );
         std::iota( perm->begin(), perm->end(), 0 );
      }
   }
   for ( int j = 0; j < q; ++j ) {
      if ( pivot ) {
         double rest = 0;
         int p = j;
         for ( int l = j; l < k; ++l ) {
            rest += cn[ l ];
            if ( cn[ l ] > cn[ p ] ) {
               p = l;
            }
         }
         if ( rest <= stop * stop ) {
            tau.resize( j );
            return j;
         }
         if ( p != j ) {
            std::swap_ranges( A + size_t( j ) * m, A + size_t( j + 1 ) * m, A + size_t( p ) * m );
            std::swap( cn[ j ], cn[ p ] );
            std::swap( cn0[ j ], cn0[ p ] );
            if ( perm ) {
               std::swap( ( *perm )[ j ], ( *perm )[ p ] );
            }
         }
      }
      double* x = A + size_t( j ) * m;
      const double alpha = x[ j ];
      double xn2 = 0;
      for ( int i = j + 1; i < m; ++i ) {
         xn2 += x[ i ] * x[ i ];
      }
      if ( xn2 > 0 ) {
         const double beta = -std::copysign( std::sqrt( alpha * alpha + xn2 ), alpha );
         tau[ j ] = ( beta - alpha ) / beta;
         const double scale = 1.0 / ( alpha - beta );
         for ( int i = j + 1; i < m; ++i ) {
            x[ i ] *= scale;
         }
         x[ j ] = beta;
         const double t = tau[ j ];
         for ( int l = j + 1; l < k; ++l ) {
            double* a = A + size_t( l ) * m;
            double w = a[ j ];
            for ( int i = j + 1; i < m; ++i ) {
               w += x[ i ] * a[ i ];
            }
            w *= t;
            a[ j ] -= w;
            for ( int i = j + 1; i < m; ++i ) {
               a[ i ] -= w * x[ i ];
            }
         }
      }
      if ( pivot ) {
         for ( int l = j + 1; l < k; ++l ) {
            const double* a = A + size_t( l ) * m;
            cn[ l ] -= a[ j ] * a[ j ];
            // Recompute a norm that has lost most of its digits to cancellation.
            if ( cn[ l ] <= 1e-6 * cn0[ l ] ) {
               double s = 0;
               for ( int i = j + 1; i < m; ++i ) {
                  s += a[ i ] * a[ i ];
               }
               cn[ l ] = s;
               cn0[ l ] = s;
            }
         }
      }
   }
   return q;
}

void apply_q( int m, int q, const double* A, const std::vector<double>& tau, int n, double* C, int ldc ) {
   for ( int j = q - 1; j >= 0; --j ) {
      const double t = tau[ j ];
      if ( t == 0.0 ) {
         continue;
      }
      const double* v = A + size_t( j ) * m;
      for ( int c = 0; c < n; ++c ) {
         double* x = C + size_t( c ) * ldc;
         double w = x[ j ];
         for ( int i = j + 1; i < m; ++i ) {
            w += v[ i ] * x[ i ];
         }
         w *= t;
         x[ j ] -= w;
         for ( int i = j + 1; i < m; ++i ) {
            x[ i ] -= w * v[ i ];
         }
      }
   }
}

void jacobi_svd( int m, int n, double* A, std::vector<double>& V, std::vector<double>& s ) {
   V.assign( size_t( n ) * n, 0.0 );
   for ( int j = 0; j < n; ++j ) {
      V[ size_t( j ) * n + j ] = 1.0;
   }
   s.assign( n, 0.0 );
   if ( n == 0 ) {
      return;
   }
   auto col = [&]( int j ) { return A + size_t( j ) * m; };
   auto norm2 = [&]( const double* a ) {
      double t = 0;
      for ( int i = 0; i < m; ++i ) {
         t += a[ i ] * a[ i ];
      }
      return t;
   };
   const double tol = std::max( 1e-15, std::sqrt( double( m ) ) * DBL_EPSILON );
   std::vector<double> d( n );
   for ( int sweep = 0; sweep < 60; ++sweep ) {
      for ( int j = 0; j < n; ++j ) {
         d[ j ] = norm2( col( j ) );
      }
      bool rotated = false;
      for ( int p = 0; p < n - 1; ++p ) {
         for ( int q = p + 1; q < n; ++q ) {
            const double alpha = d[ p ], beta = d[ q ];
            if ( alpha == 0.0 || beta == 0.0 ) {
               continue;
            }
            double* ap = col( p );
            double* aq = col( q );
            double gamma = 0;
            for ( int i = 0; i < m; ++i ) {
               gamma += ap[ i ] * aq[ i ];
            }
            if ( std::fabs( gamma ) <= tol * std::sqrt( alpha * beta ) ) {
               continue;
            }
            rotated = true;
            const double zeta = ( beta - alpha ) / ( 2.0 * gamma );
            const double t = std::copysign( 1.0, zeta ) / ( std::fabs( zeta ) + std::sqrt( 1.0 + zeta * zeta ) );
            const double c = 1.0 / std::sqrt( 1.0 + t * t ), sn = c * t;
            for ( int i = 0; i < m; ++i ) {
               const double x = ap[ i ], y = aq[ i ];
               ap[ i ] = c * x - sn * y;
               aq[ i ] = sn * x + c * y;
            }
            double* vp = V.data() + size_t( p ) * n;
            double* vq = V.data() + size_t( q ) * n;
            for ( int i = 0; i < n; ++i ) {
               const double x = vp[ i ], y = vq[ i ];
               vp[ i ] = c * x - sn * y;
               vq[ i ] = sn * x + c * y;
            }
            d[ p ] = alpha - t * gamma;
            d[ q ] = beta + t * gamma;
         }
      }
      if ( !rotated ) {
         break;
      }
   }
   for ( int j = 0; j < n; ++j ) {
      s[ j ] = std::sqrt( norm2( col( j ) ) );
   }
   // Sort by decreasing singular value.
   std::vector<int> order( n );
   std::iota( order.begin(), order.end(), 0 );
   std::stable_sort( order.begin(), order.end(), [&]( int a, int b ) { return s[ a ] > s[ b ]; } );
   bool sorted = true;
   for ( int j = 0; j < n; ++j ) {
      sorted = sorted && order[ j ] == j;
   }
   if ( sorted ) {
      return;
   }
   std::vector<double> A2( size_t( m ) * n ), V2( size_t( n ) * n ), s2( n );
   for ( int j = 0; j < n; ++j ) {
      const int o = order[ j ];
      std::copy( col( o ), col( o ) + m, A2.data() + size_t( j ) * m );
      std::copy( V.data() + size_t( o ) * n, V.data() + size_t( o + 1 ) * n, V2.data() + size_t( j ) * n );
      s2[ j ] = s[ o ];
   }
   std::copy( A2.begin(), A2.end(), A );
   V.swap( V2 );
   s.swap( s2 );
}

}  // namespace bh
