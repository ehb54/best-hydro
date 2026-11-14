// BEST H-matrix solver: small dense kernels (see dense.h).
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// gemm packs op(A) and op(B) into thread-local panels and runs an 8 x 4 register-blocked kernel;
// products with few columns (the matrix-vector products) do not pack A.  QR is blocked: unblocked
// panels of 32 columns, then their block reflectors (I - V T V^T) applied by gemm.  On x86-64 the
// inner kernels are compiled twice, for AVX2 with FMA and for the baseline, and the first is used
// when the processor has both.  Compile with -O3 so that the kernels are vectorized.
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

// op(A)[ i0 .. i0 + mc, l0 .. l0 + kc ) into strips of kMR rows: strip s, then l, then i (in double).
template <class TS>
void pack_a( char ta, const TS* A, int lda, int i0, int mc, int l0, int kc, double* p ) {
   for ( int s = 0; s < mc; s += kMR ) {
      const int mr = std::min( kMR, mc - s );
      if ( ta == 'N' ) {
         for ( int l = 0; l < kc; ++l, p += kMR ) {
            const TS* a = A + size_t( l0 + l ) * lda + i0 + s;
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

// op(B)[ l0 .. l0 + kc, j0 .. j0 + nc ) into strips of kNR columns: strip s, then l, then j (in double).
template <class TS>
void pack_b( char tb, const TS* B, int ldb, int l0, int kc, int j0, int nc, double* p ) {
   for ( int s = 0; s < nc; s += kNR ) {
      const int nr = std::min( kNR, nc - s );
      for ( int l = 0; l < kc; ++l, p += kNR ) {
         int j = 0;
         if ( tb == 'N' ) {
            for ( ; j < nr; ++j ) {
               p[ j ] = B[ size_t( j0 + s + j ) * ldb + l0 + l ];
            }
         } else {
            const TS* b = B + size_t( l0 + l ) * ldb + j0 + s;
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

// Four doubles as one vector (GCC and clang vector extensions): one AVX register in the AVX2
// clones, two SSE2 registers otherwise.  v4du is the same with 8-byte alignment, for loads and
// stores at any address.
// These helpers are always inlined, so the ABI of passing AVX vectors (-Wpsabi) never applies.
#if defined( __GNUC__ ) && !defined( __clang__ )
#pragma GCC diagnostic ignored "-Wpsabi"
#endif
typedef double v4d __attribute__( ( vector_size( 32 ) ) );
typedef double v4du __attribute__( ( vector_size( 32 ), aligned( 8 ) ) );

#define BH_INLINE inline __attribute__( ( always_inline ) )

BH_INLINE v4d load4( const double* p ) { return *reinterpret_cast<const v4du*>( p ); }
BH_INLINE void store4( double* p, v4d v ) { *reinterpret_cast<v4du*>( p ) = v; }
// The same for four floats, converted to and from double (vcvtps2pd and vcvtpd2ps with AVX).
typedef float v4fu __attribute__( ( vector_size( 16 ), aligned( 4 ) ) );
BH_INLINE v4d load4( const float* p ) {
   const v4fu f = *reinterpret_cast<const v4fu*>( p );
   return v4d{ f[ 0 ], f[ 1 ], f[ 2 ], f[ 3 ] };
}
BH_INLINE void store4( float* p, v4d v ) {
   *reinterpret_cast<v4fu*>( p ) = v4fu{ float( v[ 0 ] ), float( v[ 1 ] ), float( v[ 2 ] ), float( v[ 3 ] ) };
}
BH_INLINE v4d bcast( double x ) { return v4d{ x, x, x, x }; }
BH_INLINE double hsum( v4d v ) { return ( v[ 0 ] + v[ 1 ] ) + ( v[ 2 ] + v[ 3 ] ); }

// The 8 x 4 block acc += sum over l of a(l) b(l)^T: a(l) is 8 doubles at a + l * as, b(l) 4 doubles
// at b + 4 l.  c00 .. c31: column j, rows 0-3 and 4-7.
#define BH_KERNEL_8X4( kc, a, as, b )                                                                     \
   v4d c00 = { 0, 0, 0, 0 }, c01 = c00, c10 = c00, c11 = c00, c20 = c00, c21 = c00, c30 = c00, c31 = c00; \
   for ( int l = 0; l < ( kc ); ++l ) {                                                                   \
      const auto* al = ( a ) + size_t( l ) * ( as );                                                      \
      const double* bl = ( b ) + size_t( l ) * kNR;                                                       \
      const v4d a0 = load4( al ), a1 = load4( al + 4 );                                                   \
      v4d bj = bcast( bl[ 0 ] );                                                                          \
      c00 += a0 * bj;                                                                                     \
      c01 += a1 * bj;                                                                                     \
      bj = bcast( bl[ 1 ] );                                                                              \
      c10 += a0 * bj;                                                                                     \
      c11 += a1 * bj;                                                                                     \
      bj = bcast( bl[ 2 ] );                                                                              \
      c20 += a0 * bj;                                                                                     \
      c21 += a1 * bj;                                                                                     \
      bj = bcast( bl[ 3 ] );                                                                              \
      c30 += a0 * bj;                                                                                     \
      c31 += a1 * bj;                                                                                     \
   }                                                                                                      \
   const v4d acc[ kNR ][ 2 ] = { { c00, c01 }, { c10, c11 }, { c20, c21 }, { c30, c31 } };

// C[ 0 .. mr, 0 .. nr ) += alpha acc.
template <class TC>
BH_INLINE void add_tile( const v4d ( &acc )[ kNR ][ 2 ], int mr, int nr, double alpha, TC* C, int ldc ) {
   const v4d va = bcast( alpha );
   for ( int j = 0; j < nr; ++j ) {
      TC* c = C + size_t( j ) * ldc;
      if ( mr == kMR ) {
         store4( c, load4( c ) + va * acc[ j ][ 0 ] );
         store4( c + 4, load4( c + 4 ) + va * acc[ j ][ 1 ] );
      } else {
         for ( int i = 0; i < mr; ++i ) {
            c[ i ] += alpha * acc[ j ][ i / 4 ][ i % 4 ];
         }
      }
   }
}

// C[ mc x nc ] += alpha (packed A) (packed B).
template <class TC>
BH_INLINE void macro_body( int mc, int nc, int kc, double alpha, const double* pa, const double* pb, TC* C, int ldc ) {
   for ( int jb = 0; jb < nc; jb += kNR ) {
      const int nr = std::min( kNR, nc - jb );
      const double* bp = pb + size_t( jb / kNR ) * kc * kNR;
      for ( int ib = 0; ib < mc; ib += kMR ) {
         const int mr = std::min( kMR, mc - ib );
         const double* ap = pa + size_t( ib / kMR ) * kc * kMR;
         BH_KERNEL_8X4( kc, ap, kMR, bp )
         add_tile( acc, mr, nr, alpha, C + size_t( jb ) * ldc + ib, ldc );
      }
   }
}

// C[ m x nj ] += alpha op(A) op(B)[ :, 0 .. nj ) for nj <= kNR columns, A not packed.  bp: op(B)
// packed as one kNR-column strip (pack_b).  'N': kMR-row strips of A against all nj columns at
// once; 'T': for each column of A, its dot products with the nj columns at once.
template <class TA, class TC>
BH_INLINE void thin_body( char ta, int m, int nj, int k, double alpha, const TA* A, int lda, const double* bp, TC* C,
                          int ldc ) {
   if ( ta == 'N' ) {
      int ib = 0;
      for ( ; ib + kMR <= m; ib += kMR ) {
         BH_KERNEL_8X4( k, A + ib, lda, bp )
         add_tile( acc, kMR, nj, alpha, C + ib, ldc );
      }
      for ( ; ib < m; ++ib ) {
         v4d s = { 0, 0, 0, 0 };
         for ( int l = 0; l < k; ++l ) {
            s += bcast( A[ size_t( l ) * lda + ib ] ) * load4( bp + size_t( l ) * kNR );
         }
         for ( int j = 0; j < nj; ++j ) {
            C[ size_t( j ) * ldc + ib ] += alpha * s[ j ];
         }
      }
   } else {
      for ( int i = 0; i < m; ++i ) {
         const TA* a = A + size_t( i ) * lda;
         v4d s0 = { 0, 0, 0, 0 }, s1 = s0, s2 = s0, s3 = s0;
         int l = 0;
         for ( ; l + 4 <= k; l += 4 ) {
            const double* b = bp + size_t( l ) * kNR;
            s0 += bcast( a[ l ] ) * load4( b );
            s1 += bcast( a[ l + 1 ] ) * load4( b + 4 );
            s2 += bcast( a[ l + 2 ] ) * load4( b + 8 );
            s3 += bcast( a[ l + 3 ] ) * load4( b + 12 );
         }
         for ( ; l < k; ++l ) {
            s0 += bcast( a[ l ] ) * load4( bp + size_t( l ) * kNR );
         }
         const v4d s = ( s0 + s1 ) + ( s2 + s3 );
         for ( int j = 0; j < nj; ++j ) {
            C[ size_t( j ) * ldc + i ] += alpha * s[ j ];
         }
      }
   }
}

template <class TC>
void macro_base( int mc, int nc, int kc, double alpha, const double* pa, const double* pb, TC* C, int ldc ) {
   macro_body( mc, nc, kc, alpha, pa, pb, C, ldc );
}
template <class TA, class TC>
void thin_base( char ta, int m, int nj, int k, double alpha, const TA* A, int lda, const double* bp, TC* C, int ldc ) {
   thin_body( ta, m, nj, k, alpha, A, lda, bp, C, ldc );
}

#if defined( __x86_64__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
template <class TC>
__attribute__( ( target( "avx2,fma" ) ) ) void macro_avx2( int mc, int nc, int kc, double alpha, const double* pa,
                                                           const double* pb, TC* C, int ldc ) {
   macro_body( mc, nc, kc, alpha, pa, pb, C, ldc );
}
template <class TA, class TC>
__attribute__( ( target( "avx2,fma" ) ) ) void thin_avx2( char ta, int m, int nj, int k, double alpha, const TA* A,
                                                          int lda, const double* bp, TC* C, int ldc ) {
   thin_body( ta, m, nj, k, alpha, A, lda, bp, C, ldc );
}
bool has_avx2() {
   __builtin_cpu_init();
   return __builtin_cpu_supports( "avx2" ) && __builtin_cpu_supports( "fma" );
}
const bool g_avx2 = has_avx2();
template <class TC>
void macro( int mc, int nc, int kc, double alpha, const double* pa, const double* pb, TC* C, int ldc ) {
   ( g_avx2 ? macro_avx2<TC> : macro_base<TC> )( mc, nc, kc, alpha, pa, pb, C, ldc );
}
template <class TA, class TC>
void thin( char ta, int m, int nj, int k, double alpha, const TA* A, int lda, const double* bp, TC* C, int ldc ) {
   ( g_avx2 ? thin_avx2<TA, TC> : thin_base<TA, TC> )( ta, m, nj, k, alpha, A, lda, bp, C, ldc );
}
#else
template <class TC>
void macro( int mc, int nc, int kc, double alpha, const double* pa, const double* pb, TC* C, int ldc ) {
   macro_base( mc, nc, kc, alpha, pa, pb, C, ldc );
}
template <class TA, class TC>
void thin( char ta, int m, int nj, int k, double alpha, const TA* A, int lda, const double* bp, TC* C, int ldc ) {
   thin_base( ta, m, nj, k, alpha, A, lda, bp, C, ldc );
}
#endif

}  // namespace

template <class TA, class TB, class TC>
void gemm( char ta, char tb, int m, int n, int k, double alpha, const TA* A, int lda, const TB* B, int ldb, double beta,
           TC* C, int ldc ) {
   if ( m <= 0 || n <= 0 ) {
      return;
   }
   if ( beta != 1.0 ) {
      for ( int j = 0; j < n; ++j ) {
         TC* c = C + size_t( j ) * ldc;
         if ( beta == 0.0 ) {
            std::fill( c, c + m, TC( 0 ) );
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
         thin( ta, m, nj, k, alpha, A, lda, pb.data(), C + size_t( j0 ) * ldc, ldc );
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
            macro( mc, nc, kc, alpha, pa.data(), pb.data(), C + size_t( jc ) * ldc + ic, ldc );
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

template <class TT>
void trsm_llnu( int m, int n, const TT* L, int ldl, double* B, int ldb ) {
   for ( int kb = 0; kb < m; kb += kTB ) {
      const int nb = std::min( kTB, m - kb );
      for ( int j = 0; j < n; ++j ) {
         double* b = B + size_t( j ) * ldb;
         for ( int kk = kb; kk < kb + nb; ++kk ) {
            const double x = b[ kk ];
            if ( x == 0.0 ) {
               continue;
            }
            const TT* l = L + size_t( kk ) * ldl;
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

template <class TT>
void trsm_lunn( int m, int n, const TT* U, int ldu, double* B, int ldb ) {
   const int nblk = ( m + kTB - 1 ) / kTB;
   for ( int bk = nblk - 1; bk >= 0; --bk ) {
      const int kb = bk * kTB, nb = std::min( kTB, m - kb );
      for ( int j = 0; j < n; ++j ) {
         double* b = B + size_t( j ) * ldb;
         for ( int kk = kb + nb - 1; kk >= kb; --kk ) {
            const TT* u = U + size_t( kk ) * ldu;
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

template <class TT>
void trsm_runn( int m, int n, const TT* U, int ldu, double* B, int ldb ) {
   for ( int jb = 0; jb < n; jb += kTB ) {
      const int nb = std::min( kTB, n - jb );
      if ( jb > 0 ) {
         gemm( 'N', 'N', m, nb, jb, -1.0, B, ldb, U + size_t( jb ) * ldu, ldu, 1.0, B + size_t( jb ) * ldb, ldb );
      }
      for ( int j = jb; j < jb + nb; ++j ) {
         double* bj = B + size_t( j ) * ldb;
         const TT* u = U + size_t( j ) * ldu;
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

// Every combination of double and float operands, and the triangles in either.
#define BH_GEMM( TA, TB, TC ) \
   template void gemm<TA, TB, TC>( char, char, int, int, int, double, const TA*, int, const TB*, int, double, TC*, \
                                   int );
BH_GEMM( double, double, double )
BH_GEMM( double, double, float )
BH_GEMM( double, float, double )
BH_GEMM( double, float, float )
BH_GEMM( float, double, double )
BH_GEMM( float, double, float )
BH_GEMM( float, float, double )
BH_GEMM( float, float, float )
#undef BH_GEMM
template void trsm_llnu<double>( int, int, const double*, int, double*, int );
template void trsm_llnu<float>( int, int, const float*, int, double*, int );
template void trsm_lunn<double>( int, int, const double*, int, double*, int );
template void trsm_lunn<float>( int, int, const float*, int, double*, int );
template void trsm_runn<double>( int, int, const double*, int, double*, int );
template void trsm_runn<float>( int, int, const float*, int, double*, int );

// ---------------------------------------------------------------------------------------------
// Householder QR and one-sided Jacobi SVD.
// ---------------------------------------------------------------------------------------------

namespace {

constexpr int kQB = 32;  // QR panel width

// x . y over n entries, four partial sums.
BH_INLINE double dot_body( int n, const double* x, const double* y ) {
   v4d a0 = { 0, 0, 0, 0 }, a1 = a0;
   int i = 0;
   for ( ; i + 8 <= n; i += 8 ) {
      a0 += load4( x + i ) * load4( y + i );
      a1 += load4( x + i + 4 ) * load4( y + i + 4 );
   }
   double s = hsum( a0 + a1 );
   for ( ; i < n; ++i ) {
      s += x[ i ] * y[ i ];
   }
   return s;
}

// The Householder reflector of x[ 0 .. n ): x[ 0 ] becomes beta, x[ 1 .. n ) the vector v (v[0] = 1
// implied); returns tau (0: no reflection needed).
BH_INLINE double reflector_body( int n, double* x ) {
   const double alpha = x[ 0 ];
   const double xn2 = dot_body( n - 1, x + 1, x + 1 );
   if ( !( xn2 > 0 ) ) {
      return 0.0;
   }
   const double beta = -std::copysign( std::sqrt( alpha * alpha + xn2 ), alpha );
   const double scale = 1.0 / ( alpha - beta );
   for ( int i = 1; i < n; ++i ) {
      x[ i ] *= scale;
   }
   x[ 0 ] = beta;
   return ( beta - alpha ) / beta;
}

// a <- (I - tau v v^T) a, both n long, v[ 0 ] = 1 implied.
BH_INLINE void reflect_body( int n, double tau, const double* v, double* a ) {
   const double w = tau * ( a[ 0 ] + dot_body( n - 1, v + 1, a + 1 ) );
   a[ 0 ] -= w;
   const v4d vw = bcast( w );
   int i = 1;
   for ( ; i + 4 <= n; i += 4 ) {
      store4( a + i, load4( a + i ) - vw * load4( v + i ) );
   }
   for ( ; i < n; ++i ) {
      a[ i ] -= w * v[ i ];
   }
}

// Unblocked QR of the m x k matrix A (leading dimension lda): reflectors applied within A only.
BH_INLINE void panel_body( int m, int k, double* A, int lda, double* tau ) {
   const int q = std::min( m, k );
   for ( int j = 0; j < q; ++j ) {
      double* x = A + size_t( j ) * lda + j;
      tau[ j ] = reflector_body( m - j, x );
      if ( tau[ j ] != 0.0 ) {
         for ( int l = j + 1; l < k; ++l ) {
            reflect_body( m - j, tau[ j ], x, A + size_t( l ) * lda + j );
         }
      }
   }
}

// Pivoted unblocked QR (see householder_qr); returns the number of steps.
BH_INLINE int pivot_body( int m, int k, double* A, double* tau, double stop, int* perm,
                                                           double* cn, double* cn0 ) {
   const int q = std::min( m, k );
   for ( int l = 0; l < k; ++l ) {
      const double* a = A + size_t( l ) * m;
      cn[ l ] = cn0[ l ] = dot_body( m, a, a );
      perm[ l ] = l;
   }
   for ( int j = 0; j < q; ++j ) {
      double rest = 0;
      int p = j;
      for ( int l = j; l < k; ++l ) {
         rest += cn[ l ];
         if ( cn[ l ] > cn[ p ] ) {
            p = l;
         }
      }
      if ( rest <= stop * stop ) {
         return j;
      }
      if ( p != j ) {
         std::swap_ranges( A + size_t( j ) * m, A + size_t( j + 1 ) * m, A + size_t( p ) * m );
         std::swap( cn[ j ], cn[ p ] );
         std::swap( cn0[ j ], cn0[ p ] );
         std::swap( perm[ j ], perm[ p ] );
      }
      double* x = A + size_t( j ) * m + j;
      tau[ j ] = reflector_body( m - j, x );
      for ( int l = j + 1; l < k; ++l ) {
         double* a = A + size_t( l ) * m;
         if ( tau[ j ] != 0.0 ) {
            reflect_body( m - j, tau[ j ], x, a + j );
         }
         cn[ l ] -= a[ j ] * a[ j ];
         // Recompute a norm that has lost most of its digits to cancellation.
         if ( cn[ l ] <= 1e-6 * cn0[ l ] ) {
            cn[ l ] = cn0[ l ] = dot_body( m - j - 1, a + j + 1, a + j + 1 );
         }
      }
   }
   return q;
}

// C (m x n) <- (H_0 ... H_{q-1}) C, unblocked, reflectors in A (leading dimension lda).
BH_INLINE void applyq_body( int m, int q, const double* A, int lda, const double* tau,
                                                             int n, double* C, int ldc ) {
   for ( int j = q - 1; j >= 0; --j ) {
      if ( tau[ j ] == 0.0 ) {
         continue;
      }
      const double* v = A + size_t( j ) * lda + j;
      for ( int c = 0; c < n; ++c ) {
         reflect_body( m - j, tau[ j ], v, C + size_t( c ) * ldc + j );
      }
   }
}

using PanelFn = void ( * )( int, int, double*, int, double* );
using PivotFn = int ( * )( int, int, double*, double*, double, int*, double*, double* );
using ApplyqFn = void ( * )( int, int, const double*, int, const double*, int, double*, int );

void panel_base( int m, int k, double* A, int lda, double* tau ) { panel_body( m, k, A, lda, tau ); }
int pivot_base( int m, int k, double* A, double* tau, double stop, int* perm, double* cn, double* cn0 ) {
   return pivot_body( m, k, A, tau, stop, perm, cn, cn0 );
}
void applyq_base( int m, int q, const double* A, int lda, const double* tau, int n, double* C, int ldc ) {
   applyq_body( m, q, A, lda, tau, n, C, ldc );
}
#if defined( __x86_64__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
__attribute__( ( target( "avx2,fma" ) ) ) void panel_avx2( int m, int k, double* A, int lda, double* tau ) {
   panel_body( m, k, A, lda, tau );
}
__attribute__( ( target( "avx2,fma" ) ) ) int pivot_avx2( int m, int k, double* A, double* tau, double stop, int* perm,
                                                          double* cn, double* cn0 ) {
   return pivot_body( m, k, A, tau, stop, perm, cn, cn0 );
}
__attribute__( ( target( "avx2,fma" ) ) ) void applyq_avx2( int m, int q, const double* A, int lda, const double* tau,
                                                            int n, double* C, int ldc ) {
   applyq_body( m, q, A, lda, tau, n, C, ldc );
}
const PanelFn g_panel = g_avx2 ? panel_avx2 : panel_base;
const PivotFn g_pivot = g_avx2 ? pivot_avx2 : pivot_base;
const ApplyqFn g_applyq = g_avx2 ? applyq_avx2 : applyq_base;
#else
const PanelFn g_panel = panel_base;
const PivotFn g_pivot = pivot_base;
const ApplyqFn g_applyq = applyq_base;
#endif

thread_local std::vector<double> t_qr_v, t_qr_t, t_qr_w;

// The block reflector of panel reflectors j0 .. j0 + jb of A (m rows, leading dimension m):
// H_j0 ... H_{j0+jb-1} = I - V T V^T.  V (mm x jb, mm = m - j0) explicit, T upper triangular.
void block_reflector( int m, const double* A, const double* tau, int j0, int jb, std::vector<double>& V,
                      std::vector<double>& T ) {
   const int mm = m - j0;
   V.assign( size_t( mm ) * jb, 0.0 );
   T.assign( size_t( jb ) * jb, 0.0 );
   for ( int j = 0; j < jb; ++j ) {
      double* v = V.data() + size_t( j ) * mm;
      const double* a = A + size_t( j0 + j ) * m + j0;
      v[ j ] = 1.0;
      for ( int i = j + 1; i < mm; ++i ) {
         v[ i ] = a[ i ];
      }
   }
   // T(i, i) = tau_i; T(0:i, i) = -tau_i T(0:i, 0:i) V(:, 0:i)^T v_i.
   std::vector<double> w( jb );
   for ( int i = 0; i < jb; ++i ) {
      const double ti = tau[ j0 + i ];
      T[ size_t( i ) * jb + i ] = ti;
      if ( i == 0 || ti == 0.0 ) {
         continue;
      }
      const double* vi = V.data() + size_t( i ) * mm;
      for ( int l = 0; l < i; ++l ) {
         const double* vl = V.data() + size_t( l ) * mm;
         double s = 0;
         for ( int r = i; r < mm; ++r ) {  // vl and vi overlap from row i (vi is 0 above, 1 at i)
            s += vl[ r ] * vi[ r ];
         }
         w[ l ] = -ti * s;
      }
      for ( int r = 0; r < i; ++r ) {
         double s = 0;
         for ( int l = r; l < i; ++l ) {
            s += T[ size_t( l ) * jb + r ] * w[ l ];
         }
         T[ size_t( i ) * jb + r ] = s;
      }
   }
}

// C (mm x n) <- (I - V T V^T) C (trans = false) or (I - V T^T V^T) C (trans = true).
void apply_block( int mm, int jb, const std::vector<double>& V, const std::vector<double>& T, bool trans, int n,
                  double* C, int ldc, std::vector<double>& W ) {
   W.resize( size_t( jb ) * n );
   gemm( 'T', 'N', jb, n, mm, 1.0, V.data(), mm, C, ldc, 0.0, W.data(), jb );  // W = V^T C
   // W <- T W or T^T W (T upper triangular, jb x jb).
   for ( int c = 0; c < n; ++c ) {
      double* w = W.data() + size_t( c ) * jb;
      if ( !trans ) {
         for ( int r = 0; r < jb; ++r ) {
            double s = 0;
            for ( int l = r; l < jb; ++l ) {
               s += T[ size_t( l ) * jb + r ] * w[ l ];
            }
            w[ r ] = s;
         }
      } else {
         for ( int r = jb - 1; r >= 0; --r ) {
            double s = 0;
            for ( int l = 0; l <= r; ++l ) {
               s += T[ size_t( r ) * jb + l ] * w[ l ];
            }
            w[ r ] = s;
         }
      }
   }
   gemm( 'N', 'N', mm, n, jb, -1.0, V.data(), mm, W.data(), jb, 1.0, C, ldc );  // C -= V W
}

}  // namespace

int householder_qr( int m, int k, double* A, std::vector<double>& tau, bool pivot, double stop,
                    std::vector<int>* perm ) {
   const int q = std::min( m, k );
   tau.assign( std::max( q, 0 ), 0.0 );
   if ( q <= 0 ) {
      return 0;
   }
   if ( pivot ) {
      std::vector<int> p( k );
      std::vector<double> cn( k ), cn0( k );
      const int steps = g_pivot( m, k, A, tau.data(), stop, p.data(), cn.data(), cn0.data() );
      tau.resize( steps );
      if ( perm ) {
         perm->swap( p );
      }
      return steps;
   }
   if ( k < 2 * kQB || m < 2 * kQB ) {
      g_panel( m, k, A, m, tau.data() );
      return q;
   }
   // Blocked: factor a panel, then apply its block reflector to the columns on its right.
   for ( int j0 = 0; j0 < q; j0 += kQB ) {
      const int jb = std::min( kQB, q - j0 );
      g_panel( m - j0, jb, A + size_t( j0 ) * m + j0, m, tau.data() + j0 );
      if ( j0 + jb < k ) {
         block_reflector( m, A, tau.data(), j0, jb, t_qr_v, t_qr_t );
         apply_block( m - j0, jb, t_qr_v, t_qr_t, true, k - j0 - jb, A + size_t( j0 + jb ) * m + j0, m, t_qr_w );
      }
   }
   return q;
}

void apply_q( int m, int q, const double* A, const std::vector<double>& tau, int n, double* C, int ldc ) {
   if ( q <= 0 || n <= 0 ) {
      return;
   }
   if ( q < 2 * kQB || m < 2 * kQB ) {
      g_applyq( m, q, A, m, tau.data(), n, C, ldc );
      return;
   }
   // Q C = H_0 ... H_{q-1} C: block reflectors from the last panel to the first.
   const int nb = ( q + kQB - 1 ) / kQB;
   for ( int b = nb - 1; b >= 0; --b ) {
      const int j0 = b * kQB, jb = std::min( kQB, q - j0 );
      block_reflector( m, A, tau.data(), j0, jb, t_qr_v, t_qr_t );
      apply_block( m - j0, jb, t_qr_v, t_qr_t, false, n, C + j0, ldc, t_qr_w );
   }
}

namespace {

// The rotations of jacobi_svd: sweeps until no pair of columns of A (m x n) is far from
// orthogonal, the same rotations applied to V (n x n).
BH_INLINE void sweeps_body( int m, int n, double* A, double* V, double* d ) {
   const double tol = std::max( 1e-15, std::sqrt( double( m ) ) * DBL_EPSILON );
   for ( int sweep = 0; sweep < 60; ++sweep ) {
      for ( int j = 0; j < n; ++j ) {
         const double* a = A + size_t( j ) * m;
         d[ j ] = dot_body( m, a, a );
      }
      bool rotated = false;
      for ( int p = 0; p < n - 1; ++p ) {
         for ( int q = p + 1; q < n; ++q ) {
            const double alpha = d[ p ], beta = d[ q ];
            if ( alpha == 0.0 || beta == 0.0 ) {
               continue;
            }
            double* ap = A + size_t( p ) * m;
            double* aq = A + size_t( q ) * m;
            const double gamma = dot_body( m, ap, aq );
            if ( std::fabs( gamma ) <= tol * std::sqrt( alpha * beta ) ) {
               continue;
            }
            rotated = true;
            const double zeta = ( beta - alpha ) / ( 2.0 * gamma );
            const double t = std::copysign( 1.0, zeta ) / ( std::fabs( zeta ) + std::sqrt( 1.0 + zeta * zeta ) );
            const double c = 1.0 / std::sqrt( 1.0 + t * t ), sn = c * t;
            const v4d vc = bcast( c ), vs = bcast( sn );
            int i = 0;
            for ( ; i + 4 <= m; i += 4 ) {
               const v4d x = load4( ap + i ), y = load4( aq + i );
               store4( ap + i, vc * x - vs * y );
               store4( aq + i, vs * x + vc * y );
            }
            for ( ; i < m; ++i ) {
               const double x = ap[ i ], y = aq[ i ];
               ap[ i ] = c * x - sn * y;
               aq[ i ] = sn * x + c * y;
            }
            double* vp = V + size_t( p ) * n;
            double* vq = V + size_t( q ) * n;
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
}

using SweepsFn = void ( * )( int, int, double*, double*, double* );
void sweeps_base( int m, int n, double* A, double* V, double* d ) { sweeps_body( m, n, A, V, d ); }
#if defined( __x86_64__ ) && ( defined( __GNUC__ ) || defined( __clang__ ) )
__attribute__( ( target( "avx2,fma" ) ) ) void sweeps_avx2( int m, int n, double* A, double* V, double* d ) {
   sweeps_body( m, n, A, V, d );
}
const SweepsFn g_sweeps = g_avx2 ? sweeps_avx2 : sweeps_base;
#else
const SweepsFn g_sweeps = sweeps_base;
#endif

}  // namespace

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
   std::vector<double> d( n );
   g_sweeps( m, n, A, V.data(), d.data() );
   for ( int j = 0; j < n; ++j ) {
      const double* a = col( j );
      s[ j ] = std::sqrt( dot_body( m, a, a ) );
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
