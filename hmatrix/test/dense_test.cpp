// BEST H-matrix solver: checks of the dense kernels (dense.h) and of the low-rank truncation
// (mat.h) against plain loops, on random matrices of many shapes.  Prints one line per check and
// exits nonzero if any fails.  Run by make test.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "dense.h"
#include "mat.h"

using namespace bh;

namespace {

std::mt19937_64 rng( 2026 );

std::vector<double> random( size_t n ) {
   std::uniform_real_distribution<double> u( -1.0, 1.0 );
   std::vector<double> v( n );
   for ( double& x : v ) {
      x = u( rng );
   }
   return v;
}

int fails = 0;

void check( bool ok, const char* what, double err ) {
   std::printf( "%s - dense kernels: %s (%.1e)\n", ok ? "ok  " : "FAIL", what, err );
   fails += ok ? 0 : 1;
}

double ref( char t, const std::vector<double>& A, int ld, int i, int l ) {
   return t == 'N' ? A[ size_t( l ) * ld + i ] : A[ size_t( i ) * ld + l ];
}

// Largest |gemm - loops| / (k max|A| max|B|) over shapes that take both code paths.
double gemm_error() {
   const int shapes[][ 3 ] = { { 1, 1, 1 },   { 7, 3, 5 },     { 96, 6, 96 },   { 96, 11, 300 }, { 20, 4, 3000 },
                               { 9, 17, 33 }, { 130, 70, 260 }, { 300, 5, 17 }, { 257, 1030, 9 }, { 1, 40, 1 } };
   double worst = 0;
   for ( char ta : { 'N', 'T' } ) {
      for ( char tb : { 'N', 'T' } ) {
         for ( const auto& s : shapes ) {
            const int m = s[ 0 ], n = s[ 1 ], k = s[ 2 ];
            const int lda = ( ta == 'N' ? m : k ) + 3, ldb = ( tb == 'N' ? k : n ) + 2, ldc = m + 1;
            const std::vector<double> A = random( size_t( lda ) * ( ta == 'N' ? k : m ) );
            const std::vector<double> B = random( size_t( ldb ) * ( tb == 'N' ? n : k ) );
            std::vector<double> C = random( size_t( ldc ) * n ), C0 = C;
            const double alpha = -0.7, beta = 0.3;
            gemm( ta, tb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc );
            for ( int j = 0; j < n; ++j ) {
               for ( int i = 0; i < m; ++i ) {
                  double t = 0;
                  for ( int l = 0; l < k; ++l ) {
                     t += ref( ta, A, lda, i, l ) * ref( tb, B, ldb, l, j );
                  }
                  const double want = alpha * t + beta * C0[ size_t( j ) * ldc + i ];
                  worst = std::max( worst, std::fabs( C[ size_t( j ) * ldc + i ] - want ) / k );
               }
               // The padding row must be untouched.
               if ( C[ size_t( j ) * ldc + m ] != C0[ size_t( j ) * ldc + m ] ) {
                  worst = 1;
               }
            }
         }
      }
   }
   return worst;
}

// Largest residual |T X - B| / |B| of the three triangular solves.
double trsm_error() {
   double worst = 0;
   for ( int m : { 1, 5, 64, 65, 200 } ) {
      for ( int n : { 1, 3, 70 } ) {
         // A well-conditioned triangle: unit diagonal (L) or a diagonal near 2 (U), small entries.
         std::vector<double> T = random( size_t( m ) * m );
         for ( double& x : T ) {
            x *= 0.3 / std::sqrt( double( m ) );
         }
         for ( int i = 0; i < m; ++i ) {
            T[ size_t( i ) * m + i ] = 2.0 + 0.5 * T[ size_t( i ) * m + i ];
         }
         for ( int form = 0; form < 3; ++form ) {
            const int rows = form == 2 ? n : m, cols = form == 2 ? m : n;
            std::vector<double> B = random( size_t( rows ) * cols ), X = B;
            if ( form == 0 ) {
               trsm_llnu( m, n, T.data(), m, X.data(), rows );
            } else if ( form == 1 ) {
               trsm_lunn( m, n, T.data(), m, X.data(), rows );
            } else {
               trsm_runn( n, m, T.data(), m, X.data(), rows );
            }
            double e = 0, b = 0;
            for ( int j = 0; j < cols; ++j ) {
               for ( int i = 0; i < rows; ++i ) {
                  double t = 0;
                  if ( form == 0 ) {  // L X, unit diagonal
                     t = X[ size_t( j ) * rows + i ];
                     for ( int l = 0; l < i; ++l ) {
                        t += T[ size_t( l ) * m + i ] * X[ size_t( j ) * rows + l ];
                     }
                  } else if ( form == 1 ) {  // U X
                     for ( int l = i; l < m; ++l ) {
                        t += T[ size_t( l ) * m + i ] * X[ size_t( j ) * rows + l ];
                     }
                  } else {  // X U
                     for ( int l = 0; l <= j; ++l ) {
                        t += X[ size_t( l ) * rows + i ] * T[ size_t( j ) * m + l ];
                     }
                  }
                  e = std::max( e, std::fabs( t - B[ size_t( j ) * rows + i ] ) );
                  b = std::max( b, std::fabs( B[ size_t( j ) * rows + i ] ) );
               }
            }
            worst = std::max( worst, e / b );
         }
      }
   }
   return worst;
}

// Q R = A and Q^T Q = I, with and without pivoting, for tall, square and wide A.
double qr_error() {
   double worst = 0;
   for ( const auto& s : std::vector<std::pair<int, int>>{ { 50, 7 }, { 30, 30 }, { 8, 20 }, { 1, 1 }, { 300, 60 } } ) {
      for ( bool pivot : { false, true } ) {
         const int m = s.first, k = s.second;
         const std::vector<double> A = random( size_t( m ) * k );
         std::vector<double> F = A, tau;
         std::vector<int> perm;
         const int q = householder_qr( m, k, F.data(), tau, pivot, 0.0, &perm );
         // Q R column by column: R's column j is the upper part of F's, padded with zeros.
         std::vector<double> QR( size_t( m ) * k, 0.0 ), Q( size_t( m ) * q, 0.0 );
         for ( int j = 0; j < k; ++j ) {
            for ( int i = 0; i <= std::min( j, q - 1 ); ++i ) {
               QR[ size_t( j ) * m + i ] = F[ size_t( j ) * m + i ];
            }
         }
         apply_q( m, q, F.data(), tau, k, QR.data(), m );
         for ( int j = 0; j < q; ++j ) {
            Q[ size_t( j ) * m + j ] = 1.0;
         }
         apply_q( m, q, F.data(), tau, q, Q.data(), m );
         for ( int j = 0; j < k; ++j ) {
            const int c = pivot ? perm[ j ] : j;
            for ( int i = 0; i < m; ++i ) {
               worst = std::max( worst, std::fabs( QR[ size_t( j ) * m + i ] - A[ size_t( c ) * m + i ] ) );
            }
         }
         for ( int a = 0; a < q; ++a ) {
            for ( int b = 0; b < q; ++b ) {
               double t = 0;
               for ( int i = 0; i < m; ++i ) {
                  t += Q[ size_t( a ) * m + i ] * Q[ size_t( b ) * m + i ];
               }
               worst = std::max( worst, std::fabs( t - ( a == b ? 1.0 : 0.0 ) ) );
            }
         }
      }
   }
   return worst;
}

// A V = (U S) with V orthogonal and U S's columns orthogonal, s sorted; and the singular values
// of a matrix built from known ones.
double svd_error() {
   double worst = 0;
   for ( const auto& s : std::vector<std::pair<int, int>>{ { 40, 40 }, { 100, 12 }, { 3, 3 }, { 1, 1 }, { 96, 96 } } ) {
      const int m = s.first, n = s.second;
      const std::vector<double> A = random( size_t( m ) * n );
      std::vector<double> X = A, V, sv;
      jacobi_svd( m, n, X.data(), V, sv );
      for ( int j = 0; j < n; ++j ) {
         if ( j > 0 && sv[ j ] > sv[ j - 1 ] ) {
            worst = 1;
         }
         for ( int i = 0; i < m; ++i ) {
            double t = 0;
            for ( int l = 0; l < n; ++l ) {
               t += A[ size_t( l ) * m + i ] * V[ size_t( j ) * n + l ];
            }
            worst = std::max( worst, std::fabs( t - X[ size_t( j ) * m + i ] ) );
         }
         for ( int b = 0; b < n; ++b ) {
            double t = 0, u = 0;
            for ( int i = 0; i < n; ++i ) {
               t += V[ size_t( j ) * n + i ] * V[ size_t( b ) * n + i ];
            }
            for ( int i = 0; i < m; ++i ) {
               u += X[ size_t( j ) * m + i ] * X[ size_t( b ) * m + i ];
            }
            worst = std::max( worst, std::fabs( t - ( j == b ? 1.0 : 0.0 ) ) );
            if ( j != b ) {
               worst = std::max( worst, std::fabs( u ) / ( sv[ 0 ] * sv[ 0 ] ) );
            }
         }
      }
   }
   // Known singular values 10^-0 .. 10^-11: Q1 diag(s) Q2^T from two QR factors.
   const int m = 30, n = 12;
   std::vector<double> G1 = random( size_t( m ) * n ), G2 = random( size_t( n ) * n ), t1, t2;
   householder_qr( m, n, G1.data(), t1 );
   householder_qr( n, n, G2.data(), t2 );
   // D = Q1 [ diag(s) Q2^T ; 0 ]: build Q2 diag(s) (n x n), transpose it into D's top rows.
   std::vector<double> D( size_t( m ) * n, 0.0 );
   std::vector<double> E( size_t( n ) * n, 0.0 );
   for ( int j = 0; j < n; ++j ) {
      E[ size_t( j ) * n + j ] = std::pow( 10.0, -j );
   }
   apply_q( n, n, G2.data(), t2, n, E.data(), n );
   for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i < n; ++i ) {
         D[ size_t( j ) * m + i ] = E[ size_t( i ) * n + j ];
      }
   }
   apply_q( m, n, G1.data(), t1, n, D.data(), m );
   std::vector<double> V, sv;
   jacobi_svd( m, n, D.data(), V, sv );
   for ( int j = 0; j < n; ++j ) {
      worst = std::max( worst, std::fabs( sv[ j ] - std::pow( 10.0, -j ) ) );
   }
   return worst;
}

// truncate and dense_to_lowrank on matrices with singular values 2^-j: the error in the 2-norm,
// estimated by power iteration, over eps times the largest singular value (dense_to_lowrank may
// add a tenth of eps from the pivoted QR); the rank must be within one of the expected one.
double lowrank_error() {
   double worst = 0;
   for ( const auto& s : std::vector<std::pair<int, int>>{ { 96, 96 }, { 200, 60 }, { 30, 500 } } ) {
      const int m = s.first, n = s.second, k = std::min( m, n );
      // A = X diag(2^-j) Y^T, X and Y with orthonormal columns.
      std::vector<double> X = random( size_t( m ) * k ), Y = random( size_t( n ) * k ), tx, ty;
      householder_qr( m, k, X.data(), tx );
      householder_qr( n, k, Y.data(), ty );
      Mat Xo( m, k ), Yo( n, k );
      for ( int j = 0; j < k; ++j ) {
         Xo( j, j ) = std::pow( 2.0, -j );
         Yo( j, j ) = 1.0;
      }
      apply_q( m, k, X.data(), tx, k, Xo.data(), m );
      apply_q( n, k, Y.data(), ty, k, Yo.data(), n );
      Mat A( m, n );
      gemm( 'N', 'T', m, n, k, 1.0, Xo.data(), m, Yo.data(), n, 0.0, A.data(), m );
      for ( double eps : { 1e-1, 1e-4, 1e-8 } ) {
         for ( int how = 0; how < 2; ++how ) {
            Mat U, V;
            if ( how == 0 ) {
               U = Xo;
               V = Yo;
               truncate( U, V, eps );
            } else {
               Mat P = A;
               dense_to_lowrank( P, eps, U, V );
            }
            Mat R = A;
            gemm( 'N', 'T', m, n, U.n, -1.0, U.data(), m, V.data(), n, 1.0, R.data(), m );
            // Power iteration on R^T R for its largest singular value.
            std::vector<double> x = random( n ), y( m ), z( n );
            double sig = 0;
            for ( int it = 0; it < 60; ++it ) {
               gemm( 'N', 'N', m, 1, n, 1.0, R.data(), m, x.data(), n, 0.0, y.data(), m );
               gemm( 'T', 'N', n, 1, m, 1.0, R.data(), m, y.data(), m, 0.0, z.data(), n );
               double nz = 0;
               for ( double v : z ) {
                  nz += v * v;
               }
               nz = std::sqrt( nz );
               if ( nz == 0 ) {
                  break;
               }
               sig = std::sqrt( nz );
               for ( int i = 0; i < n; ++i ) {
                  x[ i ] = z[ i ] / nz;
               }
            }
            // Expected rank: singular values 2^-j > eps.
            int want = 0;
            while ( want < k && std::pow( 2.0, -want ) > eps ) {
               ++want;
            }
            worst = std::max( worst, sig / eps );
            if ( std::abs( U.n - want ) > 1 ) {
               worst = std::max( worst, 99.0 );
            }
         }
      }
   }
   return worst;
}

}  // namespace

int main() {
   const double eg = gemm_error();
   check( eg < 1e-14, "gemm against loops", eg );
   const double et = trsm_error();
   check( et < 1e-13, "triangular solves", et );
   const double eq = qr_error();
   check( eq < 1e-13, "Householder QR, with and without pivoting", eq );
   const double es = svd_error();
   check( es < 1e-12, "one-sided Jacobi SVD", es );
   const double el = lowrank_error();
   check( el <= 1.15, "truncation error at most about eps times the largest singular value", el );
   return fails ? 1 : 0;
}
