// BEST H-matrix solver: batched restarted GMRES with right preconditioning.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "gmres.h"

#include <algorithm>
#include <cmath>

#include "mat.h"
#include "parallel.h"

namespace bh {

namespace {

double norm( const double* x, int n ) {
   double s = 0;
   for ( int i = 0; i < n; ++i ) {
      s += x[ i ] * x[ i ];
   }
   return std::sqrt( s );
}

}  // namespace

GmresResult gmres( int n, int p, const Apply& A, const Precond& Minv, const double* B, int ldb, double* X, int ldx,
                   double tol, int m, int maxit ) {
   GmresResult res;
   res.iters.assign( p, 0 );
   res.relres.assign( p, 0.0 );
   std::vector<double> bnorm( p );
   std::vector<char> done( p, 0 );
   Mat R( n, p );
   for ( int c = 0; c < p; ++c ) {
      for ( int i = 0; i < n; ++i ) {
         X[ size_t( c ) * ldx + i ] = 0.0;
         R( i, c ) = B[ size_t( c ) * ldb + i ];
      }
      bnorm[ c ] = norm( R.col( c ), n );
      res.relres[ c ] = bnorm[ c ] > 0 ? 1.0 : 0.0;
      done[ c ] = bnorm[ c ] == 0;
   }
   int total = 0;
   while ( total < maxit ) {
      std::vector<int> act;
      for ( int c = 0; c < p; ++c ) {
         if ( !done[ c ] ) {
            act.push_back( c );
         }
      }
      if ( act.empty() ) {
         break;
      }
      const int q = int( act.size() );
      std::vector<Mat> Vk( q ), Hk( q );
      std::vector<std::vector<double>> g( q ), cs( q ), sn( q );
      std::vector<int> jlast( q, -1 );
      std::vector<char> going( q, 1 );
      for ( int a = 0; a < q; ++a ) {
         const int c = act[ a ];
         Vk[ a ].resize( n, m + 1 );
         Hk[ a ].resize( m + 1, m );
         g[ a ].assign( m + 1, 0.0 );
         cs[ a ].assign( m, 0.0 );
         sn[ a ].assign( m, 0.0 );
         const double beta = norm( R.col( c ), n );
         for ( int i = 0; i < n; ++i ) {
            Vk[ a ]( i, 0 ) = R( i, c ) / beta;
         }
         g[ a ][ 0 ] = beta;
      }
      for ( int j = 0; j < m && total < maxit; ++j ) {
         std::vector<int> cyc;
         for ( int a = 0; a < q; ++a ) {
            if ( going[ a ] ) {
               cyc.push_back( a );
            }
         }
         if ( cyc.empty() ) {
            break;
         }
         const int qc = int( cyc.size() );
         Mat Z( n, qc ), W( n, qc );
         for ( int b = 0; b < qc; ++b ) {
            std::copy( Vk[ cyc[ b ] ].col( j ), Vk[ cyc[ b ] ].col( j ) + n, Z.col( b ) );
         }
         Minv( Z.data(), n, qc );
         A( Z.data(), n, qc, W.data(), n );
         ++total;
         ++res.steps;
         parallel_for( qc, [&]( long b, int ) {
            const int a = cyc[ b ];
            Mat& V = Vk[ a ];
            Mat& H = Hk[ a ];
            double* w = W.col( int( b ) );
            for ( int i = 0; i <= j; ++i ) {
               const double* vi = V.col( i );
               double h = 0;
               for ( int r = 0; r < n; ++r ) {
                  h += w[ r ] * vi[ r ];
               }
               H( i, j ) = h;
               for ( int r = 0; r < n; ++r ) {
                  w[ r ] -= h * vi[ r ];
               }
            }
            const double hn = norm( w, n );
            H( j + 1, j ) = hn;
            if ( hn > 0 ) {
               double* vn = V.col( j + 1 );
               for ( int r = 0; r < n; ++r ) {
                  vn[ r ] = w[ r ] / hn;
               }
            }
            for ( int i = 0; i < j; ++i ) {
               const double t = cs[ a ][ i ] * H( i, j ) + sn[ a ][ i ] * H( i + 1, j );
               H( i + 1, j ) = -sn[ a ][ i ] * H( i, j ) + cs[ a ][ i ] * H( i + 1, j );
               H( i, j ) = t;
            }
            const double rr = std::hypot( H( j, j ), H( j + 1, j ) );
            cs[ a ][ j ] = rr > 0 ? H( j, j ) / rr : 1.0;
            sn[ a ][ j ] = rr > 0 ? H( j + 1, j ) / rr : 0.0;
            H( j, j ) = rr;
            H( j + 1, j ) = 0.0;
            g[ a ][ j + 1 ] = -sn[ a ][ j ] * g[ a ][ j ];
            g[ a ][ j ] = cs[ a ][ j ] * g[ a ][ j ];
            jlast[ a ] = j;
            ++res.iters[ act[ a ] ];
            if ( std::fabs( g[ a ][ j + 1 ] ) <= tol * bnorm[ act[ a ] ] || hn == 0 ) {
               going[ a ] = 0;
            }
         } );
      }
      // x += M^-1 V y for every column of this cycle.
      Mat T( n, q );
      for ( int a = 0; a < q; ++a ) {
         const int k = jlast[ a ] + 1;
         std::vector<double> y( k );
         for ( int i = k - 1; i >= 0; --i ) {
            double s = g[ a ][ i ];
            for ( int l = i + 1; l < k; ++l ) {
               s -= Hk[ a ]( i, l ) * y[ l ];
            }
            y[ i ] = s / Hk[ a ]( i, i );
         }
         gemm( 'N', 'N', n, 1, k, 1.0, Vk[ a ].data(), n, y.data(), std::max( 1, k ), 0.0, T.col( a ), n );
      }
      Vk.clear();
      Minv( T.data(), n, q );
      Mat Xa( n, q ), Ya( n, q );
      for ( int a = 0; a < q; ++a ) {
         double* xc = X + size_t( act[ a ] ) * ldx;
         for ( int i = 0; i < n; ++i ) {
            xc[ i ] += T( i, a );
            Xa( i, a ) = xc[ i ];
         }
      }
      A( Xa.data(), n, q, Ya.data(), n );
      for ( int a = 0; a < q; ++a ) {
         const int c = act[ a ];
         for ( int i = 0; i < n; ++i ) {
            R( i, c ) = B[ size_t( c ) * ldb + i ] - Ya( i, a );
         }
         res.relres[ c ] = norm( R.col( c ), n ) / bnorm[ c ];
         if ( res.relres[ c ] <= tol ) {
            done[ c ] = 1;
         }
      }
   }
   return res;
}

}  // namespace bh
