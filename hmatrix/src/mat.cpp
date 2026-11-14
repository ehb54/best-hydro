// BEST H-matrix solver: BLAS/LAPACK wrappers, QR/SVD truncation of low-rank products.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "mat.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>

#include <sys/resource.h>

extern "C" {
void dgemm_( const char*, const char*, const int*, const int*, const int*, const double*, const double*,
             const int*, const double*, const int*, const double*, double*, const int*, size_t, size_t );
void dtrsm_( const char*, const char*, const char*, const char*, const int*, const int*, const double*,
             const double*, const int*, double*, const int*, size_t, size_t, size_t, size_t );
void dgeqrf_( const int*, const int*, double*, const int*, double*, double*, const int*, int* );
void dorgqr_( const int*, const int*, const int*, double*, const int*, const double*, double*, const int*,
              int* );
void dgesdd_( const char*, const int*, const int*, double*, const int*, double*, double*, const int*, double*,
              const int*, double*, const int*, int*, int*, size_t );
void dgesvd_( const char*, const char*, const int*, const int*, double*, const int*, double*, double*,
              const int*, double*, const int*, double*, const int*, int*, size_t, size_t );
void openblas_set_num_threads( int ) __attribute__( ( weak ) );
}

namespace bh {

void gemm( char ta, char tb, int m, int n, int k, double alpha, const double* A, int lda, const double* B,
           int ldb, double beta, double* C, int ldc ) {
   if ( m <= 0 || n <= 0 ) {
      return;
   }
   if ( k <= 0 ) {
      if ( beta != 1.0 ) {
         for ( int j = 0; j < n; ++j ) {
            for ( int i = 0; i < m; ++i ) {
               C[ size_t( j ) * ldc + i ] *= beta;
            }
         }
      }
      return;
   }
   lda = std::max( 1, lda );
   ldb = std::max( 1, ldb );
   ldc = std::max( 1, ldc );
   dgemm_( &ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1 );
}

void trsm( char side, char uplo, char ta, char diag, int m, int n, double alpha, const double* A, int lda,
           double* B, int ldb ) {
   if ( m <= 0 || n <= 0 ) {
      return;
   }
   lda = std::max( 1, lda );
   ldb = std::max( 1, ldb );
   dtrsm_( &side, &uplo, &ta, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1 );
}

void qr( Mat& A, Mat& R ) {
   const int m = A.m, k = A.n, q = std::min( m, k );
   R.resize( q, k );
   if ( q == 0 ) {
      A.keep_cols( 0 );
      return;
   }
   std::vector<double> tau( q );
   int lwork = -1, info = 0;
   double wq = 0;
   const int lda = std::max( 1, m );
   dgeqrf_( &m, &k, A.data(), &lda, tau.data(), &wq, &lwork, &info );
   lwork = std::max( 1, int( wq ) );
   std::vector<double> work( lwork );
   dgeqrf_( &m, &k, A.data(), &lda, tau.data(), work.data(), &lwork, &info );
   if ( info != 0 ) {
      throw std::runtime_error( "dgeqrf failed" );
   }
   for ( int j = 0; j < k; ++j ) {
      for ( int i = 0; i <= std::min( j, q - 1 ); ++i ) {
         R( i, j ) = A( i, j );
      }
   }
   lwork = -1;
   dorgqr_( &m, &q, &q, A.data(), &lda, tau.data(), &wq, &lwork, &info );
   lwork = std::max( 1, int( wq ) );
   work.resize( lwork );
   dorgqr_( &m, &q, &q, A.data(), &lda, tau.data(), work.data(), &lwork, &info );
   if ( info != 0 ) {
      throw std::runtime_error( "dorgqr failed" );
   }
   A.keep_cols( q );
}

void svd( Mat& A, Mat& U, std::vector<double>& s, Mat& VT ) {
   const int m = A.m, n = A.n, q = std::min( m, n );
   U.resize( m, q );
   VT.resize( q, n );
   s.assign( q, 0.0 );
   if ( q == 0 ) {
      return;
   }
   Mat keep = A;  // dgesdd can fail to converge; dgesvd then works on the original
   const int lda = std::max( 1, m ), ldu = std::max( 1, m ), ldvt = std::max( 1, q );
   int lwork = -1, info = 0;
   double wq = 0;
   std::vector<int> iwork( 8 * size_t( q ) );
   dgesdd_( "S", &m, &n, A.data(), &lda, s.data(), U.data(), &ldu, VT.data(), &ldvt, &wq, &lwork, iwork.data(),
            &info, 1 );
   lwork = std::max( 1, int( wq ) );
   std::vector<double> work( lwork );
   dgesdd_( "S", &m, &n, A.data(), &lda, s.data(), U.data(), &ldu, VT.data(), &ldvt, work.data(), &lwork,
            iwork.data(), &info, 1 );
   if ( info == 0 ) {
      return;
   }
   A = keep;
   lwork = -1;
   dgesvd_( "S", "S", &m, &n, A.data(), &lda, s.data(), U.data(), &ldu, VT.data(), &ldvt, &wq, &lwork, &info, 1,
            1 );
   lwork = std::max( 1, int( wq ) );
   work.resize( lwork );
   dgesvd_( "S", "S", &m, &n, A.data(), &lda, s.data(), U.data(), &ldu, VT.data(), &ldvt, work.data(), &lwork,
            &info, 1, 1 );
   if ( info != 0 ) {
      throw std::runtime_error( "SVD failed" );
   }
}

namespace {

int kept_rank( const std::vector<double>& s, double eps ) {
   if ( s.empty() || !( s[ 0 ] > 0 ) ) {
      return 0;
   }
   const double cut = eps * s[ 0 ];
   int r = int( s.size() );
   while ( r > 0 && s[ r - 1 ] <= cut ) {
      --r;
   }
   return r;
}

}  // namespace

namespace {

std::atomic<long> g_lr_calls( 0 ), g_lr_ns( 0 );

struct LrTimer {
   std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
   ~LrTimer() {
      ++g_lr_calls;
      g_lr_ns += long( std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now() - t0 ).count() );
   }
};

}  // namespace

void lowrank_stats( long& calls, double& seconds, bool reset ) {
   calls = g_lr_calls;
   seconds = 1e-9 * double( g_lr_ns );
   if ( reset ) {
      g_lr_calls = 0;
      g_lr_ns = 0;
   }
}

double cpu_seconds() {
   rusage ru;
   getrusage( RUSAGE_SELF, &ru );
   return double( ru.ru_utime.tv_sec + ru.ru_stime.tv_sec ) + 1e-6 * double( ru.ru_utime.tv_usec + ru.ru_stime.tv_usec );
}

void truncate( Mat& U, Mat& V, double eps ) {
   LrTimer timer;
   const int m = U.m, n = V.m, k = U.n;
   if ( k == 0 ) {
      return;
   }
   if ( V.n != k ) {
      throw std::runtime_error( "truncate: rank mismatch" );
   }
   Mat Ru, Rv;
   qr( U, Ru );  // U -> Qu (m x ku), Ru ku x k
   qr( V, Rv );  // V -> Qv (n x kv), Rv kv x k
   const int ku = U.n, kv = V.n;
   Mat M( ku, kv );
   gemm( 'N', 'T', ku, kv, k, 1.0, Ru.data(), ku, Rv.data(), kv, 0.0, M.data(), ku );
   Mat W, ZT;
   std::vector<double> s;
   svd( M, W, s, ZT );
   const int r = kept_rank( s, eps );
   for ( int j = 0; j < r; ++j ) {
      for ( int i = 0; i < W.m; ++i ) {
         W( i, j ) *= s[ j ];
      }
   }
   Mat Un( m, r ), Vn( n, r );
   gemm( 'N', 'N', m, r, ku, 1.0, U.data(), m, W.data(), W.m, 0.0, Un.data(), m );
   gemm( 'N', 'T', n, r, kv, 1.0, V.data(), n, ZT.data(), ZT.m, 0.0, Vn.data(), n );
   U = std::move( Un );
   V = std::move( Vn );
}

void dense_to_lowrank( Mat& P, double eps, Mat& U, Mat& V ) {
   LrTimer timer;
   const int m = P.m, n = P.n;
   Mat W, ZT;
   std::vector<double> s;
   svd( P, W, s, ZT );
   const int r = kept_rank( s, eps );
   U.resize( m, r );
   V.resize( n, r );
   for ( int j = 0; j < r; ++j ) {
      for ( int i = 0; i < m; ++i ) {
         U( i, j ) = W( i, j ) * s[ j ];
      }
      for ( int i = 0; i < n; ++i ) {
         V( i, j ) = ZT( j, i );
      }
   }
}

void set_blas_threads( int n ) {
   if ( openblas_set_num_threads ) {
      openblas_set_num_threads( std::max( 1, n ) );
   }
}

namespace {

size_t proc_status_kb( const char* key ) {
   std::ifstream in( "/proc/self/status" );
   std::string line;
   const size_t len = std::strlen( key );
   while ( std::getline( in, line ) ) {
      if ( line.compare( 0, len, key ) == 0 ) {
         return size_t( std::strtoull( line.c_str() + len, nullptr, 10 ) );
      }
   }
   return 0;
}

}  // namespace

size_t peak_rss() { return proc_status_kb( "VmHWM:" ) * 1024; }
size_t current_rss() { return proc_status_kb( "VmRSS:" ) * 1024; }

size_t mem_available() {
   std::ifstream in( "/proc/meminfo" );
   std::string line;
   while ( std::getline( in, line ) ) {
      if ( line.compare( 0, 13, "MemAvailable:" ) == 0 ) {
         return size_t( std::strtoull( line.c_str() + 13, nullptr, 10 ) ) * 1024;
      }
   }
   return 0;
}

}  // namespace bh
