// BEST H-matrix solver: truncation of low-rank products by QR and SVD (dense.h), memory figures.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "mat.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>

#include <sys/resource.h>

namespace bh {

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

namespace {

// dense_to_lowrank without the timer, for truncate too.
void compress( Mat& P, double eps, Mat& U, Mat& V ) {
   const int m = P.m, n = P.n;
   double cmax = 0;
   for ( int j = 0; j < n; ++j ) {
      const double* c = P.col( j );
      double s = 0;
      for ( int i = 0; i < m; ++i ) {
         s += c[ i ] * c[ i ];
      }
      cmax = std::max( cmax, s );
   }
   // P Pi = Q R, stopped once the columns left are below eps / 10 of the largest column (and so of
   // the largest singular value); then R = W S Z^T from the SVD of R^T = (Z S) W^T.
   std::vector<double> tau;
   std::vector<int> perm;
   const int q = m > 0 && n > 0 ? householder_qr( m, n, P.data(), tau, true, 0.1 * eps * std::sqrt( cmax ), &perm ) : 0;
   if ( q == 0 ) {
      U.resize( m, 0 );
      V.resize( n, 0 );
      return;
   }
   Mat T( n, q );
   for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i <= std::min( j, q - 1 ); ++i ) {
         T( j, i ) = P( i, j );
      }
   }
   std::vector<double> W, s;
   jacobi_svd( n, q, T.data(), W, s );
   const int r = kept_rank( s, eps );
   U.resize( m, r );
   V.resize( n, r );
   for ( int j = 0; j < r; ++j ) {
      for ( int i = 0; i < q; ++i ) {
         U( i, j ) = W[ size_t( j ) * q + i ];
      }
      for ( int i = 0; i < n; ++i ) {
         V( perm[ i ], j ) = T( i, j );
      }
   }
   apply_q( m, q, P.data(), tau, r, U.data(), m );
}

}  // namespace

void dense_to_lowrank( Mat& P, double eps, Mat& U, Mat& V ) {
   LrTimer timer;
   compress( P, eps, U, V );
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
   // U = Qu Ru, V = Qv Rv (reflectors kept in U and V), then the SVD of Ru Rv^T.
   std::vector<double> tu, tv;
   const int qu = householder_qr( m, k, U.data(), tu );
   const int qv = householder_qr( n, k, V.data(), tv );
   Mat Ru( qu, k ), Rv( qv, k );
   for ( int j = 0; j < k; ++j ) {
      for ( int i = 0; i <= std::min( j, qu - 1 ); ++i ) {
         Ru( i, j ) = U( i, j );
      }
      for ( int i = 0; i <= std::min( j, qv - 1 ); ++i ) {
         Rv( i, j ) = V( i, j );
      }
   }
   // M = Ru Rv^T ~ Um Vm^T by the pivoted QR and SVD of dense_to_lowrank; then U = Qu Um, V = Qv Vm.
   Mat M( qu, qv ), Um, Vm;
   gemm( 'N', 'T', qu, qv, k, 1.0, Ru.data(), qu, Rv.data(), qv, 0.0, M.data(), qu );
   compress( M, eps, Um, Vm );
   const int r = Um.n;
   Mat Un( m, r ), Vn( n, r );
   for ( int j = 0; j < r; ++j ) {
      std::copy( Um.col( j ), Um.col( j ) + qu, Un.col( j ) );
      std::copy( Vm.col( j ), Vm.col( j ) + qv, Vn.col( j ) );
   }
   apply_q( m, qu, U.data(), tu, r, Un.data(), m );
   apply_q( n, qv, V.data(), tv, r, Vn.data(), n );
   U = std::move( Un );
   V = std::move( Vn );
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
