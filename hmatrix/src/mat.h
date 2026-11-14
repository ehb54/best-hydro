// BEST H-matrix solver: dense column-major matrices and the compression of low-rank products.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#pragma once

#include <cstddef>
#include <vector>

#include "dense.h"

namespace bh {

// A dense column-major m x n matrix.
struct Mat {
   int m = 0, n = 0;
   std::vector<double> a;

   Mat() = default;
   Mat( int rows, int cols ) : m( rows ), n( cols ), a( size_t( rows ) * size_t( cols ), 0.0 ) {}

   double* data() { return a.data(); }
   const double* data() const { return a.data(); }
   double& operator()( int i, int j ) { return a[ size_t( j ) * size_t( m ) + size_t( i ) ]; }
   double operator()( int i, int j ) const { return a[ size_t( j ) * size_t( m ) + size_t( i ) ]; }
   double* col( int j ) { return a.data() + size_t( j ) * size_t( m ); }
   const double* col( int j ) const { return a.data() + size_t( j ) * size_t( m ); }
   size_t bytes() const { return a.size() * sizeof( double ); }
   void resize( int rows, int cols ) {
      m = rows;
      n = cols;
      a.assign( size_t( rows ) * size_t( cols ), 0.0 );
   }
   // Keeps the first `cols` columns.
   void keep_cols( int cols ) {
      n = cols;
      a.resize( size_t( m ) * size_t( cols ) );
      a.shrink_to_fit();
   }
   void release() {
      m = n = 0;
      std::vector<double>().swap( a );
   }
};

// Truncates the low-rank product U V^T (U m x k, V n x k) to the rank r at which the singular values
// left out are all at most eps times the largest one; U and V are replaced by m x r and n x r factors.
void truncate( Mat& U, Mat& V, double eps );
// The same for a dense m x n matrix P (destroyed): P ~ U V^T.  Uses QR with column pivoting, which
// stops once the columns left are below a tenth of the tolerance, then an SVD of the triangle.
void dense_to_lowrank( Mat& P, double eps, Mat& U, Mat& V );

// Calls of truncate() and dense_to_lowrank() and the CPU seconds spent in them (all threads) since
// the last reset.
void lowrank_stats( long& calls, double& seconds, bool reset );
// CPU seconds used by this process (all threads).
double cpu_seconds();

// Peak and current resident memory of this process in bytes (Linux /proc; 0 elsewhere).
size_t peak_rss();
size_t current_rss();
// MemAvailable from /proc/meminfo in bytes (Linux; 0 elsewhere).
size_t mem_available();

}  // namespace bh
