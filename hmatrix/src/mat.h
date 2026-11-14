// BEST H-matrix solver: dense column-major matrices and the compression of low-rank products.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include "dense.h"

namespace bh {

// A dense column-major m x n matrix of double (Mat) or float (MatF, the H-matrix in single
// precision).
template <class T>
struct MatT {
   int m = 0, n = 0;
   std::vector<T> a;

   MatT() = default;
   MatT( int rows, int cols ) : m( rows ), n( cols ), a( size_t( rows ) * size_t( cols ), T( 0 ) ) {}
   // A copy in the other precision.
   template <class S, class = typename std::enable_if<!std::is_same<S, T>::value>::type>
   explicit MatT( const MatT<S>& o ) : m( o.m ), n( o.n ), a( o.a.begin(), o.a.end() ) {}

   T* data() { return a.data(); }
   const T* data() const { return a.data(); }
   T& operator()( int i, int j ) { return a[ size_t( j ) * size_t( m ) + size_t( i ) ]; }
   T operator()( int i, int j ) const { return a[ size_t( j ) * size_t( m ) + size_t( i ) ]; }
   T* col( int j ) { return a.data() + size_t( j ) * size_t( m ); }
   const T* col( int j ) const { return a.data() + size_t( j ) * size_t( m ); }
   size_t bytes() const { return a.size() * sizeof( T ); }
   void resize( int rows, int cols ) {
      m = rows;
      n = cols;
      a.assign( size_t( rows ) * size_t( cols ), T( 0 ) );
   }
   // Keeps the first `cols` columns.
   void keep_cols( int cols ) {
      n = cols;
      a.resize( size_t( m ) * size_t( cols ) );
      a.shrink_to_fit();
   }
   void release() {
      m = n = 0;
      std::vector<T>().swap( a );
   }
};

using Mat = MatT<double>;
using MatF = MatT<float>;

// d <- s, in d's precision; s is left empty.
inline void store( Mat& d, Mat&& s ) { d = std::move( s ); }
inline void store( MatF& d, Mat&& s ) {
   d.m = s.m;
   d.n = s.n;
   d.a.assign( s.a.begin(), s.a.end() );
   s.release();
}

// A double copy of s.
inline Mat to_double( const Mat& s ) { return s; }
inline Mat to_double( const MatF& s ) { return Mat( s ); }

// Runs f on M in double: on M itself, or on a double copy that is stored back into M afterwards.
template <class F>
void in_double( Mat& M, F&& f ) {
   f( M );
}
template <class F>
void in_double( MatF& M, F&& f ) {
   Mat W( M );
   f( W );
   store( M, std::move( W ) );
}

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
