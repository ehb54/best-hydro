// BEST H-matrix solver: mesh features that make GR ill-conditioned (meshcheck.h).
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "meshcheck.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace bh {

namespace {

void sub( const double* a, const double* b, double* d ) {
   for ( int c = 0; c < 3; ++c ) {
      d[ c ] = a[ c ] - b[ c ];
   }
}

double dot( const double* a, const double* b ) { return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ]; }

}  // namespace

MeshDefects mesh_defects( const Mesh& m, const std::vector<double>& x, double dclose ) {
   MeshDefects d;
   const int n = m.n;
   const double kDeg = 180.0 / 3.14159265358979323846;
   std::vector<double> nrm( 3 * size_t( n ), 0.0 );
   for ( int i = 0; i < n; ++i ) {
      const double* p[ 3 ] = { &m.v1[ 3 * size_t( i ) ], &m.v2[ 3 * size_t( i ) ], &m.v3[ 3 * size_t( i ) ] };
      double e1[ 3 ], e2[ 3 ];
      sub( p[ 1 ], p[ 0 ], e1 );
      sub( p[ 2 ], p[ 0 ], e2 );
      double* q = &nrm[ 3 * size_t( i ) ];
      q[ 0 ] = e1[ 1 ] * e2[ 2 ] - e1[ 2 ] * e2[ 1 ];
      q[ 1 ] = e1[ 2 ] * e2[ 0 ] - e1[ 0 ] * e2[ 2 ];
      q[ 2 ] = e1[ 0 ] * e2[ 1 ] - e1[ 1 ] * e2[ 0 ];
      const double l = std::sqrt( dot( q, q ) );
      for ( int c = 0; c < 3; ++c ) {
         q[ c ] = l > 0 ? q[ c ] / l : 0.0;
      }
      // The smallest angle, from the edge vectors at each corner.
      double amin = 180;
      for ( int k = 0; k < 3; ++k ) {
         double a[ 3 ], b[ 3 ];
         sub( p[ ( k + 1 ) % 3 ], p[ k ], a );
         sub( p[ ( k + 2 ) % 3 ], p[ k ], b );
         const double la = std::sqrt( dot( a, a ) ), lb = std::sqrt( dot( b, b ) );
         if ( la > 0 && lb > 0 ) {
            amin = std::min( amin, std::acos( std::max( -1.0, std::min( 1.0, dot( a, b ) / ( la * lb ) ) ) ) * kDeg );
         } else {
            amin = 0;
         }
      }
      if ( amin < 1.0 ) {
         ++d.n_sliver;
      }
      if ( amin < d.min_angle ) {
         d.min_angle = amin;
         d.min_angle_tri = i;
      }
   }
   // Pairs of collocation points closer than dclose: a grid of cells of that size.
   if ( dclose <= 0 ) {
      return d;
   }
   // Cell indices packed exactly in 21 bits each (cells of 0.05 A cover +-52,000 A).
   auto cell = [&]( const double* p, int c ) { return long( std::floor( p[ c ] / dclose ) ); };
   auto key = [&]( long a, long b, long c ) {
      const long off = 1L << 20, mask = ( 1L << 21 ) - 1;
      return ( ( ( a + off ) & mask ) << 42 ) | ( ( ( b + off ) & mask ) << 21 ) | ( ( c + off ) & mask );
   };
   std::unordered_multimap<long, int> grid;
   grid.reserve( size_t( n ) );
   for ( int i = 0; i < n; ++i ) {
      const double* p = &x[ 3 * size_t( i ) ];
      grid.emplace( key( cell( p, 0 ), cell( p, 1 ), cell( p, 2 ) ), i );
   }
   for ( int i = 0; i < n; ++i ) {
      const double* p = &x[ 3 * size_t( i ) ];
      const long c0 = cell( p, 0 ), c1 = cell( p, 1 ), c2 = cell( p, 2 );
      for ( long a = c0 - 1; a <= c0 + 1; ++a ) {
         for ( long b = c1 - 1; b <= c1 + 1; ++b ) {
            for ( long c = c2 - 1; c <= c2 + 1; ++c ) {
               auto r = grid.equal_range( key( a, b, c ) );
               for ( auto it = r.first; it != r.second; ++it ) {
                  const int j = it->second;
                  if ( j <= i ) {
                     continue;
                  }
                  const double* q = &x[ 3 * size_t( j ) ];
                  double dv[ 3 ];
                  sub( p, q, dv );
                  const double dist = std::sqrt( dot( dv, dv ) );
                  if ( dist >= dclose || dot( &nrm[ 3 * size_t( i ) ], &nrm[ 3 * size_t( j ) ] ) >= -0.5 ) {
                     continue;
                  }
                  ++d.n_close;
                  if ( d.close_a < 0 || dist < d.close_min ) {
                     d.close_min = dist;
                     d.close_a = i;
                     d.close_b = j;
                  }
               }
            }
         }
      }
   }
   return d;
}

}  // namespace bh
