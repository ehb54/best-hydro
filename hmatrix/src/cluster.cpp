// BEST H-matrix solver: cluster tree.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "cluster.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace bh {

namespace {

void split( ClusterTree& ct, const Mesh& mesh, const std::vector<double>& x, int id ) {
   Cluster cl = ct.c[ id ];
   for ( int d = 0; d < 3; ++d ) {
      cl.lo[ d ] = 1e300;
      cl.hi[ d ] = -1e300;
   }
   for ( int k = cl.beg; k < cl.beg + cl.n; ++k ) {
      const int t = ct.perm[ k ];
      for ( const std::vector<double>* v : { &mesh.v1, &mesh.v2, &mesh.v3 } ) {
         for ( int d = 0; d < 3; ++d ) {
            cl.lo[ d ] = std::min( cl.lo[ d ], ( *v )[ 3 * t + d ] );
            cl.hi[ d ] = std::max( cl.hi[ d ], ( *v )[ 3 * t + d ] );
         }
      }
   }
   if ( cl.level == ct.depth ) {
      cl.leaf_first = int( ct.leaves.size() );
      cl.leaf_count = 1;
      ct.leaves.push_back( id );
      ct.c[ id ] = cl;
      return;
   }
   // Split along the longest side of the incenters' box.
   double lo[ 3 ] = { 1e300, 1e300, 1e300 }, hi[ 3 ] = { -1e300, -1e300, -1e300 };
   for ( int k = cl.beg; k < cl.beg + cl.n; ++k ) {
      const int t = ct.perm[ k ];
      for ( int d = 0; d < 3; ++d ) {
         lo[ d ] = std::min( lo[ d ], x[ 3 * t + d ] );
         hi[ d ] = std::max( hi[ d ], x[ 3 * t + d ] );
      }
   }
   int axis = 0;
   for ( int d = 1; d < 3; ++d ) {
      if ( hi[ d ] - lo[ d ] > hi[ axis ] - lo[ axis ] ) {
         axis = d;
      }
   }
   const int half = cl.n / 2;
   int* b = ct.perm.data() + cl.beg;
   std::nth_element( b, b + half, b + cl.n, [&]( int p, int q ) {
      const double xp = x[ 3 * p + axis ], xq = x[ 3 * q + axis ];
      return xp < xq || ( xp == xq && p < q );
   } );
   Cluster s0, s1;
   s0.beg = cl.beg;
   s0.n = half;
   s1.beg = cl.beg + half;
   s1.n = cl.n - half;
   s0.level = s1.level = cl.level + 1;
   const int i0 = int( ct.c.size() );
   ct.c.push_back( s0 );
   ct.c.push_back( s1 );
   cl.son[ 0 ] = i0;
   cl.son[ 1 ] = i0 + 1;
   cl.leaf_first = int( ct.leaves.size() );
   ct.c[ id ] = cl;
   split( ct, mesh, x, i0 );
   split( ct, mesh, x, i0 + 1 );
   ct.c[ id ].leaf_count = int( ct.leaves.size() ) - ct.c[ id ].leaf_first;
}

double diam( const Cluster& t ) {
   double s = 0;
   for ( int d = 0; d < 3; ++d ) {
      s += ( t.hi[ d ] - t.lo[ d ] ) * ( t.hi[ d ] - t.lo[ d ] );
   }
   return std::sqrt( s );
}

double dist( const Cluster& t, const Cluster& u ) {
   double s = 0;
   for ( int d = 0; d < 3; ++d ) {
      const double gap = std::max( 0.0, std::max( t.lo[ d ] - u.hi[ d ], u.lo[ d ] - t.hi[ d ] ) );
      s += gap * gap;
   }
   return std::sqrt( s );
}

}  // namespace

ClusterTree build_cluster_tree( const Mesh& mesh, const std::vector<double>& x, int leaf ) {
   ClusterTree ct;
   const int n = mesh.n;
   ct.perm.resize( n );
   std::iota( ct.perm.begin(), ct.perm.end(), 0 );
   ct.depth = 0;
   while ( ( ( long( n ) + ( 1L << ct.depth ) - 1 ) >> ct.depth ) > leaf ) {
      ++ct.depth;
   }
   Cluster root;
   root.beg = 0;
   root.n = n;
   ct.c.reserve( size_t( 2 ) << ct.depth );
   ct.c.push_back( root );
   split( ct, mesh, x, 0 );
   ct.pos.resize( n );
   for ( int k = 0; k < n; ++k ) {
      ct.pos[ ct.perm[ k ] ] = k;
   }
   return ct;
}

bool admissible( const Cluster& t, const Cluster& s, double eta ) {
   return std::min( diam( t ), diam( s ) ) <= eta * dist( t, s );
}

}  // namespace bh
