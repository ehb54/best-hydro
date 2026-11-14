// BEST H-matrix solver: cluster tree over the triangles.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// Clusters are split at the median incenter along the longest side of their incenters' bounding
// box, down to a depth at which every cluster has at most `leaf` triangles.  Every cluster at one
// level is split, so all leaves are at the same depth and clusters at a level differ in size by at
// most one (a level-consistent tree; the H-LU relies on it).  A triangle's three unknowns stay
// together: cluster c holds unknowns 3 beg .. 3 (beg + n) - 1 in tree order.
#pragma once

#include <vector>

#include "kernel.h"

namespace bh {

struct Cluster {
   int beg = 0, n = 0;  // triangles perm[ beg .. beg + n )
   int level = 0;
   int son[ 2 ] = { -1, -1 };
   double lo[ 3 ] = { 0, 0, 0 }, hi[ 3 ] = { 0, 0, 0 };  // bounding box of the triangles' vertices
   int leaf_first = 0, leaf_count = 0;  // the leaf clusters below, as a range of ClusterTree::leaves
   bool is_leaf() const { return son[ 0 ] < 0; }
};

struct ClusterTree {
   std::vector<Cluster> c;   // c[ 0 ] is the root
   std::vector<int> perm;    // tree position -> triangle
   std::vector<int> pos;     // triangle -> tree position
   std::vector<int> leaves;  // leaf clusters in tree order
   int depth = 0;
};

ClusterTree build_cluster_tree( const Mesh& mesh, const std::vector<double>& x, int leaf );

// Standard admissibility: min( diam t, diam s ) <= eta dist( t, s ), with bounding boxes.
bool admissible( const Cluster& t, const Cluster& s, double eta );

}  // namespace bh
