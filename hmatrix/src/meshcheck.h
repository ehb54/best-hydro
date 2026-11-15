// BEST H-matrix solver: mesh features that make GR ill-conditioned, for the log when the
// preconditioner or GMRES has trouble.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// On probe-0 unions of atoms (and some probe-1.5 meshes with narrow tunnels) sesmesh leaves sliver
// triangles folded onto their neighbors and sheets facing each other across gaps of hundredths of
// an A.  Two collocation points that close, on triangles with opposite normals, give GR rows that
// are nearly the same: a near-null traction mode (a stress trapped in the fold or gap) that exerts
// no net force, so the results barely depend on it, but the coarse preconditioner's factorization
// breaks down on such meshes.
#pragma once

#include <vector>

#include "kernel.h"

namespace bh {

struct MeshDefects {
   double min_angle = 180;   // the smallest angle of any triangle, degrees
   int min_angle_tri = -1;   // its triangle (from 0)
   long n_sliver = 0;        // triangles with an angle below 1 degree
   long n_close = 0;         // pairs with opposite normals (dot < -0.5) and collocation points closer than dclose
   double close_min = 0;     // the closest such pair's distance (A) and triangles (from 0)
   int close_a = -1, close_b = -1;
};

// x: the collocation points (incenters), 3 per triangle.
MeshDefects mesh_defects( const Mesh& m, const std::vector<double>& x, double dclose );

}  // namespace bh
