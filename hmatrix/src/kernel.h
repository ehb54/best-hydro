// BEST H-matrix solver: the triangulation and the entries of BEST's matrix GR.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// The entries are a line-by-line port of tgetgr.f (self blocks) and triangleint.f (off-diagonal
// blocks), with the operations in the same order, so that they match the Fortran GR bit for bit
// when both are compiled without fused multiply-adds (-ffp-contract=off).
#pragma once

#include <vector>

namespace bh {

// A triangulation as BEST holds it after its set-up: centroid moved to the origin, areas with
// the area correction applied.
struct Mesh {
   int n = 0;                       // triangles (NAR)
   std::vector<double> v1, v2, v3;  // vertices, 3 coordinates per triangle
   std::vector<double> area;        // corrected triangle areas
   double afac = 1.0;               // area correction factor
};

// BEST's V1, V2, V3 (n x 3, column-major) and AREA as a Mesh.
Mesh mesh_from_best( int n, const double* v1, const double* v2, const double* v3, const double* area, double afac );

// The 3x3 blocks of GR: GR(3 n1 + a, 3 n2 + b) = g[ 3 a + b ], where n1 is the triangle whose
// incenter is the collocation point and n2 the triangle integrated over.
class Kernel {
 public:
   Kernel( const Mesh& mesh, double cg );
   int n() const { return n_; }
   // Incenters, 3 per triangle (TGETGR's X).
   const std::vector<double>& x() const { return x_; }
   void block( int n1, int n2, double* g ) const;
   // The dense block of rows of triangles rt[0..m) and columns of triangles ct[0..n): a 3m x 3n
   // column-major matrix with leading dimension lda.
   void fill( const int* rt, int m, const int* ct, int n, double* A, int lda ) const;
   // The mean of GR's diagonal entries, summed as best.f sums them.
   double diag_mean() const;

 private:
   void pair( int n1, int n2, double* g ) const;  // triangleint
   int n_;
   double cg_;
   const Mesh& mesh_;
   std::vector<double> x_, self_;
};

}  // namespace bh
