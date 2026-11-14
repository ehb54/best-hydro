// BEST H-matrix solver: restarted GMRES for several right-hand sides at once, right-preconditioned.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#pragma once

#include <functional>
#include <vector>

namespace bh {

// Y = A X for q columns (leading dimensions ldx, ldy).
using Apply = std::function<void( const double* X, int ldx, int q, double* Y, int ldy )>;
// X <- M^-1 X for q columns.
using Precond = std::function<void( double* X, int ldx, int q )>;

struct GmresResult {
   std::vector<int> iters;      // iterations per column
   std::vector<double> relres;  // final true relative residual |b - A x| / |b| per column
   int steps = 0;               // batched applications of A M^-1
};

// Solves A X = B (n x p) by GMRES(m) with right preconditioning, X starting at zero.  The columns
// run in lockstep so that A and M^-1 are applied to all unconverged columns at once.  A column stops
// when its residual estimate is at most tol |b|; the true residual is checked at each restart.
GmresResult gmres( int n, int p, const Apply& A, const Precond& Minv, const double* B, int ldb, double* X, int ldx,
                   double tol, int m, int maxit );

}  // namespace bh
