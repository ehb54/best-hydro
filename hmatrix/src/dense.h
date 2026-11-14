// BEST H-matrix solver: the small dense kernels the solver needs, without BLAS or LAPACK.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// The solver makes millions of small dense calls from every thread at once.  OpenBLAS takes a
// global lock in its buffer allocation on each call, and with 64 threads that made the
// prototype 30 to 40 times slower.  These kernels keep their scratch space per thread and take
// no locks.  All matrices are column-major.
#pragma once

#include <vector>

namespace bh {

// C = alpha op(A) op(B) + beta C, as BLAS dgemm; op is 'N' or 'T'.  op(A) is m x k, op(B) k x n.
void gemm( char ta, char tb, int m, int n, int k, double alpha, const double* A, int lda, const double* B,
           int ldb, double beta, double* C, int ldc );

// Triangular solves, as BLAS dtrsm with alpha = 1, in the three forms the solver uses:
//   trsm_llnu: B (m x n) <- L^-1 B, L m x m unit lower triangular;
//   trsm_lunn: B (m x n) <- U^-1 B, U m x m upper triangular;
//   trsm_runn: B (m x n) <- B U^-1, U n x n upper triangular.
void trsm_llnu( int m, int n, const double* L, int ldl, double* B, int ldb );
void trsm_lunn( int m, int n, const double* U, int ldu, double* B, int ldb );
void trsm_runn( int m, int n, const double* U, int ldu, double* B, int ldb );

// Householder QR of A (m x k, leading dimension m), LAPACK's compact form: R in the upper
// triangle of the first q = min(m, k) rows, the reflectors below the diagonal, their scalars in
// tau (q).  With pivot, columns are exchanged so that each step takes the largest remaining
// column (perm[ j ] = the original column now at j), and the factorization stops after q steps
// once the columns left have a Frobenius norm of at most stop; q is returned.
int householder_qr( int m, int k, double* A, std::vector<double>& tau, bool pivot = false, double stop = 0,
                    std::vector<int>* perm = nullptr );

// C (m x n, leading dimension ldc) <- Q C, Q = H_0 H_1 ... H_{q-1} from householder_qr of an m-row
// matrix A (leading dimension m).
void apply_q( int m, int q, const double* A, const std::vector<double>& tau, int n, double* C, int ldc );

// One-sided Jacobi SVD of A (m x n, m >= n, leading dimension m, destroyed): on return A holds
// the left singular vectors times the singular values, U diag(s), V (n x n) the right singular
// vectors and s the singular values, all sorted by decreasing s.
void jacobi_svd( int m, int n, double* A, std::vector<double>& V, std::vector<double>& s );

}  // namespace bh
