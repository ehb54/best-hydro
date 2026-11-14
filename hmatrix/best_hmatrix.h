/* best_hmatrix.h: a hierarchical-matrix solver for BEST's boundary-element system.
 *
 * BEST solves GR x = b, where GR is the 3N x 3N matrix of the Oseen tensor integrated over
 * N triangles (tgetgr.f, triangleint.f).  Stored dense, GR takes 72 N^2 bytes: about 58,000
 * triangles on a 256 GB node.  This library computes the same entries (a port of those two
 * routines that matches them bit for bit), compresses the blocks of GR between well-separated
 * groups of triangles to low rank (adaptive cross approximation, then QR and SVD) to a relative
 * accuracy eps, and solves with GMRES preconditioned by a coarse H-LU factorization, or with a
 * direct H-LU factorization.  Memory and time grow about as N log N.  The blocks can be stored
 * in single precision, in about half the memory; the arithmetic stays in double precision.
 *
 * BEST calls it from best/besolve.f90 through bind(C).  There is one solver per process.
 * Vectors and indices are in BEST's order: unknown 3 n + a (counting from 0) is component a of
 * the traction on triangle n.  Messages (progress, statistics, errors) are queued and read with
 * bhm_log_line; BEST writes them to its .log file.
 *
 * Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
 */
#ifndef BEST_HMATRIX_H
#define BEST_HMATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

/* 1: the solver is there (BEST built without it links stand-ins that return 0). */
int bhm_available( void );

/* Sets up the solver for a triangulation.  v1, v2, v3 are BEST's vertex arrays V1, V2, V3
 * (n x 3, column-major, after the centroid shift), area the triangle areas after the area
 * correction, afac the area correction factor and cg BEST's control value (1, or 0 with -p).
 * threads < 1: OMP_NUM_THREADS if set, otherwise every hardware thread.  Returns 0, or nonzero
 * after queuing a message. */
int bhm_create( int n, const double* v1, const double* v2, const double* v3, const double* area, double afac,
                double cg, int threads );

/* The incenters of the triangles, as TGETGR computes them: x is n x 3, column-major. */
void bhm_incenters( double* x );

/* Entries of GR (before any deflation): rows i0 .. i0 + m - 1 and columns j0 .. j0 + nc - 1
 * into a, column-major with leading dimension lda. */
void bhm_entries( int i0, int m, int j0, int nc, double* a, int lda );

/* The mean of GR's diagonal entries, summed in the order BEST uses. */
double bhm_diag_mean( void );

/* Builds the H-matrix of GR + gdef xmu xmu^T at accuracy eps and prepares the solves: with
 * direct = 0 a preconditioner for GMRES (an H-LU factorization at eps_pre, with the near field
 * compressed too), with direct = 1 the H-LU factorization at eps.  xmu has 3n entries.  With
 * single = 1 the H-matrix and its factors are stored in single precision.  Returns 0, or nonzero
 * after queuing a message (for example when the memory needed is not available). */
int bhm_build( double eps, int direct, double eps_pre, double gdef, const double* xmu, int single );

/* Solves for nrhs right-hand sides at once: b is 3n x nrhs, column-major, and holds the
 * solutions on return.  Returns 0, or nonzero after queuing a message (GMRES did not converge). */
int bhm_solve( int nrhs, double* b );

/* Copies the next queued message into buf (at most len characters, without a terminating NUL)
 * and returns its length, or returns -1 when there is none. */
int bhm_log_line( char* buf, int len );

/* Frees everything. */
void bhm_destroy( void );

#ifdef __cplusplus
}
#endif

#endif
