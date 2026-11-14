// BEST H-matrix solver: hierarchical matrix of BEST's GR, its matrix-vector product and H-LU.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// Unknowns are in cluster-tree order (3 per triangle).  The block tree is built from pairs of
// clusters at the same level: an admissible pair is a low-rank leaf (R, U V^T), a pair of leaf
// clusters is a dense leaf (D), anything else is split 2 x 2 (H).  Low-rank blocks are filled by
// adaptive cross approximation on BEST's own entries and recompressed by QR + SVD, keeping the
// singular values above eps times the block's largest.
#pragma once

#include <memory>
#include <vector>

#include "cluster.h"
#include "kernel.h"
#include "mat.h"

namespace bh {

enum class BType { H, R, D };

struct Block {
   int r0 = 0, nr = 0, c0 = 0, nc = 0;  // rows and columns (unknowns, tree order)
   int t = -1, s = -1;                  // row and column clusters
   int level = 0;
   BType type = BType::D;
   std::unique_ptr<Block> ch[ 4 ];      // H: (t0,s0) (t0,s1) (t1,s0) (t1,s1)
   Mat D;                               // D: nr x nc
   Mat U, V;                            // R: U nr x k, V nc x k; the block is U V^T
   Mat Pu, Pv;                          // H during the H-LU: a pending update Pu Pv^T not yet pushed
                                        // down to the children (the block is children + Pu Pv^T)
   int rank() const { return U.n; }
   Block* c( int i, int j ) const { return ch[ 2 * i + j ].get(); }
};

struct HOptions {
   double eps = 1e-8;        // relative accuracy of each low-rank block
   double eta = 2.0;         // admissibility parameter
   double aca_factor = 0.1;  // ACA stops at aca_factor * eps; the SVD then truncates at eps
   long dense_small = 1024;  // admissible blocks of at most this many triangle pairs: fill, then SVD
   bool near = false;        // also compress dense leaves off the diagonal at eps where that is smaller
};

struct HStats {
   long n_dense = 0, n_lowrank = 0;
   double bytes_dense = 0, bytes_lowrank = 0;
   int kmax = 0;
   double kmean = 0;                 // mean rank over low-rank blocks
   std::vector<int> kmax_level;      // largest rank per block level
   std::vector<long> nlr_level;      // low-rank blocks per level
   double bytes() const { return bytes_dense + bytes_lowrank; }
};

class HMatrix {
 public:
   HMatrix( const ClusterTree& ct, std::unique_ptr<Block> root );
   // Builds the block tree and fills it from the kernel.  kernel_calls counts 3x3 blocks evaluated.
   static HMatrix build( const ClusterTree& ct, const Kernel& K, const HOptions& opt, long* kernel_calls );
   // A deep copy; with eps > 0 every low-rank block is truncated again at eps, and with near also
   // every dense leaf off the diagonal becomes low-rank at eps where that takes less memory.
   HMatrix copy( double eps, bool near = false ) const;

   const ClusterTree& tree() const { return *ct_; }
   Block* root() const { return root_.get(); }
   int size() const { return root_->nr; }
   HStats stats() const;

   // Y = A X: X is size() x p with leading dimension ldx, Y likewise.  Parallel; the result does
   // not depend on the number of threads.
   void matvec( const double* X, int ldx, int p, double* Y, int ldy ) const;

   // In-place H-LU without pivoting: afterwards the blocks hold L (unit lower) and U, with every
   // low-rank result truncated at eps.  Returns the smallest |pivot| / largest |pivot| seen.
   double lu( double eps );
   // X <- (L U)^-1 X for p columns, after lu(); all columns at once, parallel over block rows.
   void lu_solve( double* X, int ldx, int p ) const;
   // A <- A + gamma u v^T (u, v in tree order): appended to every low-rank leaf (rank + 1) and added to
   // every dense leaf.
   void add_rank1( double gamma, const double* u, const double* v );
   // Truncation counts and CPU seconds by call site in lu() since the last reset.
   static std::string profile( bool reset );
   // lu() runs the parts of blocks of at least leaves x leaves leaf clusters on separate threads.
   static void set_fork( int leaves );

 private:
   void index();
   const ClusterTree* ct_;
   std::unique_ptr<Block> root_;
   std::vector<Block*> leaves_;
   std::vector<std::vector<int>> rows_;  // per leaf cluster: the leaves (indices) covering its rows
};

}  // namespace bh
