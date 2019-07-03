#ifndef HYPJACL2_HPP
#define HYPJACL2_HPP

#include "defines.hpp"

EXTERN_C int // 0 if OK, < 0 if invalid argument, > 0 if error
hypJacL2
(const unsigned routine,     // IN, routine ID, <= 15, (B_FS)_2
 // B: block-oriented or full-block,
 // F: full (U \Sigma V^T) or partial (U \Sigma) SVD,
 // S: innermost (32x32 blocks) sorting of eigenvalues.
 const unsigned nrow,        // IN, number of rows of G, == 0 (mod 256)
 const unsigned ncol,        // IN, number of columns of G, <= nrow, == 0 (mod 128)
 const unsigned nplus,       // IN, number of positive eigenvalues of G J G^T, <= ncol
 double *const hG,           // INOUT, ldhG x ncol host array in Fortran order,
 // IN: factor G, OUT: U \Sigma of G = U \Sigma V^T
 const unsigned ldhG,        // IN, leading dimension of G, >= nrow
 double *const hV,           // OUT, optional, ldhV x ncol host array in Fortran order,
 // V^{-T} of G = U \Sigma V^T
 const unsigned ldhV,        // IN, optional, leading dimension of V^{-T}, >= ncol
 double *const hD,           // OUT, eigenvalues of G J G^T, optionally sorted in descending order
 unsigned *const glbSwp,     // OUT, number of sweeps at the outermost level
 unsigned long long *const glb_s, // OUT, number of rotations
 unsigned long long *const glb_b, // OUT, number of ``big'' rotations
 double *const timing        // OUT, optional, in seconds, double[4] ==
 // WALL, SETUP & HOST ==> GPUs, COMPUTATION, CLEANUP & GPUs ==> HOST
) throw();

#endif // !HYPJACL2_HPP
