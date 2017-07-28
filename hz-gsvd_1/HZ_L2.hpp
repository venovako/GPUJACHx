#ifndef HZ_L2_HPP
#define HZ_L2_HPP

#include "defines.hpp"

extern int // 0 if OK, < 0 if invalid argument, > 0 if error
HZ_L2
(
 const unsigned routine,     // IN, routine ID, <= 15, (B_FS)_2
 // B: block-oriented or full-block,
 // F: full (U \Sigma V^T) or partial (U \Sigma) SVD,
 // S: innermost (32x32 blocks) sorting of eigenvalues.
 const unsigned nrow,        // IN, number of rows of G, == 0 (mod 256)
 const unsigned ncol,        // IN, number of columns of G, <= nrow, == 0 (mod 128)
 double *const hF,           // INOUT, ldhF x ncol host array in Fortran order,
 const unsigned ldhF,        // IN, leading dimension of F, >= nrow
 double *const hG,           // INOUT, ldhG x ncol host array in Fortran order,
 // IN: factor G, OUT: U \Sigma of G = U \Sigma V^T
 const unsigned ldhG,        // IN, leading dimension of G, >= nrow
 double *const hV,           // OUT, optional, ldhV x ncol host array in Fortran order,
 // V^{-T} of G = U \Sigma V^T
 const unsigned ldhV,        // IN, optional, leading dimension of V^{-T}, >= nrow
 double *const hS,           // OUT, the generalized singular values, optionally sorted in descending order
 double *const hH,           // ||F_i||_2/sqrt(||F_i||_2^2 + ||G_i||_2^2)
 double *const hK,           // ||G_i||_2/sqrt(||F_i||_2^2 + ||G_i||_2^2)
 unsigned *const glbSwp,     // OUT, number of sweeps at the outermost level
 unsigned Long *const glb_s, // OUT, number of rotations
 unsigned Long *const glb_b, // OUT, number of ``big'' rotations
 double *const timing        // OUT, optional, in seconds, double[4] ==
 // WALL, SETUP & HOST ==> GPUs, COMPUTATION, CLEANUP & GPUs ==> HOST
) throw();

#endif // !HZ_L2_HPP
