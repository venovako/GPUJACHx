#include "HypJac.hpp"
#include "HypJacL2.hpp"

#include "device_code.hpp"
#include "cuda_helper.hpp"
#include "cuda_memory_helper.hpp"
#include "my_utils.hpp"

#ifdef ANIMATE
#include "vn_lib.h"
#endif /* ANIMATE */

int // 0 if OK, < 0 if invalid argument, > 0 if error
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
) throw()
{
  long long timers[4] = { 0ll };
  stopwatch_reset(timers[0]);

  if (routine >= 16u)
    return -1;

  const bool
    cdsort = (routine & HYPJAC_CDSORT_1),
    full_svd = (routine & HYPJAC_FULL_SVD),
    blk_ori = (routine & HYPJAC_BLK_ORI);

  if (!nrow || (nrow % 64u))
    return -2;

  if (!ncol || (ncol > nrow) || (ncol % 32u))
    return -3;

  if (nplus > ncol)
    return -4;

  if (!hG)
    return -5;

  if (ldhG < nrow)
    return -6;

  if (full_svd) {
    if (!hV)
      return -7;
    if (ldhV < ncol)
      return -8;
  }

  if (!hD)
    return -9;

  if (!glbSwp)
    return -10;
  if (!glb_s)
    return -11;
  if (!glb_b)
    return -12;

  stopwatch_reset(timers[3]);

  const unsigned
    swp_max[HYPJAC_MAX_LEVELS] = { (blk_ori ? 1u : HYPJAC_NSWEEP), HYPJAC_NSWEEP };

  size_t lddG = static_cast<size_t>(nrow);
  double *const dG = allocDeviceMtx<double>(lddG, static_cast<size_t>(nrow), static_cast<size_t>(ncol), true);

  size_t lddV = static_cast<size_t>(full_svd ? ncol : 0u);
  double *const dV = (full_svd ? allocDeviceMtx<double>(lddV, static_cast<size_t>(ncol), static_cast<size_t>(ncol), true) : static_cast<double*>(NULL));

  double *const dD = allocDeviceVec<double>(static_cast<size_t>(ncol));

  volatile unsigned long long *cvg_dat = static_cast<volatile unsigned long long*>(NULL);
  CUDA_CALL(cudaHostAlloc((void**)&cvg_dat, sizeof(unsigned long long), cudaHostAllocPortable | cudaHostAllocMapped));

  CUDA_CALL(cudaMemcpy2DAsync(dG, lddG * sizeof(double), hG, ldhG * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyHostToDevice));
  if (full_svd)
    CUDA_CALL(cudaMemset2DAsync(dV, lddV * sizeof(double), 0, ncol * sizeof(double), ncol));
  initSymbols(dG, dV, cvg_dat, nrow, ncol, static_cast<unsigned>(lddG), static_cast<unsigned>(lddV), swp_max[0u]);

  CUDA_CALL(cudaDeviceSynchronize());

  timers[1] = stopwatch_lap(timers[3]);

  if (full_svd) {
    initV(dV, 0u, ncol, static_cast<unsigned>(lddV));
    CUDA_CALL(cudaDeviceSynchronize());
  }

  const int definite = (nplus ? ((nplus < ncol) ? 0 : +1) : -1);
  void (*const defJac1)(const unsigned, const int) =
    (cdsort ? (full_svd ? defJacL1sv : defJacL1s) : (full_svd ? defJacL1v : defJacL1));
  void (*const hypJac1)(const unsigned, const unsigned) =
    (cdsort ? (full_svd ? hypJacL1sv : hypJacL1s) : (full_svd ? hypJacL1v : hypJacL1));

  *glb_s = 0ull;
  *glb_b = 0ull;
  long long swp_tim = 0ll;
  stopwatch_reset(swp_tim);

  const unsigned swp = swp_max[1u];
  unsigned blk_swp = 0u;

#ifdef ANIMATE
  vn_mtxvis_ctx *ctx = static_cast<vn_mtxvis_ctx*>(NULL);
  if (ncol < 10000u) {
    char fname[8] = { '\0' };
    (void)sprintf(fname, "%c%x_%04u", (definite ? 'A' : 'H'), routine, ncol);
    SYSI_CALL(vn_mtxvis_start(&ctx, fname, (VN_MTXVIS_OP_AtA | VN_MTXVIS_FN_Lg | VN_MTXVIS_FF_Bin), nrow, ncol, 1, 1, 7));
    if (ctx)
      SYSI_CALL(vn_mtxvis_frame(ctx, hG, ldhG));
  }
#endif /* ANIMATE */

  while (blk_swp < swp) {
    *cvg_dat = 0ul;

    for (unsigned blk_stp = 0u; blk_stp < STRAT1_STEPS; ++blk_stp) {
      if (blk_stp)
        CUDA_CALL(cudaDeviceSynchronize());

      if (definite)
        defJac1(blk_stp, definite);
      else
        hypJac1(blk_stp, nplus);

#ifdef ANIMATE
      if (ctx) {
        CUDA_CALL(cudaDeviceSynchronize());
        CUDA_CALL(cudaMemcpy2DAsync(hG, ldhG * sizeof(double), dG, lddG * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaDeviceSynchronize());
        SYSI_CALL(vn_mtxvis_frame(ctx, hG, ldhG));
      }
#endif /* ANIMATE */
    }

    ++blk_swp;
    CUDA_CALL(cudaDeviceSynchronize());

    const unsigned long long cvg = *cvg_dat;
    const unsigned cvg_s = static_cast<unsigned>(cvg);
    *glb_s += cvg_s;
    const unsigned cvg_b = static_cast<unsigned>(cvg >> 32u);
    *glb_b += cvg_b;
    const double tim_s = stopwatch_lap(swp_tim) * TS2S;
    (void)fprintf(stdout, "BLK_SWP(%2u), ROT_S(%10u), ROT_B(%10u), TIME(%#12.6f s)\n", blk_swp, cvg_s, cvg_b, tim_s);
    (void)fflush(stdout);
    if (!cvg_b)
      break;
  }

#ifdef ANIMATE
  if (ctx)
    SYSI_CALL(vn_mtxvis_stop(ctx));
#endif /* ANIMATE */

  *glbSwp = blk_swp;
  initD(dG, dD, 0u, nrow, ncol, nplus, static_cast<unsigned>(lddG));
  CUDA_CALL(cudaDeviceSynchronize());

  timers[2] = stopwatch_lap(timers[3]);

  CUDA_CALL(cudaMemcpy2DAsync(hG, ldhG * sizeof(double), dG, lddG * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
  if (full_svd)
    CUDA_CALL(cudaMemcpy2DAsync(hV, ldhV * sizeof(double), dV, lddV * sizeof(double), ncol * sizeof(double), ncol, cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpyAsync(hD, dD, ncol * sizeof(double), cudaMemcpyDeviceToHost));

  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaFreeHost((void*)cvg_dat));
  CUDA_CALL(cudaFree(static_cast<void*>(dD)));
  if (dV)
    CUDA_CALL(cudaFree(static_cast<void*>(dV)));
  CUDA_CALL(cudaFree(static_cast<void*>(dG)));

  timers[3] = stopwatch_lap(timers[3]);
  timers[0] = stopwatch_lap(timers[0]);

  if (timing)
    for (unsigned i = 0u; i < 4u; ++i)
      timing[i] = timers[i] * TS2S;

  return 0;
}
