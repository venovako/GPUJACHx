#include "HZ.hpp"
#include "HZ_L2.hpp"

#include "device_code.hpp"
#include "cuda_helper.hpp"
#include "cuda_memory_helper.hpp"
#include "my_utils.hpp"

#ifdef ANIMATE
#include "vn_lib.h"
#endif // ANIMATE

int // 0 if OK, < 0 if invalid argument, > 0 if error
HZ_L2
(
 const unsigned routine,     // IN, routine ID, <= 15, (B___)_2
 // B: block-oriented or full-block
 const unsigned nrow,        // IN, number of rows of G, == 0 (mod 256)
 const unsigned ncol,        // IN, number of columns of G, <= nrow, == 0 (mod 128)
 double *const hF,           // INOUT, ldhF x ncol host array in Fortran order,
 const unsigned ldhF,        // IN, leading dimension of F, >= nrow
 double *const hG,           // INOUT, ldhG x ncol host array in Fortran order,
 // IN: factor G, OUT: U \Sigma of G = U \Sigma V^T
 const unsigned ldhG,        // IN, leading dimension of G, >= nrow
 double *const hV,           // OUT, ldhV x ncol host array in Fortran order,
 // V^{-T} of G = U \Sigma V^T
 const unsigned ldhV,        // IN, leading dimension of V^{-T}, >= nrow
 double *const hS,           // OUT, the generalized singular values, optionally sorted in descending order
 double *const hH,           // ||F_i||_2/sqrt(||F_i||_2^2 + ||G_i||_2^2)
 double *const hK,           // ||G_i||_2/sqrt(||F_i||_2^2 + ||G_i||_2^2)
 unsigned *const glbSwp,     // OUT, number of sweeps at the outermost level
 unsigned Long *const glb_s, // OUT, number of rotations
 unsigned Long *const glb_b, // OUT, number of ``big'' rotations
 double *const timing        // OUT, optional, in seconds, double[4] ==
 // WALL, SETUP & HOST ==> GPUs, COMPUTATION, CLEANUP & GPUs ==> HOST
) throw()
{
  Long timers[4] = { MkLong(0) };
  stopwatch_reset(timers[0]);

  if (routine >= 16u)
    return -1;

  const bool blk_ori = (routine & HZ_BLK_ORI);

  if (!nrow || (nrow % 64u))
    return -2;

  if (!ncol || (ncol > nrow) || (ncol % 32u))
    return -3;

  if (!hF)
    return -4;

  if (ldhF < nrow)
    return -5;

  if (!hG)
    return -6;

  if (ldhG < nrow)
    return -7;

  if (!hV)
    return -8;
  if (ldhV < nrow)
    return -9;

  if (!hS)
    return -10;
  if (!hH)
    return -11;
  if (!hK)
    return -12;

  if (!glbSwp)
    return -13;
  if (!glb_s)
    return -14;
  if (!glb_b)
    return -15;

  stopwatch_reset(timers[3]);

  const unsigned
    swp_max[HZ_MAX_LEVELS] = { (blk_ori ? 1u : HZ_NSWEEP), HZ_NSWEEP };

  const double
    alpha = +1.0,
    beta = +0.0;
  const double
    *alpha_ptr = static_cast<double*>(NULL),
    *beta_ptr = static_cast<double*>(NULL);

  size_t ldd = static_cast<size_t>(nrow);
  double *ptr = allocDeviceMtx<double>(ldd, static_cast<size_t>(nrow), static_cast<size_t>(3u * ncol + 3u), true);
  double *const dF = ptr;
  ptr += ldd * ncol;
  double *const dG = ptr;
  ptr += ldd * ncol;
  double *const dV = ptr;
  ptr += ldd * ncol;
  double *const dS = ptr;
  ptr += ldd;
  double *const dH = ptr;
  ptr += ldd;
  double *const dK = ptr;
  ptr = dF;

  volatile unsigned Long *cvg_dat = static_cast<volatile unsigned Long*>(NULL);
  CUDA_CALL(cudaHostAlloc((void**)&cvg_dat, sizeof(unsigned Long), cudaHostAllocPortable | cudaHostAllocMapped));

  CUDA_CALL(cudaMemcpy2DAsync(dF, ldd * sizeof(double), hF, ldhF * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy2DAsync(dG, ldd * sizeof(double), hG, ldhG * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemset2DAsync(dV, ldd * sizeof(double), 0, nrow * sizeof(double), ncol));
  initSymbols(dF, dG, dV, dS, dH, dK, cvg_dat, nrow, ncol, static_cast<unsigned>(ldd), static_cast<unsigned>(ldd), static_cast<unsigned>(ldd), swp_max[0u], alpha, beta, alpha_ptr, beta_ptr);

  CUDA_CALL(cudaDeviceSynchronize());

  timers[1] = stopwatch_lap(timers[3]);
  
  initV(ncol, static_cast<cudaStream_t>(NULL));
  CUDA_CALL(cudaDeviceSynchronize());

  void (*const HZ_L1)(const unsigned) = HZ_L1_sv;

  *glb_s = MkLong(0u);
  *glb_b = MkLong(0u);
  Long swp_tim = MkLong(0);
  stopwatch_reset(swp_tim);

  const unsigned swp = swp_max[1u];
  unsigned blk_swp = 0u;

#ifdef ANIMATE
#if (ANIMATE == 1)
  vn_mtxvis_ctx *ctx = static_cast<vn_mtxvis_ctx*>(NULL);
  if (ncol < 10000u) {
    char fname[8] = { '\0' };
    (void)sprintf(fname, "FG%x%04u", routine, ncol);
    SYSI_CALL(vn_mtxvis_start(&ctx, fname, (VN_MTXVIS_OP_AtA | VN_MTXVIS_FN_Lg | VN_MTXVIS_FF_Bin), nrow, ncol, 1, 1, 7));
    if (ctx) {
      SYSI_CALL(vn_mtxvis_frame(ctx, hF, ldhF));
      SYSI_CALL(vn_mtxvis_frame(ctx, hG, ldhG));
    }
  }
#elif (ANIMATE == 2)
  vn_mtxvis_ctx *ctxF = static_cast<vn_mtxvis_ctx*>(NULL);
  vn_mtxvis_ctx *ctxG = static_cast<vn_mtxvis_ctx*>(NULL);
  if (ncol < 10000u) {
    char fname[8] = { '\0' };
    (void)sprintf(fname, "%c%x_%04u", 'F', routine, ncol);
    SYSI_CALL(vn_mtxvis_start(&ctxF, fname, (VN_MTXVIS_OP_AtA | VN_MTXVIS_FN_Lg | VN_MTXVIS_FF_Bin), nrow, ncol, 1, 1, 7));
    if (ctxF)
      SYSI_CALL(vn_mtxvis_frame(ctxF, hF, ldhF));
    (void)sprintf(fname, "%c%x_%04u", 'G', routine, ncol);
    SYSI_CALL(vn_mtxvis_start(&ctxG, fname, (VN_MTXVIS_OP_AtA | VN_MTXVIS_FN_Lg | VN_MTXVIS_FF_Bin), nrow, ncol, 1, 1, 7));
    if (ctxG)
      SYSI_CALL(vn_mtxvis_frame(ctxG, hG, ldhG));
  }
#endif // ?ANIMATE
#endif // ANIMATE

  while (blk_swp < swp) {
    *cvg_dat = 0ul;

    for (unsigned blk_stp = 0u; blk_stp < STRAT1_STEPS; ++blk_stp) {
      HZ_L1(blk_stp);

#ifdef ANIMATE
#if (ANIMATE == 1)
      if (ctx) {
        CUDA_CALL(cudaDeviceSynchronize());
        CUDA_CALL(cudaMemcpy2DAsync(hF, ldhF * sizeof(double), dF, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaMemcpy2DAsync(hG, ldhG * sizeof(double), dG, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaDeviceSynchronize());
        SYSI_CALL(vn_mtxvis_frame(ctx, hF, ldhF));
        SYSI_CALL(vn_mtxvis_frame(ctx, hG, ldhG));
      }
#elif (ANIMATE == 2)
      if (ctxF) {
        CUDA_CALL(cudaDeviceSynchronize());
        CUDA_CALL(cudaMemcpy2DAsync(hF, ldhF * sizeof(double), dF, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaDeviceSynchronize());
        SYSI_CALL(vn_mtxvis_frame(ctxF, hF, ldhF));
      }
      if (ctxG) {
        CUDA_CALL(cudaDeviceSynchronize());
        CUDA_CALL(cudaMemcpy2DAsync(hG, ldhG * sizeof(double), dG, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaDeviceSynchronize());
        SYSI_CALL(vn_mtxvis_frame(ctxG, hG, ldhG));
      }
#endif // ?ANIMATE
#endif // ANIMATE
    }

    ++blk_swp;
    CUDA_CALL(cudaDeviceSynchronize());

    const unsigned Long cvg = *cvg_dat;
    const unsigned cvg_s = static_cast<unsigned>(cvg);
    *glb_s += cvg_s;
    const unsigned cvg_b = static_cast<unsigned>(cvg >> 32u);
    *glb_b += cvg_b;
    const double tim_s = stopwatch_lap(swp_tim) * TS2S;
    (void)printf("BLK_SWP(%2u), ROT_S(%10u), ROT_B(%10u), TIME(%#12.6f s)\n", blk_swp, cvg_s, cvg_b, tim_s);
    if (!cvg_b)
      break;

    initS(0, ncol, static_cast<cudaStream_t>(NULL));
    CUDA_CALL(cudaDeviceSynchronize());
  }

#ifdef ANIMATE
#if (ANIMATE == 1)
  if (ctx)
    SYSI_CALL(vn_mtxvis_stop(ctx));  
#elif (ANIMATE == 2)
  if (ctxG)
    SYSI_CALL(vn_mtxvis_stop(ctxG));
  if (ctxF)
    SYSI_CALL(vn_mtxvis_stop(ctxF));
#endif // ?ANIMATE
#endif // ANIMATE

  *glbSwp = blk_swp;
  initS(1, ncol, static_cast<cudaStream_t>(NULL));
  CUDA_CALL(cudaDeviceSynchronize());

  timers[2] = stopwatch_lap(timers[3]);

  CUDA_CALL(cudaMemcpy2DAsync(hF, ldhF * sizeof(double), dF, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy2DAsync(hG, ldhG * sizeof(double), dG, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy2DAsync(hV, ldhV * sizeof(double), dV, ldd * sizeof(double), nrow * sizeof(double), ncol, cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpyAsync(hS, dS, ncol * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpyAsync(hH, dH, ncol * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpyAsync(hK, dK, ncol * sizeof(double), cudaMemcpyDeviceToHost));

  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaFreeHost((void*)cvg_dat));
  CUDA_CALL(cudaFree(static_cast<void*>(ptr)));

  timers[3] = stopwatch_lap(timers[3]);
  timers[0] = stopwatch_lap(timers[0]);

  if (timing)
    for (unsigned i = 0u; i < 4u; ++i)
      timing[i] = timers[i] * TS2S;

  return 0;
}
