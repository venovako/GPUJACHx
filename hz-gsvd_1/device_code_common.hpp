#ifndef DEVICE_CODE_COMMON_HPP
#define DEVICE_CODE_COMMON_HPP

#ifndef MYKERN
#define MYKERN __global__ void
#else // MYKERN
#error MYKERN not definable externally
#endif // !MYKERN

#ifndef MYDEVFN
#ifdef NDEBUG
#define MYDEVFN __device__ __forceinline__
#else // DEBUG
#define MYDEVFN __device__
#endif // NDEBUG
#else // MYDEVFN
#error MYDEVFN not definable externally
#endif // !MYDEVFN

#ifndef HZ_L1_MAX_THREADS_PER_BLOCK
#define HZ_L1_MAX_THREADS_PER_BLOCK 512
#else // HZ_L1_MAX_THREADS_PER_BLOCK
#error HZ_L1_MAX_THREADS_PER_BLOCK not definable externally
#endif // !HZ_L1_MAX_THREADS_PER_BLOCK

#ifndef HZ_L1_THREADS_PER_BLOCK_X
#define HZ_L1_THREADS_PER_BLOCK_X 32u
#else // HZ_L1_THREADS_PER_BLOCK_X
#error HZ_L1_THREADS_PER_BLOCK_X not definable externally
#endif // !HZ_L1_THREADS_PER_BLOCK_X

#ifndef HZ_L1_THREADS_PER_BLOCK_Y
#define HZ_L1_THREADS_PER_BLOCK_Y 16u
#else // HZ_L1_THREADS_PER_BLOCK_Y
#error HZ_L1_THREADS_PER_BLOCK_Y not definable externally
#endif // !HZ_L1_THREADS_PER_BLOCK_Y

#ifndef HZ_L1_MIN_BLOCKS_PER_SM
#define HZ_L1_MIN_BLOCKS_PER_SM 1
#else // HZ_L1_MIN_BLOCKS_PER_SM
#error HZ_L1_MIN_BLOCKS_PER_SM not definable externally
#endif // !HZ_L1_MIN_BLOCKS_PER_SM

#if (WARP_SZ != 32u)
#error WARP_SZ not 32
#endif // WARP_SZ

#ifndef WARP_SZ_LG
#define WARP_SZ_LG 5u
#else // WARP_SZ_LG
#error WARP_SZ_LG not definable externally
#endif // !WARP_SZ_LG

#ifndef WARP_SZ_LGi
#define WARP_SZ_LGi 5
#else // WARP_SZ_LGi
#error WARP_SZ_LGi not definable externally
#endif // !WARP_SZ_LGi

#ifndef WARP_SZ_SUB1
#define WARP_SZ_SUB1 31u
#else // WARP_SZ_SUB1
#error WARP_SZ_SUB1 not definable externally
#endif // !WARP_SZ_SUB1

#ifndef F32
#define F32(A, r, c) (A)[(c) * 32u + (r)]
#else // F32
#error F32 not definable externally
#endif // !F32

#ifndef F64
#define F64(A, r, c) (A)[(c) * 64u + (r)]
#else // F64
#error F64 not definable externally
#endif // !F64

#include "device_code_globals.hpp"

// Thanks to Norbert Juffa of NVIDIA.
MYDEVFN double
my_drsqrt_rn(double a)
{
  double y, h, l, e;
  unsigned int ilo, ihi, g, f;
  int d;

  ihi = __double2hiint(a);
  ilo = __double2loint(a);
  if (((unsigned int)ihi) - 0x00100000U < 0x7fe00000U) {
    f = ihi | 0x3fe00000;
    g = f & 0x3fffffff;
    d = g - ihi;
    a = __hiloint2double(g, ilo); 
    y = rsqrt(a);
    h = __dmul_rn(y, y);
    l = __fma_rn(y, y, -h);
    e = __fma_rn(l, -a, __fma_rn(h, -a, 1.0));
    /* Round as shown in Peter Markstein, "IA-64 and Elementary Functions" */
    y = __fma_rn(__fma_rn(0.375, e, 0.5), e * y, y);
    d = d >> 1;
    a = __hiloint2double(__double2hiint(y) + d, __double2loint(y));
  } else if (a == 0.0) {
    a = __hiloint2double((ihi & 0x80000000) | 0x7ff00000, 0x00000000);
  } else if (a < 0.0) {
    a = __hiloint2double(0xfff80000, 0x00000000);
  } else if (isinf(a)) {
    a = __hiloint2double(ihi & 0x80000000, 0x00000000);
  } else if (isnan(a)) {
    a = a + a;
  } else {
    a = a * __hiloint2double(0x7fd00000, 0);
    y = rsqrt(a);
    h = __dmul_rn(y, y);
    l = __fma_rn(y, y, -h);
    e = __fma_rn(l, -a, __fma_rn(h, -a, 1.0));
    /* Round as shown in Peter Markstein, "IA-64 and Elementary Functions" */
    y = __fma_rn(__fma_rn(0.375, e, 0.5), e * y, y);
    a = __hiloint2double(__double2hiint(y) + 0x1ff00000,__double2loint(y));
  }
  return a;
}

#include "device_code_common_rotate.hpp"
#include "device_code_common_Kepler.hpp"
#include "device_code_common_Cholesky.hpp"

MYDEVFN void dMultAV
(
 double *const A0,
 double *const A1,
 volatile double *const A,
 const double *const B,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  // Cannon-like A*B
  for (unsigned i = x; i < _nRow; i += 32u) {
    F32(A, x, y0) = A0[i];
    F32(A, x, y1) = A1[i];

    __syncthreads();

    double
      Cxy0 = +0.0,
      Cxy1 = +0.0;

    // skew (mod 32)
    unsigned
      p0 = ((y0 + x) & 0x1Fu),
      p1 = ((y1 + x) & 0x1Fu);

    // mult-and-cshift (mod 32)
    #pragma unroll
    for (unsigned k = 0u; k < 32u; ++k) {
      Cxy0 = __fma_rn(F32(A, x, p0), F32(B, p0, y0), Cxy0);
      Cxy1 = __fma_rn(F32(A, x, p1), F32(B, p1, y1), Cxy1);
      p0 = (p0 + 1u) & 0x1Fu;
      p1 = (p1 + 1u) & 0x1Fu;
    }
    __syncthreads();

    A0[i] = Cxy0;
    A1[i] = Cxy1;
  }
}

MYDEVFN void dInvNrm2C
(
 const double *const bA,
 const double *const eA,
 double &ssq,
 double &inv_nrm
)
{
  double x = +0.0, y;
  for (const double *pA = bA; pA < eA; pA += WARP_SZ) {
    y = *pA;
    x = __fma_rn(y, y, x);
  }
  ssq = dSum32(x);
  inv_nrm = my_drsqrt_rn(ssq);
}

MYDEVFN void dNrm2InvC
(
 const double *const bA,
 const double *const eA,
 double &ssq,
 double &nrm,
 double &inv_nrm
)
{
  dInvNrm2C(bA, eA, ssq, inv_nrm);
  nrm = __dsqrt_rn(ssq);
}

MYDEVFN void dScalC
(
 double *const bA,
 const double *const eA,
 const double scl
)
{
  for (double *pA = bA; pA < eA; pA += WARP_SZ)
    *pA *= scl;
}

MYDEVFN void dGlobalPostScaleFast
(
 double *const F,
 double *const G,
 double *const V,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned ldF,
 const unsigned ldG,
 const unsigned ldV
)
{
  const unsigned wpb = (blockDim.x + WARP_SZ_SUB1) >> WARP_SZ_LG;
  const unsigned wid = threadIdx.x >> WARP_SZ_LG;

  const unsigned cix = blockIdx.x * wpb + wid;
  if (cix < nRank) {
    unsigned lid;
    asm ("mov.u32 %0, %%laneid;" : "=r"(lid));
    double *const bFi = F + (cix * ldF + lid);
    const double *const eFi = F + (cix * ldF + nRow);
    double *const bGi = G + (cix * ldG + lid);
    const double *const eGi = G + (cix * ldG + nRow);
    double Fi_ssq, Fi_nrm, Fi_inv_nrm;
    dNrm2InvC(bFi, eFi, Fi_ssq, Fi_nrm, Fi_inv_nrm);
    double Gi_ssq, Gi_nrm, Gi_inv_nrm;
    dNrm2InvC(bGi, eGi, Gi_ssq, Gi_nrm, Gi_inv_nrm);
    const double Rhyp = my_drsqrt_rn(Fi_ssq + Gi_ssq);
    if (Rhyp != 1.0) {
      double *const bVi = V + (cix * ldV + lid);
      const double *const eVi = V + (cix * ldV + nRank);
      dScalC(bVi, eVi, Rhyp);
    }
  }
}

MYDEVFN void dGlobalPostScaleFull
(
 double *const F,
 double *const G,
 double *const V,
 double *const S,
 double *const H,
 double *const K,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned ldF,
 const unsigned ldG,
 const unsigned ldV
)
{
  const unsigned wpb = (blockDim.x + WARP_SZ_SUB1) >> WARP_SZ_LG;
  const unsigned wid = threadIdx.x >> WARP_SZ_LG;

  const unsigned cix = blockIdx.x * wpb + wid;
  if (cix < nRank) {
    unsigned lid;
    asm ("mov.u32 %0, %%laneid;" : "=r"(lid));
    double *const bFi = F + (cix * ldF + lid);
    const double *const eFi = F + (cix * ldF + nRow);
    double *const bGi = G + (cix * ldG + lid);
    const double *const eGi = G + (cix * ldG + nRow);
    double Fi_ssq, Fi_nrm, Fi_inv_nrm;
    dNrm2InvC(bFi, eFi, Fi_ssq, Fi_nrm, Fi_inv_nrm);
    double Gi_ssq, Gi_nrm, Gi_inv_nrm;
    dNrm2InvC(bGi, eGi, Gi_ssq, Gi_nrm, Gi_inv_nrm); 
    double Sigmai = Fi_nrm;
    if (Fi_inv_nrm != 1.0)
      dScalC(bFi, eFi, Fi_inv_nrm);
    if (Gi_inv_nrm != 1.0) {
      dScalC(bGi, eGi, Gi_inv_nrm);
      Sigmai *= Gi_inv_nrm;
    }
    double Hi = Fi_nrm;
    double Ki = Gi_nrm;
    const double Rhyp = my_drsqrt_rn(Fi_ssq + Gi_ssq);
    if (Rhyp != 1.0) {
      Hi *= Rhyp;
      Ki *= Rhyp;
      double *const bVi = V + (cix * ldV + lid);
      const double *const eVi = V + (cix * ldV + nRank);
      dScalC(bVi, eVi, Rhyp);
    }
    if (!lid) {
      S[cix] = Sigmai;
      H[cix] = Hi;
      K[cix] = Ki;
    }
  }
}

MYKERN dInitS(const int full)
{
  if (full)
    dGlobalPostScaleFull(_F, _G, _V, _S, _H, _K, _nRow, _nRank, _ldF, _ldG, _ldV);
  else
    dGlobalPostScaleFast(_F, _G, _V, _nRow, _nRank, _ldF, _ldG, _ldV);
}

MYDEVFN void dGlobalInitV
(
 double *const F,
 double *const G,
 double *const V,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned ldF,
 const unsigned ldG,
 const unsigned ldV
)
{
  const unsigned wpb = (blockDim.x + WARP_SZ_SUB1) >> WARP_SZ_LG;
  const unsigned wid = threadIdx.x >> WARP_SZ_LG;

  const unsigned cix = blockIdx.x * wpb + wid;
  if (cix < nRank) {
    unsigned lid;
    asm ("mov.u32 %0, %%laneid;" : "=r"(lid));
    double *const bGi = G + (cix * ldG + lid);
    const double *const eGi = G + (cix * ldG + nRow);
    double Gi_ssq, Gi_inv_nrm;
    dInvNrm2C(bGi, eGi, Gi_ssq, Gi_inv_nrm);
    if (Gi_inv_nrm != 1.0) {
      double *const bFi = F + (cix * ldF + lid);
      const double *const eFi = F + (cix * ldF + nRow);
      dScalC(bFi, eFi, Gi_inv_nrm);
      dScalC(bGi, eGi, Gi_inv_nrm);
    }
    if (!lid)
      V[cix * ldV + cix] = Gi_inv_nrm;
  }
}

MYKERN dInitV()
{
  dGlobalInitV(_F, _G, _V, _nRow, _nRank, _ldF, _ldG, _ldV);
}

#endif // !DEVICE_CODE_COMMON_HPP
