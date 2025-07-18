#ifndef DEVICE_CODE_COMMON_HPP
#define DEVICE_CODE_COMMON_HPP

#ifndef MYKERN
#define MYKERN __global__ void
#else /* MYKERN */
#error MYKERN not definable externally
#endif /* ?MYKERN */

#ifndef MYDEVFN
#ifdef NDEBUG
#define MYDEVFN __device__ __forceinline__
#else /* DEBUG */
#define MYDEVFN __device__
#endif /* ?NDEBUG */
#else /* MYDEVFN */
#error MYDEVFN not definable externally
#endif /* ?MYDEVFN */

#ifndef HYPJACL1_MAX_THREADS_PER_BLOCK
#define HYPJACL1_MAX_THREADS_PER_BLOCK 512
#else /* HYPJACL1_MAX_THREADS_PER_BLOCK */
#error HYPJACL1_MAX_THREADS_PER_BLOCK not definable externally
#endif /* ?HYPJACL1_MAX_THREADS_PER_BLOCK */

#ifndef HYPJACL1_THREADS_PER_BLOCK_X
#define HYPJACL1_THREADS_PER_BLOCK_X 32u
#else /* HYPJACL1_THREADS_PER_BLOCK_X */
#error HYPJACL1_THREADS_PER_BLOCK_X not definable externally
#endif /* ?HYPJACL1_THREADS_PER_BLOCK_X */

#ifndef HYPJACL1_THREADS_PER_BLOCK_Y
#define HYPJACL1_THREADS_PER_BLOCK_Y 16u
#else /* HYPJACL1_THREADS_PER_BLOCK_Y */
#error HYPJACL1_THREADS_PER_BLOCK_Y not definable externally
#endif /* ?HYPJACL1_THREADS_PER_BLOCK_Y */

#ifndef HYPJACL1_MIN_BLOCKS_PER_SM
#define HYPJACL1_MIN_BLOCKS_PER_SM 1
#else /* HYPJACL1_MIN_BLOCKS_PER_SM */
#error HYPJACL1_MIN_BLOCKS_PER_SM not definable externally
#endif /* ?HYPJACL1_MIN_BLOCKS_PER_SM */

#if (WARP_SZ != 32u)
#error WARP_SZ not 32
#endif /* ?WARP_SZ */

#ifndef WARP_SZ_LG
#define WARP_SZ_LG 5u
#else /* WARP_SZ_LG */
#error WARP_SZ_LG not definable externally
#endif /* ?WARP_SZ_LG */

#ifndef WARP_SZ_LGi
#define WARP_SZ_LGi 5
#else /* WARP_SZ_LGi */
#error WARP_SZ_LGi not definable externally
#endif /* ?WARP_SZ_LGi */

#ifndef WARP_SZ_SUB1
#define WARP_SZ_SUB1 31u
#else /* WARP_SZ_SUB1 */
#error WARP_SZ_SUB1 not definable externally
#endif /* ?WARP_SZ_SUB1 */

#ifndef F32
#define F32(A, r, c) (A)[(c) * 32u + (r)]
#else /* F32 */
#error F32 not definable externally
#endif /* ?F32 */

#ifndef F64
#define F64(A, r, c) (A)[(c) * 64u + (r)]
#else /* F64 */
#error F64 not definable externally
#endif /* ?F64 */

#ifdef USE_DRMAC
#ifndef USE_QR
#define USE_QR
#endif /* !USE_QR */
#ifndef USE_RSQRT
#define USE_RSQRT
#endif /* !USE_RSQRT */
#endif /* USE_DRMAC */

#include "device_code_globals.hpp"

// Thanks to Norbert Juffa of NVIDIA.
#ifdef USE_RSQRT
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
    // Round as shown in Peter Markstein, "IA-64 and Elementary Functions"
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
    // Round as shown in Peter Markstein, "IA-64 and Elementary Functions"
    y = __fma_rn(__fma_rn(0.375, e, 0.5), e * y, y);
    a = __hiloint2double(__double2hiint(y) + 0x1ff00000,__double2loint(y));
  }
  return a;
}
#endif /* USE_RSQRT */

#include "device_code_common_rotate.hpp"
#include "device_code_common_scale.hpp"
#include "device_code_common_Kepler.hpp"
//Fermi: #include "device_code_common_Fermi.hpp"

#ifdef USE_QR
#include "device_code_common_QR.hpp"
#else /* Cholesky */
#include "device_code_common_Cholesky.hpp"
#endif /* ?USE_QR */

MYDEVFN void dMultAV
(double *const A0,
 double *const A1,
 volatile double *const A,
 volatile const double *const B,
 const unsigned x,
 const unsigned y0,
 const unsigned y1,
 const unsigned m)
{
  // Cannon-like A*B
  for (unsigned i = x; i < m; i += 32u) {
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
    __syncthreads();
  }
}

MYKERN dInitD
(double *const G,
 double *const D,
 const unsigned ifc,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned nPlus,
 const unsigned ldG)
{
  //Fermi: extern __shared__ double shPool[];
  const unsigned wpb = (blockDim.x + WARP_SZ_SUB1) >> WARP_SZ_LG;
  const unsigned wid = threadIdx.x >> WARP_SZ_LG;

  const unsigned cix = blockIdx.x * wpb + wid;
  if (cix < nRank) {
    const unsigned lid = static_cast<unsigned>(threadIdx.x) & WARP_SZ_SUB1;
    double *const bGi = G + (cix * ldG + lid);
    const double *const eGi = G + (cix * ldG + nRow);
    double x = +0.0, y;
    for (double *pGi = bGi; pGi < eGi; pGi += WARP_SZ) {
      y = *pGi;
      x = __fma_rn(y, y, x);
    }
    y = dSum32(x);
    //Fermi: volatile double *const shPtr = (volatile double*)(shPool + wid * WARP_SZ);
    //Fermi: y = dSum32(x, shPtr, lid);
    x = __dsqrt_rn(y);
    for (double *pGi = bGi; pGi < eGi; pGi += WARP_SZ)
      *pGi = __ddiv_rn(*pGi, x);
    if (!lid)
      D[cix] = (((ifc + cix) < nPlus) ? y : -y);
  }
}

MYKERN dInitV
(double *const V,
 const unsigned ifc,
 const unsigned nRank,
 const unsigned ldV)
{
  const unsigned cix = blockIdx.x * blockDim.x + threadIdx.x;
  if (cix < nRank)
    V[cix * ldV + ifc + cix] = +1.0;
}

#endif /* !DEVICE_CODE_COMMON_HPP */
