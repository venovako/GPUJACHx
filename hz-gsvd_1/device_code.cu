#include "HZ.hpp"

#include "device_code.hpp"
#include "device_code_common.hpp"
#include "device_code_accumV.hpp"
#if (CVG == 0)
#include "device_code_cdsort_0.hpp"
#elif (CVG == 1)
#include "device_code_cdsort_1.hpp"
#else // unknown CVG
#error CVG unknown
#endif // ?CVG
#include "device_code_cdsort_accumV.hpp"

#include "my_utils.hpp"

static const dim3 hzL1bD(HZ_L1_THREADS_PER_BLOCK_X, HZ_L1_THREADS_PER_BLOCK_Y, 1u);

void HZ_L1_sv(const unsigned step) throw()
{
  const dim3 hzL1gD(STRAT1_PAIRS, 1u, 1u);
  CUDA_CALL(cudaConfigureCall(hzL1gD, hzL1bD));
  CUDA_CALL(cudaSetupArgument(step, static_cast<size_t>(0u)));
  CUDA_CALL(cudaLaunch(dHZ_L1_sv));
}

void initS(const int full, const unsigned nRank, const cudaStream_t s) throw()
{
  const dim3 bD(2u * WARP_SZ, 1u, 1u);
  const dim3 gD(udiv_ceil(nRank * WARP_SZ, bD.x), 1u, 1u);
  const size_t shmD = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(gD, bD, shmD, s));
  CUDA_CALL(cudaSetupArgument(full, static_cast<size_t>(0u)));
  CUDA_CALL(cudaLaunch(dInitS));
}

void initV(const unsigned nRank, const cudaStream_t s) throw()
{
  const dim3 bD(2u * WARP_SZ, 1u, 1u);
  const dim3 gD(udiv_ceil(nRank * WARP_SZ, bD.x), 1u, 1u);
  const size_t shmD = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(gD, bD, shmD, s));
  CUDA_CALL(cudaLaunch(dInitV));
}

void initSymbols
(
 double *const F,
 double *const G,
 double *const V,
 double *const S,
 double *const H,
 double *const K,
 volatile unsigned Long *const cvg,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned ldF,
 const unsigned ldG,
 const unsigned ldV,
 const unsigned nSwp,
 const double alpha,
 const double beta,
 const double* &alpha_ptr,
 const double* &beta_ptr
) throw()
{
  CUDA_CALL(cudaMemcpyToSymbolAsync(_F, &F, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_G, &G, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_V, &V, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_S, &S, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_H, &H, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_K, &K, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_cvg, &cvg, sizeof(volatile unsigned Long*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nRow, &nRow, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nRank, &nRank, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_ldF, &ldF, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_ldG, &ldG, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_ldV, &ldV, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nSwp, &nSwp, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_STRAT0_STEPS, &STRAT0_STEPS, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_STRAT0_PAIRS, &STRAT0_PAIRS, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_alpha, &alpha, sizeof(double)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_beta, &beta, sizeof(double)));
  // copy strategy tables
  CUDA_CALL(cudaMemcpyToSymbolAsync(_strat0, strat0, sizeof(strat0)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_strat1, strat1, sizeof(strat1)));
  // get symbol addresses
  CUDA_CALL(cudaGetSymbolAddress((void**)&alpha_ptr, _alpha));
  CUDA_CALL(cudaGetSymbolAddress((void**)&beta_ptr, _beta));
}
