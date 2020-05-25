#include "HypJac.hpp"

#include "device_code.hpp"
#include "device_code_common.hpp"
#include "device_code_accumV.hpp"
#include "device_code_cdsort.hpp"
#include "device_code_nosort.hpp"
#include "device_code_cdsort_accumV.hpp"
#include "device_code_cdsort_solveV.hpp"
#include "device_code_nosort_accumV.hpp"
#include "device_code_nosort_solveV.hpp"

#include "my_utils.hpp"

static const dim3 jacL1bD(HYPJACL1_THREADS_PER_BLOCK_X, HYPJACL1_THREADS_PER_BLOCK_Y, 1u);

void defJacL1(const unsigned step, const int definite VAR_UNUSED) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dDefJacL1<<< jacL1gD, jacL1bD >>>(step);
}

void hypJacL1(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dHypJacL1<<< jacL1gD, jacL1bD >>>(step, npos);
}

void defJacL1v(const unsigned step, const int definite VAR_UNUSED) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dDefJacL1v<<< jacL1gD, jacL1bD >>>(step);
}

void hypJacL1v(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dHypJacL1v<<< jacL1gD, jacL1bD >>>(step, npos);
}

void defJacL1s(const unsigned step, const int definite) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dDefJacL1s<<< jacL1gD, jacL1bD >>>(step, definite);
}

void hypJacL1s(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dHypJacL1s<<< jacL1gD, jacL1bD >>>(step, npos);
}

void defJacL1sv(const unsigned step, const int definite) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dDefJacL1sv<<< jacL1gD, jacL1bD >>>(step, definite);
}

void hypJacL1sv(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  dHypJacL1sv<<< jacL1gD, jacL1bD >>>(step, npos);
}

void initD
(double *const G,
 double *const D,
 const unsigned ifc,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned nPlus,
 const unsigned ldG
) throw()
{
  const dim3 bD(2u * WARP_SZ, 1u, 1u);
  const dim3 gD(udiv_ceil(nRank * WARP_SZ, bD.x), 1u, 1u);
  const size_t shmD =
#if __CUDA_ARCH__ >= 300
    static_cast<size_t>(0u)
#else /* Fermi */
    bD.x * sizeof(double)
#endif /* ?__CUDA_ARCH__ */
  ;
  dInitD<<< gD, bD, shmD >>>(G, D, ifc, nRow, nRank, nPlus, ldG);
}

void initV
(double *const V,
 const unsigned ifc,
 const unsigned nRank,
 const unsigned ldV
) throw()
{
  const dim3 bD(2u * WARP_SZ, 1u, 1u);
  const dim3 gD(udiv_ceil(nRank, bD.x), 1u, 1u);
  const size_t shmD = static_cast<size_t>(0u);
  dInitV<<< gD, bD, shmD >>>(V, ifc, nRank, ldV);
}

void initSymbols
(double *const G,
 double *const V,
 volatile unsigned long long *const cvg,
 const unsigned nRow,
 const unsigned nCol,
 const unsigned ldG,
 const unsigned ldV,
 const unsigned nSwp
) throw()
{
  CUDA_CALL(cudaMemcpyToSymbolAsync(_G, &G, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_V, &V, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_cvg, &cvg, sizeof(volatile unsigned long long*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nRow, &nRow, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nCol, &nCol, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_ldG, &ldG, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_ldV, &ldV, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nSwp, &nSwp, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_STRAT0_STEPS, &STRAT0_STEPS, sizeof(unsigned)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_STRAT0_PAIRS, &STRAT0_PAIRS, sizeof(unsigned)));
  // copy strategy tables
  CUDA_CALL(cudaMemcpyToSymbolAsync(_strat0, strat0, sizeof(strat0)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_strat1, strat1, sizeof(strat1)));
}
