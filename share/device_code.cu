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
  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));
  CUDA_CALL(cudaSetupArgument(step, static_cast<size_t>(0u)));
  CUDA_CALL(cudaLaunch(dDefJacL1));
}

void hypJacL1(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));

  CUDA_CALL(cudaSetupArgument(step, off));
  off += sizeof(step);
  CUDA_CALL(cudaSetupArgument(npos, off));

  CUDA_CALL(cudaLaunch(dHypJacL1));
}

void defJacL1v(const unsigned step, const int definite VAR_UNUSED) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));
  CUDA_CALL(cudaSetupArgument(step, static_cast<size_t>(0u)));
  CUDA_CALL(cudaLaunch(dDefJacL1v));
}

void hypJacL1v(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));

  CUDA_CALL(cudaSetupArgument(step, off));
  off += sizeof(step);
  CUDA_CALL(cudaSetupArgument(npos, off));

  CUDA_CALL(cudaLaunch(dHypJacL1v));
}

void defJacL1s(const unsigned step, const int definite) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));

  CUDA_CALL(cudaSetupArgument(step, off));
  off += sizeof(step);
  CUDA_CALL(cudaSetupArgument(definite, off));

  CUDA_CALL(cudaLaunch(dDefJacL1s));
}

void hypJacL1s(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));

  CUDA_CALL(cudaSetupArgument(step, off));
  off += sizeof(step);
  CUDA_CALL(cudaSetupArgument(npos, off));

  CUDA_CALL(cudaLaunch(dHypJacL1s));
}

void defJacL1sv(const unsigned step, const int definite) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));

  CUDA_CALL(cudaSetupArgument(step, off));
  off += sizeof(step);
  CUDA_CALL(cudaSetupArgument(definite, off));

  CUDA_CALL(cudaLaunch(dDefJacL1sv));
}

void hypJacL1sv(const unsigned step, const unsigned npos) throw()
{
  const dim3 jacL1gD(STRAT1_PAIRS, 1u, 1u);
  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(jacL1gD, jacL1bD));

  CUDA_CALL(cudaSetupArgument(step, off));
  off += sizeof(step);
  CUDA_CALL(cudaSetupArgument(npos, off));

  CUDA_CALL(cudaLaunch(dHypJacL1sv));
}

void initD
(
 double *const G,
 double *const D,
 const unsigned ifc,
 const unsigned nRow,
 const unsigned nRank,
 const unsigned nPlus,
 const unsigned ldG,
 const cudaStream_t s
) throw()
{
  const dim3 bD(2u * WARP_SZ, 1u, 1u);
  const dim3 gD(udiv_ceil(nRank * WARP_SZ, bD.x), 1u, 1u);
  const size_t shmD =
#if __CUDA_ARCH__ >= 300
    static_cast<size_t>(0u)
#else // Fermi
    bD.x * sizeof(double)
#endif // __CUDA_ARCH__
  ;

  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(gD, bD, shmD, s));

  CUDA_CALL(cudaSetupArgument(G, off));
  off += sizeof(G);
  CUDA_CALL(cudaSetupArgument(D, off));
  off += sizeof(D);
  CUDA_CALL(cudaSetupArgument(ifc, off));
  off += sizeof(ifc);
  CUDA_CALL(cudaSetupArgument(nRow, off));
  off += sizeof(nRow);
  CUDA_CALL(cudaSetupArgument(nRank, off));
  off += sizeof(nRank);
  CUDA_CALL(cudaSetupArgument(nPlus, off));
  off += sizeof(nPlus);
  CUDA_CALL(cudaSetupArgument(ldG, off));

  CUDA_CALL(cudaLaunch(dInitD));
}

void initV
(
 double *const V,
 const unsigned ifc,
 const unsigned nRank,
 const unsigned ldV,
 const cudaStream_t s
) throw()
{
  const dim3 bD(2u * WARP_SZ, 1u, 1u);
  const dim3 gD(udiv_ceil(nRank, bD.x), 1u, 1u);

  size_t off = static_cast<size_t>(0u);

  CUDA_CALL(cudaConfigureCall(gD, bD, static_cast<size_t>(0u), s));

  CUDA_CALL(cudaSetupArgument(V, off));
  off += sizeof(V);
  CUDA_CALL(cudaSetupArgument(ifc, off));
  off += sizeof(ifc);
  CUDA_CALL(cudaSetupArgument(nRank, off));
  off += sizeof(nRank);
  CUDA_CALL(cudaSetupArgument(ldV, off));

  CUDA_CALL(cudaLaunch(dInitV));
}

void initSymbols
(
 double *const G,
 double *const V,
 volatile unsigned Long *const cvg,
 const unsigned nRow,
 const unsigned ldG,
 const unsigned ldV,
 const unsigned nSwp,
 const double alpha,
 const double beta,
 const double* &alpha_ptr,
 const double* &beta_ptr
) throw()
{
  CUDA_CALL(cudaMemcpyToSymbolAsync(_G, &G, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_V, &V, sizeof(double*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_cvg, &cvg, sizeof(volatile unsigned Long*)));
  CUDA_CALL(cudaMemcpyToSymbolAsync(_nRow, &nRow, sizeof(unsigned)));
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
