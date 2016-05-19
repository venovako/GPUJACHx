#ifndef DEVICE_CODE_HPP
#define DEVICE_CODE_HPP

#include "cuda_helper.hpp"

extern void defJacL1(const unsigned step, const int definite /* unused */) throw();
extern void hypJacL1(const unsigned step, const unsigned npos) throw();

extern void defJacL1v(const unsigned step, const int definite /* unused */) throw();
extern void hypJacL1v(const unsigned step, const unsigned npos) throw();

extern void defJacL1s(const unsigned step, const int definite) throw();
extern void hypJacL1s(const unsigned step, const unsigned npos) throw();

extern void defJacL1sv(const unsigned step, const int definite) throw();
extern void hypJacL1sv(const unsigned step, const unsigned npos) throw();

extern void initD(double *const G, double *const D, const unsigned ifc, const unsigned nRow, const unsigned nRank, const unsigned nPos, const unsigned ldG, const cudaStream_t s) throw();
extern void initV(double *const V, const unsigned ifc, const unsigned nRank, const unsigned ldV, const cudaStream_t s) throw();

extern void initSymbols(double *const G, double *const V, volatile unsigned Long *const cvg, const unsigned nRow, const unsigned ldG, const unsigned ldV, const unsigned nSwp, const double alpha, const double beta, const double* &alpha_ptr, const double* &beta_ptr) throw();

#endif // !DEVICE_CODE_HPP
