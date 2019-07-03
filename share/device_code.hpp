#ifndef DEVICE_CODE_HPP
#define DEVICE_CODE_HPP

#include "cuda_helper.hpp"

EXTERN_C void defJacL1(const unsigned step, const int definite /* unused */) throw();
EXTERN_C void hypJacL1(const unsigned step, const unsigned npos) throw();

EXTERN_C void defJacL1v(const unsigned step, const int definite /* unused */) throw();
EXTERN_C void hypJacL1v(const unsigned step, const unsigned npos) throw();

EXTERN_C void defJacL1s(const unsigned step, const int definite) throw();
EXTERN_C void hypJacL1s(const unsigned step, const unsigned npos) throw();

EXTERN_C void defJacL1sv(const unsigned step, const int definite) throw();
EXTERN_C void hypJacL1sv(const unsigned step, const unsigned npos) throw();

EXTERN_C void initD(double *const G, double *const D, const unsigned ifc, const unsigned nRow, const unsigned nRank, const unsigned nPos, const unsigned ldG) throw();
EXTERN_C void initV(double *const V, const unsigned ifc, const unsigned nRank, const unsigned ldV) throw();

EXTERN_C void initSymbols(double *const G, double *const V, volatile unsigned long long *const cvg, const unsigned nRow, const unsigned nCol, const unsigned ldG, const unsigned ldV, const unsigned nSwp) throw();

#endif // !DEVICE_CODE_HPP
