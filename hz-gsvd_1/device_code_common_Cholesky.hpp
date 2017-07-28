#ifndef DEVICE_CODE_COMMON_CHOLESKY_HPP
#define DEVICE_CODE_COMMON_CHOLESKY_HPP

MYDEVFN void
dAtA(const double *const A0,
     const double *const A1,
     volatile double *const A,
     const unsigned x,
     const unsigned y0,
     const unsigned y1)
{
  double
    y0x = +0.0,
    y1x = +0.0;

  const unsigned
    x32 = x + 32u;

  for (unsigned i = x; i < _nRow; i += 32u) {
    F64(A, x, y0) = A0[i];
    F64(A, x, y1) = A1[i];

    i += 32u;

    F64(A, x32, y0) = A0[i];
    F64(A, x32, y1) = A1[i];

    __syncthreads();

    #pragma unroll
    for (unsigned j = 0u; j < 64u; ++j) {
      // x_64 = (x + j) % 64u
      const unsigned x_64 = (x + j) & 0x3Fu;
      const double _x_ = F64(A, x_64, x);
      const double _y0_ = F64(A, x_64, y0);
      const double _y1_ = F64(A, x_64, y1);
      y0x = __fma_rn(_y0_, _x_, y0x);
      y1x = __fma_rn(_y1_, _x_, y1x);
    }
    __syncthreads();
  }

  // first 32 columns set (A^T A unsymmetrized)
  F32(A, x, y0) = y0x;
  F32(A, x, y1) = y1x;

  __syncthreads();
}

MYDEVFN void
dCholesky32(volatile double *const A,
            const unsigned x,
            const unsigned y0,
            const unsigned y1)
{
  //      [ L ? ? ]
  // A -> [ L L ? ]
  //      [ L L L ]

  #pragma unroll
  for (unsigned k = 0u; k < 16u; ++k) {
    // cdiv(k)
    if ((y0 == k) && (x >= k)) {
      const double Akk = F32(A, k, k);
      F32(A, x, k) = (x > k) ? (F32(A, x, k) * my_drsqrt_rn(Akk)) : __dsqrt_rn(Akk);
    }
    __syncthreads();

    unsigned j = (k + 1u) + y0;

    // cmod(j,k)
    if (x >= j) {
      const double Aij = F32(A, x, j);
      const double _Aik = -F32(A, x, k);
      const double Ajk = F32(A, j, k);
      F32(A, x, j) = __fma_rn(_Aik, Ajk, Aij);
    }
    __syncthreads();

    j += 16u;

    // cmod(j+16,k)
    if (x >= j) {
      const double Aij = F32(A, x, j);
      const double _Aik = -F32(A, x, k);
      const double Ajk = F32(A, j, k);
      F32(A, x, j) = __fma_rn(_Aik, Ajk, Aij);
    }
    __syncthreads();
  }

  #pragma unroll
  for (unsigned k = 16u; k < 32u; ++k) {
    // cdiv(k)
    if ((y1 == k) && (x >= k)) {
      const double Akk = __dsqrt_rn(F32(A, k, k));
      F32(A, x, k) = (x > k) ? __ddiv_rn(F32(A, x, k), Akk) : Akk;
    }
    __syncthreads();

    const unsigned j = (k + 1u) + y0;

    // cmod(j,k)
    if (x >= j) {
      const double Aij = F32(A, x, j);
      const double _Aik = -F32(A, x, k);
      const double Ajk = F32(A, j, k);
      F32(A, x, j) = __fma_rn(_Aik, Ajk, Aij);
    }
    __syncthreads();
  }

  //      [ U U U ]
  // A -> [ 0 U U ]
  //      [ 0 0 U ]
  
  double
    Axy0 = +0.0,
    Axy1 = +0.0;

  if (x >= y0)
    Axy0 = F32(A, x, y0);
  if (x >= y1)
    Axy1 = F32(A, x, y1);

  __syncthreads();

  F32(A, y0, x) = Axy0;
  F32(A, y1, x) = Axy1;

  __syncthreads();
}

MYDEVFN void dFactorize
(
 const double *const F0,
 const double *const F1,
 const double *const G0,
 const double *const G1,
 volatile double *const A,
 volatile double *const B,
 volatile double *const V,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  dAtA(F0, F1, A, x, y0, y1);
  dAtA(G0, G1, B, x, y0, y1);
  dCholesky32(A, x, y0, y1);
  dCholesky32(B, x, y0, y1);
}

#endif // !DEVICE_CODE_COMMON_CHOLESKY_HPP
