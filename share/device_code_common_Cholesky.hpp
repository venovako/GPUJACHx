#ifndef DEVICE_CODE_COMMON_CHOLESKY_HPP
#define DEVICE_CODE_COMMON_CHOLESKY_HPP

MYDEVFN void dGtG
(const double *const G0,
 const double *const G1,
 volatile double *const A,
 volatile double *const V,
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
    F64(A, x, y0) = G0[i];
    F64(A, x, y1) = G1[i];

    i += 32u;

    F64(A, x32, y0) = G0[i];
    F64(A, x32, y1) = G1[i];

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

  // first 32 columns set (G^T G unsymmetrized)
  F32(A, x, y0) = y0x;
  F32(A, x, y1) = y1x;

  // second 32 columns set to identity (initial V)
  y0x = ((y0 == x) ? +1.0 : +0.0);
  y1x = ((y1 == x) ? +1.0 : +0.0);

  F32(V, x, y0) = y0x;
  F32(V, x, y1) = y1x;

  __syncthreads();
}

MYDEVFN void dCholesky32
(volatile double *const A,
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
#ifdef USE_RSQRT
    const double Akk = (((y0 == k) && (x >= k)) ? F32(A, k, k) : 0.0);
#else // DIV
    const double Akk = (((y0 == k) && (x >= k)) ? __dsqrt_rn(F32(A, k, k)) : 0.0);
#endif // ?USE_RSQRT
    __syncthreads();
    if ((y0 == k) && (x >= k)) {
#ifdef USE_RSQRT
      F32(A, x, k) = (x > k) ? (F32(A, x, k) * my_drsqrt_rn(Akk)) : __dsqrt_rn(Akk);
#else // DIV
      F32(A, x, k) = (x > k) ? __ddiv_rn(F32(A, x, k), Akk) : Akk;
#endif // ?USE_RSQRT
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
#ifdef USE_RSQRT
    const double Akk = (((y1 == k) && (x >= k)) ? F32(A, k, k) : 0.0);
#else // DIV
    const double Akk = (((y1 == k) && (x >= k)) ? __dsqrt_rn(F32(A, k, k)) : 0.0);
#endif // ?USE_RSQRT
    __syncthreads();
    if ((y1 == k) && (x >= k)) {
#ifdef USE_RSQRT
      F32(A, x, k) = (x > k) ? (F32(A, x, k) * my_drsqrt_rn(Akk)) : __dsqrt_rn(Akk);
#else // DIV
      F32(A, x, k) = (x > k) ? __ddiv_rn(F32(A, x, k), Akk) : Akk;
#endif // USE_RSQRT
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
(const double *const G0,
 const double *const G1,
 volatile double *const A,
 volatile double *const V,
 const unsigned x,
 const unsigned y0,
 const unsigned y1)
{
  dGtG(G0, G1, A, V, x, y0, y1);
  dCholesky32(A, x, y0, y1);
}

#endif // !DEVICE_CODE_COMMON_CHOLESKY_HPP
