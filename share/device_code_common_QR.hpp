#ifndef DEVICE_CODE_COMMON_QR_HPP
#define DEVICE_CODE_COMMON_QR_HPP

MYDEVFN void dPeelOff
(
 volatile double *const R0,
 volatile double *const R1,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  #pragma unroll
  for (unsigned k = 0u; k < 32u; ++k) {
    const unsigned
      my0 = (k + y0),
      my1 = (k + y1);
    unsigned x_k;
    double c, s, r, a_, b_;
    bool store = false;
    const bool active = ((x >= k) && (x <= my1));
    if ((my0 < 32u) && active) {
      x_k = x - k;
      const double f = F32(R0, x, x);
      const double g = F32(R1, x_k, x);
      dGivens(f, g, c, s, r);
      if (x < my0) {
        const double a = F32(R0, x, my0);
        const double b = F32(R1, x_k, my0);
        a_ = c * a + s * b;
        b_ = c * b - s * a;
        store = true;
      }
      else if (x == my0) {
        a_ = r;
        b_ = 0.0;
        store = true;
      }
    }
    __syncthreads();
    if (store) {
      F32(R0, x, my0) = a_;
      F32(R1, x_k, my0) = b_;
      store = false;
    }
    __syncthreads();
    if ((my1 < 32u) && active) {
      if (x < my1) {
        const double a = F32(R0, x, my1);
        const double b = F32(R1, x_k, my1);
        a_ = c * a + s * b;
        b_ = c * b - s * a;
        store = true;
      }
      else if (x == my1) {
        a_ = r;
        b_ = 0.0;
        store = true;
      }
    }
    __syncthreads();
    if (store) {
      F32(R0, x, my1) = a_;
      F32(R1, x_k, my1) = b_;
    }
    __syncthreads();
  }
}

MYDEVFN void dFactorize
(
 const double *const G0,
 const double *const G1,
 volatile double *const A,
 volatile double *const V,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  F32(A, x, y0) = G0[x];
  F32(A, x, y1) = G1[x];
  __syncthreads();
  dQR32(A, x, y0, y1);

  for (unsigned i = x + 32u; i < _nRow; i += 32u) {
    F32(V, x, y0) = G0[i];
    F32(V, x, y1) = G1[i];
    __syncthreads();
    dQR32(V, x, y0, y1);
    dPeelOff(A, V, x, y0, y1);
  }

  // second 32 columns set to identity (initial V)
  F32(V, x, y0) = ((y0 == x) ? +1.0 : +0.0);
  F32(V, x, y1) = ((y1 == x) ? +1.0 : +0.0);
  __syncthreads();
}

#endif // !DEVICE_CODE_COMMON_QR_HPP
