#ifndef DEVICE_CODE_COMMON_FERMI_HPP
#define DEVICE_CODE_COMMON_FERMI_HPP

// sum x
// warp-synchornous shmem reduction
MYDEVFN double
dSum32(const double x, volatile double *const shPtr, const unsigned lid)
{
  shPtr[lid] = x;
  if (lid < 16u)
    shPtr[lid] += shPtr[lid + 16u];
  if (lid < 8u)
    shPtr[lid] += shPtr[lid + 8u];
  if (lid < 4u)
    shPtr[lid] += shPtr[lid + 4u];
  if (lid < 2u)
    shPtr[lid] += shPtr[lid + 2u];
  if (lid < 1u)
    shPtr[lid] += shPtr[lid + 1u];
  return *shPtr;
}

// max|x|
// warp-synchornous shmem reduction
MYDEVFN double
dMax32(const double x, volatile double *const shPtr, const unsigned lid)
{
  shPtr[lid] = fabs(x);
  if (lid < 16u)
    shPtr[lid] = fmax(shPtr[lid], shPtr[lid + 16u]);
  if (lid < 8u)
    shPtr[lid] = fmax(shPtr[lid], shPtr[lid + 8u]);
  if (lid < 4u)
    shPtr[lid] = fmax(shPtr[lid], shPtr[lid + 4u]);
  if (lid < 2u)
    shPtr[lid] = fmax(shPtr[lid], shPtr[lid + 2u]);
  if (lid < 1u)
    shPtr[lid] = fmax(shPtr[lid], shPtr[lid + 1u]);
  return *shPtr;
}

// min|x|, x =/= 0
// warp-synchornous shmem reduction
MYDEVFN double
dMin32(const double x, volatile double *const shPtr, const unsigned lid)
{
  shPtr[lid] = ((x == 0.0) ? DBL_MAX : fabs(x));
  if (lid < 16u)
    shPtr[lid] = fmin(shPtr[lid], shPtr[lid + 16u]);
  if (lid < 8u)
    shPtr[lid] = fmin(shPtr[lid], shPtr[lid + 8u]);
  if (lid < 4u)
    shPtr[lid] = fmin(shPtr[lid], shPtr[lid + 4u]);
  if (lid < 2u)
    shPtr[lid] = fmin(shPtr[lid], shPtr[lid + 2u]);
  if (lid < 1u)
    shPtr[lid] = fmin(shPtr[lid], shPtr[lid + 1u]);
  return *shPtr;
}

MYDEVFN double
dSsq32(const double x, volatile double *const shPtr, const unsigned lid)
{
  return dSum32(x * x, shPtr, lid);
}

MYDEVFN double
dSsq32_2(const double x, volatile double *const shPtr, const unsigned lid, const int scl)
{
  const double x_ = fabs(x);
  const double x__ = ((x_ > NU_) ? scalbn(x_, scl) : 0.0);
  return dSsq32(x__, shPtr, lid);
}

MYDEVFN double
dSsq32_1(const double x, volatile double *const shPtr, const unsigned lid)
{
  const double x_ = fabs(x);
  const double x__ = (((x_ >= MU_) && (x_ <= NU_)) ? x_ : 0.0);
  return dSsq32(x__, shPtr, lid);
}

MYDEVFN double
dSsq32_0(const double x, volatile double *const shPtr, const unsigned lid, const int scl)
{
  const double x_ = fabs(x);
  const double x__ = ((x_ < MU_) ? scalbn(x_, scl) : 0.0);
  return dSsq32(x__, shPtr, lid);
}

MYDEVFN int
dSSQ32(double &ssq, const double x, volatile double *const shPtr, const unsigned lid)
{
  const double M = dMax32(x, shPtr, lid);
  if (M == 0.0) {
    ssq = 0.0;
    return INT_MIN;
  }
  const double m = dMin32(x, shPtr, lid);
  const bool lo = (m < MU_), hi = (M > NU_);
  int sc0, sc1, sc2, scl;
  double sq0, sq1, sq2;
  if (hi) {
    sc2 = find_scl(M, NU_);
    sq2 = dSsq32_2(x, shPtr, lid, sc2);
    sc2 *= -2;
    if (m > NU_) {
      scl = sc2;
      ssq = sq2;
    }
    else {
      sc1 = 0;
      sq1 = dSsq32_1(x, shPtr, lid);
      if (sq1 > 0.0) {
        normalize2(sc1, sq1);
        normalize2(sc2, sq2);
        sum_scl(sc1, sc2, sq1, sq2, scl, ssq);
      }
      else {
        scl = sc2;
        ssq = sq2;
      }
    }
  }
  else if (lo) {
    sc0 = find_scl(m, MU_);
    sq0 = dSsq32_0(x, shPtr, lid, sc0);
    sc0 *= -2;
    if (M < MU_) {
      scl = sc0;
      ssq = sq0;
    }
    else {
      sc1 = 0;
      sq1 = dSsq32_1(x, shPtr, lid);
      if (sq1 > 0.0) {
        normalize2(sc1, sq1);
        normalize2(sc0, sq0);
        sum_scl(sc0, sc1, sq0, sq1, scl, ssq);
      }
      else {
        scl = sc0;
        ssq = sq0;
      }
    }
  }
  else {
    scl = 0;
    ssq = dSsq32(x, shPtr, lid);
  }
  return scl;
}

MYDEVFN double
dSSQ_32(const double x, volatile double *const shPtr, const unsigned lid)
{
  double ssq;
  const int scl = dSSQ32(ssq, x, shPtr, lid);
  return scalbn(ssq, scl);
}

MYDEVFN double
dNRM2_32(const double x, volatile double *const shPtr, const unsigned lid)
{
  double ssq;
  const int scl = dSSQ32(ssq, x, shPtr, lid) / 2;
  return scalbn(__dsqrt_rn(ssq), scl);
}

// Ap >= Aq
MYDEVFN double
dCOSA_(const double Ap, const double Aq, const double x, const double y, volatile double *const shPtr, const unsigned lid)
{
  double z, w;
  bool big, small;
  if (Aq >= 1.0) {
    big = ((Ap * Aq) >= DBL_MAX); // Ap*Aq may overflow, but it's safe
    small = false;
  }
  else {
    big = false;
    small = ((Ap * Aq) <= MU_EPS); // Ap*Aq may underflow, but it's safe
  }
  if (big) { // possible overflow
    z = (fabs(x) >= fabs(y)) ? (__ddiv_rn(x, Ap) * y) : (__ddiv_rn(y, Ap) * x);
    w = __ddiv_rn(dSum32(z, shPtr, lid), Aq);
  }
  else if (small) { // possible underflow
    z = (fabs(x) >= fabs(y)) ? (x * __ddiv_rn(y, Aq)) : (y * __ddiv_rn(x, Aq));
    w = __ddiv_rn(dSum32(z, shPtr, lid), Ap);
  }
  else { // no overflow nor underflow
    z = x * y;
    w = __ddiv_rn(__ddiv_rn(dSum32(z, shPtr, lid), Aq), Ap);
  }
  return w;
}

MYDEVFN double
dCOSA(const double Ap, const double Aq, const double x, const double y, volatile double *const shPtr, const unsigned lid)
{
  return (Ap >= Aq) ?
    dCOSA_(Ap, Aq, x, y, shPtr, lid) :
    dCOSA_(Aq, Ap, y, x, shPtr, lid);
}

MYDEVFN void
dQR32
(
 volatile double *const A,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  #pragma unroll
  for (unsigned k = 0u; k < 31u; ++k) {
    double nrm2 = 0.0, Axk, Axk_, beta, alpha_beta, tau;
    // compute the Householer reflector, see DLARFG
    const unsigned my0 = k + y0;
    if (my0 < 32u) {
      Axk = F32(A, x, k);
      Axk_ = ((x >= k) ? Axk : 0.0);
    }
    __syncthreads();
    if (my0 < 32u) {
      const double Axy0 = F32(A, x, y0);
      nrm2 = dNRM2_32(Axk_, &(F32(A, 0u, y0)), x);
      F32(A, x, y0) = Axy0;
    }
    __syncthreads();
    if (nrm2 > 0.0) {
      const double alpha = F32(A, k, k);
      beta = copysign(nrm2, alpha);
      alpha_beta = alpha + beta;
      tau = __ddiv_rn(alpha_beta, beta);
    }
    __syncthreads();
    if (nrm2 > 0.0) {
      // apply the Householder reflector
      if (my0 == k)
        F32(A, x, k) = ((x == k) ? -beta : ((x > k) ? 0.0 : Axk));
      else {
        const double Axy = F32(A, x, my0);
        const double Vxk = ((x == k) ? 1.0 : ((x > k) ? __ddiv_rn(Axk_, alpha_beta) : 0.0));
        const double dp = dSum32(Vxk * Axy, &(F32(A, 0u, my0)), x);
        const double _tdp = -(tau * dp);
        F32(A, x, my0) = __fma_rn(_tdp, Vxk, Axy);
      }
      const unsigned my1 = k + y1;
      if (my1 < 32u) {
        const double Axy = F32(A, x, my1);
        const double Vxk = ((x == k) ? 1.0 : ((x > k) ? __ddiv_rn(Axk_, alpha_beta) : 0.0));
        const double dp = dSum32(Vxk * Axy, &(F32(A, 0u, my1)), x);
        const double _tdp = -(tau * dp);
        F32(A, x, my1) = __fma_rn(_tdp, Vxk, Axy);
      }
    }
    __syncthreads();
  }
}

#endif // !DEVICE_CODE_COMMON_FERMI_HPP
