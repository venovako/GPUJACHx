#ifndef DEVICE_CODE_COMMON_KEPLER_HPP
#define DEVICE_CODE_COMMON_KEPLER_HPP

#ifndef _shfl_xor
#if __CUDACC_VER_MAJOR__ >= 9
#define _shfl_xor(x,y) __shfl_xor_sync(~0u, (x), (y))
#else /* __CUDACC_VER_MAJOR__ < 9 */
#define _shfl_xor(x,y) __shfl_xor((x), (y))
#endif /* ?__CUDACC_VER_MAJOR__ */
#else /* _shfl_xor */
#error _shfl_xor already defined
#endif /* ?_shfl_xor */

#ifndef _shfl
#if __CUDACC_VER_MAJOR__ >= 9
#define _shfl(x,y) __shfl_sync(~0u, (x), (y))
#else /* __CUDACC_VER_MAJOR__ < 9 */
#define _shfl(x,y) __shfl((x), (y))
#endif /* ?__CUDACC_VER_MAJOR__ */
#else /* _shfl */
#error _shfl already defined
#endif /* ?_shfl */

// sum x
// Kepler warp shuffle
MYDEVFN double
dSum32(const double x)
{
  int lo_my, hi_my, lo_his, hi_his;
  double x_my = x, x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 16);
  hi_his = _shfl_xor(hi_my, 16);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my += x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 8);
  hi_his = _shfl_xor(hi_my, 8);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my += x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 4);
  hi_his = _shfl_xor(hi_my, 4);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my += x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 2);
  hi_his = _shfl_xor(hi_my, 2);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my += x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 1);
  hi_his = _shfl_xor(hi_my, 1);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my += x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl(lo_my, 0);
  hi_his = _shfl(hi_my, 0);
  x_his = __hiloint2double(hi_his, lo_his);

  return x_his;
}

// max|x|
// Kepler warp shuffle
MYDEVFN double
dMax32(const double x)
{
  int lo_my, hi_my, lo_his, hi_his;
  double x_my = fabs(x), x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 16);
  hi_his = _shfl_xor(hi_my, 16);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmax(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 8);
  hi_his = _shfl_xor(hi_my, 8);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmax(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 4);
  hi_his = _shfl_xor(hi_my, 4);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmax(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 2);
  hi_his = _shfl_xor(hi_my, 2);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmax(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 1);
  hi_his = _shfl_xor(hi_my, 1);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmax(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl(lo_my, 0);
  hi_his = _shfl(hi_my, 0);
  x_his = __hiloint2double(hi_his, lo_his);

  return x_his;
}

// min|x|, x =/= 0
// Kepler warp shuffle
MYDEVFN double
dMin32(const double x)
{
  int lo_my, hi_my, lo_his, hi_his;
  double x_my = ((x == 0.0) ? DBL_MAX : fabs(x)), x_his;

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 16);
  hi_his = _shfl_xor(hi_my, 16);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmin(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 8);
  hi_his = _shfl_xor(hi_my, 8);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmin(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 4);
  hi_his = _shfl_xor(hi_my, 4);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmin(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 2);
  hi_his = _shfl_xor(hi_my, 2);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmin(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl_xor(lo_my, 1);
  hi_his = _shfl_xor(hi_my, 1);
  x_his = __hiloint2double(hi_his, lo_his);
  x_my = fmin(x_my, x_his);

  lo_my = __double2loint(x_my);
  hi_my = __double2hiint(x_my);
  lo_his = _shfl(lo_my, 0);
  hi_his = _shfl(hi_my, 0);
  x_his = __hiloint2double(hi_his, lo_his);

  return x_his;
}

MYDEVFN double
dSsq32(const double x)
{
  return dSum32(x * x);
}

MYDEVFN double
dSsq32_2(const double x, const int scl)
{
  const double x_ = fabs(x);
  const double x__ = ((x_ > NU_) ? scalbn(x_, scl) : 0.0);
  return dSsq32(x__);
}

MYDEVFN double
dSsq32_1(const double x)
{
  const double x_ = fabs(x);
  const double x__ = (((x_ >= MU_) && (x_ <= NU_)) ? x_ : 0.0);
  return dSsq32(x__);
}

MYDEVFN double
dSsq32_0(const double x, const int scl)
{
  const double x_ = fabs(x);
  const double x__ = ((x_ < MU_) ? scalbn(x_, scl) : 0.0);
  return dSsq32(x__);
}

MYDEVFN int
dSSQ32(double &ssq, const double x)
{
  const double M = dMax32(x);
  if (M == 0.0) {
    ssq = 0.0;
    return INT_MIN;
  }
  const double m = dMin32(x);
  const bool lo = (m < MU_), hi = (M > NU_);
  int sc0, sc1, sc2, scl;
  double sq0, sq1, sq2;
  if (hi) {
    sc2 = find_scl(M, NU_);
    sq2 = dSsq32_2(x, sc2);
    sc2 *= -2;
    if (m > NU_) {
      scl = sc2;
      ssq = sq2;
    }
    else {
      sc1 = 0;
      sq1 = dSsq32_1(x);
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
    sq0 = dSsq32_0(x, sc0);
    sc0 *= -2;
    if (M < MU_) {
      scl = sc0;
      ssq = sq0;
    }
    else {
      sc1 = 0;
      sq1 = dSsq32_1(x);
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
    ssq = dSsq32(x);
  }
  return scl;
}

MYDEVFN double
dSSQ_32(const double x)
{
  double ssq;
  const int scl = dSSQ32(ssq, x);
  return scalbn(ssq, scl);
}

MYDEVFN double
dNRM2_32(const double x)
{
  double ssq;
  const int scl = dSSQ32(ssq, x) / 2;
  return scalbn(__dsqrt_rn(ssq), scl);
}

// Ap >= Aq
MYDEVFN double
dCOSA_(const double Ap, const double Aq, const double x, const double y)
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
    w = __ddiv_rn(dSum32(z), Aq);
  }
  else if (small) { // possible underflow
    z = (fabs(x) >= fabs(y)) ? (x * __ddiv_rn(y, Aq)) : (y * __ddiv_rn(x, Aq));
    w = __ddiv_rn(dSum32(z), Ap);
  }
  else { // no overflow nor underflow
    z = x * y;
    w = __ddiv_rn(__ddiv_rn(dSum32(z), Aq), Ap);
  }
  return w;
}

MYDEVFN double
dCOSA(const double Ap, const double Aq, const double x, const double y)
{
  return (Ap >= Aq) ?
    dCOSA_(Ap, Aq, x, y) :
    dCOSA_(Aq, Ap, y, x);
}

MYDEVFN void dQR32
(volatile double *const A,
 const unsigned x,
 const unsigned y0,
 const unsigned y1)
{
  #pragma unroll
  for (unsigned k = 0u; k < 31u; ++k) {
    double nrm2 = 0.0, Axk, Axk_, beta, alpha_beta, tau;
    // compute the Householer reflector, see DLARFG
    const unsigned my0 = k + y0;
    if (my0 < 32u) {
      Axk = F32(A, x, k);
      Axk_ = ((x >= k) ? Axk : 0.0);
      nrm2 = dNRM2_32(Axk_);
      if (nrm2 > 0.0) {
        const double alpha = F32(A, k, k);
        beta = copysign(nrm2, alpha);
        alpha_beta = alpha + beta;
        tau = __ddiv_rn(alpha_beta, beta);
      }
    }
    __syncthreads();
    if (nrm2 > 0.0) {
      // apply the Householder reflector
      if (my0 == k)
        F32(A, x, k) = ((x == k) ? -beta : ((x > k) ? 0.0 : Axk));
      else {
        const double Axy = F32(A, x, my0);
        const double Vxk = ((x == k) ? 1.0 : ((x > k) ? __ddiv_rn(Axk_, alpha_beta) : 0.0));
        const double dp = dSum32(Vxk * Axy);
        const double _tdp = -(tau * dp);
        F32(A, x, my0) = __fma_rn(_tdp, Vxk, Axy);
      }
      const unsigned my1 = k + y1;
      if (my1 < 32u) {
        const double Axy = F32(A, x, my1);
        const double Vxk = ((x == k) ? 1.0 : ((x > k) ? __ddiv_rn(Axk_, alpha_beta) : 0.0));
        const double dp = dSum32(Vxk * Axy);
        const double _tdp = -(tau * dp);
        F32(A, x, my1) = __fma_rn(_tdp, Vxk, Axy);
      }
    }
    __syncthreads();
  }
}

#endif /* !DEVICE_CODE_COMMON_KEPLER_HPP */
