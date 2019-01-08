#ifndef DEVICE_CODE_COMMON_ROTATE_HPP
#define DEVICE_CODE_COMMON_ROTATE_HPP

MYDEVFN double
tCot2(const double Aq_p, const double Ap_q, const double CosA)
{
  return __ddiv_rn(Aq_p - Ap_q, scalbn(CosA, 1));
}
MYDEVFN double
hCot2(const double Aq_p, const double Ap_q, const double CosA)
{
  return __ddiv_rn(Aq_p + Ap_q, scalbn(CosA, 1));
}

#ifdef USE_RSQRT
MYDEVFN int
dRotT(const double Apq, const double Dp, const double Dq, double &c, double &t)
{
  const double h = Dq - Dp;
  const double ct2 = __ddiv_rn(h, scalbn(Apq, 1)); // h / (2 * Apq)
  const double CT2 = fabs(ct2);
  if (CT2 < SQRT_HEPS) {
    const double cot = CT2 + 1.0;
    t = copysign(__drcp_rn(cot), ct2);
    c = my_drsqrt_rn(__fma_rn(t, t, 1.0));
    return 1;
  }
  else if (CT2 >= SQRT_2_HEPS) {
    const double cot = scalbn(CT2, 1);
    c = 1.0;
    t = copysign(__drcp_rn(cot), ct2);
    return 0;
  }
  else {
    const double cot = CT2 + __dsqrt_rn(__fma_rn(CT2, CT2, 1.0));
    t = copysign(__drcp_rn(cot), ct2);
    c = my_drsqrt_rn(__fma_rn(t, t, 1.0));
    return (c < 1.0);
  }
}
MYDEVFN int
dRotH(const double Apq, const double Dp, const double Dq, double &c, double &t)
{
  const double h = Dq + Dp;
  const double ct2 = -__ddiv_rn(h, scalbn(Apq, 1)); // -h / (2 * Apq)
  const double CT2 = fabs(ct2);
  if (CT2 >= SQRT_2_HEPS) {
    const double cot = scalbn(CT2, 1);
    c = 1.0;
    t = copysign(__drcp_rn(cot), ct2);
    return 0;
  }
  else if (CT2 > 1.0) {
    const double cot = CT2 +  __dsqrt_rn(__fma_rn(CT2, CT2, -1.0));
    t = copysign(__drcp_rn(cot), ct2);
    c = my_drsqrt_rn(__fma_rn(-t, t, 1.0));
    return (c > 1.0);
  }
  else {
    c = HYPJAC_COSHI_FIX;
    t = HYPJAC_TANH1_FIX;
    return 1;
  }
}
#else // COT
MYDEVFN int
dRotT(const double Apq, const double Dp, const double Dq, double &c, double &t)
{
  const double h = Dq - Dp;
  const double ct2 = __ddiv_rn(h, scalbn(Apq, 1)); // h / (2 * Apq)
  const double CT2 = fabs(ct2);
  if (CT2 < SQRT_HEPS) {
    const double cot = CT2 + 1.0;
    c = __ddiv_rn(cot, __dsqrt_rn(__fma_rn(cot, cot, 1.0)));
    t = copysign(__drcp_rn(cot), ct2);
    return 1;
  }
  else if (CT2 >= SQRT_2_HEPS) {
    const double cot = scalbn(CT2, 1);
    c = 1.0;
    t = copysign(__drcp_rn(cot), ct2);
    return 0;
  }
  else {
    const double cot = CT2 + __dsqrt_rn(__fma_rn(CT2, CT2, 1.0));
    c = __ddiv_rn(cot, __dsqrt_rn(__fma_rn(cot, cot, 1.0)));
    t = copysign(__drcp_rn(cot), ct2);
    return (c < 1.0);
  }
}
MYDEVFN int
dRotH(const double Apq, const double Dp, const double Dq, double &c, double &t)
{
  const double h = Dq + Dp;
  const double ct2 = -__ddiv_rn(h, scalbn(Apq, 1)); // -h / (2 * Apq)
  const double CT2 = fabs(ct2);
  if (CT2 >= SQRT_2_HEPS) {
    const double cot = scalbn(CT2, 1);
    c = 1.0;
    t = copysign(__drcp_rn(cot), ct2);
    return 0;
  }
  else if (CT2 > 1.0) {
    const double cot = CT2 +  __dsqrt_rn(__fma_rn(CT2, CT2, -1.0));
    c = __ddiv_rn(cot, __dsqrt_rn(__fma_rn(cot, cot, -1.0)));
    t = copysign(__drcp_rn(cot), ct2);
    return (c > 1.0);
  }
  else {
    c = HYPJAC_COSHI_FIX;
    t = HYPJAC_TANH1_FIX;
    return 1;
  }
}
#endif // USE_RSQRT

// [  C S ] [ F ] = [ R ]
// [ -S C ] [ G ]   [ 0 ]
// R >= 0
MYDEVFN void dGivens
(const double f,
 const double g,
 double &c,
 double &s,
 double &r)
{
  const double f_ = fabs(f);
  const double g_ = fabs(g);
  if (f_ >= g_) {
    if (g_ == 0.0) {
      c = copysign(1.0, f);
      s = 0.0;
      r = f_;
    }
    else {
      const double g_f = __ddiv_rn(g_, f_);
      r = f_ * __dsqrt_rn(__fma_rn(g_f, g_f, 1.0));
      c = __ddiv_rn(f, r);
      s = __ddiv_rn(g, r);
    }
  }
  else {
    if (f_ == 0.0) {
      c = 0.0;
      s = copysign(1.0, g);
      r = g_;
    }
    else {
      const double f_g = __ddiv_rn(f_, g_);
      r = g_ * __dsqrt_rn(__fma_rn(f_g, f_g, 1.0));
      c = __ddiv_rn(f, r);
      s = __ddiv_rn(g, r);
    }
  }
}

#endif // !DEVICE_CODE_COMMON_ROTATE_HPP
