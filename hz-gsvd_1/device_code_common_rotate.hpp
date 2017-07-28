#ifndef DEVICE_CODE_COMMON_ROTATE_HPP
#define DEVICE_CODE_COMMON_ROTATE_HPP

MYDEVFN void
dRot(double &App, double &Aqq,
     const double Apq, const double Bpq, const double Bpq_,
     double &CosF, double &SinF, double &CosP, double &SinP)
{
  const double
    Bpqp = __dsqrt_rn(1.0 + Bpq),
    Bpqm = __dsqrt_rn(1.0 - Bpq),
    Xi = __ddiv_rn(Bpq, (Bpqp + Bpqm)),
    Xi_ = -Xi,
    Eta = __ddiv_rn(Bpq, ((1.0 + Bpqp) * (1.0 + Bpqm))),
    Eta_ = -Eta;

  if (Bpq_ < SQRT_HEPS) {
    const double
      Cot2T = __ddiv_rn((Aqq - App), __fma_rn(Bpq, -(App + Aqq), scalbn(Apq, 1))),
      Cot2T_ = fabs(Cot2T);
    double TanT;
    if (Cot2T_ >= SQRT_2_HEPS)
      TanT = scalbn(__drcp_rn(Cot2T), -1);
    else if (Cot2T_ < SQRT_HEPS)
      TanT = copysign(__drcp_rn(Cot2T_ + 1.0), Cot2T);
    else
      TanT = copysign(__drcp_rn(Cot2T_ + __dsqrt_rn(__fma_rn(Cot2T, Cot2T, 1.0))), Cot2T);
    if (fabs(TanT) < SQRT_HEPS) {
      CosF = __fma_rn(TanT - Eta, Xi, 1.0);
      SinF = __fma_rn(__fma_rn(TanT, Eta, 1.0), Xi_, TanT);
      CosP = __fma_rn(TanT + Eta, Xi_, 1.0);
      SinP = __fma_rn(__fma_rn(TanT, Eta_, 1.0), Xi, TanT);
    }
    else {
      const double
        CosT = my_drsqrt_rn(__fma_rn(TanT, TanT, 1.0)),
        SinT = CosT * TanT;
      CosF = __fma_rn(__fma_rn(CosT, Eta_, SinT), Xi, CosT);
      SinF = __fma_rn(__fma_rn(SinT, Eta, CosT), Xi_, SinT);
      CosP = __fma_rn(__fma_rn(CosT, Eta, SinT), Xi_, CosT);
      SinP = __fma_rn(__fma_rn(SinT, Eta_, CosT), Xi, SinT);
    }
  }
  else {
    const double
      F = my_drsqrt_rn(__fma_rn(Bpq, -Bpq, 1.0)),
      Cot2T = __ddiv_rn(Aqq - App, __fma_rn(Bpq, -(App + Aqq), scalbn(Apq, 1)) * F),
      Cot2T_ = fabs(Cot2T);
    double TanT;
    if (Cot2T_ >= SQRT_2_HEPS)
      TanT = scalbn(__drcp_rn(Cot2T), -1);
    else if (Cot2T_ < SQRT_HEPS)
      TanT = copysign(__drcp_rn(Cot2T_ + 1.0), Cot2T);
    else
      TanT = copysign(__drcp_rn(Cot2T_ + __dsqrt_rn(__fma_rn(Cot2T, Cot2T, 1.0))), Cot2T);
    if (fabs(TanT) < SQRT_HEPS) {
      CosF = __fma_rn(TanT - Eta, Xi, 1.0) * F;
      SinF = __fma_rn(__fma_rn(TanT, Eta, 1.0), Xi_, TanT) * F;
      CosP = __fma_rn(TanT + Eta, Xi_, 1.0) * F;
      SinP = __fma_rn(__fma_rn(TanT, Eta_, 1.0), Xi, TanT) * F;
    }
    else {
      const double
        CosT = my_drsqrt_rn(__fma_rn(TanT, TanT, 1.0)),
        SinT = CosT * TanT;
      CosF = __fma_rn(__fma_rn(CosT, Eta_, SinT), Xi, CosT) * F;
      SinF = __fma_rn(__fma_rn(SinT, Eta, CosT), Xi_, SinT) * F;
      CosP = __fma_rn(__fma_rn(CosT, Eta, SinT), Xi_, CosT) * F;
      SinP = __fma_rn(__fma_rn(SinT, Eta_, CosT), Xi, SinT) * F;
    }
  }

  App = CosF*CosF*App - scalbn(CosF*SinP*Apq, 1) + SinP*SinP*Aqq;
  Aqq = SinF*SinF*App + scalbn(SinF*CosP*Apq, 1) + CosP*CosP*Aqq;
}

#endif // !DEVICE_CODE_COMMON_ROTATE_HPP
