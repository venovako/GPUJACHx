#ifndef DEVICE_CODE_NOSORT_HPP
#define DEVICE_CODE_NOSORT_HPP

MYDEVFN unsigned dDefJacL0
(volatile double *const G,
 volatile double *const V,
 const unsigned x,
 const unsigned y)
{
  //Fermi: const unsigned y2 = (y << 1u);
  //Fermi: volatile double *const shPtr = &(F32(V, 0u, y2));
  unsigned
    blk_transf_s = 0u,
    blk_transf_b = 0u;

  for (unsigned swp = 0u; swp < _nSwp; ++swp) {
    int
      swp_transf_s = 0,
      swp_transf_b = 0;

    for (unsigned step = 0u; step < _STRAT0_STEPS; ++step) {
      const unsigned
        p = _strat0[step][y][0u],
        q = _strat0[step][y][1u];

      double Dp = +0.0;
      double Dq = +0.0;
      double Apq = +0.0;

      const double Gp = F32(G, x, p);
      const double Gq = F32(G, x, q);

      const double Vp = F32(V, x, p);
      const double Vq = F32(V, x, q);

      __syncthreads();

      Dp = __fma_rn(Gp, Gp, Dp);
      Dq = __fma_rn(Gq, Gq, Dq);
      Apq = __fma_rn(Gp, Gq, Apq);

      Dp = dSum32(Dp);
      Dq = dSum32(Dq);
      Apq = dSum32(Apq);
      //Fermi: Dp = dSum32(Dp, shPtr, x);
      //Fermi: Dq = dSum32(Dq, shPtr, x);
      //Fermi: Apq = dSum32(Apq, shPtr, x);

      const double
        Dp_ = __dsqrt_rn(Dp),
        Dq_ = __dsqrt_rn(Dq),
        Apq_ = fabs(Apq),
        Dpq_ = Dp_ * Dq_;

      const int transf_s = !(Apq_ < (Dpq_ * HYPJAC_MYTOL));
      swp_transf_s += (__syncthreads_count(transf_s) >> WARP_SZ_LGi);

      int transf_b = 0;
      if (transf_s) {
        double c, t;
        transf_b = dRotT(Apq, Dp, Dq, c, t);
        const double t_ = -t;

        if (transf_b) {
          F32(G, x, p) = c * __fma_rn(t_, Gq, Gp);
          F32(G, x, q) = c * __fma_rn(t, Gp, Gq);
          F32(V, x, p) = c * __fma_rn(t_, Vq, Vp);
          F32(V, x, q) = c * __fma_rn(t, Vp, Vq);
        }
        else {
          F32(G, x, p) = __fma_rn(t_, Gq, Gp);
          F32(G, x, q) = __fma_rn(t, Gp, Gq);
          F32(V, x, p) = __fma_rn(t_, Vq, Vp);
          F32(V, x, q) = __fma_rn(t, Vp, Vq);
        }
      }
      else {
        // must restore V
        F32(V, x, p) = Vp;
        F32(V, x, q) = Vq;
      }

      swp_transf_b += (__syncthreads_count(transf_b) >> WARP_SZ_LGi);
    }

    if (swp_transf_s) {
      blk_transf_s += static_cast<unsigned>(swp_transf_s);
      blk_transf_b += static_cast<unsigned>(swp_transf_b);
    }
    else
      break;
  }

  if (!y && !x && blk_transf_s) {
    if (blk_transf_b) {
      unsigned long long blk_transf = blk_transf_b;
      blk_transf <<= 32u;
      blk_transf |= blk_transf_s;
      atomicAdd((unsigned long long*)_cvg, blk_transf);
    }
    else
      atomicAdd((unsigned*)_cvg, blk_transf_s);
  }

  __syncthreads();
  return blk_transf_s;
}

MYDEVFN unsigned dHypJacL0
(volatile double *const G,
 volatile double *const V,
 const unsigned x,
 const unsigned y,
 const unsigned npos)
{
  //Fermi: const unsigned y2 = (y << 1u);
  //Fermi: volatile double *const shPtr = &(F32(V, 0u, y2));
  unsigned
    blk_transf_s = 0u,
    blk_transf_b = 0u;

  for (unsigned swp = 0u; swp < _nSwp; ++swp) {
    int
      swp_transf_s = 0,
      swp_transf_b = 0;

    for (unsigned step = 0u; step < _STRAT0_STEPS; ++step) {
      const unsigned
        p = _strat0[step][y][0u],
        q = _strat0[step][y][1u];

      double Dp = +0.0;
      double Dq = +0.0;
      double Apq = +0.0;

      const double Gp = F32(G, x, p);
      const double Gq = F32(G, x, q);

      const double Vp = F32(V, x, p);
      const double Vq = F32(V, x, q);

      __syncthreads();

      Dp = __fma_rn(Gp, Gp, Dp);
      Dq = __fma_rn(Gq, Gq, Dq);
      Apq = __fma_rn(Gp, Gq, Apq);

      Dp = dSum32(Dp);
      Dq = dSum32(Dq);
      Apq = dSum32(Apq);
      //Fermi: Dp = dSum32(Dp, shPtr, x);
      //Fermi: Dq = dSum32(Dq, shPtr, x);
      //Fermi: Apq = dSum32(Apq, shPtr, x);

      const double
        Dp_ = __dsqrt_rn(Dp),
        Dq_ = __dsqrt_rn(Dq),
        Apq_ = fabs(Apq),
        Dpq_ = Dp_ * Dq_;

      const int transf_s = !(Apq_ < (Dpq_ * HYPJAC_MYTOL));
      swp_transf_s += (__syncthreads_count(transf_s) >> WARP_SZ_LGi);

      int transf_b = 0;
      if (transf_s) {
        const bool
          trig = (((p < npos) && (q < npos)) || ((p >= npos) && (q >= npos)));

        double c, t, t_;
        if (trig) { // trig
          transf_b = dRotT(Apq, Dp, Dq, c, t);
          t_ = -t;
        }
        else { // hyp
          transf_b = dRotH(Apq, Dp, Dq, c, t);
          t_ = t;
        }

        if (transf_b) {
          F32(G, x, p) = c * __fma_rn(t_, Gq, Gp);
          F32(G, x, q) = c * __fma_rn(t, Gp, Gq);
          F32(V, x, p) = c * __fma_rn(t_, Vq, Vp);
          F32(V, x, q) = c * __fma_rn(t, Vp, Vq);
        }
        else {
          F32(G, x, p) = __fma_rn(t_, Gq, Gp);
          F32(G, x, q) = __fma_rn(t, Gp, Gq);
          F32(V, x, p) = __fma_rn(t_, Vq, Vp);
          F32(V, x, q) = __fma_rn(t, Vp, Vq);
        }
      }
      else {
        // must restore V
        F32(V, x, p) = Vp;
        F32(V, x, q) = Vq;
      }

      swp_transf_b += (__syncthreads_count(transf_b) >> WARP_SZ_LGi);
    }

    if (swp_transf_s) {
      blk_transf_s += static_cast<unsigned>(swp_transf_s);
      blk_transf_b += static_cast<unsigned>(swp_transf_b);
    }
    else
      break;
  }

  if (!y && !x && blk_transf_s) {
    if (blk_transf_b) {
      unsigned long long blk_transf = blk_transf_b;
      blk_transf <<= 32u;
      blk_transf |= blk_transf_s;
      atomicAdd((unsigned long long*)_cvg, blk_transf);
    }
    else
      atomicAdd((unsigned*)_cvg, blk_transf_s);
  }

  __syncthreads();
  return blk_transf_s;
}

#endif /* !DEVICE_CODE_NOSORT_HPP */
