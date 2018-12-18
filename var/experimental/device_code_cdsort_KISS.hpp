#ifndef DEVICE_CODE_CDSORT_KISS_HPP
#define DEVICE_CODE_CDSORT_KISS_HPP

MYDEVFN unsigned dDefJacL0posd
(
 volatile double *const G,
 volatile double *const V,
#if __CUDA_ARCH__ >= 300
#else // Fermi
 volatile double *const shPtr,
#endif // __CUDA_ARCH__
 const unsigned x,
 const unsigned y
)
{
  unsigned
    blk_transf_G = 0u,
    blk_transf_V = 0u;

  for (unsigned swp = 0u; swp < _nSwp; ++swp) {
    int
      swp_transf_G = 0,
      swp_transf_V = 0;

    for (unsigned step = 0u; step < _STRAT0_STEPS; ++step) {
      const unsigned
        p = _strat0[step][y][0u],
        q = _strat0[step][y][1u];

      const double Gp = F32(G, x, p);
      const double Gq = F32(G, x, q);

      const double Vp = F32(V, x, p);
      const double Vq = F32(V, x, q);

      __syncthreads();

#if __CUDA_ARCH__ >= 300
      double Dp = dSum32(Gp * Gp);
      double Dq = dSum32(Gq * Gq);
      const double Apq = dSum32(Gp * Gq);
#else // Fermi
      double Dp = dSum32(Gp * Gp, shPtr, x);
      double Dq = dSum32(Gq * Gq, shPtr, x);
      const double Apq = dSum32(Gp * Gq, shPtr, x);
#endif // __CUDA_ARCH__

      const double
        Dp_ = __dsqrt_rn(Dp),
        Dq_ = __dsqrt_rn(Dq),
        Apq_ = fabs(Apq),
        Dpq_ = Dp_ * Dq_;

      const int transf_G = !(Apq_ < (Dpq_ * HYPJAC_MYTOL));
      int transf_V = 0;
      swp_transf_G += (__syncthreads_count(transf_G) >> WARP_SZ_LGi);

      if (transf_G) {
        double c, t;
        transf_V = dRotT(Apq, Dp, Dq, c, t);
        const double t_ = -t;

        Dp = __fma_rn(t, Apq, Dp);
        Dq = __fma_rn(t_, Apq, Dq);

        if (transf_V) {
          if (Dp >= Dq) {
            const double Gp_ = c * __fma_rn(t_, Gq, Gp);
            F32(G, x, p) = Gp_;
            const double Gq_ = c * __fma_rn(t, Gp, Gq);
            F32(G, x, q) = Gq_;
            const double Vp_ = c * __fma_rn(t_, Vq, Vp);
            F32(V, x, p) = Vp_;
            const double Vq_ = c * __fma_rn(t, Vp, Vq);
            F32(V, x, q) = Vq_;
          }
          else {
            const double Gq_ = c * __fma_rn(t_, Gq, Gp);
            F32(G, x, q) = Gq_;
            const double Gp_ = c * __fma_rn(t, Gp, Gq);
            F32(G, x, p) = Gp_;
            const double Vq_ = c * __fma_rn(t_, Vq, Vp);
            F32(V, x, q) = Vq_;
            const double Vp_ = c * __fma_rn(t, Vp, Vq);
            F32(V, x, p) = Vp_;
          }
        }
        else {
          if (Dp >= Dq) {
            const double Gp_ = __fma_rn(t_, Gq, Gp);
            F32(G, x, p) = Gp_;
            const double Gq_ = __fma_rn(t, Gp, Gq);
            F32(G, x, q) = Gq_;
            const double Vp_ = __fma_rn(t_, Vq, Vp);
            F32(V, x, p) = Vp_;
            const double Vq_ = __fma_rn(t, Vp, Vq);
            F32(V, x, q) = Vq_;
            transf_V = (Vp != Vp_) || (Vq != Vq_);
          }
          else {
            const double Gq_ = __fma_rn(t_, Gq, Gp);
            F32(G, x, q) = Gq_;
            const double Gp_ = __fma_rn(t, Gp, Gq);
            F32(G, x, p) = Gp_;
            const double Vq_ = __fma_rn(t_, Vq, Vp);
            F32(V, x, q) = Vq_;
            const double Vp_ = __fma_rn(t, Vp, Vq);
            F32(V, x, p) = Vp_;
            transf_V = (Vq != Vp_) || (Vp != Vq_);
          }
        }
      }
      else if (Dp < Dq) {
        F32(G, x, p) = Gq;
        F32(G, x, q) = Gp;
        F32(V, x, p) = Vq;
        F32(V, x, q) = Vp;
      }
#if __CUDA_ARCH__ >= 300
#else // Fermi
      else {
        // must restore V
        F32(V, x, p) = Vp;
        F32(V, x, q) = Vq;
      }
#endif // __CUDA_ARCH__

      swp_transf_V += (__syncthreads_count(transf_V) >> WARP_SZ_LGi);
    }

    if (swp_transf_G)
      blk_transf_G += static_cast<unsigned>(swp_transf_G);
    else
      break;
    blk_transf_V += static_cast<unsigned>(swp_transf_V);
  }

  if (!y && !x) {
    if (blk_transf_V) {
      unsigned Long blk_transf = blk_transf_V;
      blk_transf <<= 32u;
      blk_transf |= blk_transf_G;
      atomicAdd((unsigned long long*)_cvg, blk_transf);
    }
    else if (blk_transf_G)
      atomicAdd((unsigned*)_cvg, blk_transf_G);
  }

  __syncthreads();
  return blk_transf_G;
}

MYDEVFN unsigned dDefJacL0negd
(
 volatile double *const G,
 volatile double *const V,
#if __CUDA_ARCH__ >= 300
#else // Fermi
 volatile double *const shPtr,
#endif // __CUDA_ARCH__
 const unsigned x,
 const unsigned y
)
{
  unsigned
    blk_transf_G = 0u,
    blk_transf_V = 0u;

  for (unsigned swp = 0u; swp < _nSwp; ++swp) {
    int
      swp_transf_G = 0,
      swp_transf_V = 0;

    for (unsigned step = 0u; step < _STRAT0_STEPS; ++step) {
      const unsigned
        p = _strat0[step][y][0u],
        q = _strat0[step][y][1u];

      const double Gp = F32(G, x, p);
      const double Gq = F32(G, x, q);

      const double Vp = F32(V, x, p);
      const double Vq = F32(V, x, q);

      __syncthreads();

#if __CUDA_ARCH__ >= 300
      double Dp = dSum32(Gp * Gp);
      double Dq = dSum32(Gq * Gq);
      const double Apq = dSum32(Gp * Gq);
#else // Fermi
      double Dp = dSum32(Gp * Gp, shPtr, x);
      double Dq = dSum32(Gq * Gq, shPtr, x);
      const double Apq = dSum32(Gp * Gq, shPtr, x);
#endif // __CUDA_ARCH__

      const double
        Dp_ = __dsqrt_rn(Dp),
        Dq_ = __dsqrt_rn(Dq),
        Apq_ = fabs(Apq),
        Dpq_ = Dp_ * Dq_;

      const int transf_G = !(Apq_ < (Dpq_ * HYPJAC_MYTOL));
      int transf_V = 0;
      swp_transf_G += (__syncthreads_count(transf_G) >> WARP_SZ_LGi);

      if (transf_G) {
        double c, t;
        transf_V = dRotT(Apq, Dp, Dq, c, t);
        const double t_ = -t;

        Dp = __fma_rn(t, Apq, Dp);
        Dq = __fma_rn(t_, Apq, Dq);

        if (transf_V) {
          if (Dp <= Dq) {
            const double Gp_ = c * __fma_rn(t_, Gq, Gp);
            F32(G, x, p) = Gp_;
            const double Gq_ = c * __fma_rn(t, Gp, Gq);
            F32(G, x, q) = Gq_;
            const double Vp_ = c * __fma_rn(t_, Vq, Vp);
            F32(V, x, p) = Vp_;
            const double Vq_ = c * __fma_rn(t, Vp, Vq);
            F32(V, x, q) = Vq_;
          }
          else {
            const double Gq_ = c * __fma_rn(t_, Gq, Gp);
            F32(G, x, q) = Gq_;
            const double Gp_ = c * __fma_rn(t, Gp, Gq);
            F32(G, x, p) = Gp_;
            const double Vq_ = c * __fma_rn(t_, Vq, Vp);
            F32(V, x, q) = Vq_;
            const double Vp_ = c * __fma_rn(t, Vp, Vq);
            F32(V, x, p) = Vp_;
          }
        }
        else {
          if (Dp <= Dq) {
            const double Gp_ = __fma_rn(t_, Gq, Gp);
            F32(G, x, p) = Gp_;
            const double Gq_ = __fma_rn(t, Gp, Gq);
            F32(G, x, q) = Gq_;
            const double Vp_ = __fma_rn(t_, Vq, Vp);
            F32(V, x, p) = Vp_;
            const double Vq_ = __fma_rn(t, Vp, Vq);
            F32(V, x, q) = Vq_;
            transf_V = (Vp != Vp_) || (Vq != Vq_);
          }
          else {
            const double Gq_ = __fma_rn(t_, Gq, Gp);
            F32(G, x, q) = Gq_;
            const double Gp_ = __fma_rn(t, Gp, Gq);
            F32(G, x, p) = Gp_;
            const double Vq_ = __fma_rn(t_, Vq, Vp);
            F32(V, x, q) = Vq_;
            const double Vp_ = __fma_rn(t, Vp, Vq);
            F32(V, x, p) = Vp_;
            transf_V = (Vq != Vp_) || (Vp != Vq_);
          }
        }
      }
      else if (Dq < Dp) {
        F32(G, x, p) = Gq;
        F32(G, x, q) = Gp;
        F32(V, x, p) = Vq;
        F32(V, x, q) = Vp;
      }
#if __CUDA_ARCH__ >= 300
#else // Fermi
      else {
        // must restore V
        F32(V, x, p) = Vp;
        F32(V, x, q) = Vq;
      }
#endif

      swp_transf_V += (__syncthreads_count(transf_V) >> WARP_SZ_LGi);
    }

    if (swp_transf_G)
      blk_transf_G += static_cast<unsigned>(swp_transf_G);
    else
      break;
    blk_transf_V += static_cast<unsigned>(swp_transf_V);
  }

  if (!y && !x) {
    if (blk_transf_V) {
      unsigned Long blk_transf = blk_transf_V;
      blk_transf <<= 32u;
      blk_transf |= blk_transf_G;
      atomicAdd((unsigned long long*)_cvg, blk_transf);
    }
    else if (blk_transf_G)
      atomicAdd((unsigned*)_cvg, blk_transf_G);
  }

  __syncthreads();
  return blk_transf_G;
}

#endif // !DEVICE_CODE_CDSORT_KISS_HPP
