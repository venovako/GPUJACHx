#ifndef DEVICE_CODE_CDSORT_DRMAC_HPP
#define DEVICE_CODE_CDSORT_DRMAC_HPP

MYDEVFN unsigned dDefJacL0posd
(volatile double *const G,
 volatile double *const V,
#if __CUDA_ARCH__ >= 300
#else /* Fermi */
 volatile double *const shPtr,
#endif /* ?__CUDA_ARCH__ */
 const unsigned x,
 const unsigned y)
{
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

      const double Gp = F32(G, x, p);
      const double Gq = F32(G, x, q);

      const double Vp = F32(V, x, p);
      const double Vq = F32(V, x, q);

      __syncthreads();

      const double
#if __CUDA_ARCH__ >= 300
        Ap = dNRM2_32(Gp),
        Aq = dNRM2_32(Gq),
        CosA = dCOSA(Ap, Aq, Gp, Gq);
#else /* Fermi */
        Ap = dNRM2_32(Gp, shPtr, x),
        Aq = dNRM2_32(Gq, shPtr, x),
        CosA = dCOSA(Ap, Aq, Gp, Gq, shPtr, x);
#endif /* ?__CUDA_ARCH__ */

      const int transf_s = (fabs(CosA) >= HYPJAC_MYTOL);
      swp_transf_s += (__syncthreads_count(transf_s) >> WARP_SZ_LGi);

      int transf_b = 0;
      if (transf_s) {
        const double
          Aq_p = __ddiv_rn(Aq, Ap),
          Ap_q = __ddiv_rn(Ap, Aq),
          ct2 = tCot2(Aq_p, Ap_q, CosA);
        if (ct2 <= NU_8) {
          const double CT2 = fabs(ct2);
          double c, t;
          if (CT2 < SQRT_HEPS) {
            const double cot = CT2 + 1.0;
            t = copysign(__drcp_rn(cot), ct2);
            c = my_drsqrt_rn(__fma_rn(t, t, 1.0));
            transf_b = 1;
          }
          else if (CT2 >= SQRT_2_HEPS) {
            const double cot = scalbn(CT2, 1);
            c = 1.0;
            t = copysign(__drcp_rn(cot), ct2);
          }
          else {
            const double cot = CT2 + __dsqrt_rn(__fma_rn(CT2, CT2, 1.0));
            t = copysign(__drcp_rn(cot), ct2);
            c = my_drsqrt_rn(__fma_rn(t, t, 1.0));
            transf_b = (c < 1.0);
          }
          const double
            t_ = -t,
            tCosA = t * CosA,
            tCosA_ = -tCosA,
            Ap_ = Ap * __dsqrt_rn(__fma_rn(tCosA_, Aq_p, 1.0)),
            Aq_ = Aq * __dsqrt_rn(__fma_rn(tCosA, Ap_q, 1.0));
          if (transf_b) {
            if (Ap_ >= Aq_) {
              F32(G, x, p) = c * __fma_rn(t_, Gq, Gp);
              F32(G, x, q) = c * __fma_rn(t, Gp, Gq);
              F32(V, x, p) = c * __fma_rn(t_, Vq, Vp);
              F32(V, x, q) = c * __fma_rn(t, Vp, Vq);
            }
            else {
              F32(G, x, q) = c * __fma_rn(t_, Gq, Gp);
              F32(G, x, p) = c * __fma_rn(t, Gp, Gq);
              F32(V, x, q) = c * __fma_rn(t_, Vq, Vp);
              F32(V, x, p) = c * __fma_rn(t, Vp, Vq);
            }
          }
          else {
            if (Ap_ >= Aq_) {
              F32(G, x, p) = __fma_rn(t_, Gq, Gp);
              F32(G, x, q) = __fma_rn(t, Gp, Gq);
              F32(V, x, p) = __fma_rn(t_, Vq, Vp);
              F32(V, x, q) = __fma_rn(t, Vp, Vq);
            }
            else {
              F32(G, x, q) = __fma_rn(t_, Gq, Gp);
              F32(G, x, p) = __fma_rn(t, Gp, Gq);
              F32(V, x, q) = __fma_rn(t_, Vq, Vp);
              F32(V, x, p) = __fma_rn(t, Vp, Vq);
            }
          }
        }
        else if (Ap >= Aq) {
          const double
            t_ = CosA,
            t = -t_;
          F32(G, x, q) = Aq * __fma_rn(t, __ddiv_rn(Gp, Ap), __ddiv_rn(Gq, Aq));
          F32(V, x, p) = __fma_rn(t_, Vq, Vp);
          F32(V, x, q) = __fma_rn(t, Vp, Vq);
        }
        else {
          const double
            t = CosA,
            t_ = -t;
          F32(G, x, q) = Ap * __fma_rn(t_, __ddiv_rn(Gq, Aq), __ddiv_rn(Gp, Ap));
          F32(G, x, p) = Gq;
          F32(V, x, q) = __fma_rn(t_, Vq, Vp);
          F32(V, x, p) = __fma_rn(t, Vp, Vq);
        }
      }
      else if (Ap < Aq) {
        F32(G, x, p) = Gq;
        F32(G, x, q) = Gp;
        F32(V, x, p) = Vq;
        F32(V, x, q) = Vp;
      }
#if __CUDA_ARCH__ >= 300
#else /* Fermi */
      else {
        // must restore V
        F32(V, x, p) = Vp;
        F32(V, x, q) = Vq;
      }
#endif /* ?__CUDA_ARCH__ */

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

MYDEVFN unsigned dDefJacL0negd
(volatile double *const G,
 volatile double *const V,
#if __CUDA_ARCH__ >= 300
#else /* Fermi */
 volatile double *const shPtr,
#endif /* ?__CUDA_ARCH__ */
 const unsigned x,
 const unsigned y)
{
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

      const double Gp = F32(G, x, p);
      const double Gq = F32(G, x, q);

      const double Vp = F32(V, x, p);
      const double Vq = F32(V, x, q);

      __syncthreads();

      const double
#if __CUDA_ARCH__ >= 300
        Ap = dNRM2_32(Gp),
        Aq = dNRM2_32(Gq),
        CosA = dCOSA(Ap, Aq, Gp, Gq);
#else /* Fermi */
        Ap = dNRM2_32(Gp, shPtr, x),
        Aq = dNRM2_32(Gq, shPtr, x),
        CosA = dCOSA(Ap, Aq, Gp, Gq, shPtr, x);
#endif /* ?__CUDA_ARCH__ */

      const int transf_s = (fabs(CosA) >= HYPJAC_MYTOL);
      swp_transf_s += (__syncthreads_count(transf_s) >> WARP_SZ_LGi);

      int transf_b = 0;
      if (transf_s) {
        const double
          Aq_p = __ddiv_rn(Aq, Ap),
          Ap_q = __ddiv_rn(Ap, Aq),
          ct2 = tCot2(Aq_p, Ap_q, CosA);
        if (ct2 <= NU_8) {
          const double CT2 = fabs(ct2);
          double c, t;
          if (CT2 < SQRT_HEPS) {
            const double cot = CT2 + 1.0;
            t = copysign(__drcp_rn(cot), ct2);
            c = my_drsqrt_rn(__fma_rn(t, t, 1.0));
            transf_b = 1;
          }
          else if (CT2 >= SQRT_2_HEPS) {
            const double cot = scalbn(CT2, 1);
            c = 1.0;
            t = copysign(__drcp_rn(cot), ct2);
          }
          else {
            const double cot = CT2 + __dsqrt_rn(__fma_rn(CT2, CT2, 1.0));
            t = copysign(__drcp_rn(cot), ct2);
            c = my_drsqrt_rn(__fma_rn(t, t, 1.0));
            transf_b = (c < 1.0);
          }
          const double
            t_ = -t,
            tCosA = t * CosA,
            tCosA_ = -tCosA,
            Ap_ = Ap * __dsqrt_rn(__fma_rn(tCosA_, Aq_p, 1.0)),
            Aq_ = Aq * __dsqrt_rn(__fma_rn(tCosA, Ap_q, 1.0));
          if (transf_b) {
            if (Ap_ <= Aq_) {
              F32(G, x, p) = c * __fma_rn(t_, Gq, Gp);
              F32(G, x, q) = c * __fma_rn(t, Gp, Gq);
              F32(V, x, p) = c * __fma_rn(t_, Vq, Vp);
              F32(V, x, q) = c * __fma_rn(t, Vp, Vq);
            }
            else {
              F32(G, x, q) = c * __fma_rn(t_, Gq, Gp);
              F32(G, x, p) = c * __fma_rn(t, Gp, Gq);
              F32(V, x, q) = c * __fma_rn(t_, Vq, Vp);
              F32(V, x, p) = c * __fma_rn(t, Vp, Vq);
            }
          }
          else {
            if (Ap_ <= Aq_) {
              F32(G, x, p) = __fma_rn(t_, Gq, Gp);
              F32(G, x, q) = __fma_rn(t, Gp, Gq);
              F32(V, x, p) = __fma_rn(t_, Vq, Vp);
              F32(V, x, q) = __fma_rn(t, Vp, Vq);
            }
            else {
              F32(G, x, q) = __fma_rn(t_, Gq, Gp);
              F32(G, x, p) = __fma_rn(t, Gp, Gq);
              F32(V, x, q) = __fma_rn(t_, Vq, Vp);
              F32(V, x, p) = __fma_rn(t, Vp, Vq);
            }
          }
        }
        else if (Ap <= Aq) {
          const double
            t_ = CosA,
            t = -t_;
          F32(G, x, q) = Aq * __fma_rn(t, __ddiv_rn(Gp, Ap), __ddiv_rn(Gq, Aq));
          F32(V, x, p) = __fma_rn(t_, Vq, Vp);
          F32(V, x, q) = __fma_rn(t, Vp, Vq);
        }
        else {
          const double
            t = CosA,
            t_ = -t;
          F32(G, x, q) = Ap * __fma_rn(t_, __ddiv_rn(Gq, Aq), __ddiv_rn(Gp, Ap));
          F32(G, x, p) = Gq;
          F32(V, x, q) = __fma_rn(t_, Vq, Vp);
          F32(V, x, p) = __fma_rn(t, Vp, Vq);
        }
      }
      else if (Aq < Ap) {
        F32(G, x, p) = Gq;
        F32(G, x, q) = Gp;
        F32(V, x, p) = Vq;
        F32(V, x, q) = Vp;
      }
#if __CUDA_ARCH__ >= 300
#else /* Fermi */
      else {
        // must restore V
        F32(V, x, p) = Vp;
        F32(V, x, q) = Vq;
      }
#endif /* ?__CUDA_ARCH__ */

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

#endif /* !DEVICE_CODE_CDSORT_DRMAC_HPP */
