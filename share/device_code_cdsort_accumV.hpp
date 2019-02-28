#ifndef DEVICE_CODE_CDSORT_ACCUMV_HPP
#define DEVICE_CODE_CDSORT_ACCUMV_HPP

MYKERN __launch_bounds__(HYPJACL1_MAX_THREADS_PER_BLOCK, HYPJACL1_MIN_BLOCKS_PER_SM)
  dDefJacL1sv(const unsigned step, const int definite)
{
  __shared__ double shMem[2048u];

  const unsigned
    x = threadIdx.x,
    y0 = threadIdx.y,
    iBlk = _strat1[step][blockIdx.x][0u],
    jBlk = _strat1[step][blockIdx.x][1u],
    y1 = y0 + 16u;

  const unsigned
    ix = iBlk * 16u + y0,
    jx = jBlk * 16u + y0;

  double
    *const G0 = _G + ix * _ldG,
    *const G1 = _G + jx * _ldG,
    *const V0 = _V + ix * _ldV,
    *const V1 = _V + jx * _ldV;
  volatile double
    *const G = shMem;
  volatile double
    *const V = shMem + 1024u;

  dFactorize(G0, G1, G, V, x, y0, y1);
  (void)dDefJacL0s(G, V, x, y0, definite);
  dMultV(G0, G1, V0, V1, G, V, x, y0, y1);
}

MYKERN __launch_bounds__(HYPJACL1_MAX_THREADS_PER_BLOCK, HYPJACL1_MIN_BLOCKS_PER_SM)
  dHypJacL1sv(const unsigned step, const unsigned nplus)
{
  __shared__ double shMem[2048u];

  const unsigned
    x = threadIdx.x,
    y0 = threadIdx.y,
    iBlk = _strat1[step][blockIdx.x][0u],
    jBlk = _strat1[step][blockIdx.x][1u],
    y1 = y0 + 16u;

  const unsigned
    ifc0 = iBlk * 16u,
    ifc1 = jBlk * 16u,
    ix = ifc0 + y0,
    jx = ifc1 + y0;

  unsigned npos;
  int definite;
  // blocks negative definite
  if (nplus <= ifc0) {
    npos = 0u;
    definite = -1;
  }
  else if (nplus <= (ifc0 + HYPJACL1_NCOLB)) {
    npos = nplus - ifc0;
    definite = 0;
  }
  else if (nplus <= ifc1) {
    npos = HYPJACL1_NCOLB;
    definite = 0;
  }
  else if (nplus <= (ifc1 + HYPJACL1_NCOLB)) {
    npos = HYPJACL1_NCOLB + (nplus - ifc1);
    definite = 0;
  }
  // blocks positive definite
  else {
    npos = 2u * HYPJACL1_NCOLB;
    definite = +1;
  }

  double
    *const G0 = _G + ix * _ldG,
    *const G1 = _G + jx * _ldG,
    *const V0 = _V + ix * _ldV,
    *const V1 = _V + jx * _ldV;
  volatile double
    *const G = shMem;
  volatile double
    *const V = shMem + 1024u;

  dFactorize(G0, G1, G, V, x, y0, y1);
  if (definite)
    (void)dDefJacL0s(G, V, x, y0, definite);
  else
    (void)dHypJacL0s(G, V, x, y0, npos);
  dMultV(G0, G1, V0, V1, G, V, x, y0, y1);
}

#endif // !DEVICE_CODE_CDSORT_ACCUMV_HPP
