#ifndef DEVICE_CODE_CDSORT_SOLVEV_HPP
#define DEVICE_CODE_CDSORT_SOLVEV_HPP

MYKERN __launch_bounds__(HZ_L1_MAX_THREADS_PER_BLOCK, HZ_L1_MIN_BLOCKS_PER_SM)
  dHZ_L1_s(const unsigned step)
{
  __shared__ double shMem[3u * 1024u];

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
    *const F0 = _F + ix * _ldF,
    *const F1 = _F + jx * _ldF,
    *const G0 = _G + ix * _ldG,
    *const G1 = _G + jx * _ldG,
    *const F = shMem,
    *const G = shMem + 1024u,
    *const V = shMem + 2048u;

  dFactorize(F0, F1, G0, G1, F, G, V, x, y0, y1);
  if (dHZ_L0_s(F, G, V, x, y0)) {
    dMultAV(F0, F1, F, V, x, y0, y1);
    dMultAV(G0, G1, G, V, x, y0, y1);
  }
}

#endif // !DEVICE_CODE_CDSORT_SOLVEV_HPP
