#ifndef DEVICE_CODE_ACCUMV_HPP
#define DEVICE_CODE_ACCUMV_HPP

MYDEVFN void dMultV
(
 double *const F0,
 double *const F1,
 double *const G0,
 double *const G1,
 double *const V0,
 double *const V1,
 volatile double *const A,
 volatile double *const B,
 const double *const C,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  dMultAV(F0, F1, A, C, x, y0, y1);
  dMultAV(G0, G1, B, C, x, y0, y1);
  dMultAV(V0, V1, A, C, x, y0, y1);
}

#endif // !DEVICE_CODE_ACCUMV_HPP
