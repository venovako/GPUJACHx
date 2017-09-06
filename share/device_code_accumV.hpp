#ifndef DEVICE_CODE_ACCUMV_HPP
#define DEVICE_CODE_ACCUMV_HPP

MYDEVFN void dMultV
(
 double *const G0,
 double *const G1,
 double *const V0,
 double *const V1,
 volatile double *const A,
 const double *const B,
 const unsigned x,
 const unsigned y0,
 const unsigned y1
)
{
  dMultAV(G0, G1, A, B, x, y0, y1);
  dMultAV(V0, V1, A, B, x, y0, y1);
}

#endif // !DEVICE_CODE_ACCUMV_HPP
