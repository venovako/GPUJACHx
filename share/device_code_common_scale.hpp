#ifndef DEVICE_CODE_COMMON_SCALE_HPP
#define DEVICE_CODE_COMMON_SCALE_HPP

// find the smallest |scl| s.t.
// if (x<t): 2**scl * x > t
// if (x>t): 2**scl * x < t
MYDEVFN int
find_scl(const double x, const double t)
{
  int x_n, t_n;
  const double x_m = frexp(x, &x_n);
  const double t_m = frexp(t, &t_n);
  return ((x <= t) ?
          ((t_n - x_n) + (x_m < t_m)) :
          ((t_n - x_n) - (x_m > t_m)));
}

MYDEVFN int
find_scl2(const double x)
{
  const int scl = ilogb(x);
  return (scl + (scl & 1));
}

MYDEVFN void
normalize2(int &scl, double &ssq)
{
  const int scl_ = find_scl2(ssq);
  scl += scl_;
  ssq = scalbn(ssq, -scl_);
}

MYDEVFN int
cmp_scl(const int s1, const int s2, const double x1, const double x2)
{
  if (s1 < s2)
    return -1;
  if (s2 < s1)
    return 1;
  if (x1 < x2)
    return -1;
  if (x2 < x1)
    return 1;
  return 0;
}

MYDEVFN void
sum_scl(const int s1, const int s2, const double x1, const double x2, int &s, double &x)
{
  if (cmp_scl(s1, s2, x1, x2) <= 0) {
    s = s2;
    x = scalbn(x1, s1 - s2) + x2;
  }
  else {
    s = s1;
    x = scalbn(x2, s2 - s1) + x1;
  }
}

#endif /* !DEVICE_CODE_COMMON_SCALE_HPP */
