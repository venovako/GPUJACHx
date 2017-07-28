#include "my_utils.hpp"

#include <dlfcn.h>

TLS char err_msg[err_msg_size] = { '\0' };

int fexist(const char *const fn) throw()
{
  struct stat buf;
  return (fn ? (0 == stat(fn, &buf)) : 0);
}

void *strat_open(const char *const sdy) throw()
{
  return dlopen(sdy, RTLD_LAZY);
}

int strat_close(void *const h) throw()
{
  return dlclose(h);
}

const void *strat_ptr(void *const h, const char *const snp, const unsigned n) throw()
{
  char arrn[12] = { '\0' };
  (void)snprintf(arrn, sizeof(arrn), "%s%05u", snp, n);
  return (arrn[11] ? NULL : dlsym(h, arrn));
}

Long timestamp() throw()
{
  struct timeval tv;
  SYSI_CALL(gettimeofday(&tv, static_cast<struct timezone*>(NULL)));
  return (tv.tv_sec * TS_S + tv.tv_usec);
}

void stopwatch_reset(Long &sw) throw()
{
  sw = timestamp();
}

Long stopwatch_lap(Long &sw) throw()
{
  const Long
    ts = timestamp(),
    lap = ts - sw;
  sw = ts;
  return lap;
}
