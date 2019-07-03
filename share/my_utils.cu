#include "my_utils.hpp"

#ifdef _WIN32
EXTERN_C void __stdcall GetSystemTimePreciseAsFileTime(long long*);
EXTERN_C void* __stdcall LoadLibraryA(const char*);
EXTERN_C void* __stdcall GetProcAddress(void*, const char*);
EXTERN_C int __stdcall FreeLibrary(void*);
#else // POSIX -ldl
#include <dlfcn.h>
#endif //? _WIN32

char err_msg[err_msg_size] = { '\0' };

int fexist(const char *const fn) throw()
{
#ifdef _WIN32
  // _stat under Windows fails
  // http://stackoverflow.com/questions/15827954/windows-with-mapped-drive-stat-fails-to-detect-existence-of-the-file
  if (FILE *const f = (fn ? fopen(fn, "rb") : static_cast<FILE*>(NULL))) {
    (void)fclose(f);
    return 1;
  }
  return 0;
#else // POSIX
  struct stat buf;
  return (fn ? (0 == stat(fn, &buf)) : 0);
#endif // ?_WIN32
}

void *strat_open(const char *const sdy) throw()
{
  return
#ifdef _WIN32
    LoadLibraryA(sdy);
#else // POSIX
    dlopen(sdy, RTLD_LAZY);
#endif // ?_WIN32
}

int strat_close(void *const h) throw()
{
  return
#ifdef _WIN32
    !FreeLibrary(h);
#else // POSIX
    dlclose(h);
#endif // ?_WIN32
}

const void *strat_ptr(void *const h, const char *const snp, const unsigned n) throw()
{
  char arrn[12] = { '\0' };
  (void)snprintf(arrn, sizeof(arrn), "%s%05u", snp, n);
  return (arrn[11] ? NULL :
#ifdef _WIN32
    GetProcAddress(h, arrn)
#else // POSIX
    dlsym(h, arrn)
#endif // ?_WIN32
  );
}

long long timestamp() throw()
{
#ifdef _WIN32
  long long ret;
  GetSystemTimePreciseAsFileTime(&ret);
  return ret;
#else // POSIX
  struct timeval tv;
  SYSI_CALL(gettimeofday(&tv, static_cast<struct timezone*>(NULL)));
  return (tv.tv_sec * TS_S + tv.tv_usec);
#endif // ?_WIN32
}

void stopwatch_reset(long long &sw) throw()
{
  sw = timestamp();
}

long long stopwatch_lap(long long &sw) throw()
{
  const long long
    ts = timestamp(),
    lap = ts - sw;
  sw = ts;
  return lap;
}
