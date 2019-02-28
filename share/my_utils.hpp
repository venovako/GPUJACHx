#ifndef MY_UTILS_HPP
#define MY_UTILS_HPP

#include "defines.hpp"

#ifndef err_msg_size
#define err_msg_size static_cast<size_t>(1024u)
#else // err_msg_size
#error err_msg_size not definable externally
#endif // !err_msg_size

EXTERN_C TLS char err_msg[err_msg_size];

#ifndef WARN
#define WARN(msg) {                                                             \
    (void)fprintf(stderr, "[WARNING] %s(%d): %s\n", __FILE__, __LINE__, (msg)); \
  }
#else // WARN
#error WARN not definable externally
#endif // !WARN

#ifndef DIE
#define DIE(msg) {                                                            \
    (void)fprintf(stderr, "[ERROR] %s(%d): %s\n", __FILE__, __LINE__, (msg)); \
    exit(EXIT_FAILURE);                                                       \
  }
#else // DIE
#error DIE not definable externally
#endif // !DIE

#ifndef SYSI_CALL
#define SYSI_CALL(call) {						\
    if (0 != static_cast<int>(call)) {					\
      (void)fprintf(stderr, "[ERROR] %s(%d): %s",                       \
		    __FILE__, __LINE__, strerror(errno));		\
      exit(EXIT_FAILURE);                                               \
    }									\
  }
#else // SYSI_CALL
#error SYSI_CALL not definable externally
#endif // !SYSI_CALL

#ifndef SYSP_CALL
#define SYSP_CALL(call) {						\
    if (NULL == static_cast<const void*>(call)) {			\
      (void)fprintf(stderr, "[ERROR] %s(%d): %s",                       \
		    __FILE__, __LINE__, strerror(errno));		\
      exit(EXIT_FAILURE);                                               \
    }									\
  }
#else
#error SYSP_CALL not definable externally
#endif // !SYSP_CALL

extern int fexist(const char *const fn) throw();

extern void *strat_open(const char *const sdy) throw();
extern int strat_close(void *const h) throw();
extern const void *strat_ptr(void *const h, const char *const snp, const unsigned n) throw();

template <typename T>
T udiv_ceil(const T a, const T b) throw()
{
  return (a + b - static_cast<T>(1u)) / b;
}

template <typename T>
T inc_mod(const T a, const T m) throw()
{
  return (m ? ((a < (m - static_cast<T>(1u))) ? (a + static_cast<T>(1u)) : static_cast<T>(0u)) : a);
}

template <typename T>
T dec_mod(const T a, const T m) throw()
{
  return (m ? (a ? (a - static_cast<T>(1u)) : (m - static_cast<T>(1u))) : a);
}

inline unsigned
bc2c(const unsigned bc, const unsigned nc_bc) throw()
{
  return (bc * nc_bc);
}

template <typename T>
inline T*
cA(T *const A, const unsigned c, const unsigned ldA) throw()
{
  return (A + c * ldA);
}

#ifndef TS2S
#ifdef _WIN32
#define TS2S 1e-7
#else // POSIX
#define TS2S 1e-6
#endif // _WIN32
#else // TS2S
#error TS2S not definable externally
#endif // !TS2S

#ifndef TS_S
#ifdef _WIN32
#define TS_S 10000000ll
#else // POSIX
#define TS_S 1000000ll
#endif // _WIN32
#else // TS_S
#error TS_S not definable externally
#endif // !TS_S

extern long long timestamp() throw();
extern void stopwatch_reset(long long &sw) throw();
extern long long stopwatch_lap(long long &sw) throw();

#endif // !MY_UTILS_HPP
