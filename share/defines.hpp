// defines.hpp: macro definitions and system header includes.

#ifndef DEFINES_HPP
#define DEFINES_HPP

#ifdef _WIN32
#ifndef _CRT_NONSTDC_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#endif // !_CRT_NONSTDC_NO_DEPRECATE
#ifndef _CRT_NONSTDC_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#endif // !_CRT_NONSTDC_NO_WARNINGS
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif // !_CRT_SECURE_NO_DEPRECATE
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif // !_CRT_SECURE_NO_WARNINGS
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES
#else // !_WIN32
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif // !_GNU_SOURCE
#endif // _WIN32

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else // NVCC host compiler
#ifdef __cplusplus
#include <cmath>
#include <complex>
#else // C99
#include <math.h>
#include <complex.h>
#endif // __cplusplus
#endif // __INTEL_COMPILER

#ifdef __cplusplus
#include <cassert>
#include <cerrno>
#include <cctype>
#include <cfloat>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#else // C
#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#endif // __cplusplus

#ifdef _WIN32
#include <io.h>
#include <malloc.h>
#ifndef alloca
#define alloca _alloca
#endif // !alloca
#else // POSIX
#include <alloca.h>
#include <sys/time.h>
#endif // _WIN32

#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#ifndef snprintf
#define snprintf _snprintf
#else // snprintf
#error snprintf not definable externally
#endif // !snprintf
#else // POSIX
#include <unistd.h>
#endif // _WIN32

// defines

#ifndef EXTERN_C
#ifdef __cplusplus
#define EXTERN_C extern "C"
#else // C
#define EXTERN_C extern
#endif // __cplusplus
#else // EXTERN_C
#error EXTERN_C not definable externally
#endif // !EXTERN_C

#ifdef TLS
#error TLS not definable externally
#endif // TLS

#ifdef USE_MULTI_GPU
#ifdef _WIN32
#define TLS __declspec(thread)
#else // POSIX
#define TLS __thread
#endif // _WIN32
#else // !USE_MULTI_GPU
#define TLS
#endif // USE_MULTI_GPU

#ifdef Long
#error Long not definable externally
#endif // Long
#ifdef Ldiv
#error Ldiv not definable externally
#endif // Ldiv
#ifdef Ldiv_t
#error Ldiv_t not definable externally
#endif // Ldiv_t
#ifdef FmtLong
#error FmtLong not definable externally
#endif // FmtLong
#ifdef MkLong
#error MkLong not definable externally
#endif // MkLong

#ifdef _WIN32
#define Long long long
#define Ldiv lldiv
#define Ldiv_t lldiv_t
#define FmtLong "ll"
#define MkLong(x) x ## ll
#else // POSIX
#define Long long
#define Ldiv ldiv
#define Ldiv_t ldiv_t
#define FmtLong "l"
#define MkLong(x) x ## l
#endif // _WIN32

#ifdef VAR_UNUSED
#error VAR_UNUSED not definable externally
#endif // VAR_UNUSED

#ifdef _WIN32
#define VAR_UNUSED
#else // POSIX
#define VAR_UNUSED __attribute__ ((unused))
#endif // _WIN32

#endif // !DEFINES_HPP
