// defines.hpp: macro definitions and system header includes.

#ifndef DEFINES_HPP
#define DEFINES_HPP

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif // !_GNU_SOURCE

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

#include <alloca.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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

#ifdef USE_MULTI_THR
#define TLS __thread
#else // !USE_MULTI_THR
#define TLS
#endif // USE_MULTI_THR

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

#define Long long
#define Ldiv ldiv
#define Ldiv_t ldiv_t
#define FmtLong "l"
#define MkLong(x) x ## l

#ifdef VAR_UNUSED
#error VAR_UNUSED not definable externally
#endif // VAR_UNUSED
#define VAR_UNUSED __attribute__ ((unused))

#endif // !DEFINES_HPP
