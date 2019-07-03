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
#endif // ?_WIN32

#ifdef USE_MKL
#ifndef MKL_Complex8
#ifdef __cplusplus
#define MKL_Complex8 std::complex<float>
#else // C99
#define MKL_Complex8 float _Complex
#endif // ?__cplusplus
#endif // !MKL_Complex8

#ifndef MKL_Complex16
#ifdef __cplusplus
#define MKL_Complex16 std::complex<double>
#else // C99
#define MKL_Complex16 double _Complex
#endif // ?__cplusplus
#endif // !MKL_Complex16
#endif // USE_MKL

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else // NVCC host compiler
#ifdef __cplusplus
#include <cmath>
#include <complex>
#else // C99
#include <math.h>
#include <complex.h>
#endif // ?__cplusplus
#endif // ?__INTEL_COMPILER

#ifdef USE_MKL
#include <mkl.h>
#endif // USE_MKL

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
#endif // ?__cplusplus

#ifdef _WIN32
#include <io.h>
#include <malloc.h>
#ifndef alloca
#define alloca _alloca
#endif // !alloca
#else // POSIX
#include <alloca.h>
#include <sys/time.h>
#endif // ?_WIN32

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
#endif // ?_WIN32

// defines

#ifndef EXTERN_C
#ifdef __cplusplus
#define EXTERN_C extern "C"
#else // C
#define EXTERN_C extern
#endif // ?__cplusplus
#else // EXTERN_C
#error EXTERN_C not definable externally
#endif // ?EXTERN_C

#ifdef VAR_UNUSED
#error VAR_UNUSED not definable externally
#endif // VAR_UNUSED

#ifdef _WIN32
#define VAR_UNUSED
#else // POSIX
#define VAR_UNUSED __attribute__ ((unused))
#endif // ?_WIN32

#endif // !DEFINES_HPP
