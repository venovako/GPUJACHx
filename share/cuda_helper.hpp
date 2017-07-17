#ifndef CUDA_HELPER_HPP
#define CUDA_HELPER_HPP

#include "defines.hpp"

#ifdef HAVE_CUBLAS
#if (HAVE_CUBLAS == 2)
#include <cublas_v2.h>
// prevent loading cublas.h
#ifndef CUBLAS_H_
#define CUBLAS_H_
#endif // !CUBLAS_H_
#elif (HAVE_CUBLAS == 1)
#include <cublas.h>
#else
#error unsupported cuBLAS version
#endif
#else // !HAVE_CUBLAS
#include <cuComplex.h>
#endif // HAVE_CUBLAS

#include <cuda_runtime.h>
#include <math_constants.h>

#ifndef CUDA_CALL
#define CUDA_CALL(call) {						\
    const cudaError_t err = (call);					\
    if (cudaSuccess != err) {						\
      (void)fprintf(stderr, "CUDA runtime error %d [%s] @ %s(%d)!\n",   \
		    static_cast<int>(err), cudaGetErrorString(err),     \
                    __FILE__, __LINE__);                                \
      exit(EXIT_FAILURE);                                               \
    }									\
}
#else // CUDA_CALL
#error CUDA_CALL not definable externally
#endif // !CUDA_CALL

#ifdef HAVE_CUBLAS
#ifndef CUBLAS_CALL
#define CUBLAS_CALL(call) {						\
    const cublasStatus_t err = (call);					\
    if (CUBLAS_STATUS_SUCCESS != err) {					\
      (void)fprintf(stderr, "CUBLAS runtime error %d @ %s(%d)!\n",      \
		    static_cast<int>(err), __FILE__, __LINE__);         \
      exit(EXIT_FAILURE);                                               \
    }									\
}
#else // CUBLAS_CALL
#error CUBLAS_CALL not definable externally
#endif // !CUBLAS_CALL
#endif // HAVE_CUBLAS

#ifndef WARP_SZ
#define WARP_SZ 32u
#else // WARP_SZ
#error WARP_SZ not definable externally
#endif // !WARP_SZ

EXTERN_C int
configureGPU(const int dev) throw();

EXTERN_C int
configureGPUex(const int dev, const unsigned maxShMemB) throw();

#endif // !CUDA_HELPER_HPP
