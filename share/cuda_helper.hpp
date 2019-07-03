#ifndef CUDA_HELPER_HPP
#define CUDA_HELPER_HPP

#include "defines.hpp"

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
#endif // ?CUDA_CALL

#ifndef WARP_SZ
#define WARP_SZ 32u
#else // WARP_SZ
#error WARP_SZ not definable externally
#endif // ?WARP_SZ

EXTERN_C int configureGPU(const int dev) throw();
EXTERN_C int configureGPUex(const int dev, const unsigned maxShMemB) throw();

#endif // !CUDA_HELPER_HPP
