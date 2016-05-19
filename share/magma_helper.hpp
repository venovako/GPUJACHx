#ifndef MAGMA_HELPER_HPP
#define MAGMA_HELPER_HPP

// Prevent MAGMA meddling with cuBLAS selection.
#include "cuda_helper.hpp"

#include <magma.h>

#ifndef MAGMA_CALL
#define MAGMA_CALL(call) {						\
    const magma_int_t err = (call);					\
    if (err) {                                                          \
      (void)fprintf(stderr, "MAGMA error %d @ %s(%d)!\n",               \
		    static_cast<int>(err), __FILE__, __LINE__);         \
      exit(EXIT_FAILURE);                                               \
    }									\
}
#else // MAGMA_CALL
#error MAGMA_CALL not definable externally
#endif // !MAGMA_CALL

#endif // !MAGMA_HELPER_HPP
