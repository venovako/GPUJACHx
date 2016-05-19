#ifndef HDF5_HELPER_HPP
#define HDF5_HELPER_HPP

#include "defines.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

#ifndef HYPJAC_IDADIM_SIZE
#define HYPJAC_IDADIM_SIZE static_cast<size_t>(4u)
#else // HYPJAC_IDADIM_SIZE
#error HYPJAC_IDADIM_SIZE not definable externally
#endif // !HYPJAC_IDADIM_SIZE

#ifndef HDF5_CALL
#define HDF5_CALL(call) {						\
    const herr_t err = static_cast<herr_t>(call);			\
    if (static_cast<herr_t>(0) > err) {					\
      (void)fprintf(stderr, "HDF5 error %d @ %s(%d)!\n",		\
		    static_cast<int>(err), __FILE__, __LINE__);         \
      exit(EXIT_FAILURE);                                               \
    }									\
  }
#else // HDF5_CALL
#error HDF5_CALL not definable externally
#endif // !HDF5_CALL

#endif // !HDF5_HELPER_HPP
