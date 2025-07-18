#include "cuda_helper.hpp"

#include "my_utils.hpp"

int configureGPUex(const int dev, const unsigned maxShMemB) throw()
{
  assert(dev >= 0);

  CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceScheduleSpin));
  CUDA_CALL(cudaSetDevice(dev));

  cudaDeviceProp cdp;
  CUDA_CALL(cudaGetDeviceProperties(&cdp, dev));
  const int dcc = cdp.major * 10 + cdp.minor;

  if (dcc < 20) {
    (void)snprintf(err_msg, err_msg_size, "CUDA Device %d Compute Capability %d.%d < 2.0", dev, cdp.major, cdp.minor);
    DIE(err_msg);
  }

  if (WARP_SZ != static_cast<unsigned>(cdp.warpSize)) {
    (void)snprintf(err_msg, err_msg_size, "CUDA Device %d has %d threads in a warp, must be %u", dev, cdp.warpSize, WARP_SZ);
    DIE(err_msg);
  }

  if (!cdp.unifiedAddressing) {
    (void)snprintf(err_msg, err_msg_size, "CUDA Device %d does not support unified addressing", dev);
    DIE(err_msg);
  }

  cudaFuncCache cacheConfig = cudaFuncCachePreferNone;
  if (maxShMemB <= 16384u) // 16 kB
    cacheConfig = cudaFuncCachePreferL1;
  else if (maxShMemB <= 32768u) // 32 kB
    cacheConfig = cudaFuncCachePreferEqual;
  else if (maxShMemB <= 49152u) // 48 kB
    cacheConfig = cudaFuncCachePreferShared;
  else { // > 48 kB
    (void)snprintf(err_msg, err_msg_size, "Maximum shared memory requested (%u B) > 48 kB", maxShMemB);
    WARN(err_msg);
  }
  CUDA_CALL(cudaDeviceSetCacheConfig(cacheConfig));
  /* deprecated:
  if (dcc >= 30)
    CUDA_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
  */
  return dcc;
}

int configureGPU(const int dev) throw()
{
  static const unsigned maxShMemB = 16384u; // 16 kB
  return configureGPUex(dev, maxShMemB);
}
