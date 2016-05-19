#include "cuda_helper.hpp"

#include "my_utils.hpp"

int configureGPU(const int dev) throw()
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

  CUDA_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
  if (dcc >= 30)
    CUDA_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));

  return dcc;
}
