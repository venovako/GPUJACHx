// main.cu: test driver.

#include "HZ.hpp"
#include "HZ_L2.hpp"

#include "cuda_memory_helper.hpp"
#include "hdf5_helper.hpp"
#include "my_utils.hpp"

struct CmdArgs {
  char *exe;
  char *dev;
  char *sdy;
  char *snp;
  char *alg;
  char *h5f;
  char *h5g;
  char *h5r;

  CmdArgs(const int argc, char *const argv[]) throw()
    : exe(static_cast<char*>(NULL)),
      dev(static_cast<char*>(NULL)),
      sdy(static_cast<char*>(NULL)),
      snp(static_cast<char*>(NULL)),
      alg(static_cast<char*>(NULL)),
      h5f(static_cast<char*>(NULL)),
      h5g(static_cast<char*>(NULL)),
      h5r(static_cast<char*>(NULL))
  {
    if (1 > argc)
      DIE("argc <= 0");
    if (!argv)
      DIE("NULL argv");

    if ((7 > argc) || (8 < argc)) {
      (void)fprintf(stderr, "%s DEV SDY SNP ALG H5F H5G [H5R]\n", argv[0]);
      exit(EXIT_FAILURE);
    }

    exe = argv[0];
    dev = argv[1];
    sdy = argv[2];
    snp = argv[3];
    alg = argv[4];
    h5f = argv[5];
    h5g = argv[6];
    h5r = ((8 == argc) ? argv[7] : static_cast<char*>(NULL));
  }
};

int main(int argc, char *argv[])
{
  int ret = EXIT_SUCCESS;

  unsigned idadim[HZ_IDADIM_SIZE] = { 0u };
  hid_t fid = static_cast<hid_t>(-1), gid = static_cast<hid_t>(-1);
  size_t ldA = static_cast<size_t>(0u), m = static_cast<size_t>(0u), n = static_cast<size_t>(0u);

  CmdArgs ca(argc, argv);

  const unsigned n0 = (HZ_L1_NCOLB << 1u);
  const unsigned n1 = (static_cast<unsigned>(atoi(ca.h5g)) + HZ_L1_NCOLB - 1u) / HZ_L1_NCOLB;
  init_strats(ca.sdy, ca.snp, n0, n1);

  const int dev = atoi(ca.dev);
  const int dcc = configureGPU(dev);
  if (dcc < 30) {
    (void)snprintf(err_msg, err_msg_size, "Device %d has CC %d < 30", dev, dcc);
    DIE(err_msg);
  }

  const unsigned routine = static_cast<unsigned>(atoi(ca.alg));

  if (!fexist(ca.h5f)) {
    (void)snprintf(err_msg, err_msg_size, "non-existent H5F(%s)", ca.h5f);
    DIE(err_msg);
  }

  HDF5_CALL(H5open());

  HDF5_CALL((fid = H5Fopen(ca.h5f, H5F_ACC_RDONLY, H5P_DEFAULT)));
  HDF5_CALL((gid = H5Gopen2(fid, ca.h5g, H5P_DEFAULT)));

  HDF5_CALL(H5LTread_dataset_int(gid, "IDADIM", reinterpret_cast<int*>(idadim)));

  unsigned
    ldhF = idadim[0],
    ldhG = idadim[0],
    ldhV = 0u,
    nrow = idadim[0], //idadim[1]
    ncol = idadim[0]; //idadim[2]

  m = static_cast<size_t>(nrow);
  n = static_cast<size_t>(ncol);

  ldA = static_cast<size_t>(ldhF);
  double *const hF = allocHostMtx<double>(ldA, m, n, true);
  SYSP_CALL(hF);
  ldhF = static_cast<unsigned>(ldA);

  HDF5_CALL(H5LTread_dataset_double(gid, "F", hF));

  ldA = static_cast<size_t>(ldhG);
  double *const hG = allocHostMtx<double>(ldA, m, n, true);
  SYSP_CALL(hG);
  ldhG = static_cast<unsigned>(ldA);

  HDF5_CALL(H5LTread_dataset_double(gid, "G", hG));

  HDF5_CALL(H5Gclose(gid));
  HDF5_CALL(H5Fclose(fid));

  double *hV = static_cast<double*>(NULL);
  ldhV = ((ldhF <= ldhG) ? ldhF : ldhG);
  ldA = static_cast<size_t>(ldhV);
  hV = allocHostMtx<double>(ldA, n, n, true);
  SYSP_CALL(hV);
  ldhV = static_cast<unsigned>(ldA);

  double *const hS = allocHostVec<double>(n);
  SYSP_CALL(hS);
  double *const hH = allocHostVec<double>(n);
  SYSP_CALL(hH);
  double *const hK = allocHostVec<double>(n);
  SYSP_CALL(hK);

  unsigned glbSwp = 0u;
  unsigned Long glb_s = MkLong(0u), glb_b = MkLong(0u);
  double timing[4] = { -0.0, -0.0, -0.0, -0.0 };
  ret = HZ_L2(routine, nrow, ncol, hF, ldhF, hG, ldhG, hV, ldhV, hS, hH, hK, &glbSwp, &glb_s, &glb_b, timing);

  if (ret)
    (void)fprintf(stderr, "%s: error %d\n", ca.exe, ret);
  else {
    (void)printf("GLB_ROT_S(%20" FmtLong "u), GLB_ROT_B(%20" FmtLong "u)\n", glb_s, glb_b);
    (void)printf("%#12.6f s %2u sweeps\n", *timing, glbSwp);
  }

  if (ca.h5r) {
    HDF5_CALL(fid = fexist(ca.h5r) ?
              H5Fopen(ca.h5r, H5F_ACC_RDWR, H5P_DEFAULT) :
              H5Fcreate(ca.h5r, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_CALL(gid = H5Gcreate2(fid, ca.h5g, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    idadim[0] = ldhG;
    hsize_t dims[2] = { static_cast<hsize_t>(HZ_IDADIM_SIZE), static_cast<hsize_t>(1u) };
    HDF5_CALL(H5LTmake_dataset_int(gid, "IDADIM", 1, dims, reinterpret_cast<int*>(idadim)));

    int istats[4] = { 0, static_cast<int>(STRAT0), static_cast<int>(STRAT1), dev };
    istats[0] = ((ret < 0) ? ret : static_cast<int>(glbSwp));
    dims[0] = static_cast<hsize_t>(4u);
    dims[1] = static_cast<hsize_t>(1u);
    HDF5_CALL(H5LTmake_dataset_int(gid, "ISTATS", 1, dims, istats));

    dims[0] = static_cast<hsize_t>(4u);
    dims[1] = static_cast<hsize_t>(1u);
    HDF5_CALL(H5LTmake_dataset_double(gid, "TIMING", 1, dims, timing));

    if (hF) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(ldhF);
      HDF5_CALL(H5LTmake_dataset_double(gid, "F", 2, dims, hF));
    }

    if (hG) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(ldhG);
      HDF5_CALL(H5LTmake_dataset_double(gid, "G", 2, dims, hG));
    }

    if (hV) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(ldhV);
      HDF5_CALL(H5LTmake_dataset_double(gid, "V", 2, dims, hV));
    }

    if (hS) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(1u);
      HDF5_CALL(H5LTmake_dataset_double(gid, "SIGMA", 1, dims, hS));
    }

    if (hH) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(1u);
      HDF5_CALL(H5LTmake_dataset_double(gid, "H", 1, dims, hH));
    }

    if (hK) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(1u);
      HDF5_CALL(H5LTmake_dataset_double(gid, "K", 1, dims, hK));
    }

    HDF5_CALL(H5Gclose(gid));
    HDF5_CALL(H5Fclose(fid));
  }

  HDF5_CALL(H5close());

  if (hK)
    CUDA_CALL(cudaFreeHost(hK));
  if (hH)
    CUDA_CALL(cudaFreeHost(hH));
  if (hS)
    CUDA_CALL(cudaFreeHost(hS));
  if (hV)
    CUDA_CALL(cudaFreeHost(hV));
  if (hG)
    CUDA_CALL(cudaFreeHost(hG));
  if (hF)
    CUDA_CALL(cudaFreeHost(hF));

  // for profiling
  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaDeviceReset());

  return ret;
}
