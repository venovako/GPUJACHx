// main.cu: test driver.

#include "HypJac.hpp"
#include "HypJacL2.hpp"

#include "cuda_memory_helper.hpp"
#include "hdf5_helper.hpp"
#include "my_utils.hpp"

struct CmdArgs {
  char *exe;
  char *dev;
  char *sdy;
  char *snp0;
  char *snp1;
  char *alg;
  char *h5f;
  char *h5g;
  char *h5r;

  CmdArgs(const int argc, char *const argv[]) throw()
    : exe(static_cast<char*>(NULL)),
      dev(static_cast<char*>(NULL)),
      sdy(static_cast<char*>(NULL)),
      snp0(static_cast<char*>(NULL)),
      snp1(static_cast<char*>(NULL)),
      alg(static_cast<char*>(NULL)),
      h5f(static_cast<char*>(NULL)),
      h5g(static_cast<char*>(NULL)),
      h5r(static_cast<char*>(NULL))
  {
    if (1 > argc)
      DIE("argc <= 0");
    if (!argv)
      DIE("NULL argv");

    if ((8 > argc) || (9 < argc)) {
      (void)fprintf(stderr, "%s DEV SDY SNP0 SNP1 ALG H5F H5G [H5R]\n", argv[0]);
      exit(EXIT_FAILURE);
    }

    exe = argv[0];
    dev = argv[1];
    sdy = argv[2];
    snp0 = argv[3];
    snp1 = argv[4];
    alg = argv[5];
    h5f = argv[6];
    h5g = argv[7];
    h5r = ((9 == argc) ? argv[8] : static_cast<char*>(NULL));
  }
};

int main(int argc, char *argv[])
{
  int ret = EXIT_SUCCESS;

  unsigned idadim[HYPJAC_IDADIM_SIZE] = { 0u };
  hid_t fid = static_cast<hid_t>(-1), gid = static_cast<hid_t>(-1);
  size_t ldA = static_cast<size_t>(0u), m = static_cast<size_t>(0u), n = static_cast<size_t>(0u);

  CmdArgs ca(argc, argv);

  const unsigned n0 = (HYPJACL1_NCOLB << 1u);
  const unsigned n1 = (static_cast<unsigned>(atoi(ca.h5g)) + HYPJACL1_NCOLB - 1u) / HYPJACL1_NCOLB;
  init_strats(ca.sdy, ca.snp0, n0, ca.snp1, n1);

  const int dev = atoi(ca.dev);
#ifndef NDEBUG
  const int dcc =
#else // NDEBUG
  (void)
#endif // !NDEBUG
    configureGPU(dev);
#ifndef NDEBUG
  (void)fprintf(stdout, "Device %d [CC %d] initialized.\n", dev, dcc);
  (void)fflush(stdout);
#endif // !NDEBUG

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
    ldhG = idadim[0],
    nrow = idadim[1],
    ncol = idadim[2],
    nplus = idadim[3];

  m = static_cast<size_t>(nrow);
  n = static_cast<size_t>(ncol);

  ldA = static_cast<size_t>(ldhG);
  double *const hG = allocHostMtx<double>(ldA, m, n, true);
  SYSP_CALL(hG);
  ldhG = static_cast<unsigned>(ldA);

  HDF5_CALL(H5LTread_dataset_double(gid, "G", hG));

  HDF5_CALL(H5Gclose(gid));
  HDF5_CALL(H5Fclose(fid));

  double *hV;
  unsigned ldhV;
  if (routine & HYPJAC_FULL_SVD) {
    hV = allocHostMtx<double>(ldA, m, n, true);
    SYSP_CALL(hV);
    ldhV = static_cast<unsigned>(ldA);
  }
  else {
    hV = static_cast<double*>(NULL);
    ldhV = 0u;
  }

  double *const hD = allocHostVec<double>(n);
  SYSP_CALL(hD);

  unsigned glbSwp = 0u;
  unsigned Long glb_s = MkLong(0u), glb_b = MkLong(0u);
  double timing[4] = { -0.0, -0.0, -0.0, -0.0 };
  ret = hypJacL2(routine, nrow, ncol, nplus, hG, ldhG, hV, ldhV, hD, &glbSwp, &glb_s, &glb_b, timing);

  if (ret)
    (void)fprintf(stderr, "%s: error %d\n", ca.exe, ret);
  else {
    (void)fprintf(stdout, "GLB_ROT_S(%20" FmtLong "u), GLB_ROT_B(%20" FmtLong "u)\n", glb_s, glb_b);
    (void)fflush(stdout);
    (void)fprintf(stdout, "%#12.6f s %2u sweeps\n", *timing, glbSwp);
    (void)fflush(stdout);
  }

  if (ca.h5r) {
    HDF5_CALL(fid = fexist(ca.h5r) ?
              H5Fopen(ca.h5r, H5F_ACC_RDWR, H5P_DEFAULT) :
              H5Fcreate(ca.h5r, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_CALL(gid = H5Gcreate2(fid, ca.h5g, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    idadim[0] = ldhG;
    hsize_t dims[2] = { static_cast<hsize_t>(HYPJAC_IDADIM_SIZE), static_cast<hsize_t>(1u) };
    HDF5_CALL(H5LTmake_dataset_int(gid, "IDADIM", 1, dims, reinterpret_cast<int*>(idadim)));

    int istats[4] = { 0, static_cast<int>(STRAT0), static_cast<int>(STRAT1), dev };
    istats[0] = ((ret < 0) ? ret : static_cast<int>(glbSwp));
    dims[0] = static_cast<hsize_t>(4u);
    dims[1] = static_cast<hsize_t>(1u);
    HDF5_CALL(H5LTmake_dataset_int(gid, "ISTATS", 1, dims, istats));

    dims[0] = static_cast<hsize_t>(4u);
    dims[1] = static_cast<hsize_t>(1u);
    HDF5_CALL(H5LTmake_dataset_double(gid, "TIMING", 1, dims, timing));
    
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

    if (hD) {
      dims[0] = static_cast<hsize_t>(ncol);
      dims[1] = static_cast<hsize_t>(1u);
      HDF5_CALL(H5LTmake_dataset_double(gid, "D", 1, dims, hD));
    }

    HDF5_CALL(H5Gclose(gid));
    HDF5_CALL(H5Fclose(fid));
  }

  HDF5_CALL(H5close());

  if (hD)
    CUDA_CALL(cudaFreeHost(hD));
  if (hV)
    CUDA_CALL(cudaFreeHost(hV));
  if (hG)
    CUDA_CALL(cudaFreeHost(hG));

  // for profiling
  CUDA_CALL(cudaDeviceSynchronize());
  CUDA_CALL(cudaDeviceReset());

  return ret;
}
