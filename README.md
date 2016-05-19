# GPUJACHx
The Jacobi-type (hyperbolic) Singular Value Decomposition in CUDA (for 1 GPU).

This software is a supplementary material for the paper
[doi:10.1137/140952429](http://dx.doi.org/10.1137/140952429 "A hierarchically blocked Jacobi SVD algorithm for single and multiple graphics processing units").
The preprint is available at [arXiv:1401.2720](http://arxiv.org/abs/1401.2720 "A hierarchically blocked Jacobi SVD algorithm for single and multiple graphics processing units").

Multi-GPU level is still under cleanup and will be added eventually.

Directories:
1     - Test program for 1 GPU,
share - Common and device code,
strat - Jacobi strategy tables.

Build strategy tables:
cd strat
./proc_all-lnx.sh         (Linux), or
./proc_all-win.bat        (Windows), or
./proc_all-mac.sh         (MacOS).

Build the test program:
cd 1
./mkLNX.sh  SM D          (Linux), or
./mkWIN.bat SM D          (Windows), or
./mkMAC.sh  SM D          (MacOS), where
SM = target GPU architecture (e.g., 30 or 35 or 37 for Kepler, 20 for Fermi),
D  = optimization level (usually 3).

Prerequisites for the test program:
HDF5 library (please specify where it is installed by editing share/Makefile.XXX).

Test run:
./x1.exe DEV SDY SNP ALG H5F H5G [H5R]
where
DEV = CUDA device number
SDY = path to strat.so (Linux), or strat.dll (Windows), or strat.dylib (MacOS)
SNP = BrentL or rowcyc or colcyc or cycwor or cycloc or mmstep strategy
ALG = algorithm number (see 1/HypJacL2.hpp)
H5F = input HDF5 file (see main.cu)
H5G = input/output HDF5 group
H5R = optional output HDF5 file
