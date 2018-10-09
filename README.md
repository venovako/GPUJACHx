# GPUJACHx
The Jacobi-type (hyperbolic) Singular Value Decomposition in CUDA (for 1 GPU).

This software is a supplementary material for the paper
[doi:10.1137/140952429](http://dx.doi.org/10.1137/140952429 "A hierarchically blocked Jacobi SVD algorithm for single and multiple graphics processing units").

The preprint is available at [arXiv:1401.2720](http://arxiv.org/abs/1401.2720 "A hierarchically blocked Jacobi SVD algorithm for single and multiple graphics processing units").

Multi-GPU level is still under cleanup and will be added eventually.

Directories:
* 1     - Test program for 1 GPU,
* share - Common and device code,
* strat - Jacobi strategy tables.

Build strategy tables:
```bash
cd strat
./proc_all-lnx.sh         (Linux), or
./proc_all-win.bat        (Windows), or
./proc_all-mac.sh         (MacOS).
```

Build the test program:
```bash
cd 1
./mkLNX.sh  SM D          (Linux), or
./mkWIN.bat SM D          (Windows), or
./mkMAC.sh  SM D          (MacOS), where
```
* SM = target GPU architecture (e.g., 30 or 35 or 37 for Kepler, 20 for Fermi),
* D  = optimization level (usually 3).

Prerequisites for the test program:
* HDF5 library (please specify where it is installed by editing share/Makefile.XXX).

Test run:
```bash
./x1.exe DEV SDY SNP ALG H5F H5G [H5R]
```
where
* DEV = CUDA device number
* SDY = path to strat.so (Linux), or strat.dll (Windows), or strat.dylib (MacOS)
* SNP = BrentL or rowcyc or colcyc or cycwor or cycloc or mmstep strategy
* ALG = algorithm number (see 1/HypJacL2.hpp)
* H5F = input HDF5 file (see main.cu)
* H5G = input/output HDF5 group
* H5R = optional output HDF5 file

## Input file format

Here, an example of the input files can be found:
http://euridika.math.hr:1846/Jacobi/data/

In essence, an HDF5 file can contain mulitple input matrices and the associated metadata, with each input in its own HDF5 group within the file.

Here is a Matlab code template to create an input dataset:
```Matlab
% G: input double precision M x N real matrix to compute (H)SVD of.
%    size(G) = [LDA N]
% LDA >= M >= N >= NPLUS >= 0;
%    for 1 GPU: (LDA = M = N) mod 64 == 0,
%    for 4 GPUs: (LDA = M = N) mod 256 == 0,
%    N = NPLUS for SVD (J = I),
%    N > NPLUS for HSVD (NPLUS positive (+1) and N-NPLUS negative (-1) signs in J).
IDADIM=int32([LDA; M; N; NPLUS]);
% Substitute 'file' and 'group' with appropriate names.
h5create('file.h5','/group/IDADIM',[4],'Datatype','int32')
h5write('file.h5','/group/IDADIM',IDADIM)
h5create('file.h5','/group/G',size(G))
h5write('file.h5','/group/G',G)
```
For Scilab it would be similar:
```Scilab
h5=h5open("file.h5","w")
h5group(h5,"/group")
IDADIM=int32([LDA;M;N;NPLUS])
h5write(h5,"/group/IDADIM",IDADIM)
h5write(h5,"/group/G",G)
h5close(h5)
```
In `file.h5` there will be a group named `group` created, with two datasets: `IDADIM` is an integer vector containing the dimensions of the matrix `G`, which has `M` rows and `N` columns, and is stored in *Fortran* array order (column by column), with the leading dimension `LDA`.

For now, please, make sure that `group` name is in fact the number of columns of `G`.
E.g., if the input matrix dataset has `2048` columns, the HDF5 group in which it is contained should be named `2048`.

Hint: if the file is viewed with, e.g., [HDFView](https://www.hdfgroup.org/downloads/hdfview/ "A free HDF5 viewer/editor."), `G` will look as if having `N` rows and `M` columns, which is perfectly normal, since the Fortran data look transposed in C, and vice versa.
Matlab and Scilab store and read data in the correct format for GPUJACHx.

For simplicity, it is expected that `LDA=M=N`, and that all dimensions are a multiple of `64`.  That is not a big issue, though, since a matrix can always be "bordered".

The bordering can be done in, e.g., Matlab or Scilab, as follows (`I` is an identity matrix):

|   N | 64-(N%64) |
| --- | --------- |
| `G` |       `0` |
| `0` |       `I` |

The output, if requested, is stored in the same-named `group` with `G` being the matrix of the left singular vectors, `V` the matrix of the right singular vectors (non-transposed), and `D` are the eigenvalues of `trans(G)*J*G`, i.e., the squares of the singular values (for the "ordinary" SVD, with `J=I`).

## Running

Note that the strategies do not support all possible input dimensions!
Please have a look at the strategies' headers in `strat` subdirectory.
It is expected that the fastest execution will be with cycwor.

Example, for the full SVD with the block oriented variant and output in `output.h5`:
```bash
./x1.exe 0 ../strat/strat.dylib cycwor 9 input.h5 group output.h5
```

This work has been supported in part by Croatian Science Foundation under the project IP--2014--09--3670.
