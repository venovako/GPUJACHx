RM=rm -rf

ifndef NVCC
NVCC=nvcc
endif # !NVCC

ifndef GPU_ARCH
GPU_ARCH=sm_30
endif # !GPU_ARCH

ifndef CUBLAS
CUBLAS=2
endif # !CUBLAS

ifndef HDF5
HDF5=$(HOME)/hdf5
endif # !HDF5

HDF5_I=-I$(HDF5)/include
HDF5_L=-L$(HDF5)/lib
HDF5_l=-l hdf5_hl,hdf5

NVPPFLAGS=-I. -I../share $(HDF5_I) -DHZ_GSVD
NVCCFLAGS=-std c++11 -arch $(GPU_ARCH) -cudart shared -restrict -Xcompiler -ffp-contract=on,-march=native,-integrated-as
NVLDFLAGS=$(HDF5_L) $(HDF5_l)

ifdef NDEBUG
NVPPFLAGS+=-DNDEBUG
NVCCFLAGS+=-O$(NDEBUG) -Xptxas -v
else # debug
NVPPFLAGS+=-UNDEBUG
NVCCFLAGS+=-O0 -g -G -src-in-ptx -Xptxas -no-bb-merge,-v
endif # NDEBUG

NVPPFLAGS+=-UUSE_MULTI_GPU
NVPPFLAGS+=-UUSE_MKL

ifdef CUBLAS
NVPPFLAGS+=-DHAVE_CUBLAS=$(CUBLAS)
NVLDFLAGS+=-l cublas
else # !CUBLAS
NVPPFLAGS+=-UHAVE_CUBLAS
endif # CUBLAS

NVLDFLAGS+=-l z,dl,m

HEADERS=				\
	HZ.hpp				\
	HZ_L.hpp			\
	HZ_L2.hpp			\
	cuda_helper.hpp			\
	cuda_memory_helper.hpp		\
	defines.hpp			\
	device_code.hpp			\
	device_code_accumV.hpp 		\
	device_code_cdsort.hpp 		\
	device_code_cdsort_accumV.hpp	\
	device_code_cdsort_solveV.hpp	\
	device_code_common.hpp		\
	device_code_globals.hpp		\
	device_code_common_rotate.hpp	\
	device_code_common_Kepler.hpp	\
	device_code_common_Cholesky.hpp	\
	hdf5_helper.hpp			\
	my_utils.hpp

SOURCES=		\
	HZ_L.cu		\
	HZ_L2.cu	\
	cuda_helper.cu	\
	device_code.cu	\
	main.cu		\
	my_utils.cu

EXE=HZ1.exe

all: $(EXE)

$(EXE): $(HEADERS) $(SOURCES) GNUmakefile
	$(NVCC) $(NVPPFLAGS) $(NVCCFLAGS) $(SOURCES) -o $@ $(NVLDFLAGS)

clean:
	-$(RM) $(EXE)
	-$(RM) $(EXE).dSYM
