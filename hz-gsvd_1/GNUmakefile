RM=rm -rf

ifndef NVCC
NVCC=nvcc
endif # !NVCC

ifndef GPU_ARCH
GPU_ARCH=sm_30
endif # !GPU_ARCH

ifndef HDF5
HDF5=$(HOME)/hdf5
endif # !HDF5

ifndef CVG
CVG=0
endif # !CVG

HDF5_I=-I$(HDF5)/include
HDF5_L=-L$(HDF5)/lib

# e.g., mpi,z
ifdef HDF5_FULL
HDF5_l=-l hdf5_hl,hdf5,$(HDF5_FULL)
else # minimal HDF5
HDF5_l=-l hdf5_hl,hdf5
endif # HDF5_FULL

NVPPFLAGS=-I. $(HDF5_I) -DHZ_GSVD -DCVG=$(CVG)
NVCCFLAGS=-std c++11 -arch $(GPU_ARCH) -cudart shared -restrict -Xcompiler -march=native
NVLDFLAGS=$(HDF5_L) $(HDF5_l)

ifdef NDEBUG
NVPPFLAGS+=-DNDEBUG
NVCCFLAGS+=-O$(NDEBUG) -Xptxas -v
else # debug
NVPPFLAGS+=-UNDEBUG
NVCCFLAGS+=-O0 -g -G -src-in-ptx -Xptxas -no-bb-merge,-v
endif # NDEBUG

# e.g. -ffp-contract=on,-integrated-as
ifdef HOST_FLAGS
NVCCFLAGS+=-Xcompiler $(HOST_FLAGS)
endif # HOST_FLAGS

NVPPFLAGS+=-UUSE_MULTI_THR
NVLDFLAGS+=-l dl,m

HEADERS=				\
	HZ.hpp				\
	HZ_L.hpp			\
	HZ_L2.hpp			\
	cuda_helper.hpp			\
	cuda_memory_helper.hpp		\
	defines.hpp			\
	device_code.hpp			\
	device_code_accumV.hpp 		\
	device_code_cdsort_$(CVG).hpp	\
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
