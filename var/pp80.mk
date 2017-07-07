RM=rm -f
FC=gfortran
ifdef NDEBUG
FOPTS=-std=f2008ts -O$(NDEBUG) -cpp -DNDEBUG -march=native -Wa,-q -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller -fopenmp -fstack-arrays -Wall -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
else
FOPTS=-std=f2008ts -Og -cpp -march=native -Wa,-q -fopenmp -fstack-arrays -fcheck=all -finit-local-zero -finit-real=snan -ffpe-trap=invalid,zero,overflow -Wall -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all # -finit-derived
endif
HDF5=$(HOME)/hdf5
H5INC=-I$(HDF5)/include
H5LIB=-Wl,-rpath $(HDF5)/lib -L$(HDF5)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
MKLIB=-Wl,-rpath ${MKLROOT}/lib -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#MKLIB=-Wl,-rpath ${MKLROOT}/lib -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LOPTS=

all: pp80 pp80_FMA

pp80: pp80J.exe pp80M.exe

pp80_FMA: pp80J_FMA.exe pp80M_FMA.exe

pp80J.exe: pp80.f90 pp80.mk
	$(FC) $(FOPTS) $(H5INC) pp80.f90 -o $@ $(LOPTS) $(H5LIB) $(MKLIB) 

pp80M.exe: pp80.f90 pp80.mk
	$(FC) $(FOPTS) $(H5INC) -DMAGMA pp80.f90 -o $@ $(LOPTS) $(H5LIB) $(MKLIB)

pp80J_FMA.exe: pp80.f90 pp80.mk
	$(FC) $(FOPTS) $(H5INC) -DUSE_FMA pp80.f90 -o $@ $(LOPTS) $(H5LIB) $(MKLIB) 

pp80M_FMA.exe: pp80.f90 pp80.mk
	$(FC) $(FOPTS) $(H5INC) -DMAGMA -DUSE_FMA pp80.f90 -o $@ $(LOPTS) $(H5LIB) $(MKLIB)

clean:
	-$(RM) pp80M.exe
	-$(RM) pp80J.exe
