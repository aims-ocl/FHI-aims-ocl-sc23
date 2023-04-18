.SUFFIXES: .inc .f .f90 .F
F90=mpiifort
#F90=ifort
FC=$(F90)
FFLAGS = -O2 -traceback
FLIB = -L$(MKL_HOME)/lib/intel64 -I$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#FLIB =  -L/opt/intel/Compiler/11.0/084/mkl/lib/em64t -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
OBJS = parameter_reader.o MODOS_FHI_aims_mpi.o
Target=MODOS_FHI_aims_mpi
$(Target):  $(OBJS)
	$(FC) $(OBJS) $(FFLAGS) $(FLIB) -o $(Target) 

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

clean:
	-@rm $(Target) *.o *.mod
