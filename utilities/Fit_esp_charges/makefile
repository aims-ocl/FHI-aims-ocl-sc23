.SUFFIXES: .inc .f .f90 .F
# Intel
FC=mpif90
FFLAGS = -O3
#-debug all -check all -warn -traceback
FLIB = -L$(MKL_HOME)/lib/intel64 -I$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

# gnu / gfortran
#FC=mpif90.openmpi
#FFLAGS = -O3 -ffree-line-length-none
#FLIB = -llapack -lblas -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -lblacsF77init-openmpi



OBJS = mpi_tasks.o Fit_esp_charges.o
Target=Fit_esp_charges
$(Target):  $(OBJS)
	$(FC) $(OBJS) $(FFLAGS) $(FLIB) -o $(Target) 

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

clean:
	-@rm $(Target) *.o *.mod
                           
                                                                                                                        
