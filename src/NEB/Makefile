VERSION          = 111907 
# serial version: thlc's anyway
COMPILER         = ifort
LAPACK           = -L/opt/intel/mkl/8.1/lib/32/ -lmkl_lapack -lmkl
COMPILER_OPTIONS = -traceback -g -check bounds -check uninit
MPIROUTINES      = mpi_dummy

# parallel version: p5's in Garching for starters
# COMPILER         = mpxlf90_r
# LAPACK           = -L/afs/rzg/@sys/lib -llapack-essl -lessl -lpesslsmp -lblacssmp
# COMPILER_OPTIONS = 
# MPIROUTINES      = mpi_tasks


SUBROUTINES = NEB_force_updates.o \
	      NEB_geometry.o      


NEB: NEB.f90 $(MPIROUTINES).o $(SUBROUTINES)
	cat NEB.f90 | sed 's/mpi_routines/'$(MPIROUTINES)'/g' > NEB_temp.f90; \
	$(COMPILER) $(COMPILER_OPTIONS) -o NEB.x NEB_temp.f90 $(SUBROUTINES) $(LAPACK) $(MPIROUTINES).o ; \
	rm NEB_temp.f90


$(MPIROUTINES).o: $(MPIROUTINES).f90
	$(COMPILER) $(COMPILER_OPTIONS) -c $(MPIROUTINES).f90

NEB_force_updates.o: NEB_force_updates.f90
	$(COMPILER) $(COMPILER_OPTIONS) -c NEB_force_updates.f90

NEB_geometry.o: NEB_geometry.f90
	$(COMPILER) $(COMPILER_OPTIONS) -c NEB_geometry.f90	

clean:
	rm NEB.x *.o *.mod
