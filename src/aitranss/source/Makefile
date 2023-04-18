############################################################
#      a i t r a n s s : ab initio transport simulations
#     (c)  2003-2012   : alexej bagrets,  andreas arnold
#                        florian weigend, ferdinand evers
#     institute of nanotechnology (int) &
#     institut fuer theorie der kondensierten materie (tkm)
#     karlsruhe institute of technology (kit)
#
#     author:         alexej.bagrets <at> kit.edu
#     date:           nov 2009
#     last revision:  jan 2012
############################################################

# This is the template Makefile for the "aitranss" module
#
# A mandatory prerequisite to build the code is a Fortran90/95
# capable compiler and a compiled version of LAPACK and BLAS
# (for example, Intel's MKL). Please, adjust variables
# FC, LD, FLAGS and LIBS appropriately.
#
# A binary ("aitranss.x") built by the Makefile will go to
# the ./bin directory of the FHI-aims.

# REMARK :  Current release of the transport code, in
#           contrast to FHI-aims, does not yet based
#           on MPI. However, your are encouraged to use a
#           fortran compiler option(s), aka "-openmp" and
#           "-O2" for Intel's ifort, to build a multithreaded
#           code based on OpenMP directives.
#
#           According to our experience, a generated code
#           can be safety executed in parallel within a
#           single compute node with multiple processors,
#           and with a significant gain in computation time.

# version format:  MMDDYY  (well, usually ...)

# comment added : Mon Feb 20 10:18:51 CET 2012
# this is to be the first release version of "aitranss"
# to work together with FHI-aims >>>
# aitranss_VERSION = 022012
aitranss_VERSION = 071813

# bin directory ../../../bin, consistent with FHI-aims
bindir = ../../../bin

###  'inttp' Linux-cluster at Inst. for Nanotechnology(INT)/KIT 
###  intel fortran compiler, no OpenMP support
#FC=ifort
#LD=ifort
#FFLAGS= -132 -implicitnone -warn all -assume byterecl -O2
#LIBS= -L/opt/software/intel/mkl/8.1/lib/32/ -lmkl_lapack -lmkl -lguide

###  'intcool' Linux-cluster at Inst. for Nanotechnology(INT)/KIT 
###  intel fortran compiler, with OpenMP support (-openmp):
###  generates multithreaded code based on OpenMP directives
FC=ifort
LD=ifort
FFLAGS= -132 -implicitnone -warn all -openmp -assume byterecl -O2
LIBS= -L/data/software/intel/Compiler/11.0/083/mkl/lib/em64t -lmkl_lapack -lmkl -lmkl_intel_thread -lguide

## for intact: ifort compiler
#FC=ifort
#LD=ifort
#FFLAGS= -132 -implicitnone -warn all -openmp -assume byterecl -O3
#LIBS= -L/data/software/intel/lib/em64t -lmkl_lapack -lmkl -lpthread

### for intact: pgf90 compiler --> best choice
#FC=pgf90
#LD=pgf90
#FFLAGS= -Minform=warn -Mextend -O3
#LIBS= -L/data/software/pgi/linux86-64/6.1/lib -lacml

### for intact: pgf90 compiler --> best choice, openMP
#FC=pgf90
#LD=pgf90
#FFLAGS= -Minform=warn -Mextend -fast -fastsse -tp=amd64 -mp -Mconcur -O3
#LIBS= -L/data/software/pgi/linux86-64/6.1/lib -lacml_mp -mp -Mconcur

###  HP XC3000 Linux-Cluster / Quad-core Intel Xeon E5540 (Nehalem, 2.53GHz)
###  at Steinbuch Centre for Computing(SCC)/KIT
###  intel fortran compiler, with OpenMP support (-openmp):
###  generates multithreaded code based on OpenMP directives
#FC=ifort
#LD=ifort
#FFLAGS= -132 -implicitnone -warn all -openmp -assume byterecl -O2 -I${MKLINCLUDE}
#LIBS= -L${MKLPATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm

modules = \
    globalvars.o \
    info.o \
    ptable.o \
    tools.o \
    rmatrix.o \
    eigenvalue_solver.o \
    cmatrix.o \
    externals.o \
    selfenergy.o \
    math.o \
    overlap.o \
    domainwall.o \
    hamiltonian.o \
    popanalysis.o \
    exthamiltonian.o \
    densmat.o \
    densmateq.o \
    int_out_core_states.o \
    hubbardu.o \
    spectralfunc.o \
    molsubsystem.o \
    ldos.o \
    ldos_rpoint.o \
    scatt_wave_func.o \
    trans.o \
    rsigma.o \
    main.o
    
.suffixes: .f.o
.f.o:
	${FC} ${FFLAGS} -c $<

aitranss.${aitranss_VERSION}.x: ${modules}	
	./generate_version_stamp.pl ${aitranss_VERSION} > output_version_stamp.f 
	${FC} ${FFLAGS} -c output_version_stamp.f
	${LD} output_version_stamp.o ${modules} ${LIBS} -o ${bindir}/$@ 

cleanmod:	
	rm *.mod

clean:	
	rm *.mod *.o
