#
#  Makefile for FHI-aims, VB / aims team, 2004-2018
#
#  This Makefile works with GNU make. NOTE that on some proprietary platforms,
#  different versions of the ancient Unix make may not work - in that case, please
#  install GNU make - if it's not already there as "gmake".
#
#  MANDATORY REQUIREMENTS:
#
#  To build FHI-aims on your own architecture, you will need to adapt at 
#  least the paths given below, possibly more. Minimal third-party requirements 
#  are: 
#
#  - A working Fortran compiler (e.g., a recent version of Intel's ifort; 
#    alternatively, the free g95 compiler, which however generates significantly 
#    slower code).
#
#  - a Lapack implementation (e.g., www.netlib.org/lapack)
#  - a working Basic Linear Algebra Subroutine (BLAS) implementation for your
#    system.
#
#  BLAS and Lapack routines alike are, e.g., distributed with Intel's Math Kernel Library (mkl). 
#  Alternatively, outstanding BLAS implementations exist and can be built using the free ATLAS
#  or Goto libraries.
#
#  OPTIONAL REQUIREMENTS:
#
#  - Simple parallelism is supported through MPI; an appropriate free implementation is
#    mpich2, which must be set up on your system. Most commercial supercomputer architectures
#    ship their own MPI implementation, which should be documented locally.
#
#  - In addition, ScaLapack support (parallel eigenvalue solver) can be enabled IF a version 
#    of ScaLapack and a version of the Basic Linear Algebra Communication subroutines (BLACS) 
#    exist on your system. It is not necessary to have ScaLapack support in order to run in 
#    parallel. However, if not already available, ScaLapack and BLACS can be built with very
#    reasonable effort from the versions distributed at www.netlib.org.
#

#############################################################################################
#
#  Following: Settings to be customized by the user - three examples provided.
#
#  In general, you will need to:
#
#  1) Adjust the AIMSVERSION version stamp
#
#  2) Adjust the Compiler / architecture specific variable list immediately following
#
#  3) Possibly, adapt one of the make targets (further below) for your needs.
#
#  All other settings (version stamp, modules, subroutines, generic make targets) are included 
#  through a file Makefile.backend. With the exception of the version stamp, these should not
#  usually concern you.
#
#  Note that the very end of Makefile.backend lists some additional make targets that are 
#  not connected to the main aims binary at all, but instead independent utilities
#  for postprocessing (vibrations, phonons).
#
#############################################################################################


# MANDATORY Compiler- / architecture-specific parameters:
# YOU MUST ADJUST ALL THESE SETTINGS TO MATCH YOUR SPECIFIC COMPUTER PLATFORM.
# Uncomment ONLY ONE of the example lists below to 
#
# FC           : Fortran compiler
# FFLAGS       : Compiler flags for f77-style fixed-format subroutines
# F90FLAGS     : Compiler flags for f90-style free-format subroutines
# ARCHITECTURE : Allows to add any absolutely unavoidable architecture-specific choices at compile time
#                by changing a few simple options. Presently, some allowed choices are:
#                * BlueGene_P
#                * Alpha
#                * PGI (for the Portland Group Compiler)
#                * Generic (default which should work for all other architectures).
#                Check Makefile.backend itself for a full list of possibilities.
#                Any other choice (including no choice at all) will always fall back to the
#                "Generic" choice. See the interim comment below for more information.
# LAPACKBLAS   : Compile path and libraries that hold Lapack and BLAS subroutines
# USE_MPI      : Either "yes" or "no", depending on whether MPI is used or not.
#                The value of this variable effects the defaults of several internal variables, including MPIF .
# MPIFC        : When compiling WITH    MPI    support, this is the MPI compiler.
# SCALAPACK    : When compiling WITH Scalapack support, these are the ScaLapack and BLACS Libraries.

# There are also some optional flags that need not be specified but that can be helpful.
# If compiling read_control.f90 takes too much time, set the
#
# F90MINFLAGS  : keyword that specifies flags without optimization for faster compilation time
#                for files that are never a speed bottleneck at runtime. Example: read_control.f90 
#                with IBM compilers.
#
# F90REDUCEDFLAGS  : Keyword that specifies flags with reduced optimization for files that cause
#                    a compiler crash. These are cases that are not within the control of FHI-aims.
#                    It can happen that perfectly legal Fortran code crashes a specific compiler
#                    version. In this case, having a workaround helps.
#                    In addition to setting F90REDUCEDFLAGS, the source file in question must be 
#                    added by hand to Makefile.backend.
#                    Example: scs_cfdm.f90 with INTEL ifort compiler 13.1.
#                    If you need to work around a compiler problem, search for this file in 
#                    Makefile.backend and do exactly the same thing for the file that is problematic
#                    for you.
#
# PRINT_STACKTRACES  : if set to "yes", this enables the output of human readable stacktraces when
#                      FHI-aims terminates abnormally, requires a Intel Fortran compiler
#
# Check the list below, or Makefile.backend itself, for more settings that are optionally possible.

# NOTE on the ARCHITECTURE keyword above:
#
#                The ARCHITECTURE tag allows to specify compile-time choices for specific architectures.
#                Please be very careful with ever adding anything architecture-specific to the code.
#                As a rule, the same FHI-aims source code should normally run on any platform, and
#                this has worked very well in the past. There are a very few cases where a specific
#                Fortran compiler will really not understand an essential construct, or where an
#                architecture-specific hack is so beneficial that we prefer it anyway (BlueGene) -
#                but again, in principle this should be the absolute exception, and all changes
#                should be kept as simple as possible.

# NOTE: Consult your compiler manual to find out whether any optimization flags beyond these listed 
#       below will be helpful for YOUR specific compiler and computer system. If so, they should be 
#       included in FFLAGS and F90FLAGS, respectively.
#
#
 FC = mpif90
#FC = swafort
# FFLAGS = -O3 -ffree-line-length-none -mieee
# FFLAGS = -O3 -g -ffree-line-length-none -mieee -fcheck=all -fbacktrace
 FFLAGS = -O3 -g -ffree-line-length-none -mieee 
 F90FLAGS = $(FFLAGS)
 ARCH = arch_generic.f90
 MPIFC = $(FC)
 USE_MPI      = yes
#------------begin xMath_20201223-----------------------------------
my_mkl = /home/export/base/shisuan/swyjs/lxh/lib/xMath_20201223_g
LAPACKBLAS = -mhybrid -static  -L$(my_mkl) -lswlapack    -L$(my_mkl) -lswblas
SCALAPACK = -L$(my_mkl) -lswscalapack
#------------end xMath 20201223-----------------------------------

SWCL_LIB = -L/home/export/base/shisuan/swyjs/cq/swcl/install/lib/swcl -lswclapi

#------------begin no_opt-----------------------------------
#my_mkl = /home/export/online1/mdt00/shisuan/swict/2020_wuxi/2020_1008_mathlib/01_no_opt/mathlib/lapack-3.8.0
#LAPACKBLAS = -mhybrid -static  -L$(my_mkl) -llapack   -L$(my_mkl) -lrefblas
#SCALAPACK =   -L/home/export/online1/mdt00/shisuan/swict/2020_wuxi/2020_1008_mathlib/01_no_opt/mathlib_beijing      -lscalapack
#------------end no_opt-----------------------------------



#ARCHITECTURE = Generic
 CC           = mpicc
 CFLAGS       = -c -std=gnu99 -O3 -g -mieee -I /home/export/base/shisuan/swyjs/cq/swcl/install/lib/swcl/host/include
 SHM          = shm.o /usr/lib64/librt.a
 # OCL_LINKER_FLAGS = -lOpenCL 


#
# First example: ifort, mkl, serial mode
#
# The corresponding make target is "make serial" .
#
# FC = ifort 
# FFLAGS = -O3 -fp-model precise
# F90FLAGS = $(FFLAGS)
# for debugging F90FLAGS = $(FFLAGS) -debug all -check all
# ARCHITECTURE = Generic
# LAPACKBLAS = -L/opt/intel/composerxe-2011.4.191/compiler/lib/intel64 \
#              -I/opt/intel/composerxe-2011.4.191/mkl/include/intel64/lp64 \
#              -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
# USE_MPI = no
#
# End first example

# Second example: ifort, mkl, MPI build without ScaLapack.
#
# The corresponding make target is "make mpi" .
#
# FC = ifort 
# FFLAGS = -O3 -module $(MODDIR) -fp-model precise
# F90FLAGS = $(FFLAGS) -check bounds -check uninit -check pointers -traceback -g -fpe0
# ARCHITECTURE = Generic
# LAPACKBLAS = -mkl
# USE_MPI = yes
# MPIFC = mpif90
#
# End second example

# Third example: ifort, mkl, MPI build WITH ScaLapack.
#
# The corresponding make target is "make scalapack.mpi" .
#
# FC = ifort 
# FFLAGS = -O3 -ip -module $(MODDIR) -fp-model precise
# ### Alternative for ifort (slow! testing only) FFLAGS = -O3 -module $(MODDIR) -check bounds -check uninit -check pointers -traceback -g -fpe0
# F90FLAGS = $(FFLAGS)
# ARCHITECTURE = Generic
# USE_MPI = yes
# MPIFC = mpiifort
# SCALAPACK = -static-mpi -L/afs/@cell/common/soft/intel/ics/2011.0.013/12.1/mkl/lib/intel64/ -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm
#
# End third example

####################################################################################################
# SPECIAL PURPOSE (low level) variables
#
# In general, you will not need to bother with the following variables
# explicitly because sensible defaults should be chosen from your above
# settings.
#
# ARCH_MODULE: Use this file instead of arch_generic.f90 for compiler-specific
#              code.
#
# ELPA_KERNEL: If different from the generic elpa2_kernels.f90, a drop-in
#              replacement for it which should be highly optimized to the
#              host system (e.g. Blue Gene).
#
# LDFLAGS    : Flags passed (only) to the linker.
#
# CBFLAGS    : Used along with FFLAGS; probably meant to be set from the
#              command line  "make CBFLAGS='-check bounds'".  Depricated.
#
# MPIF       : Defaults to "-I./external" for USE_MPI=no and is used to
#              make our mpif.h stub available to the compiler.  Depending
#              on the compiler, "-I" might need to be replaced.
#
# BINDIR     : Directory for compiled and linked binaries [../bin]
#
# BUILDDIR   : Default directory for compiled object and module files [.]
#
# OBJDIR     : Directory for compiled object files [$(BUILDDIR)]
#
# MODDIR     : Directory for compiled module headers [$(BUILDDIR)]
#
# AUTODEPEND : If active, automatically check for dependencies among
#              the modules.  This allows for parallel compiling (make -j 4)
#              and saves from having to do "make clean" after every tiny
#              change. Default: "yes".
#
# AUTODEPEND_TYPE : There are two different dependency generators, one written
#                   in perl, the other in python.  The latter does some
#                   additional checks but is less portable. Either "perl" or
#                   "python".  Default: "perl".
#
# F90MINFLAGS  : keyword that specifies flags without optimization for faster compilation time
#                for files that are never a speed bottleneck at runtime. Example: read_control.f90 
#                with IBM compilers.
#
# USE_LIBXC : Set to "yes" to compiles aims with the version of libxc distributed with aims.  Requires
#             a C compiler be specified.
# LIBXC_DIR : /path/to/libxc so as to link in the external library of XC functionals (only recommended
#             for advanced developers!)
#
################################################################################################################
#  This is where the backend should be included for all generic settings
#
  include Makefile.backend
#
################################################################################################################
#
#  Modules must always be compiled before files which use them.  In the
#  default setting, there is a hard coded sequence in which the files are
#  compiled.  After each (nontrivial) change to a module file, the only save
#  way to build is to recompile after a "make clean".  After a "make depend",
#  make is aware of the module dependencies and it should be safe to just
#  run "make <target>" even after changes to module files.  This has the
#  additional advantage that one can use parallel builds via "make -j <np>".
#
#  But do not forget to rerun "make depend" or at least "make clean" after
#  adding new files or changing "use" clauses.  And do not forget to add new
#  files to Makefile.backend, either.  You can check for the need to "make
#  depend" on every make invocation by adding "AUTODEPEND = yes" to your
#  configuration.
#
################################################################################################################


#************************************************************************
#
# All platform-independent definitions, compilation etc follow here ...
#
#************************************************************************

# Make targets - potentially, you may need to adapt one of these targets for your own needs.

serial: $(BINDIR)/aims.$(AIMSVERSION).serial.x

mpi: $(BINDIR)/aims.$(AIMSVERSION).mpi.x

libaims: $(BINDIR)/libaims.a

scalapack.libaims: $(BINDIR)/libaims.scalapack.a

scalapack.mpi: $(BINDIR)/aims.$(AIMSVERSION).scalapack.mpi.x

multi.scalapack.mpi: $(BINDIR)/aims.$(AIMSVERSION).multi.scalapack.mpi.x

libaims.mpi: $(LIBDIR)/libaims.$(AIMSVERSION).mpi.so

libaims.scalapack.mpi: $(LIBDIR)/libaims.$(AIMSVERSION).scalapack.mpi.so

# The parser binary simply stops after reading in the input files.  This is
# useful for checking control.in and geometry.in for syntax errors.
parser: $(BINDIR)/aims.$(AIMSVERSION).parser.x

# Following are the actual compilation options for the above make targets.

$(BINDIR)/aims.$(AIMSVERSION).serial.x: $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) \
   $(ANACON) $(MPI_STUBS) $(SCALAPACK_STUBS) $(METIS_STUB) $(CUBA_STUBS) \
   $(C_FILES_STUBS) $(ASM_FILES) $(SHM_STUBS) $(PLUMED_STUB) $(PARSER_STUBS) \
   $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/aims.o \
              $(OBJDIR)/main.o $(SUBROUTINES) $(MODULES) $(EXTERNAL) $(ANACON) \
              $(LAPACKBLAS) $(MPI_STUBS) $(SCALAPACK_STUBS) $(METIS_STUB) \
              $(CUDALIBS) $(CUBA_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
              $(SHM_STUBS) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) \
              $(LIBXC_STUBS) $(HDF5_STUBS) $(SOC) $(IPC_FILES) $(PERF_STUBS) \
              $(LIBXC_LINKER_FLAGS)
	@echo "Built $@."

$(BINDIR)/aims.$(AIMSVERSION).mpi.x: $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) \
   $(ANACON) $(METIS_STUB) $(SCALAPACK_STUBS) $(CUBA_STUBS) $(SHM_STUBS) \
   $(C_FILES_STUBS) $(ASM_FILES) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) \
   $(LIBXC_STUBS) $(HDF5_STUBS) $(CUDA_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS)
	@mkdir -p $(dir $@)
	$(MPIFC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/aims.o \
              $(OBJDIR)/main.o $(SUBROUTINES) $(MODULES) $(EXTERNAL) $(ANACON) \
              $(METIS_STUB) $(CUDALIBS) $(SCALAPACK_STUBS) $(LAPACKBLAS) \
              $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) \
              $(ASM_FILES) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) \
              $(LIBXC_STUBS) $(HDF5_STUBS) $(SOC) $(IPC_FILES) $(PERF_STUBS) \
              $(LIBXC_LINKER_FLAGS)
	@echo "Built $@."

$(BINDIR)/aims.$(AIMSVERSION).scalapack.mpi.x: $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) \
   $(ANACON) $(METIS_STUB) $(CUBA_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
   $(SHM_STUBS) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) \
   $(HDF5_STUBS) $(SOC) $(CUDA_STUBS) $(IPC_FILES) $(PERF_STUBS) $(SLAVE_OBJS) /home/export/online1/mdt00/shisuan/swyjs/wzk/share/lib/sw_slave_dgemm.o
	@mkdir -p $(dir $@)
	$(MPIFC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/aims.o \
              $(OBJDIR)/main.o $(SUBROUTINES) $(MODULES) $(ANACON) $(EXTERNAL) \
              $(METIS_STUB) $(CUDALIBS) $(ELSI_LIB) $(SCALAPACK) $(LAPACKBLAS) \
              $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
              $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) \
              $(HDF5_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS) $(SLAVE_OBJS) -lm_slave \
              $(LIBXC_LINKER_FLAGS) $(OCL_LINKER_FLAGS) -lrt $(SWCL_LIB) /home/export/online1/mdt00/shisuan/swyjs/wzk/share/lib/sw_slave_dgemm.o
	@echo "Built $@."

$(BINDIR)/aims.$(AIMSVERSION).multi.scalapack.mpi.x: \
   $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/multiaims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) \
   $(ANACON) $(METIS_STUB) $(CUBA_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
   $(SHM_STUBS) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) \
   $(HDF5_STUBS) $(IPC_FILES) $(PERF_STUBS)
	@mkdir -p $(dir $@)
	$(MPIFC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/multiaims.o \
              $(OBJDIR)/main.o $(SUBROUTINES) $(MODULES) $(ANACON) $(EXTERNAL) \
              $(METIS_STUB) $(CUDALIBS) $(ELSI_LIB) $(SCALAPACK) $(LAPACKBLAS) \
              $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
              $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) \
              $(HDF5_STUBS) $(IPC_FILES) $(PERF_STUBS) $(LIBXC_LINKER_FLAGS)
	@echo "Built $@."


$(BINDIR)/aims.$(AIMSVERSION).parser.x: $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) \
   $(ANACON) $(MPI_STUBS) $(SCALAPACK_STUBS) $(METIS_STUB) $(CUBA_STUBS) \
   $(C_FILES_STUBS) $(ASM_FILES) $(SHM_STUBS) $(PLUMED_STUB) $(PARSER_OBJECTS) \
   $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) $(PERF_STUBS)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/aims.o \
              $(OBJDIR)/main.o $(SUBROUTINES) $(MODULES) $(EXTERNAL) $(ANACON) \
              $(LAPACKBLAS) $(MPI_STUBS) $(SCALAPACK_STUBS) $(METIS_STUB) \
              $(CUDALIBS) $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) \
              $(ASM_FILES) $(PLUMED_STUB) $(PARSER_OBJECTS) $(IPI_STUBS) \
              $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS) \
              $(LIBXC_LINKER_FLAGS)
	@echo "Built $@."

$(BINDIR)/libaims.a: $(OBJDIR)/write_version_stamp.o $(MODULES) \
      $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) $(ANACON) \
      $(METIS_STUB) $(SCALAPACK_STUBS) $(CUBA_STUBS) $(SHM_STUBS) \
      $(C_FILES_STUBS) $(ASM_FILES) $(PLUMED_STUB) $(PARSER_STUBS) \
      $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) $(PERF_STUBS)
	ar r $(BINDIR)/libaims.a \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/main.o $(SUBROUTINES) \
              $(MODULES) $(EXTERNAL) $(ANACON) $(METIS_STUB) $(CUDALIBS) \
              $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
              $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) \
              $(HDF5_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS)
	@echo "Built $@."

$(BINDIR)/libaims.scalapack.a: $(OBJDIR)/write_version_stamp.o $(MODULES) \
      $(OBJDIR)/aims.o $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) $(ANACON) \
      $(METIS_STUB) $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) $(ASM_FILES) \
      $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) \
      $(IPC_FILES) $(PERF_STUBS)
	ar r $(BINDIR)/libaims.a \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/main.o $(SUBROUTINES) \
              $(MODULES) $(EXTERNAL) $(ANACON) $(METIS_STUB) $(CUDALIBS) \
              $(ELSI_LIB) $(CUBA_STUBS) $(SHM_STUBS) $(C_FILES_STUBS) \
              $(ASM_FILES) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) \
              $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS)
	@echo "Built $@."

# To compile aims as a shared library the -fPIC F90FLAG is required.
$(LIBDIR)/libaims.$(AIMSVERSION).mpi.so: $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) $(ANACON) \
   $(METIS_STUB) $(SCALAPACK_STUBS) $(CUBA_STUBS) $(C_FILES_STUBS) \
   $(ASM_FILES) $(SHM_STUBS) $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) \
   $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) $(SOC) $(PERF_STUBS)
	@mkdir -p $(dir $@)
	$(MPIFC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) $(SHARED_OPT) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/main.o $(SUBROUTINES) \
              $(MODULES) $(ANACON) $(EXTERNAL) $(METIS_STUB) $(CUDALIBS) \
              $(LAPACKBLAS) $(SCALAPACK_STUBS) $(CUBA_STUBS) $(SHM_STUBS) \
              $(C_FILES_STUBS) $(ASM_FILES) $(PLUMED_STUB) $(PARSER_STUBS) \
              $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) $(IPC_FILES) \
              $(PERF_STUBS) $(LIBXC_LINKER_FLAGS)
	@echo "Built $@."

# To compile aims as a shared library the -fPIC F90FLAG is required.
$(LIBDIR)/libaims.$(AIMSVERSION).scalapack.mpi.so: $(OBJDIR)/write_version_stamp.o \
   $(MODULES) $(OBJDIR)/main.o $(SUBROUTINES) $(EXTERNAL) $(ANACON) \
   $(METIS_STUB) $(CUBA_STUBS) $(C_FILES_STUBS) $(ASM_FILES) $(SHM_STUBS) \
   $(PLUMED_STUB) $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) \
   $(IPC_FILES) $(SOC) $(PERF_STUBS)
	@mkdir -p $(dir $@)
	$(MPIFC) $(FFLAGS) $(CBFLAGS) $(LDFLAGS) $(SHARED_OPT) -o $@ \
              $(OBJDIR)/write_version_stamp.o $(OBJDIR)/main.o $(SUBROUTINES) \
              $(MODULES) $(ANACON) $(EXTERNAL) $(METIS_STUB) $(CUDALIBS) \
              $(ELSI_LIB) $(SCALAPACK) $(LAPACKBLAS) $(CUBA_STUBS) \
              $(SHM_STUBS) $(C_FILES_STUBS) $(ASM_FILES) $(PLUMED_STUB) \
              $(PARSER_STUBS) $(IPI_STUBS) $(LIBXC_STUBS) $(HDF5_STUBS) \
              $(IPC_FILES) $(PERF_STUBS) $(LIBXC_LINKER_FLAGS)
	@echo "Built $@."
