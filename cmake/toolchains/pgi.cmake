set(TARGET_NAME aims.x CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

###############
# Basic Flags #
###############
set(CMAKE_Fortran_COMPILER $ENV{MPI_ROOT}/bin/mpifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-fast -O3 -Mvect=nosimd" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0" CACHE STRING "")
set(LIB_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "")
set(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(CMAKE_C_COMPILER pgcc CACHE STRING "")
set(CMAKE_C_FLAGS "-fast -O3" CACHE STRING "")
