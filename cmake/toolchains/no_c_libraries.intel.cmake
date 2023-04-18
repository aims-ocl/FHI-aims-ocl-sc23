set(TARGET_NAME aims.x CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

###############
# Basic Flags #
###############
set(ARCHITECTURE Generic CACHE STRING "")
set(CMAKE_Fortran_COMPILER $ENV{MPI_ROOT}/bin/mpifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0 -fp-model precise" CACHE STRING "")
set(LIB_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "")
set(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(USE_C_FILES OFF CACHE BOOL "")
set(USE_CXX_FILES OFF CACHE BOOL "")
