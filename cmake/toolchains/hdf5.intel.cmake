set(TARGET_NAME aims.x CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

###############
# Basic Flags #
###############
set(CMAKE_Fortran_COMPILER mpiifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0 -fp-model precise" CACHE STRING "")
set(INC_PATHS "$ENV{HDF5_HOME}/include" CACHE STRING "")
set(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 $ENV{HDF5_HOME}/lib" CACHE STRING "")
set(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core hdf5_fortran hdf5hl_fortran hdf5" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(USE_CXX_FILES ON CACHE BOOL "")
set(CMAKE_C_COMPILER icc CACHE STRING "")
set(CMAKE_C_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")

##########################
# Optional Library Flags #
##########################
set(USE_HDF5 ON CACHE BOOL "")
