set(TARGET_NAME aims.x CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

###############
# Basic Flags #
###############
set(CMAKE_Fortran_COMPILER mpiifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0 -fp-model precise" CACHE STRING "")
set(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 $ENV{ELSI_HOME}/lib" CACHE STRING "")
set(LIBS "OMM MatrixSwitch elpa NTPoly fortjson pexsi superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotcherr scotch mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(CMAKE_C_COMPILER icc CACHE STRING "")
set(CMAKE_C_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")

##########################
# External Library Flags #
##########################
set(EXTERNAL_ELSI_PATH "$ENV{ELSI_HOME}" CACHE STRING "")
