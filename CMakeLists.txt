# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

# CMAKE
cmake_minimum_required(VERSION 3.1)

# PROJECT
file(READ src/version_stamp.txt _version_stamp)
string(REGEX MATCH "AIMSVERSION *= *[0-9]+" _version ${_version_stamp})
string(REGEX REPLACE "AIMSVERSION *= *" "" _version ${_version})
# All languages are enabled later.
project(FHI-aims VERSION ${_version} LANGUAGES NONE)
# Disallow in-source builds
if (PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "
In-source builds are not allowed. Remove \
${PROJECT_BINARY_DIR}/CMakeCache.txt and run cmake again from \
a build directory. Example: mkdir build && cd build && cmake ..")
endif()
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
# Tools.cmake contains some helper functions
include(Tools)
# Stop if the Fortran compiler is not among the initial cache
# variables. Ignore if run with DEV.
if (NOT CMAKE_Fortran_COMPILER AND NOT DEV)
  message(FATAL_ERROR
    "\n !!! Error. Could not find the Fortran compiler specification (the CMAKE_Fortran_COMPILER variable). !!! \n\n !!! A Fortran compiler must be specified explicitly before proceeding. !!! \n\n !!! Also, check that you used 'cmake -C ...', not 'cmake -c ...' - this could be another reason why this error message triggers. !!!")
endif()
enable_language(Fortran)
set(MINIMUM_IFORT_VERSION 12.0.4)
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND NOT
    CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER ${MINIMUM_IFORT_VERSION} AND
    NOT CMAKE_VERSION VERSION_LESS 3.4)
  message(FATAL_ERROR "Intel compiler versions ${MINIMUM_IFORT_VERSION} and \
lower not supported. You are running version ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

# OPTIONS
#
# Do not forget to update documentation when making changes to this
# section!
#
# Base name of the compiled target.
set(TARGET_NAME "" CACHE STRING "Base name of the target")
# Whether to build multiaims or just aims
option(MULTIAIMS "Whether to build multiaims or regular aims" OFF)
# Whether to build a library or an executable
option(BUILD_LIBRARY "Build library or an executable" OFF)
# If BUILD_LIBRARY is ON, whether to build shared/static libraries
option(BUILD_SHARED_LIBS "Whether to build shared or static libraries" OFF)
# Whether to use MPI
if (DEFINED USE_MPI)
  # Auxiliary variable to detect whether MPI was initially set by the
  # user. This is used later in the FIRST RUN section.
  set(USE_MPI_DEFINED TRUE)
endif()
option(USE_MPI "Whether to compile with MPI support" OFF)
# Whether to use Scalapack
if (DEFINED USE_SCALAPACK)
  # Same logic as with USE_MPI_DEFINED above
  set(USE_SCALAPACK_DEFINED TRUE)
endif()
option(USE_SCALAPACK "Whether to compile with Scalapack support" OFF)

if (NOT DEFINED USE_RLSY)
    set(USE_RLSY ON CACHE BOOL "Whether to compile with RLSY support" FORCE)
endif()

# Architecture
set(ARCHITECTURE "" CACHE STRING "Choose AMD64_AVX/AMD64_SSE/PGI/BlueGene_Q")
if (ARCHITECTURE STREQUAL "AMD64_SSE")
  # Now that CMakeDetermineASMCompiler.cmake is read from the updated
  # CMAKE_MODULE_PATH, we can enable the assembly language. See
  # gitlab.kitware.com/cmake/cmake/merge_requests/2418 for more info.
  enable_language(ASM)
endif()
# C files
option(USE_C_FILES "Whether to compile with C source files" ON)
if (USE_C_FILES)
  enable_language(C)
endif()
# C++ files
option(USE_CXX_FILES "Whether to compile with C++ source files" OFF)
if (USE_CXX_FILES OR ENABLE_PEXSI)
  enable_language(CXX)
endif()
# LibXC
option(USE_LIBXC "Whether to compile with LibXC support" ${USE_C_FILES})
set(LIBXC_VERSION ""
  CACHE STRING "LibXC version to be downloaded and compiled into aims")
# Python 2/3 interface to FHI-aims
option(USE_CFFI "Whether to provide an interface to Python 2 or 3" OFF)
# Spglib
option(USE_SPGLIB "Whether to compile with Spglib support" ${USE_C_FILES})
# CUDA
option(USE_CUDA "Whether to enable GPU acceleration" OFF)
# ELPA_MT
option(ELPA_MT "Whether to use MPI/OpenMP version of ELPA" OFF)
# IPC
option(USE_IPC "Whether to support inter-process communication " OFF)
# i-Pi
option(USE_iPI "Whether support the i-Pi Python wrapper" ${USE_C_FILES})
# HDF5
option(USE_HDF5 "Whether to compile with HDF5 support" OFF)
# External ELSI
set(EXTERNAL_ELSI_PATH "" CACHE STRING "External ELSI directory")
set(USE_GPU_ELPA "OFF" CACHE BOOL "Whether to default to GPU-accelerated ELPA")
# Linker flags
set(CMAKE_EXE_LINKER_FLAGS "" CACHE STRING
  "Flags used by the linker when creating a target")
# Flags for .f files
set(FFLAGS "" CACHE STRING "Flags used for .f files")
# List of directories containing any external libraries
set(LIB_PATHS "" CACHE STRING
  "Any directories containing libraries to be linked against")
# List of libraries to link against
set(LIBS "" CACHE STRING
  "Any libraries to be linked against")
# List of directories with header files
set(INC_PATHS "" CACHE STRING
  "Any directories containing necessary header files")
string(REPLACE "-l" "" LIBS "${LIBS}")
# Split lists by spaces
foreach(_list LIBS LIB_PATHS INC_PATHS)
  convert_to_list(${_list})
endforeach()
# Link against the libraries
generate_library_targets(LIB_PATHS LIBS)

# FIRST RUN
# Determine the defaults for some variables here.
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  # Automatically detect USE_MPI
  if (NOT USE_MPI_DEFINED)
    try_compile(MPI_OK ${PROJECT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/cmake/tests/MPI.f90)
    if (MPI_OK)
      set(USE_MPI ON CACHE BOOL "Whether to compile with MPI support" FORCE)
    endif()
  endif()
  # Automatically detect USE_SCALAPACK
  if (NOT USE_SCALAPACK_DEFINED)
    try_compile(SCALAPACK_OK ${PROJECT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/cmake/tests/Scalapack.f90
      CMAKE_FLAGS "-DLINK_DIRECTORIES=${LIB_PATHS}"
      LINK_LIBRARIES ${LIBS})
    if (SCALAPACK_OK)
      set(USE_SCALAPACK ON CACHE BOOL
        "Whether to compile with Scalapack support" FORCE)
    endif()
  endif()
  # mpif.h for ELSI
  if(NOT USE_MPI)
    set(INC_PATHS "${INC_PATHS} ${PROJECT_SOURCE_DIR}/src/external" CACHE STRING
      "List of include directories" FORCE)
  endif()
  # Legacy
  if (LIBRARY_PATHS)
    set(LIB_PATHS ${LIBRARY_PATHS} CACHE STRING "" FORCE)
    message(WARNING "${Magenta}Use LIB_PATHS instead of LIBRARY_PATHS from now on (see documentation)${ColorReset}")
  endif()
  if (LIBRARIES)
    set(LIBS ${LIBRARIES} CACHE STRING "" FORCE)
    message(WARNING "${Magenta}Use LIBS instead of LIBRARIES from now on (see documentation)${ColorReset}")
    foreach(_list LIBS LIB_PATHS INC_PATHS)
      convert_to_list(${_list})
    endforeach()
    generate_library_targets(LIB_PATHS LIBS)
  endif()
  if (INCLUDE_PATHS)
    set(INC_PATHS ${INCLUDE_PATHS} CACHE STRING "" FORCE)
    message(WARNING "${Magenta}Use INC_PATHS instead of INCLUDE_PATHS from now on (see documentation)${ColorReset}")
  endif()
endif()

# Unless TARGET_NAME exists, set the base name of the target based on
# the current configuration.
if (NOT TARGET_NAME)
  if (USE_MPI)
    set(TARGET_NAME mpi)
  else()
    set(TARGET_NAME serial)
  endif()
  if (USE_SCALAPACK)
    set(TARGET_NAME scalapack.${TARGET_NAME})
  endif()
  if (USE_HDF5)
    set(TARGET_NAME hdf5.${TARGET_NAME})
  endif()
  if (MULTIAIMS)
    set(TARGET_NAME multi.${TARGET_NAME})
  endif()
  set(TARGET_NAME ${PROJECT_VERSION}.${TARGET_NAME})
  set(TARGET_NAME aims.${TARGET_NAME})
  if (NOT BUILD_LIBRARY)
    set(TARGET_NAME ${TARGET_NAME}.x)
  endif()
endif()

# OTHER
# We need this for including external libraries
include(ExternalProject)
# Collect all flags corresponding to the enabled languages.
get_all_flags(Fortran ALL_Fortran_FLAGS)
get_all_flags(C ALL_C_FLAGS)
get_all_flags(CXX ALL_CXX_FLAGS)
get_all_flags(CUDA ALL_CUDA_FLAGS)
# Unhide some variables
mark_as_advanced(CLEAR
  CMAKE_ASM_COMPILER
  CMAKE_Fortran_COMPILER
  CMAKE_Fortran_FLAGS
  CMAKE_C_FLAGS
  CMAKE_C_COMPILER
  CMAKE_CXX_FLAGS
  CMAKE_CXX_COMPILER)
if (NOT Fortran_MIN_FLAGS)
  set(Fortran_MIN_FLAGS ${ALL_Fortran_FLAGS} CACHE STRING "" FORCE)
endif()
# Expand the provided libraries to their full paths (for output
# purposes).
foreach(_LIB ${LIBS})
  get_target_property(_LIB_FULL ${_LIB} LOCATION)
  list(APPEND LIBS_FULL ${_LIB_FULL})
endforeach()
# Get git commit number
include(GetGitRevisionDescription)
get_git_head_revision(GIT_SHA1)
# Generate the cmake_info subroutine
if (CMAKE_Fortran_COMPILER_ID)
  string(REGEX MATCH ".*PGI.*" _is_pgi ${CMAKE_Fortran_COMPILER_ID})
endif()
if (_is_pgi)
  pgi_libs_full_fix()
endif()

# SOURCE
add_subdirectory(src)
