# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

if (USE_CUDA)
  target_sources(aims5 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/gpuError.cu
    ${CMAKE_CURRENT_LIST_DIR}/gpuInterface.cu
    ${CMAKE_CURRENT_LIST_DIR}/gpuIntegrateHam.cu
    ${CMAKE_CURRENT_LIST_DIR}/gpuDensity.cu
    ${CMAKE_CURRENT_LIST_DIR}/gpuForces.cu
    ${CMAKE_CURRENT_LIST_DIR}/gpuMR.cu
    ${CMAKE_CURRENT_LIST_DIR}/gpuMR.f90)
  target_include_directories(aims5 PRIVATE ${CMAKE_CURRENT_LIST_DIR})

  # CMake has proper CUDA support starting version 3.8
  if (${CMAKE_VERSION} VERSION_LESS 3.8)
    enable_language(CXX)
    find_package(CUDA REQUIRED)
    mark_as_advanced(CUDA_HOST_COMPILER CUDA_SDK_ROOT_DIR CUDA_TOOLKIT_ROOT_DIR
      CUDA_USE_STATIC_CUDA_RUNTIME)
    if (NOT CMAKE_CUDA_FLAGS)
      set(CMAKE_CUDA_FLAGS "" CACHE STRING "Flags for the CUDA compiler")
    endif()
  else()
    enable_language(CUDA)
    # See gitlab.kitware.com/cmake/cmake/issues/17797
    if (NOT CMAKE_GENERATOR STREQUAL Ninja AND
        ${CMAKE_VERSION} VERSION_LESS 3.11)
      set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_INCLUDES 0)
      set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_LIBRARIES 0)
      set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_OBJECTS 0)
    endif()
  endif()
  mark_as_advanced(CLEAR CMAKE_CUDA_COMPILER CMAKE_CUDA_FLAGS)
  get_all_flags(CUDA ALL_CUDA_FLAGS)
  # Check presence of cuBLAS
  if (NOT ALL_CUDA_FLAGS MATCHES "-lcublas")
    message(WARNING "${Magenta}-lcublas missing from CUDA flags${ColorReset}")
  endif()
  # If nvcc wants to link against an old version of libgfortran, then
  # this can be prevented by overriding the following variable.
  if (CUDA_LINK_DIRS)
    convert_to_list(CUDA_LINK_DIRS)
    set(CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES ${CUDA_LINK_DIRS})
  endif()
  # If using CMake version < 3.8, use the following hack.
  if (${CMAKE_VERSION} VERSION_LESS 3.8)
    if (NOT CMAKE_CUDA_COMPILER)
      message(FATAL_ERROR
        "CMAKE_CUDA_COMPILER must be set for GPU support (only for CMake <3.8)")
    endif()
    # Extract .cu files from all sources
    get_target_property(_gpu_sources aims5 SOURCES)
    if (${CMAKE_VERSION} VERSION_LESS 3.6)
      list_filter(_gpu_sources "\\.cu$")
    else()
      list(FILTER _gpu_sources INCLUDE REGEX "\\.cu$")
    endif()
    # Turn CMAKE_CUDA_FLAGS into a list
    string(REPLACE " " ";" _flags "${CMAKE_CUDA_FLAGS}")
    # Manually specify rules for compiling the .cu files
    foreach (_src ${_gpu_sources})
      get_filename_component(_obj ${_src} NAME)
      set(_obj ${PROJECT_BINARY_DIR}/${_obj}.o)
      execute_process(COMMAND
        ${CMAKE_CUDA_COMPILER} ${_flags} -c ${_src} -o ${_obj})
      _target_sources(aims5 PRIVATE ${_obj})
    endforeach()
  endif()
else()
  target_sources(aims5 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gpu_stub.f90)
endif()
