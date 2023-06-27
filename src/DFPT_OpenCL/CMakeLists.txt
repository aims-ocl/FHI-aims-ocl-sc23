# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

message(STATUS "USE_OPENCL = ${USE_OPENCL}")

if (USE_OPENCL)
  add_library_normal(aims-ocl)

  if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
  endif()

  set(OPENCL_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/sum_up.cl
    ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho.cl)

  foreach(OPENCL_SOURCE ${OPENCL_SOURCES})
    message(STATUS "${OPENCL_SOURCE}")

    file(READ ${OPENCL_SOURCE} OPENCL_CONTENTS)

    string(REPLACE "\"" "\\\"" OPENCL_CONTENTS "${OPENCL_CONTENTS}")
    string(REPLACE "\n" "\\n\"\n\"" OPENCL_CONTENTS "${OPENCL_CONTENTS}")
    set(OPENCL_CONTENTS "\"${OPENCL_CONTENTS}\\n\"")

    get_filename_component(OPENCL_FILENAME ${OPENCL_SOURCE} NAME_WE)
    configure_file("${CMAKE_CURRENT_LIST_DIR}/opencl_source.h.in" "${CMAKE_CURRENT_LIST_DIR}/opencl_source_${OPENCL_FILENAME}.h" @ONLY)

    # target_sources(aims-ocl PRIVATE "${CMAKE_CURRENT_LIST_DIR}/opencl_source_${OPENCL_FILENAME}.cpp")

  endforeach()

  configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake_help.h.in ${CMAKE_CURRENT_LIST_DIR}/cmake_help.h)

  target_sources(aims-ocl PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho.c
    ${CMAKE_CURRENT_LIST_DIR}/opencl_util_c.c
    ${CMAKE_CURRENT_LIST_DIR}/opencl_util_mpi.c
    ${CMAKE_CURRENT_LIST_DIR}/opencl_util_mpi_m.c
    ${CMAKE_CURRENT_LIST_DIR}/save_load_var.c
    ${CMAKE_CURRENT_LIST_DIR}/sum_up_whole_potential_shanghui.c)

  # target_link_libraries(aims-ocl PRIVATE fmt::fmt-header-only)
  
  find_package(OpenCL REQUIRED)
  message(STATUS "OpenCL found: ${OpenCL_FOUND}")
  message(STATUS "OpenCL version: ${OpenCL_VERSION_STRING}")
  message(STATUS "OpenCL includes: ${OpenCL_INCLUDE_DIRS}")
  message(STATUS "OpenCL libraries: ${OpenCL_LIBRARIES}")
  # INCLUDE_DIRECTORIES(${OpenCL_INCLUDE_DIRS})
  # LINK_DIRECTORIES(${OpenCL_LIBRARY})
  # -DOpenCL_FOUND=True -DOpenCL_LIBRARY= -DOpenCL_INCLUDE_DIR=
  
  target_include_directories(aims-ocl PRIVATE ${OpenCL_INCLUDE_DIRS})
  target_link_libraries(aims-ocl PRIVATE ${OpenCL_LIBRARIES}/libOpenCL.so)

  # find_package(multiarray REQUIRED)
  # target_link_libraries(aims-ocl PRIVATE multiarray::multiarray)
  # target_compile_features(aims-ocl PRIVATE cxx_std_17)

  target_compile_options(aims-ocl PRIVATE -Wall -Wno-unused-variable)
  target_include_directories(aims-ocl PRIVATE ${CMAKE_CURRENT_LIST_DIR})
  target_link_libraries(aims PRIVATE aims-ocl)
else()
  # target_sources(aims5 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/gpu_stub.f90)
endif()

target_sources(aims5 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/opencl_util.f90)