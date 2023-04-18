# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

# Copied from stackoverflow.com/a/19578320
string(ASCII 27 Esc)
set(ColorReset  "${Esc}[m")
set(ColorBold   "${Esc}[1m")
set(Red         "${Esc}[31m")
set(Green       "${Esc}[32m")
set(Yellow      "${Esc}[33m")
set(Blue        "${Esc}[34m")
set(Magenta     "${Esc}[35m")
set(Cyan        "${Esc}[36m")
set(White       "${Esc}[37m")
set(BoldRed     "${Esc}[1;31m")
set(BoldGreen   "${Esc}[1;32m")
set(BoldYellow  "${Esc}[1;33m")
set(BoldBlue    "${Esc}[1;34m")
set(BoldMagenta "${Esc}[1;35m")
set(BoldCyan    "${Esc}[1;36m")
set(BoldWhite   "${Esc}[1;37m")

# Determine current flags
macro(get_all_flags language flags)
  set(${flags} ${CMAKE_${language}_FLAGS})
  if (CMAKE_BUILD_TYPE)
    string(TOUPPER ${CMAKE_BUILD_TYPE} buildtype)
    set(${flags} "${${flags}} ${CMAKE_${language}_FLAGS_${buildtype}}")
    string(STRIP ${${flags}} ${flags})
  endif()
endmacro()

# Convert a string into a list. If the variable was already given as
# list, do nothing.
macro(convert_to_list var)
  list(LENGTH ${var} _length)
  if (_length EQUAL 1)
    string(REPLACE " " ";" ${var} ${${var}})
  endif()
endmacro()

# Go through directories listed in LIBRARY_PATHS and turn the entries
# of LIBRARIES into actual targets.
function(generate_library_targets _PATHS _LIBRARIES)
  foreach(LIB ${${_LIBRARIES}})
    # If the target already exists, skip it. This can occur when the
    # same library is set multiple times in LIBRARIES due to circular
    # dependencies.
    if (NOT TARGET ${LIB})
      find_library(LIB_FULLPATH ${LIB} HINTS ${${_PATHS}})
      if (LIB_FULLPATH)
        message(STATUS "Found ${LIB_FULLPATH}")
        add_library(${LIB} UNKNOWN IMPORTED)
        set_target_properties(${LIB} PROPERTIES
          IMPORTED_LOCATION ${LIB_FULLPATH})
        unset(LIB_FULLPATH CACHE)
      else()
        message(FATAL_ERROR "${Magenta}Could not find ${LIB}${ColorReset}")
      endif()
    endif()
  endforeach()
endfunction()

# Use this function for older CMake versions that do not support
# list(FILTER ...).
function(list_filter list regex)
  set(_tmp_list "")
  foreach(_src ${${list}})
    string(REGEX MATCH "${regex}" match ${_src})
    if (match)
      list(APPEND _tmp_list ${_src})
    endif()
  endforeach()
  set(${list} ${_tmp_list} PARENT_SCOPE)
endfunction()

# In the very specific situation of linking against ELSI, where some
# ELSI components have calls to linear algebra libraries with circular
# dependencies and the libraries are linked statically, CMake is
# unable to produce the correct order and multiplicity of the
# libraries on the linking line. This arises specifically with the MKL
# libraries, which need to be repeated at least a few times in order
# to remove all circular dependencies. This function presents a
# workaround where, if linking statically, the main target explicitly
# links against all targets that are propagated with the dependant
# target (e.g., elsi::elsi). See
# https://gitlab.kitware.com/cmake/cmake/issues/18415 for a discussion
# on this topic.
function(static_linking_fix main_target dep_target libs)
  foreach (_lib ${${libs}})
    get_target_property(LIB_FULLPATH ${_lib} IMPORTED_LOCATION)
    string(REGEX MATCH "\.a$" _result ${LIB_FULLPATH})
    if (_result)
      set(_has_static_libs TRUE)
    endif()
  endforeach()
  if (_has_static_libs)
    get_target_property(_libs ${dep_target} INTERFACE_LINK_LIBRARIES)
    target_link_libraries(${main_target} PRIVATE ${_libs})
  endif()
endfunction()

# Create library with normal flags
macro(add_library_normal target_name)
  add_library0(${target_name} ${CMAKE_Fortran_COMPILER} ${ARGN})
endmacro()

# Create library with minimal flags
macro(add_library_minimal target_name)
  set(_tmp_flags ${Fortran_MIN_FLAGS})
  add_library0(${target_name} ${CMAKE_Fortran_COMPILER} ${ARGN})
  unset(_tmp_flags)
endmacro()

# Create library with reduced flags
macro(add_library_reduced target_name)
  if (KNOWN_COMPILER_BUGS STREQUAL "ifort-14.0.1-O3")
    set(_tmp_flags -O1)
  elseif (KNOWN_COMPILER_BUGS STREQUAL "ifort-17.0.0-O3")
    set(_tmp_flags -O2)
  elseif (REDUCED_COMPILER STREQUAL "quick_develop")
    set(_tmp_flags -O0)
  endif()
  add_library0(${target_name} ${CMAKE_Fortran_COMPILER} ${ARGN})
  unset(_tmp_flags)
endmacro()

# Use FFLAGS for .f files
macro(add_library_f77 target_name)
  set(_tmp_flags ${FFLAGS})
  add_library0(${target_name} ${CMAKE_Fortran_COMPILER} ${ARGN})
  unset(_tmp_flags)
endmacro()

macro(add_library0 target_name compiler)
  if (NOT project_root_dir)
    set(project_root_dir ${PROJECT_SOURCE_DIR})
  endif()
  if (NOT EXTERNAL_ELSI_PATH AND ${target_name} STREQUAL aims1)
    add_subdirectory(external/elsi)
  endif()
  file(WRITE ${PROJECT_BINARY_DIR}/subdirectory_${target_name}/CMakeLists.txt
    "
project(project_${target_name} LANGUAGES Fortran C)
set(CMAKE_Fortran_COMPILER ${compiler})
if (_tmp_flags)
  if (CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} buildtype)
  set(CMAKE_Fortran_FLAGS_${buildtype} \"\")
  endif()
  set(CMAKE_Fortran_FLAGS \"${_tmp_flags}\")
endif()
add_library(${target_name} \"\")
if (${ARGC} EQUAL 3)
target_link_libraries(${target_name} PUBLIC ${ARGV2})
endif()
set_target_properties(${target_name} PROPERTIES
Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
target_link_libraries(${target_name} PRIVATE ${LIBS})
if (EXTERNAL_ELSI_PATH)
  find_package(elsi 2.1 REQUIRED NO_DEFAULT_PATH PATHS ${EXTERNAL_ELSI_PATH})
  target_link_libraries(${target_name} PRIVATE elsi::elsi)
  static_linking_fix(${target_name} elsi::elsi LIBS)
else()
if (${target_name} STREQUAL aims1)
  target_link_libraries(${target_name} PRIVATE elsi)
  static_linking_fix(${target_name} elsi LIBS)
endif()
endif()
if (USE_LIBXC AND LIBXC_VERSION AND ${target_name} STREQUAL aims1)
  target_link_libraries(${target_name} PUBLIC libxc-custom)
endif()
if (USE_CUDA AND NOT ${CMAKE_VERSION} VERSION_LESS 3.8)
  enable_language(CUDA)
endif()
target_include_directories(${target_name} PRIVATE ${INC_PATHS})
target_include_directories(${target_name} PRIVATE ${PROJECT_BINARY_DIR}/modules)
"
    )
  add_subdirectory(${PROJECT_BINARY_DIR}/subdirectory_${target_name}
    ${PROJECT_BINARY_DIR}/subdirectory_${target_name}/build)
endmacro()

macro(target_sources target_name mode)
  foreach(_src ${ARGN})
    if (EXISTS ${_src})
      _target_sources(${target_name} ${mode} ${_src})
    else()
      _target_sources(${target_name} ${mode} ${CMAKE_CURRENT_SOURCE_DIR}/${_src})
    endif()
  endforeach()
endmacro()

macro(pgi_libs_full_fix)
  string(REGEX REPLACE "([a-zA-Z0-9_\.]+/)" "" LIBS_FULL "${LIBS_FULL}")
  string(REGEX REPLACE "(/|.so|lib)" "" LIBS_FULL "${LIBS_FULL}")
  if (${CMAKE_VERSION} VERSION_LESS 3.2)
    string(LENGTH "${LIBS_FULL}" _length)
    if (${_length} LESS 264)
      string(SUBSTRING "${LIBS_FULL}" 0 -1 LIBS_FULL)
    else()
      string(SUBSTRING "${LIBS_FULL}" 0 264 LIBS_FULL)
    endif()
  else()
    string(SUBSTRING "${LIBS_FULL}" 0 264 LIBS_FULL)
  endif()
endmacro()
