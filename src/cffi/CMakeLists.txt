# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

enable_language(C)
find_package(PythonInterp)
target_include_directories(aims PRIVATE ${CMAKE_CURRENT_LIST_DIR})
add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/python_interface.c
  COMMAND
  ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cffi/build_python_interface.py
  DEPENDS cffi/build_python_interface.py cffi/python_interface.py
  COMMENT "Generating ${PROJECT_BINARY_DIR}/python_interface.c")
target_sources(aims3 PRIVATE
  aims_c_api.f90
  python_interface.f90
  ${PROJECT_BINARY_DIR}/python_interface.c)
