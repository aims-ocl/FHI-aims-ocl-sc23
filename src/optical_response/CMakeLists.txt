# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

add_subdirectory(libxc_stub)

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/libxc_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/optical_response.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_fxc_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/optics_module.f90)
