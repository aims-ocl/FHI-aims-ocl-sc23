# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/dimensions_soc.f90
  ${CMAKE_CURRENT_LIST_DIR}/scalapack_soc.f90
  ${CMAKE_CURRENT_LIST_DIR}/soc_utilities.f90
  )
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/calculate_second_variational_soc.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_soc_hamiltonian.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_soc_shell.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_dipolematrix_soc.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_soc_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/soc_debug.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_soc_matrix.f90)
