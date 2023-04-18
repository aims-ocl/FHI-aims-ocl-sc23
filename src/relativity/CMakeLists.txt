# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/rel_x2c_mod.f90
  ${CMAKE_CURRENT_LIST_DIR}/rel_basis_index.f90
  ${CMAKE_CURRENT_LIST_DIR}/rel_diag_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/atom_x2c.f90
  ${CMAKE_CURRENT_LIST_DIR}/scalar2spinor.f90
  ${CMAKE_CURRENT_LIST_DIR}/sigma_dot_p_int.f90
  )
target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/rel_integrations.f90
  ${CMAKE_CURRENT_LIST_DIR}/rel_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/rel_scf_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/rel_testfunctions.f90
  ${CMAKE_CURRENT_LIST_DIR}/x2c_functions.f90
  ${CMAKE_CURRENT_LIST_DIR}/x2c.f90)
