# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_aux.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_calc_w.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_diis.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_distribution.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_initial.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_intl.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_mem.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_mod.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_3d_mpi.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_calc_3d.f90)
