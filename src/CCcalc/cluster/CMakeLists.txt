# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/cc_calc_cl.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_aux.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_calc_w.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_DIIS.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_distribution.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_initial.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_intl.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_jacob.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_mem.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_mod.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_mpi.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_cl_PT.f90)
