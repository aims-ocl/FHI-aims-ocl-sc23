# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/cc_calc_aux.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_calc_E.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_calculation.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_calc_w_vect.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_initial.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_mod.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_RLE.f90
  ${CMAKE_CURRENT_LIST_DIR}/cc_util.f90)
