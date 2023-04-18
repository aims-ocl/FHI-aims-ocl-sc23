# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_right.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_right_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_right_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_right.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_right_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_right_reduce_memory.f90)
