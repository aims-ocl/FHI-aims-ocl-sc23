# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_interface.f90
  ${CMAKE_CURRENT_LIST_DIR}/test_symmetry_of_grid_and_things_on_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_index_map.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_custom_batches.f90
  ${CMAKE_CURRENT_LIST_DIR}/batch_rlsy_interface.f90
  )
