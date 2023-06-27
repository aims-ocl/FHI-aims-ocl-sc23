# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_gradient_rho_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_polar_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/prune_density_matrix_sparse_polar_reduce_memory.f90)
