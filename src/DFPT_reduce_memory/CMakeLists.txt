# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_EDM_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_gradient_rho_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_U_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_zero_order_DM_EDM_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_hellman_feynman_hessian_fixed_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_pulay_hessian_fixed_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_DM_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_EDM_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_hessian_reduce_memory.f90)
