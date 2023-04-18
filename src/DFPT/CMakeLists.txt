# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_C.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_pulay.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_moving_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_U.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_zero_order_DM_EDM.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_hellman_feynman_hessian.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_hellman_feynman_hessian_fixed.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_second_order_H_pulay.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_second_order_H_pulay_fixed.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_second_order_S.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_second_order_S_fixed.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_xc_shanghui.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_hessian.f90
  ${CMAKE_CURRENT_LIST_DIR}/sum_up_whole_potential_p2_shanghui.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_hartree_potential_p2_shanghui.f90)
