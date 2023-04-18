# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/ci_acceleration.f90
  ${CMAKE_CURRENT_LIST_DIR}/ci_acc_initial.f90
  ${CMAKE_CURRENT_LIST_DIR}/ci_calculation.f90
  ${CMAKE_CURRENT_LIST_DIR}/ci_calculation_slab.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_c_vector.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_E_ci.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_w_vector_1a1b.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_w_vector_1a.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_w_vector_1b.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_w_vector_1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_w_vector_2a.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_w_vector_2b.f90
  ${CMAKE_CURRENT_LIST_DIR}/fciqmc_module.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_4_index_integrals_slab.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_array.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_CI_amplitudes.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_FCIDUMP.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_matrix_2D.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_w_vect.f90)
