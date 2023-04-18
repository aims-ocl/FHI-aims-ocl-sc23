# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/bse_2d.f90
  ${CMAKE_CURRENT_LIST_DIR}/bse.f90
  ${CMAKE_CURRENT_LIST_DIR}/bse_serial_wrapper.f90
  ${CMAKE_CURRENT_LIST_DIR}/bse_parallel_wrapper.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_solve_bse_mat.f90
  ${CMAKE_CURRENT_LIST_DIR}/coulomb_4ks.f90
  ${CMAKE_CURRENT_LIST_DIR}/qpe_calculation_bse.f90
  ${CMAKE_CURRENT_LIST_DIR}/screened_coulomb_4ks.f90
  ${CMAKE_CURRENT_LIST_DIR}/screened_coulomb_abf.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_coulomb_integrals_mo.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_coulomb_integrals_mo_2d.f90)
