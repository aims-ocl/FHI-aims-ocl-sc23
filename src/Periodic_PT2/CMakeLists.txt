# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/mod_cpt2_blacs.f90
  )
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_cpt2_energy_kspace.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_cpt2_energy_kspace_blacs.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_cpt2_os_energy_kspace.f90
  ${CMAKE_CURRENT_LIST_DIR}/exchange_coulomb_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/exchange_ev_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/exchange_trico_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_v_multi_ovlp3fn_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/initialize_cpt2_para.f90
  ${CMAKE_CURRENT_LIST_DIR}/lvl_tricoeff_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/post_cpt2_blacs.f90
  ${CMAKE_CURRENT_LIST_DIR}/post_cpt2_lapack.f90
  ${CMAKE_CURRENT_LIST_DIR}/power_auxmat_lapack_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/power_auxmat_scalapack_complex_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/power_auxmat_scalapack_real_cpt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/pt2_calculation.f90
  ${CMAKE_CURRENT_LIST_DIR}/redistribute_cpt2_kspace.f90
  ${CMAKE_CURRENT_LIST_DIR}/restart_pt2.f90)
