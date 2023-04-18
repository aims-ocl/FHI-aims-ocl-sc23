# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_dftpt2_correlation.f90
  ${CMAKE_CURRENT_LIST_DIR}/lrc_pt2.f90
  ${CMAKE_CURRENT_LIST_DIR}/prepare_lrc_corr_energy_calc.f90)
