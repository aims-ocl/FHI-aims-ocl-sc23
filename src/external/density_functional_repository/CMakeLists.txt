# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims6 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_lyp.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_p86.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_pbe.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_pw91.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_pw92.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_pz81.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_vwn5.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_c_vwn5rpa.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_x_b3.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_x_b88.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_b3lyp.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_b97_1.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_b97_2.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_b97.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_edf1.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_ft97.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_hcth120.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_hcth147.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_hcth407.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_hcth.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_pbe.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_xc_pw91.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_x_lda.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_x_pbe.f
  ${CMAKE_CURRENT_LIST_DIR}/f77_x_pw91.f)
