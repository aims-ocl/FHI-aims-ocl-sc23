# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims8 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/anderson.f
  ${CMAKE_CURRENT_LIST_DIR}/atoaux.f
  ${CMAKE_CURRENT_LIST_DIR}/bpgcxc.f
  ${CMAKE_CURRENT_LIST_DIR}/cepvwn.f
  ${CMAKE_CURRENT_LIST_DIR}/corlyp.f
  ${CMAKE_CURRENT_LIST_DIR}/dcpv.f
  ${CMAKE_CURRENT_LIST_DIR}/dftseq.f
  ${CMAKE_CURRENT_LIST_DIR}/fxc_mgga_pk.f
  ${CMAKE_CURRENT_LIST_DIR}/gga91_sr.f
  ${CMAKE_CURRENT_LIST_DIR}/ggarad.f
  ${CMAKE_CURRENT_LIST_DIR}/ldaxc.f
  ${CMAKE_CURRENT_LIST_DIR}/moment.f
  ${CMAKE_CURRENT_LIST_DIR}/pbe.f
  ${CMAKE_CURRENT_LIST_DIR}/pberev.f
  ${CMAKE_CURRENT_LIST_DIR}/sratom_n.f
  ${CMAKE_CURRENT_LIST_DIR}/vestat.f
  ${CMAKE_CURRENT_LIST_DIR}/vexcor.f
  ${CMAKE_CURRENT_LIST_DIR}/vexcos.f)
