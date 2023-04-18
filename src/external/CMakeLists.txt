# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

if (NOT EXTERNAL_ELSI_PATH)
  add_dependencies(aims1 elsi)
  target_include_directories(aims1 PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/elsi/include)
endif()

add_subdirectory(libxc-4.0.2)
add_subdirectory(i-pi)
add_subdirectory(density_functional_repository)
add_subdirectory(spglib)
add_subdirectory(libmbd)

target_sources(aims8 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/boys.f
  ${CMAKE_CURRENT_LIST_DIR}/dffti.f
  ${CMAKE_CURRENT_LIST_DIR}/dfftb.f
  ${CMAKE_CURRENT_LIST_DIR}/dfftf.f
  ${CMAKE_CURRENT_LIST_DIR}/drc3jm.f
  ${CMAKE_CURRENT_LIST_DIR}/cdgamma.f
  ${CMAKE_CURRENT_LIST_DIR}/cor_vwn_rpa_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/cor_pbeh0_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/pbe_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/HSEFx_partial.f
  ${CMAKE_CURRENT_LIST_DIR}/ei.f
  ${CMAKE_CURRENT_LIST_DIR}/pbeint_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/am05.f
  ${CMAKE_CURRENT_LIST_DIR}/am05_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/b88_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/cepvwn-gauss-like.f
  ${CMAKE_CURRENT_LIST_DIR}/cor_vwn_deriv.f
  ${CMAKE_CURRENT_LIST_DIR}/d1mach.f
  ${CMAKE_CURRENT_LIST_DIR}/drc3jj.f
  ${CMAKE_CURRENT_LIST_DIR}/LRCLSDA.f
  ${CMAKE_CURRENT_LIST_DIR}/lyp_part_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/M06c_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/M06x_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/M08M11c_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/M08M11x_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/pw91_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/pw_lda.f
  ${CMAKE_CURRENT_LIST_DIR}/pz_lda.f
  ${CMAKE_CURRENT_LIST_DIR}/relxc.f
  ${CMAKE_CURRENT_LIST_DIR}/revpbe_x_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/rpbe_x_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/stvxc_spin.f
  ${CMAKE_CURRENT_LIST_DIR}/VS98c_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/VS98x_mod.f
  ${CMAKE_CURRENT_LIST_DIR}/xpbe_partial_derivs.f
  ${CMAKE_CURRENT_LIST_DIR}/ylm.f
  ${CMAKE_CURRENT_LIST_DIR}/ylm_real.f)

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fparser.f90
  ${CMAKE_CURRENT_LIST_DIR}/nlcorr_routines.f90
  ${CMAKE_CURRENT_LIST_DIR}/meta_scan.f90
  ${CMAKE_CURRENT_LIST_DIR}/tpss.f90
  ${CMAKE_CURRENT_LIST_DIR}/vdwkernel.f90)
