# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/elpa1_auxiliary_stub.f90)

if (ELPA_MT)
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa1_mt.f90
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_mt.f90)
else()
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa1_st.f90
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_st.f90)
endif()

if (ARCHITECTURE STREQUAL "AMD64_AVX")
  if (NOT USE_C_FILES OR NOT USE_CXX_FILES)
    message(FATAL_ERROR
      " USE_C_FILES and USE_CXX_FILES must be enabled for AMD64_AVX")
  endif()
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_real_sse-avx_2hv.c
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_complex_sse-avx_1hv.cpp
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_real_stub.f90
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_complex_stub.f90)
elseif (ARCHITECTURE STREQUAL "AMD64_SSE")
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_asm_x86_64.s)
elseif (ARCHITECTURE STREQUAL "PGI")
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_real.f90
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_complex.f90)
elseif (ARCHITECTURE STREQUAL "BlueGene_Q")
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_bgq.f90
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_complex.f90)
else()
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_real.f90
    ${CMAKE_CURRENT_LIST_DIR}/elpa2_kernels_complex.f90)
endif()
