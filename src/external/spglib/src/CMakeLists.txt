# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

if (USE_SPGLIB)
  if (NOT USE_C_FILES)
    message(FATAL_ERROR "USE_C_FILES must be enabled for Spglib")
  endif()
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/arithmetic.c
    ${CMAKE_CURRENT_LIST_DIR}/cell.c
    ${CMAKE_CURRENT_LIST_DIR}/debug.c
    ${CMAKE_CURRENT_LIST_DIR}/delaunay.c
    ${CMAKE_CURRENT_LIST_DIR}/hall_symbol.c
    ${CMAKE_CURRENT_LIST_DIR}/kgrid.c
    ${CMAKE_CURRENT_LIST_DIR}/kpoint.c
    ${CMAKE_CURRENT_LIST_DIR}/mathfunc.c
    ${CMAKE_CURRENT_LIST_DIR}/niggli.c
    ${CMAKE_CURRENT_LIST_DIR}/pointgroup.c
    ${CMAKE_CURRENT_LIST_DIR}/primitive.c
    ${CMAKE_CURRENT_LIST_DIR}/refinement.c
    ${CMAKE_CURRENT_LIST_DIR}/sitesym_database.c
    ${CMAKE_CURRENT_LIST_DIR}/site_symmetry.c
    ${CMAKE_CURRENT_LIST_DIR}/spacegroup.c
    ${CMAKE_CURRENT_LIST_DIR}/spg_database.c
    ${CMAKE_CURRENT_LIST_DIR}/spglib.c
    ${CMAKE_CURRENT_LIST_DIR}/spglib_f.c
    ${CMAKE_CURRENT_LIST_DIR}/spin.c
    ${CMAKE_CURRENT_LIST_DIR}/symmetry.c
    ${CMAKE_CURRENT_LIST_DIR}/test.c
    ${CMAKE_CURRENT_LIST_DIR}/tetrahedron_method.c)
endif()
