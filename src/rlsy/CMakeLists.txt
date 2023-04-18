# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

# OH: I have not the faintest idea if this is correct by any means. Seems to work.
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_constants.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_memtracker.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_helpers.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_mpi_helper.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_timer.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_scalapack_helper.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_sorting.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_geometry.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_linalg.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_distancetable.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_voronoi.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_verletlist.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_symmetry_helper_functions.f90
  )

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_interface.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_crystalstructure.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_spacegroup.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_basis_set.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_free_atom.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_realspace_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_extended_cluster.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_integration_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_integration_grid_helpers.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_partition_function.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_electron_density.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_kpointmesh.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_kspace_eigenproblem.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_realspace_integration_helper.f90
  )

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/test_symmetry_of_grid_and_things_on_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_custom_batches.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_index_map.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_pair_locator.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_calculate_density.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_calculate_overlap.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_hartree_potential.f90
  ${CMAKE_CURRENT_LIST_DIR}/rlsy_solve_KS_equations.f90
  ${CMAKE_CURRENT_LIST_DIR}/batch_rlsy_interface.f90
  )

