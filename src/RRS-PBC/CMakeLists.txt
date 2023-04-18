# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/read_rrs_pbc_unit_cell.f90
  ${CMAKE_CURRENT_LIST_DIR}/read_rrs_pbc_plot_band.f90
  )
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/build_rrs_pbc_k_matrices.f90
  ${CMAKE_CURRENT_LIST_DIR}/calculate_rrs_pbc_k_matrices.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagnalize_rrs_pbc.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_rrs_pbc_atom2basis.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_rrs_pbc_cell_index.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_rrs_pbc_cell_vector.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_rrs_pbc_k_point_list.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_rrs_pbc_lattice_vector.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_rrs_pbc_matrix_2D.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_rrs_pbc_matrix_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_rrs_pbc_matrix_real.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_rrs_pbc_result.f90
  ${CMAKE_CURRENT_LIST_DIR}/output_rrs_pbc_unit_cell.f90
  ${CMAKE_CURRENT_LIST_DIR}/out_rrs_pbc_band.f90
  ${CMAKE_CURRENT_LIST_DIR}/parse_rrs_pbc.f90
  ${CMAKE_CURRENT_LIST_DIR}/parse_rrs_pbc_rationality.f90
  ${CMAKE_CURRENT_LIST_DIR}/prepare_rrs_pbc_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/purge_rrs_pbc_hamiltonian.f90
  ${CMAKE_CURRENT_LIST_DIR}/read_rrs_pbc_lattice_vector.f90
  ${CMAKE_CURRENT_LIST_DIR}/rrs_pbc_complex_lapack_solver_fast.f90
  ${CMAKE_CURRENT_LIST_DIR}/rrs_pbc.f90
  ${CMAKE_CURRENT_LIST_DIR}/run_rrs_pbc.f90
  ${CMAKE_CURRENT_LIST_DIR}/store_rrs_pbc_data.f90)
