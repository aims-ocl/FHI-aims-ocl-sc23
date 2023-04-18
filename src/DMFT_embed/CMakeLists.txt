# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/dmft_para.f90
  )
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/analy_continue_green_fn_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/analy_continue_self_energy_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/check_norm_dmft_p2_cluster.f90
  ${CMAKE_CURRENT_LIST_DIR}/check_norm_dmft_p2.f90
  ${CMAKE_CURRENT_LIST_DIR}/check_norm_dmft_p2_spectrum.f90
  ${CMAKE_CURRENT_LIST_DIR}/check_the_error_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_xc_pot.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagonalize_auxmat_lapack_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagonalize_green_fn_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagonalize_green_fn_dmft_GW.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_hartree_and_exchange_embed.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_hartree_and_exchange_embed_PBE0.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_delta_density_matr.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_density_matr.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_density_matr_no_ovlp.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_density_matr_test.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_embedding_self_enrg.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_embedding_self_enrg_no_ovlp.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_embed_gf_freq.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_embed_gf_freq_no_ovlp.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_GW_self_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_NEW_embed_gf_freq.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_NEW_on_site_gf_freq.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_NEW_on_site_gf_freq_PBE0.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_on_site_gf_freq.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_on_site_gf_freq_no_ovlp.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_on_site_gf_freq_PBE0.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_particle_number.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_pbc_and_cluster_quantities.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_qp_spectrum_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_spectral_func.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_spectrum_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_spectrum_dmft_GW.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_spectrum_k_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_total_energy_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_total_energy_dmft_PBE0.f90
  ${CMAKE_CURRENT_LIST_DIR}/initialize_fitting_parameter_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/invert_overlap_matrix_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/invert_overlap_matrix_sqrt.f90
  ${CMAKE_CURRENT_LIST_DIR}/mr_min_lsq_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/multiply_ovlp_matr.f90
  ${CMAKE_CURRENT_LIST_DIR}/multiply_ovlp_self_enrg.f90
  ${CMAKE_CURRENT_LIST_DIR}/overlap_matrix_sqrt.f90
  ${CMAKE_CURRENT_LIST_DIR}/power_genmat_lapack_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/prepare_partition_tabs_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/qp_spectrum_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/self_consistent_DMFT.f90
  ${CMAKE_CURRENT_LIST_DIR}/self_consistent_DMFT_PBE0.f90
  ${CMAKE_CURRENT_LIST_DIR}/spectrum_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/spectrum_dmft_GW.f90
  ${CMAKE_CURRENT_LIST_DIR}/transform_to_KS_basis_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/transform_to_NAO_basis_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_chemical_pot_cluster.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_chemical_pot_dmft_cluster.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_chemical_pot_dmft.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_chemical_pot_dmft_spectrum.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_chemical_pot.f90
  ${CMAKE_CURRENT_LIST_DIR}/zeroin_dmft_cluster.f90
  ${CMAKE_CURRENT_LIST_DIR}/zeroin_dmft.f90)
