# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

add_subdirectory(atom_sphere)
add_subdirectory(external)
add_subdirectory(sratom)
add_subdirectory(rdirac)

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/species_data.f90
  ${CMAKE_CURRENT_LIST_DIR}/grids.f90
  ${CMAKE_CURRENT_LIST_DIR}/spline.f90
  ${CMAKE_CURRENT_LIST_DIR}/basis.f90
  ${CMAKE_CURRENT_LIST_DIR}/free_atoms.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_leb_grids.f90
  ${CMAKE_CURRENT_LIST_DIR}/divide_angular_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/divide_angular_grid_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/atomic_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/atom_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/atom_sphere_wrapper.f90
  ${CMAKE_CURRENT_LIST_DIR}/confined_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/cutoff_pot.f90
  ${CMAKE_CURRENT_LIST_DIR}/fixed_basis_potential.f90
  ${CMAKE_CURRENT_LIST_DIR}/free_atoms_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/gaussian_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_bas_indexing.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_full_bas.f90
  ${CMAKE_CURRENT_LIST_DIR}/generate_Lsp_indexing.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_angular_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_bas_dimensions.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_confined_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_dirac_hydrogenic_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_fixed_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_free_atoms.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_free_atoms_polarized.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_free_ion.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_gaussian_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_hydrogenic_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_ionic_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_species_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_sph_gaussian_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/hydrogenic_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/increment_ylm.f90
  ${CMAKE_CURRENT_LIST_DIR}/insertionsort.f90
  ${CMAKE_CURRENT_LIST_DIR}/ion_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/orthonormalize_basis_fn.f90
  ${CMAKE_CURRENT_LIST_DIR}/rel_atoms.f90
  ${CMAKE_CURRENT_LIST_DIR}/shrink_fixed_basis_phi_thresh.f90
  ${CMAKE_CURRENT_LIST_DIR}/shrink_fixed_large_basis.f90
  ${CMAKE_CURRENT_LIST_DIR}/tab_local_trigonom.f90
  ${CMAKE_CURRENT_LIST_DIR}/tab_local_ylm.f90
  ${CMAKE_CURRENT_LIST_DIR}/use_atomic_solver_basis_fns_as_species_basis_fns.f90
  ${CMAKE_CURRENT_LIST_DIR}/verify_onsite_integrals.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_sto_basis_fns.f90)
