!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  This subroutine is the last opportunity to check for any
!!  inconsistencies between keywords before the computationally
!!  expensive part starts. Because the order of keywords in control.in
!!  or geometry.in is not fixed, such checks cannot be performed while
!!  reading the input. The rule of this subroutine is that no keyword
!!  in dimensions.f90 or runtime_choices.f90 should be modified!
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
subroutine check_consistency_of_keywords()

  use dimensions,      only: calculate_perturbative_soc, n_periodic, n_spin, &
       & use_confined, use_hartree_fock, use_meta_gga, &
       & use_gga
  use MR_global,       only: active_nuclei, no_giao, j_is_active, &
       & shield_is_active, magnet_is_active, efg_is_active
  use mpi_tasks,       only: aims_stop, mpi_comm_global, myid
  use runtime_choices, only: atomic_solver, ATOMIC_SOLVER_ATOM_SPHERE, &
       & flag_rel, magnetic_response, mr_experimental, packed_matrix_format, &
       & REL_none, REL_atomic_zora, use_load_balancing, use_local_index, &
       & RI_type, RI_LVL
  use species_data,    only: include_min_basis
  use tools,           only: libxc_yes_no, newunit

  implicit none

  character(*), parameter :: THIS_SUB='check_consistency_of_keywords'
  integer :: i_status

  magnetic_response_keyword: if (magnetic_response) then
     if (.not. mr_experimental .and. myid == 0) then
        if (size(active_nuclei) <= 1 .and. j_is_active) &
             & call aims_stop('At least two atoms must be marked active in &
             &geometry.in for the calculation of J-couplings.', THIS_SUB)
        if (size(active_nuclei) == 0 .and. shield_is_active) &
             & call aims_stop('At least one atom must be marked active in &
             &geometry.in for the calculation of NMR shieldings.', THIS_SUB)
        if (size(active_nuclei) == 0 .and. efg_is_active) &
             & call aims_stop('At least one atom must be marked active in &
             &geometry.in for the calculation of the electric field &
             &gradient.', THIS_SUB)
        if (j_is_active .and. use_hartree_fock .and. RI_type >= RI_LVL) &
             & call aims_stop('LVL method of resolution of identity currently &
             &not supported for magnetic response.', THIS_SUB)
        if (shield_is_active .and. use_hartree_fock .and. .not. no_giao) &
             & call aims_stop('Hartree-fock currently not supported for &
             &the calculation of GIAO shieldings.', THIS_SUB)
        if (magnet_is_active .and. use_hartree_fock .and. .not. no_giao) &
             & call aims_stop('Hartree-fock currently not supported for &
             &the calculation of GIAO magnetizabilities.', THIS_SUB)
        if (use_meta_gga) &
             & call aims_stop('Meta functionals currently not supported for &
             &magnetic response.', THIS_SUB)
        if (n_periodic /= 0) &
             & call aims_stop('Periodic systems currently not supported for &
             &magnetic response.', THIS_SUB)
        if (packed_matrix_format /= 0 .and. .not. use_local_index) &
             & call aims_stop('Current calculation is set to use packed &
             &matrices. However, using packed matrices without use_local_index &
             &has proven to be inaccurate for magnetic response &
             &calculations.', THIS_SUB)
        if (flag_rel /= REL_none .and. flag_rel /= REL_atomic_zora) &
             & call aims_stop('flag_rel has incorrect value for magnetic &
             &response calculations.', THIS_SUB)
        if (flag_rel == REL_none .and. j_is_active .and. &
             & any(include_min_basis) .and. &
             & atomic_solver /= ATOMIC_SOLVER_ATOM_SPHERE) &
             & call aims_stop('Nonrelativistic J-coupling calculations are &
             &only available with the atom_sphere atomic solver.', THIS_SUB)
        if (flag_rel == REL_atomic_zora .and. j_is_active) &
             & call aims_stop('flag_rel=atomic_zora is currently not supported &
             &for the calculation of J-couplings.', THIS_SUB)
        if (flag_rel == REL_atomic_zora .and. shield_is_active) &
             & call aims_stop('flag_rel=atomic_zora is currently not supported &
             &for the calculation of NMR shieldings.', THIS_SUB)
        if (flag_rel == REL_atomic_zora .and. magnet_is_active) &
             & call aims_stop('flag_rel=atomic_zora is currently not supported &
             &for the calculation of magnetizability.', THIS_SUB)
        if (flag_rel == REL_atomic_zora .and. efg_is_active) &
             & call aims_stop('flag_rel=atomic_zora is currently not supported &
             &for the calculation of the electric field gradient.', THIS_SUB)
        if (n_spin == 2 .and. j_is_active .and. use_hartree_fock) &
             & call aims_stop('Open-shell systems (n_spin=2) currently not &
             &supported for the calculation of J-couplings with Hartree-Fock &
             &or a hybrid functional.', THIS_SUB)
        if (n_spin == 2 .and. efg_is_active) &
             & call aims_stop('Open-shell systems (n_spin=2) currently not &
             &supported for the calculation of the electric field gradient.', &
             & THIS_SUB)
        if (.not. libxc_yes_no() .and. use_gga) &
             & call aims_stop('This calculation only accepts LDA. Compile with &
             &LibXC for other functionals.', THIS_SUB)
        if (use_load_balancing .and. .not. use_local_index) &
             & call aims_stop('Setting load_balancing without use_local_index &
             &not supported for magnetic response calculations', THIS_SUB)
        if (flag_rel == REL_none .and. use_confined .and. j_is_active) &
             & call aims_stop('confined basis functions not suppored for &
             &nonrelativistic J-couplings', THIS_SUB)
        if (calculate_perturbative_soc) &
             & call aims_stop('Magnetic response calculations are currently &
             &unavailable with spin-orbit coupling.', THIS_SUB)
     end if
     call mpi_barrier(mpi_comm_global, i_status)
  end if magnetic_response_keyword
end subroutine check_consistency_of_keywords
