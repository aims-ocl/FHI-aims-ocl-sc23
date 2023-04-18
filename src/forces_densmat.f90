!****h* FHI-aims/forces_densmat
!  NAME
!    forces_densmat
!  SYNOPSIS
module forces_densmat
!  PURPOSE
!    The module calculates a subset of the forces using the density matrix
!    formalism.  This works with the both periodic and nonperiodic systems.
!    The contributions included are:
!      o Pulay forces
!      o GGA forces
!      o meta-GGA forces
!      o atomic ZORA correction term to forces
!      o Gnomf forces
!    Contributions which are calculated elsewhere (as of 2018 January 8th)
!    include:
!      o Hellman-Feynman forces (sum_up_whole_potential_p1.f90)
!      o Electrostatic multipole derivatives (sum_up_whole_potential_p1.f90)
!      o van der Waals forces (various implementations exist in aims)
!      o Exact exchange forces
!  USES
  implicit none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!******

contains

!-------------------------------------------------------------------------------
!****s* forces_densmat/update_forces_densmat
!  NAME
!    update_forces_densmat
!  SYNOPSIS
  subroutine update_forces_densmat &
       ( KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, &
       partition_tab, rho, rho_gradient, kinetic_density, hartree_potential, &
       basis_l_max, dens_mat, sum_forces )
!  PURPOSE
!    The subroutine updates the Pulay- and gga-forces using the density matrix
!    type of method.
!  USES
    use aims_gpu, only: gpu_forces_used
    use analytical_stress, only: AS_pulay_stress_local, &
        AS_dde_T_sa_at_di_orb_local, AS_pulay_stress, &
        AS_dde_T_sa_at_di_orb, AS_dde_potential
    use constraint, only: n_active_regions, constraint_potential
    use density_matrix_evaluation, only: evaluate_densmat
    use dimensions, only: n_k_points_task, n_full_points, n_species, n_atoms, &
        use_constraint, n_basis, n_states, n_spin, n_k_points, &
        n_hamiltonian_matrix_size, compute_heat_flux
    use elsi_wrapper, only: eh_scf,aims_elsi_get_edm
    use load_balancing, only: batch_perm, use_batch_permutation, &
        init_comm_full_local_matrix
    use localorb_io, only: use_unit, OL_norm, localorb_info
    use lpb_solver_utilities, only: Gnonmf_forces, use_Gnonmf_forces
    use mpi_tasks, only: myid
    use pbc_lists, only: k_weights
    use runtime_choices, only: AS_components, use_analytical_stress, &
        use_elsi_dm, use_symmetry_reduced_spg, AS_components, real_eigenvectors, &
        use_rlsy_symmetry
    use scalapack_wrapper, only: ham, ham_complex, eigenvec, eigenvec_complex, &
        my_k_point
    use sym_base, only: evaluate_densmat_sym
    use timing, only: gpu_forces_always_used
    use synchronize_mpi, only: sync_matrix
    use heat_flux
    use rlsy_interface, only: rlsy_h
    implicit none
!  ARGUMENTS
    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task), intent(in)  :: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task), intent(in)  :: KS_eigenvector_complex
    real*8,     dimension(n_states, n_spin, n_k_points),              intent(in)  :: occ_numbers
    real*8,     dimension(n_states, n_spin, n_k_points),              intent(in)  :: KS_eigenvalue
    real*8,     dimension(n_full_points),                             intent(in)  :: partition_tab
    real*8,     dimension(n_spin, n_full_points),                     intent(in)  :: rho
    real*8,     dimension(3, n_spin, n_full_points),                  intent(in)  :: rho_gradient
    real*8,     dimension(n_spin, n_full_points),                     intent(in)  :: kinetic_density
    real*8,     dimension(n_full_points),                             intent(in)  :: hartree_potential
    integer,    dimension(n_species),                                 intent(in)  :: basis_l_max
    real*8,     dimension(*),                                         intent(out) :: dens_mat
    real*8,     dimension(3, n_atoms),                                intent(out) :: sum_forces
!  INPUTS
!    o KS_eigenvector -- eigenvectors if real lapack eigenvectors are in use
!    o KS_eigenvector_complex  -- eigenvectors if complex lapack eigenvectors
!      are in use
!    o occ_numbers -- occupation of differenct eigenstates
!    o KS_eigenvalue -- eigenvalues
!    o partition_tab -- values of partition function
!    o rho -- electron density
!    o rho_gradient -- gradient of electron density
!    o kinetic_density -- kinetic-energy density
!    o hartree_potential -- hartree potential
!    o basis_l_max -- maximum l of basis functions
!  OUTPUT
!    o dens_mat -- work space for the density matrix.  Output is basically
!      noise and should be ignored.
!      When this routine is called, dens_mat has either the dimension
!      (n_hamiltonian_matrix_size, n_spin) or (n_local_matrix_size, n_spin)
!      so we declare it here as a 1D assumed size array.
!    o sum_forces -- the sum of Pulay, gga, meta-GGA, and gnonmf forces
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    integer :: i_term, i_k_point, i_state, i_spin
    integer :: i_coord_1, i_coord_2, i_counter
    real*8  :: occu, this_en
    logical :: dens_mat_is_en_weighted

    real*8, dimension(n_states, n_spin, n_k_points) :: aux_eigenvalue
    real*8, dimension(n_states, n_spin, n_k_points) :: scaled_occs
    real*8, dimension(1,1) :: dummy_mat

    ! Load balancing stuff
    integer :: n_bp ! Number for current batch permutation used
    integer :: ld_dens_mat ! leading dimension of density matrix in calling
                           ! routine
    integer :: i_off ! offset to account for spin channel indexing, as dens_mat
                     ! is passed around as a 1D assumed size array, whereas it
                     ! is actually a 2D array

        !unfold_stress: block
        real*8 :: AS_pulay_stress_sum(3,3)
        real*8 :: AS_dde_T_sa_at_di_orb_sum(3,3)
        integer :: i_op, i_atom, j_atom

    call localorb_info("", use_unit, '(2X,A)', OL_norm)

    n_bp = use_batch_permutation
    if (n_bp > 0) then
      call localorb_info("Evaluating density-matrix-based force terms: &
                         &batch-based integration with load balancing", &
                         use_unit,'(2X,A)',OL_norm)
      ld_dens_mat = batch_perm(n_bp)%n_local_matrix_size

      ! evaluate_densmat uses the BLACS-distributed Hamiltonian matrix as a work
      ! array to calculate the density matrix, so we need to reinitialize
      ! index arrays for local_index communication
      call init_comm_full_local_matrix( &
           batch_perm(n_bp)%n_basis_local, &
           batch_perm(n_bp)%i_basis_local )
    else
      call localorb_info("Evaluating density-matrix-based force terms: &
                         &batch-based integration", &
                         use_unit,'(2X,A)',OL_norm)
      ld_dens_mat = n_hamiltonian_matrix_size
    end if

    if(.not. use_elsi_dm) then
       call check_occs('forces_densmat', occ_numbers, .true.)
    endif

    sum_forces = 0.0d0
    ! Initialization of stress components
    if (use_analytical_stress) then
      ! CC:  We have to compute the following terms:
      !#(1)## - Non-relativistic case
      !       sum_ijl c_ijl   d/de_LM  < i | h - e_l | j >
      !     = sum_ijl c_ijl [        2 * < d/de_LM i | h - e_l       | j >
      !                         + \delta_kron_LM < i | h - e_l       | j >
      !                         +                < i | d/de_LM v_H   | j >
      !                         +                < i | d/dr_L d/dr_M | j > ]
      !  (A) The 'Pulay term' < d/de_LM i | h - e_l | j > is computed in
      !      AS_pulay_stress, whereby
      !           d/de_LM i(r) =  d/dR_L i(r) * ( R - r )_M
      !  (B) The 'Jacobian' term
      !           \delta_kron_LM < i | h - e_l | j > = < i | en_dens_xc(r) - v_xc(r) | j >
      !      is recomputed for concistency reasons together with the
      !  (C) integral over the strain derivative of the hartree potential
      !      (computed in sum_up_whole_potential)
      !           < i | d/de_LM v_H   | j >
      !  (D) and the 'hessian term'
      !           < i | d/dr_L d/dr_M | j >
      !  (B) (C) and (D) are only calculated for i_term=1
      !
      !  (E) On top of that, small radial grids lead to numerical error in the
      !      strain derivative of the kinetic energy for the on-site terms. We
      !      correct for them by computing
      !           d/de_LM < i | T | j > =
      !           2 * < d/de_LM |  T | j > + \delta_kron_LM < i | T | j > + < i | d/dr_L d/dr_M | j >
      !      Since d/de_LM < i | T | j > should be 0 for on-site i,j, we
      !      substract the numerically computed d/de_LM < i | T | j > to achieve
      !      convergence at small radial grids.
      !
      !      Accordingly,
      !      On-site, i  = j -> AS_dde_T_sa_at_sa_orb, for debug purposes
      !      On-site, i != j -> AS_dde_T_sa_at_di_orb  (incl. AS_dde_T_sa_at_sa_orb)
      !      Off-site        -> AS_dde_T_di_at_di_orb, for debug purposes
      !
      !#(2)## - Relativistic case (atomic zora)
      !  In the relativistic case, the hamiltonian h contains only ``half''
      !  the kinetic energy
      !      T_ji =  <j| \laplace |i> --> <j| t_zora_i |i>
      !  the other ``half'' of the kinetic energy, i.e.,
      !       <i| t_zora_j |j>
      !  is treated in a second step, whereby the explicit derivative
      !       kinetic_wave_deriv = d/dr (t_zora_j |j>)
      !  is computed. Accordingly, we get
      !       AS_strain_deriv_kinetic_wave = d/de (t_zora_j |j>) = kinetic_wave_deriv * (r-R_j)
      !  The contributions of [ <i| d/de (t_zora_j |j>) ]^T are directly added
      !  to AS_pulay_stress.
      !  N.B.:
      !  In this formalism, in which the derivative of the kinetic energy is
      !  computed explicitly, no hessian (D) is required anymore. This removes
      !  the respective numerical errors, so that no on-site corrections (E)
      !  seem to be required anymore, either.  As a consequence, there is no
      !  need to compute
      !       AS_jac_pot_kin = < i | en_dens_xc(r) - v_xc(r) | j > \delta_lm + < i | d/de_LM v_H   | j >
      !  at all, since the remaining quantitites are available / more rapdidly
      !  computable in sum_up_whole_potential.
      !
      !#(3)## - GGA case
      !  The GGA term involves the derivative of the xc potential with respect
      !  to the square of the density gradient: d/dgrad(rho)^2 f_xc)
      !  We have to calculate the following term:
      !       sum_ijl c_ijl < i | (d/dgrad(rho)^2 f_xc) * (d/de_LM grad(rho)^2) | j >
      !  This term splits up into two contributuions:
      !  (A) 4 * grad(< i |) (d/dgrad(rho)^2 f_xc) * grad(rho) * (d/dr_L | j >)       * (r-R_M)
      !  (B) 4 *      < i |  (d/dgrad(rho)^2 f_xc) * grad(rho) * (d/dr_L grad(| j >)) * (r-R_M)
      !  There is a scalar product between grad(rho) and the gradient of the
      !  wave function.

      AS_pulay_stress_local(:)              = 0.0d0
      AS_dde_T_sa_at_di_orb_local(:)        = 0.0d0
      ! Allocate Heat Flux
      if (compute_heat_flux) then
        if (.not.allocated(HF_stress_per_atom_PU)) allocate(HF_stress_per_atom_PU(AS_components,n_atoms))
        HF_stress_per_atom_PU(:,:)          = 0.0d0
      end if


    end if

    if(use_Gnonmf_forces) then
      ! Gnonmf forces are summed into sum_forces in this routine
      ! As a result, this variable is unused and needs zeroing.
      Gnonmf_forces = 0d0
    end if

    ! First, apply possible shift to the eigenvectors in case a
    ! fixed-spin-moment calculation was done
    if (use_constraint) then
      if (n_active_regions.eq.1) then
        do i_spin=1,n_spin,1
          aux_eigenvalue(:,i_spin,:) = &
               KS_eigenvalue(:,i_spin,:) - constraint_potential(1,i_spin)
        enddo
      end if
    else
      aux_eigenvalue = KS_eigenvalue
    end if

    ! There are two types of terms entering into the forces: those which require
    ! the energy-weighted density matrix and those requiring the "ordinary"
    ! density matrix
    ! We here calculated both, one after the other, and sum them together. The
    ! dens_mat_is_en_weighted flag is passed to
    ! integrate_force_integrals_densmat indicate which term we're presently
    ! calculating.
    do i_term = 1,2
      dens_mat(1:ld_dens_mat*n_spin) = 0.d0

      if (i_term.eq.1) then
        ! For the first term, we need the density matrix
        ! This term includes the first term for the Pulay forces (Equation 73
        ! of Blum et al., Comput Phys Comm 2009.)
        ! GGA, meta-GGA, and Gnonmf forces are also included in this term
        ! Analytic stress is also handled by this term
        dens_mat_is_en_weighted = .false.
      else
        ! For the second term, we need the energy-weighted density matrix
        ! This term is unique to the Pulay forces, where it is the second term
        ! (Equation 73 of Blum et al., Comput Phys Comm 2009.)
        dens_mat_is_en_weighted = .true.
      end if

      ! Construct the needed flavor of the density matrix (energy-weighted or
      ! not)
      do i_spin = 1, n_spin
        do i_k_point = 1, n_k_points
          do i_state = 1, n_states
            this_en = aux_eigenvalue(i_state, i_spin, i_k_point)
            occu = occ_numbers(i_state,i_spin, i_k_point)

            if (dens_mat_is_en_weighted) then
              scaled_occs(i_state, i_spin, i_k_point) = -2.0d0 * occu * this_en
            else
              scaled_occs(i_state, i_spin, i_k_point) =  2.0d0 * occu
            endif
          end do
        end do

        if(use_elsi_dm) then
          ! DM/EDM stored in eigenvec(_complex)
          if (real_eigenvectors) then
            if(dens_mat_is_en_weighted) then
              call aims_elsi_get_edm(eh_scf,eigenvec(:,:,i_spin))

              ham(:,:,i_spin) = &
                   k_weights(my_k_point)*eigenvec(:,:,i_spin)*(-2.0d0)
            else
              ham(:,:,i_spin) = &
                   k_weights(my_k_point)*eigenvec(:,:,i_spin)*2.0d0
            endif
          else
            if(dens_mat_is_en_weighted) then
              call aims_elsi_get_edm(eh_scf,eigenvec_complex(:,:,i_spin))

              ham_complex(:,:,i_spin) = &
                   k_weights(my_k_point)*eigenvec_complex(:,:,i_spin)*(-2.0d0)
            else
              ham_complex(:,:,i_spin) = &
                   k_weights(my_k_point)*eigenvec_complex(:,:,i_spin)*2.0d0
            endif
          endif
        endif

        i_off = (i_spin-1)*ld_dens_mat
        if (use_symmetry_reduced_spg) then
          call evaluate_densmat_sym(KS_eigenvector, KS_eigenvector_complex, &
               scaled_occs, dummy_mat, dens_mat(1+i_off), i_spin, .true.)
        else
          call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, &
               scaled_occs, dummy_mat, dens_mat(1+i_off), i_spin, .true.)
        endif
      end do

      !! CC: Debug routines to read/write the Pulay matrix
      !! if ( (AS_flag_write_pulay_dm) .and. (myid==0) ) then
      !!   call AS_write_pulay_dm( dens_mat, ld_dens_mat, &
      !!        n_spin, i_term )
      !! end if
      !! if (AS_flag_read_pulay_dm) then
      !!   call AS_read_pulay_dm( dens_mat, ld_dens_mat, &
      !!        n_spin, i_term )
      !! end if

      ! Now perform batch integration for the currently selected term and add
      ! the result into sum_forces
      call integrate_force_integrals_densmat(hartree_potential, rho, &
           rho_gradient, kinetic_density, partition_tab, basis_l_max, &
           dens_mat_is_en_weighted, dens_mat, sum_forces)
    end do ! i_term

    ! Map from one index (1 to 6 or 9) to 2 indices (1 to 3, 1 to 3)
    ! 1 2 3
    ! 7 4 5
    ! 8 9 6
    ! We do not symmetrize here if AS_components=6. This is only allowed if we
    ! take the total analytical stress.
    if (use_analytical_stress) then
      i_counter=0
      do i_coord_1=1,3,1
        do i_coord_2=i_coord_1,3,1
          i_counter=i_counter+1
          AS_pulay_stress( i_coord_1 , i_coord_2 ) = &
               AS_pulay_stress_local( i_counter )
          AS_dde_T_sa_at_di_orb( i_coord_1 , i_coord_2 ) = &
               AS_dde_T_sa_at_di_orb_local( i_counter )

          if (i_coord_1.ne.i_coord_2) then
            AS_pulay_stress( i_coord_2 , i_coord_1 )       = 0.0d0
            AS_dde_T_sa_at_di_orb( i_coord_2 , i_coord_1 ) = 0.0d0
          end if
        end do
      end do
      if (AS_components .eq. 9) then
        AS_pulay_stress(2,1) = AS_pulay_stress_local(7)
        AS_pulay_stress(3,1) = AS_pulay_stress_local(8)
        AS_pulay_stress(3,2) = AS_pulay_stress_local(9)
        AS_dde_T_sa_at_di_orb(2,1) = AS_dde_T_sa_at_di_orb_local(7)
        AS_dde_T_sa_at_di_orb(3,1) = AS_dde_T_sa_at_di_orb_local(8)
        AS_dde_T_sa_at_di_orb(3,2) = AS_dde_T_sa_at_di_orb_local(9)
      end if
      if (allocated(AS_dde_potential)) deallocate(AS_dde_potential)
      ! YY symmetry
      if (use_rlsy_symmetry) then
        !unfold_stress: block
        AS_pulay_stress_sum = 0.0d0
        AS_dde_T_sa_at_di_orb_sum = 0.0d0
        do i_op = 1, rlsy_h%spacegroup%n_operation
          AS_pulay_stress_sum(:,:) = AS_pulay_stress_sum(:,:) + &
              matmul( &
              matmul(rlsy_h%spacegroup%op(i_op)%m(:,:), &
                     AS_pulay_stress(:,:)), &
                     transpose(rlsy_h%spacegroup%op(i_op)%m(:,:)))
          AS_dde_T_sa_at_di_orb_sum(:,:) = AS_dde_T_sa_at_di_orb_sum(:,:) + &
              matmul( &
              matmul(rlsy_h%spacegroup%op(i_op)%m(:,:), &
                     AS_dde_T_sa_at_di_orb(:,:)), &
                     transpose(rlsy_h%spacegroup%op(i_op)%m(:,:)))
        end do
        AS_pulay_stress = AS_pulay_stress_sum
        AS_dde_T_sa_at_di_orb = AS_dde_T_sa_at_di_orb_sum
        !end block unfold_stress
      end if
      ! YY symmetry
    end if

   ! CC: Sync per atom heat flux
   if (compute_heat_flux) then
     call sync_matrix( HF_stress_per_atom_PU, AS_components, n_atoms )
   end if

  end subroutine update_forces_densmat
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/integrate_force_integrals_densmat
!  NAME
!    integrate_force_integrals_densmat
!  SYNOPSIS
  subroutine integrate_force_integrals_densmat( hartree_potential_std, rho_std,&
       rho_gradient_std, kinetic_density_std, partition_tab_std, basis_l_max, &
       dens_mat_is_en_weighted, dens_mat, sum_forces )
!  PURPOSE
!    The subroutine is called by integrate_force_integrals_densmat and it
!    calculates the Pulay- and gga-force components after
!    integrate_force_integrals_densmat is first constructed the density matrix
!    type of matrixes. For calculate whole Pulay+gga forces, this routine is
!    called 2 times. This is because the second part of Pulay forces (psi*psi)
!    needs energies included in the density matrix and the first one
!    do not. This means two different types of density matrix.
!  USES
    use aims_gpu, only: output_mem_array_gpu, gpu_forces_used
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use analytical_stress
    use basis, only: basis_deriv_ordered, basis_kinetic_ordered, &
        basis_wave_ordered
    use cartesian_ylm, only: n_max_cartesian, evaluate_cartesians, &
        evaluate_cartesian_hessian_and_gradient_terms
    use c_helper, only: fort_log_to_c_int
    use constants, only: light_speed_sq
    use dimensions
    use grids, only: batches, batch_of_points
    use load_balancing, only: use_batch_permutation, batch_perm, &
        get_batch_weights, set_batch_weights, permute_point_array
    use localorb_io, only: OL_norm, localorb_info, use_unit
    use lpb_solver_utilities, only: use_Gnonmf_forces, &
        evaluate_Gnonmf_forces_dens_mat,evaluate_Gnonmf
    use mpi_tasks, only: aims_stop, check_allocation, myid
    use pbc_lists, only: centers_basis_integrals, inv_centers_basis_integrals, &
        cbasis_to_atom
    use runtime_choices, only: AS_components, flag_rel, packed_matrix_format, &
        prune_basis_once, REL_atomic_zora, REL_none, use_analytical_stress, &
        use_gpu, use_gpu_forces, PM_none, use_rlsy_symmetry
    use timing, only: gpu_forces_always_used
    use rlsy_interface, only: rlsy_h
    use batch_rlsy_interface, only: batch_irr, &
                                    compute_batch_rlsy_full_to_irr_points, &
                                    compute_batch_rlsy_full_to_irr_points_gradient, &
                                    compute_batch_rlsy_full_to_irr
                                               
    implicit none
!  ARGUMENTS
    real*8,  target, intent(IN),    dimension(n_full_points)            :: hartree_potential_std
    real*8,  target, intent(IN),    dimension(n_spin, n_full_points)    :: rho_std
    real*8,  target, intent(IN),    dimension(3, n_spin, n_full_points) :: rho_gradient_std
    real*8,  target, intent(IN),    dimension(n_spin, n_full_points)    :: kinetic_density_std
    real*8,  target, intent(IN),    dimension(n_full_points)            :: partition_tab_std
    integer,         intent(IN),    dimension(n_species)                :: basis_l_max
    logical,         intent(IN)                                         :: dens_mat_is_en_weighted
    real*8,          intent(IN),    dimension(*)                        :: dens_mat
    real*8,          intent(INOUT), dimension(3, n_atoms)               :: sum_forces

    real*8,           dimension(3, n_atoms)                             :: sum_forces_temp

!  INPUTS
!    o hartree_potential_std -- Hartree potential
!    o rho_std -- electron density
!    o rho_gradient_std -- gradient of electron density (used if gga is used)
!    o kinetic_density_std -- gradient of kinetic energy density (used if meta-gga is used)
!    o partition_tab -- values of partition function
!    o basis_l_max -- maximum l of basis functions
!    o dens_mat -- The density matrix, possibly energy-weighted.
!      When this routine is called, dens_mat has either the dimension
!      (n_hamiltonian_matrix_size, n_spin) or (n_local_matrix_size, n_spin)
!      so we declare it here as a 1D assumed size array.
!    o dens_mat_is_en_weighted -- Whether the density matrix is energy-weighted.
!      If .false. the first part of Pulay forces, as well as various other
!      forces, are calculated; if .true., only second part (psi*psi) part of the
!      Pulay forces are calculated.
!  OUTPUT
!    o sum_forces -- the sum of Pulay, gga, meta-GGA, and gnonmf forces
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    real*8, dimension(n_spin) :: local_potential_parts

    integer :: l_ylm_max
    integer, dimension(:,:), allocatable :: index_lm
    real*8, dimension(:,:), allocatable :: ylm_tab

    real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
    real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

    real*8 coord_current(3)
    real*8 dist_tab(n_centers_integrals, n_max_batch_size)
    real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
    real*8 i_r(n_max_compute_atoms)
    real*8 dir_tab(3,n_centers_integrals, n_max_batch_size)
    real*8 trigonom_tab(4,n_max_compute_atoms)

    real*8,dimension(:)  ,allocatable:: radial_wave
    real*8,dimension(:)  ,allocatable:: radial_wave_deriv
    real*8,dimension(:)  ,allocatable:: kinetic_wave
    real*8,dimension(:)  ,allocatable:: kinetic_wave_deriv

    ! Note that, for the case where we're evaluated the term with the
    ! energy-weighted density matrix, H_times_psi just holds psi instead
    real*8,dimension(:,:,:),allocatable:: H_times_psi
    real*8,dimension(:,:,:,:),allocatable:: d_H_times_psi
    real*8,dimension(:,:)  ,allocatable:: wave
    real*8,dimension(:,:,:)  ,allocatable:: d_wave

    ! for XC
    real*8, dimension(n_max_batch_size) :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin, n_max_batch_size) :: local_xc_derivs
    real*8, dimension(3,n_spin,n_max_batch_size) :: xc_gradient_deriv
    real*8, dimension(n_spin, n_max_batch_size) :: xc_tau_deriv

    ! Auxiliary Hamiltonian matrix, to sum up contributions from only a single
    ! integration shell.  The hope is that such a separate treatment will allow
    ! to minimize numerical noise introduced through ZORA
    real*8, dimension(:), allocatable :: forces_shell

    ! optimal accounting for matrix multiplications: only use points with
    ! nonzero components
    integer :: n_points
    integer :: n_rel_points

    ! and condensed version of hamiltonian_partition_tabs on angular grids
    real*8 :: partition(n_max_batch_size)
    real*8 :: energy_partition(n_max_batch_size)

    real*8, dimension(:,:,:), allocatable :: gradient_basis_wave
    real*8, dimension(:,:,:), allocatable :: kinetic_gradient_basis_wave

    ! Following is all that is needed for the handling of ZORA scalar relativity
    real*8, dimension(n_spin) :: zora_operator
    logical :: t_zora(2)
    real*8, dimension(n_spin) :: zora_potential_parts

    real*8, dimension(:), allocatable :: dist_tab_full
    real*8, dimension(:,:), allocatable :: dir_tab_full_norm
    real*8, dimension(:), allocatable :: i_r_full

    real*8, dimension(:,:,:,:), allocatable :: zora_vector1
    real*8, dimension(:,:,:,:), allocatable :: zora_vector2

    ! This term contains contributions from the xc potential and the
    ! zora formalism (if applicable) which are summed up using Gauss' law:
    ! < grad(phi_i) | local_gradient_sum |grad(phi_j) >
    real*8, dimension(3,n_spin) :: sum_of_local_gradients

    ! for pruning of atoms, radial functions, and basis functions, to only the
    ! relevant ones ...
    integer :: n_compute_c, n_compute_a
    integer,dimension(:),allocatable :: i_basis

    integer :: n_compute_fns
    integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
    integer :: i_basis_fns_inv(n_basis_fns,n_centers)
    integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

    integer :: n_compute_atoms
    integer :: atom_index(n_centers_integrals)
    integer :: atom_index_inv(n_centers)

    integer :: spline_array_start(n_centers_integrals)
    integer :: spline_array_end(n_centers_integrals)

    ! for splitting of angular shells into "octants"
    integer division_low
    integer division_high

    ! counters
    integer :: i_basis_1
    integer :: i_basis_2
    integer :: i_compute_1
    integer :: i_atom, i_atom_2, i_dir
    integer :: i_grid
    integer :: i_index, i_l, i_m
    integer :: i_coord
    integer :: i_division
    integer :: i_species
    integer :: i_point
    integer :: i_full_points_C
    integer :: i_full_points_2C
    integer :: i_full_points_A
    integer :: i_spin, i_spin_2
    integer :: i_my_batch
    integer :: n_spin_local
    integer :: i_calculate_dimension, i_compute

    ! Stress
    integer :: AS_hessian_index, n_calculate_dimension
    integer, dimension(9,3) :: AS_index_hessian_gga

    ! GGA forces
    real*8, dimension(:,:,:), allocatable :: hessian_basis_wave
    real*8, dimension(:,:,:), allocatable :: sum_hessian
    real*8, dimension(:,:,:), allocatable :: cartesians
    real*8, dimension(:,:,:), allocatable :: sum_gradient
    real*8,dimension(:),allocatable:: radial_wave_2nd_deriv

    logical:: gga_forces_on_temp
    ! Let me check if meta_gga_forces are needed
    logical:: meta_gga_forces_on_temp
    integer :: info

    !MPB solvation
    logical :: Gnonmf_forces_on_temp
    real*8, dimension(:),allocatable :: local_Gnonmf_derivs
    real*8, dimension(:,:),allocatable :: Gnonmf_gradient_deriv

    ! CUDA
    logical :: gpu_save, index_on_gpu

    ! Variables for permuting compute basis elements so that the row/columns are
    ! ordered sequentially by associated atom.  This is needed for GPU
    ! acceleration to eliminate problems related to warp divergence.
    integer :: i_permute, my_atom
    integer, allocatable :: permute_compute_by_atom(:)
    integer, allocatable :: n_compute_for_atom(:)

    ! Load balancing stuff
    integer :: n_bp ! Number for current batch permutation used
    integer :: n_my_batches_work ! Number of batches actually used
    integer :: my_n_full_points ! Number of actual integration points
    type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches
                                                       ! actually used
    integer :: ld_dens_mat ! leading dimension of density matrix in calling
                           ! routine
    integer, dimension(:), allocatable :: ins_idx
    integer :: n_basis_local

    ! Pointers to the actually used array
    real*8, pointer :: partition_tab(:)
    real*8, pointer :: rho(:,:)
    real*8, pointer :: rho_gradient(:,:,:)
    real*8, pointer :: hartree_potential(:)
    real*8, pointer :: kinetic_density(:,:)

    real*8, pointer :: AS_dde_potential_work(:,:,:)

    ! Timing
    real*8, allocatable :: batch_times(:)
    real*8 :: time_start

    character*140 :: info_str

      !unfold_force: block
      !integer :: i_op, i_atom, j_atom
      integer :: i_op, j_atom

    sum_forces_temp = 0.0d0

    gpu_save = use_gpu

    if (use_gpu_forces .and. .not. use_gpu) then
      write(info_str,'(2X,A)')
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      write(info_str,'(2X,A)') &
           "You have request GPU acceleration for calculation of various &
           &contributions to the forces, but no GPU"
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      write(info_str,'(2X,A)') &
           "is available.  Turning off GPU acceleration for force &
           &calculations."
      call localorb_info(info_str,use_unit,'(A)',OL_norm)

      use_gpu = .false.
    end if

    if (use_gpu_forces .and. use_Gnonmf_forces .and. &
        .not.dens_mat_is_en_weighted) then
      write(info_str,'(2X,A)')
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      write(info_str,'(2X,A)') &
           "You have request GPU acceleration for calculation of Gnonmf &
           &contributions to the forces, but this"
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      write(info_str,'(2X,A)') &
           "is not supported.  Turning off GPU acceleration for force &
           &and stress calculations."
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      use_gpu = .false.
    end if

    use_gpu = use_gpu .and. use_gpu_forces

    if (use_gpu) then
      if (use_analytical_stress) then
        write(info_str,'(2X,A)') "GPU acceleration will be used when evaluting &
                                 &forces and analytical stress using the &
                                 &density matrix."
      else
        write(info_str,'(2X,A)') "GPU acceleration will be used when evaluting &
                                 &forces using the density matrix."
      end if
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      gpu_forces_used = .true.
    else
      gpu_forces_always_used = .false.
      gpu_forces_used = .false.
    end if

    !---------------------------------------------------------------------------
    ! Initialize load balancing:
    ! Set pointers either to permuted batches / arrays over integration points
    ! (for load balancing) or to standard batches / arrays (no load balancing)
    n_bp = use_batch_permutation

    if(n_bp > 0) then
      n_my_batches_work = batch_perm(n_bp)%n_my_batches
      my_n_full_points = batch_perm(n_bp)%n_full_points
      n_basis_local = batch_perm(n_bp)%n_basis_local
      ld_dens_mat = batch_perm(n_bp)%n_local_matrix_size

      allocate(ins_idx(n_basis_local))

      batches_work => batch_perm(n_bp)%batches
      partition_tab => batch_perm(n_bp)%partition_tab

      allocate(hartree_potential(my_n_full_points))
      call permute_point_array(n_bp, 1, hartree_potential_std, &
                               hartree_potential)

      allocate(rho(n_spin, my_n_full_points))
      call permute_point_array(n_bp, n_spin, rho_std, rho)

      allocate(rho_gradient(3, n_spin, my_n_full_points))
      call permute_point_array(n_bp, 3*n_spin, rho_gradient_std, rho_gradient)

      allocate(kinetic_density(n_spin, my_n_full_points))
      call permute_point_array(n_bp, n_spin, kinetic_density_std, &
                               kinetic_density)
    elseif (use_rlsy_symmetry) then
      n_my_batches_work = batch_irr%n_my_batches
      my_n_full_points = batch_irr%n_full_points
      ld_dens_mat = n_hamiltonian_matrix_size
      allocate(ins_idx(1))
      batches_work => batch_irr%batches

      !allocate(partition_tab(batch_irr%n_full_points))
      !call compute_batch_rlsy_full_to_irr_points_weighted(partition_tab_std(:), partition_tab(:))

      partition_tab =>batch_irr%partition_tab

      !call compute_batch_rlsy_full_to_irr(hartree_potential_std, rho_std, &
      !                         rho_gradient_std, kinetic_density_std)
      !hartree_potential => batch_irr%hartree_potential
      !rho => batch_irr%rho
      !rho_gradient => batch_irr%rho_gradient
      !kinetic_density => batch_irr%kinetic_density

      allocate(hartree_potential(batch_irr%n_full_points))
      call compute_batch_rlsy_full_to_irr_points(hartree_potential_std(:), hartree_potential(:))

      allocate(rho(n_spin,batch_irr%n_full_points))
      do i_spin = 1, n_spin
        call compute_batch_rlsy_full_to_irr_points(rho_std(i_spin,:),rho(i_spin,:))
      end do

      if (use_density_gradient) then
        allocate(rho_gradient(3,n_spin,batch_irr%n_full_points))
        do i_spin = 1, n_spin
          call compute_batch_rlsy_full_to_irr_points_gradient(rho_gradient_std(:,i_spin,:),rho_gradient(:,i_spin,:))
        end do
      else
        allocate(rho_gradient(3,n_spin,batch_irr%n_full_points))
      endif

      if (use_meta_gga) then
        allocate(kinetic_density(n_spin,batch_irr%n_full_points))
        do i_spin = 1, n_spin
          call compute_batch_rlsy_full_to_irr_points(kinetic_density_std(i_spin,:),kinetic_density(i_spin,:))
        end do
      else
        allocate(kinetic_density(n_spin,batch_irr%n_full_points))
      endif

      if (use_rlsy_symmetry .and. allocated(AS_dde_potential)) then
          allocate(AS_dde_potential_work(3,3,batch_irr%n_full_points))
          AS_dde_potential_work(:,:,1:batch_irr%n_full_points) = AS_dde_potential(:,:,1:batch_irr%n_full_points)
          !call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(1,1,:),AS_dde_potential_work(1,1,:))
          !call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(2,2,:),AS_dde_potential_work(2,2,:))
          !call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(3,3,:),AS_dde_potential_work(3,3,:))
          !call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(1,2,:),AS_dde_potential_work(1,2,:))
          !call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(1,3,:),AS_dde_potential_work(1,3,:))
          !call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(2,3,:),AS_dde_potential_work(2,3,:))
          !if (AS_components .eq. 9 ) then
          !  call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(2,1,:),AS_dde_potential_work(2,1,:))
          !  call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(3,1,:),AS_dde_potential_work(3,1,:))
          !  call compute_batch_rlsy_full_to_irr_points(AS_dde_potential(3,2,:),AS_dde_potential_work(3,2,:))
          !end if
      end if

    else
      n_my_batches_work = n_my_batches
      my_n_full_points = n_full_points
      n_basis_local = -1    ! Only used for load balancing
      ld_dens_mat = n_hamiltonian_matrix_size

      allocate(ins_idx(1))

      batches_work => batches
      partition_tab => partition_tab_std

      hartree_potential => hartree_potential_std
      rho => rho_std
      rho_gradient => rho_gradient_std
      kinetic_density => kinetic_density_std

      AS_dde_potential_work => AS_dde_potential
    endif

    allocate(batch_times(n_my_batches_work))
    !---------------------------------------------------------------------------

    ! Unfortunately if gga-xc is used we have to calculate gga forces in
    ! every iteration forces are calculated. This is because gga forces are
    ! summed up to pulay forces.

    ! By default, only the Pulay forces are calculated
    gga_forces_on_temp = .false.
    meta_gga_forces_on_temp = .false.
    Gnonmf_forces_on_temp = .false.

    ! When we're evaluating the term involving the non-energy-weighted density
    ! matrix, then we include the appropriate contribute to the force
    if (use_gga .and. (.not. dens_mat_is_en_weighted))then
      gga_forces_on_temp = .true.
      if (use_meta_gga) then
        meta_gga_forces_on_temp = .true.
      end if
    end if
    if (use_Gnonmf_forces.and.(.not.dens_mat_is_en_weighted)) then
        Gnonmf_forces_on_temp = .true.
    end if
    if (Gnonmf_forces_on_temp) then
      allocate( local_Gnonmf_derivs(n_max_batch_size),stat=info)
      call check_allocation(info, 'local_Gnonmf_derivs               ')
      allocate( Gnonmf_gradient_deriv(3,n_max_batch_size),stat=info)
      call check_allocation(info, 'Gnonmf_gradient_deriv             ')
    end if

    ! CC: For the stress computation, the Hessian is required even at the LDA
    !     level, but only in the non-relativistic case
    if (gga_forces_on_temp .or. &
         ( use_analytical_stress .and. (flag_rel .eq. REL_none) ).or.&
          Gnonmf_forces_on_temp ) then
      if (.not.allocated(cartesians)) then
        allocate(cartesians(n_max_cartesian, 0:l_wave_max,n_max_compute_atoms),&
             stat=info)
        call check_allocation(info, 'cartesians                    ')
      end if
      if (.not.allocated(sum_gradient)) then
        allocate(sum_gradient(3, (l_wave_max+1) ** 2, n_max_compute_atoms ),&
             stat=info)
        call check_allocation(info, 'sum_gradient                  ')
      end if
      if (.not.allocated(sum_hessian)) then
        allocate(sum_hessian(6, (l_wave_max+1) ** 2, n_max_compute_atoms ),&
             stat=info)
        call check_allocation(info, 'sum_hessian                   ')
      end if
      if (.not.allocated(hessian_basis_wave)) then
        call aims_allocate(hessian_basis_wave, n_max_compute_ham, 6, &
             n_max_batch_size, "hessian_basis_wave")
      end if
      if(.not. allocated(radial_wave_2nd_deriv))then
        allocate(radial_wave_2nd_deriv(n_max_compute_fns_ham),stat=info)
        call check_allocation(info, 'radial_wave_2nd_deriv         ')
      end if
    end if

    !  begin work

    t_zora = .false.

    if (dens_mat_is_en_weighted) then
      n_spin_local = 1
    else
      n_spin_local = n_spin
    end if

    l_ylm_max = l_wave_max
    call aims_allocate(gradient_basis_wave, n_max_compute_ham, 3,&
         n_max_batch_size, "gradient_basis_wave" )
    gradient_basis_wave = 0.
    call aims_allocate(kinetic_gradient_basis_wave, n_max_compute_ham, 3, &
         n_max_batch_size, "kinetic_gradient_basis_wave" )
    kinetic_gradient_basis_wave = 0.

    if(.not. gga_forces_on_temp.and..not.Gnonmf_forces_on_temp) then
      allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ), &
           stat=info )
      call check_allocation(info, 'dylm_dtheta_tab               ')
      dylm_dtheta_tab = 0.

      allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ), &
           stat=info )
      call check_allocation(info, 'scaled_dylm_dphi_tab          ')
      scaled_dylm_dphi_tab = 0.
    end if

    allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,stat=info)
    call check_allocation(info, 'ylm_tab                       ')
    ylm_tab = 0.

    allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) ,stat=info)
    call check_allocation(info, 'index_lm                      ')
    index_lm = 0

    call aims_allocate(forces_shell, n_max_compute_ham*n_max_compute_ham, &
         "forces_shell" )
    forces_shell = 0.
    call aims_allocate(H_times_psi, n_max_compute_ham, n_max_batch_size, &
         n_spin_local, "H_times_psi" )
    H_times_psi = 0.
    call aims_allocate(d_H_times_psi, n_max_compute_ham, n_max_batch_size, 3, &
         n_spin_local, "d_H_times_psi")
    d_H_times_psi = 0.

    allocate(radial_wave(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'radial_wave                   ')
    radial_wave = 0.

    allocate(radial_wave_deriv(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'radial_wave_deriv             ')
    radial_wave_deriv = 0.

    allocate(kinetic_wave(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'kinetic_wave                  ')
    kinetic_wave = 0.

    allocate(kinetic_wave_deriv(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'kinetic_wave_deriv            ')
    kinetic_wave_deriv = 0.

    call aims_allocate(wave, n_max_compute_ham, n_max_batch_size, &
         "wave")
    wave = 0.
    call aims_allocate(d_wave, n_max_compute_ham, n_max_batch_size, 3, &
         "d_wave" )
    d_wave = 0.

    allocate(i_basis(n_centers_basis_I),stat=info)
    call check_allocation(info, 'i_basis                       ')
    i_basis = 0

    if (flag_rel.eq.1) then
      ! allocate all arrays relevant for ZORA

      if (.not.allocated(dist_tab_full)) then
        allocate(dist_tab_full(n_centers_integrals),stat=info)
        call check_allocation(info, 'dist_tab_full                 ')
        dist_tab_full = 0.
      end if

      if (.not.allocated(dir_tab_full_norm)) then
        allocate(dir_tab_full_norm(3,n_centers_integrals),stat=info)
        call check_allocation(info, 'dir_tab_full_norm             ')
        dir_tab_full_norm = 0.
      end if

      if (.not.allocated(i_r_full)) then
        allocate(i_r_full(n_centers_integrals),stat=info)
        call check_allocation(info, 'i_r_full                      ')
        i_r_full = 0.
      end if

      if (.not.allocated(zora_vector1)) then
        call aims_allocate(zora_vector1, n_max_compute_ham, 3, &
             n_max_batch_size,n_spin_local, "zora_vector1" )
        zora_vector1 = 0.
      end if

      if (.not.allocated(zora_vector2)) then
        call aims_allocate(zora_vector2, n_max_compute_ham, 3, &
             n_max_batch_size,n_spin_local, "zora_vector2")
        zora_vector2 = 0.
      end if
    end if

    ! CC: Allocate arrays required for stress
    if (use_analytical_stress) then
      if (.not.allocated(AS_strain_deriv_wave)) &
           allocate(AS_strain_deriv_wave(1:n_max_compute_ham, &
                    1:n_max_batch_size,1:AS_components))
           AS_strain_deriv_wave = 0.
      if (.not.allocated(AS_strain_deriv_wave_shell))&
           allocate(AS_strain_deriv_wave_shell (1:n_max_compute_ham, &
                    1:n_max_compute_ham))
           AS_strain_deriv_wave_shell = 0.
      if ( use_AS_Jac_in_pulay ) then
        if (.not.allocated(AS_jac_pot_kin_times_psi)) &
           allocate(AS_jac_pot_kin_times_psi(1:n_max_compute_ham, &
                    1:n_max_batch_size,n_spin_local,1:AS_components))
           AS_jac_pot_kin_times_psi = 0.
      end if
      if (flag_rel .eq. REL_none) then
        if (.not.allocated(AS_Hessian_times_psi)) &
            allocate(AS_Hessian_times_psi(1:n_max_compute_ham, &
                     1:n_max_batch_size,1:AS_components))
            AS_Hessian_times_psi = 0.
        if (.not.allocated(AS_T_times_psi)) &
            allocate(AS_T_times_psi(1:n_max_compute_ham,1:n_max_batch_size))
            AS_T_times_psi = 0.
        if (.not.allocated(AS_corr_shell)) &
            allocate(AS_corr_shell(1:n_max_compute_ham,1:n_max_compute_ham))
            AS_corr_shell = 0.
        if (.not.allocated(AS_i_basis_on_site)) &
            allocate(AS_i_basis_on_site(1:(n_max_compute_ham)**2,2))
            ! FIXME: Find better estimate for size
            AS_i_basis_on_site = 0
      else
        if (.not.allocated(AS_strain_deriv_kinetic_wave)) &
            allocate(AS_strain_deriv_kinetic_wave (1:n_max_compute_ham, &
                     1:n_max_batch_size,1:AS_components))
            AS_strain_deriv_kinetic_wave = 0.
      end if
      if (gga_forces_on_temp) then
        if (.not.allocated(AS_hessian_times_xc_deriv_gga)) &
            allocate(AS_hessian_times_xc_deriv_gga(1:n_max_compute_ham, &
                     1:n_max_batch_size,n_spin_local,1:AS_components))
            AS_hessian_times_xc_deriv_gga = 0.

        ! For AS_eval_hessian_wave_times_xc_gradient_deriv_gga, we need a
        ! mapping AS_strain_index -> needed hessian_components
        AS_index_hessian_gga(1,1:3) = (/ 1,2,3 /)
        AS_index_hessian_gga(2,1:3) = (/ 1,2,3 /)
        AS_index_hessian_gga(3,1:3) = (/ 1,2,3 /)
        AS_index_hessian_gga(4,1:3) = (/ 2,4,5 /)
        AS_index_hessian_gga(5,1:3) = (/ 2,4,5 /)
        AS_index_hessian_gga(6,1:3) = (/ 3,5,6 /)
        AS_index_hessian_gga(7,1:3) = (/ 2,4,5 /)
        AS_index_hessian_gga(8,1:3) = (/ 3,5,6 /)
        AS_index_hessian_gga(9,1:3) = (/ 3,5,6 /)

        ! Need to declare this irresepctive of whether we are using mggas,
        ! Could change it to a dummy declaration for the GGA case?
        ! (meta_gga_forces_on_temp) then
        if (.not.allocated(AS_hessian_times_xc_deriv_mgga)) &
             allocate(AS_hessian_times_xc_deriv_mgga(1:n_max_compute_ham, &
                    1:n_max_batch_size,n_spin_local,1:AS_components))
      end if
    end if

    ! Initialization when using GPU acceleration
    if (use_gpu) then
      index_on_gpu = n_bp.gt.0

      ! This is better than the copy-paste used in the Hamiltonian integration
      ! and density update, but there must be a better way...
      ! Possibly *_create_gpu makes two passes:  one where it passes back two
      ! arrays back to the calling subroutine containing a list of all array
      ! which will be allocated.  Then we just loop over all elements of the
      ! array.  Then in the second pass, we actually allocate the array.
      write (info_str,'(A)') &
         'Reporting GPU memory usage by myid 0.  Note that it is likely not &
         &the only MPI task binding to its GPU!'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )

      call output_mem_array_gpu("partition", &
           n_max_batch_size, 8)
      call output_mem_array_gpu("dev_waveComputeA", &
           n_max_compute_ham * n_max_batch_size, 8)
      call output_mem_array_gpu("forces_shell", &
           n_max_compute_ham * n_max_compute_ham, 8)

      call output_mem_array_gpu("h_times_psi", &
           n_max_compute_ham * n_max_batch_size * n_spin, 8)
      call output_mem_array_gpu("d_wave", &
           n_max_compute_ham * n_max_batch_size * 3, 8)

      if (gga_forces_on_temp .or. flag_rel.eq.REL_atomic_zora .or. &
           use_AS_Jac_in_pulay) then
        call output_mem_array_gpu("wave", &
             n_max_compute_ham * n_max_batch_size, 8);
      end if

      if (gga_forces_on_temp) then
        call output_mem_array_gpu("xc_gradient_deriv", &
             3 * n_spin * n_max_batch_size, 8);
        call output_mem_array_gpu("hessian_basis_wave", &
             6 * n_max_compute_ham * n_max_batch_size, 8);
        call output_mem_array_gpu("gradient_basis_wave", &
             n_max_compute_ham * 3 * n_max_batch_size, 8);
        call output_mem_array_gpu("dev_matrixTmp1", &
             n_max_compute_ham * n_max_batch_size, 8);
        call output_mem_array_gpu("dev_matrixTmp2", &
             n_max_compute_ham * n_max_batch_size, 8);
      end if

      if (meta_gga_forces_on_temp) then
        call output_mem_array_gpu("xc_tau_deriv", &
             n_spin * n_max_batch_size, 8);
        if (use_analytical_stress) then
          call output_mem_array_gpu("dev_matrixTmp1MGGA", &
               n_max_compute_ham * 3 *n_max_batch_size, 8);
          call output_mem_array_gpu("dev_matrixTmp2MGGA", &
               n_max_compute_ham * 3 * n_max_batch_size, 8);
        else
          call output_mem_array_gpu("dev_matrixTmp1MGGA", &
               n_max_compute_ham * n_max_batch_size, 8);
          call output_mem_array_gpu("dev_matrixTmp2MGGA", &
               n_max_compute_ham * n_max_batch_size, 8);
        end if
      end if

      if (use_analytical_stress) then
        call output_mem_array_gpu("as_strain_deriv_wave", &
             n_max_compute_ham * n_max_batch_size * as_components, 8);
        call output_mem_array_gpu("as_strain_deriv_wave_shell", &
             n_max_compute_ham * n_max_compute_ham, 8);
        if (gga_forces_on_temp) then
          call output_mem_array_gpu("dev_matrixTmp3", &
              3 * n_max_batch_size, 8);
          call output_mem_array_gpu("as_hessian_times_xc_deriv_gga", &
              n_max_compute_ham * n_max_batch_size * n_spin * as_components, 8)
        end if
        if (meta_gga_forces_on_temp) then
          call output_mem_array_gpu("as_hessian_times_xc_deriv_mgga", &
              n_max_compute_ham * n_max_batch_size * n_spin * as_components, 8)
        end if
        if (use_AS_Jac_in_pulay) then
          call output_mem_array_gpu("as_jac_pot_kin_times_psi", &
               n_max_compute_ham * n_max_batch_size * n_spin * as_components, 8)
        end if
        if (flag_rel.eq.REL_atomic_zora) then
          call output_mem_array_gpu("as_strain_deriv_kinetic_wave", &
               n_max_compute_ham * n_max_batch_size * as_components, 8)
          call output_mem_array_gpu("as_strain_deriv_wave_shell_trans", &
               n_max_compute_ham * n_max_compute_ham, 8)
        end if
      end if

      if (flag_rel.eq.REL_atomic_zora) then
        call output_mem_array_gpu("h_times_psi", &
               n_max_compute_ham * n_max_batch_size * n_spin * 3, 8)
      end if

      if (index_on_gpu) then
        call output_mem_array_gpu("ins_idx", &
             n_basis_local, 4)
        call output_mem_array_gpu("dens_mat", &
             ld_dens_mat * n_spin, 8)
        call output_mem_array_gpu("permute_compute_by_atom", &
             n_max_compute_ham, 4)
        call output_mem_array_gpu("dev_forceValues", &
             n_max_compute_ham * n_max_compute_ham, 8)
        call output_mem_array_gpu("sum_forces", &
               3 * n_atoms, 8)
        if (use_analytical_stress) then
          call output_mem_array_gpu("as_pulay_stress_local", &
               as_components, 8);
        end if
      end if

      call forces_create_gpu( &
           n_max_compute_ham, &
           n_max_batch_size, &
           ld_dens_mat, &
           n_spin_local, &
           n_atoms, &
           n_basis, &
           n_basis_local, &
           AS_components, &
           fort_log_to_c_int(gga_forces_on_temp), &
           fort_log_to_c_int(meta_gga_forces_on_temp), &
           fort_log_to_c_int(use_analytical_stress), &
           fort_log_to_c_int(flag_rel.eq.REL_atomic_zora), &
           fort_log_to_c_int(use_AS_Jac_in_pulay), &
           fort_log_to_c_int(index_on_gpu))

      if (index_on_gpu) then
        ! Density matrix is used when indexing the force batches back into the
        ! forces, but we only do indexing on GPU for non-packed matrices or
        ! load balanced matrices
        call set_dens_mat_gpu(dens_mat, ld_dens_mat)

        ! Arrays for permuting compute basis elements so that compute basis
        ! elements are ordered sequentially by atom.  This ordering of the
        ! matrix elements allows for naive reduction on blocks of memory to
        ! generate forces on each atom.
        allocate(permute_compute_by_atom(n_max_compute_ham))
        allocate(n_compute_for_atom(n_atoms))

        ! The following set_* subroutines are needed, because we will be calling
        ! this subroutine multiple times.  If we do not reset the values on
        ! the GPU, the values from previous cycles will be wiped out, just like
        ! the day of my life spent finding this out was.
        call set_sum_forces_gpu(sum_forces, n_atoms)
        if (use_analytical_stress) then
          call set_as_pulay_stress_local_gpu(AS_pulay_stress_local)
        end if
      end if
    end if

    ! Compute basis permutation is only used when performing GPU acceleration
    ! with load balancing, allocate dummy arrays otherwise
    if (.not.allocated(permute_compute_by_atom)) then
      allocate(permute_compute_by_atom(1))
      allocate(n_compute_for_atom(1))
    end if

    ! initialize
    i_index = 0
    do i_l = 0, l_wave_max, 1
      do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
      end do
    end do

    i_full_points_C = 0
    i_full_points_2C = 0
    i_full_points_A = 0

    ! perform partitioned integration, atom by atom, and point by point
    ! This will be the outermost loop, to save evaluations of the potential.
    ! and the Y_lm functions
    do i_my_batch = 1, n_my_batches_work, 1
      n_compute_c = 0
      n_compute_a = 0
      i_basis = 0
      i_point = 0

      ! loop within the current batch
      do i_index = 1, batches_work(i_my_batch)%size, 1
        i_full_points_2C = i_full_points_2C + 1

        if (partition_tab(i_full_points_2C).gt.0.d0) then
          i_point = i_point+1

          ! get current integration point coordinate
          coord_current(:) = batches_work(i_my_batch)%points(i_index)%coords(:)

          if (n_periodic > 0) then
            call map_to_center_cell(coord_current(1:3) )
          end if

          ! compute atom-centered coordinates of current integration point,
          ! as viewed from all atoms
          call tab_atom_centered_coords_p0 &
               ( coord_current,  &
               dist_tab_sq(1,i_point),  &
               dir_tab(1,1,i_point), &
               n_centers_integrals, centers_basis_integrals )

          ! determine which basis functions are relevant at current integration
          ! points, and tabulate their indices

          ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are
          ! actually needed
          if (.not.prune_basis_once) then
            call prune_basis_p0 &
                 ( dist_tab_sq(1,i_point), &
                 n_compute_a, n_compute_c, i_basis,  &
                 n_centers_basis_I, n_centers_integrals, &
                 inv_centers_basis_integrals  )
          end if
        end if
      end do  ! end loop within one batch

      if (prune_basis_once) then
        n_compute_a = batches_work(i_my_batch)%batch_n_compute
        n_compute_c = n_compute_a
        i_basis(1:n_compute_a) = batches_work(i_my_batch)%batch_i_basis
      end if

      !! CC: Here, we construct an array of pairs AS_i_basis_on_site(:,1:2) that
      !!     includes all on-site terms. We can thereby transform
      !!          dgemm -> sum_AS_n_on_site ddot( AS_i_basis_on_site_1 , AS_i_basis_on_site_2 )
      !!     later on.
      if ( use_analytical_stress .and. (flag_rel.eq.REL_none) ) then
        call AS_compute_i_on_site_basis( n_compute_a, i_basis(1), &
             AS_n_on_site, AS_i_basis_on_site(1,1), &
             AS_i_basis_on_site(1,2) )
      end if

      n_points = i_point

      ! Perform actual integration if more than 0 basis functions
      ! are actually relevant on the present angular shell ...
      if (n_compute_c.gt.0) then
        n_rel_points = 0
        i_point = 0

        ! loop over one division of the angular grid
        do i_index = 1, batches_work(i_my_batch)%size, 1
          ! Increment the (global) counter for the grid, to access storage
          ! arrays
          i_full_points_C = i_full_points_C + 1
          i_full_points_A = i_full_points_A + 1

          if (partition_tab(i_full_points_C).gt.0.d0) then
            i_point = i_point+1

            if (flag_rel.eq.1) then
              call tab_global_geometry_p0 &
                   ( dist_tab_sq(1,i_point), &
                   dir_tab(1,1,i_point), &
                   dist_tab_full, &
                   i_r_full, &
                   dir_tab_full_norm, &
                   n_centers_integrals,  centers_basis_integrals)
            end if

            ! for all integrations
            partition(i_point) = partition_tab(i_full_points_C)
            energy_partition(i_point) = partition_tab(i_full_points_A)

            ! CC: Save strain derivative of the potential, as calculated in
            !     sum_up_...
            if ( use_analytical_stress .and. use_AS_Jac_in_pulay ) then
               ! The AS_dde_potential array is only allocated in sum_up_whole_potential_p1,
               ! i.e., potentially after the first call to pulay_forces_densmat.
               ! Only in the second s.c.f. iteration after forces are on is the
               ! AS_dde_potential actually available.
               !
               ! In principle, this term should never be reached before AS_counter has reached
               ! a value of 2 or greater. In scf_solver.f90 (and if analytical stress is used),
               ! pulay_forces_densmat is only called if AS_counter .gt. 1.
               !
               ! As an additional safeguard, we use the allocation status of AS_dde_potential
               ! (yes or no) to determine if this array is actually available.
               if (allocated(AS_dde_potential)) then
                  AS_dde_potential_local(1) = AS_dde_potential_work(1,1,i_full_points_C)
                  AS_dde_potential_local(4) = AS_dde_potential_work(2,2,i_full_points_C)
                  AS_dde_potential_local(6) = AS_dde_potential_work(3,3,i_full_points_C)
                  AS_dde_potential_local(2) = AS_dde_potential_work(1,2,i_full_points_C)
                  AS_dde_potential_local(3) = AS_dde_potential_work(1,3,i_full_points_C)
                  AS_dde_potential_local(5) = AS_dde_potential_work(2,3,i_full_points_C)
                  if (AS_components .eq. 9 ) then
                     AS_dde_potential_local(7) = AS_dde_potential_work(2,1,i_full_points_C)
                     AS_dde_potential_local(8) = AS_dde_potential_work(3,1,i_full_points_C)
                     AS_dde_potential_local(9) = AS_dde_potential_work(3,2,i_full_points_C)
                  end if
               else ! AS_dde_potential is not yet available
                  AS_dde_potential_local(:) = 0.d0  ! Omit this term while not yet available
               end if
            end if

            n_compute_atoms = 0
            n_compute_fns = 0
            i_basis_fns_inv = 0
            atom_index_inv = 0

            ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if
            ! needed) are stored in a compact spline array that can be accessed
            ! by spline_vector_waves, without any copying and without doing any
            ! unnecessary operations.  The price is that the interface is no
            ! longer explicit in terms of physical objects. See
            ! shrink_fixed_basis() for details regarding the reorganized spline
            ! arrays.
            call prune_radial_basis_p0 &
                 ( dist_tab_sq(1,i_point), &
                 dist_tab(1,i_point), &
                 dir_tab(1,1,i_point), &
                 n_compute_atoms, atom_index, atom_index_inv, &
                 n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                 i_atom_fns, spline_array_start, spline_array_end, &
                 n_centers_integrals, centers_basis_integrals)

            ! Tabulate distances, unit vectors, and inverse logarithmic grid
            ! units for all atoms which are actually relevant
            call tab_local_geometry_p0 &
                 ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                 dir_tab(1,1,i_point), dist_tab(1,i_point), i_r )


            ! compute trigonometric functions of spherical coordinate angles
            ! of current integration point, viewed from all atoms
            call tab_trigonom_p0 &
                 ( n_compute_atoms, dir_tab(1,1,i_point), trigonom_tab )

            ! tabulate distance and Ylm's w.r.t. other atoms
            call tab_wave_ylm_p0 &
                 ( n_compute_atoms, atom_index,  &
                 trigonom_tab, basis_l_max,  &
                 l_ylm_max, ylm_tab )

            ! Now evaluate radial functions
            ! from the previously stored compressed spline arrays
            call evaluate_radial_functions_p0  &
                 (   spline_array_start, spline_array_end,  &
                 n_compute_atoms, n_compute_fns,   &
                 dist_tab(1,i_point), i_r,  &
                 atom_index, i_basis_fns_inv,  &
                 basis_wave_ordered, radial_wave,  &
                 .false. , n_compute_c, n_max_compute_fns_ham )

            ! tabulate total wave function value for each basis function
            call evaluate_waves_p0  &
                 ( l_ylm_max,   &
                 ylm_tab, dist_tab(1,i_point),   &
                 index_lm, n_compute_c,   &
                 i_basis, radial_wave(1),   &
                 wave(1,i_point), n_compute_atoms,   &
                 atom_index_inv, n_compute_fns,  &
                 i_basis_fns_inv,  n_max_compute_fns_ham )

            ! in the remaining part of the subroutine, some decisions (scalar
            !  relativity) depend on the potential; must therefore evaluate the
            ! potential and derived quantities right here

            ! Local exchange-correlation parts of the potential are evaluated
            ! right here, to avoid having to store them separately elsewhere.
            ! For large systems, savings are significant

            ! This one is required for vdW-DF
            coord_current(:) = &
                 batches_work(i_my_batch)%points(i_index)%coords(:)!SAG
            if (n_periodic > 0) then
              call map_to_center_cell(coord_current(1:3) )
            end if

            call evaluate_xc  &
                 ( rho(1,i_full_points_A),   &
                 rho_gradient(1,1,i_full_points_A),  &
                 kinetic_density(1,i_full_points_A), &
                 en_density_xc(i_point), &
                 en_density_x, en_density_c,  &
                 local_xc_derivs(1,i_point),  &
                 xc_gradient_deriv(1,1,i_point), &
                 xc_tau_deriv(1,i_point), &
                 .false., &
                 coord_current(:) &
                 )
            if (Gnonmf_forces_on_temp) then
              call evaluate_Gnonmf(rho(:,i_full_points_A), &
                   rho_gradient(:,:,i_full_points_A), &
                   local_Gnonmf_derivs(i_point), &
                   Gnonmf_gradient_deriv(:,i_point), &
                   i_full_points_A)
            end if

            do i_spin = 1, n_spin_local, 1
              local_potential_parts(i_spin) =   &
                   hartree_potential(i_full_points_A) +   &
                   local_xc_derivs(i_spin,i_point)
              if (Gnonmf_forces_on_temp) then
                local_potential_parts(i_spin) = &
                     local_potential_parts(i_spin) + &
                     local_Gnonmf_derivs(i_point)
              end if
              !SR: what is the purpose of the variable
              !sum_of_local_gradients here? do my knowledge it is not
              !used in any sense??
              !CC: sum_of_local_gradients is used in add_zora_gradient_part_p0
              !    later on
              if (use_gga) then
                sum_of_local_gradients(1:3,i_spin) =   &
                     xc_gradient_deriv(1:3,i_spin,i_point)*4.d0
              else
                sum_of_local_gradients = 0.d0
              end if
              if (Gnonmf_forces_on_temp) then
                sum_of_local_gradients(1:3,i_spin) =   &
                     sum_of_local_gradients(1:3,i_spin)+&
                     Gnonmf_gradient_deriv(1:3,i_point)*2.d0
              end if

              ! CC: Save xc en_dens - pot
              if ( (use_analytical_stress).and.(use_AS_Jac_in_pulay) ) then
                AS_en_xc_minus_pot_xc(i_spin) = &
                     en_density_xc(i_point) - local_xc_derivs(i_spin,i_point)
              end if
            enddo

            ! Check whether relativistic corrections are needed at the present
            ! point.  The check is based entirely on the local parts of the
            ! potential - i.e.  in a GGA, the terms due to
            ! d(rho*exc)/d(|grad(rho|^2) is not evaluated.  Hopefully this
            ! approximation to the full ZORA energy is small.
            if (flag_rel.eq.1) then
              ! if we need ZORA, must get the _full_ local geometry in order to
              ! create the superposition of atomic potentials which is used to
              ! estimate the potential gradient for ZORA
              call evaluate_pot_superpos_p0  &
                   ( i_r_full, zora_potential_parts(1),  &
                   n_centers_integrals, centers_basis_integrals )
              do i_spin = 1, n_spin_local, 1
                zora_operator(i_spin) =  &
                     light_speed_sq /  &
                     ( 2 * light_speed_sq -  &
                     zora_potential_parts(i_spin) )
              end do
            end if

            ! tabulate radial derivatives of those radial functions
            ! which are actually non-zero at current point, using vectorized
            ! splines
            call evaluate_radial_functions_p0  &
                 ( spline_array_start, spline_array_end,  &
                 n_compute_atoms, n_compute_fns,   &
                 dist_tab(1,i_point), i_r,  &
                 atom_index, i_basis_fns_inv,  &
                 basis_deriv_ordered,   &
                 radial_wave_deriv(1), .true.,  &
                 n_compute_c, n_max_compute_fns_ham )

            ! CC: For the stress computation, the Hessian is required even at
            !     the LDA level, but only in the non-relativistic case
            if (dens_mat_is_en_weighted .or. &
                 .not.(gga_forces_on_temp .or.Gnonmf_forces_on_temp.or.&
                 (use_analytical_stress.and.(flag_rel.eq.REL_none))))&
                 then ! SR: only continue if GGA-force.and.Gnonmf-force&
                 !.and.analytical stress part are all false

              ! of current integration point, viewed from all atoms
              call tab_trigonom_p0 &
                   ( n_compute_atoms, dir_tab(1,1,i_point),  &
                   trigonom_tab(1,1) &
                   )

              ! tabulate distance and Ylm's w.r.t. other atoms
              call tab_wave_ylm_p0 &
                   ( n_compute_atoms, atom_index,  &
                   trigonom_tab(1,1), basis_l_max,  &
                   l_ylm_max, &
                   ylm_tab(1,1) )

              ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
              call tab_gradient_ylm_p0  &
                   ( trigonom_tab(1,1), basis_l_max,   &
                   l_ylm_max, n_compute_atoms, atom_index,  &
                   ylm_tab(1,1),   &
                   dylm_dtheta_tab(1,1),   &
                   scaled_dylm_dphi_tab(1,1)  )

              ! and finally, assemble the actual gradients
              call evaluate_wave_gradient_p0  &
                   ( dist_tab(1,i_point),  &
                   dir_tab(1,1,i_point),  &
                   trigonom_tab(1,1),  &
                   l_ylm_max, ylm_tab(1,1),  &
                   dylm_dtheta_tab(1,1),  &
                   scaled_dylm_dphi_tab(1,1),  &
                   index_lm, n_compute_c,  &
                   i_basis(1:n_compute_c),  &
                   radial_wave(1),  &
                   radial_wave_deriv(1),  &
                   gradient_basis_wave (1,1,i_point),  &
                   n_compute_atoms,  &
                   atom_index_inv,  &
                   n_compute_fns,  &
                   i_basis_fns_inv, n_max_compute_fns_ham   )
            else
              ! GGA-force part or Gnonmf_force part or analytical stress is
              ! calculated
              ! Derivative of the radial part of the wave functions.
              call evaluate_radial_functions_deriv_p0 &
                   (spline_array_start, spline_array_end, &
                   n_compute_atoms, n_compute_fns,  &
                   dist_tab(1,i_point), i_r(1), &
                   atom_index, i_basis_fns_inv, &
                   basis_deriv_ordered,  &
                   radial_wave_2nd_deriv(1),  &
                   .true., n_max_compute_fns_dens &
                   )

              ! Basis functions, their gradients and hessians are needed.
              ! They are evaluated on the basis of a cartesian expansion of the
              ! ylm-functions - hence, if we do need all this, we get a lot of
              ! stuff directly that we'd have to compute separately otherwise
              ! (see below)

              ! first, corresponding cartesian terms x^l_x*y^l_y*z^l_z are
              ! evaluated for the current integration point
              call evaluate_cartesians &
                   (dir_tab(1,1,i_point), basis_l_max, l_wave_max,  &
                   atom_index, n_compute_atoms, cartesians(1,0,1) )

              ! further "ingredients" of the hessian are evaluated based on the
              ! cartesians ylm's are evaluated on the fly
              call evaluate_cartesian_hessian_and_gradient_terms &
                   (basis_l_max, l_wave_max, cartesians(1,0,1), &
                   atom_index, n_compute_atoms,  &
                   ylm_tab(1,1),  &
                   sum_gradient(1,1,1),  &
                   sum_hessian(1,1,1))

              ! now, all hessians of the basis functions are evaluated based on
              ! the "ingredients" (ylm, sum_two_gradient, sum_hessian)
              ! (gradient and ylm-functions themselves would be for free, we
              ! have to check the performance)
              call evaluate_wave_gradient_cartesian_p1 &
                   (dist_tab(1,i_point), i_r(1),  &
                   dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                   n_compute_c, i_basis, atom_index_inv,  &
                   i_basis_fns_inv, radial_wave(1),  &
                   radial_wave_deriv(1), &
                   ylm_tab(1,1),  &
                   sum_gradient(1,1,1),  &
                   gradient_basis_wave(1,1,i_point),n_compute_atoms)

              ! ugh ... the hessian itself.
              call evaluate_wave_hessian_cartesian_p0 &
                   (dist_tab(1,i_point), i_r(1),  &
                   dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                   n_compute_c, i_basis, atom_index_inv,  &
                   i_basis_fns_inv, radial_wave(1),  &
                   radial_wave_deriv(1),  &
                   radial_wave_2nd_deriv(1),  &
                   ylm_tab(1,1),  &
                   sum_gradient(1,1,1),  &
                   sum_hessian(1,1,1), &
                   hessian_basis_wave(1,1,i_point), n_compute_atoms)
            end if ! gga forces/Gnonmf_forces/analytical stress end-------------

            ! CC: Whe have the gradients of the orbitals by now, so we can
            !     compute the strain derivative of orbital
            if (use_analytical_stress) then
              do AS_hessian_index=1,AS_components,1
                call AS_eval_strain_deriv_wave(          &
                    n_compute_c, n_compute_atoms,        &
                    dist_tab(1,i_point),                 &
                    dir_tab(1,1,i_point),                &
                    i_basis(1),                          &
                    i_basis_fns_inv(1,1),                &
                    atom_index_inv(1),                   &
                    gradient_basis_wave(1,1,i_point),    &
                    AS_hessian_index,                    &
                    AS_strain_deriv_wave(1,i_point,AS_hessian_index) )

                if (gga_forces_on_temp) then
                  ! GGA: Evaluate xc_gradient_deriv * hessian_basis_wave * distance.
                  do i_spin = 1, n_spin_local, 1
                    call AS_eval_hessian_wave_times_xc_gradient_deriv_gga(                       &
                         n_compute_c,                                                            &
                         n_compute_atoms,                                                        &
                         dist_tab(1,i_point),                                                    &
                         dir_tab(1,1,i_point),                                                   &
                         i_basis(1),                                                             &
                         i_basis_fns_inv(1,1),                                                   &
                         atom_index_inv(1),                                                      &
                         hessian_basis_wave(1,AS_index_hessian_gga(AS_hessian_index,1),i_point), &
                         hessian_basis_wave(1,AS_index_hessian_gga(AS_hessian_index,2),i_point), &
                         hessian_basis_wave(1,AS_index_hessian_gga(AS_hessian_index,3),i_point), &
                         xc_gradient_deriv(1,i_spin,i_point),                                    &
                         AS_hessian_index,                                                       &
                         AS_hessian_times_xc_deriv_gga(1,i_point,i_spin,AS_hessian_index)        )

                    ! Calculation of Meta-GGA hessian term: xc_tau_deriv * hessian_basis_wave * distance
                    if (meta_gga_forces_on_temp) then
                        call AS_eval_hessian_wave_times_xc_tau_deriv_mgga(                      &
                        n_compute_c,                                                            &
                        n_compute_atoms,                                                        &
                        dist_tab(1,i_point),                                                    &
                        dir_tab(1,1,i_point),                                                   &
                        i_basis(1),                                                             &
                        i_basis_fns_inv(1,1),                                                   &
                        atom_index_inv(1),                                                      &
                        hessian_basis_wave(1,AS_index_hessian_gga(AS_hessian_index,1),i_point), &
                        hessian_basis_wave(1,AS_index_hessian_gga(AS_hessian_index,2),i_point), &
                        hessian_basis_wave(1,AS_index_hessian_gga(AS_hessian_index,3),i_point), &
                        xc_tau_deriv(i_spin,i_point),                                           &
                        AS_hessian_index,                                                       &
                        AS_hessian_times_xc_deriv_mgga(1,i_point,i_spin,AS_hessian_index)        )
                    end if
                  end do ! i_spin
                end if ! gga
              end do ! AS_hessian_index
            end if ! analytic_stress

            ! d_wave(1:n_compute_c,i_point) = gradient_basis_wave(1:n_compute_c,i_calculate_dimension)
            do i_calculate_dimension = 1,3
              call copy_gradient(d_wave(1,i_point,i_calculate_dimension),&
                   gradient_basis_wave(1,1, i_point), n_compute_c, &
                   i_calculate_dimension)
            end do

            ! paula
            if (dens_mat_is_en_weighted) then
              ! In this case, we only need the orbitals and their derivatives,
              ! which we already have

              ! H_times_psi(:, :, 1) =  wave(:,:)
              H_times_psi(1:n_compute_c, i_point, 1) = &
                   wave(1:n_compute_c,i_point)
            else
              ! In this case, we evaluate vector of components H*phi(i,r) and
              ! possibly its derivative as well

              ! Local potential parts first; in the case of GGA,
              ! the real gradient parts are added further below
              !               if ( (flag_rel/=1)) then
              ! Non-relativistic treatment - simply evaluate
              ! H*phi(i,r) all in one

              ! First, obtain radial kinetic energy terms from vectorized
              ! splines
              call evaluate_radial_functions_p0  &
                   ( spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_kinetic_ordered, kinetic_wave(1),  &
                   .false., n_compute_c, n_max_compute_fns_ham )

              if (flag_rel == REL_atomic_zora) then
                kinetic_wave = kinetic_wave *0.5
              end if

              do i_spin = 1, n_spin_local, 1
                call evaluate_H_psi_p0  &
                     ( l_ylm_max,  &
                     ylm_tab(1:(l_ylm_max+1)**2,1:n_max_compute_atoms),  &
                     dist_tab(1, i_point),  &
                     index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),  &
                     H_times_psi(1, i_point, i_spin),  &
                     radial_wave(1),  &
                     local_potential_parts(i_spin),  &
                     n_compute_c,  &
                     i_basis(1:n_max_compute_ham),  &
                     n_compute_atoms, atom_index_inv,  &
                     n_compute_fns,   &
                     i_basis_fns_inv,   &
                     kinetic_wave(1),  &
                     zora_operator(i_spin), n_max_compute_fns_ham )
              end do

              ! CC: Construct all operators for the evaluation of the
              !     Jacobian, Potential, Kinetic, and Correction terms
              if (use_analytical_stress) then
                ! Non-relativistic case
                if (flag_rel.eq.REL_none) then
                  call AS_construct_T_times_psi(             &
                       n_compute_c,                               &
                       H_times_psi(1, i_point, 1),                &
                       local_potential_parts(1),                  &
                       wave(1,i_point),                           &
                       AS_T_times_psi(1,i_point) )

                  do AS_hessian_index=1,AS_components,1
                    do i_spin = 1, n_spin_local, 1
                      call AS_eval_jac_pot_kin_times_psi(                                     &
                           AS_hessian_index, n_compute_c,                                      &
                           AS_en_xc_minus_pot_xc(i_spin),                                      &
                           wave(1, i_point),                                                   &
                           AS_dde_potential_local(AS_hessian_index),                           &
                           AS_jac_pot_kin_times_psi(1,i_point,i_spin,AS_hessian_index),        &
                           (flag_rel .eq. REL_atomic_zora               ),                     &
                           hessian_basis_wave(1,AS_map_hess_to_symm(AS_hessian_index),i_point) )

                      ! Add GGA part to Jacobian
                      if ( gga_forces_on_temp .and. &
                           ((AS_hessian_index.eq.1) .or. (AS_hessian_index.eq.4) .or. (AS_hessian_index.eq.6)) ) then
                        call AS_add_gga_jac(                                             &
                             AS_jac_pot_kin_times_psi(1,i_point,i_spin,AS_hessian_index), &
                             n_compute_c,                                                 &
                             xc_gradient_deriv(1,i_spin,i_point),                         &
                             gradient_basis_wave(1,1,i_point)                             )
                      end if
                    end do

                    call AS_construct_hessian_op_times_psi(             &
                         AS_hessian_index, n_compute_c,                      &
                         hessian_basis_wave(1,AS_map_hess_to_symm(AS_hessian_index),i_point),     &
                         AS_T_times_psi(1,i_point),                          &
                         AS_Hessian_times_psi(1,i_point,AS_hessian_index) )
                  end do

                else if ( use_AS_Jac_in_pulay ) then ! Relativistic (atomic zora), JAC if requested
                  do AS_hessian_index=1,AS_components,1
                    do i_spin = 1, n_spin_local, 1
                      call AS_eval_jac_pot_kin_times_psi(                             &
                           AS_hessian_index, n_compute_c,                              &
                           AS_en_xc_minus_pot_xc(i_spin),                              &
                           wave(1, i_point),                                           &
                           AS_dde_potential_local(AS_hessian_index),                   &
                           AS_jac_pot_kin_times_psi(1,i_point,i_spin,AS_hessian_index),&
                           (flag_rel .eq. REL_atomic_zora               )              )

                      ! Add GGA part to Jacobian
                      if ( gga_forces_on_temp .and. &
                           ((AS_hessian_index.eq.1) .or. (AS_hessian_index.eq.4) .or. (AS_hessian_index.eq.6)) ) then
                        call AS_add_gga_jac(                                             &
                             AS_jac_pot_kin_times_psi(1,i_point,i_spin,AS_hessian_index), &
                             n_compute_c,                                                 &
                             xc_gradient_deriv(1,i_spin,i_point),                         &
                             gradient_basis_wave(1,1,i_point)                             )
                      end if
                    end do ! i_spin
                  end do ! AS_hessian_index
                end if ! flag_rel
              end if ! analytic_stress

              if ((flag_rel.eq.REL_atomic_zora)) then
                !paula----------------------------------------------------------
                call evaluate_radial_functions_deriv_p0 &
                     (spline_array_start, spline_array_end, &
                     n_compute_atoms, n_compute_fns,  &
                     dist_tab(1,i_point), i_r(1), &
                     atom_index, i_basis_fns_inv, &
                     basis_kinetic_ordered,  &
                     kinetic_wave_deriv(1),  &
                     .true., n_max_compute_fns_dens &
                     )
                kinetic_wave_deriv =  0.5d0*kinetic_wave_deriv

                if(.not.(gga_forces_on_temp .or.Gnonmf_forces_on_temp.or.&
                     (use_analytical_stress.and.(flag_rel.eq.REL_none)))) &
                     then ! GGA-force part is not calculated
                  ! and Gnonmf-force also not calculated
                  ! and finally, assemble the actual gradients
                  call evaluate_wave_gradient_p0  &
                       ( dist_tab(1,i_point),  &
                       dir_tab(1,1,i_point),  &
                       trigonom_tab(1,1),  &
                       l_ylm_max, ylm_tab(1,1),  &
                       dylm_dtheta_tab(1,1),  &
                       scaled_dylm_dphi_tab(1,1),  &
                       index_lm, n_compute_c,  &
                       i_basis(1:n_compute_c),  &
                       kinetic_wave(1),  &
                       kinetic_wave_deriv(1),  &
                       kinetic_gradient_basis_wave (1,1,i_point),  &
                       n_compute_atoms,  &
                       atom_index_inv,  &
                       n_compute_fns,  &
                       i_basis_fns_inv, n_max_compute_fns_ham   )
                else
                  call evaluate_wave_gradient_cartesian_p1 &
                       (dist_tab(1,i_point), i_r(1), &
                       dir_tab(1,1,i_point), index_lm,  l_wave_max, &
                       n_compute_c, i_basis, atom_index_inv, &
                       i_basis_fns_inv, kinetic_wave(1), &
                       kinetic_wave_deriv(1), &
                       ylm_tab(1,1), &
                       sum_gradient(1,1,1), &
                       kinetic_gradient_basis_wave(1,1,i_point), &
                       n_compute_atoms)
                end if

                ! CC: Now that we have the derivative of the kinetic energy, we
                !     can compute the respective strain derivative
                if (use_analytical_stress) then
                  do AS_hessian_index=1,AS_components,1
                    call AS_eval_strain_deriv_wave(          &
                         n_compute_c, n_compute_atoms,        &
                         dist_tab(1,i_point),                 &
                         dir_tab(1,1,i_point),                &
                         i_basis(1),                          &
                         i_basis_fns_inv,                     &
                         atom_index_inv(1),                   &
                         kinetic_gradient_basis_wave(1,1,i_point),    &
                         AS_hessian_index,                    &
                         AS_strain_deriv_kinetic_wave(1,i_point,AS_hessian_index))
                  end do
                end if

                i_spin = 1
                do i_calculate_dimension= 1,3
                  call copy_gradient &
                       (d_H_times_psi(1,i_point,i_calculate_dimension,1), &
                       kinetic_gradient_basis_wave(1,1, i_point), n_compute_c, &
                       i_calculate_dimension)
                end do
              end if !flag_rel == atomic_zora

              if ((flag_rel.eq.1)) then
                ! Scalar relativistic treatment.
                ! count number of "truly" relativistic points for ZORA treatment
                ! of kinetic energy ...
                do i_spin = 1, n_spin_local, 1
                  zora_operator(i_spin) =  &
                       light_speed_sq /  &
                       (2 * light_speed_sq -  &
                       zora_potential_parts(i_spin))**2
                  call add_zora_gradient_part_p0(   &
                       sum_of_local_gradients(1:3,i_spin),  &
                       i_r_full,  &
                       dir_tab_full_norm,   &
                       dist_tab_full,  &
                       zora_operator(i_spin), &
                       n_centers_integrals, centers_basis_integrals )
                end do

                do i_spin = 1, n_spin_local, 1
                  ! Evaluate difference of scalar relativistic kinetic energy
                  ! operator for the true potential and the superposition of
                  ! free atom potentials separately, and only for all
                  ! relativistic points in shell. Here, use partially integrated
                  ! version, leading to a vector:
                  !zora_operator(r)*grad(phi(r,i))
                  zora_operator(i_spin) =  &
                       light_speed_sq *  &
                       (local_potential_parts(i_spin) -  &
                       zora_potential_parts(i_spin))/  &
                       ( 2 * light_speed_sq -  &
                       local_potential_parts(i_spin))/  &
                       ( 2 * light_speed_sq -  &
                       zora_potential_parts(i_spin))

                  call evaluate_zora_vector_p1  &
                       ( zora_operator(i_spin),  &
                       partition_tab(i_full_points_C),  &
                       gradient_basis_wave(1,1,i_point),  &
                       n_compute_c,  &
                       zora_vector1(1, 1, n_rel_points+1, i_spin),  &
                       zora_vector2(1, 1, n_rel_points+1, i_spin), &
                       n_max_compute_ham, t_zora(i_spin)  )
                enddo

                if (t_zora(1) .or. t_zora(2)) n_rel_points = n_rel_points + 1
              end if  ! end ZORA preparations
            end if ! dens_mat_is_en_weighted
          end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! i_index (loop within a batch)

        ! We now have all quantities computed on individual points for the
        ! current batch.  Time to do the actual batch integration.

        if (n_bp > 0) then
          do i_compute_1 = 1, n_compute_c
            ins_idx(i_compute_1) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i_compute_1))
          enddo
        end if ! n_bp > 0

        ! When GPU acceleration is being used, shift the needed quantities
        ! to the GPU
        if (use_gpu) then
          if (index_on_gpu) then
            ! Create arrays to permute compute basis elements so that they are
            ! ordered based on the atoms they reside on.
            ! This is needed in the GPU code to suitably vectorize the final
            ! construction of forces from the forces shell.
            i_permute = 0
            permute_compute_by_atom = 0
            n_compute_for_atom = 0
            do i_atom = 1, n_atoms
              do i_compute_1 = 1, n_compute_c
                my_atom = cBasis_to_atom(i_basis(i_compute_1))

                if (my_atom .eq. i_atom) then
                  i_permute = i_permute + 1
                  permute_compute_by_atom(i_compute_1) = i_permute
                  n_compute_for_atom(i_atom) = n_compute_for_atom(i_atom) + 1
                end if
              end do ! i_compute_1
            end do ! i_atom
          end if

          ! Quantities needed for all batch integrations
          call set_partition_gpu(partition)

          ! Quantities needed for calculation of Pulay forces, and thus
          ! always applicable
          call set_h_times_psi_gpu(H_times_psi)
          call set_d_wave_gpu(d_wave)

          ! Quantities needed both for pretty much everything except
          ! non-relativistic LDA force-only calculations
          if (gga_forces_on_temp .or. flag_rel.eq.REL_atomic_zora .or. &
              use_analytical_stress) then
            call set_wave_gpu(wave)
          end if

          ! Quantities needed only for GGA forces
          if (gga_forces_on_temp) then
            call set_hessian_basis_wave_gpu(hessian_basis_wave)
            call set_gradient_basis_wave_gpu(gradient_basis_wave)
            call set_xc_gradient_deriv_gpu(xc_gradient_deriv)

            if (meta_gga_forces_on_temp) then
              call set_xc_tau_deriv_gpu(xc_tau_deriv)
            end if
          end if

          ! Quantities needed only for atomic ZORA corrections
          if (flag_rel.eq.REL_atomic_zora) then
            call set_d_h_times_psi_gpu(d_H_times_psi)
          end if

          ! Quantities needs for calculating analytical stress tensor
          if (use_analytical_stress) then
            call set_as_strain_deriv_wave_gpu(AS_strain_deriv_wave)
            if (use_AS_Jac_in_pulay) then
              call set_as_jac_pot_kin_times_psi_gpu(AS_jac_pot_kin_times_psi)
            end if
            if (gga_forces_on_temp) then
              call set_AS_hessian_times_xc_deriv_gga_gpu &
                   (AS_hessian_times_xc_deriv_gga)
            end if
            if (meta_gga_forces_on_temp) then
              call set_AS_hessian_times_xc_deriv_mgga_gpu &
                   (AS_hessian_times_xc_deriv_mgga)
            end if
            if (flag_rel.eq.REL_atomic_zora) then
              call set_as_strain_deriv_kinetic_wave_gpu &
                  (AS_strain_deriv_kinetic_wave)
            end if
          end if

          ! Quantities needed for updating forces and analytical stress tensor
          ! on GPU; dens_mat is the same for all batches and thus was
          ! communicated at the beginning of the subroutine
          if (n_bp.gt.0) then
             if (n_bp.gt.0) call set_ins_idx_gpu(ins_idx, n_compute_c)
             call set_permute_compute_by_atom_gpu(permute_compute_by_atom, &
                                                  n_compute_c)
          end if
        end if

        ! Now add all contributions to the full forces matrix, by way of matrix
        ! multiplications work separately for each spin channel

        ! CC: 3 components for forces, 6 or 9 for stress
        if (use_analytical_stress) then
          n_calculate_dimension = AS_components
        else
          n_calculate_dimension = 3
        end if

        do i_spin = 1, n_spin_local, 1
          do i_calculate_dimension = 1,n_calculate_dimension

            ! add full non-relativistic contributions and (for relativistic
            ! points) all contributions from the potential to the Hamiltonian
            ! matrix elements

            ! First, calculated the Pulay forces corresponding to the current
            ! term, as indicated by dens_mat_is_en_weighted, for the current
            ! batch.  (Remember that, when the density matrix is energy
            ! weighted, H_times_psi actually holds psi)
            if (i_calculate_dimension .le. 3) then
              if (use_gpu) then
                call eval_forces_shell_dpsi_h_psi_gpu  &
                     ( n_points, n_compute_c, n_compute_a, &
                     i_calculate_dimension, i_spin)
              else
                call evaluate_forces_shell  &
                     ( n_points, partition(1), n_compute_c, n_compute_a, &
                     H_times_psi(1,1,i_spin),  &
                     d_wave(1,1,i_calculate_dimension ),  &
                     forces_shell, n_max_compute_ham )
              end if
            end if

            ! CC: Here comes the integration for the stress part
            if (use_analytical_stress) then
              ! Orbital derivatives
              if (use_gpu) then
                call eval_AS_shell_dpsi_h_psi_gpu  &
                     ( n_points, n_compute_c, n_compute_a, &
                     i_calculate_dimension, i_spin)
              else
                call evaluate_forces_shell  &
                     ( n_points, partition(1), n_compute_c, n_compute_a,    &
                     H_times_psi(1,1,i_spin),                               &
                     AS_strain_deriv_wave(1,1,i_calculate_dimension),       &
                     AS_strain_deriv_wave_shell(1,1), n_max_compute_ham )
              end if

              if (.not.dens_mat_is_en_weighted)  then
                !! ! Jacobian, derivative of potential, and kinetic term
                !! --> Add 0.5 * jac_pot_kin to 1.0*AS_strain_deriv_wave_shell
                if (use_AS_Jac_in_pulay) then
                  if (use_gpu) then
                    call eval_AS_shell_add_psi_kin_psi_shell_gpu &
                         ( n_points, n_compute_c, n_compute_a, &
                         i_calculate_dimension, i_spin )
                  else
                    call evaluate_forces_shell_add  &
                         ( n_points, partition(1), n_compute_c, n_compute_a,    &
                         AS_jac_pot_kin_times_psi(1,1,i_spin,i_calculate_dimension),&
                         wave(1,1),  &
                         AS_strain_deriv_wave_shell(1,1), n_max_compute_ham )
                  end if

                  if (meta_gga_forces_on_temp) then
                    if (i_calculate_dimension.eq.1 .or. &
                         i_calculate_dimension.eq.4 .or. &
                         i_calculate_dimension.eq.6 ) then
                      if (use_gpu) then
                        call evaluate_forces_shell_add_mgga_gpu  &
                             ( n_points, n_compute_c, n_compute_a, i_spin )
                      else
                        call evaluate_forces_shell_add_mgga  &
                             ( n_points, partition(1), n_compute_c, n_compute_a,  &
                             xc_tau_deriv(i_spin,1), gradient_basis_wave(1,1,1),  &
                             AS_strain_deriv_wave_shell(1,1), n_max_compute_ham )
                      end if
                    end if
                  end if
                end if

                if (flag_rel.eq.REL_none) then
                  ! Kinetic energy corrections
                  ! We do not GPU accelerate this operation, as it is made up
                  ! of vector-vector operations for relatively small vectors
                  call AS_evaluate_forces_shell_unified_corrections(          &
                       n_points, partition(1), n_compute_c,                           &
                       AS_n_on_site, AS_i_basis_on_site(1,1), AS_i_basis_on_site(1,2),&
                       AS_T_times_psi(1,1),                                           &
                       AS_strain_deriv_wave(1,1,i_calculate_dimension),               &
                       AS_Hessian_times_psi(1,1,i_calculate_dimension),               &
                       wave(1,1),                                                     &
                       AS_corr_shell(1,1), n_max_compute_ham )
                end if

                if (gga_forces_on_temp) then
                  ! Add GGA part to AS_strain_deriv_wave_shell.
                  ! Updated to also include addition of mGGA.
                  if (use_gpu) then
                    call AS_evaluate_gga_stress_gpu(n_compute_c, n_points, &
                         i_calculate_dimension, i_spin)
                  else
                    call AS_evaluate_gga_stress(                                         &
                         AS_hessian_times_xc_deriv_gga(1,1,i_spin,i_calculate_dimension), &
                         AS_hessian_times_xc_deriv_mgga(1,1,i_spin,i_calculate_dimension),&
                         xc_gradient_deriv(1,1,1),                                        &
                         AS_strain_deriv_wave(1,1,i_calculate_dimension),                 &
                         gradient_basis_wave(1,1,1),                                      &
                         wave(1,1),                                                       &
                         partition(1),                                                    &
                         n_max_compute_ham,                                               &
                         n_compute_c,                                                     &
                         n_points,                                                        &
                         i_spin,                                                          &
                         meta_gga_forces_on_temp,                                         &
                         AS_strain_deriv_wave_shell(1,1)                                  )
                  end if
                end if
              end if !.not.dens_mat_is_en_weighted
            end if !analytic_stress

            ! Calculate GGA and meta-GGA forces for the current batch, if
            ! applicable
            if( (.not. dens_mat_is_en_weighted) .and. &
                  gga_forces_on_temp .and. &
                  (i_calculate_dimension .le. 3) ) then
              if (use_gpu) then
                call eval_gga_forces_dens_mat_gpu ( &
                      n_compute_c, n_points, i_calculate_dimension, i_spin)
              else
                call evaluate_gga_forces_dens_mat(forces_shell, &
                     hessian_basis_wave, &
                     xc_gradient_deriv, gradient_basis_wave, &
                     wave, n_compute_c, n_points, &
                     i_calculate_dimension, i_spin, partition(1), &
                     xc_tau_deriv, meta_gga_forces_on_temp)
              end if
            end if

            ! Calculate Gnonmf forces for the current batch, if applicable
            if ( (.not.dens_mat_is_en_weighted) .and.&
                   Gnonmf_forces_on_temp.and.&
                   (i_calculate_dimension .le. 3) ) then
              call evaluate_Gnonmf_forces_dens_mat(forces_shell, &
                   hessian_basis_wave, &
                   Gnonmf_gradient_deriv, gradient_basis_wave, &
                   wave, n_compute_c, n_points, &
                   i_calculate_dimension, i_spin, partition(1))
            end if

            ! Update the forces with the values of all previous contributions
            ! for the current batch
            if (dens_mat_is_en_weighted) then
              do i_spin_2 = 1, n_spin
                if (use_analytical_stress) then
                  ! No corrections needed of dens_mat_is_en_weighted, even in the
                  ! non-relativistic case
                  call AS_update_sum_forces_and_stress &
                       ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                       n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat, &
                       i_spin_2, i_calculate_dimension, sum_forces_temp, &
                       AS_strain_deriv_wave_shell(1,1),     &
                       dens_mat_is_en_weighted )
                else
                  call update_sum_forces &
                       ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                       n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat,&
                       i_spin_2, i_calculate_dimension, sum_forces_temp )
                end if
              end do
            else
              if (use_analytical_stress) then
                if (flag_rel.eq.REL_none) then
                  call AS_update_sum_forces_and_stress &
                       ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                       n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat, &
                       i_spin, i_calculate_dimension, sum_forces_temp, &
                       AS_strain_deriv_wave_shell(1,1),     &
                       (flag_rel.eq.REL_atomic_zora), &  ! Compute corrections
                       AS_corr_shell(1,1)                   )
                else
                  call AS_update_sum_forces_and_stress &
                       ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                       n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat, &
                       i_spin, i_calculate_dimension, sum_forces_temp, &
                       AS_strain_deriv_wave_shell(1,1),     &
                       (flag_rel.eq.REL_atomic_zora)) ! No corrections
                end if
              else
                call update_sum_forces &
                     ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                     n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat, &
                     i_spin, i_calculate_dimension, sum_forces_temp )
              end if ! use_analytical_stress
            end if ! dens_mat_is_en_weighted

            ! Finally, calculate the relativistic correction due to atomic ZORA
            ! for the current batch and add it in, if needed.  We calculate and
            ! add this contribution seperately from the other contributions,
            ! since we must transpose forces_shell in this case.
            if ((flag_rel.eq.REL_atomic_zora) .and. &
                 (.not. dens_mat_is_en_weighted)) then
              if (i_calculate_dimension .le. 3) then
                if (use_gpu) then
                  call eval_forces_shell_psi_dH_psi_gpu  &
                       ( n_points, n_compute_c, n_compute_a, &
                       i_calculate_dimension, &
                       i_spin)
                else
                  call evaluate_forces_shell  &
                       ( n_points, partition(1), n_compute_c, &
                       n_compute_a, &
                       d_H_times_psi(1,1,i_calculate_dimension, 1),  &
                       wave(1,1),  &
                       forces_shell, n_max_compute_ham )
                  call transpose_forces_shell( &
                       forces_shell, n_compute_c)
                end if
              end if

              !! CC: Now calculate the respective shell for the strain
              !!     derivative of the (rel.) kinetic energy
              if ( use_analytical_stress ) then
                if (use_gpu) then
                  call eval_as_shell_psi_dkin_psi_gpu  &
                       ( n_points, n_compute_c, n_compute_a, &
                       i_calculate_dimension, &
                       i_spin)
                  call transpose_as_shell_gpu(n_compute_c)
                else
                  call evaluate_forces_shell  &
                       ( n_points, partition(1), n_compute_c, n_compute_a, &
                       AS_strain_deriv_kinetic_wave(1,1,i_calculate_dimension), &
                       wave(1,1),  &
                       AS_strain_deriv_wave_shell(1,1), n_max_compute_ham )
                  call transpose_forces_shell( AS_strain_deriv_wave_shell(1,1), &
                       n_compute_c)
                end if

                call AS_update_sum_forces_and_stress &
                     ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                     n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat, &
                     i_spin, i_calculate_dimension, sum_forces_temp, &
                     AS_strain_deriv_wave_shell(1,1),     &
                     .true.) ! No corrections, since relativistic case
              else
                call update_sum_forces &
                     ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                     n_compute_for_atom, forces_shell, dens_mat, ld_dens_mat, &
                     i_spin, i_calculate_dimension, sum_forces_temp )
              end if ! use_analytical_stress
            end if ! i_flag_rel_eq.REL_atomic_ZORA.and..not. ...
          end do ! i_calculate_dimension
        end do ! i_spin
        ! Forces are now complete for the current batch.
      else ! n_compute_C.eq.0
        do i_index = 1, batches_work(i_my_batch)%size, 1
          i_full_points_C = i_full_points_C + 1
          i_full_points_A = i_full_points_A + 1
        enddo
      end if ! n_compute.gt.0
    end do ! end loop over batches

    ! YY symmetry
    if (use_rlsy_symmetry) then
      !unfold_force: block
      do i_op = 1, rlsy_h%spacegroup%n_operation
        do i_atom = 1, n_atoms
          j_atom = rlsy_h%spacegroup%op(i_op)%fmap(i_atom)
          sum_forces(:,j_atom) = sum_forces(:,j_atom) + &
              matmul(rlsy_h%spacegroup%op(i_op)%m, sum_forces_temp(:,i_atom))
        end do
      end do
      !end block unfold_force

      deallocate(rho)
      deallocate(rho_gradient)
      deallocate(kinetic_density)
      deallocate(ins_idx)
      deallocate(hartree_potential)
    else
        sum_forces(:,:) = sum_forces(:,:) + sum_forces_temp(:,:) 
    end if
    ! YY symmetry

    ! Finalize GPU acceleration, when applicable
    if (use_gpu) then
      ! Because we did indexing on the GPU for PM_none, now we communicate the
      ! results back to the GPU
      if (index_on_gpu) then
        call get_sum_forces_gpu(sum_forces, n_atoms)
        if (use_analytical_stress) then
          call get_as_pulay_stress_local_gpu(AS_pulay_stress_local)
        end if
      end if

      write (info_str,'(A)') &
          'Deallocating memory used on GPU for evaluation of forces.'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      call forces_destroy_gpu()
    end if
    use_gpu = gpu_save

    if (get_batch_weights) then
      call set_batch_weights(n_bp, batch_times)
    endif

    ! Allocatable arrays that are tracked
    if (allocated( hessian_basis_wave         )) call aims_deallocate( hessian_basis_wave,                   "hessian_basis_wave" )
    if (allocated( gradient_basis_wave        )) call aims_deallocate( gradient_basis_wave,                 "gradient_basis_wave" )
    if (allocated( kinetic_gradient_basis_wave)) call aims_deallocate( kinetic_gradient_basis_wave, "kinetic_gradient_basis_wave" )
    if (allocated( forces_shell               )) call aims_deallocate( forces_shell,                               "forces_shell" )
    if (allocated( H_times_psi                )) call aims_deallocate( H_times_psi,                                 "H_times_psi" )
    if (allocated( d_H_times_psi              )) call aims_deallocate( d_H_times_psi,                             "d_H_times_psi" )
    if (allocated( wave                       )) call aims_deallocate( wave,                                               "wave" )
    if (allocated( d_wave                     )) call aims_deallocate( d_wave,                                           "d_wave" )
    if (allocated( zora_vector1               )) call aims_deallocate( zora_vector1,                               "zora_vector1" )
    if (allocated( zora_vector2               )) call aims_deallocate( zora_vector2,                               "zora_vector2" )

    if (allocated( i_r_full                   )) deallocate( i_r_full                   )
    if (allocated( dir_tab_full_norm          )) deallocate( dir_tab_full_norm          )
    if (allocated( dist_tab_full              )) deallocate( dist_tab_full              )
    if (allocated( i_basis                    )) deallocate( i_basis                    )
    if (allocated( kinetic_wave_deriv         )) deallocate( kinetic_wave_deriv         )
    if (allocated( kinetic_wave               )) deallocate( kinetic_wave               )
    if (allocated( radial_wave_deriv          )) deallocate( radial_wave_deriv          )
    if (allocated( radial_wave                )) deallocate( radial_wave                )
    if (allocated( index_lm                   )) deallocate( index_lm                   )
    if (allocated( ylm_tab                    )) deallocate( ylm_tab                    )
    if (allocated( scaled_dylm_dphi_tab       )) deallocate( scaled_dylm_dphi_tab       )
    if (allocated( dylm_dtheta_tab            )) deallocate( dylm_dtheta_tab            )
    if (allocated( radial_wave_2nd_deriv      )) deallocate( radial_wave_2nd_deriv      )
    if (allocated( sum_hessian                )) deallocate( sum_hessian                )
    if (allocated( sum_gradient               )) deallocate( sum_gradient               )
    if (allocated( cartesians                 )) deallocate( cartesians                 )
    if (use_analytical_stress) then
      if (allocated(AS_strain_deriv_wave))         deallocate(AS_strain_deriv_wave)
      if (allocated(AS_strain_deriv_wave_shell))   deallocate(AS_strain_deriv_wave_shell)
      if (allocated(AS_jac_pot_kin_times_psi))     deallocate(AS_jac_pot_kin_times_psi)
      if (allocated(AS_Hessian_times_psi))         deallocate(AS_Hessian_times_psi)
      if (allocated(AS_T_times_psi))               deallocate(AS_T_times_psi)
      if (allocated(AS_corr_shell))                deallocate(AS_corr_shell)
      if (allocated(AS_i_basis_on_site))           deallocate(AS_i_basis_on_site)
      if (allocated(AS_strain_deriv_kinetic_wave)) deallocate(AS_strain_deriv_kinetic_wave)
      if (gga_forces_on_temp) then
        if (allocated(AS_hessian_times_xc_deriv_gga)) deallocate(AS_hessian_times_xc_deriv_gga)
        if (allocated(AS_hessian_times_xc_deriv_mgga)) deallocate(AS_hessian_times_xc_deriv_mgga)
      end if
    end if
    if (Gnonmf_forces_on_temp) then
      if (allocated(Gnonmf_gradient_deriv))  deallocate(Gnonmf_gradient_deriv)
      if (allocated(local_Gnonmf_derivs))  deallocate(local_Gnonmf_derivs)
    end if

  end subroutine integrate_force_integrals_densmat
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/copy_gradient
!  NAME
!    copy_gradient
!  SYNOPSIS
  subroutine copy_gradient(wave, gradient, n_compute, i_dim)
!  PURPOSE
!    Copies values from variable gradient to wave from asked dimension.
!    This is needed, because the dimensions of the gradient tables are different
!    in subroutines (n_compute,3) than in calling subroutine (n_max_compute).
!  USES
    implicit none
!  ARGUMENTS
    real*8,  intent(out) :: wave(n_compute)
    real*8,  intent(in)  :: gradient(n_compute, 3)
    integer, intent(in)  :: n_compute
    integer, intent(in)  :: i_dim
!  INPUTS
!    o  gradient -- gradient table, from where the values are copied
!    o  n_compute -- the first dimension of the gradient table
!    o  i_dim --- the dimension which is copied
!  OUTPUT
!    o  wave -- the table where the values are copied
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    wave(1:n_compute) = gradient(1:n_compute,i_dim)
  end subroutine copy_gradient
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/transpose_forces_shell
!  NAME
!    transpose_forces_shell
!  SYNOPSIS
  subroutine transpose_forces_shell( forces_shell, n_compute )
!  PURPOSE
!    Makes traspose of the table.
!    This is needed, because the dimensions of the gradient tables are different
!    in subroutines (n_compute,3)  than in calling subroutine (n_max_compute).
!  USES
    implicit none
!  ARGUMENTS
    real*8,  intent(inout) :: forces_shell( n_compute,n_compute)
    integer, intent(in)    :: n_compute
!  INPUTS
!    o forces_shell -- table before transpose
!    o n_compute -- dimension of the table
!  OUTPUT
!    o forces_shell -- table after transpose
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    forces_shell = transpose(forces_shell)
  end subroutine transpose_forces_shell
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/evaluate_forces_shell
!  NAME
!    evaluate_forces_shell
!  SYNOPSIS
  subroutine evaluate_forces_shell &
       ( n_points, partition, n_compute_c, n_compute_a,  H_times_psi, &
       wave, forces_shell, n_basis_list )
!  PURPOSE
!    Subroutine evaluates the wave * H_times_wave  integral contribution
!    of several integration points, and adds it to the forces_shell
!    This routine is simular than evaluate_hamiltonian_shell, but here
!    the matrix is not symmetric and is not symmetriced.
!  USES
    implicit none
!  ARGUMENTS
    integer, intent(in)  :: n_points
    real*8,  intent(in)  :: partition(n_points)
    integer, intent(in)  :: n_compute_c
    integer, intent(in)  :: n_compute_a
    real*8,  intent(in)  :: H_times_psi(n_basis_list, n_points)
    real*8,  intent(in)  :: wave(n_basis_list, n_points)
    real*8,  intent(out) :: forces_shell(n_compute_c, n_compute_a)
    integer, intent(in)  :: n_basis_list
!  INPUTS
!    o n_points -- number of integration grid points
!    o partition -- values of partition function
!    o n_compute_c -- the number of relevantbasis functions
!    o n_compute_a  -- the number of relevantbasis functions
!    o H_times_psi -- hamiltonian times basis function
!    o wave -- basis functions
!    o n_basis_list -- the maximum number of basis function including periodic
!      mirror images.
!  OUTPUT
!    o forces_shell -- results of the integrals over different grid points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    ! auxiliary matrices for Level 3 Blas matrix multiplications

    ! ??? Need we store partition * wave in an extra matrix? ??? need BLAS
    ! ??? descriptions to go on ...
    ! ??? Does it matter that n_points (the summation index in the matrix
    ! ??? multiplication) is the outer one, not the inner one, in wave?
    !     Yes, we need, but as long as we are inside the cache, no harm should
    !     follow...
    real*8 wave_compute_a(n_compute_a,n_points), matrix_sum

    ! counters
    integer :: i_point, i_compute

    ! First we take only the waves that are used at these points
    do i_point = 1, n_points, 1
       wave_compute_a(1:n_compute_a, i_point) =  &
            partition(i_point)*wave(1:n_compute_a, i_point)
    enddo

    ! now, integrate, initialisation of aux. ham. matrix is needed
    ! since we add to it the values in the second stage below

    ! compute wave*(H*psi) and add this to aux. ham. matrix
    call dgemm('N', 'T', n_compute_c, n_compute_a,  &
               n_points, 1.0d0,  &
               H_times_psi, n_basis_list, &
               wave_compute_a, n_compute_a, &
               0.0d0, forces_shell, n_compute_c)
  end subroutine evaluate_forces_shell
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/evaluate_forces_shell_add
!  NAME
!    evaluate_forces_shell_add
!  SYNOPSIS
  subroutine evaluate_forces_shell_add &
       ( n_points, partition, n_compute_c, n_compute_a,  H_times_psi, &
       wave, forces_shell, n_basis_list )
!  PURPOSE
!    CC: Same as evaluate_forces_shell, but
!        - we add 0.5*forces_shell to the old forces_shell
!        - hence we do not zero forces_shell
!    Subroutine evaluates the wave * H_times_wave  integral contribution
!    of several integration points, and adds it to the forces_shell
!    This routine is simular than evaluate_hamiltonian_shell, but here
!    the matrix is not symmetric and is not symmetriced.
!  USES
    implicit none
!  ARGUMENTS
    integer, intent(in)    :: n_points
    real*8,  intent(in)    :: partition(n_points)
    integer, intent(in)    :: n_compute_c
    integer, intent(in)    :: n_compute_a
    real*8,  intent(in)    :: H_times_psi(n_basis_list, n_points)
    real*8,  intent(in)    :: wave(n_basis_list, n_points)
    real*8,  intent(inout) :: forces_shell(n_compute_c, n_compute_a)
    integer, intent(in)    :: n_basis_list
!  INPUTS
!    o n_points -- number of integration grid points
!    o partition -- values of partition function
!    o n_compute_c -- the number of relevantbasis functions
!    o n_compute_a  -- the number of relevantbasis functions
!    o H_times_psi -- hamiltonian times basis function
!    o wave -- basis functions
!    o n_basis_list -- the maximum number of basis function including periodic
!      mirror images.
!  OUTPUT
!    o forces_shell -- results of the integrals over different grid points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    real*8 wave_compute_a(n_compute_a,n_points)

    !     counters
    integer :: i_point

    !     First we take only the waves that are used at these points
    do i_point = 1, n_points, 1
      wave_compute_a(1:n_compute_a, i_point) =  &
           partition(i_point)*wave(1:n_compute_a, i_point)
    enddo

    call dgemm('N', 'T', n_compute_c, n_compute_a,  &
               n_points, 0.5d0,  &
               H_times_psi, n_basis_list, &
               wave_compute_a, n_compute_a, &
               1.0d0, forces_shell, n_compute_c )
  end subroutine evaluate_forces_shell_add
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/evaluate_forces_shell_add_mgga
!  NAME
!    evaluate_forces_shell_add_mgga
!  SYNOPSIS
  subroutine evaluate_forces_shell_add_mgga &
       ( n_points, partition, n_compute_c, n_compute_a, xc_tau_deriv, &
       gradient_basis_wave, AS_strain_deriv_wave_shell, n_basis_list )
!  PURPOSE
!    CC: Same as evaluate_forces_shell_add, but
!        - we add 0.5*AS_strain_deriv_wave_shell to the old
!          AS_strain_deriv_wave_shell
!        - hence we do not zero AS_strain_deriv_wave_shell
!    Subroutine evaluates the integral of:
!    gradient_basis_wave * xc_tau_deriv * partition * gradient_basis_wave
!    of several integration points, and adds it to the !
!    AS_strain_deriv_wave_shell
!    This routine is similar to evaluate_forces_shell_add, and
!    is written as such to re-use the machinery implemented for use
!    in the SCF cycle through the subroutines of
!    evaluate_mgga_contribution_and_add_to_hamiltonian_shell
!  USES
    implicit none
!  ARGUMENTS
    integer, intent(in)    :: n_points
    real*8,  intent(in)    :: partition(n_points)
    integer, intent(in)    :: n_compute_c
    integer, intent(in)    :: n_compute_a
    real*8,  intent(in)    :: xc_tau_deriv(n_points)
    real*8,  intent(in)    :: gradient_basis_wave(n_basis_list, 3, n_points)
    real*8,  intent(inout) :: AS_strain_deriv_wave_shell(n_compute_c, n_compute_a )
    integer, intent(in)    :: n_basis_list
!  INPUTS
!    o n_points -- number of integration grid points
!    o partition -- values of partition function
!    o n_compute_c -- the number of relevantbasis functions
!    o n_compute_a  -- the number of relevantbasis functions
!    o xc_tau_deriv -- derivative of the exchange-correlation wrt tau
!    o gradient_basis_wave -- gradient of the basis functions
!    o n_basis_list -- the maximum number of basis function including periodic
!      mirror images.
!  OUTPUT
!    o AS_strain_deriv_wave_shell -- results of the integrals over different grid points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    ! local variables:
    real*8 left_side_of_mgga_dot_product(n_compute_c, 3*n_points)
    real*8 gradient_basis_wave_compute_a(n_compute_a, 3*n_points)

    ! counters
    integer :: i_point

    ! First we take need to calculate the parts of the dot product
    ! using only the appropriate basis functions at these points
    do i_point = 1, n_points, 1
      ! Includes the prefactor of -0.5, for want of a better place,
      ! similar to the GGA contribution to the hamiltonian
      call evaluate_mgga_left_side_of_dot_product &
         ( n_compute_c, n_points, i_point, &
         -0.5d0*partition(i_point), xc_tau_deriv(i_point), &
         gradient_basis_wave(1,1,i_point), &
         left_side_of_mgga_dot_product(1,1) )
      call store_gradient_basis_wave &
           ( n_compute_a, n_points, i_point, &
           gradient_basis_wave(1,1,i_point), &
           gradient_basis_wave_compute_a(1,1) )
    enddo

    !call dgemm('N', 'T', n_compute_c, n_compute_a,  &
    !     3*n_points, -0.5d0,  &
    !     left_side_of_mgga_dot_product, n_compute_c,
    !     gradient_basis_wave_compute_a, n_compute_a, 1.0d0, AS_strain_deriv_wave_shell, &
    !     n_compute_c )
    call evaluate_mgga_contribution_and_add_to_hamiltonian_shell  &
         ( n_compute_c, n_compute_a, n_points, &
         left_side_of_mgga_dot_product(1,1), &
         gradient_basis_wave_compute_a(1,1), &
         AS_strain_deriv_wave_shell )
  end subroutine evaluate_forces_shell_add_mgga
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/AS_evaluate_forces_shell_unified_corrections
!  NAME
!    AS_evaluate_forces_shell_unified_corrections
!  SYNOPSIS
  subroutine AS_evaluate_forces_shell_unified_corrections( &
         n_points, partition, n_compute_c,                         &
         n_on_site, i_basis_on_site_A, i_basis_on_site_C,          &
         H_times_psi_A, wave_A,                                    &
         H_times_psi_B, wave_B,                                    &
         AS_corr_shell, n_basis_list )
!  PURPOSE
!    Again, in principle the same as evaluate_forces_shell, but
!    - we perform two integrations at at once and add everything up with the
!      correct prefactors
!    - we restrict ourself to on-site terms, hence DGEMM -> sum_n_on_site
!      ddot(..,..)
!  USES
    implicit none
!  ARGUMENTS
    integer, intent(in)  :: n_points
    real*8,  intent(in)  :: partition(n_points)
    integer, intent(in)  :: n_compute_c
    integer, intent(in)  :: n_on_site
    integer, intent(in)  :: i_basis_on_site_A(n_on_site)
    integer, intent(in)  :: i_basis_on_site_C(n_on_site)
    real*8,  intent(in)  :: H_times_psi_A(n_basis_list, n_points)
    real*8,  intent(in)  :: wave_A(n_basis_list, n_points)
    real*8,  intent(in)  :: H_times_psi_B(n_basis_list, n_points)
    real*8,  intent(in)  :: wave_B(n_basis_list, n_points)
    real*8,  intent(out) :: AS_corr_shell(n_compute_c, n_compute_c )
    integer, intent(in)  :: n_basis_list
!  SOURCE

    real*8 wave_compute_A(n_points,n_compute_c)
    real*8 wave_compute_B(n_points,n_compute_c)
    real*8 H_times_psi_compute_A(n_points,n_compute_c)
    real*8 H_times_psi_compute_B(n_points,n_compute_c)

    !     counters
    integer :: i_point
    integer :: i_compute_a,i_compute_c,i_on_site
    real*8, external :: ddot

    do i_point = 1, n_points, 1
       wave_compute_A(i_point,1:n_compute_c) =  &
            partition(i_point)*wave_A(1:n_compute_c, i_point)
       wave_compute_B(i_point,1:n_compute_c) =  &
            partition(i_point)*wave_B(1:n_compute_c, i_point)
       H_times_psi_compute_A(i_point,1:n_compute_c) = &
            H_times_psi_A(1:n_compute_c, i_point)
       H_times_psi_compute_B(i_point,1:n_compute_c) = &
            H_times_psi_B(1:n_compute_c, i_point)
    enddo

    do i_on_site=1,n_on_site
      i_compute_a = i_basis_on_site_A(i_on_site)
      i_compute_c = i_basis_on_site_C(i_on_site)
      AS_corr_shell(i_compute_c,i_compute_a) = &
           ddot( n_points, H_times_psi_compute_A(1,i_compute_c), 1, wave_compute_A(1,i_compute_a), 1 ) + &
           0.5d0*ddot( n_points, H_times_psi_compute_B(1,i_compute_c), 1, wave_compute_B(1,i_compute_a), 1 )
    end do
  end subroutine AS_evaluate_forces_shell_unified_corrections
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/AS_evaluate_gga_stress
!  NAME
!    AS_evaluate_gga_stress
!  SYNOPSIS
  subroutine AS_evaluate_gga_stress(            &
               hessian_times_xc_gradient_deriv, &
               hessian_times_xc_tau_deriv,      &
               xc_gradient_deriv,               &
               strain_deriv_wave,               &
               gradient_basis_wave,             &
               wave,                            &
               partition_tab,                   &
               n_basis_list,                    &
               n_compute_c,                     &
               n_points,                        &
               i_spin,                          &
               meta_gga_forces_on,              &
               AS_strain_deriv_wave_shell       )
!  PURPOSE
!    Calculating and integrating of the two GGA stress terms.
!    The result is added to the variable AS_strain_deriv_wave_shell.
!  USES
    use dimensions
    implicit none
!  ARGUMENTS
    real*8, dimension(1:n_basis_list,1:n_points),    intent(in)  :: hessian_times_xc_gradient_deriv
    real*8, dimension(1:n_basis_list,1:n_points),    intent(in)  :: hessian_times_xc_tau_deriv
    real*8, dimension(1:3,1:n_spin,1:n_points),      intent(in)  :: xc_gradient_deriv
    real*8, dimension(1:n_basis_list,1:n_points),    intent(in)  :: strain_deriv_wave
    real*8, dimension(1:n_basis_list*3,1:n_points),  intent(in)  :: gradient_basis_wave
    real*8, dimension(1:n_basis_list,1:n_points),    intent(in)  :: wave
    real*8, dimension(1:n_points),                   intent(in)  :: partition_tab
    integer,                                         intent(in)  :: n_basis_list
    integer,                                         intent(in)  :: n_compute_c
    integer,                                         intent(in)  :: n_points
    integer,                                         intent(in)  :: i_spin
    logical,                                         intent(in)  :: meta_gga_forces_on
    real*8, dimension(1:n_compute_c,1:n_compute_c),  intent(out) :: AS_strain_deriv_wave_shell
!  SOURCE

    ! local
    integer :: i_point
    integer :: i_coord
    integer :: start_index
    integer :: end_index
    real*8, dimension(1:n_compute_c,1:n_points) :: matrix_term
    real*8, dimension(1:n_compute_c,1:n_points) :: matrix_term_mgga

    !! First term (A)
    ! Calculate
    !    grad(<j|) * (d/dgrad(rho)^2 f_xc) * grad(rho) * (d/dr_l |i>) * (r-R_m)
    ! There is a scalarproduct between grad(<j|) and grad(rho).
    ! The variable xc_gradient_deriv includes grad(rho).
    ! This means, we calculate:
    !      gradient_basis_wave * xc_gradient_deriv * strain_deriv_wave
    ! Only times 2 because another factor of two is inside the density matrix,
    ! see scaled_occs in forces_densmat.
    matrix_term(:,:) = 0.0d0
    matrix_term_mgga(:,:) = 0.0d0
    do i_point = 1, n_points, 1
      do i_coord = 1, 3, 1
        ! Need start and end index corresponding to i_coord of
        ! gradient_basis_wave, because gradient_basis_wave must be initialized
        ! in this routine as a 2-dim array (3*n_compute_c,n_points) instead of a
        ! intuitive 3-dim array (n_compute_c,3,n_points):
        !      i_coord=1:               1 ...   n_compute_c
        !      i_coord=2:   n_compute_c+1 ... 2*n_compute_c
        !      i_coord=3: 2*n_compute_c+1 ... 3*n_compute_c
        start_index = (i_coord-1) * n_compute_c + 1
        end_index   =  i_coord    * n_compute_c

        matrix_term(1:n_compute_c,i_point) = matrix_term(1:n_compute_c,i_point)&
             + 2 * gradient_basis_wave(start_index:end_index,i_point) &
             * xc_gradient_deriv(i_coord,i_spin,i_point)*partition_tab(i_point)

        ! Prepare right side of dot product for mgga contribution to stress
        ! tensor
        if (meta_gga_forces_on) then
          matrix_term_mgga(1:n_compute_c,i_point) = &
              matrix_term_mgga(1:n_compute_c,i_point) &
              + gradient_basis_wave(start_index:end_index,i_point) * partition_tab(i_point)
        endif
      end do
    end do

    call dgemm('N', 'T', n_compute_c, n_compute_c,  &
               n_points, 1.0d0,  &
               strain_deriv_wave, n_basis_list, &
               matrix_term, n_compute_c, &
               1.0d0, AS_strain_deriv_wave_shell, n_compute_c )

    !! Second term (B)
    ! Calculate wave * hessian_times_xc_gradient_deriv
    ! Only times 2 because another factor of two is inside the density matrix,
    ! see scaled_occs in forces_densmat.
    do i_point = 1, n_points, 1
      matrix_term(1:n_compute_c,i_point) = &
           2 * wave(1:n_compute_c,i_point) * partition_tab(i_point)
    end do

    call dgemm('N', 'T', n_compute_c, n_compute_c,  &
               n_points, 1.0d0,  &
               hessian_times_xc_gradient_deriv, n_basis_list, &
               matrix_term, n_compute_c, &
               1.0d0, AS_strain_deriv_wave_shell, n_compute_c )

    !! Third term (C)
    ! Calculate gradient_wave * hessian_times_xc_tau_deriv
    ! No need for a factor here as there is already a factor of 2 in the density
    ! matrix, and we'll re-use the dgemm structure for the above matrix
    ! multiplication
    if (meta_gga_forces_on) then
      call dgemm('N', 'T', n_compute_c, n_compute_c,  &
                  n_points, 1.0d0,  &
                  hessian_times_xc_tau_deriv, n_basis_list, &
                  matrix_term_mgga, n_compute_c, &
                  1.0d0, AS_strain_deriv_wave_shell, n_compute_c )
    end if
  end subroutine AS_evaluate_gga_stress
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/evaluate_gga_forces_dens_mat
!  NAME
!    evaluate_gga_forces_dens_mat
!  SYNOPSIS
  subroutine evaluate_gga_forces_dens_mat( &
        forces_shell, hessian_basis_wave, xc_gradient_deriv, &
        gradient_basis_wave, wave, n_compute, n_points, &
        i_dim, i_spin, partition_tab, xc_tau_deriv, meta_gga_forces_on)
!  PURPOSE
!    The subroutine evaluates gga AND meta-gga force terms for severel
!    integration points, using the density matrix formalism.
!    So this is the gga version of evaluate_forces_shell.
!  USES
    use dimensions
    implicit none
!  ARGUMENTS
    real*8,                                             intent(inout) :: forces_shell(n_compute, n_compute)
    real*8, dimension(n_max_compute_dens, 6, n_points), intent(in)    :: hessian_basis_wave
    real*8, dimension(3, n_spin, n_points),             intent(in)    :: xc_gradient_deriv
    real*8, dimension(n_max_compute_ham*3, n_points),   intent(in)    :: gradient_basis_wave
    real*8, dimension(n_max_compute_ham, n_points),     intent(in)    :: wave
    integer,                                            intent(in)    :: n_compute
    integer,                                            intent(in)    :: n_points
    integer,                                            intent(in)    :: i_dim
    integer,                                            intent(in)    :: i_spin
    real*8,                                             intent(in)    :: partition_tab(n_points)
    real*8, dimension(n_spin, n_points),                intent(in)    :: xc_tau_deriv
    logical,                                            intent(in)    :: meta_gga_forces_on
!  INPUTS
!    o hessian_basis_wave -- hessian of the basis waves
!    o xc_gradient_deriv -- xc gradient terms
!    o gradient_basis_wave -- gradients of the basis functions
!    o wave -- the basis functions
!    o n_compute -- the number of basis functions in the grid patch
!    o n_points -- the number of grid points in the grid patch
!    o i_dim -- the dimension where forces are calculated
!    o partition_tab -- values of the partition function
!    o xc_tau_deriv -- xc derivative wrt tau, for mGGA
!    o meta_gga_forces_on -- do we need to calculate the mGGA term
!  OUTPUT
!    o forces_shell -- the results of the force components are added here.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    integer :: i_index(3,3)
    integer :: i_coord_1, i_coord_2, i_basis_1, i_basis_2
    integer :: i_point, i_counter, i_compute
    real*8  :: matrix_tem1(n_compute, n_points)
    real*8  :: matrix_tem2(n_points, n_compute)

    ! mGGA temporary storage arrays
    real*8  :: matrix_tem1_mgga(n_compute, n_points)
    real*8  :: matrix_tem2_mgga(n_points, n_compute)

    i_index = 0
    i_counter = 0

    ! The Hessian has nine overall components, but since it is symmetric we only
    ! store the six independent components.  i_index is the indexing matrix
    ! between the overall and independent components.
    ! First, initialize i_index as
    ! 1 0 0
    ! 2 4 0
    ! 3 5 6
    do i_coord_1 = 1,3
      i_counter = i_counter +1
      i_index(i_coord_1, i_coord_1) = i_counter
      do i_coord_2 = i_coord_1+1,3
        i_counter = i_counter +1
        i_index(i_coord_1, i_coord_2) = i_counter
      end do
    end do
    ! Then symmetrizes i_index to
    ! 1 2 3
    ! 2 4 5
    ! 3 5 6
    do i_coord_1 = 1,3
      do i_coord_2 = 1,3
        i_index(i_coord_1, i_coord_2) = &
             max(i_index(i_coord_1, i_coord_2), i_index(i_coord_2, i_coord_1))
      end do
    end do
    ! It would have been easier to just explicitly assign them...

    do i_basis_2 = 1, n_compute
      do i_point = 1, n_points
        matrix_tem2(i_point, i_basis_2) = &
             2d0 * gradient_basis_wave(i_basis_2+(i_dim-1)*n_compute, i_point) &
             * partition_tab(i_point)
      end do
    end do

    do i_coord_2 = 1,3
      do i_basis_2 = 1, n_compute
        do i_point = 1, n_points
          matrix_tem1(i_basis_2, i_point) = &
               gradient_basis_wave(i_basis_2+(i_coord_2-1)*n_compute, i_point) &
               * xc_gradient_deriv(i_coord_2, i_spin, i_point)
        end do
      end do
      call dgemm('N', 'N', n_compute, n_compute,  &
           n_points, 1.d0,  &
           matrix_tem1, n_compute, &
           matrix_tem2, n_points, &
           1.d0, forces_shell, n_compute )
    end do

    do i_basis_2 = 1, n_compute
      do i_point = 1, n_points
        matrix_tem1(i_basis_2, i_point) = &
             2d0 * wave(i_basis_2, i_point)* partition_tab(i_point)
      end do
    end do

    do i_coord_2 = 1,3
      do i_basis_2 = 1, n_compute
        do i_point = 1, n_points
          matrix_tem2(i_point, i_basis_2) = &
               hessian_basis_wave(i_basis_2, i_index(i_dim,i_coord_2), i_point)&
               * xc_gradient_deriv(i_coord_2, i_spin, i_point)

          ! Calculate meta_gga force contributions. AJL
          ! Need to work out gradient_basis_wave*weights and
          ! hessian_basis_wave*xc_tau_deriv,, and then we need to dot product
          ! them for each required coordinate.
          ! Dot product is calculated further down:
          !      [ d(f_xc) / d(tau) ] * grad_at(grad(phi)) . grad(phi)
          ! No factor of 2 as this is in the density matrix. Same as that gga
          ! contribution has factor of 2 not 4.
          if (meta_gga_forces_on) then
            matrix_tem1_mgga(i_basis_2, i_point) = &
                 gradient_basis_wave(i_basis_2+(i_coord_2-1)*n_compute,i_point)&
                 * partition_tab(i_point)

            matrix_tem2_mgga(i_point, i_basis_2) = &
                 hessian_basis_wave(i_basis_2,i_index(i_dim,i_coord_2),i_point)&
                 * xc_tau_deriv(i_spin, i_point)
          end if
        end do
      end do

      call dgemm('N', 'N', n_compute, n_compute,  &
                 n_points, 1.d0,  &
                 matrix_tem1, n_compute, matrix_tem2, &
                 n_points, 1.0d0, forces_shell, n_compute )
      ! Work out the mGGA contribution to the matrix by calculating
      ! gradient_basis_wave.(hessian_basis_wave*xc_tau_deriv)
      ! This'll go into the hamiltonian. Forces match KS implementation! AJL
      if (meta_gga_forces_on) then
        call dgemm('N', 'N', n_compute, n_compute,  &
                   n_points, 1.d0,  &
                   matrix_tem1_mgga, n_compute, matrix_tem2_mgga, &
                   n_points, 1.0d0, forces_shell, n_compute )
      end if
    end do
  end subroutine evaluate_gga_forces_dens_mat
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/update_sum_forces
!  NAME
!    update_sum_forces
!  SYNOPSIS
  subroutine update_sum_forces &
       ( n_compute_c, n_compute_a, i_basis, ins_idx, n_compute_for_atom, &
         forces_shell, dens_mat, ld_dens_mat, i_spin, i_coord, sum_forces  )
!  PURPOSE
!    The subroutine is a wrapper around the CPU and GPU versions of
!    update_sum_forces.
!  USES
    use dimensions, only: n_atoms
    use runtime_choices, only: use_gpu
    use load_balancing, only: use_batch_permutation
    implicit none
!  ARGUMENTS
    integer, intent(in)    :: n_compute_c
    integer, intent(in)    :: n_compute_a
    integer, intent(in)    :: i_basis(n_compute_c)
    integer, intent(in)    :: ins_idx(n_compute_c)
    integer, intent(in)    :: n_compute_for_atom(n_atoms)
    real*8,  intent(inout) :: forces_shell(n_compute_c, n_compute_a)
    real*8,  intent(in)    :: dens_mat(*)
    integer, intent(in)    :: ld_dens_mat
    integer, intent(in)    :: i_spin
    integer, intent(in)    :: i_coord
    real*8,  intent(inout) :: sum_forces(3, n_atoms)
!  INPUTS
!    o n_compute_c -- the number of basis functins in the grid points
!    o n_compute_a -- the number of basis functins in the grid points
!    o i_basis -- the lists basis functions relevant in the grid patch
!    o forces_shell -- results of the integrals
!    o dens_mat -- the density matrix, possibly energy-weighted, for the current
!                  spin channel
!    o offset_densmat -- offset for spin channel of density matrix
!    o offset_densmat -- Current spin channel
!    o i_coord -- dimension of the forces calculated here
!  OUTPUT
!    o sum_forces -- the sum of Pulay, gga, meta-GGA, and gnonmf forces
!  AUTHOR
!    William Huhn (Duke University)
!  HISTORY
!    February 2018 - Created.
!  SOURCE
    integer :: i_off ! offset to account for spin channel indexing, as dens_mat
                     ! is passed around as a 1D assumed size array, whereas it
                     ! is actually a 2D array

    i_off = (i_spin-1)*ld_dens_mat

    if (use_gpu) then
      if (use_batch_permutation.gt.0) then
        ! Do indexing on the GPU
        call update_sum_forces_gpu &
             ( n_compute_c, i_coord, i_spin, n_compute_for_atom )
      else
        ! Indexing on GPU is too expensive in memory; move results
        ! back to CPU and do indexing there
        call get_forces_shell_gpu(forces_shell)
        call update_sum_forces_cpu &
             ( n_compute_c, n_compute_a, i_basis, ins_idx, &
             forces_shell, dens_mat(1+i_off), &
             i_coord, sum_forces )
      end if
    else
      call update_sum_forces_cpu &
           ( n_compute_c, n_compute_a, i_basis, ins_idx, &
           forces_shell, dens_mat(1+i_off), &
           i_coord, sum_forces )
    end if
  end subroutine update_sum_forces
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/AS_update_sum_forces_and_stress
!  NAME
!    AS_update_sum_forces_and_stress
!  SYNOPSIS
  subroutine AS_update_sum_forces_and_stress &
       ( n_compute_c, n_compute_a, i_basis, ins_idx, n_compute_for_atom, &
         forces_shell, dens_mat, ld_dens_mat, i_spin, i_coord, sum_forces, &
         AS_strain_deriv_wave_shell, relativistic, AS_corr_shell )
!  PURPOSE
!    The subroutine is a wrapper around the CPU and GPU versions of
!    AS_update_sum_forces_and_stress.
!  USES
    use analytical_stress, only: AS_pulay_stress_local
    use dimensions, only: n_atoms
    use runtime_choices, only: use_gpu
    use load_balancing, only: use_batch_permutation
    implicit none
!  ARGUMENTS
    integer, intent(in)    :: n_compute_c
    integer, intent(in)    :: n_compute_a
    integer, intent(in)    :: i_basis(1:n_compute_c)
    integer, intent(in)    :: ins_idx(n_compute_c)
    integer, intent(in)    :: n_compute_for_atom(n_atoms)
    real*8,  intent(inout) :: forces_shell(1:n_compute_c,1:n_compute_a)
    real*8,  intent(in)    :: dens_mat(*)
    integer, intent(in)    :: ld_dens_mat
    integer, intent(in)    :: i_spin
    integer, intent(in)    :: i_coord
    real*8,  intent(inout) :: sum_forces(1:3,1:n_atoms)
    real*8,  intent(in)    :: AS_strain_deriv_wave_shell(1:n_compute_c,1:n_compute_a)
    logical, intent(in)    :: relativistic
    real*8,  intent(in), optional :: AS_corr_shell(1:n_compute_c,1:n_compute_a)
!  INPUTS
!    o n_compute_c -- the number of basis functins in the grid points
!    o n_compute_a -- the number of basis functins in the grid points
!    o i_basis -- the lists basis functions relevant in the grid patch
!    o forces_shell -- results of the integrals
!    o dens_mat -- the density matrix, possibly energy-weighted, for the current
!                  spin channel
!    o offset_densmat -- offset for spin channel of density matrix
!    o offset_densmat -- Current spin channel
!    o i_coord -- dimension of the forces calculated here
!  OUTPUT
!    o sum_forces -- the sum of Pulay, gga, meta-GGA, and gnonmf forces
!  AUTHOR
!    William Huhn (Duke University)
!  HISTORY
!    February 2018 - Created.
!  SOURCE
    integer :: i_off ! offset to account for spin channel indexing, as dens_mat
                     ! is passed around as a 1D assumed size array, whereas it
                     ! is actually a 2D array

    i_off = (i_spin-1)*ld_dens_mat

    if (use_gpu) then
      if (use_batch_permutation.gt.0) then
        ! Load balancing works beautifully on GPUs, since the problem is locally
        ! dense.  So do indexing on the GPU for sum_forces and
        ! AS_pulay_stress_local.
        call AS_update_sum_forces_and_stress_gpu &
             ( n_compute_c, i_coord, i_spin, n_compute_for_atom )
        ! We construct AS_corr_shell on the CPU, because there's too little
        ! compute to make GPU acceleration worth the development time.  So we do
        ! the indexing for the non-relativistic correction on the CPU.
        if (.not. relativistic .and. present(AS_corr_shell)) then
          call AS_update_nonrel_term &
             ( n_compute_c, n_compute_a, i_basis, ins_idx, dens_mat(1+i_off), &
               i_coord, AS_corr_shell )
        end if
      else
        ! WPH: For all non-load-balanced cases, we communicate results back to
        !      the CPU and do indexing there.  The rationale is as follows:
        !      1) Non-packed matrices:  Here we can index on the GPU, because
        !      the indexing is identical to the load-balanced case.  Which we do
        !      for Hamiltonian integration.  But non-packed matrices are only
        !      supported for orbital-based density updates.
        !      2) CSR-packed matrices:  The CSR format in aims is not
        !      distributed in memory; every MPI task has a full copy of the
        !      matrix.  Add to this the fact that this matrix format has several
        !      indexing array which make vectorized approaches difficult, and
        !      this just isn't a good fit for GPUs.  I could probably get around
        !      the indexing array issue if I really tried (cuSparse does support
        !      CSR, so there is a way to vectorize this on a GPU), but the
        !      memory consumption is a critical algorithmic flaw that can't be
        !      worked around.
        !      3) Domain decomposition without load balancing:  While this code
        !      path circumvents the memory issue, it still runs into the
        !      indexing array issue inherent in sparse matrix formats.  In my
        !      opinion, if one needs GPU acceleration for an algorithm that does
        !      not support load balancing, it's better to load balance the
        !      algorithm instead and exploit the embarassing parallelism
        !      possible with dense matrix formats, rather than shoehorning in
        !      sparse matrix formats.
        call get_forces_shell_gpu(forces_shell)
        call get_AS_strain_deriv_wave_shell_gpu(AS_strain_deriv_wave_shell)
        if (present(AS_corr_shell)) then
          call AS_update_sum_forces_and_stress_cpu &
               ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                 forces_shell, dens_mat(1+i_off), &
                 i_coord, sum_forces, &
                 AS_strain_deriv_wave_shell, &
                 relativistic, AS_corr_shell )
        else
          call AS_update_sum_forces_and_stress_cpu &
               ( n_compute_c, n_compute_a, i_basis, ins_idx, &
                 forces_shell, dens_mat(1+i_off), &
                 i_coord, sum_forces, &
                 AS_strain_deriv_wave_shell, &
                 relativistic )
        end if
      end if
    else
      if (present(AS_corr_shell)) then
        call AS_update_sum_forces_and_stress_cpu &
             ( n_compute_c, n_compute_a, i_basis, ins_idx, &
               forces_shell, dens_mat(1+i_off), &
               i_coord, sum_forces, &
               AS_strain_deriv_wave_shell, &
               relativistic, AS_corr_shell )
      else
        call AS_update_sum_forces_and_stress_cpu &
             ( n_compute_c, n_compute_a, i_basis, ins_idx, &
               forces_shell, dens_mat(1+i_off), &
               i_coord, sum_forces, &
               AS_strain_deriv_wave_shell, &
               relativistic )
      end if
    end if
  end subroutine AS_update_sum_forces_and_stress
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/update_sum_forces_cpu
!  NAME
!    update_sum_forces_cpu
!  SYNOPSIS
  subroutine update_sum_forces_cpu &
       ( n_compute_c, n_compute_a, i_basis, ins_idx, forces_shell, dens_mat, &
         i_coord, sum_forces  )
!  PURPOSE
!    The subroutine calculates the Pulay and gga force components for one grid
!    patch using the density matrix. Here the actual force components are
!    evaluated from the integral results for one grid patch calculated
!    previously to forces_shell.
!
!    The link between i_compute = 1 ... n_compute
!    and i_basis = 1 ... n_basis is provided by the index array
!    i_basis(i_compute).
!  USES
    use pbc_lists
    use runtime_choices
    use dimensions, only: n_atoms, n_periodic
    use load_balancing, only: batch_perm, use_batch_permutation
    implicit none
!  ARGUMENTS
    integer, intent(in)    :: n_compute_c
    integer, intent(in)    :: n_compute_a
    integer, intent(in)    :: i_basis(n_compute_c)
    integer, intent(in)    :: ins_idx(n_compute_c)
    real*8,  intent(in)    :: forces_shell(n_compute_c, n_compute_a)
    real*8,  intent(in)    :: dens_mat(*)
    integer, intent(in)    :: i_coord
    real*8,  intent(inout) :: sum_forces(3, n_atoms)
!  INPUTS
!    o n_compute_c -- the number of basis functins in the grid points
!    o n_compute_a -- the number of basis functins in the grid points
!    o i_basis -- the lists basis functions relevant in the grid patch
!    o forces_shell -- results of the integrals
!    o dens_mat -- the density matrix, possibly energy-weighted, for the current
!                  spin channel
!    o i_coord -- dimension of the forces calculated here
!  OUTPUT
!    o sum_forces -- the sum of Pulay, gga, meta-GGA, and gnonmf forces
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    integer :: i_compute_1
    integer :: i_compute_2
    integer :: i_offset
    integer :: i_index_real
    integer :: i_offset_first_part
    integer :: i_cell_index
    integer :: i_max_basis, i_min_basis
    integer :: i_one_part
    integer :: i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell, i_cell_1
    integer :: offset(n_cells)
    integer :: offset_end(n_cells)
    integer :: help
    integer :: i_cell_old

    ! Variables specific to load balancing
    integer :: n_bp
    integer :: i_basis_loc_1, i_basis_loc_2

    ! now add the aux. ham. matrix to the actual ham. matrix
    ! this requires translating between the actually computed matrix elements
    ! and the full ham. matrix ...

    ! When basis is smaller than n_basis

    if(use_batch_permutation > 0) then
      ! If use_batch_permutation > 0 is set, the local hamiltonian is
      ! always stored in full form for the local basis functions
      ! Get position of basis functions of current batch within local
      ! matrix
      !
      ! WPH: There are three different sets of basis element indices that must
      ! be appreciated to understand (and extend) this code:
      ! 1) The set of basis elements contributing to the current batch:
      !       i_compute = 1,...,n_compute_c
      !    This is the set of basis elements with non-zero support on the
      !    current batch.  This set will change from batch to batch and varies
      !    from MPI task to MPI tasks  These are the basis elements making up
      !    the current "shell".  This set of basis elements occurs in batch
      !    integration routines and will occur in all indexing scheme.
      ! 2) The set of all basis elements for the calculation:
      !       i_basis = 1,...,n_basis
      !    This is the usual set of basis elements in aims, which is fixed
      !    across all batches and all MPI tasks.
      ! 3) The set of basis elements contributing to any batch on the current
      !    MPI task:
      !       i_basis_loc = 1,...,batch_perm(n_bp)%n_basis_local
      !    This is the set of basis elements with non-zero support on at least
      !    one batch from the current MPI task.  This set varies from MPI task
      !    to MPI task, but is fixed across all batches for a given MPI task.
      !    These are the basis elements making up matrices which have a
      !    real-space nature (Hamiltonian, density matrix, etc.)  This basis set
      !    is unique to the load-balancing indexing scheme.
      ! To convert from basis set (1) to basis set (2), one uses the i_basis
      ! array, which is generated as part of the prune_basis_* subroutine.  To
      ! convert from basis set (2) to basis set (3), one uses the
      ! batch_perm(n_bp)%i_basis_glb_to_loc array, where n_bp is the current
      ! load balancing distribution, which is generated during the load
      ! balancing initialization.  Converting from basis set (1) to (3) is done
      ! by ins_idx.
      !
      ! As long as one clearly distinguishes the three basis sets, the indexing
      ! of the load-balanced domain decomposition should be functionally
      ! identical to the non-packed case, as they both store the data in a
      ! full/dense format.  I've kept the variable names consistent to reflect
      ! this.

      do i_compute_2 = 1, n_compute_c
        i_basis_2 = i_basis(i_compute_2)
        i_basis_loc_2 = ins_idx(i_compute_2)
        i_offset = (i_basis_loc_2-1)*i_basis_loc_2/2

        do i_compute_1 = 1, i_compute_2 ! n_compute_c
          i_basis_1 = i_basis(i_compute_1)
          i_basis_loc_1 = ins_idx(i_compute_1)
          i_index_real = i_offset + i_basis_loc_1

          sum_forces(i_coord, Cbasis_to_atom( i_basis_1 )) = &
               sum_forces(i_coord, Cbasis_to_atom( i_basis_1 )) + &
               dens_mat(i_index_real) * forces_shell(i_compute_2,i_compute_1)

          if (i_basis_loc_1/= i_basis_loc_2) then
            sum_forces(i_coord, Cbasis_to_atom( i_basis_2 )) = &
                 sum_forces(i_coord, Cbasis_to_atom( i_basis_2 )) + &
                 dens_mat(i_index_real) * forces_shell(i_compute_1,i_compute_2)
          end if
        end do ! i_compute_1
      end do ! i_compute_2

      ! Finished indexing, exit
      return
    end if ! use_batch_permutation > 0

    select case(packed_matrix_format)
    case (PM_none)
      i_index_real = 0
      do i_compute_2 = 1, n_compute_c, 1
        i_basis_2 = i_basis(i_compute_2)
        i_offset = (i_basis_2-1)*i_basis_2/2

        do i_compute_1 = 1,i_compute_2 ,1
          i_basis_1 = i_basis(i_compute_1)
          i_index_real = i_offset + i_basis_1

          sum_forces(i_coord, Cbasis_to_atom( i_basis_1 )) = &
               sum_forces(i_coord, Cbasis_to_atom( i_basis_1)) + &
               dens_mat(i_index_real) * forces_shell(i_compute_2, i_compute_1)

          if (i_basis_1/= i_basis_2) then
            sum_forces(i_coord, Cbasis_to_atom( i_basis_2)) = &
                 sum_forces(i_coord, Cbasis_to_atom( i_basis_2)) + &
                 dens_mat(i_index_real) * forces_shell(i_compute_1, i_compute_2)
          end if
        enddo
      enddo
    case(PM_index)
      if(n_periodic == 0)then
        do i_compute_1 = 1, n_compute_a, 1
          i_basis_1 = i_basis(i_compute_1)
          i_start =  index_hamiltonian(1, 1, i_basis_1)
          i_end   =  index_hamiltonian(2, 1, i_basis_1)

          do i_compute_2 = 1, n_compute_a, 1
            i_basis_2 = i_basis(i_compute_2)

            place: do i_place = i_start, i_end, 1
              if( column_index_hamiltonian( i_place) == i_basis_2)then
                sum_forces(i_coord, Cbasis_to_atom( i_basis_1)) = &
                     sum_forces(i_coord, Cbasis_to_atom( i_basis_1)) + &
                     dens_mat(i_place) * forces_shell(i_compute_2, i_compute_1)
                if(i_basis_1/= i_basis_2)then
                  sum_forces(i_coord, Cbasis_to_atom( i_basis_2)) = &
                       sum_forces(i_coord, Cbasis_to_atom( i_basis_2)) + &
                       dens_mat(i_place) * forces_shell(i_compute_1, i_compute_2)
                end if
                i_index_real = i_place

                exit place
              else if(column_index_hamiltonian( i_place) > i_basis_2)then
                i_index_real = i_place

                exit place
              end if
            end do place

            i_start = i_index_real
          end do ! i_compute_2
        end do ! i_compute_1
      else ! Periodic case----------------
        ! Unfortunately the periodic systems can not use the searching routine
        ! used now in the clusters.  This is because the peridic systems have
        ! extra packing for supercell information.
        do i_compute_1 = 1, n_compute_a, 1
          i_basis_1 = Cbasis_to_basis(i_basis(i_compute_1))
          i_cell_old = center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))

          offset_end = -1
          offset = -1
          do i_cell_1 = 1, n_cells
            i_cell = position_in_hamiltonian( i_cell_old, i_cell_1)
            offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
            offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)
          end do

          do i_compute_2 = 1, n_compute_a, 1
            i_basis_2 = Cbasis_to_basis(i_basis(i_compute_2))
            if(i_basis_2 <= i_basis_1)then
              i_cell = center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))

              place_2: do i_place = offset(i_cell), offset_end(i_cell), 1
                if( column_index_hamiltonian(i_place) == i_basis_2)then
                  sum_forces(i_coord, Cbasis_to_atom( i_basis_1)) = &
                       sum_forces(i_coord, Cbasis_to_atom( i_basis_1)) + &
                       dens_mat(i_place) * forces_shell(i_compute_2, i_compute_1)

                  if(i_basis_1/= i_basis_2)then
                    sum_forces(i_coord, Cbasis_to_atom( i_basis_2)) = &
                         sum_forces(i_coord, Cbasis_to_atom( i_basis_2)) + &
                         dens_mat(i_place)*forces_shell(i_compute_1, i_compute_2)
                  end if

                  exit place_2
                else if(column_index_hamiltonian( i_place) > i_basis_2)then
                  exit place_2
                end if
              end do place_2

              offset(i_cell) = i_place
            end if
          end do ! i_compute_2
        end do ! i_compute_1
      end if ! n_periodic == 0
    end select ! packed_matrix_format
  end subroutine update_sum_forces_cpu
!******

!-------------------------------------------------------------------------------
!****s* forces_densmat/AS_update_nonrel_term
!  NAME
!    AS_update_nonrel_term
!  SYNOPSIS
  subroutine AS_update_nonrel_term &
       ( n_compute_c, n_compute_a, i_basis, ins_idx, dens_mat, &
         i_coord, AS_corr_shell )
!  PURPOSE
!    The subroutine calculates only the correction term for non-relativistic
!    calculations.
!
!    This subroutine is a reduced version of the subroutine
!    AS_update_sum_forces_and_stress_cpu.  It is used when performing GPU
!    acceleration with load balancing, because in this case the CPU continues
!    to calculate the non-relativistic correction, whereas the GPU calculates
!    both the forces and the "Pulay" terms.
!
!    The link between i_compute = 1 ... n_compute and i_basis = 1 ... n_basis is
!    provided by the index array i_basis(i_compute).
!  USES
    use pbc_lists
    use runtime_choices
    use dimensions, only: n_atoms
    use analytical_stress, only: AS_dde_T_sa_at_di_orb_local
    use load_balancing, only: batch_perm, use_batch_permutation
    implicit none
!  ARGUMENTS
    integer, intent(in)  :: n_compute_c
    integer, intent(in)  :: n_compute_a
    integer, intent(in)  :: i_basis(1:n_compute_c)
    integer, intent(in)  :: ins_idx(n_compute_c)
    real*8,  intent(in)  :: dens_mat(*)
    integer, intent(in)  :: i_coord
    real*8,  intent(in)  :: AS_corr_shell(1:n_compute_c,1:n_compute_a)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    February 2018 - Forked off of AS_update_sum_forces_and_stress_cpu.
!  SOURCE
    integer :: i_compute_1
    integer :: i_compute_2
    integer :: i_offset
    integer :: i_index_real
    integer :: i_one_part

    integer :: i_start, i_end, i_place, i_basis_1, i_basis_2, i_cell, i_cell_1
    integer :: offset(n_cells)
    integer :: offset_end(n_cells)
    integer :: i_cell_old
    integer :: i_coords_1
    integer :: i_coords_2
    integer :: i_counter

    ! Variables specific to load balancing
    integer :: i_basis_loc_1, i_basis_loc_2

    ! now add the aux. ham. matrix to the actual ham. matrix
    ! this requires translating between the actually computed matrix elements
    ! and the full ham. matrix ...

    ! When basis is smaller than n_basis

    if(use_batch_permutation > 0) then

      do i_compute_2 = 1, n_compute_c
        i_basis_2 = i_basis(i_compute_2)
        i_basis_loc_2 = ins_idx(i_compute_2)
        i_offset = (i_basis_loc_2-1)*i_basis_loc_2/2

        do i_compute_1 = 1, i_compute_2 ! n_compute_c
          i_basis_1 = i_basis(i_compute_1)
          i_basis_loc_1 = ins_idx(i_compute_1)
          i_index_real = i_offset + i_basis_loc_1

          if ( Cbasis_to_center( i_basis_1 ) &
               .eq. &
               Cbasis_to_center( i_basis_2 )) &
          then ! on-site
            AS_dde_T_sa_at_di_orb_local(i_coord) = &
                 AS_dde_T_sa_at_di_orb_local(i_coord) &
                 + dens_mat(i_index_real) &
                 * AS_corr_shell(i_compute_2, i_compute_1)
          end if

          if (i_basis_loc_1 /= i_basis_loc_2) then
            if ( Cbasis_to_center( i_basis_1 ) &
                .eq. &
                Cbasis_to_center( i_basis_2) ) &
            then ! on-site
              AS_dde_T_sa_at_di_orb_local(i_coord) = &
                   AS_dde_T_sa_at_di_orb_local(i_coord) &
                   + dens_mat(i_index_real) &
                   * AS_corr_shell(i_compute_1, i_compute_2)
            end if
          end if
        end do ! i_compute_1
      end do ! i_compute_2

      ! Finished indexing, exit
      return
    end if ! use_batch_permutation > 0

    ! CC: Just periodic case is of interest here:
    if (packed_matrix_format==PM_none) then
      i_index_real = 0
      do i_compute_2 = 1, n_compute_c, 1
        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2
        i_basis_2 = i_basis(i_compute_2)

        do i_compute_1 = 1,i_compute_2 ,1
          i_index_real = i_offset + i_basis(i_compute_1)
          i_basis_1 = i_basis(i_compute_1)

          if ( Cbasis_to_center(i_basis(i_compute_1)) &
               .eq. &
               Cbasis_to_center(i_basis(i_compute_2)) ) &
          then ! on-site
            AS_dde_T_sa_at_di_orb_local(i_coord) = &
                 AS_dde_T_sa_at_di_orb_local(i_coord) &
                 + dens_mat(i_index_real) &
                 * AS_corr_shell(i_compute_2,i_compute_1)
          end if

          if(i_basis_1/= i_basis_2)then
            if ( Cbasis_to_center(i_basis(i_compute_1)) &
                .eq. &
                Cbasis_to_center(i_basis(i_compute_2)) ) &
            then ! on-site
              AS_dde_T_sa_at_di_orb_local(i_coord) = &
                   AS_dde_T_sa_at_di_orb_local(i_coord) &
                   + dens_mat(i_index_real) &
                   * AS_corr_shell(i_compute_1,i_compute_2)
            end if
          end if
        end do
      end do
    else !packed_matrix_format
      select case(packed_matrix_format)
      ! Unfortunately the periodic systems can not use the searching routine
      ! used now in the clusters.  This is because the peridic systems have
      ! extra packing for supercell information.
      case(PM_index)
        do i_compute_1 = 1, n_compute_a, 1
          i_basis_1 = Cbasis_to_basis(i_basis(i_compute_1))
          i_cell_old = center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))

          offset_end = -1
          offset = -1

          do i_cell_1 = 1, n_cells
            i_cell = position_in_hamiltonian( i_cell_old, i_cell_1)
            offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
            offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)
          end do

          do i_compute_2 = 1, n_compute_a
            i_basis_2 = Cbasis_to_basis(i_basis(i_compute_2))

            if(i_basis_2 <= i_basis_1)then
              i_cell = center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))

              place_2_ed: do i_place = offset(i_cell), offset_end(i_cell), 1
                if( column_index_hamiltonian( i_place) == i_basis_2)then
                  if ( Cbasis_to_center(i_basis(i_compute_1)) &
                       .eq. &
                       Cbasis_to_center(i_basis(i_compute_2)) ) &
                  then ! on-site
                    AS_dde_T_sa_at_di_orb_local(i_coord) = &
                         AS_dde_T_sa_at_di_orb_local(i_coord) &
                         + dens_mat(i_place) &
                         * AS_corr_shell(i_compute_2,i_compute_1)
                  end if

                  if(i_basis_1/= i_basis_2)then
                    if ( Cbasis_to_center(i_basis(i_compute_1)) &
                         .eq. &
                         Cbasis_to_center(i_basis(i_compute_2)) )&
                    then ! on-site
                      AS_dde_T_sa_at_di_orb_local(i_coord) = &
                           AS_dde_T_sa_at_di_orb_local(i_coord) + &
                           ( dens_mat(i_place) * &
                           AS_corr_shell(i_compute_1,i_compute_2) )
                    end if
                  end if
                  exit place_2_ed
                else if(column_index_hamiltonian( i_place) > i_basis_2)then
                   exit place_2_ed
                end if
              end do place_2_ed

              offset(i_cell) = i_place
            end if ! i_basis_2 <= i_basis_1
          end do ! i_compute_2
        end do ! i_compute_1
      end select ! PM_index
    end if ! packed_matrix_format
  end subroutine AS_update_nonrel_term

!-------------------------------------------------------------------------------
!****s* forces_densmat/AS_update_sum_forces_and_stress_cpu
!  NAME
!    AS_update_sum_forces_and_stress_cpu
!  SYNOPSIS
  subroutine AS_update_sum_forces_and_stress_cpu &
       ( n_compute_c, n_compute_a, i_basis, ins_idx, forces_shell, dens_mat, &
         i_coord, sum_forces, AS_strain_deriv_wave_shell, relativistic, &
         AS_corr_shell )
!  PURPOSE
!    CC: Same as update_sum_forces, but also adds up stress
!    The subroutine calculates the Pulay and gga force components for the energy
!    density/analytical stress case for one grid patch using the density matrix.
!    Here the actual force components are evaluated from the integral results
!    for one grid patch calculated previously to forces_shell.
!
!    This subroutine is a modified version of the subroutine update_sum_forces
!    above.
!
!    The link between i_compute = 1 ... n_compute and i_basis = 1 ... n_basis is
!    provided by the index array i_basis(i_compute).
!  USES
    use pbc_lists
    use runtime_choices
    use dimensions, only: n_atoms, compute_heat_flux
    use analytical_stress, only: AS_dde_T_sa_at_di_orb_local, &
        AS_pulay_stress_local
    use heat_flux
    use load_balancing, only: batch_perm, use_batch_permutation
    implicit none
!  ARGUMENTS
    integer, intent(in)  :: n_compute_c
    integer, intent(in)  :: n_compute_a
    integer, intent(in)  :: i_basis(1:n_compute_c)
    integer, intent(in)  :: ins_idx(n_compute_c)
    real*8,  intent(in)  :: forces_shell(1:n_compute_c,1:n_compute_a)
    real*8,  intent(in)  :: dens_mat(*)
    integer, intent(in)  :: i_coord
    real*8,  intent(out) :: sum_forces(1:3,1:n_atoms)
    real*8,  intent(in)  :: AS_strain_deriv_wave_shell(1:n_compute_c,1:n_compute_a)
    logical, intent(in)  :: relativistic
    real*8,  intent(in), optional :: AS_corr_shell(1:n_compute_c,1:n_compute_a)
!  INPUTS
!    o n_compute_c -- the number of basis functins in the grid points
!    o n_compute_a -- the number of basis functins in the grid points
!    o i_basis -- the lists basis functions relevant in the grid patch
!    o forces_shell -- results of the integrals
!    o dens_mat -- the density matrix, possibly energy-weighted, for the current
!                  spin channel
!    o i_coord -- dimension of the forces calculated here
!  OUTPUT
!    o sum_forces -- the sum of Pulay, gga, meta-GGA, and gnonmf forces
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
    integer :: i_compute_1
    integer :: i_compute_2
    integer :: i_offset
    integer :: i_index_real
    integer :: i_one_part

    integer :: i_start, i_end, i_place, i_basis_1, i_basis_2, i_cell, i_cell_1
    integer :: offset(n_cells)
    integer :: offset_end(n_cells)
    integer :: i_cell_old
    integer :: i_coords_1
    integer :: i_coords_2
    integer :: i_counter
    integer :: AS_atom

    ! Variables specific to load balancing
    integer :: i_basis_loc_1, i_basis_loc_2

    ! now add the aux. ham. matrix to the actual ham. matrix
    ! this requires translating between the actually computed matrix elements
    ! and the full ham. matrix ...

    ! When basis is smaller than n_basis

    if(use_batch_permutation > 0) then
      ! If use_batch_permutation > 0 is set, the local hamiltonian is
      ! always stored in full form for the local basis functions
      ! Get position of basis functions of current batch within local
      ! matrix
      !
      ! WPH: There are three different sets of basis element indices that must
      ! be appreciated to understand (and extend) this code:
      ! 1) The set of basis elements contributing to the current batch:
      !       i_compute = 1,...,n_compute_c
      !    This is the set of basis elements with non-zero support on the
      !    current batch.  This set will change from batch to batch and varies
      !    from MPI task to MPI tasks  These are the basis elements making up
      !    the current "shell".  This set of basis elements occurs in batch
      !    integration routines and will occur in all indexing scheme.
      ! 2) The set of all basis elements for the calculation:
      !       i_basis = 1,...,n_basis
      !    This is the usual set of basis elements in aims, which is fixed
      !    across all batches and all MPI tasks.
      ! 3) The set of basis elements contributing to any batch on the current
      !    MPI task:
      !       i_basis_loc = 1,...,batch_perm(n_bp)%n_basis_local
      !    This is the set of basis elements with non-zero support on at least
      !    one batch from the current MPI task.  This set varies from MPI task
      !    to MPI task, but is fixed across all batches for a given MPI task.
      !    These are the basis elements making up matrices which have a
      !    real-space nature (Hamiltonian, density matrix, etc.)  This basis set
      !    is unique to the load-balancing indexing scheme.
      ! To convert from basis set (1) to basis set (2), one uses the i_basis
      ! array, which is generated as part of the prune_basis_* subroutine.  To
      ! convert from basis set (2) to basis set (3), one uses the
      ! batch_perm(n_bp)%i_basis_glb_to_loc array, where n_bp is the current
      ! load balancing distribution, which is generated during the load
      ! balancing initialization.  Converting from basis set (1) to (3) is done
      ! by ins_idx, which is constructed here.
      !
      ! As long as one clearly distinguishes the three basis sets, the indexing
      ! of the load-balanced domain decomposition should be functionally
      ! identical to the non-packed case, as they both store the data in a
      ! full/dense format.  I've kept the variable names consistent to reflect
      ! this.

      do i_compute_2 = 1, n_compute_c
        i_basis_2 = i_basis(i_compute_2)
        i_basis_loc_2 = ins_idx(i_compute_2)
        i_offset = (i_basis_loc_2-1)*i_basis_loc_2/2

        do i_compute_1 = 1, i_compute_2 ! n_compute_c
          i_basis_1 = i_basis(i_compute_1)
          i_basis_loc_1 = ins_idx(i_compute_1)
          i_index_real = i_offset + i_basis_loc_1

          if (i_coord.le.3) then
            sum_forces(i_coord, Cbasis_to_atom( i_basis_1 )) = &
                 sum_forces(i_coord, Cbasis_to_atom( i_basis_1 )) + &
                 dens_mat(i_index_real) * forces_shell(i_compute_2,i_compute_1)
          end if

          if (.not.relativistic) then
            if ( Cbasis_to_center( i_basis_1 ) &
                 .eq. &
                 Cbasis_to_center( i_basis_2 )) &
            then ! on-site
              AS_dde_T_sa_at_di_orb_local(i_coord) = &
                   AS_dde_T_sa_at_di_orb_local(i_coord) &
                   + dens_mat(i_index_real) &
                   * AS_corr_shell(i_compute_2, i_compute_1)
            end if
          end if
          AS_pulay_stress_local(i_coord) =  &
               AS_pulay_stress_local(i_coord) &
               + dens_mat(i_index_real) &
               * AS_strain_deriv_wave_shell(i_compute_2, i_compute_1)

          if (i_basis_loc_1 /= i_basis_loc_2) then
            if (i_coord.le.3) then
              sum_forces(i_coord, Cbasis_to_atom( i_basis_2 )) = &
                   sum_forces(i_coord, Cbasis_to_atom( i_basis_2 )) + &
                   dens_mat(i_index_real) * forces_shell(i_compute_1,i_compute_2)
            end if

            if (.not.relativistic) then
              if ( Cbasis_to_center( i_basis_1 ) &
                  .eq. &
                  Cbasis_to_center( i_basis_2) ) &
              then ! on-site
                AS_dde_T_sa_at_di_orb_local(i_coord) = &
                     AS_dde_T_sa_at_di_orb_local(i_coord) &
                     + dens_mat(i_index_real) &
                     * AS_corr_shell(i_compute_1, i_compute_2)
               end if
            end if
            AS_pulay_stress_local(i_coord) = &
                 AS_pulay_stress_local(i_coord) &
                 + dens_mat(i_index_real) &
                 * AS_strain_deriv_wave_shell(i_compute_1, i_compute_2)
          end if
        end do ! i_compute_1
      end do ! i_compute_2

      ! Finished indexing, exit
      return
    end if ! use_batch_permutation > 0

    ! CC: Just periodic case is of interest here:
    if (packed_matrix_format==PM_none) then
      i_index_real = 0
      do i_compute_2 = 1, n_compute_c, 1
        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2
        i_basis_2 = i_basis(i_compute_2)

        do i_compute_1 = 1,i_compute_2 ,1
          i_index_real = i_offset + i_basis(i_compute_1)
          i_basis_1 = i_basis(i_compute_1)

          if (i_coord.le.3) then
            sum_forces(i_coord, Cbasis_to_atom(i_basis_1 )) = &
                 sum_forces(i_coord, Cbasis_to_atom(i_basis_1)) &
                 + dens_mat(i_index_real) &
                 * forces_shell(i_compute_2, i_compute_1)
          end if

          if (.not.relativistic) then
            if ( Cbasis_to_center(i_basis(i_compute_1)) &
                 .eq. &
                 Cbasis_to_center(i_basis(i_compute_2)) ) &
            then ! on-site
              AS_dde_T_sa_at_di_orb_local(i_coord) = &
                   AS_dde_T_sa_at_di_orb_local(i_coord) &
                   + dens_mat(i_index_real) &
                   * AS_corr_shell(i_compute_2,i_compute_1)
            end if
          end if
          AS_pulay_stress_local(i_coord) =  &
               AS_pulay_stress_local(i_coord) &
               + dens_mat(i_index_real) &
               * AS_strain_deriv_wave_shell( i_compute_2,i_compute_1 )

          if(i_basis_1/= i_basis_2)then
            if (i_coord.le.3) then
              sum_forces(i_coord, Cbasis_to_atom(i_basis_2)) = &
                   sum_forces(i_coord, Cbasis_to_atom(i_basis_2)) &
                    + dens_mat(i_index_real) &
                    * forces_shell(i_compute_1, i_compute_2)
            end if

            if (.not.relativistic) then
              if ( Cbasis_to_center(i_basis(i_compute_1)) &
                  .eq. &
                  Cbasis_to_center(i_basis(i_compute_2)) ) &
              then ! on-site
                AS_dde_T_sa_at_di_orb_local(i_coord) = &
                     AS_dde_T_sa_at_di_orb_local(i_coord) &
                     + dens_mat(i_index_real) &
                     * AS_corr_shell(i_compute_1,i_compute_2)
              end if
            end if
            AS_pulay_stress_local(i_coord) = &
                 AS_pulay_stress_local(i_coord) &
                 + dens_mat(i_index_real) &
                 * AS_strain_deriv_wave_shell( i_compute_1,i_compute_2 )
          end if
        end do
      end do
    else !packed_matrix_format
      select case(packed_matrix_format)
      ! Unfortunately the periodic systems can not use the searching routine
      ! used now in the clusters.  This is because the peridic systems have
      ! extra packing for supercell information.
      case(PM_index)
        do i_compute_1 = 1, n_compute_a, 1
          i_basis_1 = Cbasis_to_basis(i_basis(i_compute_1))
          i_cell_old = center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))

          offset_end = -1
          offset = -1

          do i_cell_1 = 1, n_cells
            i_cell = position_in_hamiltonian( i_cell_old, i_cell_1)
            offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
            offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)
          end do

          do i_compute_2 = 1, n_compute_a
            i_basis_2 = Cbasis_to_basis(i_basis(i_compute_2))

            if(i_basis_2 <= i_basis_1)then
              i_cell = center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))

              place_2_ed: do i_place = offset(i_cell), offset_end(i_cell), 1
                if( column_index_hamiltonian( i_place) == i_basis_2)then
                  if (i_coord.le.3) then
                    sum_forces(i_coord, Cbasis_to_atom(i_basis_1)) = &
                         sum_forces(i_coord,Cbasis_to_atom(i_basis_1)) &
                         + dens_mat(i_place) &
                         * forces_shell(i_compute_2, i_compute_1)
                  end if

                  if (.not.relativistic) then
                    if ( Cbasis_to_center(i_basis(i_compute_1)) &
                         .eq. &
                         Cbasis_to_center(i_basis(i_compute_2)) ) &
                    then ! on-site
                      AS_dde_T_sa_at_di_orb_local(i_coord) = &
                           AS_dde_T_sa_at_di_orb_local(i_coord) &
                           + dens_mat(i_place) &
                           * AS_corr_shell(i_compute_2,i_compute_1)
                      if (compute_heat_flux) then
                        AS_atom = Cbasis_to_atom( i_basis_1)
                        HF_stress_per_atom_PU(i_coord,AS_atom)  = &
                              HF_stress_per_atom_PU(i_coord,AS_atom) &
                              + ( dens_mat(i_place) &
                              * AS_corr_shell(i_compute_2,i_compute_1) )
                      end if
                    end if
                  end if
                  AS_pulay_stress_local(i_coord) = &
                       AS_pulay_stress_local(i_coord) &
                       + dens_mat(i_place) &
                       * AS_strain_deriv_wave_shell(i_compute_2,i_compute_1)
                  if (compute_heat_flux) then
                    AS_atom = Cbasis_to_atom(i_basis_1)
                    HF_stress_per_atom_PU(i_coord,AS_atom) = &
                         HF_stress_per_atom_PU(i_coord,AS_atom) &
                         + ( 0.5d0 * dens_mat(i_place) &
                         * AS_strain_deriv_wave_shell(i_compute_2,i_compute_1) )
                    AS_atom = Cbasis_to_atom(i_basis_2)
                    HF_stress_per_atom_PU(i_coord,AS_atom) = &
                         HF_stress_per_atom_PU(i_coord,AS_atom) &
                         + ( 0.5d0 * dens_mat(i_place) &
                         * AS_strain_deriv_wave_shell(i_compute_2,i_compute_1) )
                  end if

                  if(i_basis_1/= i_basis_2)then
                    if (i_coord.le.3) then
                      sum_forces(i_coord, Cbasis_to_atom(i_basis_2)) = &
                           sum_forces(i_coord, Cbasis_to_atom(i_basis_2)) &
                           + dens_mat(i_place) &
                           * forces_shell(i_compute_1, i_compute_2)
                    end if

                    if (.not.relativistic) then
                      if ( Cbasis_to_center(i_basis(i_compute_1)) &
                           .eq. &
                           Cbasis_to_center(i_basis(i_compute_2)) )&
                      then ! on-site
                        AS_dde_T_sa_at_di_orb_local(i_coord) = &
                             AS_dde_T_sa_at_di_orb_local(i_coord) + &
                             ( dens_mat(i_place) * &
                             AS_corr_shell(i_compute_1,i_compute_2) )
                        if (compute_heat_flux) then
                          AS_atom = Cbasis_to_atom(i_basis_1)
                          HF_stress_per_atom_PU(i_coord,AS_atom) = &
                               HF_stress_per_atom_PU(i_coord,AS_atom) &
                               + ( dens_mat(i_place) &
                               * AS_corr_shell(i_compute_1,i_compute_2) )
                        end if
                      end if
                    end if
                    AS_pulay_stress_local(i_coord) = &
                         AS_pulay_stress_local(i_coord) &
                         + (dens_mat(i_place) &
                         * AS_strain_deriv_wave_shell(i_compute_1, i_compute_2))
                    if (compute_heat_flux) then
                      AS_atom = Cbasis_to_atom(i_basis_1)
                      HF_stress_per_atom_PU(i_coord,AS_atom) = &
                           HF_stress_per_atom_PU(i_coord,AS_atom) &
                           + ( 0.5d0 * dens_mat(i_place) &
                           * AS_strain_deriv_wave_shell(i_compute_1,i_compute_2) )
                      AS_atom = Cbasis_to_atom(i_basis_2)
                      HF_stress_per_atom_PU(i_coord,AS_atom) = &
                           HF_stress_per_atom_PU(i_coord,AS_atom) &
                           + ( 0.5d0 * dens_mat(i_place) &
                           * AS_strain_deriv_wave_shell(i_compute_1,i_compute_2) )
                    end if
                  end if
                  exit place_2_ed
                else if(column_index_hamiltonian( i_place) > i_basis_2)then
                   exit place_2_ed
                end if
              end do place_2_ed

              offset(i_cell) = i_place
            end if ! i_basis_2 <= i_basis_1
          end do ! i_compute_2
        end do ! i_compute_1
      end select ! PM_index
    end if ! packed_matrix_format
  end subroutine AS_update_sum_forces_and_stress_cpu
end module forces_densmat
