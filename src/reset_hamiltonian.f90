!****s* FHI-aims/reset_hamiltonian
!  NAME
!   reset_hamiltonian
!  SYNOPSIS

subroutine reset_hamiltonian(converged)

  !  PURPOSE
  !
  !  If the number of k-points in- or decreases due to new symmetries, KS-eigenvalues 
  !  and -vectors have to be reallocated. The can't be used for reinitilizing the 
  !  SCF. We have to start from scratch with an initial guess from free atoms.
  !
  !  USES
  use runtime_choices, only: use_local_index, use_scalapack, PM_none,&
                             out_aitranss, out_matrices_format_2005, out_matrices,&
			     PM_index, packed_matrix_format, first_integration,&
			     flag_KS, real_eigenvectors, collect_eigenvectors,&
			     flag_KS_k_points, condition_penalty
  use localorb_io, only: localorb_info, use_unit, OL_norm
  use vdw_correction, only: hirsh
  use dimensions, only: n_periodic, n_spin, n_hamiltonian_matrix_size, use_plus_u,&
                        n_k_points_task, n_states, n_basis, n_k_points,&
                        use_vdw_correction_hirshfeld
  use physics,only: hamiltonian, KS_eigenvector, KS_eigenvector_complex,&
                    KS_eigenvalue, overlap_matrix, hartree_potential,&
                    en_xc, en_pot_xc, en_vdw, en_pot_vdw, occ_numbers, rho,&
		    rho_gradient, kinetic_density, partition_tab,&
                    total_energy, fock_energy, entropy_correction,en_pbe_c,&
                    en_lda_c, en_ll_vdw_err, en_ll_vdw,en_density_embed,&
                    en_ion_embed, en_ion_ion,en_elec_free,hartree_energy_free,&
                    pot_ion_embed,rho_change,en_elec_delta,hartree_multipole_correction,&
                    ev_sum_shifted,ev_sum,chemical_potential_spin,chemical_potential,&
                    n_electrons, av_core_shift, hartree_delta_energy,&
                    ionic_forces, ext_charge_forces, reshape_matrices
  use species_data,only:l_shell_max
  use pbc_lists, only: remove_small_numbers_in_hamiltonian_and_ovlp
  use plus_u, only: plus_u_init_idx
  use mpi_tasks, only: check_allocation, aims_stop
  use scalapack_wrapper, only: set_sparse_local_ham_scalapack,&
                               construct_hamiltonian_scalapack,&
                               setup_hamiltonian_scalapack
  use KH_core_states, only: evaluate_KH_core_states, add_KH_states
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  implicit none

  !  ARGUMENTS
  logical :: converged 
  !  INPUTS
  !    o   converged - SCF convergence parameter
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables
  real*8, allocatable :: ham_ovlp_work(:)
  
  character(len=120)   :: info_str
  
  integer :: info
  logical :: info_t
  
  real*8 :: penalty_energy

  write(info_str,'(2X,A)') ""
  call localorb_info(info_str,use_unit,'(A)',OL_norm)    
  write(info_str,'(2X,A)') "| Recalculating the Hamiltonian because number of k-points changed."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)
  write(info_str,'(2X,A)') ""
  call localorb_info(info_str,use_unit,'(A)',OL_norm)        
     ! redo the Hamiltonian
     ! deallocate hamiltonian
        if (.not. use_local_index .and. use_scalapack) then  
            call aims_deallocate( hamiltonian,                                 "hamiltonian" )
            call aims_allocate( hamiltonian, n_hamiltonian_matrix_size,n_spin, "hamiltonian" )
        end if
      


        ! Integrate the  Hamiltonian matrix
        first_integration = .true. ! new
        if(use_local_index) then

           ! TODO: Use as calibration for load balancing
           call aims_allocate( ham_ovlp_work, n_hamiltonian_matrix_size*n_spin, "ham_ovlp_work" )

           call integrate_real_hamiltonian_matrix_p2 &
                 ( hartree_potential,  &
                  rho, rho_gradient, kinetic_density, &
                  partition_tab, l_shell_max, &
                  en_xc, en_pot_xc, ham_ovlp_work, en_vdw, en_pot_vdw )

           call set_sparse_local_ham_scalapack(ham_ovlp_work)
           call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )
        else
            call integrate_real_hamiltonian_matrix_p2 &
                 ( hartree_potential,  &
                 rho, rho_gradient, kinetic_density, &
                 partition_tab, l_shell_max, &
                 en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

        endif

      if(.not.use_local_index .and. packed_matrix_format == PM_index)then
         call remove_small_numbers_in_hamiltonian_and_ovlp(hamiltonian, overlap_matrix,info_t)
         if(info_t) call reshape_matrices
      end if

      if (.not.use_local_index .and. out_matrices) then
         if (out_matrices_format_2005) then
            call output_real_matrices &
                 ( overlap_matrix, hamiltonian )
         else
            call output_real_matrices_aij &
                 ( overlap_matrix, hamiltonian )
         end if
      end if

!     >>> AB, feb 2012
      if (.not.use_local_index .and. out_aitranss) then
        call output_ka_overlap(overlap_matrix)
      end if
!     <<< done with insert: AB, feb 2012

      if(use_plus_u)then
         call plus_u_init_idx
      endif

      ! If use_scalapack is set, allocate and set the distributed overlap matrix now.
      ! Under certain circumstances we can deallocate the (non-distributed) overlap_matrix
      ! here and save some non-scalable memory
      if (.not.use_local_index .and. use_scalapack) then
         ! If scalapack is used, the hamiltonian has to be stored
         ! into the scalapack arrays

         if(n_periodic>0 .or. packed_matrix_format /= PM_none ) then
            call construct_hamiltonian_scalapack(hamiltonian)
         else
            call setup_hamiltonian_scalapack(hamiltonian)
         endif
      endif
      !     After initial integrations, solve eigenvalue problem for the first time ...
!     ... UNLESS this is a restart run, where the initial eigenvectors are obtained from
!     a separate file anyway

      
    
 !     we obtain the Kohn-Sham eigenvalues from
            ! a Kohn_Sham solution

	if (allocated(KS_eigenvalue)) deallocate(KS_eigenvalue)
        allocate( KS_eigenvalue(n_states,n_spin,n_k_points),stat=info )
        call check_allocation(info, 'KS_eigenvalue                 ')
        
	if (allocated(KS_eigenvector)) call aims_deallocate(KS_eigenvector, "KS_eigenvector" )
	if (allocated(KS_eigenvector_complex)) call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex" )
	if(collect_eigenvectors)then
	if(real_eigenvectors)then
	      call aims_allocate( KS_eigenvector, n_basis,n_states,n_spin,n_k_points_task, "KS_eigenvector" )
	else
	      call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin,n_k_points_task, "KS_eigenvector_complex" )
	end if
	end if
	if (.not.allocated(KS_eigenvector)) call aims_allocate(KS_eigenvector, 1,1,1,1, "KS_eigenvector") ! allocate dummies
	if (.not.allocated(KS_eigenvector_complex)) call aims_allocate(KS_eigenvector_complex, 1,1,1,1, "KS_eigenvector_complex")
 
        if (allocated(occ_numbers)) deallocate(occ_numbers)
        allocate (occ_numbers(n_states,n_spin,n_k_points),stat=info)
        call check_allocation(info, 'occ_numbers                   ') 
        
	if (allocated(flag_KS_k_points)) deallocate (flag_KS_k_points)       
        allocate (flag_KS_k_points(n_k_points),stat=info)
        call check_allocation(info, 'flag_KS_k_points              ')
        flag_KS_k_points = flag_KS
        
 !  1st reset values : new
            KS_eigenvalue(:,:,:)=0.d0
            KS_eigenvector(:,:,:,:)=0.d0
            KS_eigenvector_complex(:,:,:,:)=0.d0
            occ_numbers(:,:,:) = 0.0d0
            chemical_potential=0.d0           
            chemical_potential_spin(:)=0.d0 ! new

! use newer version of get_KS_orbitals_p0
            call advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, &
                    KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                    occ_numbers, chemical_potential, chemical_potential_spin)

!            call get_KS_orbitals_p0 &
!            ( overlap_matrix, hamiltonian, n_electrons, &
!              KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
!              occ_numbers, chemical_potential, chemical_potential_spin &
!            )

      call evaluate_KH_core_states(   overlap_matrix, rho, &
           rho_gradient, kinetic_density, &
           partition_tab, hartree_potential ) !new

      call add_KH_states(KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex)!new
      
      call output_real_eigenfunctions &
      ( KS_eigenvalue, KS_eigenvector, occ_numbers)

!     determine sum of eigenvalues
      call get_ev_sum_p0 &
           ( occ_numbers, KS_eigenvalue, av_core_shift, ev_sum, &
           ev_sum_shifted )
      
           
       !continue to prepare for the actual scf cycle

      if(use_vdw_correction_hirshfeld) then
         write(info_str,'(2X,A)') &
              "Evaluating Hirshfeld Density Update"
         call localorb_info(info_str, use_unit, '(A)', OL_norm)
         call hirsh()
      endif

        converged = .false.
        hartree_delta_energy = 0.d0
        hartree_multipole_correction = 0.d0
        en_elec_delta = 0.d0

!       initialize self-consistency criteria
          deallocate (rho_change)
          allocate(rho_change(n_spin),stat=info)
          call check_allocation(info, 'rho_change                    ')
          
       call get_free_superpos_energy_p1 &   ! new
          ( partition_tab, &
            rho, rho_gradient, kinetic_density, pot_ion_embed, &
            hartree_energy_free, en_elec_free, en_xc, en_pot_xc, en_ion_ion, &
            en_ion_embed, en_density_embed, en_vdw, en_ll_vdw, en_ll_vdw_err, &
            en_lda_c, en_pbe_c, ionic_forces, ext_charge_forces &
          )                
                
        ! Compute total energy of first iteration

!       if finite-width smearing was used, must get entropy correction for T=0 total energy
        call get_entropy_correction_p1(KS_eigenvalue, occ_numbers, &
             chemical_potential, chemical_potential_spin, entropy_correction)

        call get_penalty_energy(KS_eigenvector, KS_eigenvector_complex, &
        &                       occ_numbers, condition_penalty, penalty_energy)
        
        fock_energy=0.d0 ! new
        call get_total_energy &
        ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
          en_ion_embed, en_density_embed, &
          en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, en_elec_free, &
          en_elec_delta, fock_energy, &
          hartree_energy_free, hartree_delta_energy, &
          hartree_multipole_correction, &
          entropy_correction, penalty_energy, total_energy &
        )
        
        first_integration = .false.

  ! Allocatable arrays that are tracked
  if (allocated(ham_ovlp_work)) call aims_deallocate( ham_ovlp_work, "ham_ovlp_work" )
                
   endsubroutine reset_hamiltonian      

