!****s* FHI-aims/prepare_fine_k_quantities
!  NAME
!    prepare_fine_k_quantities
!  SYNOPSIS

  subroutine prepare_fine_k_quantities ()
!  USES
    use physics
    use species_data
    use dimensions
    use constants
    use pbc_lists
    use runtime_choices
    use scaled_zora_transform
    use lapack_wrapper
    use hartree_fock
    use lvl_triples
    use tight_binding_auxmat
    use prodbas
    use synchronize_mpi_basic
    use aims_memory_tracking, only : aims_allocate, aims_deallocate
    use elsi_wrapper, only : aims_elsi_occ
    use synchronize_mpi, only: sync_eigenvalues
    use ks_wrapper, only: solve_KS_elsi_serial

!  PURPOSE
!   The subroutine compute the LVL triple coefficients and KS eigenvalues and eigenvectors on
!   a fine k grid mesh.
!
    implicit none

!  INPUTS
!    none
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8,    dimension(:,:),allocatable :: hamiltonian_w
    real*8,    dimension(:),  allocatable :: overlap_matrix_w
    complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
    complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex

    real*8, dimension(:,:,:),allocatable :: work_ham
    real*8, dimension(:,:),allocatable :: work_ovl

    real*8:: diff_electrons
    real*8:: en_xc_tmp, en_pot_xc_tmp
    integer :: n_k_points_xyz_fine(3)

    logical :: t_out = .true.

!  counter
    integer :: i_spin
    integer :: i_counter
    integer :: i_k_point, i_k_point_local, i_state, i_k, i_q_point
    integer :: i_task
    integer :: i_x, i_y, i_z
    integer :: i_cell_n

    ! for possible output of Kohn-Sham eigenvectors

    ! file name infrastructure etc.

    integer :: output_priority_old
    character*200 :: info_str

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "Computing the LVL triple coefficients and KS band structure on a finer k-point mesh:"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

!    use_hf_kspace = .false.
    call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, partition_tab, &
         l_shell_max, en_xc_tmp, en_pot_xc_tmp, hamiltonian, en_vdw, en_pot_vdw )

    if(real_eigenvectors) then
       allocate(hamiltonian_w_complex (1,1))
       allocate(overlap_matrix_w_complex(1))
       allocate(hamiltonian_w (n_basis*(n_basis+1)/2,n_spin) )
       allocate(overlap_matrix_w (n_basis*(n_basis+1)/2) )
    else
       allocate(hamiltonian_w_complex (n_basis*(n_basis+1)/2,n_spin))
       allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2))
       allocate(hamiltonian_w (1,1) )
       allocate(overlap_matrix_w (1) )
    endif

    if(packed_matrix_format == PM_none)then
       allocate(work_ham(n_centers_basis_I, n_centers_basis_I, n_spin))
       allocate(work_ovl(n_centers_basis_I, n_centers_basis_I))
    else
       ! dummy only, never touched
       allocate(work_ham( 1, 1, 1))
       allocate(work_ovl( 1, 1))
    end if

    if (flag_rel.eq.1) then
       call allocate_scaled_zora_transform
       call integrate_scaled_zora_transf_p2( &
            rho, rho_gradient, kinetic_density, hartree_potential,    &
            partition_tab, l_shell_max)
    end if

    n_k_points_xyz_fine(1) = n_k_points_xyz(1)*k_enhancement_factor(1)
    n_k_points_xyz_fine(2) = n_k_points_xyz(2)*k_enhancement_factor(2)
    n_k_points_xyz_fine(3) = n_k_points_xyz(3)*k_enhancement_factor(3)

!    n_k_points = n_k_points_xyz_fine(1)*n_k_points_xyz_fine(2)*n_k_points_xyz_fine(3)

    if(allocated(k_point_list)) then
      deallocate(k_point_list)
    endif
    allocate(k_point_list(n_k_points,3))

    if(allocated(k_phase)) then
       deallocate(k_phase)
    endif
    allocate(k_phase(n_cells,n_k_points))

    i_k_point = 0
    do i_x = 1, n_k_points_xyz_fine(1)
      do i_y = 1, n_k_points_xyz_fine(2)
        do i_z = 1, n_k_points_xyz_fine(3)

             i_k_point = i_k_point + 1

             k_point_list(i_k_point,1) = dble(i_x-1) / dble(n_k_points_xyz_fine(1)) + k_points_offset(1)
             k_point_list(i_k_point,2) = dble(i_y-1) / dble(n_k_points_xyz_fine(2)) + k_points_offset(2)
             k_point_list(i_k_point,3) = dble(i_z-1) / dble(n_k_points_xyz_fine(3)) + k_points_offset(3)

        enddo
      enddo
    enddo

    do i_k_point = 1, n_k_points, 1
       do i_cell_n = 1, n_cells, 1
          k_phase( i_cell_n, i_k_point) = exp((0,2)*pi*sum(k_point_list(i_k_point,:)*cell_index(i_cell_n,:)))
       end do
    end do


    if(allocated(k_weights)) then
       deallocate(k_weights)
    endif
    allocate(k_weights(n_k_points))
    k_weights = 1.0d0 /n_k_points

    n_k_points_task = 0
    do i_k_point = 1, n_k_points, 1
      if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then
          n_k_points_task = n_k_points_task + 1
      end if
    end do

    if(allocated(lvl_tricoeff_recip1))then
      deallocate(lvl_tricoeff_recip1)
    endif

    allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_task))
! NOTE: no need to reallocate "lvl_tricoeff_recip2", which does dot depend on k point mesh.

    call cleanup_lvl_triples()
    call deallocate_tb_auxmat()
    call initialize_lvl_triples(OVLP_TYPE_COULOMB)
    call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
    call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)
    call cleanup_lvl_triples()
    call deallocate_tb_auxmat()

!    write(use_unit,*) "KS_eigenvector"
!    do i_k = 1, n_k_points_task
!      write(use_unit,*) "ik-1", i_k
!      write(use_unit,'(100f18.8)') KS_eigenvector_complex(:,:,:,i_k)
!    enddo
!    write(use_unit,*) "lvl_tricoeff_recip2", lvl_tricoeff_recip2
!    use_hf_kspace = .false.


    if (allocated(occ_numbers))            deallocate(occ_numbers)
    if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
    if (allocated(KS_eigenvector))         call aims_deallocate( KS_eigenvector, "KS_eigenvector" )
    if (allocated(KS_eigenvector_complex)) call aims_deallocate( KS_eigenvector_complex, "KS_eigenvector_complex" )
    allocate (occ_numbers(n_states,n_spin,n_k_points))
    allocate( KS_eigenvalue(n_states,n_spin,n_k_points) )
    if(real_eigenvectors) then
       call aims_allocate( KS_eigenvector, n_basis,n_states,n_spin,n_k_points_task, "KS_eigenvector" )
       call aims_allocate( KS_eigenvector_complex, 1,1,1,1, "KS_eigenvector_complex" )
    else
       call aims_allocate( KS_eigenvector, 1,1,1,1, "KS_eigenvector" )
       call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin,n_k_points_task, "KS_eigenvector_complex" )
    endif

    i_k = 0
    do i_k_point = 1, n_k_points, 1
       if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then
          i_k = i_k + 1

          call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                hamiltonian_w, overlap_matrix_w, &
                hamiltonian_w_complex, overlap_matrix_w_complex, i_k_point, work_ham, work_ovl)

          call clear_lapack_overlap_matrix()

          ! solve KS-equations for both spins
          do i_spin = 1, n_spin, 1

             flag_KS_k_points(1) = BASIS_SINGULAR_NOT_TESTED
             ! calculate the eigenvalues and eigenvectors

             output_priority_old = output_priority
             output_priority = OL_high

             if(use_elsi .and. .not. use_scalapack) then
                if(real_eigenvectors) then
                   call solve_KS_elsi_serial(hamiltonian_w(:,i_spin),&
                        overlap_matrix_w,KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector(:,:,i_spin,i_k),i_spin,1)
                else
                   call solve_KS_elsi_serial(hamiltonian_w_complex(:,i_spin),&
                        overlap_matrix_w_complex,&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector_complex(:,:,i_spin,i_k),i_spin,1)
                end if
             else
                if(real_eigenvectors) then
                   call improve_real_eigenfunctions(overlap_matrix_w,&
                        hamiltonian_w(:,i_spin),t_out,&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector(:,:,i_spin,i_k),1)
                else
                   call improve_complex_eigenfunctions(&
                        overlap_matrix_w_complex,&
                        hamiltonian_w_complex(:,i_spin),&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector_complex(:,:,i_spin,i_k),1)
                end if
             end if

             output_priority = output_priority_old

          enddo
       else
          KS_eigenvalue(:,:,i_k_point) = 0.0d0
       end if
    end do

    call sync_eigenvalues(KS_eigenvalue)

!    do i_k = 1, n_k_points
!      write(use_unit,*) "ik-2", i_k
!      do i_state = 1, n_states
!        write(use_unit,'(100f18.8)') KS_eigenvalue(i_state,1,i_k)
!      enddo
!    enddo

    if (flag_rel.eq.1) then

        write(info_str,'()')
        call localorb_info(info_str,use_unit,'(A)')
        ! This routine performs the scaled ZORA correction for every k-point on the
        ! present MPI task, and it also re-synchronizes the eigenvalues after it is done.
        ! So we have the correct eigenvalues for all k-points on MPI task 0.
        call evaluate_scaled_zora_transf_p1(KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue)

    end if

    if(.not. fixed_spin_moment) then
       if(mu_method == 0) then ! zeroin
          call get_occupation_numbers_p0(KS_eigenvalue,n_electrons,t_out,&
                  occ_numbers,chemical_potential)

          call check_norm_p0(chemical_potential,KS_eigenvalue,n_electrons,&
                  occ_numbers,diff_electrons,i_counter)
       else ! bisection
          call aims_elsi_occ(n_electrons,n_states,n_spin,n_k_points,k_weights,&
                  KS_eigenvalue,occ_numbers,chemical_potential)
       endif

       chemical_potential_spin(:) = chemical_potential
    else ! fixed_spin_moment
       do i_spin = 1,n_spin
          if(mu_method == 0) then ! zeroin
             call get_occupation_numbers_single_channel(&
                     KS_eigenvalue(1:n_states,i_spin,1:n_k_points),&
                     fixed_spin_moment_electrons(i_spin),t_out,&
                     occ_numbers(1:n_states,i_spin,1:n_k_points),&
                     chemical_potential_spin(i_spin),i_spin)

             call check_norm_periodic_v2(chemical_potential_spin(i_spin),&
                     KS_eigenvalue(:,i_spin,:),&
                     fixed_spin_moment_electrons(i_spin),&
                     occ_numbers(:,i_spin,:),diff_electrons,i_counter,i_spin)
          else ! bisection
             call aims_elsi_occ(fixed_spin_moment_electrons(i_spin),n_states,1,&
                     n_k_points,k_weights,&
                     KS_eigenvalue(1:n_states,i_spin,1:n_k_points),&
                     occ_numbers(1:n_states,i_spin,1:n_k_points),&
                     chemical_potential_spin(i_spin))
          endif
       enddo
    endif

!    do i_k = 1, n_k_points
!      write(use_unit,*) "ik-2", i_k
!      do i_state = 1, n_states
!        write(use_unit,'(2I4, f18.8)') i_k, i_state, occ_numbers(i_state,1,i_k)
!      enddo
!    enddo

    deallocate(work_ham)
    deallocate(work_ovl)

    deallocate(hamiltonian_w_complex)
    deallocate(overlap_matrix_w_complex)
    deallocate(hamiltonian_w)
    deallocate(overlap_matrix_w)

    if (flag_rel.eq.1) then
       call deallocate_scaled_zora_transform
    end if

  end subroutine prepare_fine_k_quantities
!******

