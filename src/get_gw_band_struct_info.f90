!****s* FHI-aims/get_gw_band_struct_info
!  NAME
!    out_plot_band
!  SYNOPSIS

  subroutine get_gw_band_struct_info &
             (xc_realspace, gw_selfe_band)
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
    use hartree_fock_p0, only : hf_exchange_matr_complex
    use lvl_triples
    use tight_binding_auxmat
    use prodbas
    use gw_para
    use localorb_io
    use synchronize_mpi_basic
    use ex_mat_ksk_p0
    use lvl_tricoeff
    use synchronize_mpi, only: sync_eigenvalues
    use ks_wrapper, only: solve_KS_elsi_serial

!  PURPOSE
!   The subroutine produces data files for GW band structure plotting, modified from
!   routine "out_plot_band_hf_k_space.f90".
!
    implicit none

!    integer :: n_low_state
!    integer :: n_high_state
    real*8, dimension(n_hamiltonian_matrix_size,n_spin) ::  xc_realspace
    complex*16, dimension(n_freq,n_low_state:n_high_state,n_spin,n_band_kpoints_task) :: gw_selfe_band
!  INPUTS
! o  n_low_state  --  the lowest KS/HF eigenstate taken into account in the GW self-energy calculations
! o  n_high_state --  the highest KS/HF eigenstate taken into account
! o  xc_realspace -- the local/semilocal exchange-correlation matrix in realspace
! o  gw_selfe_band -- gw self energy on a set of k points for band plotting
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

    integer :: n_k_points_min
    real*8,    dimension(:,:),allocatable :: hamiltonian_w
    real*8,    dimension(:),  allocatable :: overlap_matrix_w
    complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
    complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex

    real*8, dimension(:,:,:),allocatable :: work_ham
    real*8, dimension(:,:),allocatable :: work_ovl

    complex*16 :: k_phase_exx_old
    real*8, dimension(:), allocatable :: k_weights_old
    real*8, dimension(:,:,:), allocatable :: occ_numbers_old
    real*8, dimension(:,:), allocatable :: fock_matrix
    complex*16, dimension(:,:), allocatable :: fock_matrix_complex
    logical :: real_eigenvectors_old

    real*8, dimension(:,:,:,:), allocatable :: KS_eigenvector_old
    complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex_old

    complex*16, dimension(:,:,:), allocatable :: lvl_tricoeff_mod_r_k
    real*8, allocatable :: exact_x_kspace(:,:,:)
    real*8, allocatable :: xc_kspace(:,:,:)
    real*8, allocatable :: qp_energy(:,:,:)
    real*8, allocatable :: qp_energy_tmp(:,:)
    complex*16, allocatable :: sigma_par_band(:,:,:,:)
    complex*16, allocatable :: sigma_par_tmp(:,:,:)


    real*8:: k(3) !, k_lowest_un_occ(3), k_highest_occ(3),  k_lowest_un_occ_whole_system(3), k_highest_occ_whole_system(3)
    real*8:: k_cbm(3), k_vbm(3)
    real*8:: lowest_un_occ,highest_occ, lowest_un_occ_whole_system, highest_occ_whole_system
    real*8:: diff_electrons
    real*8:: e_diff, en, mu, delta_mu
    integer :: n_homo_global

    complex*16 selfe, dselfe

    real*8, parameter ::  qp_energy_thr = 1.d-5

    ! for possible output of Kohn-Sham eigenvectors
    complex*16, dimension(:,:,:),allocatable :: KS_eigenvector_tmp

    ! file name infrastructure etc.

    character*10 :: num_char
    integer :: output_priority_old
    character*50 :: file_name, file_name2
    character*200 :: info_str

    character(*), parameter :: func='get_gw_band_struct_info'

! counter
    integer:: i_band,    i_cell_n, i_spin
    integer:: i_cell(3), i_counter
    integer:: i_k_point, i_k_point_local, i_state, i_k, i_q_point
    integer :: i_k_point_old, i_k_old, i_basis_1, i_basis_2
    integer :: i_index, i_index_local
    integer :: i_index_bas
    integer :: i_task, id_send
    integer :: i_count
    integer :: info
    integer:: n_k_points_band,  n_k_points_band_min, n_k_points_band_task, n_ks_points_band_task

    call perfon('getgwband')
    output_priority = OL_high

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "Writing the requested band structure output:"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

! first analytically continue the self energy
    allocate(sigma_par_band(n_max_par,n_low_state:n_high_state,n_spin,n_band_kpoints_task), stat=info)
    call check_allocation(info,'sigma_par_band',func)
    allocate(sigma_par_tmp(n_max_par,n_low_state:n_high_state,n_spin),stat=info)
    call check_allocation(info,'sigma_par_tmp',func)
    allocate(qp_energy_tmp(n_low_state:n_high_state,n_spin),stat=info)
    call check_allocation(info,'qp_energy_tmp',func)

    call analy_continue_self_energy_p0 &
        (anacon_type, n_band_kpoints, n_band_kpoints_task, &
         n_max_par, n_low_state,n_high_state,n_freq, &
         omega_grid, gw_selfe_band, &
         sigma_par_band)

    call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, partition_tab, &
         l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

    allocate(hamiltonian_w_complex   (n_basis*(n_basis+1)/2,n_spin))
    allocate(overlap_matrix_w_complex(n_basis*(n_basis+1)/2))
    allocate(fock_matrix_complex(n_basis*(n_basis+1)/2,n_spin))

    ! dummy allocations for the real arrays.
    ! These trigger a compiler warning for -check pointers otherwise
    ! as they are included (but never touched!) in construct_hamiltonian_and_ovl below.
    allocate(hamiltonian_w   (1,1) )
    allocate(overlap_matrix_w (1) )

    if(packed_matrix_format == PM_none)then
       allocate(work_ham(n_centers_basis_I, n_centers_basis_I, n_spin))
       allocate(work_ovl(n_centers_basis_I, n_centers_basis_I))
    else
       ! dummy only, never touched
       allocate(work_ham( 1, 1, 1))
       allocate(work_ovl( 1, 1))
    end if

    if (out_eigenvec) then
      allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))
    end if

    if (flag_rel.eq.1) then
       call allocate_scaled_zora_transform
       call integrate_scaled_zora_transf_p2( &
            rho, rho_gradient, kinetic_density, hartree_potential,    &
            partition_tab, l_shell_max)
    end if

    allocate(k_weights_old(n_k_points))
    k_weights_old = k_weights
    allocate(occ_numbers_old(n_states,n_spin,n_k_points))
    occ_numbers_old = occ_numbers

    real_eigenvectors_old = real_eigenvectors

    use_hf_kspace = .true.

    if(real_eigenvectors_old)then
       allocate(KS_eigenvector_old(n_basis,n_states,n_spin,n_q_points_task))
       allocate(KS_eigenvector_complex_old(1,1,1,1))
       KS_eigenvector_old = KS_eigenvector
    else
       allocate(KS_eigenvector_complex_old(n_basis,n_states,n_spin,n_q_points_task))
       allocate(KS_eigenvector_old(1,1,1,1))
       KS_eigenvector_complex_old = KS_eigenvector_complex
    endif


    lowest_un_occ_whole_system = 1d100
    highest_occ_whole_system    = -1d100

    deallocate(coulomb_matr_blacs)
    allocate(coulomb_matr_blacs(lbb_row:ubb_row,lbb_col:ubb_col, n_ks_points_task))

!    call perfon('gwbl1')
    i_index = 0
    do i_band = 1, n_plot_band

       n_k_points_band =  n_points_in_band(i_band)
       n_k_points_band_task = 0
       do i_k_point = 1, n_k_points_band, 1
          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points_band )then
             n_k_points_band_task = n_k_points_band_task + 1
          end if
       end do
       n_k_points_band_min = max(n_k_points_band_task,1)

       n_ks_points_band_task = 0
       do i_k_point = 1, n_k_points_band
          if(myid_irkq.eq.mod(i_k_point-1,n_tasks_irkq) .and. myid_irkq .le. n_k_points_band-1) then
             n_ks_points_band_task = n_ks_points_band_task + 1
          endif
       enddo

       !for the moment, n_k_points has to be redefined because it is used as global variable in some subroutines
       n_k_points=n_k_points_band
       n_k_points_task=n_k_points_band_task

       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,I4,A,I4,A)') "Treating all ",n_k_points_band," k-points in band plot segment #", i_band, ":"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')

       deallocate(k_weights)
       allocate(k_weights(n_k_points_band))
       k_weights = 1.0d0 /n_k_points_band

       deallocate(k_phase)
       allocate(k_phase(n_cells,n_k_points_band))

       do  i_k_point = 1, n_k_points_band
          k(:) = band_begin(i_band,:) +  real(i_k_point-1)/real(n_k_points_band-1) &
                *( band_end(i_band,:) -  band_begin(i_band,:))

          do i_cell_n = 1, n_cells
             k_phase( i_cell_n, i_k_point) = exp((0,2)*pi*sum(k(:)*cell_index(i_cell_n,:)))
          end do
       end do

       if(abs(sum(k_weights)-1) >1e-10)then
          write(use_unit,*) 'Error: sum of k-vector weights is not one!', sum(k_weights)
          stop
       end if

       allocate(lvl_tricoeff_mod_r_k(lbb_row:ubb_row,max_n_basis_sp,n_basis),stat=info)
       call check_allocation(info, 'lvl_tricoeff_mod_r_k', func)

       if(use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse) then
           call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
       else
           call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
       endif

       real_eigenvectors = .false.
       if (allocated(occ_numbers))            deallocate(occ_numbers)
       if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
       if (allocated(KS_eigenvector_complex)) deallocate(KS_eigenvector_complex)
       if (allocated(hf_exchange_matr_complex)) deallocate(hf_exchange_matr_complex)
       allocate (occ_numbers(n_states,n_spin,n_k_points_band))
       allocate( KS_eigenvalue(n_states,n_spin,n_k_points_band) )
       allocate( hf_exchange_matr_complex(n_basis,n_basis,n_k_points_band_min,n_spin) )

       allocate( KS_eigenvector_complex(n_basis,n_states,n_spin, n_k_points_band_min))

!       call perfon('gwbl2')
       i_k = 0
       do i_k_point = 1, n_k_points_band
!          call perfon('gwbl211')
          call gw_get_single_lvl_tricoeff_noev(n_cells, k_phase(:,i_k_point), lvl_tricoeff_mod_r_k)
!          call perfoff
! SVL add Fock matrix
!          call perfon('gwbl212')
          fock_matrix_complex = (0d0,0d0)
          k(:) = band_begin(i_band,:) +  real(i_k_point-1)/real(n_k_points_band-1) &
                  *( band_end(i_band,:) -  band_begin(i_band,:))
             ! Calculate q-dependent Coulomb matrix for k-q points
          do i_q_point = 1, n_q_points
              k_minus_q_point_list(i_q_point,:) = k(:) - k_point_list(i_q_point,:)
          enddo

          call get_full_coulomb_matr_blacs(coulomb_matr_blacs,0)
!          call perfoff
          call perfon('gwbl213')
          call my_evaluate_exchange_matr_kspace_single_kpoint_p0 &
               (n_k_points_band, n_ks_points_band_task, i_k_point,KS_eigenvector_old,KS_eigenvector_complex_old,&
               real_eigenvectors_old,occ_numbers_old,k_weights_old,&
               lvl_tricoeff_mod_r_k, fock_matrix_complex)
          if ((i_k_point.eq.1).and.(myid.eq.1)) then
             print*,'FMOUT: ', sum(abs(fock_matrix_complex)),sum(fock_matrix_complex)
          end if
          call perfoff
 !         call perfon('gwbl22')
          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points_band )then

             i_k = i_k + 1

             i_index_bas = 0
             do i_basis_1=1, n_basis, 1
               do i_basis_2=1, i_basis_1, 1
                i_index_bas = i_index_bas + 1
                hf_exchange_matr_complex(i_basis_2,i_basis_1,i_k,:) = fock_matrix_complex(i_index_bas,:)
                hf_exchange_matr_complex(i_basis_1,i_basis_2,i_k,:) = conjg(fock_matrix_complex(i_index_bas,:))
               enddo
             enddo

             call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                  hamiltonian_w, overlap_matrix_w, &
                  hamiltonian_w_complex, overlap_matrix_w_complex, i_k_point, work_ham, work_ovl)

!             if(use_periodic_hf)&
!                  hamiltonian_w_complex = hamiltonian_w_complex - hybrid_coeff*fock_matrix_complex

             call clear_lapack_overlap_matrix()

             ! solve KS-equations for both spins
             do i_spin = 1, n_spin, 1

                flag_KS_k_points(1) = BASIS_SINGULAR_NOT_TESTED
                ! calculate the eigenvalues and eigenvectors

                output_priority_old = output_priority
                output_priority = OL_high

                if(use_elsi .and. .not. use_scalapack) then
                   call solve_KS_elsi_serial(hamiltonian_w_complex(:,i_spin),&
                        overlap_matrix_w_complex,&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector_complex(:,:,i_spin,i_k),i_spin,1)
                else
                   call improve_complex_eigenfunctions(&
                        overlap_matrix_w_complex,&
                        hamiltonian_w_complex(:,i_spin),&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector_complex(:,:,i_spin,i_k),1)
                end if

                output_priority = output_priority_old

             enddo
          else
             KS_eigenvalue(:,:,i_k_point) = 0.0d0
          end if
!          call perfoff
       end do
!       call perfoff
       deallocate(lvl_tricoeff_mod_r_k)

       call sync_eigenvalues(KS_eigenvalue)

       if (flag_rel.eq.1) then

             write(info_str,'()')
             call localorb_info(info_str,use_unit,'(A)')
             ! This routine performs the scaled ZORA correction for every k-point on the
             ! present MPI task, and it also re-synchronizes the eigenvalues after it is done.
             ! So we have the correct eigenvalues for all k-points on MPI task 0.
             call evaluate_scaled_zora_transf_p1(KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue)

       end if

       if (.not.fixed_spin_moment) then

         call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)

       else ! if fixed_spin_moment

         do i_spin = 1, n_spin, 1

           call check_norm_periodic_v2(chemical_potential_spin(i_spin), KS_eigenvalue(:,i_spin,:), &
               fixed_spin_moment_electrons(i_spin), occ_numbers(:,i_spin,:), diff_electrons,&
               i_counter,i_spin)

         enddo

       end if

       allocate(exact_x_kspace(n_low_state:n_high_state,n_spin,n_k_points_band_min),stat=info)
       call check_allocation(info,'exact_x_kspace',func)
       allocate(xc_kspace(n_low_state:n_high_state,n_spin,n_k_points_band_min),stat=info)
       call check_allocation(info,'xc_kspace',func)

       use_inv_symmetry = .false.
       call evaluate_ex_and_xc_matr_kspace &
           ( n_k_points_band, n_k_points_band_min, n_low_state, n_high_state, &
             KS_eigenvector, KS_eigenvector_complex, &
             xc_realspace, exact_x_kspace, xc_kspace &
           )
!       call evaluate_exx_matr_kspace &
!           ( n_k_points_band, n_k_points_band_min, n_low_state, n_high_state, &
!             KS_eigenvector, KS_eigenvector_complex, &
!             exact_x_kspace &
!           )
!       call evaluate_xc_matr_kspace &
!            (n_k_points_band, n_k_points_band_min, n_low_state, n_high_state, &
!             KS_eigenvector, KS_eigenvector_complex, &
!             xc_realspace, xc_kspace &
!            )

       allocate(qp_energy(n_low_state:n_high_state,n_spin,n_k_points_band),stat=info)
       call check_allocation(info,'qp_energy',func)
! qusiparticle energy calculation
       do i_k_point = 1, n_k_points_band, 1
         i_index = i_index + 1
         i_index_local = (i_index-1)/n_tasks + 1
         i_task = mod(i_k_point, n_tasks)
         id_send = mod(i_index, n_tasks)

         if(myid .eq. id_send) then
            sigma_par_tmp(:,:,:) = sigma_par_band(:,:,:,i_index_local)
         endif

         if(id_send .ne. i_task) then
             if(myid.eq.id_send) then
               call send_complex_vector(sigma_par_tmp, n_max_par*(n_high_state-n_low_state+1)*n_spin, i_task)
             elseif(myid.eq.i_task) then
               call receive_complex_vector(sigma_par_tmp, n_max_par*(n_high_state-n_low_state+1)*n_spin, id_send)
             endif
         endif

         if (myid .eq. mod(i_k_point, n_tasks) ) then

            i_k_point_local = (i_k_point-1)/n_tasks + 1

            qp_energy(n_low_state:n_high_state,:,i_k_point)= &
                   KS_eigenvalue(n_low_state:n_high_state,:,i_k_point)

            qp_energy_tmp(n_low_state:n_high_state,:)= &
            KS_eigenvalue(n_low_state:n_high_state,:,i_k_point)

            do i_spin = 1, n_spin, 1
              do i_state = n_low_state, n_high_state, 1

                e_diff = 1.d-3
                i_count =0
                do while (abs(e_diff).gt.qp_energy_thr)
                  i_count = i_count +1
                  qp_energy(i_state,i_spin,i_k_point) = &
                        qp_energy_tmp(i_state,i_spin) + 0.5d0* e_diff
                  qp_energy_tmp(i_state,i_spin) = qp_energy(i_state,i_spin,i_k_point)

                  mu =  chemical_potential_spin(i_spin)

                  en = qp_energy(i_state,i_spin,i_k_point) - mu

                  call get_real_selfenergy(anacon_type,n_freq,omega_grid, &
                            dcmplx(en,0.d0), n_max_par, &
                            sigma_par_tmp(1:n_max_par,i_state,i_spin), selfe)


                  qp_energy(i_state,i_spin,i_k_point) = &
                             KS_eigenvalue(i_state,i_spin,i_k_point) &
                           + real(selfe) &
                           + exact_x_kspace(i_state,i_spin,i_k_point_local) &
                           - xc_kspace(i_state,i_spin,i_k_point_local)

!                  if(i_k_point .eq. 7) then
!                    write(use_unit,'(I4,5f16.8)') i_state,  KS_eigenvalue(i_state,i_spin,i_k_point), &
!                               exact_x_kspace(i_state,i_spin,i_k_point_local), &
!                               xc_kspace(i_state,i_spin,i_k_point_local), &
!                               real(selfe), &
!                               qp_energy(i_state,i_spin,i_k_point)
!                  endif

                  e_diff =  qp_energy(i_state,i_spin,i_k_point) &
                           - qp_energy_tmp(i_state,i_spin)

                 if(i_count .gt. 100) then
                      write(use_unit,'(2X,3A,I4,A,I4 )') &
                     " * Error: QUASI_PARTILCE_ENERGY: self-consistent", &
                     " quasiparticle solution can not be found for", &
                     " i_state = ",  i_state, "  i_spin = ", i_spin
                  exit
                 endif

   ! end of do while
              enddo

  !             write(use_unit,'(4I4,4f18.6)')myid, i_state, i_band, i_k_point, KS_eigenvalue(i_state,i_spin,i_k_point), &
  !                    exact_x_kspace(i_state,i_spin,i_k_point_local), &
  !                    xc_kspace(i_state,i_spin,i_k_point_local), &
  !                    qp_energy(i_state,i_spin,i_k_point)
   ! end of do i_state
             enddo
   ! end of do i_spin
            enddo

           else
             qp_energy(:,:,i_k_point) = 0.d0
! end of if myid == mod(i_k_point, n_tasks)
           endif
   ! end of loop over i_k_point
          enddo
          call sync_vector(qp_energy(n_low_state:n_high_state, :, :), (n_high_state-n_low_state+1)*n_spin*n_k_points_band)

!       do i_k_point = 1, n_k_points_band, 1
!         i_k_point_local = (i_k_point-1)/n_tasks +1
!         if(myid.eq.mod(i_k_point, n_tasks)) then
!           write(use_unit,*)i_band, i_k_point
!           do i_state = n_low_state, n_high_state, 1
!             write(use_unit,'(I4,3f18.8)')i_state, KS_eigenvalue(i_state,1, i_k_point), &
!                     exact_x_kspace(i_state, 1, i_k_point_local), &
!                     xc_kspace(i_state, 1, i_k_point_local)
!           enddo
!         endif
!       enddo

       if(myid==0)then

          if(i_band < 10)then
             write(file_name, '(A10,I1,A4)') 'GW_band100', i_band,'.out'
          else if(i_band < 100)then
             write(file_name, '(A9,I2,A4)') 'GW_band10', i_band,'.out'
          else if(i_band < 1000)then
             write(file_name, '(A8,I3,A4)') 'GW_band1', i_band,'.out'
          else
             write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
             stop
          end if

          open(88, file=file_name)

          if(n_spin >1)then

             if(i_band < 10)then
                write(file_name2, '(A10,I1,A4)') 'GW_band200', i_band,'.out'
             else if(i_band < 100)then
                write(file_name2, '(A9,I2,A4)') 'GW_band20', i_band,'.out'
             else if(i_band < 1000)then
                write(file_name2, '(A8,I3,A4)') 'GW_band2', i_band,'.out'
             else
                write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
                stop
             end if
             open(89, file=file_name2)

          end if


          do  i_k_point = 1,  n_k_points_band

             k(:) = band_begin(i_band,:) +  real(i_k_point-1)/real(n_k_points_band-1) &
                   *( band_end(i_band,:) -  band_begin(i_band,:))
!	     k(:) = modulo(k(:),1.d0)

             write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

             do  i_state = n_low_state, n_high_state

                write(88,'(F12.5,F15.5)',ADVANCE='NO') occ_numbers(i_state,1,i_k_point), &! *n_k_points_band, &
                     (qp_energy(i_state,1,i_k_point)-chemical_potential)* hartree
             end do
             write(88,'()')

             if(n_spin == 2) then

                write(89,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

                do  i_state = n_low_state, n_high_state

                   write(89,'(F12.5,F15.5)',ADVANCE='NO') occ_numbers(i_state,2,i_k_point), &!*n_k_points_band, &
                        (qp_energy(i_state,2,i_k_point)-chemical_potential)* hartree

                end do
                write(89,'()')
             end if

          end do

          close(88)

          if(n_spin==2) close(89)

       end if

       lowest_un_occ = 1d100
       highest_occ    = -1d100

       n_homo_global = nint(n_electrons/2)
!       write(use_unit,*) "n_electrons, n_homo_global", n_electrons, n_homo_global
       do  i_spin = 1,  n_spin
          do  i_k_point = 1,  n_k_points_band

             k(:) = band_begin(i_band,:) +  real(i_k_point-1)/real(n_k_points_band-1) &
                   *( band_end(i_band,:) -  band_begin(i_band,:))

             do  i_state = n_low_state,  n_high_state

!                if( qp_energy(i_state,i_spin,i_k_point) >  chemical_potential)then
                if( i_state .gt.  n_homo_global)then
                   if(lowest_un_occ .gt. qp_energy(i_state,i_spin,i_k_point) .and.  &
                      lowest_un_occ_whole_system .gt. qp_energy(i_state,i_spin,i_k_point) ) then
                      k_cbm(:)=k(:)
                   endif
                   lowest_un_occ = min(lowest_un_occ,  qp_energy(i_state,i_spin,i_k_point))
                else
                   if(highest_occ .lt. qp_energy(i_state,i_spin,i_k_point) .and.  &
                      highest_occ_whole_system .lt. qp_energy(i_state,i_spin,i_k_point)) then
                      k_vbm(:)=k(:)
                   endif
                   highest_occ = max(highest_occ, qp_energy(i_state,i_spin,i_k_point))
                end if

             end do
          end do
       end do

       if(lowest_un_occ < lowest_un_occ_whole_system) then
          lowest_un_occ_whole_system = lowest_un_occ
       end if


       if(highest_occ > highest_occ_whole_system)      highest_occ_whole_system  = highest_occ

       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')

       write(info_str,'(2X,A,I4)')  '"GW Band gap" along reciprocal space direction number: ',i_band
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F18.6,A)') '| Lowest unoccupied state:',lowest_un_occ* hartree, ' eV'

       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F18.6,A)') '| Highest occupied state :',highest_occ* hartree,  ' eV'
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F18.6,A)') '| Energy difference      :', (lowest_un_occ- highest_occ) * hartree,  ' eV'
       call localorb_info ( info_str )

       write(info_str,'(A)') ''
       call localorb_info ( info_str )

       call deallocate_tb_auxmat()

       deallocate (occ_numbers)
       deallocate( KS_eigenvalue)
       deallocate( KS_eigenvector_complex)
       deallocate( hf_exchange_matr_complex)
       deallocate(exact_x_kspace)
       deallocate(xc_kspace)
       deallocate(qp_energy)

    end do ! i_band
!    call perfoff

    write(info_str,'(2X,A)')  '"GW Band gap" of total set of bands: '
    call localorb_info ( info_str )

    write(info_str,'(2X,A,1X,F16.8,A,3F13.5)') '| Lowest unoccupied state:', &
       lowest_un_occ_whole_system* hartree, ' eV at k point ', k_cbm(:)
    call localorb_info ( info_str )

    write(info_str,'(2X,A,1X,F16.8,A,3F13.5)') '| Highest occupied state :', &
       highest_occ_whole_system* hartree, ' eV at k point ', k_vbm(:)
    call localorb_info ( info_str )

    write(info_str,'(2X,A,1X,F16.8,A)') '| Energy difference      :', &
      (lowest_un_occ_whole_system- highest_occ_whole_system) * hartree,' eV'
    call localorb_info ( info_str )

    write(info_str,'(A)') ''
    call localorb_info ( info_str )


    deallocate(work_ham)
    deallocate(work_ovl)

    if (allocated(KS_eigenvector_tmp)) then
      deallocate(KS_eigenvector_tmp)
    end if

    if(allocated(KS_eigenvector_old)) deallocate(KS_eigenvector_old)
    if(allocated(KS_eigenvector_complex_old)) deallocate(KS_eigenvector_complex_old)

    deallocate(hamiltonian_w_complex)
    deallocate(overlap_matrix_w_complex)
    deallocate(hamiltonian_w)
    deallocate(overlap_matrix_w)
    deallocate(k_weights_old)
    deallocate(occ_numbers_old)
    deallocate(fock_matrix_complex)
    deallocate(sigma_par_band)
    deallocate(sigma_par_tmp)
    deallocate(qp_energy_tmp)

    if (flag_rel.eq.1) then

       call deallocate_scaled_zora_transform
    end if

    call perfoff

  end subroutine get_gw_band_struct_info
!******

