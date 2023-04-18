!****s* FHI-aims/out_plot_band_hf_k_space
!  NAME
!    out_plot_band
!  SYNOPSIS

  subroutine out_plot_band_hf_k_space ( )
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
    use load_balancing, only : use_batch_permutation, batch_perm, n_bp_integ
    use soc_utilities, only : perform_soc_perturbation, &
                              convert_sr_to_soc_environment, &
                              revert_soc_to_sr_environment
    use aims_memory_tracking, only : aims_allocate, aims_deallocate
    use hdf5_output, only: output_complex_eigenvector
    use dimensions_soc, only : n_states_soc
    use timing_core, only: get_timestamps, get_times
    use timing, only: tot_clock_time_band_dos, tot_time_band_dos
    use synchronize_mpi, only: sync_eigenvalues, sync_eigenvector_complex
    use synchronize_mpi_basic, only: sync_vector
    use ks_wrapper, only: solve_KS_elsi_serial
!  PURPOSE
!   The subroutine plots band structure. This works only with lapack type of eigenvectors.
!   The routine can be called only after self consistant iterations, because it destrois
!   the original k-point information.
!
    implicit none
!  INPUTS
!
!   None
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

    integer:: i_band,    i_cell_n, i_spin
    integer:: i_cell(3), i_counter
    integer:: i_k_point, i_state, j_state, i_k, i_q_point, i_k_index
    integer :: i_k_point_old, n_k_points_old, i_index, i_k_old, i_basis_1, i_basis_2
    complex*16 :: k_phase_exx_old
    real*8, dimension(:), allocatable :: k_weights_old
    real*8, dimension(:,:,:), allocatable :: occ_numbers_old
    real*8, dimension(:,:), allocatable :: fock_matrix
    complex*16, dimension(:,:), allocatable :: fock_matrix_complex
    logical :: real_eigenvectors_old
    complex*16, dimension(:,:,:,:), allocatable :: lvl_tricoeff_recip1_new
    real*8, dimension(:,:,:), allocatable :: lvl_tricoeff_recip2_new
    real*8, dimension(:,:,:,:), allocatable :: KS_eigenvector_old
    complex*16, dimension(:,:,:,:), allocatable :: KS_eigenvector_complex_old

    real*8:: k(3) !, k_lowest_un_occ(3), k_highest_occ(3),  k_lowest_un_occ_whole_system(3), k_highest_occ_whole_system(3)
    real*8:: lowest_un_occ,highest_occ, lowest_un_occ_whole_system, highest_occ_whole_system
    real*8:: lowest_un_occ_whole_system_soc, highest_occ_whole_system_soc
    real*8:: diff_electrons

    ! for possible output of Kohn-Sham eigenvectors
    complex*16, dimension(:,:,:),allocatable :: KS_eigenvector_tmp

    ! file name infrastructure etc.

    character*10 :: num_char
    integer :: output_priority_old
    character*50 :: file_name, file_name2, sr_suffix
    character*100 :: info_str

    interface
       subroutine evaluate_exchange_matr_kspace_single_kpoint_p0 &
            (k_point,KS_egnv,KS_egnv_complex,real_eigen,occ_numbers,&
            q_weights,lvl_tricoeff_recip1_new,lvl_tricoeff_recip2_new,fock_m)
         integer :: k_point
         real*8, dimension(:,:,:,:) :: KS_egnv
         complex*16, dimension(:,:,:,:) :: KS_egnv_complex
         logical :: real_eigen
         real*8, dimension(:,:,:) :: occ_numbers
         real*8, dimension(:) :: q_weights
         complex*16, intent(IN)  :: lvl_tricoeff_recip1_new(:,:,:,:)
         real*8, intent(IN)      :: lvl_tricoeff_recip2_new(:,:,:)
         complex*16, intent(OUT) :: fock_m(:,:)
       end subroutine evaluate_exchange_matr_kspace_single_kpoint_p0
    end interface

    ! Variables for periodic SOC
    complex*16, dimension(:,:), allocatable :: SOC_Hamiltonian
    real*8,dimension(:,:), allocatable :: soc_matrix
    complex*16, dimension(:,:), allocatable :: eigenvec_soc_wf_basis
    real*8 :: dummy
    integer :: i, dummy_int, this_k_point, info
    ! my_k_points is an array for converting between local k-point
    ! indexing (i.e. the k-point indexing of KS_eigenvector) and
    ! global/shared k-point indexing (i.e. the k-point indexing of
    ! KS_eigenvalue)
    integer, dimension(:), allocatable :: my_k_points

    real*8 :: band_time = 0.d0
    real*8 :: clock_band_time = 0.d0

    call get_timestamps(band_time, clock_band_time)

    if (use_local_index) then
      call aims_stop_coll("* Error:  use_local_index not supported for exx_band_structure_version 2.  Exiting.")
    end if

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "Writing the requested band structure output (EXX k-space version):"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

    call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, &
         partition_tab, l_shell_max, &
         en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw)

    ! If SOC enabled, set up initial matrices
    if (calculate_perturbative_soc) then
      allocate( SOC_Hamiltonian( n_states_soc, n_states_soc ),stat=info )
      call check_allocation(info, 'SOC_Hamiltonian               ')

      if (use_local_index.and.use_load_balancing) then
        allocate( soc_matrix(batch_perm(n_bp_integ)%n_local_matrix_size, 3) )
        use_batch_permutation = n_bp_integ
      else
        allocate( soc_matrix(n_hamiltonian_matrix_size, 3) )
      end if
      call check_allocation(info, 'soc_matrix                    ')

      allocate( eigenvec_soc_wf_basis( n_states_soc, n_states_soc ), stat=info )
      call check_allocation(info, 'eigenvec_soc_wf_basis         ')

      call integrate_soc_matrix (rho, hartree_potential, partition_tab, soc_matrix )
   else ! otherwise SOC code will never be touched, set up dummy indices
      allocate( SOC_Hamiltonian(1,1),stat=info )
      call check_allocation(info, 'SOC_Hamiltonian               ')

      allocate( soc_matrix(1,1),stat=info )
      call check_allocation(info, 'soc_matrix                    ')

      allocate( eigenvec_soc_wf_basis( 1, 1 ), stat=info )
      call check_allocation(info, 'eigenvec_soc_wf_basis         ')

      allocate( my_k_points(1),stat=info )
      call check_allocation(info, 'my_k_points                   ')
    end if

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
    n_k_points_old = n_k_points ! SVL store for periodic exx bands
    allocate(k_weights_old(n_k_points_old))
    k_weights_old = k_weights
    allocate(occ_numbers_old(n_states,n_spin,n_k_points_old))
    occ_numbers_old = occ_numbers

    real_eigenvectors_old = real_eigenvectors

    ! k-space HF implementation will be used
    if(use_periodic_hf)then
       use_hf_kspace = .true.
       call allocate_hartree_fock()
       call initialize_lvl_triples(OVLP_TYPE_COULOMB)
       ! extend the cut-Coulomb operator
!       cutCb_rcut = cutCb_rcut*2
       call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
       call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)
       call cleanup_lvl_triples()
       call deallocate_tb_auxmat()
    endif
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

    do i_band = 1, n_plot_band

       n_k_points =  n_points_in_band(i_band)

       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,I4,A,I4,A)') "Treating all ",n_k_points," k-points in band plot segment #", i_band, ":"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')

       deallocate(k_weights)
       allocate(k_weights(n_k_points))
       k_weights = 1.0d0 /n_k_points

       deallocate(k_phase)
       allocate(k_phase(n_cells,n_k_points))

       do  i_k_point = 1, n_k_points
          k(:) = band_k_frac(i_k_point, i_band)

          do i_cell_n = 1, n_cells
             k_phase( i_cell_n, i_k_point) = exp((0,2)*pi*sum(k(:)*cell_index(i_cell_n,:)))
          end do
       end do

       if(abs(sum(k_weights)-1) >1e-10)then
          write(use_unit,*) 'Error: sum of k-vector weights is not one!', sum(k_weights)
          stop
       end if

       n_k_points_task = 0
       do i_k_point = 1, n_k_points, 1
          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then
             n_k_points_task = n_k_points_task + 1
          end if
       end do

       if (use_periodic_hf) then
          allocate(lvl_tricoeff_recip1_new(max_n_basbas_sp,n_basis,n_basis,n_k_points_task))
          allocate(lvl_tricoeff_recip2_new(max_n_basbas_sp,n_basis,n_basis))
          call initialize_lvl_triples(OVLP_TYPE_COULOMB)
          call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
          call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1_new,lvl_tricoeff_recip2_new)
          call cleanup_lvl_triples()
          call deallocate_tb_auxmat()
          if(use_hse .and. hse_omega_hf /= 0.d0 .and. .not. use_gw_and_hse) then
             call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
          else
             call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
          endif
       endif



       real_eigenvectors = .false.
       if (allocated(occ_numbers))            deallocate(occ_numbers)
       if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
       if (allocated(KS_eigenvector_complex)) call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex")
       allocate (occ_numbers(n_states,n_spin,n_k_points))
       allocate( KS_eigenvalue(n_states,n_spin,n_k_points) )

       if ( (flag_rel.eq.1) .or. (out_eigenvec) .or. (calculate_perturbative_soc) ) then
          call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, n_k_points_task, "KS_eigenvector_complex" )
       else
          call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, 1, "KS_eigenvector_complex" )
       end if




       i_k = 0
       do i_k_point = 1, n_k_points

! SVL add Fock matrix
          if(use_periodic_hf)then

             fock_matrix_complex = (0d0,0d0)
             k(:) = band_k_frac(i_k_point, i_band)

             ! Calculate q-dependent Coulomb matrix for k-q points
             ! n_q_points is equal to n_k_points_old in current implementation
             do i_q_point = 1, n_q_points
                k_minus_q_point_list(i_q_point,:) = k(:) - k_point_list(i_q_point,:)
             enddo
             call get_coulomb_matr_recip(coulomb_matr_recip,0)
             call evaluate_exchange_matr_kspace_single_kpoint_p0 &
                  (i_k_point,KS_eigenvector_old,KS_eigenvector_complex_old,&
                  real_eigenvectors_old,occ_numbers_old,k_weights_old,&
                  lvl_tricoeff_recip1_new,lvl_tricoeff_recip2_new,fock_matrix_complex)
          endif

          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then

             i_k = i_k + 1

             call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                  hamiltonian_w, overlap_matrix_w, &
                  hamiltonian_w_complex, overlap_matrix_w_complex, i_k_point, work_ham, work_ovl)

             if(use_periodic_hf)&
                  hamiltonian_w_complex = hamiltonian_w_complex - hybrid_coeff*fock_matrix_complex

             call clear_lapack_overlap_matrix()

             ! solve KS-equations for both spins
             do i_spin = 1, n_spin, 1

                flag_KS_k_points(1) = BASIS_SINGULAR_NOT_TESTED
                ! calculate the eigenvalues and eigenvectors

                output_priority_old = output_priority
                output_priority = OL_high

                if (flag_rel == 1 .or. out_eigenvec) then
                   i_k_index = i_k
                else
                   i_k_index = 1
                end if

                if (use_elsi .and. .not. use_scalapack) then
                   call solve_KS_elsi_serial(hamiltonian_w_complex(:,i_spin),&
                        overlap_matrix_w_complex,&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector_complex(:,:,i_spin,i_k_index),i_spin,1)
                else
                   call improve_complex_eigenfunctions(&
                        overlap_matrix_w_complex,&
                        hamiltonian_w_complex(:,i_spin),&
                        KS_eigenvalue(:,i_spin,i_k_point),&
                        KS_eigenvector_complex(:,:,i_spin,i_k_index),1)
                end if

                output_priority = output_priority_old

             enddo
          else
             KS_eigenvalue(:,:,i_k_point) = 0.0d0
          end if
       end do

       call sync_eigenvalues( KS_eigenvalue)

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

       ! First, we output the scalar-relativistic band structure

       if(myid==0)then

          ! For a scalar-relativistic calculation, the scalar relativistic band
          ! structures will have file names of form band####.out
          ! For a spin-orbit-coupled calculation, the scalar relativistic band
          ! structure will have file names of form band####.out.no_soc, to
          ! distinguish them from the SOC band structure
          if (calculate_perturbative_soc) then
            sr_suffix = ".out.no_soc"
          else
            sr_suffix = ".out"
          end if

          if(i_band < 10)then
             write(file_name, '(A7,I1,A11)') 'band100', i_band, adjustl(sr_suffix)
          else if(i_band < 100)then
             write(file_name, '(A6,I2,A11)') 'band10', i_band, adjustl(sr_suffix)
          else if(i_band < 1000)then
             write(file_name, '(A5,I3,A11)') 'band1', i_band, adjustl(sr_suffix)
          else
             write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
             stop
          end if


          open(88, file=file_name)

          if(n_spin >1)then

             if(i_band < 10)then
                write(file_name2, '(A7,I1,A11)') 'band200', i_band, adjustl(sr_suffix)
             else if(i_band < 100)then
                write(file_name2, '(A6,I2,A11)') 'band20', i_band, adjustl(sr_suffix)
             else if(i_band < 1000)then
                write(file_name2, '(A5,I3,A11)') 'band2', i_band, adjustl(sr_suffix)
             else
                write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
                stop
             end if
             open(89, file=file_name2)

          end if

          do  i_k_point = 1,  n_k_points

             k(:) = band_k_frac(i_k_point, i_band)

             write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

             do  i_state = 1,  n_states

                write(88,'(F12.5,F15.5)',ADVANCE='NO') occ_numbers(i_state,1,i_k_point), &! *n_k_points, &
                     (KS_eigenvalue(i_state,1,i_k_point)-chemical_potential)* hartree
             end do
             write(88,'()')

             if(n_spin ==2)then

                write(89,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

                do  i_state = 1,  n_states

                   write(89,'(F12.5,F15.5)',ADVANCE='NO') occ_numbers(i_state,2,i_k_point), &!*n_k_points, &
                        (KS_eigenvalue(i_state,2,i_k_point)-chemical_potential)* hartree

                end do
                write(89,'()')
             end if

          end do

          close(88)

          if(n_spin==2) close(89)

       end if

       lowest_un_occ = 1d100
       highest_occ    = -1d100

       do  i_spin = 1,  n_spin
          do  i_k_point = 1,  n_k_points

             do  i_state = 1,  n_states

                if( KS_eigenvalue(i_state,i_spin,i_k_point) >  chemical_potential)then
                   lowest_un_occ = min(lowest_un_occ,  KS_eigenvalue(i_state,i_spin,i_k_point))
                else
                   highest_occ = max(highest_occ, KS_eigenvalue(i_state,i_spin,i_k_point))
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

       if (calculate_perturbative_soc) then
         write(info_str,'(2X,A,I4)')  'Scalar-relativistic "band gap" along reciprocal space direction number: ',i_band
       else
         write(info_str,'(2X,A,I4)')  '"Band gap" along reciprocal space direction number: ',i_band
       end if
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F16.8,A)') '| Lowest unoccupied state:',lowest_un_occ* hartree, ' eV'
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F16.8,A)') '| Highest occupied state :',highest_occ* hartree,  ' eV'
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F16.8,A)') '| Energy difference      :', (lowest_un_occ- highest_occ) * hartree,  ' eV'
       call localorb_info ( info_str )

       write(info_str,'(A)') ''
       call localorb_info ( info_str )

       if (out_eigenvec) then

         write(info_str,'(2X,A,I4,A)') 'Kohn-Sham eigenvectors for band number', i_band, ' will be written to files.'
         call localorb_info ( info_str )

       ! For each k-point
         ! Pull KS_eigenvector_complex to temporary array on thread 0
         ! For each spin channel
           !
           ! Create KS_eigenvector file name
           ! Only on thread 0, call subroutine that writes into the desired file:
           !
           !                   eigenvalue_1  eigenvalue_2  eigenvalue_3  ... eigenvalue_n
           !                   occ_number_1  occ_number_2  occ_number_3  ... occ_number_n
           !
           !  basis_fn_info 1  state_1       state_2       state_3       ... state_n
           !  basis_fn_info 2  state_1       state_2       state_3       ... state_n
           !  etc.

       ! Note: In principle, every task that has a k-point could write, without any temporary
       !       synchronization to task 0.
       !       In the lapack version (present routine), I am preventing this, simply to make sure that we do not
       !       suddenly have 200 writing tasks. In the scalapack version, we should
       !       probably have every "mother" task of each k-point write its own information, period.

         i_k = 0
         do i_k_point = 1,  n_k_points

         ! synchronize the copied Kohn-Sham eigenvector such that a
         ! copy ends up (among others) on the output task, myid=0
           if (myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
           ! k-point stored on current mpi task
             i_k = i_k + 1
             KS_eigenvector_tmp(:,:,:) = KS_eigenvector_complex(:,:,:,i_k)
           else
           ! zero temp. KS eigenvector on all other threads, prior to allreduce
             KS_eigenvector_tmp(:,:,:) = 0.d0
           end if
           call sync_eigenvector_complex(KS_eigenvector_tmp)

         ! The rest is the output operation, which involves taks 0 only
           if (myid.eq.0) then

             k(:) = band_k_frac(i_k_point, i_band)

             do i_spin = 1, n_spin, 1

             call output_complex_eigenvector &
             ( KS_eigenvector_tmp(:,:,i_spin), KS_eigenvalue(:,i_spin,i_k_point), chemical_potential, &
               occ_numbers(:,i_spin,i_k_point), i_band, i_k_point, i_spin, k(1), k(2), k(3) )

             enddo ! end spin loop

           end if ! end restriction of operations to task number 0

         enddo ! end loop over k-points in current band

       end if

       if (calculate_perturbative_soc) then

         ! Because we will never store the SOC-perturbed eigenvectors, here we
         ! compute and save all eigenvalues in the second-variational window,
         ! since they're cheap
         if (allocated(KS_eigenvalue_soc_perturbed)) &
              call aims_deallocate(KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed")
         call aims_allocate( KS_eigenvalue_soc_perturbed, n_states_soc, 1, n_k_points, "KS_eigenvalue_soc_perturbed" )
         KS_eigenvalue_soc_perturbed = 0.0d0

         if (allocated(occ_numbers_soc)) &
              call aims_deallocate(occ_numbers_soc, "occ_numbers_soc")
         call aims_allocate( occ_numbers_soc, n_states_soc, 1, n_k_points, "occ_numbers_soc")

         if (allocated(my_k_points)) deallocate(my_k_points)
         allocate( my_k_points(n_k_points_task),stat=info )
         call check_allocation(info, 'my_k_points                   ')

         ! Undo the round-robin allocation
         i = 1
         do i_k_point = 1, n_k_points, 1
           if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
             my_k_points(i) = i_k_point
             i = i + 1
           end if
         end do

         ! Calculate the second-variational step on each k-point.  This block is the
         ! core functionality of calculate_second_variational_soc,
         ! with the exception of no outputting and lack of determination of
         ! Fermi level, as this has already been set in that main SOC routine.
         do  i_k_point = 1,  n_k_points_task
           this_k_point = my_k_points(i_k_point)

           call construct_SOC_Hamiltonian(n_hamiltonian_matrix_size, soc_matrix, &
                n_basis, n_states, KS_eigenvector(:,:,:,1), KS_eigenvector_complex(:,:,:,i_k_point), this_k_point, &
                n_states_soc, n_states_soc, SOC_Hamiltonian)

           call perform_soc_perturbation( n_states_soc, n_states_soc, SOC_Hamiltonian, KS_eigenvalue(1,1,this_k_point),&
                KS_eigenvalue_soc_perturbed(1, 1, this_k_point), dummy, &
                n_states_soc, n_states_soc, eigenvec_soc_wf_basis )
         end do
         if (use_scalapack .and. my_scalapack_id .ne. 0) then
           KS_eigenvalue_soc_perturbed = 0.0d0
         end if
         call sync_vector(KS_eigenvalue_soc_perturbed, n_states_soc*n_k_points, mpi_comm_global)

         ! Get occupation numbers
         call convert_sr_to_soc_environment ()
         call check_norm_p0(chemical_potential_soc, KS_eigenvalue_soc_perturbed, n_electrons, occ_numbers_soc, dummy, &
              dummy_int)
         call revert_soc_to_sr_environment ()

         if(myid==0)then
            ! No spin quantum number means only one output file, but using old
            ! spin-dependent nomenclature to allow usage of aims plotting scripts
            if(i_band < 10)then
               write(file_name, '(A7,I1,A4)') 'band100', i_band,'.out'
            else if(i_band < 100)then
               write(file_name, '(A6,I2,A4)') 'band10', i_band,'.out'
            else if(i_band < 1000)then
               write(file_name, '(A5,I3,A4)') 'band1', i_band,'.out'
            else
               write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
               stop
            end if

            open(88, file=file_name)

            do  i_k_point = 1,  n_k_points
              ! Having calculated
              k(:) = band_k_frac(i_k_point, i_band)

              write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

              do  i_state = 1,  n_states_soc

                write(88,'(F12.5,F15.5)',ADVANCE='NO') occ_numbers_soc(i_state,1,i_k_point), &! *n_k_points, &
                     (KS_eigenvalue_soc_perturbed(i_state,1,i_k_point)-chemical_potential_soc)* hartree
              end do
              write(88,'()')

            end do

            close(88)

         end if

         lowest_un_occ = 1d100
         highest_occ    = -1d100

         do  i_k_point = 1,  n_k_points

           do  i_state = 1,  n_states_soc

             if( KS_eigenvalue_soc_perturbed(i_state,1,i_k_point) >  chemical_potential_soc)then
                 lowest_un_occ = min(lowest_un_occ,  KS_eigenvalue_soc_perturbed(i_state, 1, i_k_point))
             else
                 highest_occ = max(highest_occ, KS_eigenvalue_soc_perturbed(i_state, 1, i_k_point))
             end if

           end do
         end do

         if(lowest_un_occ < lowest_un_occ_whole_system_soc) then
           lowest_un_occ_whole_system_soc = lowest_un_occ
         end if

         if(highest_occ > highest_occ_whole_system_soc)      highest_occ_whole_system_soc  = highest_occ


         write(info_str,'()')
         call localorb_info(info_str, use_unit,'(A)')

         write(info_str,'(2X,A,I4)')  'Spin-orbit-coupled "band gap" along reciprocal space direction number: ',i_band
         call localorb_info ( info_str )

         write(info_str,'(2X,A,1X,F16.8,A)') '| Lowest unoccupied state:',lowest_un_occ* hartree, ' eV'
         call localorb_info ( info_str )

         write(info_str,'(2X,A,1X,F16.8,A)') '| Highest occupied state :',highest_occ* hartree,  ' eV'
         call localorb_info ( info_str )

         write(info_str,'(2X,A,1X,F16.8,A)') '| Energy difference      :', (lowest_un_occ- highest_occ) * hartree,  ' eV'
         call localorb_info ( info_str )

         write(info_str,'(A)') ''
         call localorb_info ( info_str )

         if (out_eigenvec) then

           write(info_str,'(A)') 'Not writing out SOC-perturbed eigenvectors for band structures.'
           call localorb_info ( info_str )

         end if

      end if ! calculate_perturbative_soc

       if (use_periodic_hf) call deallocate_tb_auxmat()
       if(allocated(lvl_tricoeff_recip1_new)) deallocate(lvl_tricoeff_recip1_new)
       if(allocated(lvl_tricoeff_recip2_new)) deallocate(lvl_tricoeff_recip2_new)
       deallocate (occ_numbers)
       deallocate( KS_eigenvalue)
       call aims_deallocate( KS_eigenvector_complex, "KS_eigenvector_complex" )

    end do ! i_band

    if (allocated(coulomb_matr_recip)) deallocate(coulomb_matr_recip)

    if (calculate_perturbative_soc) then
      write(info_str,'(2X,A)')  'Scalar-relativistic "band gap" of total set of bands: '
    else
      write(info_str,'(2X,A)')  '"Band gap" of total set of bands: '
    end if
    call localorb_info ( info_str )

    write(info_str,'(2X,A,1X,F16.8,A)') '| Lowest unoccupied state:', &
       lowest_un_occ_whole_system* hartree, ' eV'
    call localorb_info ( info_str )

    write(info_str,'(2X,A,1X,F16.8,A)') '| Highest occupied state :', &
       highest_occ_whole_system* hartree,  ' eV'
    call localorb_info ( info_str )

    write(info_str,'(2X,A,1X,F16.8,A)') '| Energy difference      :', &
      (lowest_un_occ_whole_system- highest_occ_whole_system) * hartree,' eV'
    call localorb_info ( info_str )

    ! If requested, output the spin-orbit-coupled band gap predicted by the band
    ! structure
    if (calculate_perturbative_soc) then
      write(info_str,'(2X,A)')
      call localorb_info ( info_str )
      write(info_str,'(2X,A)')  'Spin-orbit-coupled "band gap" of total set of bands: '
      call localorb_info ( info_str )

      write(info_str,'(2X,A,1X,F16.8,A)') '| Lowest unoccupied state:', &
         lowest_un_occ_whole_system_soc* hartree, ' eV'
      call localorb_info ( info_str )

      write(info_str,'(2X,A,1X,F16.8,A)') '| Highest occupied state :', &
        highest_occ_whole_system_soc* hartree,  ' eV'
      call localorb_info ( info_str )

      write(info_str,'(2X,A,1X,F16.8,A)') '| Energy difference      :', &
        (lowest_un_occ_whole_system_soc- highest_occ_whole_system_soc) * hartree,' eV'
      call localorb_info ( info_str )
    end if


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
    if (calculate_perturbative_soc) then
      if (allocated(KS_eigenvalue_soc_perturbed)) &
           call aims_deallocate(KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed")
      if (allocated(occ_numbers_soc)) &
           call aims_deallocate(occ_numbers_soc, "occ_numbers_soc")
    end if

    if (flag_rel.eq.1) then

       call deallocate_scaled_zora_transform
    end if

    call get_times(band_time, clock_band_time, tot_time_band_dos, tot_clock_time_band_dos, .true.)

    write(info_str,'(A)') ''
    call localorb_info(info_str)
    write(info_str,'(2X,A)') &
         "Band Structure                                          :  max(cpu_time) wall_clock(cpu1)"
    call localorb_info(info_str)
    write(info_str, "(2X, A,F15.3,F17.3)")  "| Total Time                                            :", &
         band_time, clock_band_time
    call localorb_info(info_str)
    write(info_str,'(A)') "------------------------------------------------------------"
    call localorb_info(info_str)

  end subroutine out_plot_band_hf_k_space
!******

