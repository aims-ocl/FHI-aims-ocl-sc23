!****s* FHI-aims/out_plot_band_mulliken_scalapack
!  NAME
!    out_plot_band_mulliken_scalapack
!  SYNOPSIS

subroutine out_plot_band_mulliken_scalapack ( )
!  USES
    use physics,               only : occ_numbers, occ_numbers_soc, &
                                      KS_eigenvalue, chemical_potential, &
                                      KS_eigenvalue_soc_perturbed, &
                                      KS_eigenvector, KS_eigenvector_complex, &
                                      KS_eigenvector_soc_perturbed, &
                                      rho, rho_gradient, kinetic_density, &
                                      hamiltonian, overlap_matrix, partition_tab, &
                                      en_vdw, en_pot_vdw, en_xc, en_pot_xc, &
                                      hartree_potential, n_electrons, &
                                      chemical_potential_soc, chemical_potential_spin
    use species_data,          only : l_shell_max, species_name, &
                                      species_pseudoized ! add
    use dimensions,            only : calculate_perturbative_soc, n_atoms, &
                                      n_basis, n_states, n_spin, n_k_points, &
                                      n_hamiltonian_matrix_size, use_periodic_hf, &
                                      n_plot_band, use_lc_wpbeh, l_wave_max !
!add l_wave_max
    use constants,             only : bohr, hartree, pi
    use geometry,              only : coords, lattice_vector, species
    use pbc_lists,             only : k_point_list, k_weights, k_phase, cell_index_bvk,&
                                      n_cells, n_cells_bvk, cell_index, band_k_frac
    use mpi_tasks,             only : myid
    use runtime_choices,       only : collect_eigenvectors, fixed_spin_moment, &
                                      flag_rel, hybrid_coeff, out_eigenvec, &
                                      out_hamiltonian, out_eigenvec, out_overlap, &
                                      real_eigenvectors, use_load_balancing, &
                                      use_local_index, n_points_in_band, &
                                      band_mulliken_orbit_num, & ! add cl
                                      fixed_spin_moment_electrons, use_elsi
    use localorb_io,           only : use_unit, localorb_info
    use scaled_zora_transform, only : allocate_scaled_zora_transform, &
                                      integrate_scaled_zora_transf_p2, &
                                      evaluate_scaled_zora_transf_p1, &
                                      deallocate_scaled_zora_transform
    use scalapack_wrapper,     only : ham_complex, eigenvec, eigenvec_complex, &
                                      mxld, mxcol, l_row, l_col, my_scalapack_comm_all, &
                                      my_scalapack_id, my_k_point, &
                                      construct_overlap_scalapack, &
                                      construct_hamiltonian_scalapack, &
                                      set_sparse_local_ham_scalapack, &
                                      set_sparse_local_ovlp_scalapack, &
                                      construct_mulliken_decomp_scalapack, &
                                      save_overlap_scalapack, &
                                      solve_evp_scalapack_complex
    use synchronize_mpi,       only : sync_eigenvalues
    use synchronize_mpi_basic, only : sync_vector, sync_vector_complex
    use hartree_fock_p0,       only : hf_exchange_matr_complex, &
                                      hf_exchange_matr_complex_SR
    use soc_utilities,         only : perform_soc_perturbation, &
                                      convert_sr_to_soc_environment, &
                                      revert_soc_to_sr_environment
    use load_balancing,        only : use_batch_permutation, batch_perm, n_bp_integ
    use scalapack_soc,         only : mxld_soc, mxcol_soc
    use timing,                only : tot_clock_time_band_dos, tot_time_band_dos, &
                                      get_times, get_timestamps
    use aims_memory_tracking,  only : aims_allocate, aims_deallocate
    use hdf5_output,           only : output_overlap_scalapack, &
                                      output_complex_EV_scalapack, &
                                      output_complex_KS_OCC_scalapack
    use dimensions_soc,        only : n_states_soc
    use ks_wrapper,            only : solve_KS_elsi_parallel
    use elsi_wrapper,          only : aims_elsi_reinit_scf
!  PURPOSE
!   The subroutine output mulliken charge analysis of any k point in plots band
!   structure.
!   Developed by Chi Liu (garnett.liu@duke.edu)
!
    implicit none

!  SOURCE

    complex*16,dimension(:,:), allocatable:: k_phase_band
    complex*16, dimension(:), allocatable :: k_phase_exx_old
    complex*16, dimension(:), allocatable :: k_phase_exx_new
    complex*16, dimension(:,:,:), allocatable :: fock_matr_elem
    complex*16, dimension(:), allocatable :: fock_matr_elem_new_k
    complex*16, dimension(:,:,:), allocatable :: fock_matr_elem_SR
    complex*16, dimension(:), allocatable :: fock_matr_elem_new_k_SR


    integer:: i_band,    i_cell_n, i_spin, i_periodic, i_atom
    integer:: i_cell_1,  i_cell_2, i_cell_3, i_counter
    integer:: i_k_point, i_state, j_state, n_k_points_temp, i
    integer:: n_k_points_band, my_k_point_band
    integer :: i_k_point_old, i_index, i_k_old, i_basis_1, i_basis_2
    integer :: i_cell_bvk

    real*8:: i_x, i_y, i_z
    real*8:: T1, T2, T3
    real*8:: temp_chemical_potential
    real*8:: lowest_un_occ,highest_occ
    real*8:: diff_electrons
    real*8, dimension(:), allocatable :: k_weights_old
    real*8 :: k(3)


    character*30:: file_name, file_name2, sr_suffix
    character*100 :: info_str


    logical :: info

    ! Variables for periodic SOC
    ! Note:  This subroutine is tricky to adapt to SOC, because it throws away
    ! information about eigenvectors, whereas SOC requires the eigenvectors to
    ! generate SOC_Hamiltonian.  If some bug shows up down the line, the LAPACK
    ! variant should be more stable.
    complex*16, dimension(:,:), allocatable :: soc_ham
    real*8,dimension(:,:), allocatable :: soc_matrix
    real*8 :: dummy
    real*8, dimension(:,:), allocatable :: KS_eigenvalue_temp
    complex*16, dimension(:,:), allocatable :: eigenvec_soc_wf_basis
    integer :: dummy_int, this_k_point, info_int
    integer :: ld_soc_matrix
    ! my_k_points is an array for converting between local k-point
    ! indexing (i.e. the k-point indexing of KS_eigenvector) and
    ! global/shared k-point indexing (i.e. the k-point indexing of
    ! KS_eigenvalue)
    real*8:: lowest_un_occ_whole_system, highest_occ_whole_system

    real*8 :: band_time = 0.d0
    real*8 :: clock_band_time = 0.d0

    ! Variables for local indexing and load balancing
    real*8, dimension(:,:), allocatable :: my_hamiltonian_matrix
    real*8, dimension(:),   allocatable :: ovlp_work

! Variables for Mulliken
  real*8, dimension(:, :, :, :, :), allocatable :: mulliken_decomp
  character(LEN=80) :: uuid_str
  integer :: i_size, i_l
  complex*16:: mul_temp
  character*40 :: filename
  ! darn, need an extra aux array of characters
  character*3,dimension(0:l_wave_max) :: l_channel
  integer :: homo_id
  integer :: num_orbit_output_half

    call get_timestamps(band_time, clock_band_time)
     homo_id = int(n_electrons) / 2
     num_orbit_output_half = band_mulliken_orbit_num / 2
     do i_l = 0, l_wave_max, 1
        write(l_channel(i_l),'(A,I1)') "l=",i_l
     enddo

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "Writing the requested band structure output (ScaLAPACK version):"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

! skip real_eigenvectors
! skip use_local_index, use_load_balancing

    ! Integrate the real-space Hamiltonian.
!    call aims_allocate( my_hamiltonian_matrix, 1,1, "my_hamiltonian_matrix" )
!    call aims_allocate( ovlp_work, 1,                           "ovlp_work" )
    call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, &
           partition_tab, l_shell_max, &
           en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

! skip collect_eigenvectors, but since we need eigenvectors in mulliken
! analysis, so allocate KS_eigenvector_complex in i_band loop !

! skip SOC

    lowest_un_occ_whole_system = 1d100
    highest_occ_whole_system    = -1d100

    if(use_periodic_hf)then
       call aims_allocate( k_phase_exx_old, n_cells_bvk,                                              "k_phase_exx_old" )
       do i_cell_bvk = 1, n_cells_bvk
          k_phase_exx_old(i_cell_bvk) = &
               exp(-(0d0,2d0)*pi*sum(k_point_list(my_k_point,:)*cell_index_bvk(i_cell_bvk,:)))
       enddo
       call aims_allocate( k_phase_exx_new, n_cells_bvk,                                              "k_phase_exx_new" )
       call aims_allocate( fock_matr_elem, n_cells_bvk,n_basis,n_spin,                                 "fock_matr_elem" )
       call aims_allocate( fock_matr_elem_new_k, n_spin,                                         "fock_matr_elem_new_k" )
       if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
         call aims_allocate( fock_matr_elem_SR, n_cells_bvk,n_basis,n_spin,                         "fock_matr_elem_SR" )
         call aims_allocate( fock_matr_elem_new_k_SR, n_spin,                                 "fock_matr_elem_new_k_SR" )
       end if
       call aims_allocate( k_weights_old, n_k_points,                                                   "k_weights_old" )
       k_weights_old = k_weights
    endif

    do i_band = 1, n_plot_band
       n_k_points_band =  n_points_in_band(i_band)

       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,I4,A,I4,A)') "Treating all k-points in band plot segment #", i_band, ", one-by-one:"
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')

       ! Need new k-weights later for the occupation numbers

       if (allocated(k_weights))        deallocate(k_weights)
       allocate(k_weights(n_k_points_band))
       k_weights = 1.0d0 /n_k_points_band

       if (allocated(k_phase_band))     call aims_deallocate(k_phase_band,                               "k_phase_band" )
       call aims_allocate( k_phase_band, n_cells,n_k_points_band,                                        "k_phase_band" )

       k_phase_band = 1.0

       do  i_k_point = 1,  n_k_points_band
          k(:) = band_k_frac(i_k_point, i_band)
          i_x = k(1)
          i_y = k(2)
          i_z = k(3)

          if(i_x < 0)then
             i_x = 1 + i_x
          end if
          if(i_y < 0)then
             i_y = 1 + i_y
          end if
          if(i_z < 0)then
             i_z = 1 + i_z
          end if


          do i_cell_n = 1, n_cells, 1

             i_cell_1 = cell_index(i_cell_n,1)
             i_cell_2 = cell_index(i_cell_n,2)
             i_cell_3 = cell_index(i_cell_n,3)


             T1 = 2*pi*(i_x) ! / (n_k_points_band_xyz(1))
             T2 = 2*pi*(i_y) !/ (n_k_points_band_xyz(2))
             T3 = 2*pi*(i_z) !/ (n_k_points_band_xyz(3))

             k_phase_band( i_cell_n, i_k_point) &
                  = ((exp((0,1)*T1*i_cell_1))* (exp((0,1)*T2*i_cell_2)) * (exp((0,1)*T3*i_cell_3)))

          end do
       end do ! i_k_point

       if(abs(sum(k_weights)-1) >1e-10)then
          write(use_unit,*) 'Error: sum of k-vector weights is not one!', sum(k_weights)
          stop
       end if

       ! Allocate those arrays that change with the bands

       if (allocated(occ_numbers))            deallocate(occ_numbers)
       if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
       allocate (occ_numbers(n_states,n_spin,n_k_points_band))
       allocate (KS_eigenvalue(n_states,n_spin,n_k_points_band))
! add
      if (allocated( KS_eigenvector_complex )) call aims_deallocate(KS_eigenvector_complex,   "KS_eigenvector_complex" )
      call aims_allocate( KS_eigenvector_complex ,n_basis, n_states, n_spin, 1, "KS_eigenvector_complex" )
!      call aims_allocate( KS_eigenvector_complex ,n_basis, n_states, n_spin, n_k_points_band,  "KS_eigenvector_complex" )
!      if (allocated( eigenvec_complex )) call aims_deallocate( eigenvec_complex, "eigenvec_complex" )
!      call aims_allocate( eigenvec_complex, mxld, mxcol, n_spin, "eigenvec_complex" )

      if (allocated(mulliken_decomp)) call aims_deallocate(mulliken_decomp, "mulliken_decomp" )
      call aims_allocate(mulliken_decomp, 0,l_wave_max, 1,n_atoms, 1,n_states, 1,n_spin, 1,n_k_points_band, "+mulliken_decomp")
!      if (allocated (mulliken_decomp)) deallocate(mulliken_decomp)
!      allocate(mulliken_decomp( 0:l_wave_max, n_atoms, n_states, n_spin,n_k_points_band ))

!      write(info_str,'(2X,A)')  'after allocate mulliken_decomp'
!      call localorb_info(info_str)
      mulliken_decomp = 0.d0

! skip: Here we need to allocate the SOC-specific arrays

       ! Set KS_eigenvalue and occ_numbers to zero. Otherwise they are filled with random numbers
       ! at indices which are not written for a given (partial)  eigenwert solution of lower dimension n_k_points.
       ! This  may disturb sync.
       KS_eigenvalue=0.0
       occ_numbers=0.0
       KS_eigenvector_complex = 0.0
!       eigenvec_complex = 0.0

       ! Solve the eigenwert problem

       ! The parallelization for this section is as follows.  As we are using ScaLAPACK,
       ! each CPU has been assigned a single k-point that it calculates during
       ! the SCF cycle, the index of which is stored in the global variable my_k_point,
       ! and various communicators exist between all CPUs assigned to the same k-point.
       ! In order to use ScaLAPACK routines and not have to reinitialize the communicator
       ! infrastructure, we want to respect this original assignment, but we now
       ! have a different number of k-points (n_k_points_band) than what the original
       ! infrastructure was constructed for (n_k_points).  Thus we process the
       ! k-points of the band in batches of size n_k_points.  Because each batch
       ! has the same number of k-points as the original k-grid, which is what
       ! ScaLAPACK was initialized for, we can now use the original infrastructure
       ! within the current batch to solve the KS equations for all k-points within this
       ! batch.  We then move onto the next batch, and continue doing so until
       ! we've process all k-points in the band.  (For the last batch, which
       ! in general will not contain n_k_points unique k-points, we simply redo
       ! one k-point continuously on the excess CPUs then zero it out our effort
       ! on all but one CPU.)
       do i_k_point = 1, n_k_points_band, n_k_points

           write(info_str,'(2X,A,I4,A,I4,A)') "Calculating k-point ", i_k_point, " in band plot segment #", i_band, ":"
           call localorb_info(info_str,use_unit,'(A)')

           ! One iteration of this loop solves the KS equations for
           ! the k-points:   i_k_point ... i_k_point+n_k_points-1

           ! Set k-phase - if there are less k-points left than
           ! n_k_points, calculate last one repeatedly
           do i = 1, n_k_points
              k_phase(:,i) = k_phase_band(:,MIN(i_k_point+i-1,n_k_points_band))
           end do

           ! My k-point for this iteration
           ! Here my_k_point is the k-point assigned by the original ScaLAPACK
           ! infrastructure for the original scf cycle, and i_k_point is the
           ! offset for the beginning of the current batch.  So my_k_point_band
           ! functions as the k-point index within the band, and my_k_point as the
           ! k-point index within the batch.
           my_k_point_band = MIN (i_k_point+my_k_point-1, n_k_points_band)

           if(use_periodic_hf)then
              k(:) = band_k_frac(my_k_point_band, i_band)
              do i_cell_bvk = 1, n_cells_bvk
                 k_phase_exx_new(i_cell_bvk) = exp((0d0,2d0)*pi*sum(k(:)*cell_index_bvk(i_cell_bvk,:)))
              enddo
           endif

              ! Construct the Hamiltonian and overlap matrices in BLACS storage format
! skip use_local_index, use_load_balancing, out_overlap
           call construct_overlap_scalapack(overlap_matrix)
           call construct_hamiltonian_scalapack(hamiltonian)
! SAVE here !!!!!!!!
           call save_overlap_scalapack

              if(use_periodic_hf)then

                 do i_basis_1 = 1, n_basis
                    fock_matr_elem = (0d0,0d0)
                    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                       fock_matr_elem_SR = (0d0,0d0)
                    end if
                    if(l_row(i_basis_1)>0)then
                       do i_basis_2 = 1, n_basis
                          if(l_col(i_basis_2)>0)then
                             do i_spin = 1, n_spin
                                fock_matr_elem(:,i_basis_2,i_spin) = &
                                     hf_exchange_matr_complex(l_row(i_basis_1),l_col(i_basis_2),1,i_spin)*&
                                     k_phase_exx_old(:)*k_weights_old(my_k_point)
                                if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                   fock_matr_elem_SR(:,i_basis_2,i_spin) = &
                                        hf_exchange_matr_complex_SR(l_row(i_basis_1),l_col(i_basis_2),1,i_spin)*&
                                        k_phase_exx_old(:)*k_weights_old(my_k_point)
                                end if
                             enddo
                          endif
                       enddo
                    endif
                    call sync_vector_complex(fock_matr_elem,n_basis*n_cells_bvk*n_spin)
                    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                       call sync_vector_complex(fock_matr_elem_SR,n_basis*n_cells_bvk*n_spin)
                    end if
                    !debug
!                    do i_cell_1 = 1, n_cells_bvk
!                       do i_basis_2 = 1, n_basis
!                          do i_spin = 1, n_spin
!                             if(myid.eq.0)&
!                                  write(use_unit,*) "fmeb", cell_index_bvk(i_cell_1,:), i_basis_1, i_basis_2, fock_matr_elem(i_cell_1,i_basis_2,i_spin)
!                          enddo
!                       enddo
!                    enddo
                    !enddebug
                    if(l_row(i_basis_1)>0)then
                       do i_basis_2 = 1, n_basis
                          if(l_col(i_basis_2)>0)then
                             do i_spin = 1, n_spin
                               fock_matr_elem_new_k(i_spin) = sum(fock_matr_elem(:,i_basis_2,i_spin)*k_phase_exx_new(:))
                               if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                  fock_matr_elem_new_k_SR(i_spin) = sum(fock_matr_elem_SR(:,i_basis_2,i_spin)*k_phase_exx_new(:))
                               end if
                             enddo
                             if (use_lc_wpbeh) then
                               if (hybrid_coeff /= 0.d0) then
                                 ham_complex(l_row(i_basis_1),l_col(i_basis_2),:) = &
                                      ham_complex(l_row(i_basis_1),l_col(i_basis_2),:) - &
                                      (fock_matr_elem_new_k(:) + hybrid_coeff*fock_matr_elem_new_k_SR(:))
                               else
                                 ham_complex(l_row(i_basis_1),l_col(i_basis_2),:) = &
                                      ham_complex(l_row(i_basis_1),l_col(i_basis_2),:) - &
                                      fock_matr_elem_new_k(:)
                               end if
                             else
                               ham_complex(l_row(i_basis_1),l_col(i_basis_2),:) = &
                                   ham_complex(l_row(i_basis_1),l_col(i_basis_2),:) - &
                                   hybrid_coeff*fock_matr_elem_new_k(:)
                             end if
!                             write(use_unit,*) "fmenk", myid, k, i_basis_1, i_basis_2, fock_matr_elem_new_k
                          endif
                       enddo
                    endif
                 enddo

              endif! if(use_periodic_hf)

              do i_spin = 1, n_spin
                 if (use_elsi) then
                    ! Reinit for new overlap
                    call aims_elsi_reinit_scf()

                    call solve_KS_elsi_parallel(&
                         KS_eigenvalue(:,i_spin,my_k_point_band),&
                         KS_eigenvector_complex(:,:,i_spin,1),i_spin)
                 else
                    call solve_evp_scalapack_complex(&
                         KS_eigenvalue(:,i_spin,my_k_point_band),&
                         KS_eigenvector_complex(:,:,i_spin,1),i_spin)
                 end if
              end do

              ! If relativistic correction is required, calculate scaled Zora corrections
              if (flag_rel.eq.1) then
                 call allocate_scaled_zora_transform
                 call integrate_scaled_zora_transf_p2( &
                          rho, rho_gradient, kinetic_density, hartree_potential,&
                          partition_tab, l_shell_max )
                 call evaluate_scaled_zora_transf_p1( &
                          KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue)
                 call deallocate_scaled_zora_transform
              end if

              ! If the last point is calculated more than once
              ! throw away the additional eigenvalues
              if (i_k_point+my_k_point-1 > n_k_points_band) then
                  KS_eigenvalue(:,:,my_k_point_band)=0
!                  KS_eigenvector_complex(:,:,:,my_k_point_band) = 0
! skip soc
              end if

! skip out_eigenvec
          write(info_str,'(A)') ""
          call localorb_info ( info_str )

!       call save_overlap_scalapack, not here
       call construct_mulliken_decomp_scalapack(mulliken_decomp(0,1,1,1,my_k_point_band))
       if(myid == 0) print*, 'end construct_mulliken_decomp_scalapack'
! not a problem of sync, since the actual numbers are not correct
!       print*, "myid, mulliken_decomp(0:l_wave_max,1,171,1,1)", myid, mulliken_decomp(0:l_wave_max,1,171,1,1)
       if (i_k_point+my_k_point-1 > n_k_points_band) then
         mulliken_decomp(:, :, :, :, my_k_point_band) = 0
       end if


      end do ! i_k_point
!       if (myid == 0) print*, "after end do i_k_point"

       ! Temporary adapt n_k_points to the given k_point distribution
       n_k_points_temp = n_k_points
       n_k_points      = n_k_points_band

       ! Synchronize eignevalues
       call sync_eigenvalues(KS_eigenvalue)

!       if (myid == 0) print*, "after sync_eigenvalues"

       ! Obtain occupation numbers
       if (.not.fixed_spin_moment) then

         call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)

       else ! if fixed_spin_moment

         do i_spin = 1, n_spin, 1

           call check_norm_periodic_v2(chemical_potential_spin(i_spin), KS_eigenvalue(:,i_spin,:), &
               fixed_spin_moment_electrons(i_spin), occ_numbers(:,i_spin,:), diff_electrons,&
               i_counter,i_spin)

         enddo

       end if

       ! Restore initial k_point information , comment out !!
!       n_k_points = n_k_points_temp


! not a problem of sync, since the actual numbers are not correct
       call sync_vector(mulliken_decomp(0,1,1,1,1),(1+l_wave_max)*n_atoms*n_states*n_spin*n_k_points)
!       if (myid == 0) print*, "after sync_vector mulliken_decomp"

!  keep n_k_points = n_k_points_band
!       n_k_points = n_k_points_temp
! skip output including: out_overlap .or. out_hamiltonian .or. out_eigenvec

       if (myid==0) then
         ! file names for band plotting
! skip calculate_perturbative_soc
! output band.out
            ! file names for band plotting
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

            k(:) = band_k_frac(i_k_point, i_band)

            write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

            do  i_state = 1,  n_states

              write(88,'(F12.5,2X,F15.5)',ADVANCE='NO') occ_numbers(i_state,1,i_k_point), &! *n_k_points, &
                   (KS_eigenvalue(i_state,1,i_k_point)-chemical_potential)* hartree
            end do
            write(88,'()')
          end do
            close(88)
! output bandmlk.out


           sr_suffix = ".out"

         if(i_band < 10)then
            write(file_name, '(A10,I1,A11)') 'bandmlk100', i_band, adjustl(sr_suffix)
         else if(i_band < 100)then
           write(file_name, '(A9,I2,A11)') 'bandmlk10', i_band, adjustl(sr_suffix)
         else if(i_band < 1000)then
           write(file_name, '(A8,I3,A11)') 'bandmlk1', i_band, adjustl(sr_suffix)
         else
           write(use_unit,*) 'Error: automatic file name do not work with more than 999 bands!'
           stop
         end if

         open(88, file=file_name)

! skip n_spin > 1

        !write bands

         do  i_k_point = 1,  n_k_points_band
            k(:) = band_k_frac(i_k_point, i_band)
            i_x = k(1)
            i_y = k(2)
            i_z = k(3)

           write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, i_x, i_y, i_z

           if (out_eigenvec) then
             do i_spin = 1, n_spin, 1
               call output_complex_KS_OCC_scalapack &
                      ( KS_eigenvalue(:,i_spin,i_k_point), chemical_potential, occ_numbers(:,i_spin,i_k_point), &
                        i_band, i_k_point, i_spin)
             end do
           end if

           do  i_state = 1,  n_states
             write(88,'(F12.5,2X,F15.5)',ADVANCE='NO') occ_numbers(i_state,1,i_k_point), &! *n_k_points_band, &
                   (KS_eigenvalue(i_state,1,i_k_point)-chemical_potential)* hartree
           end do

           write(88,'()')

! skip n_spin == 2

         end do

         close(88)

       end if ! myid .eq. 0

       lowest_un_occ = 1d100
       highest_occ    = -1d100

       do  i_spin = 1,  n_spin
          do  i_k_point = 1,  n_k_points_band
             do  i_state = 1,  n_states
                if (  occ_numbers(i_state,i_spin,i_k_point) < 0.5)then
                  lowest_un_occ = min(lowest_un_occ,  KS_eigenvalue(i_state,i_spin,i_k_point))
                else
                  highest_occ = max(highest_occ, KS_eigenvalue(i_state,i_spin,i_k_point))
                end if
             end do
          end do
       end do

       if (calculate_perturbative_soc) then
         write(info_str,'(2X,A,I4)')  'Scalar-relativistic "band gap" along reciprocal space direction number: ',i_band
       else
         write(info_str,'(2X,A,I4)')  '"Band gap" along reciprocal space direction number: ',i_band
       end if
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F16.8,A)')  '|Lowest unoccupied state:',lowest_un_occ* hartree, 'eV'
       call localorb_info ( info_str )

       write(info_str,'(2X,A,1X,F16.8,A)') '|Highest occupied state:  ',highest_occ* hartree,  'eV'
       call localorb_info ( info_str )

       write(info_str,'(A)') ''
       call localorb_info ( info_str )


!
       if (myid.eq.0) then
!           print*, 'n_k_points', n_k_points
           open(50, file=file_name)

           do i_k_point = 1, n_k_points
              k(:) = band_k_frac(i_k_point, i_band)
              write(50,'(A,I5,A,F12.8,1X,F12.8,1X,F12.8,A)') &
                "k point number: ",i_k_point, ": ( ", k(1), k(2), k(3), " )"
              write(50,'(4X,A,7X,A,2X,A,1X,A, 7X,A, 10(9X,A3))') &
                "State", "eigenvalue", "occ.number", "atom", "total"

              do i_state=1,n_states,1
                if (i_state <= homo_id + num_orbit_output_half .and. &
                   i_state > homo_id - num_orbit_output_half) then
                   write(50,'(4X,A,2X,I7)') "State", i_state
                end if

                do i_atom = 1, n_atoms, 1
                  if(species_pseudoized(species(i_atom))) cycle

                  do i_spin = 1, n_spin, 1

                    if (n_spin.gt.1) then
                       write(50,*)
                       if (i_spin.eq.1) then
                          write(50,'(2X,A,A)') "Spin channel: ", "up"
                       else
                          write(50,'(2X,A,A)') "Spin channel: ", "down"
                       end if
                       write(50,*)
                    end if

                    if (i_state <= homo_id + num_orbit_output_half .and. i_state > &
                        homo_id - num_orbit_output_half) then
                      write(50,'(2X,I7,2X,F15.5,2X,F10.7,2X, I3,2X, F10.5,10(2X,F10.5))') &
                          i_state, (KS_eigenvalue(i_state,i_spin,i_k_point)  - chemical_potential)* hartree, &
                          occ_numbers(i_state,i_spin,i_k_point), i_atom, &
                          sum(mulliken_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,i_k_point)), &
                          ( mulliken_decomp(i_l,i_atom,i_state,i_spin,i_k_point), i_l=0,l_shell_max(species(i_atom)) )
                    end if
                  enddo
               end do
            enddo

         enddo

         close(50)
       end if ! myid == 0
! end out_mulliken

       n_k_points = n_k_points_temp
! Skip       ! Now, output the SOC band structure if it's been requested
    end do   ! loop over i_band


    ! Deallocate

   if (allocated (occ_numbers))                 deallocate (occ_numbers)
   if (allocated (KS_eigenvalue))               deallocate (KS_eigenvalue)

   ! Allocatable variables from this subroutine
   if (allocated (k_phase_band))               call aims_deallocate( k_phase_band,                               "k_phase_band" )
   if (allocated (k_phase_exx_old))            call aims_deallocate( k_phase_exx_old,                         "k_phase_exx_old" )
   if (allocated (k_phase_exx_new))            call aims_deallocate( k_phase_exx_new,                         "k_phase_exx_new" )
   if (allocated (fock_matr_elem))             call aims_deallocate( fock_matr_elem,                           "fock_matr_elem" )
   if (allocated (fock_matr_elem_new_k))       call aims_deallocate( fock_matr_elem_new_k,               "fock_matr_elem_new_k" )
   if (allocated (fock_matr_elem_SR))          call aims_deallocate( fock_matr_elem_SR,                     "fock_matr_elem_SR" )
   if (allocated (fock_matr_elem_new_k_SR))    call aims_deallocate( fock_matr_elem_new_k_SR,         "fock_matr_elem_new_k_SR" )
   if (allocated (k_weights_old))              call aims_deallocate( k_weights_old,                             "k_weights_old" )
   if (allocated (soc_ham))                    call aims_deallocate( soc_ham,                                         "soc_ham" )
   if (allocated (soc_matrix))                 call aims_deallocate( soc_matrix,                                   "soc_matrix" )
   if (allocated (KS_eigenvalue_temp))         call aims_deallocate( KS_eigenvalue_temp,                   "KS_eigenvalue_temp" )
   if (allocated (eigenvec_soc_wf_basis))      call aims_deallocate( eigenvec_soc_wf_basis,             "eigenvec_soc_wf_basis" )
   if (allocated (my_hamiltonian_matrix))      call aims_deallocate( my_hamiltonian_matrix,             "my_hamiltonian_matrix" )
   if (allocated (ovlp_work))                  call aims_deallocate( ovlp_work,                                     "ovlp_work" )
   if (allocated (mulliken_decomp))            call aims_deallocate(mulliken_decomp, "mulliken_decomp")

   ! Allocatable variables from other subroutines (which are tracked)
   if (allocated(KS_eigenvector))              call aims_deallocate( KS_eigenvector,                           "KS_eigenvector" )
   if (allocated(KS_eigenvector_complex))      call aims_deallocate( KS_eigenvector_complex,           "KS_eigenvector_complex" )
   if (allocated(KS_eigenvalue_soc_perturbed)) call aims_deallocate( KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed" )
   if (allocated(occ_numbers_soc))             call aims_deallocate( occ_numbers_soc,                         "occ_numbers_soc" )


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

  end subroutine out_plot_band_mulliken_scalapack
!******

