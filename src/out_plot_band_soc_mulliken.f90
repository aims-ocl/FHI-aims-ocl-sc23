!****s* FHI-aims/out_plot_band
!  NAME
!    out_plot_band
!  SYNOPSIS

  subroutine out_plot_band_soc_mulliken ( )
!  USES
    use physics
    use species_data
    use dimensions
    use constants
    use pbc_lists
    use runtime_choices
    use scaled_zora_transform
    use lapack_wrapper
    use hartree_fock_p0
    use prodbas
    use soc_utilities, only : perform_soc_perturbation, &
                              convert_sr_to_soc_environment, &
                              revert_soc_to_sr_environment, &
                              convert_wf_basis_to_compute_basis
    use timing
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use generate_aims_uuid, only: write_aims_uuid
    use hdf5_output, only: output_complex_overlap_matrix, &
                           output_complex_eigenvector, &
                           output_complex_hamiltonian_matrix
    use dimensions_soc, only : n_states_soc, n_saved_states_soc, n_basis_soc, &
                               n_basis_soc_coll, soc_saved_state_start
    use basis, only: basis_l
    use geometry, only: coords, lattice_vector, species
    use synchronize_mpi
    use ks_wrapper, only: solve_KS_elsi_serial
! add n_saved_states_soc, n_basis_soc, convert_wf_basis_to_compute_basis

!   The subroutine plots mulliken charge decomposition of any k point in the
!   band k path.
!   This works only with lapack type of eigenvectors.
!   This routine combines the out_plot_band and output_mulliken(lapack part).
!   Deal with SOC case.
!   Developed by Chi Liu (garnett.liu@duke.edu)
!
    implicit none

    real*8,    dimension(:,:),allocatable :: hamiltonian_w
    real*8,    dimension(:),  allocatable :: overlap_matrix_w
    complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
    complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex

    real*8, dimension(:,:,:),allocatable :: work_ham
    real*8, dimension(:,:),allocatable :: work_ovl

    integer:: i_band,    i_cell_n, i_spin
    integer:: i_cell(3), i_counter
    integer:: i_k_point, i_state, j_state, i_k, i_k_index
    integer :: i_k_point_old, n_k_points_old, i_index, i_k_old, i_basis_1, i_basis_2
    complex*16 :: k_phase_exx_old
    real*8, dimension(:), allocatable :: k_weights_old
    complex*16, dimension(:,:), allocatable :: fock_matrix_complex
    complex*16, dimension(:,:), allocatable :: fock_matrix_complex_SR
    logical :: real_eigenvectors_old

    real*8:: k(3) !, k_lowest_un_occ(3), k_highest_occ(3),  k_lowest_un_occ_whole_system(3), k_highest_occ_whole_system(3)
    real*8:: lowest_un_occ,highest_occ, lowest_un_occ_whole_system, highest_occ_whole_system
    real*8:: lowest_un_occ_whole_system_soc, highest_occ_whole_system_soc
    real*8:: diff_electrons

    ! for possible output of Kohn-Sham eigenvectors
    complex*16, dimension(:,:,:),allocatable :: KS_eigenvector_tmp

    ! file name infrastructure etc.

    character*10 :: num_char
    integer :: output_priority_old
    integer :: i_atom, i_periodic
    character*50 :: file_name, file_name2, sr_suffix
    character*100 :: info_str

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
! add
  real*8 :: i_x, i_y, i_z
  real*8, dimension(:, :, :, :, :), allocatable :: mulliken_decomp
  integer :: i_size, i_l, basis_offset, i_cell_i
  complex*16:: mul_temp
  character*40 :: filename
  character(LEN=80) :: uuid_str
  ! darn, need an extra aux array of characters
  character*3,dimension(0:l_wave_max) :: l_channel
  integer :: homo_id
  integer :: num_orbit_output_half

     homo_id = int(n_electrons)
     num_orbit_output_half = band_mulliken_orbit_num
     do i_l = 0, l_wave_max, 1
        write(l_channel(i_l),'(A,I1)') "l=",i_l
     enddo

! end add
    call get_timestamps(band_time, clock_band_time)

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "Writing the requested band structure output (LAPACK version):"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

    call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, &
         partition_tab, l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

    ! If SOC enabled, set up initial matrices
    if (calculate_perturbative_soc) then
      call aims_allocate( SOC_Hamiltonian, n_states_soc, n_states_soc, "SOC_Hamiltonian" )
      call aims_allocate( soc_matrix, n_hamiltonian_matrix_size,3,          "soc_matrix" )

      call integrate_soc_matrix (rho, hartree_potential, partition_tab, soc_matrix )
   else ! otherwise SOC code will never be touched, set up dummy indices
      call aims_allocate( SOC_Hamiltonian, 1,1, "SOC_Hamiltonian" )
      call aims_allocate( soc_matrix, 1,1,           "soc_matrix" )
      call aims_allocate( my_k_points, 1,           "my_k_points" )
    end if

    call aims_allocate(hamiltonian_w_complex, n_basis*(n_basis+1)/2,n_spin,      "hamiltonian_w_complex")
    call aims_allocate(overlap_matrix_w_complex, n_basis*(n_basis+1)/2,       "overlap_matrix_w_complex")
    if(use_periodic_hf)then
      call aims_allocate(fock_matrix_complex,   n_basis*(n_basis+1)/2,n_spin,      "fock_matrix_complex")
    end if
    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
      call aims_allocate(fock_matrix_complex_SR, n_basis*(n_basis+1)/2, n_spin, "fock_matrix_complex_SR")
    end if

    ! dummy allocations for the real arrays.
    ! These trigger a compiler warning for -check pointers otherwise
    ! as they are included (but never touched!) in construct_hamiltonian_and_ovl below.
    call aims_allocate(hamiltonian_w, 1,1,    "hamiltonian_w" )
    call aims_allocate(overlap_matrix_w,1, "overlap_matrix_w" )

    if(packed_matrix_format == PM_none)then
       call aims_allocate(work_ham, n_centers_basis_I, n_centers_basis_I, n_spin, "work_ham")
       call aims_allocate(work_ovl, n_centers_basis_I, n_centers_basis_I,         "work_ovl")
    else
       ! dummy only, never touched
       call aims_allocate(work_ham, 1, 1, 1, "work_ham")
       call aims_allocate(work_ovl, 1, 1,    "work_ovl")
    end if

    if (out_eigenvec) then
      call aims_allocate(KS_eigenvector_tmp, n_basis,n_states,n_spin,    "KS_eigenvector_tmp")
    end if

    if (flag_rel.eq.1) then
       call allocate_scaled_zora_transform
       call integrate_scaled_zora_transf_p2( &
            rho, rho_gradient, kinetic_density, hartree_potential,    &
            partition_tab, l_shell_max)
    end if
    n_k_points_old = n_k_points ! SVL store for periodic exx bands
    call aims_allocate(k_weights_old, n_k_points_old, "n_k_points_old")
    k_weights_old = k_weights
    real_eigenvectors_old = real_eigenvectors
    lowest_un_occ_whole_system = 1d100
    highest_occ_whole_system    = -1d100
    lowest_un_occ_whole_system_soc = 1d100
    highest_occ_whole_system_soc    = -1d100

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
       k_weights = 1.0d0 / dble(n_k_points)

       deallocate(k_phase)
       allocate(k_phase(n_cells,n_k_points))

       do  i_k_point = 1, n_k_points
          k(:) = band_k_frac(i_k_point, i_band)

          do i_cell_n = 1, n_cells
             k_phase( i_cell_n, i_k_point) = exp((0.d0,2.d0)*pi*sum(k(:)*dble(cell_index(i_cell_n,:))))
          end do
       end do

       real_eigenvectors = .false.

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



       if (allocated(occ_numbers)) deallocate(occ_numbers)
       if (allocated(KS_eigenvalue)) deallocate(KS_eigenvalue)
       if (allocated(KS_eigenvector_complex)) call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex,")
       allocate(occ_numbers(n_states,n_spin,n_k_points))
       allocate(KS_eigenvalue(n_states,n_spin,n_k_points))

       if ( (flag_rel.eq.1) .or. (out_eigenvec) .or. (calculate_perturbative_soc) ) then
          call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, n_k_points_task, "KS_eigenvector_complex" )
       else
          call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, 1, "KS_eigenvector_complex" )
       end if

! add
       if (allocated(KS_eigenvector_soc_perturbed)) deallocate(KS_eigenvector_soc_perturbed)
       allocate(KS_eigenvector_soc_perturbed(n_basis_soc, n_saved_states_soc, 1, n_k_points))
! end add

       if (out_overlap .or. out_hamiltonian) then
       !write geometry when out_overlap or out_hamiltonian has been specified
          if(myid==0)then

             open(88, file='geometry.in_plot_band')
             write(88,'(A)') '# '
             write(88,'(A)') '# This is the geometry file that corresponds to the geometry input'
             write(88,'(A)') '# The geometry is given in its center cell.'
             call write_aims_uuid(info_str)
             write(88,'(A,2X,A)') '#', info_str
             write(88,'(A)') '# '


             do i_periodic = 1, 3
                write(88,'(A,3F16.8)') &
                & 'lattice_vector', lattice_vector(:,i_periodic)*bohr
             end do

             do i_atom = 1, n_atoms
                write(88,'(A,3F16.8,A,A)') &
                & 'atom ', coords(:,i_atom)*bohr, ' ',trim(species_name(species(i_atom)))
             end do

             close(88)
          end if
       end if !(out_overlap)

       i_k = 0
       do i_k_point = 1, n_k_points

! SVL add Fock matrix
          if(use_periodic_hf)then
             fock_matrix_complex = (0d0,0d0)
             if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                fock_matrix_complex_SR = (0d0,0d0)
             end if
             k(:) = band_k_frac(i_k_point, i_band)
             do i_spin = 1, n_spin
                  i_k_old = 0
                  do i_k_point_old = 1, n_k_points_old
                    if (myid.eq.  MOD(i_k_point_old, n_tasks) .and. myid<= n_k_points_old ) then
                       i_k_old = i_k_old + 1
                       do i_cell_n = 1, n_cells_bvk
                          k_phase_exx_old = exp((0d0,2d0)*pi*(sum(k(:)* dble(cell_index_bvk(i_cell_n,:)))-&
                          sum(k_point_list(i_k_point_old,:)* dble(cell_index_bvk(i_cell_n,:)))))

                          i_index = 0
                          do i_basis_2 = 1,n_basis, 1
                           do i_basis_1 = 1,i_basis_2,1
                               i_index = i_index + 1
                               if(real_eigenvectors_old)then
                                  fock_matrix_complex(i_index,i_spin) = fock_matrix_complex(i_index,i_spin)+&
                                       hf_exchange_matr_real(i_basis_1,i_basis_2,i_k_old,i_spin)*&
                                       dble(k_phase_exx_old)*k_weights_old(i_k_point_old)
                                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                     fock_matrix_complex_SR(i_index,i_spin) = fock_matrix_complex_SR(i_index,i_spin)+&
                                          hf_exchange_matr_real_SR(i_basis_1,i_basis_2,i_k_old,i_spin)*&
                                          dble(k_phase_exx_old)*k_weights_old(i_k_point_old)
                                  end if
                               else
                                  fock_matrix_complex(i_index,i_spin) = fock_matrix_complex(i_index,i_spin)+&
                                       hf_exchange_matr_complex(i_basis_1,i_basis_2,i_k_old,i_spin)*&
                                       k_phase_exx_old*k_weights_old(i_k_point_old)
                                  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                                     fock_matrix_complex_SR(i_index,i_spin) = fock_matrix_complex_SR(i_index,i_spin)+&
                                          hf_exchange_matr_complex_SR(i_basis_1,i_basis_2,i_k_old,i_spin)*&
                                          k_phase_exx_old*k_weights_old(i_k_point_old)
                                  end if
                               end if
                           enddo
                          enddo
                       enddo
                    endif
                  end do
             enddo
             call sync_vector_complex(fock_matrix_complex,n_basis*(n_basis+1)/2*n_spin)
             if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                call sync_vector_complex(fock_matrix_complex_SR,n_basis*(n_basis+1)/2*n_spin)
             end if
          endif ! use_periodic_hf

          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then

             i_k = i_k + 1

             call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                  hamiltonian_w, overlap_matrix_w, &
                  hamiltonian_w_complex, overlap_matrix_w_complex, i_k_point, work_ham, work_ovl)

             if(use_periodic_hf) then
                if (use_lc_wpbeh) then
                   if (hybrid_coeff /= 0.d0) then
                      hamiltonian_w_complex = hamiltonian_w_complex - (fock_matrix_complex + hybrid_coeff*fock_matrix_complex_SR)
                   else
                      hamiltonian_w_complex = hamiltonian_w_complex - fock_matrix_complex
                   end if
                else
                   hamiltonian_w_complex = hamiltonian_w_complex - hybrid_coeff*fock_matrix_complex
                end if
             end if

             if (out_overlap) then
                ! This is the place to do it. Since this is Lapack, each k-point lives
                ! on a separate task. All tasks are involved in the writing.

                ! Recreate k after all (for information only)
                k(:) = band_k_frac(i_k_point, i_band)

                call output_complex_overlap_matrix &
                  ( overlap_matrix_w_complex, i_band, i_k_point, k(1), k(2), k(3) )


             end if

             if (out_hamiltonian) then
                ! Recreate k after all (for information only)
                k(:) = band_k_frac(i_k_point, i_band)

                call output_complex_hamiltonian_matrix &
                  ( hamiltonian_w_complex, i_band, i_k_point, k(1), k(2), k(3) )

             end if

             ! This routine clears, deallocates and reallocates an internal (temporary) array
             ! that holds the overlap matrix inside lapack - and that may need to be resized.
             call clear_lapack_overlap_matrix()

             ! solve KS-equations for both spins
             do i_spin = 1, n_spin, 1

                flag_KS_k_points(1) = BASIS_SINGULAR_NOT_TESTED

                ! calculate the eigenvalues and eigenvectors

                output_priority_old = output_priority
                output_priority = OL_high

                if (flag_rel == 1 .or. out_eigenvec .or. calculate_perturbative_soc) then
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
       end do ! i_k_point

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

          ! file names for band plotting

          ! For a scalar-relativistic calculation, the scalar relativistic band
          ! structures will have file names of form band####.out
          ! For a spin-orbit-coupled calculation, the scalar relativistic band
          ! structure will have file names of form band####.out.no_soc, to
          ! distinguish them from the SOC band structure
          if (calculate_perturbative_soc) then
            sr_suffix = ".out"
          else
            sr_suffix = ".out"
          end if

          if(i_band < 10)then
             write(file_name, '(A10,I1,A11)') 'bandmlk100', i_band, adjustl(sr_suffix)
          else if(i_band < 100)then
             write(file_name, '(A9,I2,A11)') 'bandmlk10', i_band, adjustl(sr_suffix)
          else if(i_band < 1000)then
             write(file_name, '(A8,I3,A11)') 'bandmlk1', i_band, adjustl(sr_suffix)
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

              write(88,'(F12.5,2X,F15.5)',ADVANCE='NO') occ_numbers(i_state,1,i_k_point), &! *n_k_points, &
                   (KS_eigenvalue(i_state,1,i_k_point)-chemical_potential)* hartree
            end do
            write(88,'()')

            if(n_spin ==2)then

              write(89,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)

              do  i_state = 1,  n_states

                 write(89,'(F12.5,2X,F15.5)',ADVANCE='NO') occ_numbers(i_state,2,i_k_point), &!*n_k_points, &
                      (KS_eigenvalue(i_state,2,i_k_point)-chemical_potential)* hartree

              end do
              write(89,'()')
            end if

          end do

          close(88)

          if(n_spin==2) close(89)

       end if ! myid == 0

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
       call localorb_info(info_str, use_unit,'(A)')

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

           ! The rest is the output operation, which involves tasks 0 only
           if (myid.eq.0) then

             k(:) = band_k_frac(i_k_point, i_band)

             do i_spin = 1, n_spin, 1

               call output_complex_eigenvector &
               ( KS_eigenvector_tmp(:,:,i_spin), KS_eigenvalue(:,i_spin,i_k_point), chemical_potential, &
                 occ_numbers(:,i_spin,i_k_point), i_band, i_k_point, i_spin, k(1), k(2), k(3))

             enddo ! end spin loop

           end if ! end restriction of operations to task number 0

       enddo ! end loop over k-points in current band


       end if

       ! Now, output the SOC band structure if it's been requested

       if (calculate_perturbative_soc) then
         if (allocated(KS_eigenvalue_soc_perturbed)) then
           call aims_deallocate(KS_eigenvalue_soc_perturbed,&
                "KS_eigenvalue_soc_perturbed")
         end if
         ! Because we will never store the SOC-perturbed eigenvectors, here we
         ! compute and save all eigenvalues in the second-variational window,
         ! since they're cheap
         call aims_allocate( KS_eigenvalue_soc_perturbed, n_states_soc, 1, n_k_points, &
              "KS_eigenvalue_soc_perturbed" )
         KS_eigenvalue_soc_perturbed = 0.0d0

         if (allocated(eigenvec_soc_wf_basis)) then
           call aims_deallocate( eigenvec_soc_wf_basis, "eigenvec_soc_wf_basis" )
         end if
         call aims_allocate( eigenvec_soc_wf_basis, n_states_soc, n_states_soc, "eigenvec_soc_wf_basis" )

         if (allocated(occ_numbers_soc)) then
           call aims_deallocate(occ_numbers_soc, "occ_numbers_soc")
         end if
         call aims_allocate( occ_numbers_soc, n_states_soc, 1, n_k_points, &
              "occ_numbers_soc")

         if (allocated(my_k_points)) call aims_deallocate(my_k_points, "my_k_points")
         call aims_allocate( my_k_points, n_k_points_task,             "my_k_points")

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
         ! with the exception of no output and lack of determination of
         ! Fermi level, as this has already been set in that main SOC routine.
         do  i_k_point = 1,  n_k_points_task
           this_k_point = my_k_points(i_k_point)

           call construct_SOC_Hamiltonian( n_hamiltonian_matrix_size, soc_matrix, &
                n_basis, n_states, KS_eigenvector(:,:,:,1), KS_eigenvector_complex(:,:,:,i_k_point), this_k_point, &
                n_states_soc, n_states_soc, SOC_Hamiltonian )

           call perform_soc_perturbation( n_states_soc, n_states_soc, SOC_Hamiltonian, KS_eigenvalue(1,1,this_k_point),&
                KS_eigenvalue_soc_perturbed(1, 1, this_k_point), dummy , &
                n_states_soc, n_states_soc, eigenvec_soc_wf_basis )
! add
           call convert_wf_basis_to_compute_basis( &
                n_basis,      n_states,           KS_eigenvector(1,1,1,1), KS_eigenvector_complex(1,1,1,i_k_point),&
                n_states_soc, n_states_soc,       eigenvec_soc_wf_basis, &
                n_basis_soc,  n_saved_states_soc, KS_eigenvector_soc_perturbed(1,1,1,i_k_point) )
! end
         end do
! add
       if (allocated(mulliken_decomp)) call aims_deallocate(mulliken_decomp, "mulliken_decomp" )
       call aims_allocate(mulliken_decomp, 0,l_wave_max, 1,n_atoms, 1,n_saved_states_soc, 1,2, 1,n_k_points, "+mulliken_decomp")

!       if (allocated(mulliken_decomp)) call aims_deallocate(mulliken_decomp, "mulliken_decomp")
!       call aims_allocate(mulliken_decomp, l_wave_max, n_atoms, n_saved_states_soc, 2, n_k_points, "mulliken_decomp")
       mulliken_decomp = 0.d0

       do i_spin = 1, 2, 1
         i_k = 0
         if (i_spin .eq. 1) then
           ! Spin-up components are requested
           basis_offset = 0
         else
           ! Spin-dn components are requested
           basis_offset = n_basis_soc_coll/2
         end if

         do i_k_point = 1, n_k_points
           ! write(use_unit,*) i_k_point
           if  ( myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
             i_k = i_k + 1
             do i_state = 1, n_saved_states_soc, 1
               do i_cell_i = 1,n_cells_in_hamiltonian-1
                 do i_basis_2 = 1, n_basis_soc_coll/2
                   if( index_hamiltonian(1,i_cell_i, i_basis_2) > 0 )then
                     i_index = index_hamiltonian(1,i_cell_i, i_basis_2)-1
                     do i_size = index_hamiltonian(1,i_cell_i, i_basis_2),index_hamiltonian(2,i_cell_i, i_basis_2)
                       i_index = i_index + 1
                       i_basis_1 =  column_index_hamiltonian(i_index)
                       mul_temp =  KS_eigenvector_soc_perturbed(basis_offset+i_basis_1, i_state,1,i_k) * &
                            dconjg(KS_eigenvector_soc_perturbed(basis_offset+i_basis_2, i_state,1,i_k)) &
                            * dconjg(k_phase(i_cell_i,i_k_point)) &
                            * overlap_matrix(i_index)
                       mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) = &
                            mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) + &
                            dble(mul_temp)
                       ! 2nd pass: must average all off-diagonal matrix elements
                       ! (but not the diagonal)
                       if (i_basis_1.ne.i_basis_2) then
                         mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) = &
                              mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) + &
                              dble(mul_temp)

                       end if
                     end do ! i_size
                   end if
                 end do ! i_basis_2
               end do ! i_cell
             end do ! i_state
           end if
         end do ! i_k_point
       end do ! i_spin
! end add

         if (use_scalapack .and. my_scalapack_id .ne. 0) then
           KS_eigenvalue_soc_perturbed = 0.0d0
         end if

         call sync_vector(KS_eigenvalue_soc_perturbed, n_states_soc*n_k_points, mpi_comm_global)
! add
         call sync_vector(mulliken_decomp(0,1,1,1,1),(1+l_wave_max)*n_atoms*n_saved_states_soc*2*n_k_points)
! end
         ! Get occupation numbers
         call convert_sr_to_soc_environment ()
         call check_norm_p0(chemical_potential_soc, KS_eigenvalue_soc_perturbed, n_electrons, occ_numbers_soc, dummy, &
              dummy_int)
         call revert_soc_to_sr_environment ()

if(myid==0)then
            sr_suffix = ".out"
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

           !write bands
           do i_k_point = 1,  n_k_points
             k(:) = band_k_frac(i_k_point, i_band)
             i_x = k(1)
             i_y = k(2)
             i_z = k(3)

             write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, i_x, i_y, i_z

             do i_state = 1,  n_states_soc
               write(88,'(F12.5,F15.5)',ADVANCE='NO') occ_numbers_soc(i_state,1,i_k_point), &! *n_k_points_band, &
                    (KS_eigenvalue_soc_perturbed(i_state,1,i_k_point)-chemical_potential_soc)* hartree
             end do

             write(88,'()')
            end do

            close(88)

! add mlk
! format change write
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

           open(50, file=file_name)

           do i_k_point = 1, n_k_points
              k(:) = band_k_frac(i_k_point, i_band)
              write(50,'(A,I5,A,F12.8,1X,F12.8,1X,F12.8,A)') &
                "k point number: ",i_k_point, ": ( ", k(1), k(2), k(3), " )"
              write(50,'(4X,A,7X,A,2X,A,1X,A, 7X,A, 7X,A, 10(9X,A3))') &
                "State", "eigenvalue", "occ.number", "atom", "spin", "total"
              do i_state=1,n_saved_states_soc,1
                if (i_state <= homo_id + num_orbit_output_half .and. &
                   i_state > homo_id - num_orbit_output_half) then
                   write(50,'(4X,A,2X,I7)') "State", i_state
                end if
                do i_atom = 1, n_atoms, 1
                  if(species_pseudoized(species(i_atom))) cycle

                  do i_spin = 1, 2, 1
!                  do i_spin = 1, n_spin, 1

                    if (i_state <= homo_id + num_orbit_output_half .and. i_state > &
                        homo_id - num_orbit_output_half) then
                      write(50,'(2X,I7,2X,F15.5,2X,F10.7,2X, I3,2X, I3,2X,F10.5,10(2X,F10.5))') &
                          i_state,(KS_eigenvalue_soc_perturbed(i_state,1,i_k_point) - chemical_potential)*hartree, &
                          occ_numbers_soc(i_state,1,i_k_point), i_atom, i_spin, &
                          sum(mulliken_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,i_k_point)),&
                          (mulliken_decomp(i_l,i_atom,i_state,i_spin,i_k_point), i_l=0,l_shell_max(species(i_atom)) )
                    end if
                  enddo
               end do
            enddo

         enddo

         close(50)
         end if ! myid == 0 write out

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

           write(info_str,'(2X,A)') 'Not writing out SOC-perturbed eigenvectors for band structures.'
           call localorb_info ( info_str )

         end if

    end if ! calculate_perturbative_soc

    deallocate(occ_numbers)
    deallocate(KS_eigenvalue)
    call aims_deallocate( KS_eigenvector_complex, "KS_eigenvector_complex" )

    end do ! i_band

    ! Output the (scalar-relativistic) band gap predicted by the band structure
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

    write(info_str,'(A)') ''
    call localorb_info ( info_str )

    ! Allocatable variables from this subroutine
    if (allocated(hamiltonian_w))               call aims_deallocate(hamiltonian_w,                             "hamiltonian_w")
    if (allocated(overlap_matrix_w))            call aims_deallocate(overlap_matrix_w,                       "overlap_matrix_w")
    if (allocated(hamiltonian_w_complex))       call aims_deallocate(hamiltonian_w_complex,             "hamiltonian_w_complex")
    if (allocated(overlap_matrix_w_complex))    call aims_deallocate(overlap_matrix_w_complex,       "overlap_matrix_w_complex")
    if (allocated(work_ham))                    call aims_deallocate(work_ham,                                       "work_ham")
    if (allocated(work_ovl))                    call aims_deallocate(work_ovl,                                       "work_ovl")
    if (allocated(k_weights_old))               call aims_deallocate(k_weights_old,                             "k_weights_old")
    if (allocated(fock_matrix_complex))         call aims_deallocate(fock_matrix_complex,                 "fock_matrix_complex")
    if (allocated(fock_matrix_complex_SR))      call aims_deallocate(fock_matrix_complex_SR,           "fock_matrix_complex_SR")
    if (allocated(KS_eigenvector_tmp))          call aims_deallocate(KS_eigenvector_tmp,                   "KS_eigenvector_tmp")
    if (allocated(SOC_Hamiltonian))             call aims_deallocate(SOC_Hamiltonian,                         "SOC_Hamiltonian")
    if (allocated(soc_matrix))                  call aims_deallocate(soc_matrix,                                   "soc_matrix")
    if (allocated(eigenvec_soc_wf_basis))       call aims_deallocate(eigenvec_soc_wf_basis,             "eigenvec_soc_wf_basis")
    if (allocated(my_k_points))                 call aims_deallocate(my_k_points,                                 "my_k_points")
    if (allocated(mulliken_decomp))             call aims_deallocate(mulliken_decomp,                         "mulliken_decomp")

    ! Allocatable variables from other subroutines (which are tracked)
    if (allocated(KS_eigenvector_complex))      call aims_deallocate(KS_eigenvector_complex,           "KS_eigenvector_complex")
    if (allocated(KS_eigenvalue_soc_perturbed)) call aims_deallocate(KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed")
    if (allocated(occ_numbers_soc))             call aims_deallocate(occ_numbers_soc,                         "occ_numbers_soc")

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

  end subroutine out_plot_band_soc_mulliken
!******

