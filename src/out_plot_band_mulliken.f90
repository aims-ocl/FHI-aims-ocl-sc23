subroutine out_plot_band_mulliken()
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
                              revert_soc_to_sr_environment
    use timing
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use generate_aims_uuid, only: write_aims_uuid
    use hdf5_output, only: output_complex_overlap_matrix, &
                           output_complex_eigenvector, &
                           output_complex_hamiltonian_matrix
    use dimensions_soc, only : n_states_soc
    use basis, only: basis_l
    use geometry, only: species
    use synchronize_mpi
    use ks_wrapper, only: solve_KS_elsi_serial
!  PURPOSE
!   The subroutine plots mulliken charge decomposition of any k point in the band k path.
!   This works only with lapack type of eigenvectors.
!   This routine combines the out_plot_band and output_mulliken(lapack part).
!   Developed by Chi Liu (garnett.liu@duke.edu)

    implicit none
    real*8,    dimension(:,:),allocatable :: hamiltonian_w
    real*8,    dimension(:),  allocatable :: overlap_matrix_w
    complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
    complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex

    real*8, dimension(:,:,:),allocatable :: work_ham
    real*8, dimension(:,:),allocatable :: work_ovl

    integer:: i_band,    i_cell_n, i_spin
    integer:: i_cell(3), i_counter
    integer:: i_k_point, i_state, j_state, i_k
    integer :: i_k_point_old, n_k_points_old, i_index, i_k_old, i_basis_1, i_basis_2
    complex*16 :: k_phase_exx_old
    real*8, dimension(:), allocatable :: k_weights_old
    complex*16, dimension(:,:), allocatable :: fock_matrix_complex
    complex*16, dimension(:,:), allocatable :: fock_matrix_complex_SR
    logical :: real_eigenvectors_old

    real*8:: k(3) !, k_lowest_un_occ(3), k_highest_occ(3), k_lowest_un_occ_whole_system(3), k_highest_occ_whole_system(3)
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

  ! We will do a Mulliken charge decomposition by the following criteria:
  ! Number of electrons in each KS state, per atom, spin channel, and
  ! angular momentum component
  ! mulliken_decomp will contain this decomposition. All derived quantities
  ! are then accessible by way of appropriate sums.
!  real*8, dimension( 0:l_wave_max, n_atoms, n_states, n_spin,n_k_points ) :: mulliken_decomp
  real*8:: i_x, i_y, i_z
  real*8, dimension(:, :, :, :, :), allocatable :: mulliken_decomp
  integer :: i_size, i_l
  complex*16:: mul_temp
  character*40 :: filename
  character(LEN=80) :: uuid_str
  ! darn, need an extra aux array of characters
  character*3,dimension(0:l_wave_max) :: l_channel
  integer :: homo_id
  integer :: num_orbit_output_half

     homo_id = int(n_electrons) / 2
     num_orbit_output_half = band_mulliken_orbit_num / 2
     do i_l = 0, l_wave_max, 1
        write(l_channel(i_l),'(A,I1)') "l=",i_l
     enddo

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "Writing the requested band Mulliken decomposition output (LAPACK version):"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

! Integrates the matrix elements for the Hamiltonian matrix,
    call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, &
         partition_tab, l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

! SOC allocation here

    call aims_allocate(hamiltonian_w_complex, n_basis*(n_basis+1)/2,n_spin, "hamiltonian_w_complex")
    call aims_allocate(overlap_matrix_w_complex, n_basis*(n_basis+1)/2, "overlap_matrix_w_complex")
    if(use_periodic_hf) then
      call aims_allocate(fock_matrix_complex,   n_basis*(n_basis+1)/2,n_spin, "fock_matrix_complex")
    end if
    if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
      call aims_allocate(fock_matrix_complex_SR, n_basis*(n_basis+1)/2, n_spin, "fock_matrix_complex_SR")
    end if

    ! dummy allocations for the real arrays.
    ! These trigger a compiler warning for -check pointers otherwise
    ! as they are included (but never touched!) in construct_hamiltonian_and_ovl below.
    call aims_allocate(hamiltonian_w, 1,1,    "hamiltonian_w" )
    call aims_allocate(overlap_matrix_w,1, "overlap_matrix_w" )

!  if (out_eigenvec) and if (flag_rel.eq.1)
    call aims_allocate(KS_eigenvector_tmp, n_basis,n_states,n_spin,"KS_eigenvector_tmp")

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

       if (allocated(occ_numbers))            deallocate(occ_numbers)
       if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
       if (allocated(KS_eigenvector_complex)) call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex,")
       allocate (occ_numbers(n_states,n_spin,n_k_points))
       allocate( KS_eigenvalue(n_states,n_spin,n_k_points) )
! if ( (flag_rel.eq.1) .or. (out_eigenvec) .or. (calculate_perturbative_soc) )
       call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, n_k_points_task, "KS_eigenvector_complex" )

!        if (out_overlap .or. out_hamiltonian)
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
                    if (myid.eq.  MOD(i_k_point_old, n_tasks) .and. myid<=n_k_points_old ) then
                       i_k_old = i_k_old + 1
                       do i_cell_n = 1, n_cells_bvk
                          k_phase_exx_old = exp((0d0,2d0)*pi*(sum(k(:)*dble(cell_index_bvk(i_cell_n,:)))-&
                          sum(k_point_list(i_k_point_old,:)*dble(cell_index_bvk(i_cell_n,:)))))

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
!             if (out_overlap) then and out_hamiltonian
             ! This routine clears, deallocates and reallocates an internal (temporary) array
             ! that holds the overlap matrix inside lapack - and that may need
             ! to be resized.
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
          end if  !myid ==  MOD(i_k_point, n_tasks)

       end do ! i_k_point
       call sync_eigenvalues( KS_eigenvalue)

! .not.fixed_spin_moment
       call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)

       ! First, we output the scalar-relativistic band structure, here, mulliken
       ! decomposition
       if(myid==0)then

          ! file names for band plotting

          ! For a scalar-relativistic calculation, the scalar relativistic band
          ! structures will have file names of form band####.out
          ! For a spin-orbit-coupled calculation, the scalar relativistic band
          ! structure will have file names of form band####.out.no_soc, to
          ! distinguish them from the SOC band structure
!          if (calculate_perturbative_soc) then
!            sr_suffix = ".out.no_soc"
!          else
            sr_suffix = ".out"
!          end if

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
! n_spin >1, bandmlk200 ???

!          open(88, file=file_name)
!
!          do  i_k_point = 1,  n_k_points
!            k(:) = band_k_frac(i_k_point, i_band)
!            write(88,'(I4,2X,3F15.7)',ADVANCE='NO') i_k_point, k(1), k(2), k(3)
!            do  i_state = 1,  n_states
!              write(88,'(F12.5,2X,F15.5)',ADVANCE='NO') occ_numbers(i_state,1,i_k_point), &! *n_k_points, &
!                   (KS_eigenvalue(i_state,1,i_k_point)-chemical_potential)*hartree
!            end do
!            write(88,'()')
!
!          end do
!
!          close(88)
       end if ! myid == 0
! determine lowest_un_occ and highest_occ
! Mulliken charge analysis
       write(info_str,'(2X,A)') &
         'Performing scalar-relativistic Mulliken charge analysis on all atoms.'
       call localorb_info(info_str)

       write(info_str, *) "n_atoms, n_states, n_spin,n_k_points", n_atoms, n_states, n_spin,n_k_points
       call localorb_info(info_str)
       if (allocated(mulliken_decomp)) call aims_deallocate(mulliken_decomp, "mulliken_decomp" )
       call aims_allocate(mulliken_decomp, 0,l_wave_max, 1,n_atoms, 1,n_states, 1,n_spin, 1,n_k_points, "+mulliken_decomp")
!       allocate(mulliken_decomp( 0:l_wave_max, n_atoms, n_states, n_spin,n_k_points ))
       mulliken_decomp = 0.d0
! not real_eigenvector ???
       do i_spin = 1, n_spin, 1
          i_k = 0
           do i_k_point = 1, n_k_points
!             write(use_unit,*) "i_k_point", i_k_point
              if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                 i_k = i_k + 1

                 do i_state = 1, n_states, 1
                    do i_cell_n = 1,n_cells_in_hamiltonian-1
                       do i_basis_2 = 1, n_basis
                          if( index_hamiltonian(1,i_cell_n, i_basis_2) > 0 )then
                             i_index = index_hamiltonian(1,i_cell_n, i_basis_2)-1
                             do i_size = index_hamiltonian(1,i_cell_n, i_basis_2),index_hamiltonian(2,i_cell_n, i_basis_2)
                                i_index = i_index + 1
                                i_basis_1 =  column_index_hamiltonian(i_index)

                                mul_temp =  KS_eigenvector_complex(i_basis_1, i_state, i_spin,i_k) * &
                                     dconjg(KS_eigenvector_complex(i_basis_2, i_state, i_spin,i_k)) &
                                     * dconjg(k_phase(i_cell_n,i_k_point)) &
                                     * overlap_matrix(i_index)

                                mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) = &
                                mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k_point ) + &
                                dble(mul_temp)

                                ! 2nd pass: must average all off-diagonal matrix
                                ! elements (but not the diagonal)
                                if (i_basis_1.ne.i_basis_2) then
                                   mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) = &
                                   mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k_point ) + &
                                   dble(mul_temp)
                                end if
                             end do ! i_size
                          end if ! index_hamiltonian
                       end do ! i_basis_2
                    end do ! i_cell_n
                 end do ! i_state
              end if  ! myid
           end do ! i_k_point
       end do ! i_spin

  call sync_vector(mulliken_decomp(0,1,1,1,1),(1+l_wave_max)*n_atoms*n_states*n_spin*n_k_points)

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
                          sum(mulliken_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,i_k_point)),&
                          (mulliken_decomp(i_l,i_atom,i_state,i_spin,i_k_point),i_l=0,l_shell_max(species(i_atom)) )
                    end if
                  enddo
               end do
            enddo

         enddo

         close(50)

       end if ! myid .eq. 0


     if (allocated (mulliken_decomp))            call aims_deallocate(mulliken_decomp,                           "mulliken_decomp")
!     deallocate(mulliken_decomp)

    end do ! i_band


end subroutine
