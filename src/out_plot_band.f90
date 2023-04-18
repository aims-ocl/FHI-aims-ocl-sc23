!****s* FHI-aims/out_plot_band
!  NAME
!    out_plot_band
!  SYNOPSIS

  subroutine out_plot_band ( WCC_plane_index )
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
                              revert_soc_to_sr_environment
    use timing
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use generate_aims_uuid, only: write_aims_uuid
    use hdf5_output, only: output_complex_overlap_matrix, &
                           output_complex_eigenvector, &
                           output_complex_hamiltonian_matrix
    use dimensions_soc, only : n_states_soc
    use geometry, only: coords, lattice_vector, species
    use synchronize_mpi
    use ks_wrapper, only: solve_KS_elsi_serial
    use WannierCenters
    use rel_x2c_mod
!  PURPOSE
!   The subroutine plots band structure. This works only with lapack type of eigenvectors.
!   The routine can be called only after self consistant iterations, because it destrois
!   the original k-point information.
!
    implicit none
!  INPUTS
    ! CC: Optional input that allows to hijack the band structure
    !     subroutine for WCC calculations
    integer,intent(in),optional                          :: WCC_plane_index
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

    call get_timestamps(band_time, clock_band_time)

    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    if (.not.WCC_calc) then
      write(info_str,'(2X,A)') "Writing the requested band structure output (LAPACK version):"
    else
      write(info_str,'(2X,A)') "Computing Wannier centers of charge         (LAPACK version):"
    end if
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') "-------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()')
    call localorb_info(info_str,use_unit,'(A)')

    if(flag_rel.ne.REL_x2c.and.flag_rel.ne.REL_4c_dks)then
       call integrate_real_hamiltonian_matrix_p2 &
         ( hartree_potential,   rho, rho_gradient, kinetic_density, &
         partition_tab, l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
    endif

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

    ! CC: Safe-checks:
    if ( (.not.WCC_calc) .and. present(WCC_plane_index) ) then
      write(info_str,'(2x,A)') " *** wcc_plane_index present in out_plot_band call, but not WCC_calc! "
      call aims_stop(info_str)
    end if
    if ( (WCC_calc) .and. (.not.present(WCC_plane_index))) then
      write(info_str,'(2x,A)') " *** wcc_plane_index not present in out_plot_band call! "
      call aims_stop(info_str)
    end if

    ! CC: Now setup  everything for WCC
    if (WCC_calc) then
      ! CC: We assume that the n_electron lowest
      !     states are fully occupied. Safe checks
      !     are done during the run
      WCC_highest_occ = n_electrons

      ! CC: Get basic dimensions as function of WCC_plane_index
      call localorb_info('',use_unit)
      call WCC_get_dimensions( WCC_plane_index, WCC_nr_of_bands, n_k_points, n_k_points_task )

      ! CC: Hijack band_* infrastructure:
      if (allocated(band_begin)) deallocate(band_begin)
      if (allocated(band_end  )) deallocate(band_end  )
      if (allocated(n_points_in_band)) deallocate(n_points_in_band)
      n_plot_band = WCC_nr_of_bands
      allocate(band_begin(WCC_nr_of_bands,3))
      allocate(band_end  (WCC_nr_of_bands,3))
      allocate(n_points_in_band(WCC_nr_of_bands))
      call WCC_get_paths( WCC_plane_index, WCC_nr_of_bands, n_k_points, band_begin, band_end, n_points_in_band )


      ! CC: Now allocate the necessary arrays
      if (allocated(WCC_imag))               call aims_deallocate(WCC_imag       ,"WCC_imag")
      if (allocated(overlap_for_WCC))        call aims_deallocate(overlap_for_WCC,"overlap_for_WCC")
      if (allocated(KS_eigenvector_for_WCC)) call aims_deallocate(KS_eigenvector_for_WCC,"KS_eigenvector_for_WCC")
      call aims_allocate( WCC_imag, WCC_highest_occ, WCC_nr_of_bands,              "WCC_imag")
      call aims_allocate(overlap_for_WCC, n_basis*(n_basis+1)/2, n_k_points_task, "overlap_for_WCC")
      if ( calculate_perturbative_soc) then
        call aims_allocate( KS_eigenvector_for_WCC, n_basis, n_states_soc, 2, n_k_points_task,"KS_eigenvector_for_WCC")
      else
        !CC: Ugly, but pragmatic: We use the SOC dimensions also in all other cases
        call aims_allocate( KS_eigenvector_for_WCC, n_basis, 2*n_states,   2, n_k_points_task,"KS_eigenvector_for_WCC")
      end if 
    end if

    do i_band = 1, n_plot_band

       n_k_points =  n_points_in_band(i_band)

       write(info_str,'()')
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(2X,A,I4,A,I4,A)') "Treating all ",n_k_points," k-points in band plot segment #", i_band, ":"
       call localorb_info(info_str,use_unit,'(A)')
       if (WCC_calc) then
         write(info_str,'(2X,A)') " Sampling over :"
         call localorb_info(info_str,use_unit,'(A)')
         do  i_k_point = 1, n_k_points
           write(info_str,'(2X,A,3F7.3)') " - k = ", band_k_frac(i_k_point, i_band)
           call localorb_info(info_str,use_unit,'(A)')
         end do
       end if
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

       if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
         ! For Q4C, I need to do the real space integration here, with the new k_phase
          scf_iteration = 0 ! Set scf_iteration=0 so that the integrate_real_hamiltonian_matrix_p2 subroutine
                            ! can generate v, w, ss, ssc2, instead of only v
          call rel_matrix_initialize()
          call integrate_real_hamiltonian_matrix_p2 &
            ( hartree_potential,   rho, rho_gradient, kinetic_density, &
            partition_tab, l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
          call integrate_ovlp_matrix_p2 ( partition_tab, l_shell_max, overlap_matrix )
       endif

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

       if ( (flag_rel.eq.1) .or. (out_eigenvec) .or. (calculate_perturbative_soc) .or. (WCC_calc) ) then
          call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, n_k_points_task, "KS_eigenvector_complex" )
       else
          call aims_allocate( KS_eigenvector_complex, n_basis,n_states,n_spin, 1, "KS_eigenvector_complex" )
       end if

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
          endif

          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then

             i_k = i_k + 1

             if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
               call rel_band_plot_k(i_k_point,KS_eigenvalue(1,1,i_k_point))
             else

               call construct_hamiltonian_and_ovl(hamiltonian, overlap_matrix, &
                    hamiltonian_w, overlap_matrix_w, &
                    hamiltonian_w_complex, overlap_matrix_w_complex, i_k_point, work_ham, work_ovl)
               
               ! CC: Need to store overlap_matrix_w_complex for WCC calculation
               if ( WCC_calc) then
                  overlap_for_WCC(:,i_k) = overlap_matrix_w_complex
               end if 

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

                  if (flag_rel == 1 .or. out_eigenvec .or. calculate_perturbative_soc .or. WCC_calc) then
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

             endif ! rel / nonrel
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

       ! CC: Store KS_evecs already here if no SOC requested
       !     Also, check occupation numbers:
       if ( WCC_calc .and. (.not. calculate_perturbative_soc)) then
         call WCC_convert_and_store_KS_evec_no_SOC( &
           n_basis, n_states, n_spin, n_k_points, n_k_points_task, & 
           occ_numbers, KS_eigenvalue, KS_eigenvector_complex , &
           KS_eigenvector_for_WCC, WCC_highest_occ )
       end if

       ! First, we output the scalar-relativistic band structure

       if(myid==0 .and. (.not.WCC_calc) ) then

          ! file names for band plotting

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
       call localorb_info(info_str, use_unit,'(A)')

       if (.not. WCC_calc) then
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

           ! CC: Convert KS_evec from SOC(n_states_soc) to standard (n_basis)
           !     representation and store it in KS_eigenvector_for_WCC
           if (WCC_calc) then
             call WCC_convert_and_store_KS_evec( i_k_point, &
               n_basis, n_states, n_spin, n_k_points_task, KS_eigenvector_complex , &
               n_states_soc, eigenvec_soc_wf_basis, &
               KS_eigenvector_for_WCC )
           end if 

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

         ! CC: Check that no fractional occupations occur,
         !     since metals are non-sensical in this context.
         if (WCC_calc) then
           call WCC_check_occ_numbers( n_states_soc, 1,  n_k_points, occ_numbers_soc, WCC_highest_occ )
         end if

         if(myid==0 .and. (.not.WCC_calc) )then
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


         if(.not.WCC_calc) then

           write(info_str,'()') 
           call localorb_info(info_str, use_unit,'(A)')

           write(info_str,'(A,I4)')  'Spin-orbit-coupled "band gap" along reciprocal space direction number: ',i_band
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
         end if

         if (allocated(occ_numbers_soc))                   call aims_deallocate(occ_numbers_soc,"occ_numbers_soc")
         if (allocated(my_k_points))                       call aims_deallocate(my_k_points,"my_k_points")
    end if ! calculate_perturbative_soc

    if ( WCC_calc ) then 
      write(info_str,'(2X,A,I6)') "|-> Computing Wannier center evolution at k-point number: ", i_band
      call localorb_info(info_str,use_unit,'(A)')
      call compute_WCC(i_band,n_states,n_basis,n_k_points,n_k_points_task,WCC_highest_occ,WCC_nr_of_bands,KS_eigenvector_for_WCC, overlap_for_WCC, WCC_imag)
    end if

    if (allocated(occ_numbers))            deallocate(occ_numbers)
    if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
    if (allocated(KS_eigenvector_complex)) call aims_deallocate( KS_eigenvector_complex, "KS_eigenvector_complex" )

    end do ! i_band

    ! Output the (scalar-relativistic) band gap predicted by the band structure
    if (.not.WCC_calc ) then  
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
    end if

    ! CC: Output of WCC for plane WCC_plane_index:
    if (WCC_calc) then
      write(info_str,'(2X,A,I6)') "*-> Finished: Writing Wannier center evolution to file."
      call localorb_info ( info_str )
      if (myid.eq.0) then     
        call ouput_WCC_to_file(WCC_plane_index, WCC_highest_occ, WCC_nr_of_bands, WCC_imag )
      endif
      if(allocated(KS_eigenvector_for_WCC)) call aims_deallocate(KS_eigenvector_for_WCC,"KS_eigenvector_for_WCC")
      if(allocated(overlap_for_WCC))        call aims_deallocate(overlap_for_WCC,"overlap_for_WCC")
      if(allocated(WCC_imag))               call aims_deallocate(WCC_imag,"WCC_imag")
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
    if (WCC_calc) then
      write(info_str,'(2X,A)') &
           "Wannier center evolution                                :  max(cpu_time) wall_clock(cpu1)"
    else
      write(info_str,'(2X,A)') &
           "Band Structure                                          :  max(cpu_time) wall_clock(cpu1)"
    end if
    call localorb_info(info_str)
    write(info_str, "(2X, A,F15.3,F17.3)")  "| Total Time                                            :", &
         band_time, clock_band_time
    call localorb_info(info_str)
    write(info_str,'(A)') "------------------------------------------------------------"
    call localorb_info(info_str)

  end subroutine out_plot_band
!******

