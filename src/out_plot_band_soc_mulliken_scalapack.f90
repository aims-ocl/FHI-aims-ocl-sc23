!****s* FHI-aims/out_plot_band_scalapack
!  NAME
!    out_plot_band_scalapack
!  SYNOPSIS
subroutine out_plot_band_soc_mulliken_scalapack ()
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
! add mlk
                                      species_pseudoized
! end mlk
    use dimensions,            only : calculate_perturbative_soc, n_atoms, &
                                      n_basis, n_states, n_spin, n_k_points, &
                                      n_hamiltonian_matrix_size, use_periodic_hf, &
                                      n_plot_band, use_lc_wpbeh, &
! add mlk
                                      l_wave_max
! end mlk
    use constants,             only : bohr, hartree, pi
    use geometry,              only : coords, lattice_vector, species
    use pbc_lists,             only : k_point_list, k_weights, k_phase, cell_index_bvk,&
                                      n_cells, n_cells_bvk, cell_index, band_k_frac, &
! add mlk
                                      Cbasis_to_atom
! end mlk
    use mpi_tasks,             only : myid, n_tasks, mpi_comm_global, mpi_integer, aims_stop
    use runtime_choices,       only : collect_eigenvectors, fixed_spin_moment, &
                                      flag_rel, hybrid_coeff, out_eigenvec, &
                                      out_hamiltonian, out_eigenvec, out_overlap, &
                                      real_eigenvectors, use_load_balancing, &
                                      use_local_index, n_points_in_band, &
                                      fixed_spin_moment_electrons, &
                                      band_mulliken_orbit_num, use_elsi
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
                                      solve_evp_scalapack_complex, &
! add mlk
                                      sc_desc, ovlp_complex_stored, &
                                      save_overlap_scalapack, nprow, npcol, &
                                      my_scalapack_comm_work
! end mlk

    use synchronize_mpi,       only : sync_eigenvalues
    use synchronize_mpi_basic, only : sync_vector, sync_vector_complex, &
                                      bcast_integer, bcast_integer_vector, &
                                      send_real_vector, receive_real_vector
    use hartree_fock_p0,       only : hf_exchange_matr_complex, &
                                      hf_exchange_matr_complex_SR
    use soc_utilities,         only : perform_soc_perturbation, &
                                      convert_sr_to_soc_environment, &
                                      revert_soc_to_sr_environment, &
! add mlk
                                      convert_wf_basis_to_compute_basis
! end mlk
    use load_balancing,        only : use_batch_permutation, batch_perm, &
                                      n_bp_integ, &
                                      set_full_local_ham, &
                                      set_full_local_ovlp
    use scalapack_soc,         only : mxld_soc, mxcol_soc, &
! add mlk
                                      mxld_soc_vec, mxcol_soc_vec, sc_desc_soc_vec, &
                                      l_row_soc_vec, l_col_soc_vec
! end mlk
    use timing,                only : tot_clock_time_band_dos, tot_time_band_dos, &
                                      get_times, get_timestamps
    use aims_memory_tracking,  only : aims_allocate, aims_deallocate
    use hdf5_output,           only : output_overlap_scalapack, &
                                      output_complex_EV_scalapack, &
                                      output_complex_KS_OCC_scalapack
    use dimensions_soc,        only : n_states_soc, &
! add mlk
                                      n_saved_states_soc, soc_saved_state_start, &
                                      n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll
    use basis,                 only : basis_l
    use ks_wrapper,            only : solve_KS_elsi_parallel
    use elsi_wrapper,          only : aims_elsi_reinit_scf

! end mlk

!   The subroutine plots mulliken charge decomposition of any k point in the
!   band k path.
!   This is a memory optimized version works with scalapack type of eigenvectors.
!   This routine combines the out_plot_band and output_mulliken(lapack part).
!   Deal with SOC case.
!   Developed by Ruyi Song (ruyi.song@duke.edu), Chi Liu (garnett.liu@duke.edu)

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
!  FIXME
!   This routine uses the infrastructure of the scf cycle. The eigenvalue equation is
!   solved in chuncks (of size n_k_points) as in the scf. Since scalapack is not allocated new
!   this only works of the scf already uses complex eigenvectors (for k_grid other than 1 1 1, 2,1,1, ... 2,2,2).
!   Need to extend to the real case as well. (Note that the lapack version works for both cases).

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
    integer:: i_k_point, i_state, j_state, n_k_points_temp, i, j
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

! add mlk
!    complex*16, dimension(:, :, :, :), allocatable :: eigenvec_soc ! actually the KS_eigenvector_soc_perturbed
    real*8, dimension(:, :, :, :, :), allocatable :: mulliken_decomp
!  real*8, dimension( 0:l_wave_max, n_atoms, n_saved_states_soc, 2, n_points_in_band(1) ) :: mulliken_decomp

    complex*16, dimension(:,:),     allocatable :: tmp_complex
    character(LEN=80) :: uuid_str
    integer :: i_size, i_l
    complex*16:: mul_temp
    character*40 :: filename
    ! darn, need an extra aux array of characters
    character*3,dimension(0:l_wave_max) :: l_channel
    integer :: homo_id
    integer :: num_orbit_output_half
    integer :: i_l_wave_max

    integer :: mpierr
    integer :: i_k, k_temp1, k_temp2, i_core, cyclenumber
    integer, dimension(:), allocatable :: k_band_each, id_each

     homo_id = int(n_electrons)
     num_orbit_output_half = band_mulliken_orbit_num
     do i_l = 0, l_wave_max, 1
        write(l_channel(i_l),'(A,I1)') "l=",i_l
     enddo
! end mlk

    call get_timestamps(band_time, clock_band_time)

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


if (real_eigenvectors ) then
     if (myid.eq.0) then
       write(use_unit,'(1X,A)') "* Warning: Bands plotting currently not supported for real eigenvectors. "
       write(use_unit,'(1X,A)') "* Real eigenvectors are used for a  k_grid  of 1 1 1, 2 1 1, 1 2 1, 1 1 2,..., 2 2 2."
       write(use_unit,'(1X,A)') "* You could:"
       write(use_unit,'(1X,2A)') "*  - use a k_grid other than above ", &
          "(i.e. request in at least one direction 3 or more k-points)"
       write(use_unit,'(1X,A)') "*  - shift by a VERY tiny k_offset ( 0.00001 0 0 or so) "
       write(use_unit,'(1X,A)') "*  - use lapack"
       write(use_unit,'(1X,A)') "* "
       write(use_unit,'(1X,A)') "* The real fix is to deallocate the entire (real) scalapack infrastructure,"
       write(use_unit,'(1X,A)') "* and reinitialize the same infrastructure but for complex matrices."
       write(use_unit,'(1X,A)') "* The possible catch is that complex matrices will take twice as much memory"
       write(use_unit,'(1X,A)') "* going into each k-point."
       write(use_unit,'(1X,A)') "* "
       write(use_unit,'(1X,A)') "* ... quitting bands plotting routine."
    end if
else

    ! Integrate the real-space Hamiltonian.  When local indexing is enabled, also
    ! integrate the real-space overlap matrix.
    if (use_local_index) then
      if (use_load_balancing) then
        use_batch_permutation = n_bp_integ
        call aims_allocate( my_hamiltonian_matrix, batch_perm(n_bp_integ)%n_local_matrix_size,n_spin, "my_hamiltonian_matrix")
        call aims_allocate( ovlp_work,             batch_perm(n_bp_integ)%n_local_matrix_size,                   "ovlp_work" )
      else
        call aims_allocate( my_hamiltonian_matrix, n_hamiltonian_matrix_size,n_spin, "my_hamiltonian_matrix" )
        call aims_allocate( ovlp_work,             n_hamiltonian_matrix_size,                    "ovlp_work" )
      end if

      call integrate_real_hamiltonian_matrix_p2 &
           ( hartree_potential,   rho, rho_gradient, kinetic_density, &
             partition_tab, l_shell_max, &
             en_xc, en_pot_xc, my_hamiltonian_matrix, en_vdw, en_pot_vdw )
      call integrate_ovlp_matrix_p2 &
           ( partition_tab, l_shell_max, ovlp_work )
    else
      call aims_allocate( my_hamiltonian_matrix, 1,1, "my_hamiltonian_matrix" )
      call aims_allocate( ovlp_work, 1,                           "ovlp_work" )
      call integrate_real_hamiltonian_matrix_p2 &
           ( hartree_potential,   rho, rho_gradient, kinetic_density, &
             partition_tab, l_shell_max, &
             en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
    end if

    ! WPH:  I don't know why this is only causing problems now, but these arrays should be allocated when
    !       collect_eigenvectors.eq..true. anyway.
    if (collect_eigenvectors) then
      if (allocated( KS_eigenvector )) call aims_deallocate( KS_eigenvector,                           "KS_eigenvector" )
      call aims_allocate( KS_eigenvector, 1, 1, 1, 1,                                                  "KS_eigenvector" )

      if (allocated( KS_eigenvector_complex )) call aims_deallocate( KS_eigenvector_complex,   "KS_eigenvector_complex" )
      call aims_allocate( KS_eigenvector_complex ,n_basis, n_states, n_spin, 1,                "KS_eigenvector_complex" )
    else
      if (allocated( KS_eigenvector )) call aims_deallocate( KS_eigenvector,                           "KS_eigenvector" )
      call aims_allocate( KS_eigenvector, 1, 1, 1, 1,                                                  "KS_eigenvector" )

      if (allocated( KS_eigenvector_complex )) call aims_deallocate( KS_eigenvector_complex,   "KS_eigenvector_complex" )
      call aims_allocate( KS_eigenvector_complex ,1, 1, n_spin, 1,                             "KS_eigenvector_complex" )
    end if

    ! If SOC enabled, set up initial matrices
    if (calculate_perturbative_soc) then
      ! The eigenvector arrays may have been deallocated by now, and we need
      ! them
      if (allocated( eigenvec )) call aims_deallocate( eigenvec,                                             "eigenvec" )
      call aims_allocate( eigenvec, 1, 1, 1,                                                                 "eigenvec" )

      if (allocated( eigenvec_complex )) call aims_deallocate( eigenvec_complex,                     "eigenvec_complex" )
      call aims_allocate( eigenvec_complex, mxld, mxcol, n_spin,                                     "eigenvec_complex" )

      call aims_allocate( eigenvec_soc_wf_basis, mxld_soc, mxcol_soc,                           "eigenvec_soc_wf_basis" )
      call aims_allocate( soc_ham, mxld_soc, mxcol_soc,                                                       "soc_ham" )

      call aims_allocate(k_band_each, n_tasks, "k_band_each")

      if (use_local_index.and.use_load_balancing) then
        call aims_allocate(soc_matrix, batch_perm(n_bp_integ)%n_local_matrix_size, 3, "+soc_matrix")
        ld_soc_matrix = batch_perm(n_bp_integ)%n_local_matrix_size
      else
        call aims_allocate(soc_matrix, n_hamiltonian_matrix_size, 3,                  "+soc_matrix")
        ld_soc_matrix = n_hamiltonian_matrix_size
      end if

      call integrate_soc_matrix (rho, hartree_potential, partition_tab, soc_matrix )
   else ! otherwise set up dummy arrays
      call aims_allocate( eigenvec_soc_wf_basis, 1, 1,                                          "eigenvec_soc_wf_basis" )
      call aims_allocate( soc_ham, 1,1,                                                                       "soc_ham" )
      call aims_allocate( soc_matrix, 1,1,                                                                 "soc_matrix" )
    end if

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

       call aims_allocate(id_each, n_k_points_band , "id_each") ! allocate the core id for calculating each k_point in band

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
       end do ! i_k_point, k phase and coordinates
! add mlk
!      if (myid == 0) print*, "after loop to get k_phase_band and k coordinates"
! end mlk
       if(abs(sum(k_weights)-1) >1e-10)then
          write(use_unit,*) 'Error: sum of k-vector weights is not one!', sum(k_weights)
          stop
       end if

       ! Allocate those arrays that change with the bands

       if (allocated(occ_numbers))            deallocate(occ_numbers)
       if (allocated(KS_eigenvalue))          deallocate(KS_eigenvalue)
       allocate (occ_numbers(n_states,n_spin,n_k_points_band))
       allocate (KS_eigenvalue(n_states,n_spin,n_k_points_band))

       ! Here we need to allocate the SOC-specific arrays
       if (calculate_perturbative_soc) then
         ! Because we will never store the SOC-perturbed eigenvectors, here we
         ! compute and save all eigenvalues in the second-variational window,
         ! since they're cheap
         if (allocated(KS_eigenvalue_soc_perturbed)) call aims_deallocate(KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed" )
         call aims_allocate( KS_eigenvalue_soc_perturbed, n_states_soc, 1, n_k_points_band,            "KS_eigenvalue_soc_perturbed" )
         if (allocated(occ_numbers_soc)) call aims_deallocate(occ_numbers_soc,                                     "occ_numbers_soc" )
         call aims_allocate( occ_numbers_soc, n_states_soc, 1, n_k_points_band,                                    "occ_numbers_soc" )

         if (allocated(KS_eigenvalue_temp)) call aims_deallocate(KS_eigenvalue_temp,                            "KS_eigenvalue_temp" )
         call aims_allocate (KS_eigenvalue_temp, n_states,n_spin,                                               "KS_eigenvalue_temp" )

! add mlk, here? or up, does it change with bands? yes, though dimension doesn't
!         if (myid == 0) print*, "before re allocate KS_eigenvector_soc_perturbed"
!         if (allocated(KS_eigenvector_soc_perturbed)) deallocate(KS_eigenvector_soc_perturbed)
!         allocate(KS_eigenvector_soc_perturbed(mxld_soc_vec, mxcol_soc_vec, 1, 1))
         if (allocated(KS_eigenvector_soc_perturbed)) call aims_deallocate(KS_eigenvector_soc_perturbed, "KS_eigenvector_soc_perturbed")
         call aims_allocate(KS_eigenvector_soc_perturbed, mxld_soc_vec, mxcol_soc_vec, 1, 1, "KS_eigenvector_soc_perturbed")

!         if (myid == 0) print*, "after re allocate KS_eigenvector_soc_perturbed"
!         if (allocated(eigenvec_soc)) call aims_deallocate(eigenvec_soc, "eigenvec_soc")
!         call aims_allocate (eigenvec_soc, mxld_soc_vec, mxcol_soc_vec, 1, 1, "eigenvec_soc")

!         if (allocated (mulliken_decomp)) deallocate(mulliken_decomp)
!         if (allocated (mulliken_decomp)) call aims_deallocate(mulliken_decomp, "mulliken_decomp")
!         if (myid == 0) print*, "before re allocate mulliken_decomp"

!         allocate(mulliken_decomp(0:l_wave_max, n_atoms, n_saved_states_soc, 2, n_k_points_band))

!         if (myid == 0) print*, "after re allocate mulliken_decomp"
         KS_eigenvector_soc_perturbed(:, :, :, :) = 0.d0
!         if (myid == 0) print*, "after KS_eigenvector_soc_perturbed = 0"
!         KS_eigenvector_soc_perturbed = (0.0d0, 0.0d0)
!         eigenvec_soc = (0.0d0, 0.0d0)
!         do i_k_point = 1, n_k_points_band
!           mulliken_decomp(:, :, :, :, i_k_point) = 0.0d0
!         end do
!         do i_k_point = 1, n_k_points_band
!           do i_spin = 1, 2
!             do i_state = 1, n_saved_states_soc
!               do i_atom = 1, n_atoms
!                 do i_l_wave_max = 0, l_wave_max
!           mulliken_decomp(i_l_wave_max, i_atom, i_state, i_spin, i_k_point) = 0.0d0
!                 end do
!               end do
!             end do
!           end do
!         end do
         if (allocated(mulliken_decomp)) call aims_deallocate(mulliken_decomp, "mulliken_decomp" )
         cyclenumber = (n_k_points_band/n_k_points)+1
         call aims_allocate(mulliken_decomp, 0,l_wave_max, 1,n_atoms, 1,n_saved_states_soc, 1,2, 1,cyclenumber, "+mulliken_decomp")
         mulliken_decomp(:, :, :, :, :) = 0.0d0
!         if (myid == 0) print*, "after mlk = 0"
! end mlk

         occ_numbers_soc = 0.0d0
         KS_eigenvalue_soc_perturbed = 0.0d0
       else
         if (allocated(KS_eigenvalue_soc_perturbed)) call aims_deallocate(KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed" )
         call aims_allocate( KS_eigenvalue_soc_perturbed, 1, 1, 1,                                     "KS_eigenvalue_soc_perturbed" )
         if (allocated(occ_numbers_soc)) call aims_deallocate(occ_numbers_soc,                                     "occ_numbers_soc" )
         call aims_allocate( occ_numbers_soc, 1, 1, 1,                                                             "occ_numbers_soc" )

         if (allocated(KS_eigenvalue_temp))          call aims_deallocate(KS_eigenvalue_temp,                   "KS_eigenvalue_temp" )
         call aims_allocate (KS_eigenvalue_temp, 1,1,                                                           "KS_eigenvalue_temp" )
       end if
! add mlk
!       if (myid == 0) print*, "after calculate_perturbative_soc"
! end mlk
       ! Set KS_eigenvalue and occ_numbers to zero. Otherwise they are filled with random numbers
       ! at indices which are not written for a given (partial)  eigenwert solution of lower dimension n_k_points.
       ! This  may disturb sync.
       KS_eigenvalue=0.0
       occ_numbers=0.0
! add mlk
!       if (myid == 0) print*, "before i_k_point loop"
! end mlk
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

       cyclenumber = 0
       do i_k_point = 1, n_k_points_band, n_k_points

           cyclenumber = cyclenumber + 1

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
           if (.not.use_local_index) then
             call construct_overlap_scalapack(overlap_matrix)
             if (out_overlap) then
               ! This is the place to do it. Since this is Lapack, each k-point lives
               ! on a separate task. All tasks are involved in the writing.
               call output_overlap_scalapack( i_band, i_k_point )
             end if
             call construct_hamiltonian_scalapack(hamiltonian)
! add mlk
             call save_overlap_scalapack
! end mlk
           else
             if (use_load_balancing) then
               call set_full_local_ovlp(ovlp_work)
               call set_full_local_ham(my_hamiltonian_matrix)
             else
               call set_sparse_local_ovlp_scalapack(ovlp_work)
               call set_sparse_local_ham_scalapack(my_hamiltonian_matrix)
             end if
           end if

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
!                 do i_cell_1 = 1, n_cells_bvk
!                    do i_basis_2 = 1, n_basis
!                       do i_spin = 1, n_spin
!                          if(myid.eq.0)&
!                               write(use_unit,*) "fmeb", cell_index_bvk(i_cell_1,:), i_basis_1, i_basis_2, fock_matr_elem(i_cell_1,i_basis_2,i_spin)
!                       enddo
!                    enddo
!                 enddo
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
!                          write(use_unit,*) "fmenk", myid, k, i_basis_1, i_basis_2, fock_matr_elem_new_k
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

! add mlk
!       if (myid == 0) print*, "before calculate_perturbative_soc"
! end mlk

           ! Normally in SOC code, we perform the second-variational step after
           ! the eigenvalue problem has been fully solved.  Here, however, the
           ! eigenvectors for the current batch are thrown away (more
           ! precisely, overwritten) during the next batch.  So we need
           ! to do the second-variational step now and store the results.
           if (calculate_perturbative_soc) then
              ! At this point, only my_scalapack_id.eq.0 is expected to have the full set of eigenvalues
              ! We do not want to work with KS_eigenvalue directly, since we don't want to modify the
              ! scalar-relativistic results, but we also want each task to have the eigenvalues for the
              ! current k-point, hence the usage of a temporary variable
              KS_eigenvalue_temp = KS_eigenvalue(:,:,my_k_point_band)
              if (my_scalapack_id .ne. 0) then
                KS_eigenvalue_temp = 0.0d0
              end if
              call sync_vector( KS_eigenvalue_temp, n_states*n_spin, my_scalapack_comm_all )

! add mlk
!       if (myid == 0) print*, "after sync KS_eigenvalue_temp "
! end mlk
              ! A questionable design choice was made to support only LAPACK-style eigenvectors
              ! in this subroutine
              ! Here, I undo that design choice for SOC
              if (collect_eigenvectors) then
                 do i_spin = 1, n_spin
                   call convert_global_to_local_matrix_complex( n_basis, n_states, l_row, l_col, &
                        KS_eigenvector_complex(:,:,i_spin,1), mxld, mxcol, eigenvec_complex(:,:,i_spin) )
                 end do
              endif
! add mlk
!       if (myid == 0) print*, "after sync convert_global_to_local_matrix_complex"
! end mlk

              call construct_SOC_Hamiltonian(ld_soc_matrix, soc_matrix, &
                   mxld, mxcol, eigenvec, eigenvec_complex, my_k_point, &
                   mxld_soc, mxcol_soc, soc_ham)
! add mlk
!       if (myid == 0) print*, "after construct_SOC_Hamiltonian"
! end mlk

              call perform_soc_perturbation( mxld_soc, mxcol_soc, soc_ham, KS_eigenvalue_temp,&
                   KS_eigenvalue_soc_perturbed(1, 1, my_k_point_band), dummy, &
                   mxld_soc, mxcol_soc, eigenvec_soc_wf_basis )

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !mulliken decomposition starts here!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call convert_wf_basis_to_compute_basis( &
                  mxld,         mxcol,         eigenvec, eigenvec_complex,&
                  mxld_soc,     mxld_soc,      eigenvec_soc_wf_basis, &
                  mxld_soc_vec, mxcol_soc_vec, KS_eigenvector_soc_perturbed)
!                 mxld_soc_vec, mxcol_soc_vec, eigenvec_soc)

! add mlk
!       if (myid == 0) print*, "after convert_wf_basis_to_compute_basis"
! end mlk
              if(my_scalapack_id>=npcol*nprow) return

              call aims_allocate( tmp_complex, mxld_soc_vec, mxcol_soc_vec, "tmp_complx" )
              tmp_complex = (0.0d0, 0.0d0)

              call pzgemm('N','N',n_basis,n_saved_states_soc,n_basis,&
                        (1.d0,0.d0),&
                        ovlp_complex_stored, 1,1,sc_desc, &
                        KS_eigenvector_soc_perturbed, 1,1,sc_desc_soc_vec,&
                        (0.d0,0.d0),&
                        tmp_complex, 1,1,sc_desc_soc_vec)
              call pzgemm('N','N',n_basis,n_saved_states_soc,n_basis,&
                        (1.d0,0.d0),&
                        ovlp_complex_stored, 1,1,sc_desc, &
                        KS_eigenvector_soc_perturbed, 1+n_basis_soc_coll/2,1,sc_desc_soc_vec,&
                        (1.d0,0.d0),&
                        tmp_complex, 1+n_basis_soc_coll/2,1,sc_desc_soc_vec)

              do i_state = 1, n_saved_states_soc
                 if(l_col_soc_vec(i_state) == 0) cycle
                 do i_basis_1 = 1, n_basis_soc
                    if(l_row_soc_vec(i_basis_1) == 0) cycle

                    ! The SOC basis elements explicitly include spinors, but the Mulliken
                    ! analysis and basis indexing arrays do not, so conversions are needed
                    if (i_basis_1 .gt. n_basis_soc_coll/2) then
                       i_spin    = 2
                       i_basis_2 = i_basis_1 - n_basis_soc_coll/2
                    else
                       i_spin    = 1
                       i_basis_2 = i_basis_1
                    endif

                    mul_temp = dble(conjg(KS_eigenvector_soc_perturbed(l_row_soc_vec(i_basis_1),l_col_soc_vec(i_state),1,1)) * &
                      tmp_complex(l_row_soc_vec(i_basis_1),l_col_soc_vec(i_state)))
                    mulliken_decomp(basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin, cyclenumber) = &
                      mulliken_decomp(basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin, cyclenumber) + mul_temp
                 enddo
              enddo

!              if ((myid==0).and.(my_k_point_band==1).and.(i_band==1)) then
!                 open(122, file="mulliken_decomp")
!                 do j = 1, n_saved_states_soc
!                    do i = 1,10
!                       write(122,"(2F25.13)") mulliken_decomp(1,i,j,1,1)
!                    enddo
!                 enddo
!                 close(122)
!              endif

              call aims_deallocate( tmp_complex, "tmp_complex" )
! add mlk
!       if (myid == 0) print*, "after mulliken_decomp calculation"
! end mlk

              ! If the last point is calculated more than once
              ! throw away the other ones
              if (i_k_point+my_k_point-1 > n_k_points_band) then
                 KS_eigenvalue(:,:,my_k_point_band)=0
                 if(calculate_perturbative_soc) then
                   KS_eigenvalue_soc_perturbed(:, 1, my_k_point_band) = 0.0d0
                   mulliken_decomp(:, :, :, :, cyclenumber) = 0.0d0
                 end if
              end if
! add mlk
!       if (myid == 0) print*, "after perform_soc_perturbation"
! end mlk

              if (out_eigenvec) then
                 do i_spin = 1, n_spin
                    call output_complex_EV_scalapack (i_band, i_k_point, i_spin)
                 end do
              end if
              write(info_str,'(A)') ""
              call localorb_info ( info_str )

           endif ! for calculate_perturbative_soc flag, soc is calculated at this point

       enddo ! for i_k_point, outter wrapper soc and mulliken calculation


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! normal soc band output begins here!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Temporary adapt n_k_points to the given k_point distribution
       n_k_points_temp = n_k_points
       n_k_points      = n_k_points_band

       ! Synchronize eignevalues
       call sync_eigenvalues(KS_eigenvalue)

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

       ! Restore initial k_point information
       n_k_points = n_k_points_temp

       ! Now, output the SOC band structure if it's been requested
       if (calculate_perturbative_soc) then

          n_k_points_temp = n_k_points
          n_k_points      = n_k_points_band

          ! All MPI ranks within a given ScaLAPACK descriptor will have a full copy of the eigenvalues
          ! at the assigned k-point of the band, but only one can be allowed to contribute to the sync
          if ( my_scalapack_id .ne. 0 ) then
            KS_eigenvalue_soc_perturbed = 0.0d0
          end if
          !  Unlike the LAPACK version, we've already done the second-variational step when
          !  solving the eigenvalue problem, so here we need only sync, calculate occupations,
          !  then print
          call sync_vector(KS_eigenvalue_soc_perturbed, n_states_soc*n_k_points)
          ! Get occupation numbers
          call convert_sr_to_soc_environment ()
          call check_norm_p0(chemical_potential_soc, KS_eigenvalue_soc_perturbed, n_electrons, occ_numbers_soc, dummy, &
               dummy_int)
          call revert_soc_to_sr_environment ()

          !Restore initial k_point information
          n_k_points = n_k_points_temp

          ! file names for band plotting
          if(myid==0)then
             sr_suffix = ".out"
             if(i_band < 10)then
                write(file_name, '(A7,I1,A4)') 'band100', i_band,'.out'
             else if(i_band < 100)then
                write(file_name, '(A6,I2,A4)') 'band10', i_band,'.out'
             else if(i_band < 1000)then
                write(file_name, '(A5,I3,A4)') 'band1', i_band,'.out'
             else
                write(use_unit,*) 'Band output error: automatic file name does  not work with more than 999 bands!'
                stop
             end if

             open(88, file=file_name)

             !write bands
             do i_k_point = 1,  n_k_points_band
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
          endif ! for myid=0, here the normal band output ends

          lowest_un_occ = 1d100
          highest_occ    = -1d100

          do i_k_point = 1,  n_k_points_band
             do i_state = 1,  n_states_soc
               if(  occ_numbers_soc(i_state,1,i_k_point) < 0.5)then
                 lowest_un_occ = min(lowest_un_occ,  KS_eigenvalue_soc_perturbed(i_state,1,i_k_point))
               else
                 highest_occ = max(highest_occ, KS_eigenvalue_soc_perturbed(i_state,1,i_k_point))
               end if
             end do
          end do

          write(info_str,'(2X,A,I4)')  'Spin-orbit-coupled "band gap" along reciprocal space direction number: ',i_band
          call localorb_info ( info_str )

          write(info_str,'(2X,A,1X,F16.8,A)')  '|Lowest unoccupied state:',lowest_un_occ* hartree, 'eV'
          call localorb_info ( info_str )

          write(info_str,'(2X,A,1X,F16.8,A)') '|Highest occupied state:  ',highest_occ* hartree,  'eV'
          call localorb_info ( info_str )

          write(info_str,'(A)') ''
          call localorb_info ( info_str )
       endif ! for calculate_perturbative_soc, here the normal band output outter wrapper ends


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! soc-mulliken band output begins here!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !outputing the soc-mulliken needs the occupation number generated by the previous
       !part, so, it follows

       if (calculate_perturbative_soc) then ! mulliken decomp output outter wrapper

          ! prepare output file bandmlk1xxx.out for soc-mlk band output
          if(myid==0)then
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
          end if ! myid == 0
          !file preparasiion ends here

          cyclenumber = 0
          do i_k_point = 1, n_k_points_band, n_k_points
             cyclenumber = cyclenumber + 1

             !some parameter initialization
             do i = 1, n_k_points
                k_phase(:,i) = k_phase_band(:,MIN(i_k_point+i-1,n_k_points_band))
             end do
             my_k_point_band = MIN (i_k_point+my_k_point-1, n_k_points_band)
             if(use_periodic_hf)then
                k(:) = band_k_frac(my_k_point_band, i_band)
                do i_cell_bvk = 1, n_cells_bvk
                   k_phase_exx_new(i_cell_bvk) = exp((0d0,2d0)*pi*sum(k(:)*cell_index_bvk(i_cell_bvk,:)))
                enddo
             endif
             !parameter initialization ends
! add mlk
             call sync_vector(mulliken_decomp(0,1,1,1,cyclenumber),(1+l_wave_max)*n_atoms*n_saved_states_soc*2*1, &
                             my_scalapack_comm_work)
             k_band_each(:) = 0
             id_each(:) = 0

             call MPI_Gather(my_k_point_band, 1, mpi_integer, k_band_each, 1, mpi_integer, &
                             0, mpi_comm_global, mpierr) ! let the 0th core know the k points calculated on each core
             if (mpierr /= 0) then
                write(info_str,'(A)') "***ERROR:MPI_Gather Failed"
                call aims_stop(info_str)
             endif

             k_temp1 = my_k_point_band
             k_temp2 = my_k_point_band

             if (myid==0) then
                !print *, k_band_each(:)
                do i_core = 1, (n_tasks-1)
                   !print *, k_band_each(i_core+1), k_temp2
                   if (k_band_each(i_core+1) == (k_temp2 + 1)) then
                      k_temp2 = k_band_each(i_core+1)
                      id_each(k_temp2) = i_core
                   endif
                enddo
                !print *, "0th", k_temp1, k_temp2
                !print *, "0th", id_each(:)
             endif


             call bcast_integer_vector(id_each, n_k_points_band, 0)
             call bcast_integer(k_temp1, 0)
             call bcast_integer(k_temp2, 0)

             !print *, "after bcast", myid, k_temp1, k_temp2

             ! firstly, let the 0th core write out its own data
             if (myid==0) then
                !print *, "my_k_point_band", my_k_point_band
                k(:) = band_k_frac(my_k_point_band, i_band)
                write(50,'(A,I5,A,F12.8,1X,F12.8,1X,F12.8,A)') &
                     "k point number: ",my_k_point_band , ": ( ", k(1), k(2), k(3), " )"
                write(50,'(4X,A,7X,A,2X,A,1X,A, 7X,A, 7X,A, 10(9X,A3))') &
                     "State", "eigenvalue", "occ.number", "atom", "spin", "total"
                do i_state=1,n_saved_states_soc,1
                   if (i_state <= homo_id + num_orbit_output_half .and. &
                      i_state > homo_id - num_orbit_output_half) then
                      write(50,'(4X,A,2X,I7)') "State", i_state
                   endif
                   do i_atom = 1, n_atoms, 1
                      if(species_pseudoized(species(i_atom))) cycle
                      do i_spin = 1, 2, 1
                         if (i_state <= homo_id + num_orbit_output_half .and. i_state > &
                            homo_id - num_orbit_output_half) then
                            write(50,'(2X,I7,2X,F15.5,2X,F10.7,2X, I3,2X, I3,2X, F10.5, 10(2X,F10.5))') &
                            i_state, (KS_eigenvalue_soc_perturbed(i_state,1,my_k_point_band) - chemical_potential_soc)* hartree, &
                            occ_numbers_soc(i_state,1,my_k_point_band), i_atom, i_spin, &
                            sum(mulliken_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,cyclenumber)),&
                            (mulliken_decomp(i_l,i_atom,i_state,i_spin,cyclenumber), i_l=0,l_shell_max(species(i_atom)) )
                         endif
                      enddo
                   enddo
                enddo
                ! here the 0th core's own data writing out ends
             endif ! myid == 0


             ! here, transfer the data for each k point to the 0th core, one by one
             do i_k = k_temp1+1, k_temp2
                ! here, let the 0th core receive data fro each k point, one by one and write out
                if (myid==0) then

                   call receive_real_vector(mulliken_decomp(0,1,1,1,cyclenumber), &
                                      (1+l_wave_max)*n_atoms*n_saved_states_soc*2*1, id_each(i_k))

                   k(:) = band_k_frac(i_k, i_band)
                   write(50,'(A,I5,A,F12.8,1X,F12.8,1X,F12.8,A)') &
                        "k point number: ",i_k , ": ( ", k(1), k(2), k(3), " )"
                   write(50,'(4X,A,7X,A,2X,A,1X,A, 7X,A, 7X,A, 10(9X,A3))') &
                        "State", "eigenvalue", "occ.number", "atom", "spin", "total"
                   do i_state=1,n_saved_states_soc,1
                      if (i_state <= homo_id + num_orbit_output_half .and. &
                         i_state > homo_id - num_orbit_output_half) then
                         write(50,'(4X,A,2X,I7)') "State", i_state
                      endif
                      do i_atom = 1, n_atoms, 1
                         if(species_pseudoized(species(i_atom))) cycle
                         do i_spin = 1, 2, 1
                            if (i_state <= homo_id + num_orbit_output_half .and. i_state > &
                               homo_id - num_orbit_output_half) then
                               write(50,'(2X,I7,2X,F15.5,2X,F10.7,2X, I3,2X, I3,2X, F10.5,10(2X,F10.5))') &
                               i_state, (KS_eigenvalue_soc_perturbed(i_state,1,i_k) - chemical_potential_soc)* hartree, &
                               occ_numbers_soc(i_state,1,i_k), i_atom, i_spin, &
                               sum(mulliken_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,cyclenumber)),&
                               (mulliken_decomp(i_l,i_atom,i_state,i_spin,cyclenumber), i_l=0,l_shell_max(species(i_atom)) )
                            endif
                         enddo
                      enddo
                   enddo
                endif ! myid=0
                ! here, the data of i_taskth core writing out through 0th core ends
                if (myid==id_each(i_k) .and. myid /= 0) then
                   call send_real_vector(mulliken_decomp(0,1,1,1,cyclenumber),(1+l_wave_max)*n_atoms*n_saved_states_soc*2*1, 0)
                endif
             enddo
!end mlk
          enddo ! for i_k_point loop
       endif ! for calculate_perturbative_soc, here the mulliken band output outter wrapper ends
       close(50)


       if (allocated(id_each))                    call aims_deallocate(id_each,                                           "id_each")
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

   ! Allocatable variables from other subroutines (which are tracked)
   if (allocated(KS_eigenvector))              call aims_deallocate( KS_eigenvector,                           "KS_eigenvector" )
   if (allocated(KS_eigenvector_complex))      call aims_deallocate( KS_eigenvector_complex,           "KS_eigenvector_complex" )
   if (allocated(KS_eigenvalue_soc_perturbed)) call aims_deallocate( KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed" )
   if (allocated(occ_numbers_soc))             call aims_deallocate( occ_numbers_soc,                         "occ_numbers_soc" )
   if (allocated(KS_eigenvector_soc_perturbed)) call aims_deallocate( KS_eigenvector_soc_perturbed, "KS_eigenvector_soc_perturbed")
   if (allocated(mulliken_decomp))            call aims_deallocate(mulliken_decomp,                           "mulliken_decomp")
   if (allocated(k_band_each))                call aims_deallocate(k_band_each,                                   "k_band_each")

end if ! .not. real_eigenvectors

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

  end subroutine out_plot_band_soc_mulliken_scalapack
!******

