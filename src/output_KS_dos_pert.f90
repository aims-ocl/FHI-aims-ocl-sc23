!****s* FHI-aims/output_KS_dos_pert
!  NAME
!   output_KS_dos_pert
!  SYNOPSIS

subroutine output_KS_dos_pert()

  ! PURPOSE
  !
  !   The subroutine calculates the density of states perturbatively on a
  !   denser grid after the scf cycle. The dimensions of the new k-point grid
  !   are
  !
  !     n_k_points_xyz(1)*n_k_points_dos_xyz(1)
  !     n_k_points_xyz(2)*n_k_points_dos_xyz(2)
  !     n_k_points_xyz(3)*n_k_points_dos_xyz(3)
  !
  !   where n_k_points_xyz are dimensions of the old k-point grid while
  !   n_k_points_dos_xyz are the factors specified by dos_kgrid_factors
  !   keyword.
  !
  !   The new k-points are divided into batches of n_k_points k-points. Within
  !   each batch, k-points are parallelized as in the main routine.
  !
  !   The following global arrays are modified (and/or reallocated):
  !      KS_eigenvalue, KS_eigenvector_complex, k_weights, k_phase,
  !      flag_KS_k_points, and even n_k_points (if symmetry_reduced_k_points).
  !
  !   Note that the present version, same as the partial density of states,
  !   does NOT include support scalapack. This means scalapack type of
  !   eigenvectors. If normal form of eigenvectors are created from scalapack
  !   routine then this this can be called.
  !
  !   WPH:  This function could seriously use some touching up.  n_k_points_full
  !   is NOT the number of all k-points.  Rather, it is the size of the k-grid
  !   before interpolation.  The total number of k-points after interpolation,
  !   the true "full" grid, is n_k_points_dos.  Most importantly, it circumvents
  !   the ScaLAPACK architecture entirely by doing its own round-robin
  !   distribution of k-points, though it does coexist with
  !   use_scalapack.eq..true.
  !
  !   11 September 2017, WPH:  I have chosen not to update this subroutine for
  !   local indexing and load balancing, because it would require a code
  !   overhaul, and I believe having such a subroutine is the wrong approach
  !   in the first place.
  !   What should be done for k-point interpolation is that we should completely
  !   recreate the k-grid at some point in main() and re-solve the KS equation
  !   on the new k-grid using the SCF density.  Not only does this eliminate
  !   code forking, but more importantly *all* post-processed functionality uses
  !   the same k-grid; we don't run into the issue that our total DOS and
  !   projected DOS do not match up because one uses an interpolated k-grid
  !   and the other uses the SCF k-grid.  And so on.
  !   This could also be done by using the proposed restart functionality:
  !   run a calculation to self-consistency and k-point convergence, write out
  !   the density matrix, restart the calculation using a denser k-grid, and
  !   then run a "postprocessing only" mode that sets up all the necessary
  !   arrays.  This may require that we do one initialization step as well
  !   as one SCF step.
  !
  !  USES

  use dimensions,            only : calculate_perturbative_soc, n_basis, n_centers_basis_I, &
                                    n_hamiltonian_matrix_size, n_k_points, n_k_points_task, &
                                    n_states, n_spin, n_periodic, use_periodic_hf, &
                                    use_lc_wpbeh, n_states_k
  use constants,             only : bohr, hartree, pi
  use mpi_tasks,             only : myid, n_tasks, aims_stop, check_allocation
  use localorb_io,           only : localorb_info, use_unit
  use species_data,          only : l_shell_max
  use runtime_choices,       only : flag_KS_k_points, dos_alpha, dos_high_energy, dos_low_energy, &
                                    dos_n_en_points, fixed_spin_moment, flag_rel, hybrid_coeff, &
                                    out_dos, packed_matrix_format, PM_none, postscf_eigenvalues, &
                                    read_k_points_from_file, real_eigenvectors, REL_zora, &
                                    use_load_balancing, use_local_index, use_scalapack, &
                                    use_symmetry_reduced_k_grid, n_k_points_xyz, k_points_offset, &
                                    n_k_points_dos_xyz, fixed_spin_moment_electrons, use_elsi
  use pbc_lists,             only : k_phase, k_weights, k_point_list, n_cells, n_cells_bvk, &
                                    cell_index_bvk, cell_index
  use geometry,              only : recip_lattice_vector
  use mpi_tasks,             only : SYNC_OR
  use synchronize_mpi_basic, only : sync_logical, sync_vector, sync_vector_complex
  use physics,               only : KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                                    chemical_potential, chemical_potential_spin, &
                                    KS_eigenvalue_soc_perturbed, chemical_potential_soc, &
                                    rho, rho_gradient, partition_tab, hartree_potential, &
                                    occ_numbers, occ_numbers_soc, n_electrons, &
                                    hamiltonian, overlap_matrix, kinetic_density, &
                                    en_vdw, en_pot_vdw, en_xc, en_pot_xc
  use scaled_zora_transform, only : allocate_scaled_zora_transform, &
                                    integrate_scaled_zora_transf_p2, &
                                    evaluate_scaled_zora_transf_p1
  use lapack_wrapper,        only : BASIS_NON_SINGULAR, BASIS_SINGULAR, &
                                    BASIS_SINGULAR_NOT_TESTED, &
                                    clear_lapack_overlap_matrix
  use hartree_fock_p0,       only : hf_exchange_matr_real, hf_exchange_matr_real_SR, &
                                    hf_exchange_matr_complex, hf_exchange_matr_complex_SR
  use load_balancing,        only : use_batch_permutation, batch_perm, n_bp_integ
  use soc_utilities,         only : perform_soc_perturbation, &
                                    convert_sr_to_soc_environment, &
                                    revert_soc_to_sr_environment
  use timing,                only : tot_time_band_dos, tot_clock_time_band_dos, &
                                    get_times, get_timestamps
  use aims_memory_tracking,  only : aims_allocate, aims_deallocate
  use generate_aims_uuid,    only : write_aims_uuid
  use dimensions_soc,        only : n_states_soc
  use ks_wrapper,            only : solve_KS_elsi_serial
  implicit none

  !  ARGUMENTS
  !    none
  !  INPUTS
  !   o KS_eigenvalue -- Kohn-Sham eigenvalues
  !   o occ_numbers -- occupations weights of the Kohn-Sham eigenstates
  !   o n_electrons -- total number of electrons in the system
  !  OUTPUT
  !   o chemical_potential -- chemical potential
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

  real*8  :: diff_electrons
  integer :: i_counter

  real*8  :: de, en
  real*8, dimension(:,:), allocatable :: KS_dos
  real*8, dimension(:), allocatable :: KS_dos_soc
  real*8,    dimension(:,:),allocatable :: hamiltonian_w
  real*8,    dimension(:),  allocatable :: overlap_matrix_w
  complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex
  complex*16,dimension(:),  allocatable :: overlap_matrix_w_complex
  complex*16, dimension(:,:), allocatable :: fock_matrix_complex
  complex*16, dimension(:,:), allocatable :: fock_matrix_complex_SR
  real*8, dimension(:,:), allocatable :: k_point_list_old

  real*8, dimension(:,:,:),allocatable :: work_ham
  real*8, dimension(:,:),allocatable :: work_ovl

  integer :: i_lattice, i_coord
  integer :: i_full_k_points
  real*8  :: k_x, k_y, k_z

  integer :: i_k_index, n_k_index

  integer:: i_band,    i_cell_n, i_spin
  integer:: i_cell_1,  i_cell_2, i_cell_3
  integer:: i_k_point, i_state, i_e
  integer:: i_batch_x, i_batch_y, i_batch_z
  integer:: i_k_x, i_k_y, i_k_z, i_k

  real*8:: i_x, i_y, i_z
  real*8:: T1, T2, T3
  real*8 :: DE_eV

  real*8:: k(3)
  complex*16 :: k_phase_exx_old
  integer :: i_k_point_old, n_k_points_old, i_index, i_k_old, i_basis_1, i_basis_2
  logical :: real_eigenvectors_old
  real*8, dimension(:), allocatable :: k_weights_old

  logical :: singularity_check !look out for this keyword for changes VB

  integer:: i_paral,i_paral_0
  integer :: n_k_points_full, n_k_points_dos
  integer, allocatable :: parallel_index(:)

  integer :: info
  character*100 :: info_str

  character(*), parameter :: func = 'output_KS_dos_pert'

  ! Variables for periodic SOC
  complex*16, dimension(:,:), allocatable :: SOC_Hamiltonian
  real*8,dimension(:,:), allocatable :: soc_matrix
  complex*16, dimension(:,:), allocatable :: eigenvec_soc_wf_basis
  real*8 :: dummy
  integer :: i, dummy_int
  logical :: use_scalapack_save
  integer :: ld_soc_matrix

  real*8 :: dos_time = 0.d0
  real*8 :: clock_dos_time = 0.d0

  call get_timestamps(dos_time, clock_dos_time)

  i_counter = 0

  !!!VB: I removed the call to check_norm , for now.
  !      * For one thing, I do not see how occ_numbers are used in this routine.
  !      * Then, the Fermi level may be spin dependent if a spin constraint is
  !        used - the results below are then wrong.
  !      * The occupation numbers must be updated whenever the KS eigenvalues
  !        change, right away, not here
  !      * Finally, if the KS eigenvalues changed, the Fermi level should
  !        itself be updated. The routine below only recomputes occupation
  !        numbers for a fixed Fermi level, the Fermi level does not change.
  ! calculate position of Fermi level, among other things ...
  ! call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)

  ! --- Initial output
  write(info_str,'(2X,A)') &
  & 'Post-scf processing of Kohn-Sham eigenvalues on a denser k-point grid.'
  call localorb_info(info_str)

  ! One of the following conditions must be true:
  if (out_dos) then
     write(info_str,'(2X,A)') &
     & '| Using dos_kgrid_factors to calculate post-scf total density of states.'
     call localorb_info(info_str)
  end if
  if (postscf_eigenvalues) then
     write(info_str,'(2X,A)') &
     & '| Using dos_kgrid_factors to write post-scf eigenvalues to a separate file.'
     call localorb_info(info_str)
  end if

  if (read_k_points_from_file) then
     write(info_str,'(2X,A)') &
     & '| Warning: k-points are read from external file.'
     call localorb_info(info_str)
     write(info_str,'(2X,A)') &
     & '| The dimensions of the new k-point grid are equal to factors'
     call localorb_info(info_str)
     write(info_str,'(2X,A)') '| specified by the dos_kgrid_factors keyword'
     call localorb_info(info_str)
  endif

  write(info_str,'(2X,A,F18.12,A)') &
  & '| Electron chemical potential: ', chemical_potential* hartree,' eV.'
  call localorb_info(info_str)

  if (fixed_spin_moment) then
     ! In this case we work with different chemical potentials in the
     ! two spin channels to enforce a certain spin state
     do i_spin = 1, n_spin, 1
        ! safeguard in case n_spin = 1 anyway
        if (i_spin.eq.1) then
           write(info_str,'(2X,A,F18.12,A)') &
           & '| Spin-up chemical potential: ', chemical_potential_spin(i_spin)* hartree,' eV.'
           call localorb_info(info_str)
        else
           write(info_str,'(2X,A,F18.12,A)') &
           & '| Spin-dn chemical potential: ', chemical_potential_spin(i_spin)* hartree,' eV.'
           call localorb_info(info_str)
        end if
     enddo
  end if

  ! If we are writing the eigenvalues at all k-points, we need to prepare the necessary files.
  if (postscf_eigenvalues) then
     if (calculate_perturbative_soc) then
       open(86, file='Final_KS_eigenvalues.dat.no_soc')
       open(87, file='Final_KS_eigenvalues.dat')
     else
       open(86, file='Final_KS_eigenvalues.dat')
     end if

     call write_aims_uuid(info_str)
     write(86,'(A,2X,A)') '#', info_str
     write(86,'(A)') '# Plain ASCII output of the Kohn-Sham eigenvalues on a k-space grid'
     write(86,'(A,F18.12,A)') '# Output is written in electron Volts. 1 Ha = ',hartree,' eV.'
     write(86,'(A)') '# Output of reciprocal lattice vectors in Angstrom^{-1} :'
     do i_lattice = 1, n_periodic,1
        write(86,'(A,2X,F18.12,2X,F18.12,2X,F18.12)') &
          '# ',(recip_lattice_vector(i_coord,i_lattice)/bohr, i_coord = 1,3,1)
     enddo
     write(86,'(A)') '# For each k-point, we write: '
     write(86,'(A)') '#   - k-point (in units of the reciprocal lattice vectors) '
     write(86,'(A)') '#   - k-point (in units of cartesian coordinates, Angstrom^{-1}) '
     write(86,'(A)') '#   - Column format for:'
     if (n_spin.eq.1) then
       write(86,'(A)') '#     Eigenvalue number, occupation number, eigenvalue '
     else if (n_spin.eq.2) then
       write(86,'(A)') '#     Eigenvalue number, occupation number (up), eigenvalue (up),',&
                            ' occupation number (dn), eigenvalue (dn). '
     end if
     write(86,*)

     if (calculate_perturbative_soc) then
       call write_aims_uuid(info_str)
       write(87,'(A,2X,A)') '#', info_str
       write(87,'(A)') '# Plain ASCII output of the Kohn-Sham eigenvalues on a k-space grid'
       write(87,'(A,F18.12,A)') '# Output is written in electron Volts. 1 Ha = ',hartree,' eV.'
       write(87,'(A)') '# Output of reciprocal lattice vectors in Angstrom^{-1} :'
       do i_lattice = 1, n_periodic,1
          write(87,'(A,2X,F18.12,2X,F18.12,2X,F18.12)') &
            '# ',(recip_lattice_vector(i_coord,i_lattice)/bohr, i_coord = 1,3,1)
       enddo
       write(87,'(A)') '# For each k-point, we write: '
       write(87,'(A)') '#   - k-point (in units of the reciprocal lattice vectors) '
       write(87,'(A)') '#   - k-point (in units of cartesian coordinates, Angstrom^{-1}) '
       write(87,'(A)') '#   - Column format for:'
       write(87,'(A)') '#     Eigenvalue number, occupation number, eigenvalue '
       write(87,*)
     end if

  end if

  ! --- Check for singularity in original k-points

  ! before going through the additional k-points, see if any of the original
  ! k points had a singular overlap matrix. If so, we do a singularity check
  ! for all new k points below. Else, we don't bother (takes time).
  singularity_check = any(flag_KS_k_points .eq. BASIS_SINGULAR)
  call sync_logical(singularity_check, SYNC_OR)


  n_k_points_old = n_k_points
  call aims_allocate( k_weights_old, n_k_points_old, "k_weights_old" )
  k_weights_old = k_weights
  real_eigenvectors_old = real_eigenvectors

  ! --- Get k-point sizes

  n_k_points_full = product(n_k_points_xyz)
  n_k_points_dos = n_k_points_full * product(n_k_points_dos_xyz)
  call aims_allocate( k_point_list_old, n_k_points,3, "k_point_list_old" )
  k_point_list_old = k_point_list

  ! --- parallel_index

  !  This array is a mapping of k-points to the CPU assigned to that
  !  given k-point.
  call aims_allocate( parallel_index, n_k_points_dos, "parallel_index" )

  !  For regular runs, once we've calculated the interpolated
  !  (incorrectly called "perturbed") eigenvalues, we can throw
  !  away the eigenvector information, as a DOS is an integral
  !  over eigenvalues.  However, in the case of scaled ZORA or
  !  perturbative SOC (actual perturbative methods), we need to
  !  save the eigenvectors to perform the perturbation to get
  !  the perturbed interpolated (see what I mean?) eigenvalues.
  !  Thus there are forks in the logic of k-point indexing
  !  depending on whether we need to save the eigenvectors or not.

  !  We will eventually loop over all enhancement factors (these
  !  are the BATCH's), and then for a given enhancement factor
  !  loop over all n_k_points_full k-points assigned to it.  To
  !  minimize memory usage, when we need the eigenvectors we'll
  !  only deal with the k-points within the current batch.  Thus
  !  our parallelization will be over this restricted set.
  if((flag_rel.eq.REL_zora).or.calculate_perturbative_soc) then
     ! Distribute original k-points
     ! (Actually, "distribute k-points within a batch in the same
     ! way that the original k-points were distributed")
     i_k = 0
     parallel_index = -1
     do i_k_point = 1, n_k_points_full,1
        i_k = i_k +1
        parallel_index(i_k) =  MOD(i_k, n_tasks)
     end do
  !  When we are not saving eigenvectors, then we can freely
  !  parallelize over *all* k-points.
  else
     ! Distribute all k-points
     i_k = 0
     parallel_index = -1
     do i_batch_x = 1, n_k_points_dos_xyz(1)
        do i_batch_y = 1, n_k_points_dos_xyz(2)
           do i_batch_z = 1, n_k_points_dos_xyz(3)

              do i_k_point = 1, n_k_points_full
                 i_k = i_k +1

                 parallel_index(i_k) =  MOD(i_k, n_tasks)
              end do

           end do
        end do
     end do
  end if
  n_k_points = n_k_points_full
  n_k_points_task = count(myid == parallel_index)
  if (flag_rel == REL_zora.or.calculate_perturbative_soc) then
     n_k_index = n_k_points_task
  else
     n_k_index = 1
  end if

  ! --- Get hamiltonian

  call integrate_real_hamiltonian_matrix_p2(hartree_potential, &
  & rho, rho_gradient, kinetic_density, partition_tab, l_shell_max, &
  & en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

  if (calculate_perturbative_soc) then
    ! Because this particular subroutine circumvents the ScaLAPACK branch of the code
    ! entirely, we must always take the LAPACK branch of the SOC code
    ! May the circle of questionable design choices be unbroken
    if (use_scalapack) then
      write(info_str,'(2X,A)')
      call localorb_info(info_str)
      write(info_str,'(2X,A)') "** Perturbed DOS calculations will create a LAPACK-style distribution of k-points."
      call localorb_info(info_str)
      write(info_str,'(2X,A)') "** Accordingly, the SOC code will temporarily convert from ScaLAPACK (currently used) to LAPACK"
      call localorb_info(info_str)
      write(info_str,'(2X,A)') "** and will revert back to ScaLAPACK at the end of this calculation."
      call localorb_info(info_str)
      write(info_str,'(2X,A)') "** For large systems, you will run into memory issues!"
      call localorb_info(info_str)
      write(info_str,'(2X,A)')
      call localorb_info(info_str)
    end if

    ! If periodic SOC enabled, set up initial matrices
    call aims_allocate( SOC_Hamiltonian, n_states_soc, n_states_soc,            "SOC_Hamiltonian" )
    call aims_allocate( eigenvec_soc_wf_basis, n_states_soc,n_states_soc, "eigenvec_soc_wf_basis" )

    if (use_local_index.and.use_load_balancing) then
      ! Not supported as of this writing, but putting this here for future proofings
      call aims_allocate(soc_matrix, batch_perm(n_bp_integ)%n_local_matrix_size, 3, "+soc_matrix")
      use_batch_permutation = n_bp_integ
      ld_soc_matrix = batch_perm(n_bp_integ)%n_local_matrix_size
    else
      call aims_allocate(soc_matrix, n_hamiltonian_matrix_size, 3,                  "+soc_matrix")
      ld_soc_matrix = n_hamiltonian_matrix_size
    end if
    call integrate_soc_matrix (rho, hartree_potential, partition_tab, soc_matrix )

  else ! otherwise periodic SOC code will never be touched, set up dummy indices
    call aims_allocate( SOC_Hamiltonian, 1,1,             "SOC_Hamiltonian" )
    call aims_allocate( eigenvec_soc_wf_basis, 1,1, "eigenvec_soc_wf_basis" )
    call aims_allocate( soc_matrix, 1,1,                       "soc_matrix" )
  end if

  if (flag_rel.eq.REL_zora) then
     call allocate_scaled_zora_transform
     call integrate_scaled_zora_transf_p2( &
     & rho, rho_gradient, kinetic_density, hartree_potential, partition_tab, l_shell_max)
  end if

  ! --- Allocate to-be-diagonalized (complex) arrays

  call aims_allocate( hamiltonian_w_complex, n_basis*(n_basis+1)/2,n_spin,    "hamiltonian_w_complex" )
  call aims_allocate( overlap_matrix_w_complex, n_basis*(n_basis+1)/2,    "overlap_matrix_w_complex"  )
  if (use_periodic_hf) then
    call aims_allocate( fock_matrix_complex, n_basis*(n_basis+1)/2,n_spin,      "fock_matrix_complex" )
  end if
  if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
    call aims_allocate( fock_matrix_complex_SR, n_basis*(n_basis+1)/2,n_spin,"fock_matrix_complex_SR" )
  end if

  ! dummy allocations to allow -check pointers to function: hamiltonian_w is never needed here,
  ! only as part of a subroutine interface
  call aims_allocate( hamiltonian_w, 1,1,     "hamiltonian_w" )
  call aims_allocate( overlap_matrix_w, 1, "overlap_matrix_w" )


  if(packed_matrix_format == PM_none)then
     call aims_allocate( work_ham, n_centers_basis_I, n_centers_basis_I, n_spin, "work_ham" )
     call aims_allocate( work_ovl, n_centers_basis_I, n_centers_basis_I,         "work_ovl" )
  else
     ! arrays are never touched. Only dummy allocation statements.
     call aims_allocate( work_ham, 1, 1, 1, "work_ham" )
     call aims_allocate( work_ovl, 1, 1,    "work_ovl" )
  end if


  real_eigenvectors = .false.
  if (allocated(KS_eigenvector_complex)) call aims_deallocate( KS_eigenvector_complex, "KS_eigenvector_complex" )
  call aims_allocate(KS_eigenvector_complex, n_basis, n_states, n_spin, n_k_index,     "KS_eigenvector_complex" )
!  write(use_unit,*) "n_k_index ", n_k_index, " myid ", myid

  if (use_symmetry_reduced_k_grid) then
     if (allocated(k_phase)) deallocate(k_phase)
     allocate(k_phase(n_cells, n_k_points_full), stat=info)
     call check_allocation(info, 'k_phase', func)

     if (allocated(k_weights)) deallocate(k_weights)
     allocate(k_weights(n_k_points_full), stat=info)
     call check_allocation(info, 'k_weights', func)

     if (allocated(n_states_k)) deallocate(n_states_k)
     allocate(n_states_k(n_k_points_full), stat=info)
     call check_allocation(info, 'n_states_k', func)

     if (allocated(k_point_list)) deallocate(k_point_list)
     allocate(k_point_list(n_k_points_full,3))
     call check_allocation(info, 'k_point_list', func)

     if (allocated(flag_KS_k_points)) deallocate(flag_KS_k_points)
     allocate(flag_KS_k_points(n_k_points_full), stat=info)
     call check_allocation(info, 'flag_KS_k_points', func)
     flag_KS_k_points = BASIS_NON_SINGULAR  ! Corrected if needed.

     if (allocated(KS_eigenvalue)) deallocate(KS_eigenvalue)
     allocate(KS_eigenvalue(n_states, n_spin, n_k_points_full), stat=info)
     call check_allocation(info, 'KS_eigenvalue', func)

     if (allocated(occ_numbers)) deallocate(occ_numbers)
     allocate(occ_numbers(n_states, n_spin, n_k_points_full), stat=info)
     call check_allocation(info, 'occ_numbers', func)
  end if

  k_weights = 1.0d0 / dble(n_k_points_dos)

  call aims_allocate( KS_dos, n_spin, dos_n_en_points, "KS_dos" )
  if (calculate_perturbative_soc) then
    call aims_allocate( KS_dos_soc, dos_n_en_points, "KS_dos_soc" )
  else
    ! Spin is not a good quantum number in SOC-perturbed version, so only
    ! one "spin" channel exists
    call aims_allocate( KS_dos_soc, 1, "KS_dos_soc" )
  end if
  KS_dos = 0d0
  KS_dos_soc = 0d0
  de = (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)

  ! --- Main loop

  ! counter for total number of k-points (for postscf_eigenvalues output)
  i_full_k_points = 0

  i_paral = 0
  !   The main do loop's of the function (the BATCH loops) are over the enhancement
  !   factors.  For each enhancement factor, we assign a number of k-points equal to
  !   the size of original k-grid, and they're done as a nested set of do-loops within
  !   the main do loop's.  Note we haven't done any parallelization yet.
  BATCH_X: do i_batch_x = 1, n_k_points_dos_xyz(1),1
     BATCH_Y: do i_batch_y = 1, n_k_points_dos_xyz(2),1
        BATCH_Z: do i_batch_z = 1, n_k_points_dos_xyz(3),1

           write(info_str,'(2X,A,1X,I5,1X,I5,1X,I5)') &
           & 'Processing batch of k-points:', i_batch_x, i_batch_y, i_batch_z
           call localorb_info(info_str)

           ! Restart if zora.
           ! For ZORA and periodic SOC, we're saving eigenvectors, so ignore the
           ! "global" indexing by resetting i_paral.  For the normal case,
           ! however, do not reset this, so that i_paral will continue to be
           ! incremented across batches.
           if (flag_rel.eq.REL_zora.or.calculate_perturbative_soc) i_paral = 0
           if (singularity_check) flag_KS_k_points = BASIS_SINGULAR_NOT_TESTED

           ! --- k_phase
           !  Here we calculate the phase factors for each of the k-points
           !  assigned to the present batch (the given enhancement factor)
           !  Still no parallelization.
           i_k_point = 0
           do i_k_x = 1, n_k_points_xyz(1)
              i_x = (i_batch_x - 1)*n_k_points_xyz(1) + i_k_x
              do i_k_y = 1, n_k_points_xyz(2)
                 i_y = (i_batch_y - 1)*n_k_points_xyz(2) + i_k_y
                 do i_k_z = 1, n_k_points_xyz(3)
                    i_z = (i_batch_z - 1)*n_k_points_xyz(3) + i_k_z

                    i_k_point = i_k_point + 1
                    if (i_k_point > n_k_points_full) then
                       call aims_stop('Invalid i_k_point', func)
                    end if

                    k_point_list(i_k_point,1) = (i_x-1)/(n_k_points_xyz(1)*n_k_points_dos_xyz(1)) + k_points_offset(1)
                    k_point_list(i_k_point,2) = (i_y-1)/(n_k_points_xyz(2)*n_k_points_dos_xyz(2)) + k_points_offset(2)
                    k_point_list(i_k_point,3) = (i_z-1)/(n_k_points_xyz(3)*n_k_points_dos_xyz(3)) + k_points_offset(3)

                    T1 = 2*pi*k_point_list(i_k_point,1)
                    T2 = 2*pi*k_point_list(i_k_point,2)
                    T3 = 2*pi*k_point_list(i_k_point,3)

                    do i_cell_n = 1, n_cells, 1

                       i_cell_1 = cell_index(i_cell_n,1)
                       i_cell_2 = cell_index(i_cell_n,2)
                       i_cell_3 = cell_index(i_cell_n,3)

                       k_phase( i_cell_n, i_k_point) &
                            = ((exp((0,1)*T1*i_cell_1))* (exp((0,1)*T2*i_cell_2)) * (exp((0,1)*T3*i_cell_3)))

                    enddo
                 enddo
              enddo
           enddo

           ! --- Diagonalization

           !  We need to save this variable because we're going to use it twice.
           !  First, we're going to use it as a running index to globally index the KS
           !  eigenvalues for each k-point within this batch.  We are then going
           !  to run through all those k-points within the batch a second time, this
           !  time calculating their contribution to the DOS.  (This will
           !  require us to reset the global index to where it was before we
           !  calculated all the KS eigenvalues, hence the need to save it.)
           i_paral_0 = i_paral

           !  Now we're going to solve the KS equations for each k-point (on a
           !  grid the size of the original k-grid) that has been mapped to the
           !  present batch.
           !  What makes this hard to follow is that there are THREE
           !  parallelization schemes (that haven't been documented, naturally).
           !  The first is completely internal to this loop.  Because this loop contains
           !  the same number of k-points as the original SCF cycle, we can use
           !  that parallelization scheme if we want to parallelize within this
           !  loop.  In this scheme, the current k-point is indexed by
           !  i_k_point, and the parallelization scheme is the familiar
           !  myid.eq.MOD... that you should have encountered repeatedly in AIMS
           !  by now.  This parallelization is only used when construction the
           !  Fock matrix, and is used in this way regardless of the usage of
           !  i_paral.
           !  However, we also may need to be able to parallelize over *all*
           !  k-points and need a way to identify the current k-point relative
           !  to all other k-points.   Alternatively, we might not want this,
           !  and instead want a way to continue to parallelize only within the
           !  batch (i.e. when saving eigenvectors.)  i_paral allows us to do
           !  both.  When we do not reset i_paral between batches, it's a running
           !  counter to indicate which k-point this is relative to
           !  *all* k-points, and when we do, it's an indexing only within this
           !  batch.   In both cases, the parallel_index matrix indicates which
           !  CPU the current k-point is mapped to.
           i_k = 0
           KS_eigenvalue = 0.0d0
           do i_k_point = 1, n_k_points_full
              i_paral = i_paral + 1

              !  Here parallelize according to original SCF's parallelization
              !  scheme
              if(use_periodic_hf)then
                 fock_matrix_complex = (0d0,0d0)
                 if (use_lc_wpbeh .and. hybrid_coeff /= 0.d0) then
                   fock_matrix_complex_SR = (0d0,0d0)
                 end if
                 k(:) = k_point_list(i_k_point,:)
                 do i_spin = 1, n_spin
                    i_k_old = 0
                    do i_k_point_old = 1, n_k_points_old
                       if (myid.eq.  MOD(i_k_point_old, n_tasks) .and. myid<= n_k_points_old ) then
                          i_k_old = i_k_old + 1
                          do i_cell_n = 1, n_cells_bvk
                             k_phase_exx_old = exp((0d0,2d0)*pi*(sum(k(:)* dble(cell_index_bvk(i_cell_n,:)))-&
                                  sum(k_point_list_old(i_k_point_old,:)* dble(cell_index_bvk(i_cell_n,:)))))

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

              !  From here, parallelize according to new dual scheme
              if (myid /= parallel_index(i_paral)) cycle

              !  This indexes the eigenvector within the given CPU, for
              !  cases where we need the eigenvector
              i_k = i_k + 1

              ! Make sure we always recreate the overlap matrix
              call clear_lapack_overlap_matrix


              call construct_hamiltonian_and_ovl(&
              & hamiltonian, overlap_matrix, &
              & hamiltonian_w, overlap_matrix_w, &
              & hamiltonian_w_complex, overlap_matrix_w_complex, &
              & i_k_point, work_ham, work_ovl)

              !Add Fock matrix to Hamiltonian
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

              ! solve KS-equations for both spins
              do i_spin = 1, n_spin, 1
                 ! calculate the eigenvalues and eigenvectors

                 !  We set i_k_index = 1 in the "regular" case, because
                 !  we don't care about the eigenvectors, and thus want
                 !  to continuously overwrite the same portion of memory
                 !  in order to reduce memory overhead.
                 if(flag_rel.eq.REL_zora.or.calculate_perturbative_soc) then
                    i_k_index = i_k
                 else
                    i_k_index = 1
                 end if

                 if (use_elsi .and. .not. use_scalapack) then
                    call solve_KS_elsi_serial(hamiltonian_w_complex(:,i_spin),&
                         overlap_matrix_w_complex,&
                         KS_eigenvalue(:,i_spin,i_k_point),&
                         KS_eigenvector_complex(:,:,i_spin,i_k_index),i_spin,&
                         i_k_point)
                 else
                    call improve_complex_eigenfunctions(&
                         overlap_matrix_w_complex,&
                         hamiltonian_w_complex(:,i_spin),&
                         KS_eigenvalue(:,i_spin,i_k_point),&
                         KS_eigenvector_complex(:,:,i_spin,i_k_index),i_k_point)
                 end if
              end do
           end do

           call sync_vector(KS_eigenvalue, n_states*n_spin*n_k_points_full)
           if (flag_rel.eq.REL_zora) then
              call evaluate_scaled_zora_transf_p1(KS_eigenvector, &
              &                          KS_eigenvector_complex, KS_eigenvalue)
           end if

           ! First, calculate the scalar-relativistic eigenvalues

           ! Update occupation numbers of the eigenstates on present k-points only.
           ! Note that the chemical_potential value is not changed by this routine.
           ! The diff_electrons and i_counter values have no meaning (or relevance)
           ! in the present context.
           if (.not.fixed_spin_moment) then
              call check_norm_p0( chemical_potential, KS_eigenvalue, &
              &      n_electrons, occ_numbers, diff_electrons, i_counter)
           else
              ! fixed_spin_moment is a special case - two Fermi levels used during the enforcement
              do i_spin = 1, n_spin,1
                 call check_norm_periodic_v2 &
                 ( chemical_potential_spin(i_spin), KS_eigenvalue(1:n_states,i_spin,1:n_k_points), &
                   fixed_spin_moment_electrons(i_spin), occ_numbers(1:n_states,i_spin,1:n_k_points), &
                   diff_electrons, i_counter, i_spin )
              enddo
           end if

           !  We've now calculated all the eigenvalues and occupation numbers
           !  for the k-points assigned to this batch.  Now we're going to go
           !  through them all over again, this time assigning them to the
           !  DOS.  Thus we need to restart the indexing to point to the first
           !  k-point in the batch.
           i_paral = i_paral_0

           if (out_dos) then
           ! Only compute the DOS if it was actually requested with proper parameters
              do i_k_point = 1, n_k_points_full
                 i_paral = i_paral + 1
                 if (myid /= parallel_index(i_paral)) cycle

                 do i_spin = 1, n_spin
                   do i_e = 1, dos_n_en_points
                     en= dos_low_energy + dble(i_e-1)*de
                     do i_state = 1, n_states, 1
                        DE_eV = &
                        & en - KS_eigenvalue(i_state, i_spin, i_k_point)*hartree
                        KS_dos(i_spin, i_e) = KS_dos(i_spin, i_e) &
                        & + k_weights(i_k_point) &
                        & *exp(-DE_eV**2/2.0/dos_alpha**2) /sqrt(2*pi)/dos_alpha
                     enddo
                   enddo
                 end do
              enddo
           end if ! out_dos

           if (postscf_eigenvalues) then
           ! Only output all eigenvalues to a separate file if this was actually requested.
           ! Note that this is a lot of output. We write out everything on task 0 (eigenvalues
           ! were synced further up) and we write out in a somewhat odd order since the
           ! eigenvalue are computed on separate batches of k-points.
              do i_k_point = 1, n_k_points_full
                 i_full_k_points = i_full_k_points + 1

                 write(86,'(A,2X,I8)') 'k-point number: ', i_full_k_points
                 write(86,'(A,3(2X,F18.12))') 'k-point in recip. lattice units: ', &
                    (k_point_list(i_k_point,i_coord),i_coord= 1,3,1)
                 k_x = k_point_list(i_k_point,1) * recip_lattice_vector(1,1) + &
                       k_point_list(i_k_point,2) * recip_lattice_vector(1,2) + &
                       k_point_list(i_k_point,3) * recip_lattice_vector(1,3)
                 k_y = k_point_list(i_k_point,1) * recip_lattice_vector(2,1) + &
                       k_point_list(i_k_point,2) * recip_lattice_vector(2,2) + &
                       k_point_list(i_k_point,3) * recip_lattice_vector(2,3)
                 k_z = k_point_list(i_k_point,1) * recip_lattice_vector(3,1) + &
                       k_point_list(i_k_point,2) * recip_lattice_vector(3,2) + &
                       k_point_list(i_k_point,3) * recip_lattice_vector(3,3)
                 write(86,'(A,3(2X,F18.12))') 'k-point in cartesian units (A^{-1}): ', &
                    k_x/bohr, k_y/bohr, k_z/bohr
                 do i_state = 1, n_states, 1
                    if (n_spin.eq.1) then
                      write(86,'(I8,2X,F15.12,2X,F22.12)') &
                      i_state, occ_numbers(i_state,1,i_k_point),KS_eigenvalue(i_state,1,i_k_point)*hartree
                    else if (n_spin.eq.2) then
                      write(86,'(I8,2(2X,F15.12,2X,F22.12))') &
                      i_state, occ_numbers(i_state,1,i_k_point),KS_eigenvalue(i_state,1,i_k_point)*hartree,&
                      occ_numbers(i_state,2,i_k_point),KS_eigenvalue(i_state,2,i_k_point)*hartree
                    end if
                 enddo ! i_states
                 write(86,*)
              enddo
           end if ! postscf_eigenvalues

           ! Then, if requested, the spin-orbit-coupled eigenvalues

           if (calculate_perturbative_soc) then
             ! Because we will never store the SOC-perturbed eigenvectors, here we
             ! compute and save all eigenvalues in the second-variational window,
             ! since they're cheap
             if (allocated(KS_eigenvalue_soc_perturbed)) &
                  call aims_deallocate(KS_eigenvalue_soc_perturbed,                             "KS_eigenvalue_soc_perturbed" )
             call aims_allocate( KS_eigenvalue_soc_perturbed, n_states_soc, 1, n_k_points_full, "KS_eigenvalue_soc_perturbed" )
             KS_eigenvalue_soc_perturbed = 0.0d0

             if (allocated(occ_numbers_soc)) call aims_deallocate(occ_numbers_soc,                          "occ_numbers_soc" )
             call aims_allocate( occ_numbers_soc, n_states_soc, 1, n_k_points_full,                         "occ_numbers_soc" )

             ! Calculate the SOC perturbation on each k-point.  This block is the
             ! core functionality of calculate_second_variational_soc,
             ! with the exception of no outputting and lack of determination of
             ! Fermi level, as this has already been set in that main SOC routine.

             ! In other SOC calculations, this is where the setup of
             ! the k-point indexing would go (my_k_point).  However, here we are
             ! using a completely different k-point indexing scheme, and this
             ! boilerplate code must be tweaked to fit with it.

             ! We need to go through the k-points within this batch a second time,
             ! this time perturbing them, and so we reset the counter back to the
             ! beginning of the batch.
             i_paral = i_paral_0

             i_k = 0

             ! Because this particular subroutine circumvents the ScaLAPACK branch of the code
             ! entirely, we must always take the LAPACK branch of the SOC code
             ! May the circle of questionable design choices be unbroken
             use_scalapack_save = use_scalapack
             use_scalapack = .false.

             do  i_k_point = 1,  n_k_points_full
               i_paral = i_paral + 1
               if (myid /= parallel_index(i_paral)) cycle
               i_k = i_k + 1

               ! Added to make code more symmetric with ZORA case, but not needed
               i_k_index = i_k

               ! Apply perturbative theory at each k-point
               ! Eigenvectors are assumed to be imaginary, so no switch
               ! statement.

               call construct_SOC_Hamiltonian(ld_soc_matrix, soc_matrix, &
                    n_basis, n_states, KS_eigenvector(:,:,:,1), KS_eigenvector_complex(:,:,:,i_k_index), i_k_point, &
                    n_states_soc, n_states_soc, SOC_Hamiltonian)

               call perform_soc_perturbation( n_states_soc, n_states_soc, SOC_Hamiltonian, KS_eigenvalue(1,1,i_k_point),&
                    KS_eigenvalue_soc_perturbed(1, 1, i_k_point), dummy, &
                    n_states_soc, n_states_soc, eigenvec_soc_wf_basis )
             end do

             ! Now back to our regularly scheduled programming
             use_scalapack = use_scalapack_save

             call sync_vector(KS_eigenvalue_soc_perturbed, n_states_soc*n_k_points_full)


             ! Get occupation numbers
             call convert_sr_to_soc_environment ()
             call check_norm_p0(chemical_potential_soc, KS_eigenvalue_soc_perturbed, &
                  n_electrons, occ_numbers_soc, dummy, dummy_int)
             call revert_soc_to_sr_environment ()

             !  We've now calculated all the perturbed eigenvalues and occupation numbers
             !  for the k-points assigned to this batch.  Now we're going to go
             !  through the batch a third (and final) time, this time assigning the energies
             !  to the DOS.
             i_paral = i_paral_0

             if (out_dos) then
             ! Only compute the DOS if it was actually requested with proper parameters
                do i_k_point = 1, n_k_points_full
                   i_paral = i_paral + 1
                   if (myid /= parallel_index(i_paral)) cycle

                   do i_e = 1, dos_n_en_points
                     en= dos_low_energy + dble(i_e-1)*de
                     do i_state = 1, n_states_soc, 1
                        DE_eV = &
                        & en - KS_eigenvalue_soc_perturbed(i_state, 1, i_k_point)*hartree
                        KS_dos_soc(i_e) = KS_dos_soc(i_e) &
                        & + k_weights(i_k_point) &
                        & *exp(-DE_eV**2/2.0/dos_alpha**2) /sqrt(2*pi)/dos_alpha
                     enddo
                   enddo
                enddo
             end if ! out_dos

             if (postscf_eigenvalues) then
             ! Only output all eigenvalues to a separate file if this was actually requested.
             ! Note that this is a lot of output. We write out everything on task 0 (eigenvalues
             ! were synced further up) and we write out in a somewhat odd order since the
             ! eigenvalue are computed on separate batches of k-points.
                do i_k_point = 1, n_k_points_full
                   i_full_k_points = i_full_k_points + 1

                   write(87,'(A,2X,I8)') 'k-point number: ', i_full_k_points
                   write(87,'(A,3(2X,F18.12))') 'k-point in recip. lattice units: ', &
                      (k_point_list(i_k_point,i_coord),i_coord= 1,3,1)
                   k_x = k_point_list(i_k_point,1) * recip_lattice_vector(1,1) + &
                         k_point_list(i_k_point,2) * recip_lattice_vector(1,2) + &
                         k_point_list(i_k_point,3) * recip_lattice_vector(1,3)
                   k_y = k_point_list(i_k_point,1) * recip_lattice_vector(2,1) + &
                         k_point_list(i_k_point,2) * recip_lattice_vector(2,2) + &
                         k_point_list(i_k_point,3) * recip_lattice_vector(2,3)
                   k_z = k_point_list(i_k_point,1) * recip_lattice_vector(3,1) + &
                         k_point_list(i_k_point,2) * recip_lattice_vector(3,2) + &
                         k_point_list(i_k_point,3) * recip_lattice_vector(3,3)
                   write(87,'(A,3(2X,F18.12))') 'k-point in cartesian units (A^{-1}): ', &
                      k_x/bohr, k_y/bohr, k_z/bohr
                   do i_state = 1, n_states_soc, 1
                     write(87,'(I8,2X,F15.12,2X,F22.12)') &
                          i_state, occ_numbers_soc(i_state,1,i_k_point), &
                          KS_eigenvalue_soc_perturbed(i_state,1,i_k_point)*hartree
                   enddo ! i_states
                   write(87,*)
                enddo
             end if ! postscf_eigenvalues

          end if

        enddo BATCH_Z
     enddo BATCH_Y
  enddo BATCH_X

  if (out_dos) then
  ! Only write the DOS if it was actually requested.
     call sync_vector(KS_dos, n_spin*dos_n_en_points)
     if (calculate_perturbative_soc) then
       call sync_vector(KS_dos_soc, dos_n_en_points)
     end if

     if(myid.eq.0)then
        ! arrange output for DOS: different for spin polarized and normal
        ! calculations

        if (calculate_perturbative_soc) then
          write(use_unit,'(2X,A)') '| writing scalar-relativistic perturbative DOS (shifted by electron chemical potential) &
               &to file KS_dos_total.dat.no_soc'
          write(use_unit,'(2X,A)') '| writing scalar-relativistic perturbative DOS (raw data) to file KS_dos_total_raw.dat.no_soc'
          open(88, file='KS_DOS_total.dat.no_soc')
          open(89, file='KS_DOS_total_raw.dat.no_soc')
        else
          write(use_unit,'(2X,A)') '| writing perturbative DOS (shifted by electron chemical potential) to file KS_dos_total.dat'
          write(use_unit,'(2X,A)') '| writing perturbative DOS (raw data) to file KS_dos_total_raw.dat'
          open(88, file='KS_DOS_total.dat')
          open(89, file='KS_DOS_total_raw.dat')
        end if

        if (n_spin.eq.2) then
           write(88,'(2A)') '# total density of states output by FHI-aims, ',&
                'for a spin polarized system'
           write(88,'(A,F15.6,A)') '# The energy reference for this output is the average chemical potential, mu = ', &
                chemical_potential*hartree, ' eV'
           call write_aims_uuid(info_str)
           write(88,'(A,2X,A)') '#', info_str
           write(88,'(A)') '#    Energy (eV)      DOS(spin up)   DOS(spin down)'
           write(89,'(2A)') '# total density of states output by FHI-aims, ',&
                'for a spin polarized system'
           write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
           call write_aims_uuid(info_str)
           write(89,'(A,2X,A)') '#', info_str
           write(89,'(A)') '#    Energy (eV)      DOS(spin up)   DOS(spin down)'
           do i_e = 1, dos_n_en_points ,1
              en = dos_low_energy + dble(i_e-1)*de
              write(88,'(3F16.8)') en-chemical_potential*hartree, KS_dos(1,i_e), KS_dos(2,i_e)
              write(89,'(3F16.8)') en, KS_dos(1,i_e), KS_dos(2,i_e)
           enddo
        else
           write(88,'(2A)') '# total density of states output by FHI-aims '
           write(88,'(A,F15.6,A)') '# The energy reference for this output is the average chemical potential, mu = ', &
                chemical_potential*hartree, ' eV'
           write(88,'(A)') '#    Energy (eV)      DOS '
           write(89,'(2A)') '# total density of states output by FHI-aims '
           write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
           write(89,'(A)') '#    Energy (eV)      DOS '
           do i_e = 1, dos_n_en_points ,1
              en = dos_low_energy + dble(i_e-1)*de
              ! WPH:  Added factors of two to account for both spin channels
              write(88,'(3F20.8)') en-chemical_potential*hartree, 2*KS_dos(1,i_e)
              write(89,'(3F20.8)') en, 2*KS_dos(1,i_e)
           enddo
        end if

        close(88)
        close(89)

        if (calculate_perturbative_soc) then
          write(use_unit,'(2X,A)') '| writing spin-orbit-coupled perturbative DOS (shifted by electron chemical potential) &
               &to file KS_dos_total.dat'
          write(use_unit,'(2X,A)') '| writing spin-orbit-coupled perturbative DOS (raw data) to file KS_dos_total_raw.dat'
          open(88, file='KS_DOS_total.dat')
          open(89, file='KS_DOS_total_raw.dat')

          write(88,'(2A)') '# total density of states output by FHI-aims ',&
               'with perturbative spin-orbit coupling'
          write(88,'(A,F15.6,A)') '# The energy reference for this output is the average chemical potential, mu = ', &
               chemical_potential_soc*hartree, ' eV'
          write(88,'(A)') '#    Energy (eV)      DOS '
          write(89,'(2A)') '# total density of states output by FHI-aims '
          write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
          write(89,'(A)') '#    Energy (eV)      DOS '
          do i_e = 1, dos_n_en_points ,1
            en = dos_low_energy + dble(i_e-1)*de
            write(88,'(3F20.8)') en-chemical_potential_soc*hartree, KS_dos_soc(i_e)
            write(89,'(3F20.8)') en, KS_dos_soc(i_e)
          enddo

          close(88)
          close(89)
        end if
     end if ! (myid.eq.0)
  end if ! (out_dos==.true.)

  ! If we are writing the eigenvalues at all k-points, we need to prepare the necessary files.
  if (postscf_eigenvalues) then
     close(86)
     if (calculate_perturbative_soc) then
       close(87)
     end if
  end if

  ! NR Old k-grid (without multipliers) is for example needed for a subsequent band structure calculation.
  k_point_list = k_point_list_old
  n_k_points = n_k_points_old
  k_weights = k_weights_old
  real_eigenvectors = real_eigenvectors_old

  ! Allocatable variables from this subroutine
  if (allocated(KS_dos))                      call aims_deallocate(KS_dos,                                            "KS_dos" )
  if (allocated(KS_dos_soc))                  call aims_deallocate(KS_dos_soc,                                    "KS_dos_soc" )
  if (allocated(hamiltonian_w))               call aims_deallocate(hamiltonian_w,                              "hamiltonian_w" )
  if (allocated(overlap_matrix_w))            call aims_deallocate(overlap_matrix_w,                        "overlap_matrix_w" )
  if (allocated(hamiltonian_w_complex))       call aims_deallocate(hamiltonian_w_complex,              "hamiltonian_w_complex" )
  if (allocated(overlap_matrix_w_complex))    call aims_deallocate(overlap_matrix_w_complex,        "overlap_matrix_w_complex" )
  if (allocated(fock_matrix_complex))         call aims_deallocate(fock_matrix_complex,                  "fock_matrix_complex" )
  if (allocated(fock_matrix_complex_SR))      call aims_deallocate(fock_matrix_complex_SR,            "fock_matrix_complex_SR" )
  if (allocated(k_point_list_old))            call aims_deallocate(k_point_list_old,                        "k_point_list_old" )
  if (allocated(work_ham))                    call aims_deallocate(work_ham,                                        "work_ham" )
  if (allocated(work_ovl))                    call aims_deallocate(work_ovl,                                        "work_ovl" )
  if (allocated(k_weights_old))               call aims_deallocate(k_weights_old,                              "k_weights_old" )
  if (allocated(parallel_index))              call aims_deallocate(parallel_index,                            "parallel_index" )
  if (allocated(SOC_Hamiltonian))             call aims_deallocate(SOC_Hamiltonian,                          "SOC_Hamiltonian" )
  if (allocated(soc_matrix))                  call aims_deallocate(soc_matrix,                                    "soc_matrix" )
  if (allocated(eigenvec_soc_wf_basis))       call aims_deallocate(eigenvec_soc_wf_basis,              "eigenvec_soc_wf_basis" )

  ! Allocatable variables from other subroutines (which are tracked)
  if (allocated(KS_eigenvector_complex))      call aims_deallocate( KS_eigenvector_complex,           "KS_eigenvector_complex" )
  if (allocated(KS_eigenvalue_soc_perturbed)) call aims_deallocate( KS_eigenvalue_soc_perturbed, "KS_eigenvalue_soc_perturbed" )
  if (allocated(occ_numbers_soc))             call aims_deallocate( occ_numbers_soc,                         "occ_numbers_soc" )

  call get_times(dos_time, clock_dos_time, tot_time_band_dos, tot_clock_time_band_dos, .true.)

  write(info_str,'(A)')
  call localorb_info(info_str)
  write(info_str,'(2X,A)') &
       "Perturbative DOS                                        :  max(cpu_time) wall_clock(cpu1)"
  call localorb_info(info_str)
  write(info_str, "(2X, A,F15.3,F17.3)")  "| Total Time                                            :", &
       dos_time, clock_dos_time
  call localorb_info(info_str)
  write(info_str,'(2X,A)') "------------------------------------------------------------"
  call localorb_info(info_str)

end subroutine output_KS_dos_pert
!******
