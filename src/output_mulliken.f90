!****h* FHI-aims/mulliken
!*  NAME
!*    mulliken
!*  SYNOPSIS
module mulliken
!*  PURPOSE
!*    This module calculates and outputs the projected charge and DOS based on a Mulliken
!*    decomposition.
!*  USES
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University), based on the output_mulliken subroutine
!*  NOTES
!*    o Despite the name of this module being "mulliken" (for legacy reasons),
!*      almost nothing in this module depends on a Mulliken decomposition having
!*      been done, and we could extend the subroutines in this module to any method 
!*      which projects states onto atoms or atom-centered l-channels.
!*    o What should be done IMO is this module should be split into a "Mmlliken" module
!*      which contains only the Mulliken decomposition construction and output_mulliken,
!*      and a "proj_output" module which contains the output subroutines slightly 
!*      rewritten to eliminate explicit mention of a Mulliken decomposition, though with 
!*      a comment line passed in with the type of decomposition used so that people 
!*      reading the output file know which type of decomposition was used.  Right now 
!*      this is a bit anal even for my tastes, since we only do Mulliken, but if we ever 
!*      get around to implementing Bader analysis or decide to weave Voronoi polyhedra 
!*      through output_proj_charge, this should be done.
!*    o write_species_proj_dos and write_atom_proj_dos are very similar.  The
!*      core logic may be fusable into a third subroutine, but all the baggage
!*      re: filenaming prevents them from being cleanly fused into one clean
!*      subroutine
!*    o Throughout the commenting in this module, I've been a bit reckless in mixing
!*      "collinear" with "scalar-relativistic" and "non-collinear" with "spin-orbit-coupled".
!*      My gut instinct is that the former designation holds and that the "scalar-relativistic"
!*      case holds for any collinear calculation and the "spin-orbit-coupled" case holds for
!*      any (two-component) non-collinear calculation, but I haven't proven this yet.
!*  TODO
!*    o White space should be systematized.  output_mulliken had a lot of previous authors, 
!*      and it shows.
!*    o I am of the opinion that output subroutines shouldn't be responsible for opening files.
!*      They should be passed a unit and write to that unit, which may be a file or a system
!*      buffer.  The subroutine shouldn't care.  This is a minor quibble that doesn't matter 
!*      right now, so I've kept the I/O behavior of the output subroutines unaltered.
!*    o construct_mulliken_scalapack should be moved back into the module.  I in general feel 
!*      that scalapack_wrapper does too much.
!*  HISTORY
!*    November 2017     - Created; largely the same as the original output_mulliken subroutine,
!*                        but split up into module subroutines to better differentiate what is 
!*                        being done where and aid in reusability
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!*  SOURCE
  public

  real*8, dimension(:,:,:,:), allocatable :: mulliken_decomp_total ! Mulliken decomposition, summed over l channels
                                                                   ! Only used for locally constrained DFT
  ! Options for the "mode" variable
  ! As of this writing (24 November 2017), the two option specify whether scalar-relativistic or
  ! spin-orbit-coupled outputs are desired.
  integer, parameter                      :: MULLIKEN_SR  = 0
  integer, parameter                      :: MULLIKEN_SOC = 1

  ! Main driver subroutine
  public :: output_mulliken
  ! Mulliken decomposition construction subroutines
  public :: construct_mulliken_decomp
  public :: construct_mulliken_decomp_soc
  public :: construct_mulliken_decomp_total
  ! Output subroutines
  public :: output_proj_charge
  public :: write_decomp_file
  public :: write_species_proj_dos
  public :: write_atom_proj_dos

contains

  !****f* mulliken/output_mulliken
  !*  NAME
  !*    output_mulliken
  !*  SYNOPSIS
  subroutine output_mulliken &
       ( n_rows, n_cols, n_spin_states, KS_eigenvector, KS_eigenvector_complex, &
         overlap_matrix, n_states_in, KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
         filename, mode )
  !*  PURPOSE
  !*    The subroutine makes and prints out Mulliken analysis and l-projected density of states.
  !*  USES
    use dimensions,            only : n_hamiltonian_matrix_size, n_atoms, calculate_perturbative_soc, &
                                      l_wave_max, n_k_points, n_k_points_task, use_constraint
    use mpi_tasks,             only : aims_stop
    use localorb_io,           only : localorb_info
    use runtime_choices,       only : out_atom_dos, out_l_proj_dos, out_mulliken, &
                                      out_l_proj_dos_tetrahedron, out_atom_dos_tetrahedron
    use aims_memory_tracking,  only : aims_allocate, aims_deallocate
    implicit none
  !*  ARGUMENTS
    integer,                                                               intent(in) :: n_rows
    integer,                                                               intent(in) :: n_cols
    integer,                                                               intent(in) :: n_spin_states
    real*8, dimension(n_rows, n_cols, n_spin_states, n_k_points_task),     intent(in) :: KS_eigenvector
    complex*16, dimension(n_rows, n_cols, n_spin_states, n_k_points_task), intent(in) :: KS_eigenvector_complex
    real*8, dimension( n_hamiltonian_matrix_size ),                        intent(in) :: overlap_matrix
    integer,                                                               intent(in) :: n_states_in
    real*8, dimension(n_states_in, n_spin_states, n_k_points),             intent(in) :: KS_eigenvalue 
    real*8, dimension(n_states_in, n_spin_states, n_k_points),             intent(in) :: occ_numbers
    real*8,                                                                intent(in) :: chemical_potential
    real*8,                                                                intent(in) :: n_electrons
    character*40,                                                          intent(in) :: filename
    integer, optional,                                                     intent(in) :: mode
  !*  INPUTS
  !*    o n_rows                  -- Number of rows in eigenvector array (i.e. basis function index)
  !*                                 May differ from n_basis for different storage formats.
  !*    o n_cols                  -- Number of columns in eigenvector array (i.e. states index)
  !*                                 May differ from n_states for different storage formats.
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvector          -- Eigenvector, real case
  !*                                 In the SOC cases, this is a dummy array.
  !*    o KS_eigenvector_complex  -- Eigenvector, complex case
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o KS_eigenvalue           -- Eigenvalues array
  !*    o occ_numbers             -- Occupation numbers array
  !*    o chemical_potential      -- Fermi level/chemical potential of electrons
  !*    o n_electrons             -- Number of electrons (per computational cell) in the system
  !*    o overlap_matrix          -- Overlap matrix for non-packed and packed Hamiltonians in the LAPACK case
  !*                                 In the ScaLAPACK case, we use the overlap matrices directly
  !*                                 from scalapack_wrapper (allowing support for local indexing)
  !*    o filename                -- The filename for the Mulliken decomposition, if desired
  !*                                 (I'm pretty sure this is a legacy argument and should be removed)
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (called subroutines write to screen and disk)
  !*  AUTHOR
  !*    William Huhn (Duke University)
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  NOTES
  !*    o In the case of scalar-relativistic ScaLAPACK calculations, the construction of the Mulliken
  !*      decomposition ignores the provided arguments and pulls directly from scalapack_wrapper.  I
  !*      personally don't like this, and in the spin-orbit-coupled case, we use the provided eigenvectors
  !*      in the proper BLACS format.  Though, as I write this, I do realize that we also pull the overlap
  !*      matrix from scalapack_wrapper.  Oops.
  !*  HISTORY
  !*    November 2017 - Rewritten basd on an older output_mulliken subroutine
  !*                    to have a more clear, modular design
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
  
    ! local variables

    ! We will do a Mulliken charge decomposition by the following criteria:
    ! Number of electrons in each KS state, per atom, spin channel, and 
    ! angular momentum component
    ! mulliken_decomp will contain this decomposition. All derived quantities
    ! are then accessible by way of appropriate sums.
    real*8, dimension(:, :, :, :, :), allocatable :: mulliken_decomp
 
    integer       :: n_spin_proj
    character*120 :: info_str
    character*11  :: sr_suffix
  
    character(*), parameter :: func = 'output_mulliken'
  
    ! begin work

    ! In the case that we're using SOC in this calculation but are currently requesting
    ! the scalar-relativistic value, modify the output to distinguish SR and SOC values 
    write(info_str,'(2X,A)') ''
    call localorb_info(info_str)

    if (calculate_perturbative_soc) then
      if (mode .eq. MULLIKEN_SR) then
        write(info_str,'(2X,A)') &
             'Performing scalar-relativistic Mulliken charge analysis on all atoms.'
      else if (mode .eq. MULLIKEN_SOC) then
        write(info_str,'(2X,A)') &
             'Performing spin-orbit-coupled Mulliken charge analysis on all atoms.'
      else 
        call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
      end if
    else
      write(info_str,'(2X,A)') &
           'Performing Mulliken charge analysis on all atoms.'
    end if
    call localorb_info(info_str)
  
    if (calculate_perturbative_soc .and. mode .eq. MULLIKEN_SR) then
      sr_suffix = ".dat.no_soc"
    else
      sr_suffix = ".dat"
    end if

    if (mode .eq. MULLIKEN_SR) then
       n_spin_proj = n_spin_states
    else if (mode .eq. MULLIKEN_SOC) then
       n_spin_proj = 2
    else
       call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
    end if

    ! Allocate and construct the Mulliken decomposition
    call aims_allocate( mulliken_decomp, 0, l_wave_max, 1, n_atoms, 1, n_states_in, 1, n_spin_proj, &
         1, n_k_points_task, "mulliken_decomp" )
    if (.not.present(mode) .or. mode .eq. MULLIKEN_SR) then
      call construct_mulliken_decomp( n_rows, n_cols, n_spin_states, KS_eigenvector, KS_eigenvector_complex, &
                                      overlap_matrix, n_states_in, n_spin_proj, mulliken_decomp )
    else if (mode .eq. MULLIKEN_SOC) then
      if (.not.calculate_perturbative_soc) then
        call aims_stop('Spin-orbit-coupled Mulliken was requested, but SOC was not used in this calculation, &
                       &exiting', func)
      end if
      call construct_mulliken_decomp_soc( n_rows, n_cols, n_spin_states, KS_eigenvector, KS_eigenvector_complex, &
                                          overlap_matrix, n_states_in, n_spin_proj, mulliken_decomp )
    else
      call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
    end if
 
    ! When doing locally constrained DFT, calculate the Mulliken decomposition summed over angular channels
    ! This shouldn't be done in the middle of an output subroutine, and the "total" Mulliken decomposition 
    ! should be distributed over MPI tasks, but I (WPH) don't want to break the locally constrained DFT code
    if (use_constraint .and. mode .eq. MULLIKEN_SR) then
      call construct_mulliken_decomp_total( n_states_in, n_spin_proj, mulliken_decomp )
    end if
    
    call check_occs('output_mulliken', occ_numbers, .false.)
 
    if( out_mulliken )then
      if (calculate_perturbative_soc) then
        write(info_str,'(2X,A)') &
           "Full analysis will be written to separate file 'Mulliken.out.no_soc'."
      else
        write(info_str,'(2X,A)') &
           "Full analysis will be written to separate file 'Mulliken.out'."
      end if 
      call localorb_info(info_str)
    else
      if (calculate_perturbative_soc) then
        write(info_str,'(2X,A)') &
           "Full analysis (per state, per k-point, etc.) will NOT be written to separate file 'Mulliken.out.no_soc'."
      else
        write(info_str,'(2X,A)') &
           "Full analysis (per state, per k-point, etc.) will NOT be written to separate file 'Mulliken.out'."
      end if
      call localorb_info(info_str)
      write(info_str,'(2X,A)') &
         "This file can be requested by stating 'output mulliken' explicitly."
      call localorb_info(info_str)
    end if
 
    ! Always write the projected charge per species to screen
    call output_proj_charge( n_states_in, n_spin_states, occ_numbers, n_spin_proj, mulliken_decomp, mode )
  
    ! Now perform the various output to disk requested by user 
    if( out_mulliken )then
      call write_decomp_file( n_states_in, n_spin_states, KS_eigenvalue, occ_numbers, &
                              n_spin_proj, mulliken_decomp, filename, mode )
    end if 
    if(out_l_proj_dos)then
      call write_species_proj_dos( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                   n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
    end if
    if(out_l_proj_dos_tetrahedron)then
      call write_species_proj_dos_tetrahedron( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                   n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
    end if
    if(out_atom_dos) then
      call write_atom_proj_dos( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
    end if 
    if(out_atom_dos_tetrahedron) then
      call write_atom_proj_dos_tetrahedron( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
    end if 
  
    if (allocated (mulliken_decomp) ) call aims_deallocate(mulliken_decomp, "mulliken_decomp")
  
  end subroutine output_mulliken
  !******

  !****f* mulliken/construct_mulliken_decomp
  !*  NAME
  !*    construct_mulliken_decomp
  !*  SYNOPSIS
  subroutine construct_mulliken_decomp( n_rows, n_cols, n_spin_states, KS_eigenvector, KS_eigenvector_complex, &
                                        overlap_matrix, &
                                        n_states_in, n_spin_proj, mulliken_decomp )
  !*  PURPOSE
  !*    The subroutine performs the Mulliken decomposition for the provided scalar-relativistic eigenvectors
  !*  USES
     use dimensions,            only : n_states, n_basis, n_spin, n_hamiltonian_matrix_size, &
                                       n_atoms, l_wave_max, n_k_points, n_k_points_task, &
                                       n_periodic
     use localorb_io,           only : use_unit
     use mpi_tasks,             only : myid, n_tasks, aims_stop
     use basis,                 only : basis_l, basis_atom
     use runtime_choices,       only : use_scalapack, packed_matrix_format, & 
                                       PM_index, PM_none, real_eigenvectors
     use pbc_lists,             only : n_cells_in_hamiltonian, index_hamiltonian, &
                                       column_index_hamiltonian, k_phase, &
                                       cbasis_to_basis, cbasis_to_atom
     use synchronize_mpi_basic, only : sync_vector
     use scalapack_wrapper,     only : my_scalapack_comm_all, construct_mulliken_decomp_scalapack
     implicit none
  !*  ARGUMENTS
     integer,                                                               intent(in)  :: n_rows
     integer,                                                               intent(in)  :: n_cols
     integer,                                                               intent(in)  :: n_spin_states
     real*8, dimension(n_rows, n_cols, n_spin_states, n_k_points_task),     intent(in)  :: KS_eigenvector
     complex*16, dimension(n_rows, n_cols, n_spin_states, n_k_points_task), intent(in)  :: KS_eigenvector_complex
     real*8, dimension( n_hamiltonian_matrix_size ),                        intent(in)  :: overlap_matrix
     integer,                                                               intent(in)  :: n_states_in
     integer,                                                               intent(in)  :: n_spin_proj
     real*8, dimension( 0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task) , &
                                                                            intent(out) :: mulliken_decomp
  !*  INPUT
  !*    o n_rows                  -- Number of rows in eigenvector array (i.e. basis function index)
  !*                                 May differ from n_basis for different storage formats.
  !*    o n_cols                  -- Number of columns in eigenvector array (i.e. states index)
  !*                                 May differ from n_states for different storage formats.
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvector          -- Eigenvector, real case
  !*                                 In the SOC cases, this is a dummy array.
  !*    o KS_eigenvector_complex  -- Eigenvector, complex case
  !*    o overlap_matrix          -- Overlap matrix for non-packed and packed Hamiltonians in the LAPACK case
  !*                                 In the ScaLAPACK case, we use the overlap matrices directly
  !*                                 from scalapack_wrapper (allowing support for local indexing)
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*  OUTPUT
  !*    o mulliken_decomp         -- The scalar-relativistic Mulliken decomposition onto atoms, angular 
  !*                                 channels, and spin channels
  !*  AUTHOR
  !*    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  NOTES
  !*    o For ScaLAPACK calculations, this subroutine is a wrapper around construct_mulliken_decomp_scalapack
  !*      in scalapack_wrapper.  I don't like this behavior at all, but I don't have the time to change it.
  !*  HISTORY
  !*    November 2017 - Forked out of the output_mulliken and output_mulliken_soc to
  !*                    allow for reusability
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer     :: i_basis_1, i_basis_2, i_cell, i_state, i_size, i_index, i_spin
    integer     :: i_k, i_k_point
    complex*16  :: mul_temp

    character(*), parameter :: func = 'construct_mulliken_decomp'

    ! Make sure the input array dimensions correspond to the "normal" (scalar-relativistic)
    ! indices for arrays
    if (n_rows .ne. n_basis) then
      call aims_stop('Incorrect array dimensions for number of rows in eigenvectors when calculating&
                     & Mulliken decomposition, exiting.', func) 
    end if
    if (n_cols .ne. n_states) then
      call aims_stop('Incorrect array dimensions for number of columns in eigenvectors when calculating&
                     & Mulliken decomposition, exiting.', func) 
    end if
    if (n_spin_states .ne. n_spin) then
      call aims_stop('Incorrect array dimensions for number of spin channels in eigenvectors when calculating&
                     & Mulliken decomposition, exiting.', func) 
    end if
    if (n_spin_proj .ne. n_spin) then
      call aims_stop('Incorrect array dimensions for number of spin channels in Mulliken decomposition when&
                     & calculating Mulliken decomposition, exiting.', func) 
    end if
    if (n_states_in .ne. n_states) then
      call aims_stop('Incorrect array dimensions for number of states in Mulliken decomposition when&
                     & calculating  Mulliken decomposition, exiting.', func) 
    end if 
    ! From here on, we use the "normal" SR indices for KS_eigenvector and other arrays 

    mulliken_decomp = 0.d0

    ! Notice we do not sum up the density matrix because want a per-state 
    ! projection analysis first - Mulliken charges are summed up only thereafter.
    if(use_scalapack) then
      call construct_mulliken_decomp_scalapack(mulliken_decomp(0,1,1,1,1))
      call sync_vector( mulliken_decomp(0,1,1,1,1), (l_wave_max+1)*n_atoms*n_states*n_spin*1, &
           my_scalapack_comm_all )
    else
      ! Note the format of overlap_matrix is packed - hence we have to set up the 
      ! necessary matrix multiplications in a really awkward way, without use of
      ! BLAS, and with an if statement in a really bad place ...
      if(n_periodic == 0 .and. packed_matrix_format==PM_none)then
         do i_spin = 1, n_spin, 1
            do i_state = 1, n_states, 1
 
               i_index = 0
               do i_basis_2 = 1, n_basis, 1
                  do i_basis_1 = 1, i_basis_2, 1
                     i_index = i_index+1
 
                     ! 1st pass over all matrix elements
                     mulliken_decomp( basis_l(i_basis_1), basis_atom(i_basis_1), i_state, i_spin,1 ) = &
                          mulliken_decomp( basis_l(i_basis_1), basis_atom(i_basis_1), i_state, i_spin,1 ) + &
                          KS_eigenvector(i_basis_1, i_state, i_spin,1) * &
                          overlap_matrix(i_index) * &
                          KS_eigenvector(i_basis_2, i_state, i_spin,1)
 
                     ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                     if (i_basis_1.ne.i_basis_2) then
                        mulliken_decomp( basis_l(i_basis_2), basis_atom(i_basis_2), i_state, i_spin,1 ) = &
                             mulliken_decomp( basis_l(i_basis_2), basis_atom(i_basis_2), i_state, i_spin,1 ) + &
                             KS_eigenvector(i_basis_2, i_state, i_spin,1) * &
                             overlap_matrix(i_index) * &
                             KS_eigenvector(i_basis_1, i_state, i_spin,1) 
                     end if
 
                  enddo
               enddo
 
            enddo
         enddo
      else
    
         if(packed_matrix_format /= PM_index)then
            write(use_unit,*) 'Error: periodic Mulliken supports only packed matrix format index'
            return
         end if
    
    
         if(real_eigenvectors)then
    
            do i_spin = 1, n_spin, 1
               i_k = 0           
               do i_k_point = 1, n_k_points
                  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                     i_k = i_k + 1
    
                     do i_state = 1, n_states, 1
                        do i_cell = 1,n_cells_in_hamiltonian-1
    
                           do i_basis_2 = 1, n_basis
    
                              if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
    
                                 i_index = index_hamiltonian(1,i_cell, i_basis_2)-1
    
                                 do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
    
                                    i_index = i_index + 1
                                    i_basis_1 =  column_index_hamiltonian(i_index)
    
    
               ! 1st pass over all matrix elements
               mulliken_decomp( basis_l(Cbasis_to_basis(i_basis_1)), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k ) = &
                    mulliken_decomp( basis_l(Cbasis_to_basis(i_basis_1)), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k ) + &
                    KS_eigenvector(Cbasis_to_basis(i_basis_1), i_state, i_spin,i_k) * &
                    overlap_matrix(i_index) * &
                    KS_eigenvector(Cbasis_to_basis(i_basis_2), i_state, i_spin,i_k) * dble(k_phase(i_cell,i_k_point))
    
               ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
               if (i_basis_1.ne.i_basis_2) then
                  mulliken_decomp( basis_l(Cbasis_to_basis(i_basis_2)), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k ) = &
                     mulliken_decomp( basis_l(Cbasis_to_basis(i_basis_2)), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k ) + &
                     KS_eigenvector(Cbasis_to_basis(i_basis_2), i_state, i_spin,i_k) * &
                     overlap_matrix(i_index) * &
                     KS_eigenvector(Cbasis_to_basis(i_basis_1), i_state, i_spin,i_k) * dble(k_phase(i_cell,i_k_point))
    
               end if
                                 end do
                              end if
                           end do
                        end do
                     end do
                  end if
               end do
            end do
    
    
         else
    
    
            do i_spin = 1, n_spin, 1
               i_k = 0           
    
               do i_k_point = 1, n_k_points
                 ! write(use_unit,*) i_k_point
    
                  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                     i_k = i_k + 1
    
                     do i_state = 1, n_states, 1
                        do i_cell = 1,n_cells_in_hamiltonian-1
    
                           do i_basis_2 = 1, n_basis
    
                              if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
    
                                 i_index = index_hamiltonian(1,i_cell, i_basis_2)-1
    
                                 do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
    
                                    i_index = i_index + 1
                                    i_basis_1 =  column_index_hamiltonian(i_index)
    
    
                                    mul_temp =  KS_eigenvector_complex(i_basis_1, i_state, i_spin,i_k) * &
                                         dconjg(KS_eigenvector_complex(i_basis_2, i_state, i_spin,i_k)) &
                                         * dconjg(k_phase(i_cell,i_k_point)) &
                                         * overlap_matrix(i_index)
    
                                    mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k ) = &
                                    mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k ) + &
                                    dble(mul_temp)
    
                                    ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                                    if (i_basis_1.ne.i_basis_2) then
                                       mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k ) = &
                                       mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin,i_k ) + &
                                       dble(mul_temp)
    
                                    end if
                                 end do
                              end if
                           end do
                        end do
                     end do
    
    
    
                  end if
               end do
    
            end do
    
         end if
      end if
    
    end if ! use_scalapack
  end subroutine construct_mulliken_decomp   
  !******

  !****f* mulliken/construct_mulliken_decomp_total
  !*  NAME
  !*    construct_mulliken_decomp_total
  !*  SYNOPSIS
  subroutine construct_mulliken_decomp_total( n_states_in, n_spin_proj, mulliken_decomp )
  !*  PURPOSE
  !*    This subroutine sums the Mulliken decomposition over l-channels.  Used for
  !*    locally constrained DFT.
  !*  USES
    use dimensions,            only : n_states, n_spin, n_atoms, l_wave_max, n_k_points, &
                                      n_k_points_task
    use mpi_tasks,             only : myid, n_tasks, aims_stop
    use geometry,              only : species
    use species_data,          only : l_shell_max, species_pseudoized
    use synchronize_mpi_basic, only : sync_vector
    implicit none
  !*  ARGUMENTS
    integer, intent(in) :: n_states_in
    integer, intent(in) :: n_spin_proj
    real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
             intent(in) :: mulliken_decomp
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*  OUTPUT
  !*    none (modifies module variable)
  !*  AUTHOR
  !*    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
  
    ! local variables
    integer :: i_atom, i_spin, i_k, i_k_point, i_state

    character(*), parameter :: func = 'construct_mulliken_decomp'

    ! Make sure the input array dimensions correspond to the "normal" (scalar-relativistic)
    ! indices for arrays
    if (n_states_in .ne. n_states) then
      call aims_stop('Incorrect array dimensions for number of states when&
                     & calculating total Mulliken decomposition, exiting.', func) 
    end if 
    if (n_spin_proj .ne. n_spin) then
      call aims_stop('Incorrect array dimensions for number of spin channels when calculating&
                     & total Mulliken decomposition, exiting.', func) 
    end if
    ! From here on, we use the "normal" SR indices for KS_eigenvector and other arrays 

    if (allocated(mulliken_decomp_total)) deallocate( mulliken_decomp_total )
    allocate( mulliken_decomp_total(n_atoms, n_states, n_spin, n_k_points) )

    do i_atom = 1, n_atoms
       if(species_pseudoized(species(i_atom))) cycle
       do i_spin = 1, n_spin, 1
          i_k = 0
          do i_k_point = 1, n_k_points
             if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                i_k = i_k + 1
                do i_state=1,n_states,1
                   mulliken_decomp_total(i_atom, i_state, i_spin, i_k_point) = &
                        sum( mulliken_decomp(0:l_shell_max(species(i_atom)),i_atom,i_state,i_spin,i_k)) 
                end do
             end if
          end do
       end do
    end do
    call sync_vector( mulliken_decomp_total, n_atoms*n_states*n_spin*n_k_points )

  end subroutine construct_mulliken_decomp_total
  !******

  !****f* mulliken/construct_mulliken_decomp_soc
  !*  NAME
  !*    construct_mulliken_decomp_soc
  !*  SYNOPSIS
  subroutine construct_mulliken_decomp_soc( n_rows, n_cols, n_spin_states, eigenvec_soc_dummy, eigenvec_soc, &
                                        overlap_matrix, &
                                        n_states_in, n_spin_proj, mulliken_decomp )
  !*  PURPOSE
  !*    The subroutine performs the Mulliken decomposition for the provided spin-orbit-coupled eigenvectors
  !*  USES
     use dimensions,            only : n_basis, n_hamiltonian_matrix_size, n_atoms, l_wave_max, &
                                       n_k_points, n_k_points_task, n_periodic
     use dimensions_soc,        only : n_saved_states_soc, n_basis_soc, n_basis_soc_coll
     use localorb_io,           only : use_unit
     use mpi_tasks,             only : myid, n_tasks, aims_stop
     use basis,                 only : basis_l, basis_atom
     use runtime_choices,       only : use_scalapack, packed_matrix_format, & 
                                       PM_index, PM_none, real_eigenvectors
     use pbc_lists,             only : n_cells_in_hamiltonian, index_hamiltonian, &
                                       column_index_hamiltonian, k_phase, &
                                       cbasis_to_basis, cbasis_to_atom
     use synchronize_mpi_basic, only : sync_vector
     use scalapack_wrapper,     only : npcol, nprow, sc_desc, mxld, mxcol, my_scalapack_id, &
                                       my_scalapack_comm_all, ovlp_stored, ovlp_complex_stored
     use scalapack_soc,         only : sc_desc_soc_vec, mxld_soc_vec, mxcol_soc_vec, l_row_soc_vec, &
                                       l_col_soc_vec
     use aims_memory_tracking,  only : aims_allocate, aims_deallocate
     implicit none
  !*  ARGUMENTS
     integer,                                                               intent(in)  :: n_rows
     integer,                                                               intent(in)  :: n_cols
     integer,                                                               intent(in)  :: n_spin_states
     real*8, dimension(n_rows, n_cols, n_spin_states, n_k_points_task),     intent(in)  :: eigenvec_soc_dummy
     complex*16, dimension(n_rows, n_cols, n_spin_states, n_k_points_task), intent(in)  :: eigenvec_soc
     real*8, dimension( n_hamiltonian_matrix_size ),                        intent(in)  :: overlap_matrix
     integer,                                                               intent(in)  :: n_states_in
     integer,                                                               intent(in)  :: n_spin_proj
     real*8, dimension( 0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task) , &
                                                                            intent(out) :: mulliken_decomp
  !*  INPUT
  !*    o n_rows                  -- Number of rows in eigenvector array (i.e. basis function index)
  !*                                 May differ from n_basis for different storage formats.
  !*    o n_cols                  -- Number of columns in eigenvector array (i.e. states index)
  !*                                 May differ from n_states for different storage formats.
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o eigenvec_soc_dummy      -- Eigenvector, real case
  !*                                 In the SOC cases, this is a dummy array.
  !*    o eigenvec_soc            -- Eigenvector, complex case
  !*    o overlap_matrix          -- Overlap matrix for non-packed and packed Hamiltonians in the LAPACK case
  !*                                 In the ScaLAPACK case, we use the overlap matrices directly
  !*                                 from scalapack_wrapper (allowing support for local indexing)
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*  OUTPUT
  !*    o mulliken_decomp         -- The spin-orbit-coupled Mulliken decomposition onto atoms, angular channels, 
  !*                                 and spin channels
  !*  AUTHOR
  !*    William Huhn (Duke University)
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  NOTES
  !*    o In the case of scalar-relativistic ScaLAPACK calculations, the construction of the Mulliken
  !*      decomposition ignores the provided arguments and pulls directly from scalapack_wrapper.  I
  !*      personally don't like this, and in the spin-orbit-coupled case, we use the provided eigenvectors
  !*      in the proper BLACS format.  Though, as I write this, I do realize that we also pull the overlap
  !*      matrix from scalapack_wrapper.  Oops.
  !*  HISTORY
  !*    November 2017 - Forked out of the output_mulliken and output_mulliken_soc to
  !*                    allow for reusability
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer                                 :: i_basis_1, i_basis_2, i_cell, i_state, i_size, i_index, i_spin
    integer                                 :: i_k, i_k_point
    integer                                 :: basis_offset
    complex*16                              :: mul_temp
    complex*16, dimension(:,:), allocatable :: tmp_complex
    complex*16, dimension(:,:), allocatable :: temp_ovlp_complex_stored

    character(*), parameter :: func = 'construct_mulliken_decomp_soc'

    ! Make sure the input array dimensions correspond to the "normal" spin-orbit-coupled
    ! indices for arrays
    if (use_scalapack) then
      if (n_rows .ne. mxld_soc_vec) then
        call aims_stop('Incorrect array dimensions for number of rows in eigenvectors when calculating&
                       & SOC-perturbed Mulliken decomposition, exiting.', func) 
      end if
      if (n_cols .ne. mxcol_soc_vec) then
        call aims_stop('Incorrect array dimensions for number of columns in eigenvectors when calculating&
                       & SOC-perturbed Mulliken decomposition, exiting.', func) 
      end if
    else 
      if (n_rows .ne. n_basis_soc) then
        call aims_stop('Incorrect array dimensions for number of rows in eigenvectors when calculating&
                       & SOC-perturbed Mulliken decomposition, exiting.', func) 
      end if
      if (n_cols .ne. n_saved_states_soc) then
        call aims_stop('Incorrect array dimensions for number of columns in eigenvectors when calculating&
                       & SOC-perturbed Mulliken decomposition, exiting.', func) 
      end if
    end if
    if (n_spin_states .ne. 1) then
      call aims_stop('Incorrect array dimensions for number of spin channels in eigenvectors when calculating&
                     & SOC-perturbed Mulliken decomposition, exiting.', func) 
    end if
    if (n_spin_proj .ne. 2) then
      call aims_stop('Incorrect array dimensions for number of spin channels in Mulliken decomposition when&
                     & calculating SOC-perturbed Mulliken decomposition, exiting.', func) 
    end if
    if (n_states_in .ne. n_saved_states_soc) then
      call aims_stop('Incorrect array dimensions for number of states in Mulliken decomposition when&
                     & calculating SOC-perturbed Mulliken decomposition, exiting.', func) 
    end if
    ! From here on, we use the "normal" SOC indices for KS_eigenvector and other arrays 

    mulliken_decomp = 0.d0

    ! Notice we do not sum up the density matrix because want a per-state 
    ! projection analysis first - Mulliken charges are summed up only thereafter.

    if(use_scalapack) then
      ! The work in this routine must be done only on the working set
      if(my_scalapack_id>=npcol*nprow) return
  
      call aims_allocate( tmp_complex, mxld_soc_vec, mxcol_soc_vec, "tmp_complx" )
      tmp_complex = (0.0d0, 0.0d0)
  
      ! Multiply the spin-up and spin-down blocks of the eigenvectors independently by the overlap matrix
      ! (this is equivalent to multiplying by the "full" overlap matrix, where the off-diagonal blocks are 
      ! zero because the spinors for basis functions are orthonormal)
      if (real_eigenvectors) then
        call aims_allocate( temp_ovlp_complex_stored, mxld, mxcol, "tmp_ovlp_complex_stored" )
        temp_ovlp_complex_stored = dcmplx(ovlp_stored)
  
        call pzgemm('N','N',n_basis,n_saved_states_soc,n_basis,&
                    (1.d0,0.d0),&
                    temp_ovlp_complex_stored, 1,1,sc_desc, &
                    eigenvec_soc, 1,1,sc_desc_soc_vec,&
                    (0.d0,0.d0),&
                    tmp_complex, 1,1,sc_desc_soc_vec)
        call pzgemm('N','N',n_basis,n_saved_states_soc,n_basis,&
                    (1.d0,0.d0),&
                    temp_ovlp_complex_stored, 1,1,sc_desc, &
                    eigenvec_soc, 1+n_basis_soc_coll/2,1,sc_desc_soc_vec,&
                    (1.d0,0.d0),&
                    tmp_complex, 1+n_basis_soc_coll/2,1,sc_desc_soc_vec)
  
       call aims_deallocate( temp_ovlp_complex_stored, "temp_ovlp_complex_stored" )
      else
        call pzgemm('N','N',n_basis,n_saved_states_soc,n_basis,&
                    (1.d0,0.d0),&
                    ovlp_complex_stored, 1,1,sc_desc, &
                    eigenvec_soc, 1,1,sc_desc_soc_vec,&
                    (0.d0,0.d0),&
                    tmp_complex, 1,1,sc_desc_soc_vec)
        call pzgemm('N','N',n_basis,n_saved_states_soc,n_basis,&
                    (1.d0,0.d0),&
                    ovlp_complex_stored, 1,1,sc_desc, &
                    eigenvec_soc, 1+n_basis_soc_coll/2,1,sc_desc_soc_vec,&
                    (1.d0,0.d0),&
                    tmp_complex, 1+n_basis_soc_coll/2,1,sc_desc_soc_vec)
      end if
  
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
          end if
  
          mul_temp = dble(conjg(eigenvec_soc(l_row_soc_vec(i_basis_1),l_col_soc_vec(i_state),1,1)) * &
                                tmp_complex (l_row_soc_vec(i_basis_1),l_col_soc_vec(i_state)))
          mulliken_decomp(basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin, 1) = &
               mulliken_decomp(basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin, 1) + mul_temp
        enddo
      enddo
  
      call aims_deallocate( tmp_complex, "tmp_complex" )
      call sync_vector( mulliken_decomp(0,1,1,1,1), (l_wave_max+1)*n_atoms*n_saved_states_soc*2*1, &
           my_scalapack_comm_all )
    else
    ! Note the format of overlap_matrix is packed - hence we have to set up the 
    ! necessary matrix multiplications in a really awkward way, without use of
    ! BLAS, and with an if statement in a really bad place ...
      if(n_periodic == 0 .and. packed_matrix_format==PM_none)then
         do i_spin = 1, 2, 1
           if (i_spin .eq. 1) then
             ! Spin-up components are requested
             basis_offset = 0
           else
             ! Spin-dn components are requested
             basis_offset = n_basis_soc_coll/2
           end if
 
           do i_state = 1, n_saved_states_soc, 1
             i_index = 0
             do i_basis_2 = 1, n_basis_soc_coll/2, 1
               do i_basis_1 = 1, i_basis_2, 1
                 i_index = i_index+1
                 ! 1st pass over all matrix elements
                 mulliken_decomp( basis_l(i_basis_1), basis_atom(i_basis_1), i_state, i_spin,1 ) = &
                      mulliken_decomp( basis_l(i_basis_1), basis_atom(i_basis_1), i_state, i_spin,1 ) + &
                      real(eigenvec_soc(basis_offset+i_basis_1, i_state, 1, 1) * &
                      overlap_matrix(i_index) * &
                      dconjg(eigenvec_soc(basis_offset+i_basis_2, i_state, 1, 1)))
                 ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                 if (i_basis_1.ne.i_basis_2) then
                   mulliken_decomp( basis_l(i_basis_2), basis_atom(i_basis_2), i_state, i_spin,1 ) = &
                        mulliken_decomp( basis_l(i_basis_2), basis_atom(i_basis_2), i_state, i_spin,1 ) + &
                        real(eigenvec_soc(basis_offset+i_basis_2, i_state, 1, 1) * &
                        overlap_matrix(i_index) * &
                        dconjg(eigenvec_soc(basis_offset+i_basis_1, i_state, 1, 1)))
                 end if
               enddo
             enddo
           enddo
         enddo
      else
         if(packed_matrix_format /= PM_index)then
           write(use_unit,*) 'Error: periodic Mulliken supports only packed matrix format index'
           return
         end if
    
         do i_spin = 1, 2, 1
           if (i_spin .eq. 1) then
             ! Spin-up components of basis set
             basis_offset = 0
           else
             ! Spin-dn components of basis set
             basis_offset = n_basis_soc_coll/2
           end if
    
           i_k = 0           
           do i_k_point = 1, n_k_points
             ! write(use_unit,*) i_k_point
             if  ( myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
               i_k = i_k + 1
               do i_state = 1, n_saved_states_soc, 1
                 do i_cell = 1,n_cells_in_hamiltonian-1
                   do i_basis_2 = 1, n_basis_soc_coll/2
                     if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then
                       i_index = index_hamiltonian(1,i_cell, i_basis_2)-1
                       do i_size = index_hamiltonian(1,i_cell, i_basis_2),index_hamiltonian(2,i_cell, i_basis_2)
                         i_index = i_index + 1
                         i_basis_1 =  column_index_hamiltonian(i_index)
                         mul_temp =  eigenvec_soc(basis_offset+i_basis_1, i_state,1,i_k) * &
                              dconjg(eigenvec_soc(basis_offset+i_basis_2, i_state,1,i_k)) &
                              * dconjg(k_phase(i_cell,i_k_point)) &
                              * overlap_matrix(i_index)
                         mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k ) = &
                              mulliken_decomp( basis_l(i_basis_1), Cbasis_to_atom(i_basis_1), i_state, i_spin, i_k ) + &
                              dble(mul_temp)
                         ! 2nd pass: must average all off-diagonal matrix elements (but not the diagonal)
                         if (i_basis_1.ne.i_basis_2) then
                           mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin, i_k ) = &
                                mulliken_decomp( basis_l(i_basis_2), Cbasis_to_atom(i_basis_2), i_state, i_spin, i_k ) + &
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
      end if ! n_periodic == 0
    end if ! use_scalapack
  end subroutine construct_mulliken_decomp_soc 
  !******

  !****f* mulliken/output_proj_charge
  !*  NAME
  !*    output_proj_charge
  !*  SYNOPSIS
  subroutine output_proj_charge( n_states_in, n_spin_states, occ_numbers, n_spin_proj, mulliken_decomp, mode)
  !*  PURPOSE
  !*    This subroutine calculates the electron occupation predicted by the decomposition to be associated with 
  !*    atoms and angular channels, as well as the difference from a neutral ion, and outputs the result to 
  !*    screen.
  !*  USES
    use dimensions,            only : n_atoms, l_wave_max, n_k_points, n_k_points_task
    use localorb_io,           only : use_unit, localorb_info
    use mpi_tasks,             only : myid, n_tasks, aims_stop
    use geometry,              only : species
    use species_data,          only : l_shell_max, species_name, species_pseudoized, species_z
    use runtime_choices,       only : spin_treatment, out_aims_json_log
    use pbc_lists,             only : k_weights
    use synchronize_mpi_basic, only : sync_vector
    use json_output,           only : write_mulliken_to_json_log
    implicit none
  !*  ARGUMENTS
    integer,                                                          intent(in) :: n_states_in
    integer,                                                          intent(in) :: n_spin_states
    real*8, dimension(n_states_in, n_spin_states, n_k_points),        intent(in) :: occ_numbers
    integer,                                                          intent(in) :: n_spin_proj
    real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
                                                                      intent(in) :: mulliken_decomp
    integer,                                                          intent(in) :: mode
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o occ_numbers             -- Occupation numbers array
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (writes to screen)
  !*  AUTHOR
  !*    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    real*8, dimension( 0:l_wave_max, n_atoms, n_spin_proj ) :: l_projected_charge
    real*8, dimension( n_atoms, n_spin_proj )               :: at_projected_charge
    real*8, dimension( n_atoms )                             :: spin_per_atom
    real*8, dimension( n_spin_proj )                        :: total_charge
    real*8                                                   :: charge_difference
 
    integer, dimension( n_spin_states, n_k_points )          :: max_occ_number

    ! darn, need an extra aux array of characters
    character*3,dimension(0:l_wave_max) :: l_channel
  
    integer       :: i_atom, i_l, i_k, i_k_point, i_spin_state, i_spin_proj, i_state
    character*120 :: info_str
 
    character(*), parameter :: func = 'output_proj_charge'
 
    l_projected_charge  = 0.d0
    at_projected_charge = 0.d0
    spin_per_atom       = 0.d0
  
    total_charge        = 0.d0

    write(info_str,'(2X,A)') "Summary of the per-atom charge analysis:"
    call localorb_info(info_str)
 
    do i_spin_proj = 1, n_spin_proj, 1

       if ( mode .eq. MULLIKEN_SR ) then 
          i_spin_state = i_spin_proj
       else if ( mode .eq. MULLIKEN_SOC ) then
          i_spin_state = 1
       else
          call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
       end if

       i_k = 0           
       do i_k_point = 1, n_k_points
  
          if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
             i_k = i_k + 1
             i_state = n_states_in
             do while ( (i_state.gt.1) .and. (occ_numbers(i_state,i_spin_state,i_k_point).eq.0.d0) )
                i_state = i_state - 1
             enddo
  
             if (i_state.eq.1) then
               ! workaround to avoid accessing a non-existing element occ_numbers(i_state=0,i_spin,i_k_point)
               ! which some compiler flags do not like
               if (occ_numbers(1,i_spin_state,i_k_point).eq.0.d0) then
                  max_occ_number(i_spin_state,i_k_point) = 0
               else
                  max_occ_number(i_spin_state,i_k_point) = 1
               end if
             else
               max_occ_number(i_spin_state,i_k_point) = i_state
             end if 
  
             do i_atom = 1, n_atoms, 1
                if(species_pseudoized(species(i_atom))) cycle
                ! sum up angular-momentum resolved charge projected on atoms in each spin channel
                ! sum up total charge projected on atoms in each spin channel
                do i_l = 0, l_shell_max(species(i_atom))
  
                   do i_state = 1, max_occ_number(i_spin_state,i_k_point), 1
                      l_projected_charge( i_l, i_atom, i_spin_proj) = &
                           l_projected_charge( i_l, i_atom, i_spin_proj) + &
                           occ_numbers(i_state,i_spin_state,i_k_point) & 
                           * mulliken_decomp(i_l, i_atom, i_state, i_spin_proj,i_k)*k_weights(i_k_point)
  
                   enddo
  
                enddo
  
             enddo
  
          end if
       end do
    enddo
    call sync_vector( l_projected_charge(0,1,1),(1+l_wave_max)*n_atoms*n_spin_proj )
   
    do i_spin_proj = 1, n_spin_proj, 1
       do i_atom = 1, n_atoms, 1
         if(species_pseudoized(species(i_atom))) cycle
          at_projected_charge(i_atom, i_spin_proj) =  sum(l_projected_charge( :, i_atom, i_spin_proj ))
       end do
       total_charge(i_spin_proj) = sum(at_projected_charge(:, i_spin_proj))
    end do
  
    if (myid.eq.0) then
  
       charge_difference = 0.d0
       do i_spin_proj = 1, n_spin_proj, 1
          charge_difference = charge_difference + total_charge(i_spin_proj)
       enddo
       do i_atom = 1, n_atoms, 1
          if(species_pseudoized(species(i_atom))) cycle
          if (spin_treatment .eq. 1) then
             spin_per_atom(i_atom) = at_projected_charge(i_atom, 1) - at_projected_charge(i_atom, 2)
          end if
          charge_difference = charge_difference - species_z(species(i_atom))
       enddo
  
       ! At this point, write detailed Mulliken analysis output
  
       ! prepare angular momentum output strings
       do i_l = 0, l_wave_max, 1
          write(l_channel(i_l),'(A,I1)') "l=",i_l
       enddo
  
       ! First, write summary output into standard output file
       write(use_unit,'(2X,A)') "|"
       write(use_unit,'(2X,A,A,7X,A,10X,A,10(13X,A))') &
            "| "," atom","electrons","charge",( l_channel(i_l), i_l=0,l_wave_max,1 )
  
       do i_atom = 1, n_atoms, 1
          if(species_pseudoized(species(i_atom))) cycle
          write(use_unit,'(2X,A,I5,2X,F14.6,2X,F14.6,10(2X,F14.6))') &
               "| ", i_atom, sum(at_projected_charge(i_atom,1:n_spin_proj)), &
               -( sum(at_projected_charge(i_atom,1:n_spin_proj))-species_z(species(i_atom)) ), &
               (sum(l_projected_charge(i_l,i_atom,1:n_spin_proj)), i_l=0,l_shell_max(species(i_atom)),1 )
  
       enddo
  
       write(use_unit,'(2X,A)') "|"
  
       write(use_unit,'(2X,A,2X,F14.6,2X,F14.6)') "| Total", sum(total_charge(1:n_spin_proj)), -charge_difference
  
       write(use_unit,*)
  
       if (spin_treatment .eq. 1) then
  
          write(info_str,'(2X,A)') &
               "Summary of the per-atom spin analysis:"
          call localorb_info(info_str)
  
          write(use_unit,'(2X,A)') "|"
          write(use_unit,'(2X,A,A,3X,A,6X,A,10(9X,A))') &
               "| "," atom","spin",( l_channel(i_l), i_l=0,l_wave_max,1 )
  
          do i_atom = 1, n_atoms, 1
            if(species_pseudoized(species(i_atom))) cycle
             write(use_unit,'(2X,A,I5,2X,F10.6,2X,F10.6,10(2X,F10.6))') &
                  "| ", i_atom, spin_per_atom(i_atom), &
                  (l_projected_charge(i_l,i_atom,1) - l_projected_charge(i_l,i_atom,2), i_l=0,l_shell_max(species(i_atom)),1 )
  
          enddo
  
          write(use_unit,'(2X,A)') "|"
  
          write(use_unit,'(2X,A,2X,F10.6)') "| Total", sum(spin_per_atom(1:n_atoms))
  
          write(use_unit,*)
       end if
  
       if (out_aims_json_log) then
         call write_mulliken_to_json_log(n_spin_proj, at_projected_charge, &
              l_projected_charge, spin_per_atom, total_charge, charge_difference)
       end if
    end if ! end exclusion of all threads other than #0
  end subroutine output_proj_charge
  !******

  !****f* mulliken/write_decomp_file
  !*  NAME
  !*    write_decomp_file
  !*  SYNOPSIS
  subroutine write_decomp_file( n_states_in, n_spin_states, KS_eigenvalue, occ_numbers, &
                                n_spin_proj, mulliken_decomp, filename, mode )
  !*  PURPOSE
  !*    This subroutine writes out the decomposition at every k-point to disk.
  !*  USES
      use dimensions,            only : n_atoms, l_wave_max, n_k_points, n_k_points_task
      use constants,             only : hartree
      use localorb_io,           only : localorb_info
      use mpi_tasks,             only : myid, n_tasks, aims_stop, mpi_comm_global, &
                                        MPI_STATUS_SIZE, MPI_DOUBLE_PRECISION
      use geometry,              only : species
      use species_data,          only : l_shell_max, species_pseudoized
      use pbc_lists,             only : k_point_list, k_weights
      use synchronize_mpi_basic, only : sync_vector
      use generate_aims_uuid,    only : write_aims_uuid
  !*  ARGUMENTS
      implicit none

      integer,                                                          intent(in) :: n_states_in
      integer,                                                          intent(in) :: n_spin_states
      real*8, dimension(n_states_in, n_spin_states, n_k_points),        intent(in) :: KS_eigenvalue
      real*8, dimension(n_states_in, n_spin_states, n_k_points),        intent(in) :: occ_numbers
      integer,                                                          intent(in) :: n_spin_proj
      real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
                                                                        intent(in) :: mulliken_decomp
      character*40,                                                     intent(in) :: filename
      integer,                                                          intent(in) :: mode
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvalue           -- Eigenvalues array
  !*    o occ_numbers             -- Occupation numbers array
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*    o filename                -- The filename for the Mulliken decomposition
  !*                                 (I'm pretty sure this is a legacy argument and should be removed)
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (writes to disk)
  !*  AUTHOR
  !*    William Huhn (Duke University), based on an older non-parallel version
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

      integer :: i_atom, i_k, i_k_point, i_l, i_spin_state, i_spin_proj, i_state

      real*8, dimension(:, :),          allocatable :: mulliken_decomp_buffer

      character*120     :: info_str
      character(LEN=80) :: uuid_str

      ! darn, need an extra aux array of characters
      character*3,dimension(0:l_wave_max) :: l_channel

      ! Variables to handle parallel Mulliken output
      integer :: sending_task ! In all cases, myid.eq.0 receives
      integer :: mpi_status(MPI_STATUS_SIZE)
      integer :: mpi_err
    
      character(*), parameter :: func = 'write_decomp_file'

      ! In the periodic case, Mulliken.out is only written if requested.
      ! It can be quite large.
       
      ! Now write spin-/state-resolved projection to file
  
      ! Since the k-points are distributed across MPI tasks, but we want the output to occur in
      ! a well-formatted pattern which, annoyingly, loops over atoms and then k-points (ugh), we work 
      ! using a master-slave dynamic where task #0 requests data from other tasks when needed.  
      ! Accordingly, we break the logic into two cases:  all MPI tasks execept #0, and MPI task #0
  
      ! This method has n_atoms*n_k_points many communications, which is a bit much for my taste,
      ! but they're point-to-point communications.  The other ways I can think of to do this is 
      ! parallel I/O or communicating all k-point information to myid #0 for each i_atom cycle.
      ! Neither of these are pleasant alternatives.  The latter in particular re-introduces the
      ! memory issue that we're trying solve by having all k-point information on one rank (though 
      ! only on task 0, to be fair) and has the same number of communications, only masking them
      ! via (even worse) all-to-one calls.  

      write(info_str,'(2X,A)') "Writing Mulliken decomposition to disk ..."
      call localorb_info(info_str)
      write(info_str,'(2X,A)') ""
      call localorb_info(info_str)
  
      allocate( mulliken_decomp_buffer( 0:l_wave_max, n_states_in ) )
 
      ! prepare angular momentum output strings
      do i_l = 0, l_wave_max, 1
         write(l_channel(i_l),'(A,I1)') "l=",i_l
      enddo
  
      ! For all MPI tasks that are not MPI task #0, cycle through all k-points.  When we reach a 
      ! k-point for which we're responsible, patiently wait for MPI task #0 for request information from
      ! us, and send it when it does.
      if (myid.ne.0) then
         do i_atom = 1, n_atoms, 1
            if (species_pseudoized(species(i_atom))) cycle
            do i_spin_proj = 1, n_spin_proj, 1
               i_k = 0
               do i_k_point = 1, n_k_points
                  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                     i_k = i_k + 1
                     mulliken_decomp_buffer = mulliken_decomp(:,i_atom,:,i_spin_proj,i_k)
                     call MPI_send( mulliken_decomp_buffer(0,1), (l_wave_max+1)*n_states_in, MPI_DOUBLE_PRECISION, &
                                    0, i_k_point, &
                                    mpi_comm_global, mpi_status, mpi_err)
                  end if
               end do
            end do
         end do
      ! For MPI task #0, do all writing.  When we reach a k-point that we do not own, request it 
      ! from the MPI tasks which owns it
      else ! myid.eq.0
         open(50, file=filename)
  
         call write_aims_uuid(uuid_str)
         write(50,'(A,2X,A)') '#', uuid_str
  
         do i_atom = 1, n_atoms, 1
            if(species_pseudoized(species(i_atom))) cycle
            write(50,*)
            write(50,'(A,I5,A)') "Atom number ",i_atom, ":"
            write(50,*)
  
            do i_spin_proj = 1, n_spin_proj, 1
               if ( mode .eq. MULLIKEN_SR ) then 
                  i_spin_state = i_spin_proj
               else if ( mode .eq. MULLIKEN_SOC ) then
                  i_spin_state = 1
               else
                  call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
               end if

               if (n_spin_proj.gt.1) then
                  write(50,*)
                  if (i_spin_proj.eq.1) then
                     write(50,'(2X,A,A)') "Spin channel: ", "up"
                  else 
                     write(50,'(2X,A,A)') "Spin channel: ", "down"
                  end if
                  write(50,*)
               end if
  
               write(50,'(4X,A,7X,A,2X,A,7X,A,10(9X,A3))') & 
                    "State", "eigenvalue", "occ.number", "total", (l_channel(i_l),i_l=0,l_shell_max(species(i_atom)))
  
               i_k = 0
               do i_k_point = 1, n_k_points
                  ! Here, we fill the Mulliken decomposition buffer pertaining to the current k-point
   
                 ! We've reached a k-point that we own, populate the output buffer ourselves
                  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
                     i_k = i_k + 1
                     mulliken_decomp_buffer = mulliken_decomp(:,i_atom,:,i_spin_proj,i_k)
                  ! We've reached a k-point that we do not own, request the buffer from the MPI
                  ! task that owns it
                  else
                     sending_task = MOD(i_k_point, n_tasks) ! The task that has the current k-point
                     mulliken_decomp_buffer = 0.0d0
                     call MPI_recv( mulliken_decomp_buffer(0,1), (l_wave_max+1)*n_states_in, MPI_DOUBLE_PRECISION, &
                                    sending_task, i_k_point, &
                                    mpi_comm_global, mpi_status, mpi_err )
                  end if 
  
                  ! Now we've received the needed Mulliken decomposition at this k-point, write to disk
                  write(50,*)
                  write(50,'(A,I5,A,F12.8,1X,F12.8,1X,F12.8,A,F12.8)') & 
                       "k point number: ",i_k_point, ": ( ", k_point_list(i_k_point,1), & 
                       k_point_list(i_k_point,2), k_point_list(i_k_point,3), " ); weight: ", k_weights(i_k_point)
                  write(50,*)
  
                  do i_state=1,n_states_in,1
                     write(50,'(2X,I7,2X,F15.5,2X,F10.7,2X,F10.5,10(2X,F10.5))') &
                          i_state, KS_eigenvalue(i_state,i_spin_state,i_k_point)* hartree, &
                          occ_numbers(i_state,i_spin_state,i_k_point), &
                          sum( mulliken_decomp_buffer(0:l_shell_max(species(i_atom)),i_state) ), &
                          ( mulliken_decomp_buffer(i_l,i_state), i_l=0,l_shell_max(species(i_atom)) )
                  end do
               end do
            end do
         end do
  
         close(50)
      end if ! myid.ne.0
  end subroutine write_decomp_file
  !******
 
  !****f* mulliken/write_species_proj_dos
  !*  NAME
  !*    write_species_proj_dos
  !*  SYNOPSIS
  subroutine write_species_proj_dos( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                          n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
  !*  PURPOSE
  !*    This subroutine calculates the projected DOS predicted by the decomposition to be associated with 
  !*    angular channels, summed up over all atoms of a given species, and outputs the results to disk.
  !*  USES
     use dimensions,            only : n_atoms, l_wave_max, n_k_points, n_k_points_task, n_species, &
                                       spin_degeneracy
     use constants,             only : hartree, one_over_sqrt2
     use localorb_io,           only : use_unit, localorb_info
     use mpi_tasks,             only : myid, n_tasks, aims_stop, check_allocation
     use geometry,              only : species
     use species_data,          only : l_shell_max, species_name, species_pseudoized
     use runtime_choices,       only : l_proj_dos_high_energy, l_proj_dos_low_energy,  &
                                       l_proj_dos_n_en_points, l_proj_dos_alpha
     use pbc_lists,             only : k_weights
     use arch_specific,         only : arch_erf
     use synchronize_mpi_basic, only : sync_vector
     use generate_aims_uuid,    only : write_aims_uuid
  !*  ARGUMENTS
     implicit none

     integer,                                                   intent(in) :: n_states_in
     integer,                                                   intent(in) :: n_spin_states
     real*8, dimension(n_states_in, n_spin_states, n_k_points), intent(in) :: KS_eigenvalue
     real*8,                                                    intent(in) :: chemical_potential
     integer,                                                   intent(in) :: n_spin_proj
     real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
                                                                intent(in) :: mulliken_decomp
     character*11,                                              intent(in) :: sr_suffix
     integer,                                                   intent(in) :: mode
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvalue           -- Eigenvalues array
  !*    o chemical_potential      -- Fermi level/chemical potential of electrons
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*    o sr_suffix               -- Suffix to put on output file names for scalar-relativistic calculations
  !*                                 This is a legacy argument that should absorbed into the body of this
  !*                                 subroutine
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (writes to disk)
  !*  AUTHOR
  !*    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

     integer :: i_atom, i_e, i_k, i_k_point, i_l, i_spin_state, i_spin_proj, i_state, i_species 
     real*8  :: de, en, E1, E2

     character*120     :: info_str, outputformat
     character*50      :: proj_dos_filename
     character*120     :: outputformat_header
     character(LEN=80) :: uuid_str

     real*8, dimension(:,:,:,:), allocatable :: KS_dos
     real*8, dimension(:,:,:),   allocatable :: KS_sum_dos
     real*8, dimension(:),       allocatable :: proj_dos_erfs
    
     character(*), parameter :: func = 'write_species_proj_dos'

     !   VB: In my view, this call to check_norm is incorrect - at the very least,
     !   it will now no longer work when two chemical potentials (spin constraint) are employed.
     !   call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_spin)
  
     if(.not. allocated (KS_dos)) then
        allocate(KS_dos(0:l_wave_max,n_spin_states,l_proj_dos_n_en_points, n_species),stat=i_state) 
        call check_allocation(i_state, 'KS_dos                        ')
     endif
     if(.not. allocated (KS_sum_dos)) then
        allocate(KS_sum_dos(n_spin_states,l_proj_dos_n_en_points, n_species),stat=i_state) 
        call check_allocation(i_state, 'KS_sum_dos                    ')
     endif
     if (.not.allocated(proj_dos_erfs)) then
        allocate(proj_dos_erfs(l_proj_dos_n_en_points),stat=i_state)
        call check_allocation(i_state, 'proj_dos_erfs                 ')
     end if
     
     write(info_str,'(2X,A)') 'Calculating angular momentum projected density of states ...'
     call localorb_info(info_str)
     write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
     call localorb_info(info_str)
     
     de= (l_proj_dos_high_energy - l_proj_dos_low_energy)/dble(l_proj_dos_n_en_points-1)
     KS_dos =0.d0 
     KS_sum_dos =0.d0 
     
     ! compute DOS for angular momentum projection. Make sure that the mulliken decomposition is properly normalized!!!
     i_k = 0
     do i_k_point = 1, n_k_points
        if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
           i_k = i_k + 1
           do i_spin_state = 1, n_spin_states
              do i_state = 1, n_states_in
                 ! calculate state- and k-dependent prefactors
                 do i_e = 1, l_proj_dos_n_en_points
                    en= l_proj_dos_low_energy + dble(i_e-1)*de   
                    E1 = en - de/2d0
                    E2 = en + de/2d0
                    proj_dos_erfs(i_e) = &
                         arch_erf((E2-(KS_eigenvalue(i_state,i_spin_state,i_k_point))*hartree)*one_over_sqrt2/l_proj_dos_alpha)  &
                         -arch_erf((E1-(KS_eigenvalue(i_state,i_spin_state,i_k_point))*hartree)*one_over_sqrt2/l_proj_dos_alpha)
                 end do
                 proj_dos_erfs(:) = proj_dos_erfs(:)*k_weights(i_k_point)/(2d0*dE)
                 ! calculate dos projection 
                 do i_e = 1, l_proj_dos_n_en_points
                    do i_atom = 1, n_atoms, 1
                       if(species_pseudoized(species(i_atom))) cycle
                       if (mode.eq.MULLIKEN_SR) then
                          i_spin_proj = i_spin_state
                          do i_l = 0, l_shell_max(species(i_atom))
                             KS_dos(i_l, i_spin_state,i_e,species(i_atom)) =  KS_dos (i_l,i_spin_state,i_e,species(i_atom)) + &
                                  spin_degeneracy * mulliken_decomp(i_l,i_atom,i_state,i_spin_proj,i_k) * proj_dos_erfs(i_e)
                          end do
                       else if (mode .eq. MULLIKEN_SOC .and. i_spin_state .eq. 1) then
                          ! The second condition is paranoia; n_spin_states should equal 1 for SOC calculations
                          do i_l = 0, l_shell_max(species(i_atom))
                             KS_dos(i_l, i_spin_state,i_e,species(i_atom)) = KS_dos (i_l,i_spin_state,i_e,species(i_atom)) + &
                                  (mulliken_decomp(i_l,i_atom,i_state,1,i_k) + mulliken_decomp(i_l,i_atom,i_state,2,i_k)) &
                                  * proj_dos_erfs(i_e)
                          end do
                       else
                          call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
                       end if
                    end do
                 end do
              end do
           end do
        end if
     end do
     call sync_vector( KS_dos(0,1,1,1),(1+l_wave_max)*n_spin_states*l_proj_dos_n_en_points*n_species)
     call sync_vector( KS_sum_dos(1,1,1),n_spin_states*l_proj_dos_n_en_points*n_species)
  
     do i_spin_state = 1, n_spin_states, 1
        do i_e = 1, l_proj_dos_n_en_points, 1
           do i_species = 1, n_species
              do i_l = 0, l_shell_max(i_species)
                 KS_sum_dos(i_spin_state,i_e,i_species) =  KS_sum_dos(i_spin_state,i_e,i_species) + &
                      KS_dos(i_l, i_spin_state,i_e,i_species)
              end do
           end do
        end do
     end do
  
     ! output on thread 0 only 
     if (myid.eq.0) then
        if (n_spin_states.eq.1) then
           do i_species = 1, n_species
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos', trim(sr_suffix)
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(i_species)+3
              write(use_unit,'(2X,5A)') '| writing projected DOS (shifted by the chemical potential) for species ',& 
                   trim(species_name(i_species)), &
                   ' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
              write(88,'(3A)') '# Angular momentum resolved density for species ',& 
                trim(species_name(i_species)), 'as calculated by FHI-aims'
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                chemical_potential*hartree, ' eV'
              
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_raw', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing projected DOS (raw data) for species ',& 
                trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
              write(89,'(3A)') '# Angular momentum resolved density for species ',& 
                trim(species_name(i_species)), 'as calculated by FHI-aims'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level '
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(i_species)+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
              
              do i_e = 1, l_proj_dos_n_en_points ,1
                 en = l_proj_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_species), &
                      (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(89,outputformat) & 
                      en, KS_sum_dos(1,i_e,i_species), (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
              enddo
              close(88)
              close(89)
           end do
        else 
           do i_species = 1, n_species
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_up', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_down', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)           
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
              write(88,'(3A)') '# spin-up angular momentum projected density of states for species ',& 
                   species_name(i_species), ' as calculated by FHI-aims '
              write(89,'(3A)') '# spin-down angular momentum projected density of states for species ',& 
                   species_name(i_species), ' as calculated by FHI-aims '
              
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_up_raw', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (raw data) for species ',trim(species_name(i_species)),&
                   ' to file ',trim(proj_dos_filename),'.'
              open(90, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(90,'(A,2X,A)') '#', uuid_str
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_down_raw', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (raw data) for species ',trim(species_name(i_species)),&
                   ' to file ',trim(proj_dos_filename),'.'
              open(91, file=proj_dos_filename)                      
              call write_aims_uuid(uuid_str)
              write(91,'(A,2X,A)') '#', uuid_str
  
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(i_species)+3
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ',& 
                   chemical_potential*hartree,' eV'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ',&
                   chemical_potential*hartree,' eV'
              write(90,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              write(91,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(i_species)+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
              write(90,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
              write(91,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
              
              do i_e = 1, l_proj_dos_n_en_points ,1
                 en = l_proj_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_species), &
                      (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(89,outputformat) en-chemical_potential*hartree, KS_sum_dos(2,i_e,i_species), &
                      (KS_dos(i_l, 2,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(90,outputformat) en, KS_sum_dos(1,i_e,i_species), (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(91,outputformat) en, KS_sum_dos(2,i_e,i_species), (KS_dos(i_l, 2,i_e,i_species), i_l=0,l_shell_max(i_species))
              enddo
              close(89)
              close(88)  
              close(90)
              close(91)
           end do
        end if
        
     end if
     call localorb_info(" ")
     if (allocated(KS_dos)       ) deallocate(KS_dos)        
     if (allocated(KS_sum_dos)   ) deallocate(KS_sum_dos)    
     if (allocated(proj_dos_erfs)) deallocate(proj_dos_erfs)
  end subroutine write_species_proj_dos
  !******

  !****f* mulliken/write_atom_proj_dos
  !*  NAME
  !*    write_atom_proj_dos
  !*  SYNOPSIS
  subroutine write_atom_proj_dos( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                          n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
  !*  PURPOSE
  !*    This subroutine calculates the projected DOS predicted by the decomposition to be associated with 
  !*    angular channels for each individual atom and outputs the results to disk.
  !*  USES
     use dimensions,            only : n_atoms, l_wave_max, n_k_points, n_k_points_task, n_species, &
                                       spin_degeneracy
     use constants,             only : hartree, one_over_sqrt2
     use localorb_io,           only : use_unit, localorb_info
     use mpi_tasks,             only : myid, n_tasks, aims_stop, check_allocation
     use geometry,              only : species
     use species_data,          only : l_shell_max, species_name, species_pseudoized
     use runtime_choices,       only : atom_dos_high_energy, atom_dos_low_energy,  &
                                       atom_dos_n_en_points, atom_dos_alpha
     use pbc_lists,             only : k_weights
     use arch_specific,         only : arch_erf
     use synchronize_mpi_basic, only : sync_vector
     use generate_aims_uuid,    only : write_aims_uuid
  !*  ARGUMENTS
     implicit none

     integer,                                                          intent(in) :: n_states_in
     integer,                                                          intent(in) :: n_spin_states
     real*8, dimension(n_states_in, n_spin_states, n_k_points),        intent(in) :: KS_eigenvalue
     real*8,                                                           intent(in) :: chemical_potential
     integer,                                                          intent(in) :: n_spin_proj
     real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
                                                                       intent(in) :: mulliken_decomp
     character*11,                                                     intent(in) :: sr_suffix
     integer,                                                          intent(in) :: mode
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvalue           -- Eigenvalues array
  !*    o chemical_potential      -- Fermi level/chemical potential of electrons
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*    o sr_suffix               -- Suffix to put on output file names for scalar-relativistic calculations
  !*                                 This is a legacy argument that should absorbed into the body of this
  !*                                 subroutine
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (writes to disk)
  !*  AUTHOR
  !*    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

     integer :: i_atom, i_e, i_k, i_k_point, i_l, i_spin_state, i_spin_proj, i_state 
     real*8  :: de, en, E1, E2

     character*120     :: info_str, outputformat
     character*50      :: proj_dos_filename
     character*120     :: outputformat_header
     character(LEN=80) :: uuid_str

     real*8, dimension(:,:,:,:), allocatable :: KS_dos
     real*8, dimension(:,:,:),   allocatable :: KS_sum_dos
     real*8, dimension(:),       allocatable :: proj_dos_erfs

     character(*), parameter :: func = 'write_atom_proj_dos'

     !   VB: In my view, this call to check_norm is incorrect - at the very least,
     !   it will now no longer work when two chemical potentials (spin constraint) are employed.
     !   call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_spin)

     call localorb_info(" ")
     write(info_str,'(2X,A)') 'Calculating atom-projected density of states ...'
     call localorb_info(info_str)
     write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
     call localorb_info(info_str)
     
     if(.not. allocated (KS_dos)) then
        allocate(KS_dos(0:l_wave_max,n_spin_states,atom_dos_n_en_points, n_atoms),stat=i_state) 
        call check_allocation(i_state, 'KS_dos                        ')
     endif
     if(.not. allocated (KS_sum_dos)) then
        allocate(KS_sum_dos(n_spin_states,atom_dos_n_en_points, n_atoms),stat=i_state) 
        call check_allocation(i_state, 'KS_sum_dos                    ')
     endif
     if (.not. allocated (proj_dos_erfs)) then
        allocate(proj_dos_erfs(atom_dos_n_en_points),stat=i_state)
        call check_allocation(i_state, 'proj_dos_erfs                 ')
     end if
     
     de= (atom_dos_high_energy - atom_dos_low_energy)/dble(atom_dos_n_en_points-1)
     KS_dos =0.d0 
     i_k = 0
     do i_k_point = 1, n_k_points
        if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
           i_k = i_k + 1
           do i_spin_state = 1, n_spin_states
              do i_state = 1, n_states_in
                 do i_e = 1, atom_dos_n_en_points
                    en= atom_dos_low_energy + dble(i_e-1)*de   
                    E1 = en - de/2d0
                    E2 = en + de/2d0
                    proj_dos_erfs(i_e) = &
                          arch_erf((E2-(KS_eigenvalue(i_state,i_spin_state,i_k_point))*hartree)*one_over_sqrt2/atom_dos_alpha)  &
                         -arch_erf((E1-(KS_eigenvalue(i_state,i_spin_state,i_k_point))*hartree)*one_over_sqrt2/atom_dos_alpha)   
                 end do
                 proj_dos_erfs(:) = proj_dos_erfs(:)*k_weights(i_k_point)/(2d0*dE)
                 do i_e = 1, atom_dos_n_en_points
                    do i_atom = 1, n_atoms
                       if(species_pseudoized(species(i_atom))) cycle
                       if (mode.eq.MULLIKEN_SR) then
                          i_spin_proj = i_spin_state
                          do i_l = 0, l_shell_max(species(i_atom))
                             KS_dos(i_l, i_spin_state,i_e,i_atom) =  KS_dos (i_l,i_spin_state,i_e,i_atom) + &
                                  spin_degeneracy * mulliken_decomp(i_l,i_atom,i_state,i_spin_proj,i_k) * proj_dos_erfs(i_e)
                          end do
                       else if (mode .eq. MULLIKEN_SOC .and. i_spin_state .eq. 1) then
                          ! The second condition is paranoia; n_spin_states should equal 1 for SOC calculations
                          do i_l = 0, l_shell_max(species(i_atom))
                             KS_dos(i_l, i_spin_state,i_e,i_atom) =  KS_dos (i_l,i_spin_state,i_e,i_atom) + &
                                  (mulliken_decomp(i_l,i_atom,i_state,1,i_k) + mulliken_decomp(i_l,i_atom,i_state,1,i_k)) & 
                                  * proj_dos_erfs(i_e)
                          end do
                       else
                          call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
                       end if
                    end do
                 end do
              end do
           end do
        end if
     end do
     call sync_vector( KS_dos(0,1,1,1),(1+l_wave_max)*n_spin_states*atom_dos_n_en_points*n_atoms)
     call sync_vector( KS_sum_dos(1,1,1),n_spin_states*atom_dos_n_en_points*n_atoms)
  
     KS_sum_dos =0.d0 
     do i_spin_state = 1, n_spin_states
        do i_e = 1, atom_dos_n_en_points
           do i_atom = 1, n_atoms
            if(species_pseudoized(species(i_atom))) cycle
              do i_l = 0, l_shell_max(species(i_atom))
                 KS_sum_dos(i_spin_state,i_e,i_atom) =  KS_sum_dos(i_spin_state,i_e,i_atom) + KS_dos(i_l, i_spin_state,i_e,i_atom)
              end do
           end do
        end do
     end do
     
     if (myid.eq.0) then
        if (n_spin_states.eq.1) then
           do i_atom = 1, n_atoms
            if(species_pseudoized(species(i_atom))) cycle            
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A11)') 'atom_projected_dos_',trim(species_name(species(i_atom))),'000',&
                      i_atom,adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A11)') 'atom_projected_dos_',trim(species_name(species(i_atom))),'00',&
                      i_atom,adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A11)') 'atom_projected_dos_',trim(species_name(species(i_atom))),'0',&
                      i_atom,adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A11)') 'atom_projected_dos_',trim(species_name(species(i_atom))),&
                      i_atom,adjustl(sr_suffix)
              end if
              
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(species(i_atom))+3
              write(use_unit,'(2X,5A)') '| writing projected DOS (shifted by the chemical potential) for species ',& 
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
              write(88,'(3A)') '# Angular momentum resolved density for species ',& 
                   trim(species_name(species(i_atom))), 'as calculated by FHI-aims'
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                   chemical_potential*hartree, ' eV'
              
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A4,A11)') 'atom_proj_dos_',trim(species_name(species(i_atom))),'000',i_atom,&
                   '_raw', adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A4,A11)') 'atom_proj_dos_',trim(species_name(species(i_atom))),'00',i_atom,&
                   '_raw', adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A4,A11)') 'atom_proj_dos_',trim(species_name(species(i_atom))),'0',i_atom,&
                   '_raw', adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A4,A11)') 'atom_proj_dos_',trim(species_name(species(i_atom))),i_atom,&
                   '_raw', adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') '| writing projected DOS (raw data) for species ',& 
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
              write(89,'(3A)') '# Angular momentum resolved density for species ',& 
                   trim(species_name(species(i_atom))), 'as calculated by FHI-aims'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level '
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(species(i_atom))+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
              
              do i_e = 1, atom_dos_n_en_points ,1
                 en = atom_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_atom), &
                      (KS_dos(i_l, 1,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
                 write(89,outputformat) en, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), & 
                      i_l=0,l_shell_max(species(i_atom)))
              enddo
              close(88)
              close(89)
           end do
        else 
           do i_atom = 1, n_atoms
  
             if(species_pseudoized(species(i_atom))) cycle
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'000',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'00',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'0',i_atom, adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),i_atom, adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
  
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'000',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'00',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'0',i_atom, adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),i_atom, adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)           
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
  
              write(88,'(3A)') '# spin-up angular momentum projected density of states for species ',& 
                   species_name(species(i_atom)), ' as calculated by FHI-aims '
              write(89,'(3A)') '# spin-down angular momentum projected density of states for species ',& 
                   species_name(species(i_atom)), ' as calculated by FHI-aims '
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A4,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'000',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A4,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'00',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A4,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),'0',i_atom,'_raw', adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A4,A11)') & 
                      'atom_proj_dos_spin_up',trim(species_name(species(i_atom))),i_atom,'_raw', adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') &
                 '| writing spin-up projected DOS (raw data) for species ', &
                 trim(species_name(species(i_atom))),&
                 ' to file ',trim(proj_dos_filename),'.'
              open(90, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(90,'(A,2X,A)') '#', uuid_str
  
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A4,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'000',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A4,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'00',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A4,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),'0',i_atom,'_raw', adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A4,A11)') & 
                      'atom_proj_dos_spin_dn',trim(species_name(species(i_atom))),i_atom,'_raw', adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') &
                 '| writing spin-down projected DOS (raw data) for species ', &
                 trim(species_name(species(i_atom))), &
                 ' to file ',trim(proj_dos_filename),'.'
              open(91, file=proj_dos_filename)                      
              call write_aims_uuid(uuid_str)
              write(91,'(A,2X,A)') '#', uuid_str
  
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(species(i_atom))+3
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                   chemical_potential*hartree, ' eV'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                   chemical_potential*hartree, ' eV'
              write(90,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              write(91,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(species(i_atom))+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
              write(90,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
              write(91,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
              
              do i_e = 1, atom_dos_n_en_points ,1
                 en = atom_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), &
                      i_l=0,l_shell_max(species(i_atom)))
                 write(89,outputformat) en-chemical_potential*hartree, KS_sum_dos(2,i_e,i_atom), (KS_dos(i_l, 2,i_e,i_atom), &
                      i_l=0,l_shell_max(species(i_atom)))
                 write(90,outputformat) en, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
                 write(91,outputformat) en, KS_sum_dos(2,i_e,i_atom), (KS_dos(i_l, 2,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
              enddo
              close(89)
              close(88) 
              close(90)
              close(91)
           end do
        end if
        
     end if  ! output on master thread
    
     if (allocated (KS_dos)        ) deallocate(KS_dos)
     if (allocated (KS_sum_dos)    ) deallocate(KS_sum_dos)
     if (allocated (proj_dos_erfs) ) deallocate(proj_dos_erfs)
  end subroutine write_atom_proj_dos
  !******

  !****f* mulliken/write_species_proj_dos_tetrahedron
  !*  NAME
  !*    write_species_proj_dos_tetrahedron
  !*  SYNOPSIS
  subroutine write_species_proj_dos_tetrahedron( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                          n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
  !*  PURPOSE
  !*    This subroutine calculates the projected DOS predicted by the decomposition to be associated with 
  !*    angular channels, summed up over all atoms of a given species, and outputs the results to disk.
  !*  USES
     use dimensions,            only : n_atoms, l_wave_max, n_k_points, n_k_points_task, n_species, &
                                       spin_degeneracy, ik2irred_map
     use constants,             only : hartree, one_over_sqrt2
     use localorb_io,           only : use_unit, localorb_info
     use mpi_tasks,             only : myid, n_tasks, aims_stop, check_allocation, mpi_comm_global
     use geometry,              only : species, lattice_vector
     use species_data,          only : l_shell_max, species_name, species_pseudoized
     use runtime_choices,       only : l_proj_dos_high_energy, l_proj_dos_low_energy,  &
                                       l_proj_dos_n_en_points, l_proj_dos_alpha, &
                                       n_k_points_xyz_nosym
     use pbc_lists,             only : k_weights
     use arch_specific,         only : arch_erf
     use synchronize_mpi_basic, only : sync_vector
     use generate_aims_uuid,    only : write_aims_uuid
     use tetrahedron_integration, only : ltispectral
  !*  ARGUMENTS
     implicit none

     integer,                                                   intent(in) :: n_states_in
     integer,                                                   intent(in) :: n_spin_states
     real*8, dimension(n_states_in, n_spin_states, n_k_points), intent(in) :: KS_eigenvalue
     real*8,                                                    intent(in) :: chemical_potential
     integer,                                                   intent(in) :: n_spin_proj
     real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
                                                                intent(in) :: mulliken_decomp
     character*11,                                              intent(in) :: sr_suffix
     integer,                                                   intent(in) :: mode
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvalue           -- Eigenvalues array
  !*    o chemical_potential      -- Fermi level/chemical potential of electrons
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*    o sr_suffix               -- Suffix to put on output file names for scalar-relativistic calculations
  !*                                 This is a legacy argument that should absorbed into the body of this
  !*                                 subroutine
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (writes to disk)
  !*  AUTHOR
  !*    Yi Yao
  !*    modified from  write_species_proj_dos
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

     integer :: i_atom, i_e, i_k, i_k_point, i_l, i_spin_state, i_spin_proj, i_state, i_species 
     integer :: i_x, i_y, i_z
     real*8  :: de, en, E1, E2
     real*8  :: energies(l_proj_dos_n_en_points)
     real*8  :: eigs(n_states_in,n_k_points_xyz_nosym(1),n_k_points_xyz_nosym(2),n_k_points_xyz_nosym(3),n_spin_states)
     real*8  :: funcs(n_states_in,n_k_points_xyz_nosym(1),n_k_points_xyz_nosym(2),n_k_points_xyz_nosym(3),n_spin_states,n_species)

     character*120     :: info_str, outputformat
     character*50      :: proj_dos_filename
     character*120     :: outputformat_header
     character(LEN=80) :: uuid_str

     real*8, dimension(:,:,:,:), allocatable :: KS_dos
     real*8, dimension(:,:,:),   allocatable :: KS_sum_dos
     real*8, dimension(:),       allocatable :: proj_dos_erfs

     real*8, dimension(:,:,:,:,:), allocatable :: state_mulliken_decomp
     real*8, dimension(:,:,:,:),   allocatable :: state_mulliken_decomp_sum
    
     character(*), parameter :: func = 'write_species_proj_dos_tetrahedron'

     !   VB: In my view, this call to check_norm is incorrect - at the very least,
     !   it will now no longer work when two chemical potentials (spin constraint) are employed.
     !   call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_spin)
  
     if(.not. allocated (KS_dos)) then
        allocate(KS_dos(0:l_wave_max,n_spin_states,l_proj_dos_n_en_points, n_species),stat=i_state) 
        call check_allocation(i_state, 'KS_dos                        ')
     endif
     if(.not. allocated (KS_sum_dos)) then
        allocate(KS_sum_dos(n_spin_states,l_proj_dos_n_en_points, n_species),stat=i_state) 
        call check_allocation(i_state, 'KS_sum_dos                    ')
     endif
     if (.not.allocated(proj_dos_erfs)) then
        allocate(proj_dos_erfs(l_proj_dos_n_en_points),stat=i_state)
        call check_allocation(i_state, 'proj_dos_erfs                 ')
     end if

     if(.not. allocated (state_mulliken_decomp)) then
        allocate(state_mulliken_decomp(n_states_in,0:l_wave_max,n_spin_states, n_species, n_k_points),stat=i_state) 
        call check_allocation(i_state, 'state_mulliken_decomp         ')
     endif
     if(.not. allocated (state_mulliken_decomp_sum)) then
        allocate(state_mulliken_decomp_sum(n_states_in, n_spin_states, n_species, n_k_points),stat=i_state) 
        call check_allocation(i_state, 'state_mulliken_decomp_sum     ')
     endif
     
     write(info_str,'(2X,A)') 'Calculating angular momentum projected density of states ...'
     call localorb_info(info_str)
     write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
     call localorb_info(info_str)
     
     de= (l_proj_dos_high_energy - l_proj_dos_low_energy)/dble(l_proj_dos_n_en_points-1)
     do i_e = 1, l_proj_dos_n_en_points
       energies(i_e) = l_proj_dos_low_energy + de * (i_e - 1)
     enddo
     KS_dos =0.d0 
     KS_sum_dos =0.d0 
     state_mulliken_decomp = 0.d0
     state_mulliken_decomp_sum = 0.d0
     
     ! compute DOS for angular momentum projection. Make sure that the mulliken decomposition is properly normalized!!!
     i_k = 0
     do i_k_point = 1, n_k_points
        ! YY: why myid <= n_k_points ( I guess related to how mulliken parallelized)
        if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
           i_k = i_k + 1
           do i_spin_state = 1, n_spin_states
              do i_state = 1, n_states_in
                    do i_atom = 1, n_atoms, 1
                       if(species_pseudoized(species(i_atom))) cycle
                       if (mode.eq.MULLIKEN_SR) then
                          i_spin_proj = i_spin_state
                          do i_l = 0, l_shell_max(species(i_atom))
                             state_mulliken_decomp(i_state, i_l, i_spin_state,species(i_atom),i_k_point) =  state_mulliken_decomp(i_state,i_l,i_spin_state,species(i_atom),i_k_point) + &
                                  mulliken_decomp(i_l,i_atom,i_state,i_spin_proj,i_k) 
                          end do
                       !else if (mode .eq. MULLIKEN_SOC .and. i_spin_state .eq. 1) then
                       !   ! The second condition is paranoia; n_spin_states should equal 1 for SOC calculations
                       !   do i_l = 0, l_shell_max(species(i_atom))
                       !      KS_dos(i_l, i_spin_state,i_e,species(i_atom)) = KS_dos (i_l,i_spin_state,i_e,species(i_atom)) + &
                       !           (mulliken_decomp(i_l,i_atom,i_state,1,i_k) + mulliken_decomp(i_l,i_atom,i_state,2,i_k)) &
                       !           * proj_dos_erfs(i_e)
                       !   end do
                       else
                          call aims_stop('tetrahedron pdos with SOC not implemented.', func) 
                          call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
                       end if
                    end do
              end do
           end do
        end if
     end do
     call sync_vector( KS_dos(0,1,1,1),(1+l_wave_max)*n_spin_states*l_proj_dos_n_en_points*n_species)
     call sync_vector( KS_sum_dos(1,1,1),n_spin_states*l_proj_dos_n_en_points*n_species)
     call sync_vector( state_mulliken_decomp(1,0,1,1,1),n_states_in*(1+l_wave_max)*n_spin_states*n_species*n_k_points)
     call sync_vector( state_mulliken_decomp_sum(1,1,1,1),n_states_in*n_spin_states*n_species*n_k_points)
  
     do i_spin_state = 1, n_spin_states, 1
        do i_state = 1, n_states_in, 1
           do i_species = 1, n_species
              do i_l = 0, l_shell_max(i_species)
                 do i_k = 1, n_k_points
                 state_mulliken_decomp_sum(i_state,i_spin_state,i_species,i_k) =  &
                    state_mulliken_decomp_sum(i_state,i_spin_state,i_species,i_k) + &
                    state_mulliken_decomp(i_state,i_l,i_spin_state,i_species,i_k)
                 end do
              end do
           end do
        end do
     end do

     do i_x = 1, n_k_points_xyz_nosym(1)
       do i_y = 1, n_k_points_xyz_nosym(2)
           do i_z = 1, n_k_points_xyz_nosym(3)
             i_k_point = ik2irred_map(i_x,i_y,i_z)
             do i_spin_state = 1, n_spin_states
               eigs(:,i_x,i_y,i_z,i_spin_state) = KS_eigenvalue(:,i_spin_state,i_k_point) * hartree
               do i_species = 1, n_species
                 funcs(:,i_x,i_y,i_z,i_spin_state,i_species) = state_mulliken_decomp_sum(:,i_spin_state,i_species,i_k_point)
               end do
             end do
           enddo
       enddo
     enddo
     do i_species = 1, n_species
       do i_spin_state = 1, n_spin_states
         call ltispectral(lattice_vector, n_k_points_xyz_nosym(1), n_k_points_xyz_nosym(2), &
                     n_k_points_xyz_nosym(3), n_states_in, &
                     eigs(:,:,:,:,i_spin_state), funcs(:,:,:,:,i_spin_state,i_species), l_proj_dos_n_en_points,&
                     energies, KS_sum_dos(i_spin_state,:,i_species), mpi_comm_global)
         
       end do
     end do
     KS_sum_dos = KS_sum_dos * spin_degeneracy
  
     ! output on thread 0 only 
     if (myid.eq.0) then
        if (n_spin_states.eq.1) then
           do i_species = 1, n_species
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_tetrahedron', trim(sr_suffix)
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(i_species)+3
              write(use_unit,'(2X,5A)') '| writing projected DOS (shifted by the chemical potential) for species ',& 
                   trim(species_name(i_species)), &
                   ' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
              write(88,'(3A)') '# Angular momentum resolved density for species ',& 
                trim(species_name(i_species)), 'as calculated by FHI-aims'
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                chemical_potential*hartree, ' eV'
              
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_raw_tetrahedron', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing projected DOS (raw data) for species ',& 
                trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
              write(89,'(3A)') '# Angular momentum resolved density for species ',& 
                trim(species_name(i_species)), 'as calculated by FHI-aims'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level '
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(i_species)+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
              
              do i_e = 1, l_proj_dos_n_en_points ,1
                 en = l_proj_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_species), &
                      (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(89,outputformat) & 
                      en, KS_sum_dos(1,i_e,i_species), (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
              enddo
              close(88)
              close(89)
           end do
        else 
           do i_species = 1, n_species
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_up_tetrahedron', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_down_tetrahedron', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(i_species)),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)           
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
              write(88,'(3A)') '# spin-up angular momentum projected density of states for species ',& 
                   species_name(i_species), ' as calculated by FHI-aims '
              write(89,'(3A)') '# spin-down angular momentum projected density of states for species ',& 
                   species_name(i_species), ' as calculated by FHI-aims '
              
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_up_raw_tetrahedron', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (raw data) for species ',trim(species_name(i_species)),&
                   ' to file ',trim(proj_dos_filename),'.'
              open(90, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(90,'(A,2X,A)') '#', uuid_str
              write(proj_dos_filename,'(3A)') trim(species_name(i_species)), '_l_proj_dos_spin_down_raw_tetrahedron', trim(sr_suffix)
              write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (raw data) for species ',trim(species_name(i_species)),&
                   ' to file ',trim(proj_dos_filename),'.'
              open(91, file=proj_dos_filename)                      
              call write_aims_uuid(uuid_str)
              write(91,'(A,2X,A)') '#', uuid_str
  
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(i_species)+3
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ',& 
                   chemical_potential*hartree,' eV'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ',&
                   chemical_potential*hartree,' eV'
              write(90,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              write(91,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(i_species)+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
              write(90,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))
              write(91,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(i_species))   
              
              do i_e = 1, l_proj_dos_n_en_points ,1
                 en = l_proj_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_species), &
                      (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(89,outputformat) en-chemical_potential*hartree, KS_sum_dos(2,i_e,i_species), &
                      (KS_dos(i_l, 2,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(90,outputformat) en, KS_sum_dos(1,i_e,i_species), (KS_dos(i_l, 1,i_e,i_species), i_l=0,l_shell_max(i_species))
                 write(91,outputformat) en, KS_sum_dos(2,i_e,i_species), (KS_dos(i_l, 2,i_e,i_species), i_l=0,l_shell_max(i_species))
              enddo
              close(89)
              close(88)  
              close(90)
              close(91)
           end do
        end if
        
     end if
     call localorb_info(" ")
     if (allocated(KS_dos)       ) deallocate(KS_dos)        
     if (allocated(KS_sum_dos)   ) deallocate(KS_sum_dos)    
     if (allocated(proj_dos_erfs)) deallocate(proj_dos_erfs)
     if (allocated(state_mulliken_decomp)       ) deallocate(state_mulliken_decomp)        
     if (allocated(state_mulliken_decomp_sum)       ) deallocate(state_mulliken_decomp_sum)        
  end subroutine write_species_proj_dos_tetrahedron
  !******

  !****f* mulliken/write_atom_proj_dos_tetrahedron
  !*  NAME
  !*    write_atom_proj_dos_tetrahedron
  !*  SYNOPSIS
  subroutine write_atom_proj_dos_tetrahedron( n_states_in, n_spin_states, KS_eigenvalue, chemical_potential, &
                                          n_spin_proj, mulliken_decomp, sr_suffix, mode ) 
  !*  PURPOSE
  !*    This subroutine calculates the projected DOS predicted by the decomposition to be associated with 
  !*    angular channels, summed up over all atoms of a given species, and outputs the results to disk.
  !*  USES
     use dimensions,            only : n_atoms, l_wave_max, n_k_points, n_k_points_task, n_species, &
                                       spin_degeneracy, ik2irred_map
     use constants,             only : hartree, one_over_sqrt2
     use localorb_io,           only : use_unit, localorb_info
     use mpi_tasks,             only : myid, n_tasks, aims_stop, check_allocation, mpi_comm_global
     use geometry,              only : species, lattice_vector
     use species_data,          only : l_shell_max, species_name, species_pseudoized
     use runtime_choices,       only : atom_dos_high_energy, atom_dos_low_energy,  &
                                       atom_dos_n_en_points, atom_dos_alpha, &
                                       n_k_points_xyz_nosym
     use pbc_lists,             only : k_weights
     use arch_specific,         only : arch_erf
     use synchronize_mpi_basic, only : sync_vector
     use generate_aims_uuid,    only : write_aims_uuid
     use tetrahedron_integration, only : ltispectral
  !*  ARGUMENTS
     implicit none

     integer,                                                   intent(in) :: n_states_in
     integer,                                                   intent(in) :: n_spin_states
     real*8, dimension(n_states_in, n_spin_states, n_k_points), intent(in) :: KS_eigenvalue
     real*8,                                                    intent(in) :: chemical_potential
     integer,                                                   intent(in) :: n_spin_proj
     real*8, dimension(0:l_wave_max, n_atoms, n_states_in, n_spin_proj, n_k_points_task), &
                                                                intent(in) :: mulliken_decomp
     character*11,                                              intent(in) :: sr_suffix
     integer,                                                   intent(in) :: mode
  !*  INPUT
  !*    o n_states_in             -- Number of states in eigenvalue and occupation numbers array
  !*                                 Since eigenvalues are stored in full for every process, this will be n_states
  !*                                 or some variant thereof
  !*    o n_spin_states           -- Number of spin indices in eigenvector and eigenvalue arrays
  !*                                 May differ from n_spin for different physical objects (e.g. SOC eigenvectors)
  !*    o KS_eigenvalue           -- Eigenvalues array
  !*    o chemical_potential      -- Fermi level/chemical potential of electrons
  !*    o n_spin_proj             -- The number of spin channels in the Mulliken decomposition
  !*                                 For collinear calculations, this will be identical to n_spin_states
  !*                                 For non-collinear calculations, this will be 2 (maybe 4 one day...)
  !*    o mulliken_decomp         -- The Mulliken decomposition onto atoms, angular channels, and spin channels
  !*    o sr_suffix               -- Suffix to put on output file names for scalar-relativistic calculations
  !*                                 This is a legacy argument that should absorbed into the body of this
  !*                                 subroutine
  !*    o mode                    -- A generic flag to alter the details of the calculations.  Options are
  !*                                 module parameters (see module header for more information).
  !*  OUTPUT
  !*    none (writes to disk)
  !*  AUTHOR
  !*    Yi Yao
  !*    modified from  write_atom_proj_dos
  !*  SEE ALSO
  !*    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !*    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !*    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !*    Computer Physics Communications (2008), submitted.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

     integer :: i_atom, i_e, i_k, i_k_point, i_l, i_spin_state, i_spin_proj, i_state, i_species 
     integer :: i_x, i_y, i_z
     real*8  :: de, en, E1, E2
     real*8  :: energies(atom_dos_n_en_points)
     real*8  :: eigs(n_states_in,n_k_points_xyz_nosym(1),n_k_points_xyz_nosym(2),n_k_points_xyz_nosym(3),n_spin_states)
     real*8  :: funcs(n_states_in,n_k_points_xyz_nosym(1),n_k_points_xyz_nosym(2),n_k_points_xyz_nosym(3),n_spin_states,n_atoms)

     character*120     :: info_str, outputformat
     character*50      :: proj_dos_filename
     character*120     :: outputformat_header
     character(LEN=80) :: uuid_str

     real*8, dimension(:,:,:,:), allocatable :: KS_dos
     real*8, dimension(:,:,:),   allocatable :: KS_sum_dos
     real*8, dimension(:),       allocatable :: proj_dos_erfs

     real*8, dimension(:,:,:,:,:), allocatable :: state_mulliken_decomp
     real*8, dimension(:,:,:,:),   allocatable :: state_mulliken_decomp_sum
    
     character(*), parameter :: func = 'write_atom_proj_dos_tetrahedron'

     !   VB: In my view, this call to check_norm is incorrect - at the very least,
     !   it will now no longer work when two chemical potentials (spin constraint) are employed.
     !   call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_spin)
  
     if(.not. allocated (KS_dos)) then
        allocate(KS_dos(0:l_wave_max,n_spin_states,atom_dos_n_en_points, n_atoms),stat=i_state) 
        call check_allocation(i_state, 'KS_dos                        ')
     endif
     if(.not. allocated (KS_sum_dos)) then
        allocate(KS_sum_dos(n_spin_states,atom_dos_n_en_points, n_atoms),stat=i_state) 
        call check_allocation(i_state, 'KS_sum_dos                    ')
     endif
     if (.not.allocated(proj_dos_erfs)) then
        allocate(proj_dos_erfs(atom_dos_n_en_points),stat=i_state)
        call check_allocation(i_state, 'proj_dos_erfs                 ')
     end if

     if(.not. allocated (state_mulliken_decomp)) then
        allocate(state_mulliken_decomp(n_states_in,0:l_wave_max,n_spin_states, n_atoms, n_k_points),stat=i_state) 
        call check_allocation(i_state, 'state_mulliken_decomp         ')
     endif
     if(.not. allocated (state_mulliken_decomp_sum)) then
        allocate(state_mulliken_decomp_sum(n_states_in, n_spin_states, n_atoms, n_k_points),stat=i_state) 
        call check_allocation(i_state, 'state_mulliken_decomp_sum     ')
     endif
     
     write(info_str,'(2X,A)') 'Calculating angular momentum projected density of states ...'
     call localorb_info(info_str)
     write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
     call localorb_info(info_str)
     
     de= (atom_dos_high_energy - atom_dos_low_energy)/dble(atom_dos_n_en_points-1)
     do i_e = 1, atom_dos_n_en_points
       energies(i_e) = atom_dos_low_energy + de * (i_e - 1)
     enddo
     KS_dos =0.d0 
     KS_sum_dos =0.d0 
     state_mulliken_decomp = 0.d0
     state_mulliken_decomp_sum = 0.d0
     
     ! compute DOS for angular momentum projection. Make sure that the mulliken decomposition is properly normalized!!!
     i_k = 0
     do i_k_point = 1, n_k_points
        ! YY: why myid <= n_k_points ( I guess related to how mulliken parallelized)
        if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
           i_k = i_k + 1
           do i_spin_state = 1, n_spin_states
              do i_state = 1, n_states_in
                    do i_atom = 1, n_atoms, 1
                       if(species_pseudoized(species(i_atom))) cycle
                       if (mode.eq.MULLIKEN_SR) then
                          i_spin_proj = i_spin_state
                          do i_l = 0, l_shell_max(species(i_atom))
                             state_mulliken_decomp(i_state, i_l, i_spin_state,i_atom,i_k_point) =  state_mulliken_decomp(i_state,i_l,i_spin_state,i_atom,i_k_point) + &
                                  mulliken_decomp(i_l,i_atom,i_state,i_spin_proj,i_k) 
                          end do
                       !else if (mode .eq. MULLIKEN_SOC .and. i_spin_state .eq. 1) then
                       !   ! The second condition is paranoia; n_spin_states should equal 1 for SOC calculations
                       !   do i_l = 0, l_shell_max(species(i_atom))
                       !      KS_dos(i_l, i_spin_state,i_e,species(i_atom)) = KS_dos (i_l,i_spin_state,i_e,species(i_atom)) + &
                       !           (mulliken_decomp(i_l,i_atom,i_state,1,i_k) + mulliken_decomp(i_l,i_atom,i_state,2,i_k)) &
                       !           * proj_dos_erfs(i_e)
                       !   end do
                       else
                          call aims_stop('tetrahedron pdos with SOC not implemented.', func) 
                          call aims_stop('An incorrect value was specified for the mode of operation, exiting.', func) 
                       end if
                    end do
              end do
           end do
        end if
     end do
     call sync_vector( KS_dos(0,1,1,1),(1+l_wave_max)*n_spin_states*atom_dos_n_en_points*n_atoms)
     call sync_vector( KS_sum_dos(1,1,1),n_spin_states*atom_dos_n_en_points*n_atoms)
     call sync_vector( state_mulliken_decomp(1,0,1,1,1),n_states_in*(1+l_wave_max)*n_spin_states*n_atoms*n_k_points)
     call sync_vector( state_mulliken_decomp_sum(1,1,1,1),n_states_in*n_spin_states*n_atoms*n_k_points)
  
     do i_spin_state = 1, n_spin_states, 1
        do i_state = 1, n_states_in, 1
           do i_atom = 1, n_atoms
              do i_l = 0, l_shell_max(species(i_atom))
                 do i_k = 1, n_k_points
                 state_mulliken_decomp_sum(i_state,i_spin_state,i_atom,i_k) =  &
                    state_mulliken_decomp_sum(i_state,i_spin_state,i_atom,i_k) + &
                    state_mulliken_decomp(i_state,i_l,i_spin_state,i_atom,i_k)
                 end do
              end do
           end do
        end do
     end do

     do i_x = 1, n_k_points_xyz_nosym(1)
       do i_y = 1, n_k_points_xyz_nosym(2)
           do i_z = 1, n_k_points_xyz_nosym(3)
             i_k_point = ik2irred_map(i_x,i_y,i_z)
             do i_spin_state = 1, n_spin_states
               eigs(:,i_x,i_y,i_z,i_spin_state) = KS_eigenvalue(:,i_spin_state,i_k_point) * hartree
               do i_atom = 1, n_atoms
                 funcs(:,i_x,i_y,i_z,i_spin_state,i_atom) = state_mulliken_decomp_sum(:,i_spin_state,i_atom,i_k_point)
               end do
             end do
           enddo
       enddo
     enddo
     do i_atom = 1, n_atoms
       do i_spin_state = 1, n_spin_states
         call ltispectral(lattice_vector, n_k_points_xyz_nosym(1), n_k_points_xyz_nosym(2), &
                     n_k_points_xyz_nosym(3), n_states_in, &
                     eigs(:,:,:,:,i_spin_state), funcs(:,:,:,:,i_spin_state,i_atom), atom_dos_n_en_points,&
                     energies, KS_sum_dos(i_spin_state,:,i_atom), mpi_comm_global)
         
       end do
     end do
     KS_sum_dos = KS_sum_dos * spin_degeneracy
  

     if (myid.eq.0) then
        if (n_spin_states.eq.1) then
           do i_atom = 1, n_atoms
            if(species_pseudoized(species(i_atom))) cycle            
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A11)') 'atom_projected_dos_tetrahedron',trim(species_name(species(i_atom))),'000',&
                      i_atom,adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A11)') 'atom_projected_dos_tetrahedron',trim(species_name(species(i_atom))),'00',&
                      i_atom,adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A11)') 'atom_projected_dos_tetrahedron',trim(species_name(species(i_atom))),'0',&
                      i_atom,adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A11)') 'atom_projected_dos_tetrahedron',trim(species_name(species(i_atom))),&
                      i_atom,adjustl(sr_suffix)
              end if
              
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(species(i_atom))+3
              write(use_unit,'(2X,5A)') '| writing projected DOS (shifted by the chemical potential) for species ',& 
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
              write(88,'(3A)') '# Angular momentum resolved density for species ',& 
                   trim(species_name(species(i_atom))), 'as calculated by FHI-aims'
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                   chemical_potential*hartree, ' eV'
              
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A4,A11)') 'atom_proj_dos_tetrahedron',trim(species_name(species(i_atom))),'000',i_atom,&
                   '_raw', adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A4,A11)') 'atom_proj_dos_tetrahedron',trim(species_name(species(i_atom))),'00',i_atom,&
                   '_raw', adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A4,A11)') 'atom_proj_dos_tetrahedron',trim(species_name(species(i_atom))),'0',i_atom,&
                   '_raw', adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A4,A11)') 'atom_proj_dos_tetrahedron',trim(species_name(species(i_atom))),i_atom,&
                   '_raw', adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') '| writing projected DOS (raw data) for species ',& 
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
              write(89,'(3A)') '# Angular momentum resolved density for species ',& 
                   trim(species_name(species(i_atom))), 'as calculated by FHI-aims'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level '
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(species(i_atom))+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
              
              do i_e = 1, atom_dos_n_en_points ,1
                 en = atom_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_atom), &
                      (KS_dos(i_l, 1,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
                 write(89,outputformat) en, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), & 
                      i_l=0,l_shell_max(species(i_atom)))
              enddo
              close(88)
              close(89)
           end do
        else 
           do i_atom = 1, n_atoms
  
             if(species_pseudoized(species(i_atom))) cycle
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),'000',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),'00',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),'0',i_atom, adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),i_atom, adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') '| writing spin-up projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(88, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(88,'(A,2X,A)') '#', uuid_str
  
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),'000',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),'00',i_atom, adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),'0',i_atom, adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),i_atom, adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') '| writing spin-down projected DOS (shifted by the chemical potential) for species ',&
                   trim(species_name(species(i_atom))),' to file ',trim(proj_dos_filename),'.'
              open(89, file=proj_dos_filename)           
              call write_aims_uuid(uuid_str)
              write(89,'(A,2X,A)') '#', uuid_str
  
              write(88,'(3A)') '# spin-up angular momentum projected density of states for species ',& 
                   species_name(species(i_atom)), ' as calculated by FHI-aims '
              write(89,'(3A)') '# spin-down angular momentum projected density of states for species ',& 
                   species_name(species(i_atom)), ' as calculated by FHI-aims '
              
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A4,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),'000',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A4,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),'00',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A4,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),'0',i_atom,'_raw', adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A4,A11)') & 
                      'atom_proj_dos_spin_up_tetrahedron',trim(species_name(species(i_atom))),i_atom,'_raw', adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') &
                 '| writing spin-up projected DOS (raw data) for species ', &
                 trim(species_name(species(i_atom))),&
                 ' to file ',trim(proj_dos_filename),'.'
              open(90, file=proj_dos_filename)
              call write_aims_uuid(uuid_str)
              write(90,'(A,2X,A)') '#', uuid_str
  
              if (i_atom.lt.10) then
                 write(proj_dos_filename,'(3A,I1,A4,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),'000',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.100) then
                 write(proj_dos_filename,'(3A,I2,A4,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),'00',i_atom,'_raw', adjustl(sr_suffix)
              else if (i_atom.lt.1000) then
                 write(proj_dos_filename,'(3A,I3,A4,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),'0',i_atom,'_raw', adjustl(sr_suffix)
              else
                 write(proj_dos_filename,'(2A,I4,A4,A11)') & 
                      'atom_proj_dos_spin_dn_tetrahedron',trim(species_name(species(i_atom))),i_atom,'_raw', adjustl(sr_suffix)
              end if
              
              write(use_unit,'(2X,5A)') &
                 '| writing spin-down projected DOS (raw data) for species ', &
                 trim(species_name(species(i_atom))), &
                 ' to file ',trim(proj_dos_filename),'.'
              open(91, file=proj_dos_filename)                      
              call write_aims_uuid(uuid_str)
              write(91,'(A,2X,A)') '#', uuid_str
  
              write(outputformat,'("(",I2,"(2X,F16.8))")') l_shell_max(species(i_atom))+3
              write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                   chemical_potential*hartree, ' eV'
              write(89,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', & 
                   chemical_potential*hartree, ' eV'
              write(90,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              write(91,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
              
              write(outputformat_header,'("(A18,2X,A16,",I2,"(A7,I2,9X))")') l_shell_max(species(i_atom))+1
              write(88,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
              write(89,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
              write(90,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))
              write(91,outputformat_header) "#      Energy (eV)","total dos",("     l=",i_l, i_l = 0, l_shell_max(species(i_atom)))   
              
              do i_e = 1, atom_dos_n_en_points ,1
                 en = atom_dos_low_energy + dble(i_e -1) *de 
                 write(88,outputformat) en-chemical_potential*hartree, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), &
                      i_l=0,l_shell_max(species(i_atom)))
                 write(89,outputformat) en-chemical_potential*hartree, KS_sum_dos(2,i_e,i_atom), (KS_dos(i_l, 2,i_e,i_atom), &
                      i_l=0,l_shell_max(species(i_atom)))
                 write(90,outputformat) en, KS_sum_dos(1,i_e,i_atom), (KS_dos(i_l, 1,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
                 write(91,outputformat) en, KS_sum_dos(2,i_e,i_atom), (KS_dos(i_l, 2,i_e,i_atom), i_l=0,l_shell_max(species(i_atom)))
              enddo
              close(89)
              close(88) 
              close(90)
              close(91)
           end do
        end if
        
     end if  ! output on master thread
     call localorb_info(" ")
     if (allocated(KS_dos)       ) deallocate(KS_dos)        
     if (allocated(KS_sum_dos)   ) deallocate(KS_sum_dos)    
     if (allocated(proj_dos_erfs)) deallocate(proj_dos_erfs)
     if (allocated(state_mulliken_decomp)       ) deallocate(state_mulliken_decomp)        
     if (allocated(state_mulliken_decomp_sum)       ) deallocate(state_mulliken_decomp_sum)        
  end subroutine write_atom_proj_dos_tetrahedron
  !******
 
end module mulliken
!******
