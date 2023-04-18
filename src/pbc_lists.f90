!****h*  FHI-aims/pbc_lists
!  NAME
!    pbc_lists
!  SYNOPSIS

module pbc_lists

  !  PURPOSE
  !
  !    Provides data needed to handle long-range sums and integrations for
  !    periodic boundary conditions OR (possibly) large molecules efficiently:
  !    * Lists of Hartree potential "centers"
  !    * Lists of integration centers
  !    * Lists of periodic images of basis functions
  !    * ...
  !    * k-point lists
  !    * ...
  !
  !  DESIGN
  !
  !    * Periodicity
  !
  !      In a periodic system, there are besides the main, "central" unit cell
  !      a bunch of (n_cells - 1) additional unit cells.  In practice, they are
  !      arranged in a parallel epipede sometimes called supercell.
  !
  !    * Atomic positions:
  !
  !      While "atom" generally refers to an atom within the central unit cell
  !      (0-cell) and is specified in coords(3, n_atoms), species(n_atoms), a
  !      general atom within the whole supercell is called "center" and saved
  !      in coords_center(3, n_centers), species_center(n_centers) with
  !      n_centers = n_atoms * n_cells.  One can obtain the corresponding atom
  !      in the central unit cell by center_to_atom(i_center) and the id of
  !      the unit cell by center_to_cell(i_center).  Please note that "atoms"
  !      are always in the central parallel epipede in the sense of
  !      map_to_center_cell().
  !
  !    * k-points
  !
  !      All what is needed from the k-points is the phases k_phases(n_cells,
  !      n_k_points) they generate in all the cells and the corresponding
  !      integration weight in k_weights(n_k_points).  Please note that the
  !      convention is to take the phase only from the cell vectors, ignoring
  !      the actual atomic positions.
  !
  !      For convenience, the explicit k-points are supplied in
  !      k_point_list(n_k_points,3).
  !
  !    * Integrations
  !
  !      All integrations over atomic quantities are formulated as
  !      integrations over the central, parallel epipede shaped, unit cell.
  !      The grid in this unit cell is constructed from all "atoms", and each
  !      single grid point is mapped back to the 0-cell by
  !      map_to_center_cell().
  !
  !      Depending on the radius of the quantities to integrate, different
  !      centers may touch the 0-cell and therefore grid points.  There are
  !      several index lists (centers_...) for efficient access to significant
  !      centers for a given type of interaction radius (see below).  There is
  !      also a list of all basis functions and there periodic replica
  !      (Cbasis) which touch some point in the 0-cell.  
  !
  !      JW: A note to those who (like me) are used to integrations over whole
  !      space like is more usually done for analytically given orbitals: the
  !      restriction to the 0-cell comes at the price that one has to really
  !      loop over all n_Cbasis functions or all (n_Cbasis, n_Cbasis) pairs
  !      instead of only n_basis or only (n_basis, n_Cbasis), respectively.
  !
  ! USES

  implicit none

  !                         k-points without symmetries, BB: should be defined here?
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
  !
  ! SOURCE

  ! VB: centers_lists and k-point lists moved here from geometry.f90 
  !     since all these are now handled by a separate module (initialize_lists) anyway,
  !     there is no point in keeping them in geometry.f90 - all these lists are strictly
  !     only needed for periodic boundary conditions ...

  !How many super cell different directions is maximally needed in periodic systems ?
  integer,dimension(3) ::       number_of_super_cells

  ! * Atomic positions (centers)

  ! Number of super cells.
  ! WARNING: While this is product(2*number_of_super_cells+1) during most
  ! parts of initialize_bc_dependent_lists(), it gets re-set to some smaller
  ! value in initialize_position_in_hamiltonian().
  integer:: n_cells


  integer :: n_cells_pairs ! SVL It is needed because the set of
                           ! active supercell can be different
                           ! in periodic EXX and in other calculations

  integer :: n_cells_bvk ! SVL A set of Born-von Karman cells

  ! Definition of centers:
  ! (center_to_atom(i_center), center_to_cell(i_center)) == (i_atom, i_cell)
  integer,  dimension(:),      allocatable:: center_to_atom  ! (n_centers)
  integer,  dimension(:),      allocatable:: center_to_cell  ! (n_centers)

  ! Convenience arrays:
  real*8, dimension(:,:), allocatable:: coords_center    ! coords_center(3,n_centers)
  integer, dimension(:), allocatable:: species_center    ! species_center(n_centers)

  integer,  dimension(:),      allocatable:: center_to_cell_pairs !SVL index of supercell to which a center belongs
  ! can be different from center_to_cell, needed for the exact exchange calculations


  ! * k-point related stuff

  ! Construct compelete phase factor from:
  !    k_phase(i_cell, i_k_point) &
  !    & = product(k_phase_base(:, i_k_point)**cell_index(i_cell,:))
  complex*16, allocatable :: k_phase_base(:,:)
  complex*16,dimension(:,:),   allocatable:: k_phase           ! k-points phase factor
  real*8,   dimension(:),      allocatable:: k_weights         ! weights to the k-points

  complex*16,dimension(:,:),   allocatable:: k_phase_exx       ! SVL k_phase for exact exchange
                                                               ! It is needed because the set of
                                                               ! active supercell can be different 
                                                               ! in periodic EXX and in other calculations
  real*8, dimension(:,:), allocatable :: k_point_list          ! For periodic EXX and for output
                                                               ! This information has been missing for a long time: 
                                                               ! List of k_point coordinates in relative units of reciprocal 
                                                               ! lattice vectors

  real*8, dimension(:,:), allocatable :: k_minus_q_point_list  ! An independent list of k-q points, required for band structure calculations with exact exchange
  real*8, dimension(:,:), allocatable :: irk_point_list  ! irreducible k point list due to symmetry
  real*8, dimension(:), allocatable :: irk_weight  ! irreducible k point weight
  integer, dimension(:), allocatable :: irk_point_mapping  ! mapping the k point list to the irreducile k point list 
  integer, dimension(:), allocatable :: inv_irk_point_mapping  ! mapping the irreducible k point list to full k point list 
  logical, dimension(:), allocatable :: irk_point_included  ! whether a given k point included in the irreducible set
  integer, dimension(:,:), allocatable:: k_point_loc        ! mapping between i_k_point and (i_task,i_k_point_local) location for the eigenvectors

  ! The index i_Cbasis is the generalization to i_basis much like i_center to
  ! i_atom.

  ! Cbasis_to_basis(i_Cbasis) == i_basis
  integer, dimension(:), allocatable:: Cbasis_to_basis    ! (n_Cbasis)
  integer, dimension(:), allocatable:: Cbasis_to_basis_s  ! For relativistic small comp. basis
  ! Cbasis_to_center(i_Cbasis) == i_center
  integer, dimension(:), allocatable:: Cbasis_to_center   ! (n_Cbasis)
  integer, dimension(:), allocatable:: Cbasis_to_center_s ! For relativistic small comp. basis
  ! Cbasis_to_center(i_Cbasis) == center_to_atom(Cbasis_to_center(i_Cbasis))
  integer, dimension(:), allocatable:: Cbasis_to_atom     ! (n_Cbasis)
  integer, dimension(:), allocatable:: Cbasis_to_atom_s   ! For relativistic small comp. basis

  ! * Different lists of centers

  ! '... touches atom in 0-cell ' means that it might have overlap to some
  ! basis function of an atom within the central unit cell.  All lists start
  ! with the atoms of the 0-cell.

  ! List of centers whose (screened) Hartree potential touches atom in 0-cell.
  integer, dimension(:), allocatable:: centers_hartree_potential
  ! List of centers whose partitioned multipole charge touches atom in 0-cell.
  integer, dimension(:), allocatable:: centers_hartree_multipole
  ! List of centers whose basis functions touch atom in 0-cell.
  integer, dimension(:), allocatable:: centers_ele_summation

  ! The next array actually saves two lists:
  ! List of centers whose free-atom charge touches point in 0-cell:
  !   centers_integrals := centers_basis_integrals(1:n_centers_integrals)
  ! List of centers whose free-atom charge touches point /or atom/ in 0-cell.
  integer, dimension(:), allocatable:: centers_basis_integrals
  ! inv_centers_basis_integrals(centers_basis_integrals(i_cbi)) == i_cbi
  ! if (i_center not in centers_basis_integrals) then
  !   inv_centers_basis_integrals(i_center) == 0
  integer, dimension(:), allocatable:: inv_centers_basis_integrals
  ! Array for "empty" atom sites (ghosts for BSSE correction)
  integer, dimension(:), allocatable:: new_centers_basis_integrals
 


  ! * Hamiltonian indexing

  ! See comments on hamiltonian(:,:) in physics.f90 for the following:
  integer:: n_cells_in_hamiltonian
  integer:: column_index_hamiltonian_size
  integer:: index_hamiltonian_dim2
  integer,dimension(:, :, :), allocatable:: index_hamiltonian
  integer,dimension(:), allocatable:: column_index_hamiltonian
  logical, private :: sparse_matrix_needs_reduce

  integer,dimension(:, :, :), allocatable:: index_hamiltonian_no_symmetry
  integer,dimension(:), allocatable:: column_index_hamiltonian_no_symmetry


  ! * Radii

  ! Range of the short-range (erfc) part of the Hartree pot.
  real*8, dimension(:), allocatable, private :: atom_radius_hartree  ! (n_species)
  real*8, dimension(:), allocatable :: atom_radius_hartree_sq        ! (n_species)

  ! * Cells

  integer, dimension(:,:), allocatable:: cell_index

  integer, dimension(:,:), allocatable :: cell_index_pairs !SVL It is needed because the set of
                                                           ! active supercell can be different
                                                           ! in periodic EXX and in other calculations

  integer, dimension(:,:), allocatable :: cell_index_bvk   ! SVL Index for Born-von Karman cells

  logical, dimension(:), allocatable :: ovlp3fn_cells !SVL array of dimension n_cells_pairs,
! contains .true. if there are non-zero three-center integrals belonging to the given cell,
! and .false. otherwise; filled up in integrate_ovlp3fn_p0.f90

  integer, dimension(:,:), allocatable :: basis_cell_to_all ! SVL array that for each basis function
                                                            ! from n_basis and a unit cell index gives 
                                                            ! the index of that basis function
                                                            ! among n_centers_basis_T functions; 
                                                            ! allocated in condense_basis_pairs_v2.f90

  integer, dimension(:,:), allocatable:: position_in_hamiltonian  ! i_cell_1, i_cell_2 -> i_cell_diff (_1 - _2)
  integer :: position_in_hamiltonian_dim1
  integer :: position_in_hamiltonian_dim2

  integer:: n_trash_start

  ! The following set of arrays are a set of auxiliary arrays for use only in subroutine initialize_centers_list .
  ! They are used to pass variables to subroutine get_include_centers, but in their original version
  ! are not intended for any other use.
  real*8, allocatable, private :: dist_sq_basis_integrals(:,:)
  real*8, allocatable, private :: dist_sq_ele_summation(:,:)
  real*8, allocatable, private :: dist_sq_hartree_potential(:,:)
  real*8, allocatable, private :: dist_sq_hartree_multipole(:,:)
  real*8, private              :: max_atom_diameter_sq

!=================shanghui add for DFPT_phonon====================================

  !------------------supercell basis function information--------
  integer, dimension (:,:),    allocatable :: cell_and_basis_to_basis_supercell 
  !-------todo : basis_diff_PBC(io,jo,juo) = io - (jo - juo)  = iuo 
  !              to replace position_in_hamiltonian


  !-----------------supercell center information---------------
  !                cell+atom ==> n_center 
  integer,  dimension(:,:),    allocatable:: cell_and_atom_to_center 
  ! n_centers_in_hamiltonian ==> n_center 
  integer, dimension(:), allocatable:: centers_in_hamiltonian        
  !                n_center  ==> n_centers_in_hamiltonian
  integer, dimension(:), allocatable:: inv_centers_in_hamiltonian  


  !-----------------supercell translation information---------------
  ! i_cell_1, i_cell_2 -> i_cell_diff (_1 - _2) in PBC  ---> 
  integer, dimension(:,:), allocatable:: position_in_hamiltonian_PBC  

  ! i_cell_1, i_cell_2 -> i_cell_add (_1 + _2) in PBC ---> 
  integer, dimension(:,:), allocatable:: position_in_hamiltonian_PBC_add  


  !-----------------begin sc_DFPT----------------------------------------
   integer,dimension(3) ::       min_nsc_DFPT
   integer,dimension(3) ::       max_nsc_DFPT
   integer,dimension(3) ::       nsc_DFPT
   integer:: n_cells_in_sc_DFPT
   integer, dimension(:,:), allocatable:: cell_index_sc_DFPT


   integer, dimension(:), allocatable:: center_in_sc_DFPT_to_center  ! n_centers_in_sc_DFPT ==> n_centers
   integer, dimension(:), allocatable:: center_to_center_in_sc_DFPT ! n_centers  ==> n_centers_in_sc_DFPT

   integer, dimension(:), allocatable:: center_in_sc_DFPT_to_cell_in_sc_DFPT
   integer, dimension(:), allocatable:: center_in_sc_DFPT_to_atom

  ! i_cell_1, i_cell_2 -> i_cell_diff ( 1 -  2) in  sc_DFPT under PBC  ---> 
   integer, dimension(:,:), allocatable:: cell_diff_sc_DFPT

  ! i_cell_1, i_cell_2 -> i_cell_add ( 1 +  2) in sc_DFPT under PBC ---> 
   integer, dimension(:,:), allocatable:: cell_add_sc_DFPT 

   integer, dimension(:,:), allocatable :: cell_and_atom_to_center_sc_DFPT
   integer, dimension (:,:), allocatable :: cell_and_basis_to_basis_sc_DFPT

   integer, dimension (:), allocatable :: cell_in_hamiltonian_to_cell_in_sc_DFPT
  !-----------------end sc_DFPT----------------------------------------

!=================shanghui end add for DFPT_phonon====================================


  ! -------

  ! Fortran 90 does not allow arrays to pointers, so we have to do it the following way

  type :: ptr_to_array
     integer, dimension(:), allocatable :: arr
     integer :: size
  end type ptr_to_array
  type(ptr_to_array), allocatable, private:: column_index(:,:)
  type(ptr_to_array), allocatable, private:: column_index_no_symmetry(:,:)

  integer :: max_matrix_size ! Maximum size of all local hamiltonians

! simple check to see if initialize_k_points was called previously with a different number
! of k points. Presently, only a minor safeguard for periodic Hartree-Fock
  integer:: prev_n_k_points

  !******
contains

  !-----------------------------------------------------------------------
  !****s* pbc_lists/initialize_bc_dependent_lists
  !  NAME
  !    initialize_bc_dependent_lists
  !  SYNOPSIS

  subroutine initialize_bc_dependent_lists 

    !  PURPOSE
    !    Initializes the periodic boundary condition lists.  
    !    Calls mainly other subroutines.
    !  USES
    use dimensions
    use runtime_choices
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use species_data
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer i
    character*100 :: info_str

    write(info_str,'(2X,A)')' '
    call localorb_info(info_str,use_unit,'(A)',OL_norm)    
    write(info_str,'(2X,A)') "Initializing index lists of integration centers etc. from given atomic structure:"
    call localorb_info(info_str,use_unit,'(A)',OL_norm)    

    atoms_in_supercell = atoms_in_structure !SVL needed for construction of auxiliary basis

    prev_n_k_points = 0

    if(n_periodic.eq.0)then
       n_k_points = 1
    end if

    call initialize_atom_radius_hartree()

    if(n_periodic > 0)then
       call evaluate_number_of_supercells_needed_in_periodic_systems()
    else
       n_cells = 1
       number_of_super_cells = 0  ! Only 0-cell
    end if

    call check_geometry_for_sanity

    call initialize_centers_data()

    ! dummy allocations performed even in the non-periodic case
    call initialize_k_points()

    call initialize_basis_lists

    call initialize_centers_list

    ! (Rundong) For fully-relativistic cases, we currently use PM_none for benchmark calculations.
    ! However, the following subroutine initialize_position_in_hamiltonian contains other functions
    ! that are needed, thus, we set packed_matrix_format /= PM_none here; but later, use PM_none 
    ! to pack the matrix in integrate_hamiltonian_matrix_p2.f90.
    if (packed_matrix_format /= PM_none .and. n_periodic > 0) then

       if(use_DFPT_phonon .or. use_friction) then 
       call initialize_sc_DFPT_phonon                 ! used in 1D ~ 3D
       call initialize_position_in_hamiltonian_PBC    ! used in 1D   
       endif   
 
       call initialize_position_in_hamiltonian
    end if

  end subroutine initialize_bc_dependent_lists
  !******
  !------------------------------------------------------------------------
  !****s* pbc_lists/allocate_center_data
  !  NAME
  !    allocate_center_data
  !  SYNOPSIS

  subroutine allocate_center_data 

    !  PURPOSE
    !    Allocates variables for centers information storage.
    use dimensions, only: n_centers
    use mpi_tasks, only: check_allocation
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: info
    character(*), parameter :: func = 'allocate_center_data'

    allocate(coords_center(3,n_centers),stat=info)
    call check_allocation(info, 'coords_center', func)

    allocate(species_center(n_centers),stat=info)
    call check_allocation(info, 'species_center', func)

    allocate(center_to_atom(n_centers),stat=info)
    call check_allocation(info, 'center_to_atom', func)

    allocate(center_to_cell(n_centers),stat=info)
    call check_allocation(info, 'center_to_cell', func)


! SVL Needed for periodic EXX
    allocate(center_to_cell_pairs(n_centers),stat=info)
    call check_allocation(info, 'center_to_cell_pairs', func)

  end subroutine allocate_center_data
  !******
  !------------------------------------------------------------------------
  !****s* pbc_lists/allocate_basis_lists
  !  NAME
  !    allocate_basis_lists
  !  SYNOPSIS

  subroutine allocate_basis_lists 

    !  PURPOSE
    !    Allocated variables for basis list information.
    use dimensions, only: n_centers_basis_T
    use mpi_tasks, only: check_allocation
    use runtime_choices, only: flag_rel, REL_4c_dks, REL_x2c
    use rel_x2c_mod, only: n_centers_basis_I_small
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE

    integer:: info
    character(*), parameter :: func = 'allocate_basis_lists'

    allocate( Cbasis_to_basis(n_centers_basis_T),stat=info)
    call check_allocation(info, 'Cbasis_to_basis', func)

    allocate( Cbasis_to_center (n_centers_basis_T),stat=info)
    call check_allocation(info, 'Cbasis_to_center', func)

    allocate( Cbasis_to_atom   (n_centers_basis_T),stat=info)
    call check_allocation(info, 'Cbasis_to_atom', func)

    if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       allocate( Cbasis_to_basis_s(n_centers_basis_I_small),stat=info)
       call check_allocation(info, 'Cbasis_to_basis_s', func)

       allocate( Cbasis_to_center_s (n_centers_basis_I_small),stat=info)
       call check_allocation(info, 'Cbasis_to_center_s', func)

       allocate( Cbasis_to_atom_s   (n_centers_basis_I_small),stat=info)
       call check_allocation(info, 'Cbasis_to_atom_s', func)
    endif

  end subroutine allocate_basis_lists
  !******
  !------------------------------------------------------------------------
  !****s* pbc_lists/allocate_center_lists
  !  NAME
  !     allocate_center_lists
  !  SYNOPSIS

  subroutine allocate_center_lists 

    !  PURPOSE
    !    Allocates variables for centers list.
    use dimensions
    use mpi_tasks, only: check_allocation
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer :: info
    character(*), parameter :: func = 'allocate_centers_lists'

    allocate(centers_hartree_potential(n_centers_hartree_potential),stat=info)
    call check_allocation(info, 'centers_hartree_potential', func, &
                          n_centers_hartree_potential)

    allocate(centers_hartree_multipole(n_centers_hartree_multipole),stat=info)
    call check_allocation(info, 'centers_hartree_multipole', func, &
          n_centers_hartree_multipole)

    allocate(centers_ele_summation( n_centers_ele_summation ),stat=info)
    call check_allocation(info, 'centers_ele_summation', func, &
          n_centers_ele_summation)

    allocate(centers_basis_integrals( n_centers_basis_integrals ),stat=info)
    call check_allocation(info, 'centers_basis_integrals', func, &
         n_centers_basis_integrals )

    allocate(new_centers_basis_integrals( n_centers_basis_integrals ),stat=info)
    call check_allocation(info, 'new_centers_basis_integrals', func, &
         n_centers_basis_integrals )

    allocate(inv_centers_basis_integrals( n_centers ),stat=info)
    call check_allocation(info, 'inv_centers_basis_integrals', func, &
          n_centers )


!---------shanghui add for DFPT_phonon------------------ 
    allocate(centers_in_hamiltonian( n_centers_in_hamiltonian ),stat=info)
    call check_allocation(info, 'centers_in_hamiltonian', func, &
          n_centers_in_hamiltonian )
    allocate(inv_centers_in_hamiltonian( n_centers ),stat=info)
    call check_allocation(info, 'inv_centers_in_hamiltonian', func, &
          n_centers )

!---------shanghui end add for DFPT_phonon------------------ 



  end subroutine allocate_center_lists
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/initialize_atom_radius_hartree
  !  NAME
  !    initialize_atom_radius_hartree
  !  SYNOPSIS

  subroutine initialize_atom_radius_hartree

    !  PURPOSE
    !    Calculates the outer radius of atoms and the outer radius of the real
    !    part of Hartree potential
    use dimensions
    use runtime_choices
    use Hartree_F_p_functions, only: F_erfc_single
    use mpi_tasks, only: check_allocation
    use basis, only: atom_radius
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8:: hartree_radius, pot
    integer :: info
    character(*), parameter :: func = 'initialize_atom_radius_hartree'

    ! This has already been done well before the current stage in prepare_scf:
    ! call  initialize_F_p_functions

    if (allocated(atom_radius_hartree)) then
       deallocate(atom_radius_hartree)
    endif
    allocate(atom_radius_hartree(n_species),stat=info)
    call check_allocation(info, 'atom_radius_hartree', func)

    if (allocated(atom_radius_hartree_sq)) then
       deallocate(atom_radius_hartree_sq)
    end if
    allocate(atom_radius_hartree_sq(n_species),stat=info)
    call check_allocation(info, 'atom_radius_hartree_sq', func)

    if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then

       hartree_radius = 12
       pot = 1e9
       do while (abs(pot) > far_distance_hartree_radius_threshold)
          hartree_radius = hartree_radius + 0.1
          pot =  F_erfc_single(hartree_radius, 0) !* 1000
       end do
       atom_radius_hartree = max(atom_radius, hartree_radius)
       atom_radius_hartree_sq = atom_radius_hartree**2

    else  ! cluster case without Ewald decomposition

       ! unscreened Coulomb tails => essentially infinite
       atom_radius_hartree =  1e10
       atom_radius_hartree_sq =  1e10

    end if

  end subroutine initialize_atom_radius_hartree
  !******
  !--------------------------------------------------------------------------------------------------------------------------  
  !****s* pbc_lists/evaluate_number_of_supercells_needed_in_periodic_systems
  !  NAME
  !     evaluate_number_of_supercells_needed_in_periodic_systems
  !  SYNOPSIS

  subroutine evaluate_number_of_supercells_needed_in_periodic_systems

    !  PURPOSE
    ! How many supercell in every direction is maximally needed in periodic systems ?
    ! Before going here, do calculate the atom_radius_hartree, we need
    ! to know how far a way from atom center interesting regions locate.    
    !  USES

    use constants
    use dimensions
    use runtime_choices
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: check_allocation
    use species_data
    use basis, only: atom_radius
    use geometry, only: periodic_unit_cell_translations, coords, frac_coords, &
        lattice_vector, map_to_center_cell_matrix, species, cart2frac
    implicit none

    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8,dimension(3):: coord_current
    real*8,dimension(3):: coefficients, coefficients_min, coefficients_max

    real*8:: radius_of_atom
    real*8:: angle_1
    real*8:: z_mid

    integer:: info

    integer:: i_atom, i_latt,  i_latt_2
    integer:: i_angular_1,i_angular_2
    character*100 :: info_str
    real*8,allocatable, dimension(:,:):: coords_temp
    integer:: n_coords_temp

    integer, parameter :: nsteps = 400
    !         real*8, dimension(2*nsteps) :: sintab, costab
    real*8, dimension(:),allocatable :: sintab, costab
    character(*), parameter :: func = 'evaluate_number_of_supercells_needed_in_periodic_systems'

    allocate(sintab(2*nsteps), stat=i_atom)
    call check_allocation(i_atom, 'sintab', func)
    allocate(costab(2*nsteps), stat=i_atom)
    call check_allocation(i_atom, 'costab', func)

    ! --- Map to center cell

    write(info_str,'(2X,A)') 'Mapping all atomic coordinates to central unit cell.'
    call localorb_info(info_str)
    
    ! VA: need to map periodic unit cell translations in fractional coordinates
    !     in particular when the lattice changes.
    !     if the lattice stays fix, one could also use cartesian coords right away
    call cart2frac(lattice_vector, coords,  frac_coords)
    periodic_unit_cell_translations(:,:) = periodic_unit_cell_translations(:,:) + frac_coords(:,:)

    if( transport_lead_calculation .or. transport_calculation .or. &
       flag_out_dipmat .or. out_esp .or. out_esp_full)then

       ! Shift coords to get the median point to the origin.
       do i_latt = 1,3
          z_mid = 0.5d0 * (minval(coords(i_latt,:)) + maxval(coords(i_latt,:)))
          do i_atom = 1, n_atoms
             coords(i_latt, i_atom) = coords(i_latt, i_atom) - z_mid
          end do
       end do

       ! Make sure that all atoms are within the Wigner-Seitz cell, now.
       do i_atom = 1,n_atoms
          coord_current= coords(:,i_atom)

          call map_to_center_cell(coords(1:3,i_atom))

          if( abs(coord_current(1) -coords(1,i_atom))> 0.1 .or. &
          abs(coord_current(2) -coords(2,i_atom))> 0.1 .or. &
          abs(coord_current(3) -coords(3,i_atom))> 0.1 )then
             write(use_unit,*)' Atom coordinates are not around origin!'
             write(use_unit,*)' old ', coord_current*bohr
             write(use_unit,*) 'new ',  coords(:,i_atom)*bohr
             stop
          end if

       end do
    end if

    !CC: Take care in case of TDI:
    if(use_thermodynamic_integration) then
      do i_atom = 1,n_atoms
         call TDI_map_to_center_cell(coords(1:3,i_atom),i_atom)
      end do
    else
      do i_atom = 1,n_atoms
         call map_to_center_cell(coords(1:3,i_atom))
      end do
    end if

    call cart2frac(lattice_vector, coords,  frac_coords)
    periodic_unit_cell_translations(:,:) = periodic_unit_cell_translations(:,:) - frac_coords(:,:)


    ! --- Prepare arrays

    ! coords_temp is coords + corners of unit cell
    n_coords_temp = n_atoms + 6
    allocate(coords_temp(3,n_coords_temp), stat=info)
    call check_allocation(info, 'coords_temp', func)
    coords_temp(1:3, 1:n_atoms) = coords(1:3,1:n_atoms)
    i_atom = n_atoms
    do i_latt = 1,3
       do i_latt_2 = 1,3
          if(i_latt /= i_latt_2)then
             i_atom = i_atom + 1
             coords_temp(1:3, i_atom) = lattice_vector(1:3,i_latt)*0.5d0 + lattice_vector(1:3,i_latt_2)*0.5d0
          end if
       end do
    end do


    ! Tabulate sin and cos over a full circle
    do i_angular_1 = 1,2*nsteps
       angle_1 = i_angular_1*pi/nsteps
       sintab(i_angular_1) = sin(angle_1)
       costab(i_angular_1) = cos(angle_1)
    enddo

    ! --- Figure out actual number of supercells.

    ! Take into account every unit cell where some atom might have a nonzero
    ! (erfc-screened) Coulomb interaction with an atom in the 0-cell.

    ! JW: I guess that this is "only Coulomb" because this has the farthest
    !     range.

    ! Go through a surface of sphere. The radius of sphere is the farthest distance
    ! we are interested in periodic systems.
    ! Calculate how many supercells in every directions we need in order to go there from center.
    ! Save maximum values. Do this for every atom.

    coefficients_max = 0.0
    coefficients_min = 0.0
    do i_atom = 1,n_coords_temp

       ! The farthest distance we are interested in periodic systems
       if(i_atom <= n_atoms)then
          radius_of_atom =  atom_radius(species(i_atom)) + &
               max(maxval(atom_radius_hartree),maxval(multipole_radius_free))  &
               + extra_adding_to_hartree_potential_distance
       else
          radius_of_atom =  maxval(atom_radius) + &
               max(maxval(atom_radius_hartree),maxval(multipole_radius_free)) &
               + extra_adding_to_hartree_potential_distance
       end if

       do i_angular_1 = 1,2*nsteps
          do i_angular_2 = 1,nsteps

             ! The coordinated of the sphere is discretized by the above steps. There
             ! is no good reason to choose 800 or 400 or something else. I thought that would be
             ! enough. -Paula
             coord_current(1) = coords_temp(1, i_atom) + radius_of_atom *sintab(i_angular_1)* costab(i_angular_2)
             coord_current(2) = coords_temp(2, i_atom) + radius_of_atom *sintab(i_angular_1)* sintab(i_angular_2)
             coord_current(3) = coords_temp(3, i_atom) + radius_of_atom *costab(i_angular_1)

             coefficients = matmul(map_to_center_cell_matrix, coord_current)

             ! Save largers and smallest values. Note: smallest are negative.
             do i_latt = 1,n_periodic
                coefficients_max(i_latt) =  max(coefficients_max(i_latt), coefficients(i_latt))
                coefficients_min(i_latt) =  min(coefficients_min(i_latt), coefficients(i_latt))
             end do
          end do
       end do
    end do

    coefficients(:) = max( abs(coefficients_max(:)), abs(coefficients_min(:)))

    number_of_super_cells(:) = ceiling(coefficients(:))

    ! Directions which are not periodic do not need supercells.
    do i_latt = n_periodic+1,3
       number_of_super_cells(i_latt) = 0
    end do

    n_cells = product(2*number_of_super_cells+1)


    deallocate(coords_temp)
    deallocate(sintab)
    deallocate(costab)

  end subroutine evaluate_number_of_supercells_needed_in_periodic_systems
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/check_geometry_for_sanity
  !  NAME
  !    check_geometry_for_sanity
  !  SYNOPSIS

  subroutine check_geometry_for_sanity

    !  PURPOSE
    !
    !    Checks sanity of geometry (atoms too close to other atoms or vacuum).
    !    Please note that both the Bravais lattice
    !    (initialize_bravais_quantities) and the outer_radii must already be
    !    initialized.
    !
    !  USES
    use dimensions
    use runtime_choices
    use constants, only: bohr
    use localorb_io, only: localorb_info
    use mpi_tasks, only: aims_stop
    use timing, only: warn_surface_dipole
    use basis, only: atom_radius
    use geometry, only: coords, lattice_vector, species
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_atom_1
    real*8:: dis

    real*8 :: cell_length, shifted_vacuum_level
    integer :: mul_int
    integer :: i_shift

    character*100 :: info_str
    character(*), parameter :: func = 'check_geometry_for_sanity'

    call check_for_close_encounters(.false.)

    if (.not. allocated(atom_radius)) then
       call aims_stop('atom_radius not allocated', func)
    end if

    ! Also, check validity of z distance of atoms from vacuum level, if 
    ! a surface work function or dipole calculation was requested.
    ! Sufficient to check z coordinates
    if (calculate_work_function .or. use_dipole_correction) then
      ! first, shift the requested vacuum level into the first cell
      cell_length =  abs((maxval(lattice_vector(3,:)) - minval(lattice_vector(3,:))))
      mul_int = int(floor((vacuum_z_level)/cell_length))
      shifted_vacuum_level = vacuum_z_level - mul_int * cell_length    

      do i_atom_1 = 1, n_atoms, 1

        dis = 10000.
        do i_shift = -1, 1, 1
          dis = min( dis, abs( coords(3,i_atom_1)-(i_shift*cell_length+shifted_vacuum_level) ) )
        enddo

        if (dis .lt. 6.d0/bohr) then
          ! Hard stop. A work function calculation or dipole correction less that 6 Angstrom from
          ! the nearest atom will not work.

             write(info_str,'(1X,A,I8,A)') & 
               '* ERROR: atom ',i_atom_1,' is closer than 6 Angstroms to the requested vacuum level'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* for the surface dipole correction or work function calculation.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* A surface dipole determined close to the atoms of a slab will not be accurate.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* As a result, any numerical values that are affected by the dipole (including'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* total energies) may be subject to noticeable inaccuracies.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* FHI-aims allows one to use rather large vacuum regions at very low cost - '
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* easily 40 AA or 80 AA if keyword "hartree_convergence_parameter" is set correctly.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* To ensure that your calculated surface dipole is accurate, '
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* please use a larger vacuum spacing in your slab calculation.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* Stopping the calculation as a safety precaution.'
             call localorb_info(info_str)
             stop
          
        else if ( dis.lt.atom_radius(species(i_atom_1)) ) then
          ! Give a loud warning but let the code run. 

             write(info_str,'(1X,A,I8,A)') & 
               '* Warning: Atom ', i_atom_1, ' is close to the z vacuum level requested for the'
             call localorb_info(info_str)
             write(info_str,'(1X,A,I8,A)') & 
               '* surface dipole correction or work function calculation.'
             call localorb_info(info_str)
             write(info_str,'(1X,A,F12.8,A)') & 
               '* Distance: ', dis*bohr, ' Angstroms.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* A surface dipole determined close to the atoms of a slab will not be accurate.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* As a result, any numerical values that are affected by the dipole (including'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* total energies) may be subject to noticeable inaccuracies.'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* The calculation will be allowed to continue, but we strongly recommend'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* to repeat it with a larger vacuum layer between slabs, large enough to make'
             call localorb_info(info_str)
             write(info_str,'(1X,A)') & 
               '* this warning disappear.'
             call localorb_info(info_str)

             warn_surface_dipole = .true.

        end if

      enddo

    end if

  end subroutine check_geometry_for_sanity
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/check_for_close_encounters
  !  NAME
  !    check_for_close_encounters
  !  SYNOPSIS

  subroutine check_for_close_encounters(output_mindist)

    !  PURPOSE
    !
    !    Checks that minimum distance of atoms is larger than 0.4 bohr.
    !    Please note that for this two work, the Bravais lattice must already
    !    be initialized via initialize_bravais_quantities().
    !
    !  USES

    use bravais
    use constants, only: bohr
    use dimensions
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use mpi_tasks, only: aims_stop
    use geometry, only: coords, lattice_vector, min_grid_centre_dist
    ! AJL Sept2017: No longer needed here. Moved down to min_atomic_dist
!    use species_data, only: species_pseudoized, no_basis
    implicit none

    !  ARGUMENTS

    logical, intent(IN) :: output_mindist

    !  INPUTS
    !    o output_mindist -- Output minimal distance?
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !    Edited Sept 2017 by AJL so that distances beteen physically meaningful centres (e.g. atoms)
    !    can be differentiates from non-physical centres (e.g. empty sites), and the checks made on
    !    distance made more appropriate. Does not consider distances between multipoles and grid centres.
    !  SOURCE

    real*8 :: dis
    integer :: i_grid_1, i_grid_2
    integer :: cell_grid_1, cell_grid_2, cell_grid_3
    character*150 :: info_str
    character(*), parameter :: func = 'check_for_close_encounters'

    ! AJL, Sept2017
    real*8 :: dis_physical
    integer :: i_physical_1, i_physical_2
    integer :: cell_physical_1, cell_physical_2, cell_physical_3

    ! AJL, Sept2017. Adapted min_atomic_dist so it does the checks for both all grid centres (as previously)
    ! and also for just physically meaningful centres. Thus, sanity checks can be made a bit more focused.
    !call min_atomic_dist(lattice_vector, coords, dis, i_atom_1, i_atom_2, cell_1, cell_2, cell_3)
    call min_grid_centre_dist(lattice_vector, coords, dis, i_grid_1, i_grid_2, cell_grid_1, cell_grid_2, cell_grid_3, &
                         dis_physical, i_physical_1, i_physical_2, cell_physical_1, cell_physical_2, cell_physical_3)

    if (output_mindist) then
       if ((n_atoms.eq.1).and.(n_periodic.eq.0)) then
          write(info_str, '(2X,A)') '| Isolated single atom, no quantum-mechanical atoms close by.' 
          call localorb_info(info_str, use_unit, '(A)', OL_norm)
       else
          ! AJL, Sept2017: Should this now be updated to distance between just physical centres?
          ! I presume so and so have taken the libery to make the changes appropriately.
          write(info_str, "(2X,'| The smallest distance between any two atoms is ',F18.8,' AA.')")&
          & dis_physical*bohr
          call localorb_info(info_str, use_unit, '(A)', OL_norm)
          write(info_str, "(2X,'| The first atom of this pair is atom number     ',I18,' .')")&
          & i_physical_1
!          & i_atom_1
          call localorb_info(info_str, use_unit, '(A)', OL_norm)
          write(info_str, "(2X,'| The second atom of this pair is atom number    ',I18,' .')")&
          & i_physical_2
!          & i_atom_2
          call localorb_info(info_str, use_unit, '(A)', OL_norm)
          if (n_periodic.eq.3) then
            write(info_str, "(2X,'| Wigner-Seitz cell of the first atom image      ',1X,I5,1X,I5,1X,I5,' .')") &
            & cell_physical_1, cell_physical_2, cell_physical_3
!            & cell_1, cell_2, cell_3
            call localorb_info(info_str, use_unit, '(A)', OL_norm)
            write(info_str, "(2X,'| (The Wigner-Seitz cell of the second atom is ',I1,1X,I1,1X,I1,1X,A)")&
            & 0, 0, 0, ' by definition.) '
            call localorb_info(info_str, use_unit, '(A)', OL_norm)
          end if
       end if
    end if

    ! AJL, Sept2017
    ! Now check sanity for all grid centres, i.e. including empty sites. Threshold arbitrarily set to
    ! 0.001 but no idea if this is reasonable for all cases: it needs to be quite small because sometimes
    ! in QM/MM we have empty sites on point charges very close to, but not on top of, pseudocores. In fact,
    ! the reason I've implemented these different categories of neighbour distances is because we get bad errors
    ! when the grid centres are absolutely on top of each other, i.e. dist = 0.0. So in theory I would be happy making
    ! this number much smaller, but definitely not bigger. Perhaps an input variable could be used to define the distance
    ! in these cases?
    if (dis < 0.001d0) then
          ! We are now looking at all grid centres, so no need for cumbersome if statement to differentiate centres.
          call localorb_info(' ')
          write(info_str,'(1X,A,I8,A,I8,A)') &
          '* ERROR: Grid centres ',i_grid_1,' and ',i_grid_2,' are too close to each other!'
          call localorb_info(info_str)
          write(info_str,'(1X,A,F12.8,A)') &
          & '* Distance: ', dis*bohr, &
          & ' Angstroms sounds clearly unphysical. Please check your input geometry or, if you think '
          call localorb_info(info_str)
          write(info_str,'(A)') &
          & ' that you need grid centres this close together, contact the Development team to discuss changing this error message. '
          call localorb_info(info_str)
          write(info_str,'(1X,A,I8,A,I8,A)') &
          '* ERROR: Grid centres ',i_grid_1,' and ',i_grid_2,' too close to each other!'
          call aims_stop(info_str, func)
    end if

    ! AJL, Sept2017: Now check for sanity of atomic structure
    if (dis_physical < 0.4d0) then
! AJL, Sept2017: This is no longer needed as its been shifted into min_atomic_dist.
! Kept here as a reference in case of outstanding errors.
!       if ( (.not.(empty(i_atom_1).and.species_pseudoized(species(i_atom_2)))) .and. & 
!            (.not.(species_pseudoized(species(i_atom_1)).and.empty(i_atom_2))) .and. &
!            (.not.(no_basis(species(i_atom_1)).and.(no_basis(species(i_atom_2)))) )) then

          call localorb_info(' ')
          write(info_str,'(1X,A,I8,A,I8,A)') & 
          '* ERROR: atoms ',i_physical_1,' and ',i_physical_2,' are too close to each other, see above!'
          call localorb_info(info_str)
          write(info_str,'(1X,A,F12.8,A)') & 
          & '* Distance: ', dis_physical*bohr, &
          & ' Angstroms sounds clearly unphysical. Please check your input geometry.'
          call localorb_info(info_str)
          write(info_str,'(1X,A,I8,A,I8,A)') & 
          '* ERROR: atoms ',i_physical_1,' and ',i_physical_2,' too close to each other!'
          call aims_stop(info_str, func)
!       end if
    end if

    ! test
    ! call aims_stop ('Stop for testing purposes.', func)

  end subroutine check_for_close_encounters
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/check_multipole_distances
  !  NAME
  !    check_for_close_encounters
  !  SYNOPSIS

  subroutine check_multipole_distances()

    !  PURPOSE
    !
    !  In the presence of external (embedding) Coulomb charges or multipoles,
    !  we must also make sure that there
    !  are no Coulomb singularities without their own integration grids in the distance 
    !  range in which the quantum system's density is non-zero.
    !  Integration grids can be placed for instance by using empty sites.
    !
    !  USES

    use constants, only: bohr
    use bravais
    use dimensions
    use localorb_io, only: localorb_info
    use basis, only: atom_radius_sq
    use geometry, only: coords, empty_coords, lattice_vector, multipole_coords, &
        species, min_multipole_dist
    implicit none

    !  ARGUMENTS

    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2012).
    !  SOURCE

    real*8 :: dis
    integer :: i_atom, i_multipole
    integer, dimension(n_multipoles) :: offending_multipole
    character*150 :: info_str
    character(*), parameter :: func = 'check_multipole_distances'

    call min_multipole_dist( lattice_vector, coords, species, atom_radius_sq, multipole_coords, & 
    &                        empty_coords, dis, i_atom, i_multipole, offending_multipole)

    ! offending_multipole is initialised 0, so check if it has been changed
    if (offending_multipole(1).ne.0) then
       write(info_str,'(1X,A,I8,A,I12,A)') & 
       '* ERROR: Atom ',i_atom,' and multipole number ',i_multipole,' are too close to each other!'
       call localorb_info(info_str)
       write(info_str,'(1X,A,F12.8,A)') & 
       & '* Distance: ', dis*bohr, '.'
       call localorb_info(info_str)
       write(info_str,'(1X,A,A)') & 
         '* Multipoles inside or close to the basis function radius of an atom should',&
       & ' have integration grids of their own.'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') & 
         '* Grids can be placed on a multipole by using empty sites without basis functions.'
       call localorb_info(info_str)

       ! write out all offending multipoles
       write(info_str,'(1X,A)') &
       '* The following multipole sites require empty sites:'
       call localorb_info(info_str)

       do i_multipole = 1, n_multipoles
          if (offending_multipole(i_multipole).ne.0) then
             write(info_str,'(1X,I8)') offending_multipole(i_multipole)
             call localorb_info(info_str)
          end if
       end do

       ! This *should* be a critical error.
       ! However, if we implement this check then currently there is no
       ! safety net to check implementations changes which might effect the
       ! the bq implementations in the NaCl_qmm regression test.
       
       ! Really what this needs is separate regression tests for (a) bqs and (b) psps
       ! To be implemented
       ! call aims_stop('Multipoles without integration grids found.', func)
    end if

  end subroutine check_multipole_distances
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/initialize_centers_data
  !  NAME
  !    initialize_centers_data
  !  SYNOPSIS

  subroutine initialize_centers_data

    !  PURPOSE
    ! Initializes:
    ! o  species_center
    ! o  coords_center
    ! o  center_to_atom
    ! o  center_to_cell
    ! These tables are created for periodic boundary conditions, but in order to use same 
    ! subroutines also in clusters case, we have to allocate and initialize them.
    use dimensions
    use runtime_choices
    use mpi_tasks, only: check_allocation
    use geometry, only: coords, species, lattice_vector
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE


    integer:: i_cell_1, i_cell_2, i_cell_3, i_cell_n
    integer:: i_atom,   i_center
    integer:: i_coord
    integer :: info
    character(*), parameter :: func = 'initialize_centers_data'

    if(n_periodic > 0)then

       ! Estimate to the number of basis includet to the periodic system
       ! This estimate gives far too large values!

       n_centers = n_atoms*n_cells

       !            write(info_str,'(2X,A,I10)') '| Number of centers                              :', n_centers
       !            call localorb_info(info_str,use_unit,'(A)')    
    else

       ! In nonperiodic cluster case, centers are always same than atoms.
       n_centers = n_atoms

    end if

    ! Copy to get sum over BvK cells and decouple from number of k-points
    n_k_points_xyz_nosym = n_k_points_xyz
       
    ! can now allocate actual data for each center 
    call allocate_center_data

    ! These are really for periodic systems but we need to dummy-initialize them also for clusters,
    ! because we want to use same routines to periodic and nonperiodic systems.

    ! In nonperiodic cluster case, centers are always same than atoms.
    ! and coords_center == coords,  species_center == species
    do i_atom = 1, n_atoms
       coords_center(1:3,i_atom) = coords(1:3,i_atom)
    end do

    do i_atom = 1, n_atoms
       species_center(i_atom) = species(i_atom)
    end do

    ! In the center cell and clusters center_to_atom is pointing to itself and
    ! cell index is zero.

    i_center  = 0
    i_cell_n = 1
    do i_atom = 1,n_atoms

       i_center = i_center + 1
       center_to_atom(i_center) = i_atom
       center_to_cell(i_center)  =  i_cell_n
       center_to_cell_pairs(i_center) = i_cell_n

    end do

! SVL Basis pair related variables
    n_cells_pairs = (2*number_of_super_cells(1)+1)*(2*number_of_super_cells(2)+1)*&
         (2*number_of_super_cells(3)+1)
    n_cells_bvk = n_k_points_xyz_nosym(1)*n_k_points_xyz_nosym(2)*n_k_points_xyz_nosym(3)
    if (.not.allocated(cell_index_bvk)) then
       allocate(cell_index_bvk(n_cells_bvk,3), stat=info)
       call check_allocation(info, 'cell_index_bvk        ')
    endif
    cell_index_bvk = 0
    if (.not.allocated(cell_index_pairs)) then
       allocate(cell_index_pairs(n_cells_pairs,3), stat=info)
       call check_allocation(info, 'cell_index_pairs', func)
    endif
    cell_index_pairs = 0
    if (.not.allocated(ovlp3fn_cells)) then
       allocate(ovlp3fn_cells(n_cells_pairs), stat=info)
       call check_allocation(info, 'ovlp3fn_cells', func)
    endif
    ovlp3fn_cells = .false.
    
    ! In periodic case we have some more information in these tables.
    if(n_periodic .gt. 0)then

       ! These have to always done in this order, so that basis function initialozation
       ! goes correctly.

       do i_cell_1 = -number_of_super_cells(1),  number_of_super_cells(1)
          do i_cell_2 = -number_of_super_cells(2),  number_of_super_cells(2)
             do i_cell_3 = -number_of_super_cells(3),  number_of_super_cells(3) 

                if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then

                   i_cell_n = i_cell_n + 1

                   do i_atom = 1,n_atoms

                      i_center = i_center + 1

                      do i_coord = 1,3
                         coords_center(i_coord,i_center) =  coords(i_coord, i_atom) &
                              + lattice_vector(i_coord,1)*i_cell_1 &
                              + lattice_vector(i_coord,2)*i_cell_2 &
                              + lattice_vector(i_coord,3)*i_cell_3


                      end do


                      center_to_atom(i_center) = i_atom

                      center_to_cell(i_center)  =  i_cell_n


 
! SVL need for periodic EXX
                      center_to_cell_pairs(i_center) = i_cell_n

                      species_center(i_center) = species(i_atom)

! SVL information for pairs of basis functions
                      cell_index_pairs(i_cell_n,1) = i_cell_1
                      cell_index_pairs(i_cell_n,2) = i_cell_2
                      cell_index_pairs(i_cell_n,3) = i_cell_3

                   end do
                end if
             end do
          end do
       end do
       i_cell_n = 1
       do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
          do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
             do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2

                if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then

                   i_cell_n = i_cell_n + 1
                   cell_index_bvk(i_cell_n,1) = i_cell_1
                   cell_index_bvk(i_cell_n,2) = i_cell_2
                   cell_index_bvk(i_cell_n,3) = i_cell_3
                endif
             enddo
          enddo
       enddo

    end if


  end subroutine initialize_centers_data
  !----------------------------------------------------------------------------
  !****s* pbc_lists/band_k_frac
  !  NAME
  !    band_k_frac
  !  SYNOPSIS
  
  ! These function calculates the fractional coordinate of the k-point
  ! used in the band structure calculation.
  
   pure function band_k_frac(my_k, i_band)
   ! Calculate the reciprocal lattice vector
   ! in fractional coordinates

      use runtime_choices, only: band_begin, band_end, n_points_in_band
      implicit none
       
      integer, intent(in) :: my_k
      integer, intent(in) :: i_band
      real(kind=8), dimension(3)  :: band_k_frac

      ! Calculate reciprocal coordinates of current k-point
      band_k_frac(:) =  band_begin(i_band,:) +  real(my_k-1)/real(n_points_in_band(i_band)-1) &
                     *( band_end(i_band,:) -  band_begin(i_band,:))

   end function band_k_frac

  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/initialize_k_points
  !  NAME
  !   initialize_k_points
  !  SYNOPSIS

  subroutine initialize_k_points()

    !  PURPOSE
    !    Initializes the k-point list: 
    !    o n_k_points
    !    o k_phase
    !    o  k_weights. 
    !    Either reads lists from file k_list.in or constructs them from the
    !    control.in settings.
    !  USES

    use dimensions
    use runtime_choices
    use constants, only: pi
    use localorb_io, only: localorb_info, localorb_allinfo, use_unit, OL_low, &
        OL_norm
    use mpi_tasks
    use spglib_symmetry, only: k_point_symmetry_check_spg
    use timing, only : warn_idle_cpus
    implicit none

    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer:: i_cell_1, i_cell_2, i_cell_3, i_cell_n
    integer:: i_k_point
    integer:: i_x, i_y, i_z,info
    integer:: i_normal ! YY for ik2irred_map
    integer, allocatable, dimension(:) :: normal_to_stable ! YY for ik2irred_map
    !! YY: sorry if i_normal, normal_to_stable feels odd.
    !!     but they are related to the stabiliser in symmetry_reduce_k_points.f90
    !!     just as "JW : Why not?" in that file. I am also confused about why 
    !!     we need it. However, I don't what to modify it in case it is important.
    !!     In order to make the ik2irred_map the correct mapping I have to remapping
    !!     it from normal irred kpt to stabiliser irred kpt.
    !!     I use them for the tetrahedron method.
    !! YY
    real*8:: r_x, r_y, r_z
    character*100 :: info_str
    integer, allocatable,dimension(:,:,:):: k_number
    character(*), parameter :: func = 'initialize_k_points'

    if(n_periodic.eq.0)then
       ! Non-periodic systems - only dummy allocations

       n_k_points = 1
       allocate(k_weights(n_k_points),stat=info)
       call check_allocation(info, 'k_weights', func)
       allocate(k_phase_base(3, n_k_points), stat=info)
       call check_allocation(info, 'k_phase_base', func)
       allocate(k_phase(1,n_k_points),stat=info)
       call check_allocation(info, 'k_phase', func)

! SVL k_phase for exact exchange
       if(use_periodic_hf .or. (use_full_spectrum .and. n_periodic .gt. 0) )then
          allocate(k_phase_exx(1,n_k_points),stat=info)
          call check_allocation(info, 'k_phase_exx', func)
       endif

       allocate(k_point_list(n_k_points,3),stat=info)
       call check_allocation(info, 'k_point_list', func)
       k_point_list(:,:) = 0.d0
       allocate(k_minus_q_point_list(n_k_points,3),stat=info)
       call check_allocation(info, 'k_minus_q_point_list', func)
       k_minus_q_point_list(:,:) = 0.d0

       k_weights = 1
       k_phase_base = (1.d0, 0.d0)
       k_phase = 1
       if(use_periodic_hf .or. (use_full_spectrum .and. n_periodic .gt. 0) ) then
           k_phase_exx = (1.d0,0.d0)
       endif

       n_k_points_group = 1

    else

       write(info_str,'(2X,A,I2,A,I10)') ' '
       call localorb_info(info_str,use_unit,'(A)') 

       write(info_str,'(2X,A,I2,A,I10)')'Initializing the k-points'
       call localorb_info(info_str,use_unit,'(A)') 

       
       if(read_k_points_from_file)then


          write(info_str,'(2X,A)')  'External k-points from file k_list.in:'
          call localorb_info(info_str,use_unit,'(A)') 

          open(88,file='k_list.in')

          read(88,*)  n_k_points_xyz(1), n_k_points_xyz(2), n_k_points_xyz(3)
          n_k_points_xyz_nosym = n_k_points_xyz
          write(info_str,'(2X,A,3I10)')'| k-points in every direction (for information only):', &
             n_k_points_xyz(1), n_k_points_xyz(2), n_k_points_xyz(3)
          call localorb_info(info_str,use_unit,'(A)') 

          read(88,*)  n_k_points

          write(info_str,'(2X,A,I10)') '| Number of k-points:', n_k_points
          call localorb_info(info_str,use_unit,'(A)') 
          
          n_k_points_nosym = n_k_points
          
          allocate(k_weights(n_k_points),stat=info)
          call check_allocation(info, 'k_weights', func)

          if(packed_matrix_format == PM_none) then
             ! The full k_phase has to be calculated for PM_none since initialize_position_in_hamiltonian
             ! is not called in this case.
             write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells*n_k_points*16)/2**10, ' KiB for k_phase.'
             call localorb_info(info_str,use_unit,'(A)')
             allocate(k_phase( n_cells, n_k_points),stat=info)
             call check_allocation(info, 'k_phase', func)
          endif



          allocate(k_phase_base(3, n_k_points),stat=info)
          call check_allocation(info, 'k_phase_base', func)

! SVL k_phase for exact exchange
          if(use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt. 0))then
                write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells_bvk*n_k_points*16)/2**10, ' KiB for k_phase_exx.'
                call localorb_info(info_str,use_unit,'(A)')
                allocate(k_phase_exx( n_cells_bvk, n_k_points),stat=info)
                call check_allocation(info, 'k_phase_exx', func)
          endif

          allocate(k_point_list(n_k_points,3),stat=info) ! only for output purposes!
          call check_allocation(info, 'k_point_list', func)
          k_point_list(:,:) = 0.d0
          allocate(k_minus_q_point_list(n_k_points,3),stat=info)
          call check_allocation(info, 'k_minus_q_point_list', func)
          k_minus_q_point_list(:,:) = 0.d0
          allocate(irk_point_mapping(n_k_points),stat=info)
          call check_allocation(info, 'irk_point_mapping', func)

          ! This array is now kept past the final scf cycle. But, in a reinitialization
          ! (for a new geometry) it has to be renewed, too - so deallocate first.
          if (allocated(cell_index)) then
            deallocate (cell_index)
          end if

          allocate(cell_index(n_cells,3),stat=info)
          call check_allocation(info, 'cell_index', func)

          cell_index = 0



          i_k_point = 0
          if(packed_matrix_format == PM_none) k_phase = 1.0

          if((use_periodic_hf .and. use_hf_realspace) .or. (use_full_spectrum .and. n_periodic .gt. 0))&
               k_phase_exx = (1.d0,0.d0)

          do  i_k_point = 1,  n_k_points

             read(88,*)  r_x, r_y, r_z,  k_weights(i_k_point)

             k_point_list (i_k_point,1) = r_x
             k_point_list (i_k_point,2) = r_y
             k_point_list (i_k_point,3) = r_z

             k_phase_base(1, i_k_point) = exp((0.d0,1.d0) * 2*pi * r_x)
             k_phase_base(2, i_k_point) = exp((0.d0,1.d0) * 2*pi * r_y)
             k_phase_base(3, i_k_point) = exp((0.d0,1.d0) * 2*pi * r_z)

             if(r_x < 0)then
                r_x = 1 + r_x
             end if
             if(r_y < 0)then
                r_y = 1 + r_y
             end if
             if(r_z < 0)then
                r_z = 1 + r_z
             end if

             i_cell_n = 1

             do i_cell_1 = -number_of_super_cells(1),  number_of_super_cells(1)
                do i_cell_2 = -number_of_super_cells(2),  number_of_super_cells(2)
                   do i_cell_3 = -number_of_super_cells(3),  number_of_super_cells(3) 

                      if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then

                         i_cell_n = i_cell_n + 1                         

                         if(packed_matrix_format == PM_none) then
                            k_phase(i_cell_n, i_k_point) &
                            & = k_phase_base(1, i_k_point)**i_cell_1 &
                            & * k_phase_base(2, i_k_point)**i_cell_2 &
                            & * k_phase_base(3, i_k_point)**i_cell_3
                         end if


                         cell_index(i_cell_n, 1) = i_cell_1
                         cell_index(i_cell_n, 2) = i_cell_2
                         cell_index(i_cell_n, 3) = i_cell_3
! SVL calculate the phase factors for the exact exchange
                         !if (use_hf_realspace) then
                         !   k_phase_exx( i_cell_n, i_k_point) &
                         !   & = k_phase_base(1, i_k_point)**i_cell_1 &
                         !   & * k_phase_base(2, i_k_point)**i_cell_2 &
                         !   & * k_phase_base(3, i_k_point)**i_cell_3
                         !end if
                      end if
                   end do  ! i_cell_3
                end do  ! i_cell_2
             end do  ! i_cell_1
             ! SVL calculate the phase factors for the exact exchange
             if(use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt. 0))then

                i_cell_n = 1

                do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
                   do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                      do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2

                         if (i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0) then
                            i_cell_n = i_cell_n + 1
                            k_phase_exx( i_cell_n, i_k_point) &
                            & = k_phase_base(1, i_k_point)**i_cell_1 &
                            & * k_phase_base(2, i_k_point)**i_cell_2 &
                            & * k_phase_base(3, i_k_point)**i_cell_3
                         endif
                      enddo
                   enddo
                enddo

             endif

          end do  ! i_k_point

          k_minus_q_point_list = k_point_list

          close(88)
          real_eigenvectors = .false.

          if(abs(sum(k_weights)-1) >1e-5)then
             write(use_unit,*) 'Error: sum of k-vector weights is not one!', sum(k_weights)
             stop
          end if

          ! There can be numerical inaccuracy. This takes care of them.
          k_weights = k_weights/  sum(k_weights)


       else ! Not reading k-points from the file.

          n_k_points = n_k_points_xyz(1) *  n_k_points_xyz(2) * n_k_points_xyz(3)
          allocate(k_number(n_k_points_xyz(1), n_k_points_xyz(2), n_k_points_xyz(3)),stat=info)
          k_number = 1
          call check_allocation(info, 'k_number', func)
          
          ! Copy n_k_points to get BvK related quantities
          n_k_points_nosym = n_k_points
          
          ! allocate mapping from nosym kpoints to reduced kpoints
          allocate(ik2irred_map(n_k_points_xyz(1), n_k_points_xyz(2), n_k_points_xyz(3)), stat=info)
          call check_allocation(info, 'ik2irred_map', func)
          i_k_point = 0
          do i_x = 1, n_k_points_xyz(1)
             do i_y = 1, n_k_points_xyz(2)
                do i_z = 1, n_k_points_xyz(3)
                   i_k_point = i_k_point + 1
                   ik2irred_map(i_x,i_y,i_z) = i_k_point
                enddo
             enddo
          enddo
          !
          
          if(use_symmetry_reduced_k_grid)then  
             ! Calculating the k-grid symmetries
             if (use_symmetry_reduced_spg)then
               call k_point_symmetry_check_spg(n_k_points_xyz, k_points_offset,&
               &                           k_number)
             else
               call k_point_symmetry_check(n_k_points_xyz, k_points_offset, &
               &                           k_number)
             endif
             n_k_points = count(k_number > 0)
          end if



          write(info_str,'(2X,A,I10)')'| Number of k-points                             :',n_k_points
          call localorb_info(info_str,use_unit,'(A)') 


          if(packed_matrix_format == PM_none) then
             ! The full k_phase has to be calculated for PM_none since initialize_position_in_hamiltonian
             ! is not called in this case.
             write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells*n_k_points*16)/2**10, ' KiB for k_phase.'
             call localorb_info(info_str,use_unit,'(A)')
             allocate(k_phase(n_cells,n_k_points),stat=info)
             call check_allocation(info, 'k_phase', func)
          endif



          allocate(k_phase_base(3, n_k_points),stat=info)
          call check_allocation(info, 'k_phase_base', func)

! SVL k_phase for the exact exchange
          if((use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt.0)).and..not.use_symmetry_reduced_spg )then
                write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells_bvk*n_k_points*16)/2**10, ' KiB for k_phase_exx.'
                call localorb_info(info_str,use_unit,'(A)')
                allocate(k_phase_exx(n_cells_bvk,n_k_points),stat=info)
                call check_allocation(info, 'k_phase_exx', func)
          endif

          allocate(k_point_list(n_k_points,3),stat=info) ! only for output purposes!
          call check_allocation(info, 'k_point_list', func)
          k_point_list(:,:) = 0.d0
          if(.not.use_symmetry_reduced_spg )then
             allocate(k_minus_q_point_list(n_k_points,3),stat=info)
             call check_allocation(info, 'k_minus_q_point_list', func)
             k_minus_q_point_list(:,:) = 0.d0
             allocate(irk_point_mapping(n_k_points),stat=info)
             call check_allocation(info, 'irk_point_mapping', func)
          endif
          allocate(k_weights(n_k_points),stat=info)
          call check_allocation(info, 'k_weights', func)

          allocate(cell_index(n_cells,3),stat=info)
          call check_allocation(info, 'cell_index', func)
          cell_index = 0

          i_k_point = 0
          if(packed_matrix_format == PM_none) k_phase = 1.0
          if(((use_periodic_hf .and. use_hf_realspace) .or. (use_full_spectrum .and. n_periodic .gt. 0).and..not.use_symmetry_reduced_spg ))&
               k_phase_exx = (1.d0,0.d0)


          if(n_periodic > 0 )then
            allocate(normal_to_stable(n_k_points),stat=info)
            call check_allocation(info, 'normal_to_stable', func)

             do i_x = 1, n_k_points_xyz(1)
                do i_y = 1, n_k_points_xyz(2)
                   do i_z = 1, n_k_points_xyz(3)

                      if(k_number(i_x, i_y, i_z)> 0)then
                         i_k_point = i_k_point + 1
                         ! YY for ik2irred_map
                         if (.not. use_symmetry_reduced_spg) then
                           i_normal = ik2irred_map(i_x, i_y, i_z)
                           normal_to_stable(i_normal) = i_k_point
                         end if
                         ! YY
                         i_cell_n = 1
                         k_weights(i_k_point) = dble(k_number(i_x,i_y,i_z))

                         !write(use_unit,*) i_x,n_k_points_xyz(1),r_x
                         r_x = dble(i_x-1) / dble(n_k_points_xyz(1)) + k_points_offset(1)
                         r_y = dble(i_y-1) / dble(n_k_points_xyz(2)) + k_points_offset(2)
                         r_z = dble(i_z-1) / dble(n_k_points_xyz(3)) + k_points_offset(3)
                         k_point_list(i_k_point,1) = r_x
                         k_point_list(i_k_point,2) = r_y
                         k_point_list(i_k_point,3) = r_z

                         k_phase_base(1, i_k_point) = exp((0.d0,1.d0) * 2*pi * r_x)
                         k_phase_base(2, i_k_point) = exp((0.d0,1.d0) * 2*pi * r_y)
                         k_phase_base(3, i_k_point) = exp((0.d0,1.d0) * 2*pi * r_z)


                         do i_cell_1 = -number_of_super_cells(1),  number_of_super_cells(1)
                            do i_cell_2 = -number_of_super_cells(2),  number_of_super_cells(2)
                               do i_cell_3 = -number_of_super_cells(3),  number_of_super_cells(3) 

                                  if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then

                                     i_cell_n = i_cell_n + 1 


                                     if(packed_matrix_format == PM_none) then
                                        k_phase(i_cell_n, i_k_point) &
                                        & = k_phase_base(1, i_k_point)**i_cell_1 &
                                        & * k_phase_base(2, i_k_point)**i_cell_2 &
                                        & * k_phase_base(3, i_k_point)**i_cell_3
                                     end if


! SVL k-phase for exact exchange
! XR: not needed anymore
!                                     if(use_hf_kspace) then
!                                        k_phase_exx(i_cell_n, i_k_point) &
!                                        & = k_phase_base(1, i_k_point)**i_cell_1 &
!                                        & * k_phase_base(2, i_k_point)**i_cell_2 &
!                                        & * k_phase_base(3, i_k_point)**i_cell_3
!                                     end if

                                     cell_index(i_cell_n, 1) = i_cell_1
                                     cell_index(i_cell_n, 2) = i_cell_2
                                     cell_index(i_cell_n, 3) = i_cell_3

                                  end if
                               end do
                            end do
                         end do ! i_cell_1
                         ! SVL k-phase for exact exchange
                         if((use_hf_realspace .or. (use_full_spectrum .and. n_periodic .gt. 0)).and..not.use_symmetry_reduced_spg)then
                            i_cell_n = 1
                            do i_cell_1 = -(n_k_points_xyz_nosym(1)-1)/2, n_k_points_xyz_nosym(1)/2
                               do i_cell_2 = -(n_k_points_xyz_nosym(2)-1)/2, n_k_points_xyz_nosym(2)/2
                                  do i_cell_3 = -(n_k_points_xyz_nosym(3)-1)/2, n_k_points_xyz_nosym(3)/2
                                     
                                     if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then
                                        i_cell_n = i_cell_n + 1
                                        k_phase_exx(i_cell_n, i_k_point) &
                                        & = k_phase_base(1, i_k_point)**i_cell_1 &
                                        & * k_phase_base(2, i_k_point)**i_cell_2 &
                                        & * k_phase_base(3, i_k_point)**i_cell_3
                                     endif
                                  enddo
                               enddo
                            enddo

                         endif

                      end if ! k_number(i_x, i_y, i_z)> 0
                   end do
                end do  ! i_y
             end do  ! i_x
             do i_x = 1, n_k_points_xyz(1)
                do i_y = 1, n_k_points_xyz(2)
                   do i_z = 1, n_k_points_xyz(3)
                     if (.not. use_symmetry_reduced_spg)  then
                       i_normal = ik2irred_map(i_x, i_y, i_z)
                       ik2irred_map(i_x, i_y, i_z) = normal_to_stable(i_normal)
                     end if
                   end do
                end do  ! i_y
             end do  ! i_x
             if (allocated(normal_to_stable)) deallocate(normal_to_stable)
          end if  ! n_periodic > 0

          if(.not.use_symmetry_reduced_spg )then
             k_minus_q_point_list = k_point_list
          endif

          do i_x = 1,3,1
             if(n_k_points_xyz(i_x)>2 .or. k_points_offset(i_x)/=0 )then

                real_eigenvectors = .false.

             end if
          end do

          k_weights(:) = k_weights(:) / ( n_k_points_xyz(1)* n_k_points_xyz(2)* n_k_points_xyz(3))

       end if

       ! CC: Force eigenvectors to be complex even if this
       !      is mathematically not necessary (Im(evec) = 0)
       if ( flag_force_complex_eigenvectors ) then
         real_eigenvectors = .false.
         write(info_str,'(2X,A)') 'Forcing the eigenvectors in the calculation to be COMPLEX.'
         call localorb_info(info_str,use_unit,'(A)') 
       end if

       if(real_eigenvectors)then
          write(info_str,'(2X,A)') 'The eigenvectors in the calculations are REAL.'
          call localorb_info(info_str,use_unit,'(A)') 
       else
          write(info_str,'(2X,A)') 'The eigenvectors in the calculations are COMPLEX.'
          call localorb_info(info_str,use_unit,'(A)') 

       end if

       ! see runtime_choices - this is a global setting. 
       n_k_points_group = min(n_k_points, n_k_points_group)
       if (allocated(k_number)) deallocate(k_number)

    end if ! n_periodic

    !array containing the mapping between k index and (task, local k) indices
    allocate(k_point_loc(2,n_k_points))
    k_point_loc=-1

    if ((n_periodic.eq.0).or.(use_scalapack)) then
       ! dummy value only - no k points
       n_k_points_task = 1
       n_q_points_task = 1
    else 
       n_k_points_task = 0
       n_q_points_task = 0
       do i_k_point = 1, n_k_points, 1

          if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then
             n_k_points_task = n_k_points_task + 1
          end if
          k_point_loc(1,i_k_point)=MOD(i_k_point,n_tasks)
          k_point_loc(2,i_k_point)=(i_k_point-1)/n_tasks+1

       end do
       n_q_points_task = n_k_points_task

       if (n_k_points_task == 0) warn_idle_cpus = .true.

       write(info_str,'(2X,A,I4,A,I10)') '| K-points in task', myid, ':', n_k_points_task
       call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
 
       ! The next option loop the complete list of k-points in a periodic calculation.
       ! For now, ONLY do this for output_level full please. (or use a separate flag)
       ! Adding 10000 lines to the standard output on a slow screen can really be unpleasant ...
       ! We must query the output_level explicitly instead of just using OL_low, as in the
       ! initial_scf routine, the output_priority is already set to OL_low.
       if ( (output_level.eq.'full') .or. out_k_point_list ) then
         do i_k_point=1,n_k_points,1
           write(info_str,'(2X,A,I0,A,3(2X,F10.6),2X,A,F12.8)') '| k-point: ', i_k_point, ' at ', &
             &  k_point_list(i_k_point,1),  k_point_list(i_k_point,2),  k_point_list(i_k_point,3), &
             &  ', weight: ', k_weights(i_k_point)
           call localorb_info(info_str, use_unit, '(A)', OL_low)
         enddo
       end if

    end if

! SVL Create communicator for processes containing k-points

    if( (n_periodic.gt.0) .and. use_mpi .and. use_periodic_hf ) then

      if (prev_n_k_points.ne.n_k_points) then ! only set up first time, or if number of k points is changed
         call setup_kpoint_communicator ( prev_n_k_points, n_k_points )
         !     The resulting communicator is 'kpt_comm' and is declared in mpi_tasks.f90 .
         !     Its size is n_k_points (from dimensions.f90), only created if n_k_points < n_tasks .
      end if 
      prev_n_k_points = n_k_points ! remember the communicator size
 
   endif

   ! periodic post-DFT calculations involving unoccupied states
!   if (use_full_spectrum) then
     ! check allocation in case geometry is ever relaxed along the way
     if (allocated(n_states_k)) then
        deallocate(n_states_k)
      end if
     if ( n_periodic .gt. 0 ) then
        allocate(n_states_k(n_k_points))
     else
        ! need dummy allocation - apparently this variable is also touched in non-periodic cases
        allocate(n_states_k(1))
     end if
!   endif
   n_q_points = n_k_points

  end subroutine initialize_k_points
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/initialize_basis_lists
  !  NAME
  !    initialize_basis_lists
  !  SYNOPSIS

  subroutine initialize_basis_lists

    !  PURPOSE
    ! Initializes: 
    ! o  basis_atom
    ! o  Cbasis_to_center
    ! o  n_hamiltonian_matrix_size
    ! o  Cbasis_to_basis
    ! These tables are created for periodic boundary conditions, but in order to use same 
    ! subroutines also in clusters case, we have to allocate and initialize them.    
     use dimensions
     use localorb_io, only: localorb_info, use_unit
     use mpi_tasks, only: myid, n_tasks, check_allocation, aims_stop, SYNC_OR
     use synchronize_mpi_basic, only: sync_logical_vector
     use basis, only: basis_atom, basis_small_atom, basis_fn, basis_small_fn, outer_radius
     use runtime_choices, only: flag_rel, REL_4c_dks, REL_x2c
     use rel_x2c_mod, only: n_basis_small, n_centers_basis_I_small
     implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE




    integer:: i_cell_1, i_cell_2, i_cell_3, i_cell_n
    integer:: i_basis, i_basis_fn, i_Cbasis, i_Cbasis_full, i_center
    integer :: n_centers_basis_full, n_centers_basis_full_s

    real*8:: this_radius
    logical,dimension(:),allocatable:: include_Cbasis, include_Cbasis_s
    character*100 :: info_str
    real*8:: center(3)
    integer, allocatable :: cell_no(:)
    integer:: n_cell_no, info
    integer*8 :: matrix_size_i8

    character(*), parameter :: func = 'initialize_basis_list'

    ! --- Cell permutation

    ! JW: This is a permutation array where (roughly) the first half consists
    !     of representatives of equivalence classes with respect to inversion.
    !     Needed for the statement function Cbasis_to_center_full().
    !
    ! QUESTION: Does this permutation actually any good?

    allocate(cell_no(0:n_cells-1),stat=info)
    call check_allocation(info, 'cell_no', func)
    n_cell_no = 0
    cell_no(0) = 0   ! Central unit cell

    if(n_periodic .gt. 0)then

       ! basis which are included to the integrations
       i_cell_n = 0
       do i_cell_1 = -number_of_super_cells(1),  number_of_super_cells(1)
          do i_cell_2 = -number_of_super_cells(2),  number_of_super_cells(2)
             do i_cell_3 = -number_of_super_cells(3),  number_of_super_cells(3) 
                if(i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0) then
                   i_cell_n = i_cell_n + 1

                   if ((i_cell_1 > 0) .or. &
                   &   (i_cell_1 == 0 .and. i_cell_2 > 0) .or. &
                   &   (i_cell_1 == 0 .and. i_cell_2==0 .and. i_cell_3 > 0)) then
                      n_cell_no = n_cell_no+1
                      cell_no(n_cell_no) = i_cell_n
                   end if
                end if
             end do
          end do
       end do

       ! And the rest of them
       i_cell_n = 0
       do i_cell_1 = -number_of_super_cells(1),  number_of_super_cells(1)
          do i_cell_2 = -number_of_super_cells(2),  number_of_super_cells(2)
             do i_cell_3 = -number_of_super_cells(3),  number_of_super_cells(3) 
                if (i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0) then
                   i_cell_n = i_cell_n + 1

                   if (.not. ((i_cell_1 > 0) .or. &
                   &          (i_cell_1 == 0 .and. i_cell_2 > 0) .or. &
                   &          (i_cell_1 == 0 .and. i_cell_2==0 .and. i_cell_3 > 0))) then
                      n_cell_no = n_cell_no+1
                      cell_no(n_cell_no) = i_cell_n
                   end if
                end if
             end do
          end do
       end do

       ! --- include_Cbasis

       n_centers_basis_full = n_cells * n_basis
       ! n_centers_basis_I_full = ((n_cells-1)/2 + 1) * n_basis
       allocate(include_Cbasis(n_centers_basis_full),stat=info)
       call check_allocation(info, 'include_Cbasis', func)

       include_Cbasis = .false.

       do i_Cbasis_full = n_basis+1, n_centers_basis_full
          
          ! Distribute work over MPI tasks
          if(mod(i_Cbasis_full, n_tasks) /= myid) cycle

          i_center = Cbasis_to_center_full(i_Cbasis_full)
          center = coords_center(:,i_center)
          i_basis_fn = basis_fn(Cbasis_to_basis_full(i_Cbasis_full))
          this_radius = outer_radius(i_basis_fn)

          ! The delta term 0.25 is here to create a basis list that is safely
          ! larger than the real amount of basis functions ever touched by the
          ! grid.  Note that this is not performance critical - at each grid
          ! point, only non-zero basis functions are evaluated anyway.
          call ball_touches_zero_cell(center, this_radius, 0.25d0, &
          &                           include_Cbasis(i_Cbasis_full))

       end do

       call sync_logical_vector(include_Cbasis, n_centers_basis_full, SYNC_OR)

       ! Actually remove basis not needed.
       include_Cbasis(1:n_basis) = .true.
       n_centers_basis_T = count(include_Cbasis(:))
       ! n_centers_basis_I = count(include_Cbasis(1:n_centers_basis_I_full))

       if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then

          n_centers_basis_full_s = n_cells * n_basis_small
          allocate(include_Cbasis_s(n_centers_basis_full_s),stat=info)
          call check_allocation(info, 'include_Cbasis_s', func)

          include_Cbasis_s = .false.
  
          do i_Cbasis_full = n_basis_small+1, n_centers_basis_full_s
             
             ! Distribute work over MPI tasks
             if(mod(i_Cbasis_full, n_tasks) /= myid) cycle
  
             i_center = Cbasis_to_center_full_s(i_Cbasis_full)
             center = coords_center(:,i_center)
             i_basis_fn = basis_small_fn(Cbasis_to_basis_full_s(i_Cbasis_full))
             this_radius = outer_radius(i_basis_fn)

             call ball_touches_zero_cell(center, this_radius, 0.25d0, include_Cbasis_s(i_Cbasis_full))
  
          end do

          call sync_logical_vector(include_Cbasis_s, n_centers_basis_full_s, SYNC_OR) ! ask Victor
          include_Cbasis_s(1:n_basis_small) = .true.
          n_centers_basis_I_small = count(include_Cbasis_s(:))

       endif

    else

       n_centers_basis_T = n_basis
       n_centers_basis_I = n_basis
       if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) n_centers_basis_I_small = n_basis_small

    end if   ! n_periodic


    ! --- Prepare actual arrays

    call allocate_basis_lists

    ! Prepare 0-cell: Simply copy, i_center==i_atom
    do i_basis = 1, n_basis
       Cbasis_to_basis(i_basis) = i_basis
       Cbasis_to_atom(i_basis) = basis_atom(i_basis)
       Cbasis_to_center(i_basis) = basis_atom(i_basis)
    end do
    if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
       do i_basis = 1, n_basis_small
          Cbasis_to_basis_s(i_basis) = i_basis
          Cbasis_to_atom_s(i_basis) = basis_small_atom(i_basis)
          Cbasis_to_center_s(i_basis) = basis_small_atom(i_basis)
       end do
    endif

    if (n_periodic > 0) then

       ! Other cells, only included needed.
       i_Cbasis = n_basis
       do i_Cbasis_full = n_basis+1, n_centers_basis_full
          if(include_Cbasis(i_Cbasis_full))then
             i_Cbasis = i_Cbasis + 1
             Cbasis_to_basis(i_Cbasis) = Cbasis_to_basis_full(i_Cbasis_full)
             Cbasis_to_atom(i_Cbasis) = Cbasis_to_atom_full(i_Cbasis_full) 
             Cbasis_to_center(i_Cbasis) = Cbasis_to_center_full(i_Cbasis_full)
          end if
       end do
       n_centers_basis_I = i_Cbasis

       write(info_str,'(2X,A,I10)') &
       & '| Number of basis functions in the Hamiltonian integrals :', &
       & n_centers_basis_I
       call localorb_info(info_str,use_unit,'(A)')

       write(info_str,'(2X,A,I10)') &
       & '| Number of basis functions in a single unit cell        :', n_basis
       call localorb_info(info_str,use_unit,'(A)')

       ! VB - I thought n_centers_basis_T was no longer needed?
       !       write(info_str,'(2X,A,I10)') &
       !       & '| Number of basis functions totally - VB: What is this??   :', n_centers_basis_T
       !       call localorb_info(info_str,use_unit,'(A)') 

       deallocate(include_Cbasis)

       ! For relativistic basis, repeat the procedure for small comp. counterpart
       if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then

          i_Cbasis = n_basis_small
          do i_Cbasis_full = n_basis_small+1, n_centers_basis_full_s
             if(include_Cbasis_s(i_Cbasis_full))then
                i_Cbasis = i_Cbasis + 1
                Cbasis_to_basis_s(i_Cbasis) = Cbasis_to_basis_full_s(i_Cbasis_full)
                Cbasis_to_atom_s(i_Cbasis) = Cbasis_to_atom_full_s(i_Cbasis_full) 
                Cbasis_to_center_s(i_Cbasis) = Cbasis_to_center_full_s(i_Cbasis_full)
             end if
          end do
          n_centers_basis_I_small = i_Cbasis

          write(info_str,'(2X,A,I10)') &
          & '| Number of small component basis functions in the Hamiltonian integrals :', &
          & n_centers_basis_I_small
          call localorb_info(info_str,use_unit,'(A)')

          write(info_str,'(2X,A,I10)') &
          & '| Number of small component basis functions in a single unit cell        :', n_basis_small
          call localorb_info(info_str,use_unit,'(A)')
        
          deallocate(include_Cbasis_s)

       endif

    end if

    matrix_size_i8 = ((n_centers_basis_I+1)*n_centers_basis_I)/2
    !BL: Check if we are integer conform, otherwise stop here 
    if ( matrix_size_i8 .le. huge(n_hamiltonian_matrix_size)) then
       n_hamiltonian_matrix_size = int(matrix_size_i8)
    else
       call aims_stop('Integer overflow detected for n_hamiltonian_matrix_size. Stop here.', func)
    end if

    if (n_centers_basis_I /= n_centers_basis_T) then
       call aims_stop('Assertion failed: n_centers_basis_I /= n_centers_basis_T', func)
    end if
    n_Cbasis = n_centers_basis_I

  contains

    integer function Cbasis_to_basis_full(i_basis)
        integer :: i_basis
        Cbasis_to_basis_full = mod(i_basis-1,n_basis) + 1
    end function
    integer function Cbasis_to_atom_full(i_basis)
        integer :: i_basis
        Cbasis_to_atom_full = basis_atom(Cbasis_to_basis_full(i_basis))
    end function
    integer function Cbasis_to_center_full(i_basis)
        integer :: i_basis
        Cbasis_to_center_full = Cbasis_to_atom_full(i_basis) + cell_no((i_basis-1)/n_basis)*n_atoms
    end function

    integer function Cbasis_to_basis_full_s(i_basis)
        integer :: i_basis
        Cbasis_to_basis_full_s = mod(i_basis-1,n_basis_small) + 1
    end function
    integer function Cbasis_to_atom_full_s(i_basis)
        integer :: i_basis
        Cbasis_to_atom_full_s = basis_small_atom(Cbasis_to_basis_full_s(i_basis))
    end function
    integer function Cbasis_to_center_full_s(i_basis)
        integer :: i_basis
        Cbasis_to_center_full_s = Cbasis_to_atom_full_s(i_basis) + cell_no((i_basis-1)/n_basis_small)*n_atoms
    end function
  end subroutine initialize_basis_lists
    !******
  !-----------------------------------------------------------
  !****s* pbc_lists/initialize_centers_list
  !  NAME
  !      initialize_centers_list
  !  SYNOPSIS

  subroutine initialize_centers_list

    !  PURPOSE
    !  Calculates, how many centers is needed in different routines.
    !  Initializes:
    !  o  centers_ele_summation 
    !  o  centers_hartree_potential
    !  o  centers_hartree_multipole
    !  o  centers_basis_integrals
    !  o  centers_integrals
    !  o  centers_in_hamiltonian
    !  USES
    use dimensions
    use runtime_choices, only: extra_adding_to_hartree_potential_distance
    use localorb_io, only: localorb_info ,use_unit, OL_norm
    use mpi_tasks, only: myid, n_tasks, check_allocation, SYNC_OR
    use species_data
    use synchronize_mpi_basic, only: sync_logical_vector
    use basis, only: atom_radius
    use geometry, only: species, empty
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_center, i_center_1, i_center_2
    integer :: i_atom, i_atom_1, i_atom_2, i_centers
    integer :: i_species, i_species_2
    integer:: i_hartree_pot, i_ele_sum, i_basis_integrals, i_hartree_multipole

    real*8:: this_radius, dist_min, dist, center(3)

    logical:: include_hartree_potential, include_ele_summation
    logical:: include_basis_integrals, include_hartree_multipole


    integer :: centers_basis_integrals_other(n_centers)
    logical:: include_center_in_integral(n_centers)

    character*100 :: info_str
    integer :: current_atom, current_center
    integer :: atom_min
    integer :: center_swap
    integer :: n_centers_basis_integrals_other
    integer :: center_temp, i_centers_integrals, i_centers_basis_integrals_other

    logical :: ghost_atom, pseudo_atom

    integer :: i_occ_atom, i_empty_atom, i_pp_atom
    character(*), parameter :: func = 'initialize_centers_list'
 
   !---------shanghui add for DFPT_phonon------------------
    integer  :: i_centers_in_hamiltonian
    logical  :: include_in_hamiltonian
   !---------shanghui end add for DFPT_phonon--------------



    ! First, find dimensions of lists of "centers", to be allocated later

    n_centers_ele_summation             = n_atoms
    n_centers_hartree_potential         = n_atoms
    n_centers_hartree_multipole         = n_atoms
    ! Make sure that in centers_basis_integrals, those centers whose ball does
    ! not touch the 0-cell come last (..._other).
    n_centers_integrals                 = n_atoms
    n_centers_basis_integrals_other     = 0
    n_centers_in_hamiltonian            = n_atoms

    if(n_periodic > 0)then

       ! --- Get array sizes
       include_center_in_integral = .false.

       ! Here we go over surfaces and calculate all centers needed for splines
       ! after mapping to center cell.
       do i_center = n_atoms+1, n_centers

          if(mod(i_center,n_tasks) /= myid) cycle ! distribute work over MPI tasks
          center = coords_center(:, i_center)
          this_radius = multipole_radius_free(species_center(i_center))
          call ball_touches_zero_cell(center, this_radius, 0.25d0, &
          &                           include_center_in_integral(i_center))

       end do
       call sync_logical_vector(include_center_in_integral,n_centers, SYNC_OR)

       ! Next, create tables of species-dependent limiting distances for use in get_include_centers.
       ! This tables should save some recomputation in that routine. 
       !
       ! Issues that are not yet addressed in get_include_centers are:
       ! - the execution is not parallel (far too expensive)
       ! - the main loop in get_include atoms depends on the order of the
       !   atoms tabulated from 1 to n_atoms. The meaning of the line
       !
       !        dist_min_sq = min(dist_sq, dist_min_sq)
       !
       ! is unclear. It looks like this like should either be left away, or the 
       ! absolute minimum distance of a given center should be determined completely
       ! before any decisions are made based on this value.

       ! Allocate the aforementioned tables.
       if (.not.allocated(dist_sq_basis_integrals)) then
          allocate(dist_sq_basis_integrals(n_species,n_species))
       end if
       if (.not.allocated(dist_sq_ele_summation)) then
          allocate(dist_sq_ele_summation(n_species,n_species))
       end if
       if (.not.allocated(dist_sq_hartree_potential)) then
          allocate(dist_sq_hartree_potential(n_species,n_species))
       end if
       if (.not.allocated(dist_sq_hartree_multipole)) then
          allocate(dist_sq_hartree_multipole(n_species,n_species))
       end if

       ! Populate the aforementioned tables.
       do i_species = 1, n_species, 1
          do i_species_2 = 1, n_species, 1

             dist_sq_basis_integrals(i_species, i_species_2) = & 
               ( atom_radius(i_species_2) + multipole_radius_free(i_species) )**2

             dist_sq_ele_summation(i_species, i_species_2) = &
               ( atom_radius(i_species_2) + atom_radius(i_species) )**2 
            
             dist_sq_hartree_potential(i_species,i_species_2) = &
               ( atom_radius(i_species_2) & 
                 + max(atom_radius_hartree(i_species), multipole_radius_free(i_species)) & 
                 + extra_adding_to_hartree_potential_distance )**2

             dist_sq_hartree_multipole(i_species,i_species_2) = &
               ( atom_radius(i_species_2) &
                 + max(atom_radius(i_species), multipole_radius_free(i_species)) & 
                 + extra_adding_to_hartree_potential_distance )**2

          enddo
       enddo
       max_atom_diameter_sq = ( 2.0d0* maxval(atom_radius) )**2

       do i_center = n_atoms+1, n_centers,1
          i_species = species_center(i_center)
          call get_include_centers(coords_center(:,i_center), i_species, &
          &                        include_basis_integrals, include_ele_summation, &
          &                        include_hartree_potential, include_hartree_multipole, &
          &                        include_in_hamiltonian)

          if (include_hartree_potential) then
             n_centers_hartree_potential = n_centers_hartree_potential + 1
          end if
          if (include_hartree_multipole) then
             n_centers_hartree_multipole = n_centers_hartree_multipole + 1
          end if
          if (include_ele_summation) then
             n_centers_ele_summation = n_centers_ele_summation + 1
             atoms_in_supercell(i_species) = atoms_in_supercell(i_species)+1  ! SVL for periodic EXX
          end if
          if (include_center_in_integral(i_center)) then
             n_centers_integrals = n_centers_integrals + 1
          else if (include_basis_integrals) then
             n_centers_basis_integrals_other = n_centers_basis_integrals_other + 1
          end if

          if (include_in_hamiltonian) then 
             n_centers_in_hamiltonian = n_centers_in_hamiltonian + 1
          endif 

       end do

    end if   ! n_periodic

    n_centers_basis_integrals = n_centers_integrals + n_centers_basis_integrals_other

    ! --- Output array sizes

    write(info_str,'(2X,A,I10)') '| Number of centers in hartree potential         :',n_centers_hartree_potential
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,I10)') '| Number of centers in hartree multipole         :',n_centers_hartree_multipole
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,I10)') '| Number of centers in electron density summation:',n_centers_ele_summation
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,I10)') '| Number of centers in basis integrals           :',n_centers_basis_integrals
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,I10)') '| Number of centers in integrals                 :',n_centers_integrals
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A,I10)') '| Number of centers in hamiltonian               :',n_centers_in_hamiltonian
    call localorb_info(info_str,use_unit,'(A)',OL_norm)

    ! Can now allocate lists of integration centers, Hartree potential centers, etc.    
    call allocate_center_lists ( )

    ! Initialize number to the list
    do i_atom = 1,n_atoms
       centers_hartree_potential(i_atom)         = i_atom
       centers_hartree_multipole(i_atom)         = i_atom
       centers_ele_summation(i_atom)             = i_atom
       centers_basis_integrals(i_atom)           = i_atom

       centers_in_hamiltonian(i_atom)            = i_atom
       inv_centers_in_hamiltonian(i_atom)        = i_atom
    end do


    if(n_periodic .gt. 0)then

       i_centers = n_atoms
       i_hartree_pot = n_atoms
       i_hartree_multipole = n_atoms
       i_ele_sum = n_atoms
       i_basis_integrals = n_atoms
       i_centers_integrals = n_atoms
       i_centers_basis_integrals_other = 0

       i_centers_in_hamiltonian = n_atoms

       ! --- Actually populate arrays

       do i_center = n_atoms+1, n_centers,1
          i_species = species_center(i_center)
          call get_include_centers(coords_center(:,i_center), i_species, &
          &                        include_basis_integrals, include_ele_summation, &
          &                        include_hartree_potential, include_hartree_multipole, &
          &                        include_in_hamiltonian)


          if(include_hartree_potential)then
             i_hartree_pot = i_hartree_pot + 1
             centers_hartree_potential(i_hartree_pot) =  i_center
          end if
          if(include_hartree_multipole)then
             i_hartree_multipole = i_hartree_multipole + 1
             centers_hartree_multipole(i_hartree_multipole) =  i_center
          end if
          if(include_ele_summation)then
             i_ele_sum = i_ele_sum + 1
             centers_ele_summation(i_ele_sum) = i_center
          end if
          if (include_center_in_integral(i_center)) then
             i_centers_integrals = i_centers_integrals + 1
             centers_basis_integrals(i_centers_integrals) = i_center
          else if (include_basis_integrals) then
             i_centers_basis_integrals_other = i_centers_basis_integrals_other + 1
             centers_basis_integrals_other(i_centers_basis_integrals_other) = i_center
          end if
 
          if (include_in_hamiltonian) then
             i_centers_in_hamiltonian = i_centers_in_hamiltonian + 1
             centers_in_hamiltonian(i_centers_in_hamiltonian) = i_center
             inv_centers_in_hamiltonian(i_center) = i_centers_in_hamiltonian
          endif
 
       end do


       ! --- Reorder centers_hartree_potential

       ! Here we order the array centers_hartree_potential so that centers relating to
       ! a same atom follow each other, this will help the distributed spline storage routines

       ! The forces need the real atoms to be first.

       do i_center_1 = n_atoms+1, n_centers_hartree_potential, 1

          atom_min = n_atoms
          center_swap = i_center_1

          do i_center_2 = i_center_1, n_centers_hartree_potential, 1

             current_center = centers_hartree_potential(i_center_2)
             current_atom = center_to_atom(current_center)
             if (current_atom < atom_min) then
                atom_min = current_atom
                center_swap = i_center_2
             end if
          end do

          if (center_swap.ne.i_center_1) then
             center_temp = centers_hartree_potential(i_center_1)
             centers_hartree_potential(i_center_1) = centers_hartree_potential(center_swap)
             centers_hartree_potential(center_swap) = center_temp
          end if

       end do

       ! deallocate any special lists only created for use in 
       ! get_include_centers .
       if (allocated(dist_sq_basis_integrals)) then
          deallocate(dist_sq_basis_integrals)
       end if
       if (allocated(dist_sq_ele_summation)) then
          deallocate(dist_sq_ele_summation)
       end if
       if (allocated(dist_sq_hartree_potential)) then
          deallocate(dist_sq_hartree_potential)
       end if
       if (allocated(dist_sq_hartree_multipole)) then
          deallocate(dist_sq_hartree_multipole)
       end if

    end if  ! n_periodic


    ! --- Finalize center_basis_integrals

    ! centers_basis_integrals
    ! = [1 .. n_atoms]
    !   + [i_center : include_center_in_integral(i_center)]
    !   + [i_center : include_basis_integrals .and.
    !                 .not. include_center_in_integral(i_center)]
    centers_basis_integrals(n_centers_integrals+1:n_centers_basis_integrals) &
    & = centers_basis_integrals_other(1:n_centers_basis_integrals_other)

    ! --- Invert center_basis_integrals permutation

    inv_centers_basis_integrals = 0
    do i_center_1= 1, n_centers
       do i_center_2 =  1, n_centers_basis_integrals
          if (centers_basis_integrals(i_center_2) == i_center_1) then
             inv_centers_basis_integrals(i_center_1) = i_center_2
          end if
       end do
    end do


    ! --- occ/empty atoms

    !   Here I reorganize the whole list, because of ghost atoms and real ones (Andrea Sanfilippo)
    n_occ_centers_integr                 = n_occ_atoms
    n_occ_centers_basis_integr           = n_occ_atoms 


    if(n_periodic .gt. 0) then


       i_centers = n_atoms
       i_hartree_pot = n_atoms
       i_hartree_multipole = n_atoms
       i_ele_sum = n_atoms
       i_basis_integrals = n_atoms
       i_centers_integrals = n_atoms


       do i_atom_1 = n_atoms+1, n_centers,1
          dist_min = 10.0e8      
          include_hartree_potential = .false.
          include_hartree_multipole = .false.
          include_ele_summation = .false.
          include_basis_integrals = .false.


          do i_atom_2 = 1,n_atoms


             dist = sqrt( (coords_center(1,i_atom_1) - coords_center(1,i_atom_2))**2 + &
                  (coords_center(2,i_atom_1) - coords_center(2,i_atom_2))**2 + &
                  (coords_center(3,i_atom_1) - coords_center(3,i_atom_2))**2)

             dist_min = min(dist,dist_min)

             if( dist_min <=  atom_radius(species(i_atom_2))+ multipole_radius_free(species_center(i_atom_1)))then
                include_basis_integrals   = .true.
             end if


          end do


          if( include_center_in_integral(i_atom_1).and..not.(empty(center_to_atom(i_atom_1)))) then
             n_occ_centers_integr = n_occ_centers_integr + 1
          end if



          if (include_basis_integrals .and.  &
              (.not. include_center_in_integral(i_atom_1)).and. &
               .not.(empty(center_to_atom(i_atom_1)))) then
             n_occ_centers_basis_integr = n_occ_centers_basis_integr + 1
          end if

       end do


       n_occ_centers_basis_integr = n_occ_centers_basis_integr + n_occ_centers_integr - n_occ_atoms
    endif

    i_occ_atom=0
    i_empty_atom=n_occ_centers_basis_integr
    i_pp_atom = n_centers_basis_integrals - n_pp_in_qm !pseudocores are always listed after empty atoms 

    do i_center=1,n_centers_basis_integrals

       ghost_atom = empty(center_to_atom(centers_basis_integrals(i_center)))
       pseudo_atom = species_pseudoized(species(center_to_atom(centers_basis_integrals(i_center))))
       if (ghost_atom) then
          i_empty_atom=i_empty_atom+1 
          new_centers_basis_integrals(i_empty_atom)= centers_basis_integrals(i_center)
       else if (pseudo_atom) then
          i_pp_atom=i_pp_atom+1
          new_centers_basis_integrals(i_pp_atom)= centers_basis_integrals(i_center)
       else
          i_occ_atom=i_occ_atom+1
          new_centers_basis_integrals(i_occ_atom)= centers_basis_integrals(i_center)
       endif
    enddo


    return

  contains

    !----------------------------------------------------------------------------
    !****s* pbc_lists/get_include_centers
    !  NAME
    !    get_include_centers
    !  SYNOPSIS

    subroutine get_include_centers(center, i_species, &
    &                              include_basis_integrals, include_ele_summation, &
    &                              include_hartree_potential, include_hartree_multipole, & 
                                   include_in_hamiltonian)

      !  PURPOSE
      !    Figure out if a given center is needed.
      !  USES

      use geometry, only: species
      implicit none

      !  ARGUMENTS

      real*8, intent(IN) :: center(3)
      integer, intent(IN) :: i_species
      logical, intent(OUT) :: include_basis_integrals
      logical, intent(OUT) :: include_ele_summation
      logical, intent(OUT) :: include_hartree_potential
      logical, intent(OUT) :: include_hartree_multipole
      logical, intent(OUT) :: include_in_hamiltonian
      
      !  INPUTS
      !    o center_1 -- Position of center
      !    o i_species -- Atomic species of center
      !  OUTPUTS  ["ovlp..." == "overlap with some 0-cell atom"]
      !    o include_basis_integrals -- Center's atomic charge ovlp...
      !    o include_ele_summation -- One-particle overlap with some atom
      !    o include_hartree_potential -- Center's Hartree potential ovlp...
      !    o include_hartree_multipole -- Center's multipole charge ovlp...
      !    o include_in_hamiltonian -- include_centers_in_hamiltonian
      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2011).
      !  SOURCE

      integer :: i_atom_2, i_species_2
      real*8 :: dist_min_sq, distvec(3), dist_sq
      integer :: info
      character(*), parameter :: func = 'get_include_centers'

      dist_min_sq = 1.0d18
      include_hartree_potential = .false.
      include_hartree_multipole = .false.
      include_ele_summation = .false.
      include_basis_integrals = .false.
      include_in_hamiltonian = .false.

      do i_atom_2 = 1,n_atoms
         i_species_2 = species(i_atom_2)
         distvec = center - coords_center(:,i_atom_2)
         dist_sq = sum(distvec**2)
         dist_min_sq = min(dist_sq, dist_min_sq)

         if (dist_min_sq <= dist_sq_basis_integrals(i_species,i_species_2)) then
            include_basis_integrals = .true.
         end if
         if (dist_min_sq <= dist_sq_ele_summation (i_species, i_species_2) ) then
            include_ele_summation = .true.
         end if
         if (dist_min_sq <= dist_sq_hartree_potential(i_species,i_species_2) ) then
            include_hartree_potential = .true.
         end if
         if (dist_min_sq <= dist_sq_hartree_multipole(i_species,i_species_2) ) then
            include_hartree_multipole = .true.
         end if
         !----------shanghui add for DFPT_phonon, ref(initialize_position_in_hamiltonian)---------
         if (dist_min_sq <= max_atom_diameter_sq) then 
            include_in_hamiltonian = .true.
         end if
         !----------shanghui end add for DFPT_phonon, ref(initialize_position_in_hamiltonian)---------
      end do

    end subroutine get_include_centers
    !******
  end subroutine initialize_centers_list


  !******
  !****s* pbc_lists/initialize_sc_DFPT_phonon
  !  NAME
  !     initialize_sc_DFPT_phonon
  !  SYNOPSIS
      ! initialize_sc_DFPT_phonon
  subroutine initialize_sc_DFPT_phonon 
    !  PURPOSE
    !    Get the supercell using in DFPT_phonon 
    ! 
    use dimensions
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: check_allocation, aims_stop
    use basis, only: atom_radius
    use geometry, only: coords, lattice_vector
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    logical,dimension(:), allocatable:: include_cell
    integer,dimension(:), allocatable:: new_index
    integer,dimension(:,:),allocatable:: cell_index_hamiltonian

    integer,dimension(:,:,:), allocatable:: cell_index_to_i_cell_in_sc_DFPT

    integer :: i_coord, delta_cell(3),cell_index_PBC(3)
    integer :: i_cell_1, i_cell_2, i_cell_3
    integer :: i_atom_1, i_atom_2, i_atom
    integer :: i_cell,i_cell_new , i_cell_in_sc_DFPT, i_cell_in_hamiltonian
    integer :: i_center,  i_center_in_sc_DFPT 
    real*8  :: dist_vec(3), dist_atom_vec(3), dist, min_dist
    integer :: io,i_basis


    character*100 :: info_str
    integer :: info
    character(*), parameter :: func = 'initialize_sc_DFPT_phonon'


!------------------------(1) calculate number_of_sc_DFPT------------------
    allocate(include_cell(n_cells),stat=info)
    call check_allocation(info, 'include_cell', func)

    allocate(new_index(n_cells),stat=info)
    call check_allocation(info, 'new_index', func)

    include_cell = .false.
    do i_cell = 1,n_cells,1
         dist_vec = cell_index(i_cell,1) * lattice_vector(1:3,1) &
                  + cell_index(i_cell,2) * lattice_vector(1:3,2) &
                  + cell_index(i_cell,3) * lattice_vector(1:3,3)
         min_dist = 100000.d0
         do i_atom_1 = 1,n_atoms
            do i_atom_2 = 1,n_atoms
               dist_atom_vec = coords(1:3, i_atom_1) - coords(1:3,i_atom_2) + dist_vec
               dist = dot_product(dist_atom_vec, dist_atom_vec)
               min_dist = min( min_dist, dist)
            end do
         end do

         dist = sqrt( min_dist )
         if(dist < 2* maxval(atom_radius))then
            include_cell(i_cell) = .true.
         end if
     end do

     i_cell_new = 0
     new_index = 0
     do i_cell = 1,n_cells,1
       if(include_cell(i_cell))then
          i_cell_new = i_cell_new + 1
          new_index(i_cell) = i_cell_new
       end if
     end do
     n_cells_in_hamiltonian = i_cell_new + 1

     write(info_str,*) " |-----begin in initialize_sc_DFPT_phonon----------"
     call localorb_info(info_str,use_unit,'(A)')
     write(info_str,*) " | [n_cells_in_hamiltonian]:", n_cells_in_hamiltonian
     call localorb_info(info_str,use_unit,'(A)')

     allocate(cell_index_hamiltonian(n_cells_in_hamiltonian,3),stat=info)
     call check_allocation(info, 'cell_index_hamiltonian',func)

     do i_cell = 1,n_cells,1
       if(new_index(i_cell)>0)then
          cell_index_hamiltonian(new_index(i_cell),1:3) = cell_index(i_cell,1:3)
       end if
     end do

     do i_coord = 1,3
     min_nsc_DFPT(i_coord) = minval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,i_coord))
     max_nsc_DFPT(i_coord) = maxval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,i_coord))
     enddo       

    !------------begin siesta method for sc-----------------------------
    !do i=1,3
    !   veclen = sqrt(lattice_vector(1,i)**2+lattice_vector(2,i)**2+lattice_vector(3,i)**2)
    !   write(use_unit,*) 'how many cells in half:   ',i,  2 * maxval(atom_radius) / veclen   
    !   write(use_unit,*) 'how many cells in seista: ',i,  ceiling( 2* 2 * maxval(atom_radius) / veclen )  
    !enddo
    !------------end siesta method for sc-----------------------------
 
    write(info_str,*) " | Number of origin super-cells(1) :", '[',-number_of_super_cells(1),':', &
                                                           number_of_super_cells(1), ']' 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | Number of origin super-cells(2) :", '[',-number_of_super_cells(2),':', &
                                                           number_of_super_cells(2), ']' 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | Number of origin super-cells(3) :", '[',-number_of_super_cells(3),':', &
                                                           number_of_super_cells(3), ']' 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | Number of  DFPT  super-cells(1) :", '[',min_nsc_DFPT(1),':', &
                                                           max_nsc_DFPT(1), ']' 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | Number of  DFPT  super-cells(2) :", '[',min_nsc_DFPT(2),':', &
                                                           max_nsc_DFPT(2), ']' 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | Number of  DFPT  super-cells(3) :", '[',min_nsc_DFPT(3),':', &
                                                           max_nsc_DFPT(3), ']' 
    call localorb_info(info_str,use_unit,'(A)')

 
!------------------------(2) calculate n_cells_in_sc_DFPT------------------

    do i_coord = 1,3 
       nsc_DFPT(i_coord) = max_nsc_DFPT(i_coord)-min_nsc_DFPT(i_coord) + 1
    enddo
 
    n_cells_in_sc_DFPT = nsc_DFPT(1)*nsc_DFPT(2)*nsc_DFPT(3)


!------------------------(3) calculate n_centers_in_sc_DFPT  <==> n_centers ------------------
    n_centers_in_sc_DFPT = n_atoms * n_cells_in_sc_DFPT

    write(info_str,*) " | nsc_DFPT                        :", nsc_DFPT(1:3)
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | n_cells_in_sc_DFPT              :", n_cells_in_sc_DFPT 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | n_center_in_sc_DFPT             :", n_centers_in_sc_DFPT 
    call localorb_info(info_str,use_unit,'(A)')

    allocate(center_to_center_in_sc_DFPT(n_centers),stat=info)
    call check_allocation(info, 'centers_to_centers_in_sc_DFPT',func)

    allocate(center_in_sc_DFPT_to_center(n_centers_in_sc_DFPT),stat=info)
    call check_allocation(info, 'center_in_sc_DFPT_to_centers',func)

    center_in_sc_DFPT_to_center(1:n_centers_in_sc_DFPT) = 0  ! n_centers_in_sc_DFPT ==> n_centers 
    center_to_center_in_sc_DFPT(1:n_centers) =0              ! n_centers            ==> n_centers_in_sc_DFPT

    i_center_in_sc_DFPT = 0
    do i_center = 1, n_centers 
        !note:  (1) center_to_cell is given in initialize_centers_data(), 
        !           cell_index is given in initialize_k_points() 
        !note:  (2) center_to_cell and cell_index are both changed in initialize_position_in_hamiltonian
        !note:  (3) our initialize_sc_DFPT_phonon is after (1) and before (2), so it is the origin_cell_index
        i_cell_1 = cell_index(center_to_cell(i_center),1)  
        i_cell_2 = cell_index(center_to_cell(i_center),2) 
        i_cell_3 = cell_index(center_to_cell(i_center),3) 

          if(i_cell_1.ge.min_nsc_DFPT(1).and. i_cell_1.le.max_nsc_DFPT(1) .and. &  
             i_cell_2.ge.min_nsc_DFPT(2).and. i_cell_2.le.max_nsc_DFPT(2) .and. &
             i_cell_3.ge.min_nsc_DFPT(3).and. i_cell_3.le.max_nsc_DFPT(3) ) then

             i_center_in_sc_DFPT = i_center_in_sc_DFPT + 1

             center_in_sc_DFPT_to_center(i_center_in_sc_DFPT) = i_center 

             center_to_center_in_sc_DFPT(i_center) = i_center_in_sc_DFPT
           endif  

    enddo 

!-------------------(4) n_center_in_sc_DFPT <==> n_cells_in_sc_DFPT--------------------
!                       cell_index_sc_DFPT  <==> i_cell_in_sc_DFPT 

    allocate(cell_index_sc_DFPT(n_cells_in_sc_DFPT,3),stat=info)
    call check_allocation(info, 'cell_index_sc_DFPT',func)
    allocate(cell_index_to_i_cell_in_sc_DFPT(min_nsc_DFPT(1):max_nsc_DFPT(1)+1, &
                                             min_nsc_DFPT(2):max_nsc_DFPT(2)+1, &
                                             min_nsc_DFPT(3):max_nsc_DFPT(3)+1) &
                                             ,stat=info)
    call check_allocation(info, 'cell_index_to_i_cell_in_sc_DFPT', func)


    allocate(center_in_sc_DFPT_to_atom(n_centers_in_sc_DFPT),stat=info)
    call check_allocation(info, 'center_in_sc_DFPT_to_atom', func)
    allocate(center_in_sc_DFPT_to_cell_in_sc_DFPT(n_centers_in_sc_DFPT),stat=info)
    call check_allocation(info, 'center_in_sc_DFPT_to_cell_in_sc_DFPT', func)
    allocate(cell_and_atom_to_center_sc_DFPT(n_cells_in_sc_DFPT,n_atoms),stat=info)
    call check_allocation(info, 'cell_and_atom_to_center_sc_DFPT', func)


 
    center_in_sc_DFPT_to_atom = 0
    center_in_sc_DFPT_to_cell_in_sc_DFPT = 0 
    cell_and_atom_to_center_sc_DFPT = 0
    
    i_center_in_sc_DFPT  = 0
    i_cell_in_sc_DFPT = 1

    cell_index_sc_DFPT(i_cell_in_sc_DFPT,1) = 0
    cell_index_sc_DFPT(i_cell_in_sc_DFPT,2) = 0
    cell_index_sc_DFPT(i_cell_in_sc_DFPT,3) = 0
    cell_index_to_i_cell_in_sc_DFPT(0,0,0)                        = i_cell_in_sc_DFPT
    do i_atom = 1,n_atoms
       i_center_in_sc_DFPT = i_center_in_sc_DFPT + 1
       center_in_sc_DFPT_to_atom(i_center_in_sc_DFPT)             =  i_atom
       center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center_in_sc_DFPT)  =  i_cell_in_sc_DFPT
       cell_and_atom_to_center_sc_DFPT(i_cell_in_sc_DFPT,i_atom)  =  i_center_in_sc_DFPT
    end do

    
    do i_cell_1 =  min_nsc_DFPT(1), max_nsc_DFPT(1)
    do i_cell_2 =  min_nsc_DFPT(2), max_nsc_DFPT(2)
    do i_cell_3 =  min_nsc_DFPT(3), max_nsc_DFPT(3)
       if(  i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0)then

          i_cell_in_sc_DFPT = i_cell_in_sc_DFPT + 1

          cell_index_sc_DFPT(i_cell_in_sc_DFPT,1) = i_cell_1
          cell_index_sc_DFPT(i_cell_in_sc_DFPT,2) = i_cell_2
          cell_index_sc_DFPT(i_cell_in_sc_DFPT,3) = i_cell_3

          cell_index_to_i_cell_in_sc_DFPT(i_cell_1,i_cell_2,i_cell_3)   =  i_cell_in_sc_DFPT


          do i_atom = 1,n_atoms
             i_center_in_sc_DFPT = i_center_in_sc_DFPT + 1
             center_in_sc_DFPT_to_atom(i_center_in_sc_DFPT)             =  i_atom
             center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center_in_sc_DFPT)  =  i_cell_in_sc_DFPT
             cell_and_atom_to_center_sc_DFPT(i_cell_in_sc_DFPT,i_atom)  =  i_center_in_sc_DFPT 
          enddo 

       endif 
    enddo 
    enddo 
    enddo 



!------------------------(5) calculate cell_diff_sc_DFPT ------------------
    allocate(cell_diff_sc_DFPT(n_cells_in_sc_DFPT, n_cells_in_sc_DFPT),stat=info)
    call check_allocation(info, 'cell_diff_sc_DFPT', func)
   
    cell_diff_sc_DFPT(1:n_cells_in_sc_DFPT,1:n_cells_in_sc_DFPT) = 0
   

    do i_cell_1 = 1,n_cells_in_sc_DFPT
    do i_cell_2 = 1,n_cells_in_sc_DFPT

      do i_coord = 1, 3, 1

         delta_cell(i_coord) = cell_index_sc_DFPT(i_cell_1,i_coord) -  &  ! <-----minus 
                               cell_index_sc_DFPT(i_cell_2,i_coord)
       if(mod(nsc_DFPT(i_coord),2).ne.1) then
          call aims_stop('nsc_DFPT not odd, please check', func)
       endif

       if(abs(delta_cell(i_coord)).le.(nsc_DFPT(i_coord)-1)/2) then
         cell_index_PBC(i_coord) = delta_cell(i_coord)
       else
         cell_index_PBC(i_coord) = &
         - mod( nsc_DFPT(i_coord), abs(delta_cell(i_coord)) )* &
         delta_cell(i_coord)/abs(delta_cell(i_coord))
       endif

      enddo

       cell_diff_sc_DFPT(i_cell_1, i_cell_2) = &
       cell_index_to_i_cell_in_sc_DFPT(cell_index_PBC(1),cell_index_PBC(2),cell_index_PBC(3))

    end do
    end do



!------------------------(6) calculate cell_add_sc_DFPT ------------------
    allocate(cell_add_sc_DFPT(n_cells_in_sc_DFPT, n_cells_in_sc_DFPT),stat=info)
    call check_allocation(info, 'cell_add_sc_DFPT', func)
   
    cell_add_sc_DFPT(1:n_cells_in_sc_DFPT,1:n_cells_in_sc_DFPT) = 0

    do i_cell_1 = 1,n_cells_in_sc_DFPT
    do i_cell_2 = 1,n_cells_in_sc_DFPT

      do i_coord = 1, 3, 1

         delta_cell(i_coord) = cell_index_sc_DFPT(i_cell_1,i_coord) +  &  ! <-----add
                               cell_index_sc_DFPT(i_cell_2,i_coord)
       if(mod(nsc_DFPT(i_coord),2).ne.1) then
          call aims_stop('nsc_DFPT not odd, please check', func)
       endif

       if(abs(delta_cell(i_coord)).le.(nsc_DFPT(i_coord)-1)/2) then
         cell_index_PBC(i_coord) = delta_cell(i_coord)
       else
         cell_index_PBC(i_coord) = &
         - mod( nsc_DFPT(i_coord), abs(delta_cell(i_coord)) )* &
         delta_cell(i_coord)/abs(delta_cell(i_coord))
       endif
      enddo

       cell_add_sc_DFPT(i_cell_1, i_cell_2) = &
       cell_index_to_i_cell_in_sc_DFPT(cell_index_PBC(1),cell_index_PBC(2),cell_index_PBC(3))

    end do
    end do



!------------------------------(7) calculate basis_in_supercell---------------------
    allocate(cell_and_basis_to_basis_sc_DFPT(n_cells_in_sc_DFPT,n_basis),stat=info)
    call check_allocation(info, 'cell_and_basis_to_basis_sc_DFPT', func)

    cell_and_basis_to_basis_sc_DFPT = 0
    io = 0
    do i_cell_in_sc_DFPT = 1,n_cells_in_sc_DFPT,1
    do i_basis = 1,n_basis
       io = io + 1
       cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,i_basis) = io
    enddo
    enddo
    n_basis_sc_DFPT = io

    write(info_str,*) " | n_basis in unit cell           :", n_basis
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " | n_basis in sc_DFPT             :", n_basis_sc_DFPT 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,*) " |-----end in initialize_sc_DFPT_phonon---"
    call localorb_info(info_str,use_unit,'(A)')


!---------------------------(8) n_cells_in_hamiltonian ====> n_cells_in_sc_DFPT----------------
    allocate(cell_in_hamiltonian_to_cell_in_sc_DFPT(n_cells_in_hamiltonian),stat=info)
    call check_allocation(info, 'cell_in_hamiltonian_to_cell_in_sc_DFPT', func)

    cell_in_hamiltonian_to_cell_in_sc_DFPT(1:n_cells_in_hamiltonian) = 0
    
    do i_cell_in_hamiltonian  = 1,n_cells_in_hamiltonian -1 

       do i_cell_in_sc_DFPT   = 1,n_cells_in_sc_DFPT 
 
       if(     cell_index_sc_DFPT(i_cell_in_sc_DFPT,1).eq. & 
           cell_index_hamiltonian(i_cell_in_hamiltonian,1) .and. & 
               cell_index_sc_DFPT(i_cell_in_sc_DFPT,2).eq. & 
           cell_index_hamiltonian(i_cell_in_hamiltonian,2) .and. & 
               cell_index_sc_DFPT(i_cell_in_sc_DFPT,3).eq. & 
           cell_index_hamiltonian(i_cell_in_hamiltonian,3) ) then  

           cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian) = i_cell_in_sc_DFPT 
           
        endif 
        
        enddo 
     enddo 




    deallocate(include_cell)
    deallocate(new_index)
    deallocate(cell_index_hamiltonian)
    deallocate(cell_index_to_i_cell_in_sc_DFPT)


  end subroutine initialize_sc_DFPT_phonon

  !******
  !****s* pbc_lists/initialize_position_in_hamiltonian_PBC
  !  NAME
  !     initialize_position_in_hamiltonian_PBC
  !  SYNOPSIS
      ! initialize_position_in_hamiltonian_PBC
  subroutine initialize_position_in_hamiltonian_PBC !DFPT_phonon_PBC
    !  PURPOSE
    !    Get the position_in_hamiltonian_PBC. 
    !    Here we have a more neat version code for n_cells_in_hamiltonian, you can compare with the one 
    !    get in initialize_position_in_hamiltonian.
    ! 
    use dimensions
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: check_allocation, aims_stop
    use basis, only: atom_radius
    use geometry, only: lattice_vector, coords
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    logical,dimension(:), allocatable:: include_cell
    integer,dimension(:), allocatable:: new_index
    integer,dimension(:,:),allocatable:: cell_index_hamiltonian
    integer,dimension(:,:,:), allocatable:: cell_index_to_i_cell


    integer:: i_coord, delta_cell(3),cell_index_PBC(3), n_cells_in_hamiltonian_array(3)

    integer:: i_cell,i_cell_new,i_cell_1,i_cell_2
    character*100 :: info_str
    integer :: info

    real*8:: dist_vec(3), dist_atom_vec(3), dist, min_dist
    real*8 :: dir_cell(3), dir_cell_norm(3), dir_cell_PBC(3), dist_cell, inv_lattice_vector(3,3) 
    real*8 :: temp(3)
 
    integer,    dimension(:,:),allocatable:: cell_index_tmp
    integer :: i_atom_1, i_atom_2

    character(*), parameter :: func = 'initialize_position_in_hamiltonian_PBC'
    integer:: io,i_basis

!------------------------(1) calculate n_cells_in_hamiltonian------------------
    allocate(include_cell(n_cells),stat=info)
    call check_allocation(info, 'include_cell', func)

    allocate(new_index(n_cells),stat=info)
    call check_allocation(info, 'new_index', func)

    include_cell = .false.

    do i_cell = 1,n_cells,1
         dist_vec = cell_index(i_cell,1) * lattice_vector(1:3,1) &
                  + cell_index(i_cell,2) * lattice_vector(1:3,2) &
                  + cell_index(i_cell,3) * lattice_vector(1:3,3) 

         min_dist = 100000.d0
         do i_atom_1 = 1,n_atoms
            do i_atom_2 = 1,n_atoms

               dist_atom_vec = coords(1:3, i_atom_1) - coords(1:3,i_atom_2) + dist_vec
               dist = dot_product(dist_atom_vec, dist_atom_vec)
               min_dist = min( min_dist, dist)

            end do
         end do

         dist = sqrt( min_dist )

         if(dist < 2* maxval(atom_radius))then

            include_cell(i_cell) = .true.

         end if

    end do

    ! new index

    i_cell_new = 0
    new_index = 0
  
    do i_cell = 1,n_cells,1
       if(include_cell(i_cell))then

          i_cell_new = i_cell_new + 1               
          new_index(i_cell) = i_cell_new  

       end if
    end do

    n_cells_in_hamiltonian = i_cell_new + 1

    write(info_str,*) " |-----begin in position_in_hamiltonian_PBC---"
    call localorb_info(info_str,use_unit,'(A)')
 
   write(info_str,*) " | Number of super-cells in hamiltonian [n_cells_in_hamiltonian]:", n_cells_in_hamiltonian
    call localorb_info(info_str,use_unit,'(A)')

    write(info_str,*) " |-----end in position_in_hamiltonian_PBC---"
    call localorb_info(info_str,use_unit,'(A)')



!-----------------------------(2) calculate cell_index_hamiltonian and cell_index_to_i_cell----------------- 

    allocate(cell_index_hamiltonian(n_cells_in_hamiltonian,3),stat=info)
    call check_allocation(info, 'cell_index_hamiltonian',func)

    allocate(cell_index_to_i_cell(-number_of_super_cells(1):number_of_super_cells(1), &  
                                  -number_of_super_cells(2):number_of_super_cells(2), & 
                                  -number_of_super_cells(3):number_of_super_cells(3)) & 
                                  ,stat=info)
    call check_allocation(info, 'cell_index_to_i_cell', func)

    cell_index_to_i_cell=0
    do i_cell = 1,n_cells,1
       if(new_index(i_cell)>0)then
          cell_index_hamiltonian(new_index(i_cell),1:3) = cell_index(i_cell,1:3) 

          cell_index_to_i_cell(   &  
          cell_index(i_cell,1),cell_index(i_cell,2),cell_index(i_cell,3))= & 
          new_index(i_cell)

       end if
    end do
    
 
    n_cells_in_hamiltonian_array(1)=maxval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,1)) - & 
                                    minval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,1)) + 1
    n_cells_in_hamiltonian_array(2)=maxval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,2)) - & 
                                    minval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,2)) + 1
    n_cells_in_hamiltonian_array(3)=maxval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,3)) - & 
                                    minval(cell_index_hamiltonian(1:n_cells_in_hamiltonian-1,3)) + 1


!------------------------------(3) calculate position_in_hamiltonian_PBC---------------------
   !! Note about lattice_vector:  
   !! lattice_vector = ( a1 a2 a3) 
   !! 
   !! inv_lattice_vector = ( b1
   !!                        b2 
   !!                        b3 ) 
   !!
   !! recip_lattic_vector = 2pi ( b1 b2 b3 ) = trans( inv_lattice_vector)                      
   !!
   !! so (1) lattice_vector     * frac_coords = real_coords
   !!    (2) inv_lattice_vector * real_coords = frac_coords
   !! call get_map_to_center_cell_matrix(n_periodic,   &
   !!                                    lattice_vector,           &
   !!                                    inv_lattice_vector)
   !!  inv_lattice_vector = map_to_center_cell_matrix

    !------------(3.1) R1-R2 ----------------------------------------
    allocate(position_in_hamiltonian_PBC(n_cells_in_hamiltonian,n_cells_in_hamiltonian),stat=info)
    call check_allocation(info, 'position_in_hamiltonian_PBC', func)

   position_in_hamiltonian_PBC = n_cells_in_hamiltonian

!---------------------------------------------------------------------------
!note: this code is right for 1d, but wrong for 3d
! 
!wrong_code_for_3d  do i_cell_1 = 1,n_cells_in_hamiltonian-1,1
!wrong_code_for_3d  do i_cell_2 = 1,n_cells_in_hamiltonian-1,1
!wrong_code_for_3d
!wrong_code_for_3d       delta_cell(1:3) =  cell_index_hamiltonian(i_cell_1,1:3) -  &  ! <-----minus 
!wrong_code_for_3d                          cell_index_hamiltonian(i_cell_2,1:3)
!wrong_code_for_3d
!wrong_code_for_3d       dir_cell = matmul(lattice_vector, delta_cell)
!wrong_code_for_3d
!wrong_code_for_3d       dist_cell = dsqrt(dot_product(dir_cell, dir_cell))
!wrong_code_for_3d
!wrong_code_for_3d       dir_cell_norm(1:3) = dir_cell(1:3)/dist_cell
!wrong_code_for_3d
!wrong_code_for_3d    if(dist_cell.lt.2* maxval(atom_radius)) then 
!wrong_code_for_3d       cell_index_PBC(1:3)  = delta_cell(1:3)  
!wrong_code_for_3d    else 
!wrong_code_for_3d       dir_cell_PBC(1:3) = (dist_cell - 2*2* maxval(atom_radius))*dir_cell_norm(1:3)      
!wrong_code_for_3d       
!wrong_code_for_3d       temp = matmul(inv_lattice_vector, dir_cell_PBC)
!wrong_code_for_3d       write(use_unit,*) ' ' 
!wrong_code_for_3d       write(use_unit,*) '=========================' 
!wrong_code_for_3d       write(use_unit,*) 'i_cell_1',cell_index_hamiltonian(i_cell_1,1:3)
!wrong_code_for_3d       write(use_unit,*) 'i_cell_2',cell_index_hamiltonian(i_cell_2,1:3)
!wrong_code_for_3d       write(use_unit,*) 'temp_cell_index_PBC:', temp(1:3) 
!wrong_code_for_3d       write(use_unit,*) '=========================' 
!wrong_code_for_3d       do i_coord = 1,3 
!wrong_code_for_3d       if(temp(i_coord).ge.0.0d0) then 
!wrong_code_for_3d         cell_index_PBC(i_coord) = ceiling(temp(i_coord)) 
!wrong_code_for_3d       else
!wrong_code_for_3d         cell_index_PBC(i_coord) = floor(temp(i_coord))
!wrong_code_for_3d       endif
!wrong_code_for_3d       enddo
!wrong_code_for_3d
!wrong_code_for_3d    endif       
!wrong_code_for_3d    write(use_unit,*) 'radius:', 2* maxval(atom_radius)/lattice_vector(1,1)  
!wrong_code_for_3d    write(use_unit,*) 'i_cell_1:',i_cell_1, cell_index_hamiltonian(i_cell_1,1:3)
!wrong_code_for_3d    write(use_unit,*) 'i_cell_2:',i_cell_2, cell_index_hamiltonian(i_cell_2,1:3) 
!wrong_code_for_3d    write(use_unit,*) 'delta_cell:   ',delta_cell(1:3)
!wrong_code_for_3d    write(use_unit,*) 'cell_index_PBC', cell_index_PBC(1:3)
!wrong_code_for_3d    position_in_hamiltonian_PBC(i_cell_1, i_cell_2) = & 
!wrong_code_for_3d    cell_index_to_i_cell(cell_index_PBC(1),cell_index_PBC(2),cell_index_PBC(3))     
!wrong_code_for_3d    write(use_unit,*) 'position_in_hamiltonian_PBC:', position_in_hamiltonian_PBC(i_cell_1,i_cell_2) 
!---------------------------------------------------------------------------


  
!---------------------------------------------------------------------------------------
!note: this code is right for 1d, will also be used for 3d after I change other code
   do i_cell_1 = 1,n_cells_in_hamiltonian-1,1
   do i_cell_2 = 1,n_cells_in_hamiltonian-1,1
 
      do i_coord = 1, 3, 1

         delta_cell(i_coord) = cell_index_hamiltonian(i_cell_1,i_coord) -  &  ! <-----minus 
                               cell_index_hamiltonian(i_cell_2,i_coord)
       if(mod(n_cells_in_hamiltonian_array(i_coord),2).ne.1) then 
          call aims_stop('n_cells_in_hamiltonian_array is not odd, please check', func)
       endif 
 
       if(abs(delta_cell(i_coord)).le.(n_cells_in_hamiltonian_array(i_coord)-1)/2) then
         cell_index_PBC(i_coord) = delta_cell(i_coord)
       else 
         cell_index_PBC(i_coord) = &  
         - mod( n_cells_in_hamiltonian_array(i_coord), abs(delta_cell(i_coord)) )* &
         delta_cell(i_coord)/abs(delta_cell(i_coord))
       endif
      enddo  
 
       position_in_hamiltonian_PBC(i_cell_1, i_cell_2) = & 
       cell_index_to_i_cell(cell_index_PBC(1),cell_index_PBC(2),cell_index_PBC(3))     

    end do
    end do
!---------------------------------------------------------------------------------------

    !------------(3.2) R1+R2 ----------------------------------------
    allocate(position_in_hamiltonian_PBC_add(n_cells_in_hamiltonian,n_cells_in_hamiltonian),stat=info)
    call check_allocation(info, 'position_in_hamiltonian_PBC_add', func)

    position_in_hamiltonian_PBC_add = n_cells_in_hamiltonian
    do i_cell_1 = 1,n_cells_in_hamiltonian-1,1
    do i_cell_2 = 1,n_cells_in_hamiltonian-1,1
 
      do i_coord = 1, 3, 1
       delta_cell(i_coord) = cell_index_hamiltonian(i_cell_1,i_coord) +  & ! <------ add  
                             cell_index_hamiltonian(i_cell_2,i_coord)
 
       if(mod(n_cells_in_hamiltonian_array(i_coord),2).ne.1) then 
          call aims_stop('n_cells_in_hamiltonian_array is not odd, please check', func)
       endif 
 
       if(abs(delta_cell(i_coord)).le.(n_cells_in_hamiltonian_array(i_coord)-1)/2) then
         cell_index_PBC(i_coord) = delta_cell(i_coord)
       else 
         cell_index_PBC(i_coord) = &  
         - mod( n_cells_in_hamiltonian_array(i_coord), abs(delta_cell(i_coord)) )* &
         delta_cell(i_coord)/abs(delta_cell(i_coord))
       endif
 
      enddo  
 
       position_in_hamiltonian_PBC_add(i_cell_1, i_cell_2) = & 
       cell_index_to_i_cell(cell_index_PBC(1),cell_index_PBC(2),cell_index_PBC(3))     

    end do
    end do


!------------------------------(4) calculate basis_in_supercell---------------------
    allocate(cell_and_basis_to_basis_supercell(n_cells_in_hamiltonian,n_basis),stat=info)

    cell_and_basis_to_basis_supercell = 0
    io = 0
    do i_cell = 1,n_cells_in_hamiltonian-1,1
    do i_basis = 1,n_basis 
       io = io + 1
       cell_and_basis_to_basis_supercell(i_cell,i_basis) = io 
    enddo 
    enddo
    n_basis_supercell = io 

    deallocate(include_cell)
    deallocate(new_index)
    deallocate(cell_index_hamiltonian)
    deallocate(cell_index_to_i_cell)

  end subroutine initialize_position_in_hamiltonian_PBC


  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/initialize_position_in_hamiltonian
  !  NAME
  !    initialize_position_in_hamiltonian
  !  SYNOPSIS

  subroutine initialize_position_in_hamiltonian

    !  PURPOSE
    !    Initializes the packed matrix structure for periodic systems
    !    Packed matrixes are hamiltonian, ovl_matrix and density matrix
    use dimensions
    use runtime_choices
    use aims_memory_tracking, only: aims_allocate
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use localorb_io, only: localorb_info, use_unit
    use synchronize_mpi_basic, only: sync_integer, sync_int_vector, &
        sync_integer_vector
    use basis, only: atom_radius, outer_radius, basis_fn
    use geometry, only: lattice_vector, coords
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE






    logical,dimension(:), allocatable:: include_cell
    logical,dimension(:), allocatable:: include_cell2
    logical,dimension(:), allocatable:: include_cell_hartree
    integer,dimension(:), allocatable:: new_index, new_center_to_cell
    integer,dimension(:), allocatable:: column_index_temp
    integer,dimension(:), allocatable:: column_index_no_symmetry_temp

    integer:: i_cell_1,i_cell_2,i_cell_3, i_center,i_basis,i_cell_new, i_cell_L1, i_cell_L2, i_k_point
    character*100 :: info_str
    integer:: n_cells_new, n_cells_new_basis
    real*8:: dist_vec(3), dist_atom_vec(3), dist, min_dist
    integer,    dimension(:,:),allocatable:: cell_index_tmp

    integer :: i_basis_1, i_basis_2,info
    integer :: number_of_non_zeros
    integer :: i_atom_1, i_atom_2
    integer :: my_index, nnz_count
    character(*), parameter :: func = 'initialize_position_in_hamiltonian'

   !------------shanghui add for DFPT_phonon----------------------- 
    integer,dimension(:), allocatable:: inv_cells_in_hamiltonian
   !------------shanghui end add for DFPT_phonon----------------------- 


    allocate(new_index(n_cells),stat=i_center)
    call check_allocation(i_center, 'new_index                     ')

    allocate(include_cell(n_cells),stat=i_center)
    call check_allocation(i_center, 'include_cell                  ')

    allocate(include_cell2(n_cells),stat=i_center)
    call check_allocation(i_center, 'include_cell2                 ')

    allocate(include_cell_hartree(n_cells),stat=i_center)
    call check_allocation(i_center, 'include_cell_hartree          ')

    allocate(new_center_to_cell(n_centers),stat=i_center)
    call check_allocation(i_center, 'new_center_to_cell            ')

    allocate( cell_index_tmp(n_cells,3),stat=i_center)
    call check_allocation(i_center, 'cell_index_tmp                ')

   !-----shanghui add for DFPT_phonon-----------------------------------------
    allocate(inv_cells_in_hamiltonian(n_cells),stat=i_center)
    call check_allocation(i_center, 'inv_cells_in_hamiltonian                     ')
   !-----shanghui end add for DFPT_phonon-----------------------------------------


    include_cell = .false.


    do i_basis = 1, n_centers_basis_I

       include_cell(center_to_cell(Cbasis_to_center( i_basis))) =.true.

    end do


    i_cell_new = 0
    do i_cell_1 = 1,n_cells,1
       if(include_cell(i_cell_1))then

          i_cell_new = i_cell_new + 1
          new_index(i_cell_new) =  i_cell_1

       end if
    end do


    n_cells_new_basis = i_cell_new

    include_cell2 =.false.



    do i_cell_L1 = 1,n_cells_new_basis,1

       i_cell_1 = new_index(i_cell_L1)


       do i_cell_L2 = 1,n_cells_new_basis,1

          i_cell_2 = new_index(i_cell_L2)


          do i_cell_3 = 1,n_cells,1


             if( cell_index(i_cell_1,1) - cell_index(i_cell_2,1) == cell_index(i_cell_3,1) .and. &
                 cell_index(i_cell_1,2) - cell_index(i_cell_2,2) == cell_index(i_cell_3,2) .and. &
                 cell_index(i_cell_1,3) - cell_index(i_cell_2,3) == cell_index(i_cell_3,3)) then


                dist_vec = (cell_index(i_cell_1,1) - cell_index(i_cell_2,1) ) * lattice_vector(1:3,1) &
                         + (cell_index(i_cell_1,2) - cell_index(i_cell_2,2) ) * lattice_vector(1:3,2) &
                         + (cell_index(i_cell_1,3) - cell_index(i_cell_2,3) ) * lattice_vector(1:3,3) 


                min_dist = 100000.d0
                do i_atom_1 = 1,n_atoms
                   do i_atom_2 = 1,n_atoms

                      dist_atom_vec = coords(1:3, i_atom_1) - coords(1:3,i_atom_2) + dist_vec
                      dist = dot_product(dist_atom_vec, dist_atom_vec)
                      min_dist = min( min_dist, dist)

                   end do
                end do

                if(min_dist < ( 2* maxval(atom_radius))**2 )then

                   include_cell(i_cell_3) = .true.
                   include_cell2(i_cell_3) = .true.

                end if
                exit

             end if
          end do
       end do
    end do ! i_cell_L1


    ! new index


    i_cell_new = 0
    new_index = 0
    inv_cells_in_hamiltonian = 0 
  
    do i_cell_1 = 1,n_cells,1
       if(include_cell2(i_cell_1))then

          i_cell_new = i_cell_new + 1               
          new_index(i_cell_1) = i_cell_new  

         !-----shanghui add for check-----------------------------------------
          !new_index() is  cells_in_hamiltonian(i_cell_1)=i_cell_new
          inv_cells_in_hamiltonian(i_cell_new) = i_cell_1
         !-----shanghui end add for check-------------------------------------
       end if
    end do

    n_cells_in_hamiltonian = i_cell_new + 1

    do i_cell_1 = 1,n_cells,1
       if(include_cell(i_cell_1) .and. .not. include_cell2(i_cell_1) )then

          i_cell_new = i_cell_new + 1         
          new_index(i_cell_1) = i_cell_new

       end if
    end do

    n_cells_new = i_cell_new + 1

    cell_index_tmp(:,:) = cell_index(:,:)

    deallocate(cell_index)
    allocate(cell_index(n_cells_new,3),stat=info)
    call check_allocation(info, 'cell_index', func)


    do i_cell_1 = 1,n_cells,1
       if(new_index(i_cell_1)>0)then
          cell_index(new_index(i_cell_1),1:3) = cell_index_tmp(i_cell_1,1:3) 
       end if
    end do
    
    ! The last cell is not set in the above loop since it is a dummy cell
    cell_index(n_cells_new,1:3) = 999999999 ! force an out-of-bound error if this should ever be used!
 

    write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells_new*n_k_points*16)/2**10, ' KiB for k_phase.'
    call localorb_info(info_str,use_unit,'(A)')
    allocate(k_phase(n_cells_new,n_k_points),stat=i_center)
    call check_allocation(i_center, 'k_phase')

    do i_cell_1 = 1, n_cells_new-1
       do i_k_point = 1, n_k_points
          k_phase(i_cell_1, i_k_point) &
          & = product(k_phase_base(:, i_k_point)**cell_index(i_cell_1,:))
       enddo
    end do
    ! The last cell is not set in the above loops since it is a dummy cell
    k_phase(n_cells_new,:) = 0

    new_center_to_cell = 0
    do i_center = 1, n_centers,1
       new_center_to_cell(i_center)  = new_index(center_to_cell(i_center))
    end do

    center_to_cell =  new_center_to_cell

 
  !shanghui add note: center_to_cell and k_phase are changed after subroutine initialize_position_in_hamiltonian.
  ! center_to_cell(n_centers)                  ---> n_cells_in_hamiltonian
  ! k_phase(n_cells_in_hamiltonian,n_k_points) ---> exp(ik*R)



   !---------shanghui add for DFPT_phonon------------
    allocate(cell_and_atom_to_center(n_cells_new,n_atoms),stat=info)
    call check_allocation(info, 'cell_and_atom_to_center', func)

    do i_center = 1, n_centers,1
       if(center_to_cell(i_center).ne.0) then
       cell_and_atom_to_center(center_to_cell(i_center),center_to_atom(i_center))=i_center 
       endif
    enddo
   !---------shanghui end add for DFPT_phonon--------

    allocate(position_in_hamiltonian(n_cells_new,n_cells_new),stat=info)
    call check_allocation(info, 'position_in_hamiltonian', func)

    position_in_hamiltonian = n_cells_in_hamiltonian


    do i_cell_1 = 1,n_cells_new-1,1

       do i_cell_2 = 1,n_cells_new-1,1

          do i_cell_3 = 1,n_cells_in_hamiltonian-1,1


             if( cell_index(i_cell_1,1) - cell_index(i_cell_2,1) == cell_index(i_cell_3,1) .and. &
                 cell_index(i_cell_1,2) - cell_index(i_cell_2,2) == cell_index(i_cell_3,2) .and. &
                 cell_index(i_cell_1,3) - cell_index(i_cell_2,3) == cell_index(i_cell_3,3)) then


                position_in_hamiltonian(i_cell_1, i_cell_2) = i_cell_3


                exit
             end if

          end do
       end do
    end do



    if(transport_lead_calculation)then
    
       i_cell_2 = 0

       do i_cell_3 = 1,n_cells_in_hamiltonian-1,1

          i_cell_2 = max( i_cell_2, abs(cell_index(i_cell_3,3)))

       end do
          
       if(i_cell_2 > 1)then
             
          write(use_unit,*) 'Error: Too narrow lead region, overlap: ',i_cell_2
          stop
       end if
       if(i_cell_2 == 0)then
             
          write(use_unit,*) 'Error: The lead do not overlap at all: ',i_cell_2
          stop
       end if
    end if

    !--------------------------3----------------------------------------------

    select case(packed_matrix_format)


    case( PM_index) !---------------------------------------------------------------------------------------


       ! Memory allocations

       allocate(index_hamiltonian(2, n_cells_in_hamiltonian, n_basis),stat=info) 
       call check_allocation(info, 'index_hamiltonian', func)
       index_hamiltonian_dim2 = n_cells_in_hamiltonian

       ! column_index is distributed over all tasks, allocate only the needed length
       allocate(column_index(n_cells_in_hamiltonian-1,(n_basis-1)/n_tasks+1),stat=info)
       call check_allocation(info, 'column_index', func)

       column_index(:,:)%size = 0

       ! auxiliary array to store the column indices
       allocate(column_index_temp(n_basis),stat=info)
       call check_allocation(info, 'column_index_temp', func)

       ! Store the global column index in the distributed array column_index
       ! We avoid to store it as a whole on 1 processor because it is rather big

       i_cell_2 = 1
       nnz_count = 0 ! counts total nonzero entries

       do i_cell_1 = 1,n_cells_in_hamiltonian-1

          dist_vec = (cell_index(i_cell_1,1) - cell_index(i_cell_2,1) ) * lattice_vector(1:3,1) &
               + (cell_index(i_cell_1,2) - cell_index(i_cell_2,2) ) * lattice_vector(1:3,2) &
               + (cell_index(i_cell_1,3) - cell_index(i_cell_2,3) ) * lattice_vector(1:3,3)

          do i_basis_1 = 1,n_basis

             if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! distribute work and storage over all tasks

             number_of_non_zeros =  0

             do i_basis_2 = 1,i_basis_1

                dist_atom_vec = coords(1:3, Cbasis_to_center(i_basis_1)) - coords(1:3,Cbasis_to_center(i_basis_2)) + dist_vec(1:3)
                dist = dot_product(dist_atom_vec, dist_atom_vec)

                if(dist < ( outer_radius(basis_fn(i_basis_1)) + outer_radius(basis_fn(i_basis_2)) )**2 )then
                   number_of_non_zeros =  number_of_non_zeros + 1
                   column_index_temp(number_of_non_zeros) = i_basis_2
                end if

             end do

             my_index = (i_basis_1-1)/n_tasks + 1

             column_index(i_cell_1,my_index)%size = number_of_non_zeros
             if(number_of_non_zeros>0) then
                allocate(column_index(i_cell_1,my_index)%arr(number_of_non_zeros),stat=info)
                call check_allocation(info, 'column_index', func)   
                column_index(i_cell_1,my_index)%arr = column_index_temp(1:number_of_non_zeros)
             endif

             nnz_count = nnz_count + number_of_non_zeros

          end do

       end do

       deallocate(column_index_temp)

       if(use_local_index) then

          ! The index arrays can only be set after the grid is partitioned
          ! This is done in set_index_hamiltonian()
          ! We just get the size of the total hamiltonian (only for output)
          ! and do some safety settings here (although not necessary)

          call sync_integer(nnz_count)
          n_hamiltonian_matrix_size = nnz_count + 1
          index_hamiltonian(1, :, :) =  0
          index_hamiltonian(2, :, :) = -1

       else

          ! Gather the size of every column to be set in index_hamiltonian(2,.,.)

          index_hamiltonian = 0

          do i_cell_1 = 1,n_cells_in_hamiltonian-1
             do i_basis_1 = 1,n_basis

                if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

                my_index = (i_basis_1-1)/n_tasks + 1
                index_hamiltonian(2,i_cell_1,i_basis_1) = column_index(i_cell_1,my_index)%size

             end do
          end do

          call sync_int_vector(index_hamiltonian, 2*n_cells_in_hamiltonian*n_basis)

          ! Set index_hamiltonian

          number_of_non_zeros = 0
          do i_cell_1 = 1,n_cells_in_hamiltonian-1
             do i_basis_1 = 1,n_basis
                if(index_hamiltonian(2,i_cell_1,i_basis_1)>0) then
                   index_hamiltonian(1,i_cell_1,i_basis_1) = number_of_non_zeros+1
                   number_of_non_zeros = number_of_non_zeros + index_hamiltonian(2,i_cell_1,i_basis_1)
                   index_hamiltonian(2,i_cell_1,i_basis_1) = number_of_non_zeros
                else
                   ! empty column
                   index_hamiltonian(1,i_cell_1,i_basis_1) = 0
                   index_hamiltonian(2,i_cell_1,i_basis_1) =-1 
                endif
             end do
          end do
          index_hamiltonian(1,n_cells_in_hamiltonian,:) = 0
          index_hamiltonian(2,n_cells_in_hamiltonian,:) =-1

          ! Set column_index_hamiltonian

          call aims_allocate( column_index_hamiltonian, number_of_non_zeros+1, "column_index_hamiltonian" )

          column_index_hamiltonian = 0

          !i_want =0
          do i_cell_1 = 1,n_cells_in_hamiltonian-1
             do i_basis_1 = 1,n_basis

                if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

                if(index_hamiltonian(1,i_cell_1,i_basis_1)>0) then
                   my_index = (i_basis_1-1)/n_tasks + 1
                   column_index_hamiltonian(index_hamiltonian(1,i_cell_1,i_basis_1):    &
                        index_hamiltonian(2,i_cell_1,i_basis_1)) =  &
                        column_index(i_cell_1,my_index)%arr
                endif
             end do
          end do

          call sync_integer_vector(column_index_hamiltonian,number_of_non_zeros)

          n_trash_start =  number_of_non_zeros + 1

          n_hamiltonian_matrix_size = number_of_non_zeros + 1 ! RJ: why "+ 1" ???

       endif ! use_local_index

    end select ! packed_matrix_format-----------------------------------------------------------------------


    write(info_str,*) " | Number of super-cells (origin) [n_cells]                     :", n_cells
    call localorb_info(info_str,use_unit,'(A)')

    n_cells = n_cells_new

    write(info_str,*) " | Number of super-cells (after PM_index) [n_cells]             :", n_cells
    call localorb_info(info_str,use_unit,'(A)')

    write(info_str,*) " | Number of super-cells in hamiltonian [n_cells_in_hamiltonian]:", n_cells_in_hamiltonian
    call localorb_info(info_str,use_unit,'(A)')

    select case(packed_matrix_format)

    case(PM_index)

       write(info_str,*) ' | Size of matrix packed + index [n_hamiltonian_matrix_size] :', n_hamiltonian_matrix_size
       call localorb_info(info_str,use_unit,'(A)')
       sparse_matrix_needs_reduce =.true.

    end select


    if(use_DFPT_dielectric.or.use_DFPT_phonon_reduce_memory) then
!-------------(4) For index_hamiltonian_no_symmetry-------------------
       ! Memory allocations

       if (allocated(index_hamiltonian_no_symmetry)) deallocate(index_hamiltonian_no_symmetry)
       allocate(index_hamiltonian_no_symmetry(2, n_cells_in_hamiltonian, n_basis),stat=info) 
       call check_allocation(info, 'index_hamiltonian_no_symmetry', func)

       ! column_index is distributed over all tasks, allocate only the needed length
       if (allocated(column_index_no_symmetry)) deallocate(column_index_no_symmetry)
       allocate(column_index_no_symmetry(n_cells_in_hamiltonian-1,(n_basis-1)/n_tasks+1),stat=info)
       call check_allocation(info, 'column_index_no_symmetry', func)

       column_index_no_symmetry(:,:)%size = 0

       ! auxiliary array to store the column indices
       if (allocated(column_index_no_symmetry_temp)) deallocate(column_index_no_symmetry_temp)
       allocate(column_index_no_symmetry_temp(n_basis),stat=info)
       call check_allocation(info, 'column_index_no_symmetry_temp', func)

       ! Store the global column index in the distributed array column_index
       ! We avoid to store it as a whole on 1 processor because it is rather big

       i_cell_2 = 1
       nnz_count = 0 ! counts total nonzero entries

       do i_cell_1 = 1,n_cells_in_hamiltonian-1

          dist_vec = (cell_index(i_cell_1,1) - cell_index(i_cell_2,1) ) * lattice_vector(1:3,1) &
               + (cell_index(i_cell_1,2) - cell_index(i_cell_2,2) ) * lattice_vector(1:3,2) &
               + (cell_index(i_cell_1,3) - cell_index(i_cell_2,3) ) * lattice_vector(1:3,3)
          

          do i_basis_1 = 1,n_basis

             if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! distribute work and storage over all tasks

             number_of_non_zeros =  0

             do i_basis_2 = 1, n_basis ! shanghui change here for no_symmetry

                dist_atom_vec = coords(1:3, Cbasis_to_center(i_basis_1)) - coords(1:3,Cbasis_to_center(i_basis_2)) + dist_vec(1:3)
                dist = dot_product(dist_atom_vec, dist_atom_vec)

                if(dist < ( outer_radius(basis_fn(i_basis_1)) + outer_radius(basis_fn(i_basis_2)) )**2 )then
                   number_of_non_zeros =  number_of_non_zeros + 1
                   column_index_no_symmetry_temp(number_of_non_zeros) = i_basis_2
                end if

             end do

             my_index = (i_basis_1-1)/n_tasks + 1

             column_index_no_symmetry(i_cell_1,my_index)%size = number_of_non_zeros

             if(number_of_non_zeros>0) then
                allocate(column_index_no_symmetry(i_cell_1,my_index)%arr(number_of_non_zeros),stat=info)
                call check_allocation(info, 'column_index_no_symmetry', func)   
                column_index_no_symmetry(i_cell_1,my_index)%arr = column_index_no_symmetry_temp(1:number_of_non_zeros)
             endif

             nnz_count = nnz_count + number_of_non_zeros

          end do

       end do

       deallocate(column_index_no_symmetry_temp)

          ! Gather the size of every column to be set in index_hamiltonian(2,.,.)

          index_hamiltonian_no_symmetry = 0

          do i_cell_1 = 1,n_cells_in_hamiltonian-1
             do i_basis_1 = 1,n_basis

                if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

                my_index = (i_basis_1-1)/n_tasks + 1
                index_hamiltonian_no_symmetry(2,i_cell_1,i_basis_1) = column_index_no_symmetry(i_cell_1,my_index)%size
             end do
          end do

          call sync_int_vector(index_hamiltonian_no_symmetry, 2*n_cells_in_hamiltonian*n_basis)

          ! Set index_hamiltonian

          number_of_non_zeros = 0
          do i_cell_1 = 1,n_cells_in_hamiltonian-1
             do i_basis_1 = 1,n_basis
                if(index_hamiltonian_no_symmetry(2,i_cell_1,i_basis_1)>0) then
                   index_hamiltonian_no_symmetry(1,i_cell_1,i_basis_1) = number_of_non_zeros+1
                   number_of_non_zeros = number_of_non_zeros + index_hamiltonian_no_symmetry(2,i_cell_1,i_basis_1)
                   index_hamiltonian_no_symmetry(2,i_cell_1,i_basis_1) = number_of_non_zeros
                else
                   ! empty column
                   index_hamiltonian_no_symmetry(1,i_cell_1,i_basis_1) = 0
                   index_hamiltonian_no_symmetry(2,i_cell_1,i_basis_1) =-1 
                endif
             end do
          end do
          index_hamiltonian_no_symmetry(1,n_cells_in_hamiltonian,:) = 0
          index_hamiltonian_no_symmetry(2,n_cells_in_hamiltonian,:) =-1

          ! Set column_index_hamiltonian

          if (allocated(column_index_hamiltonian_no_symmetry)) deallocate(column_index_hamiltonian_no_symmetry)
          allocate(column_index_hamiltonian_no_symmetry(number_of_non_zeros+1),stat=info)
          call check_allocation(info, 'column_index_hamiltonian_no_symmetry', func)

          column_index_hamiltonian_no_symmetry = 0

          do i_cell_1 = 1,n_cells_in_hamiltonian-1
             do i_basis_1 = 1,n_basis

                if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

                if(index_hamiltonian_no_symmetry(1,i_cell_1,i_basis_1)>0) then
                   my_index = (i_basis_1-1)/n_tasks + 1
                   column_index_hamiltonian_no_symmetry(index_hamiltonian_no_symmetry(1,i_cell_1,i_basis_1):    &
                        index_hamiltonian_no_symmetry(2,i_cell_1,i_basis_1)) =  &
                        column_index_no_symmetry(i_cell_1,my_index)%arr

                endif
             end do
          end do

          call sync_integer_vector(column_index_hamiltonian_no_symmetry,number_of_non_zeros)

          n_hamiltonian_matrix_size_no_symmetry = number_of_non_zeros + 1 ! RJ: why "+ 1" ???



    select case(packed_matrix_format)

    case(PM_index)

       write(info_str,*) ' | Size of matrix packed + index (no_symmetry):', n_hamiltonian_matrix_size_no_symmetry
       call localorb_info(info_str,use_unit,'(A)')
       sparse_matrix_needs_reduce =.true.

    end select
    endif ! use_DFPT_dielectric


    deallocate(cell_index_tmp)
    deallocate(new_center_to_cell)
    deallocate(include_cell_hartree)
    deallocate(include_cell2)
    deallocate(include_cell)
    deallocate(new_index)

  end subroutine initialize_position_in_hamiltonian
  !******
  !----------------------------------------------------------------------------
  !****s* pbc_lists/initialize_packed_matrix_cluster
  !  NAME
  !    initialize_packed_matrix_cluster
  !  SYNOPSIS

  subroutine initialize_packed_matrix_cluster

    !  PURPOSE
    !    Initializes sparse matrix packing of Hamiltonian/Overlap matrix for non-periodic systems
    use dimensions
    use runtime_choices
    use aims_memory_tracking, only: aims_allocate
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: n_tasks, myid, check_allocation, aims_stop
    use synchronize_mpi_basic, only: sync_integer, sync_int_vector, &
        sync_integer_vector
    use basis, only: basis_fn, outer_radius
    use geometry, only: coords
    implicit none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    character*100 :: info_str
    real*8:: dist_atom_vec(3), dist
    integer:: i_basis_1, i_basis_2,info
    integer :: number_of_non_zeros
    integer,dimension(:), allocatable:: column_index_temp
    integer :: my_index, nnz_count
    character(*), parameter :: func = 'initialize_packed_matrix_cluster'

    n_cells_in_hamiltonian = 2
    select case(packed_matrix_format)

    case( PM_index) !---------------------------------------------------------------------------------------

       ! Memory allocations

       allocate(index_hamiltonian(2, 1, n_basis),stat=info)
       call check_allocation(info, 'index_hamiltonian', func)
       index_hamiltonian_dim2 = 1

       ! column_index is distributed over all tasks, allocate only the needed length
       allocate(column_index(1,(n_basis-1)/n_tasks+1),stat=info)
       call check_allocation(info, 'column_index', func)
       column_index(:,:)%size = 0

       ! auxiliary array to store the column indices
       allocate(column_index_temp(n_basis),stat=info)
       call check_allocation(info, 'column_index_temp', func)


       ! Store the global column index in the distributed array column_index
       ! We avoid to store it as a whole on 1 processor because it is rather big

       nnz_count = 0 ! counts total nonzero entries

       do i_basis_1 = 1,n_basis,1

          if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! distribute work and storage over all tasks

          number_of_non_zeros =  0

          do i_basis_2 = 1,i_basis_1, 1

             dist_atom_vec = coords(1:3,Cbasis_to_center(i_basis_1)) - coords(1:3,Cbasis_to_center(i_basis_2))

             dist = sqrt(dot_product(dist_atom_vec, dist_atom_vec))

             if(dist <= (outer_radius(basis_fn(i_basis_1)) + outer_radius(basis_fn(i_basis_2)))  )then
                number_of_non_zeros =  number_of_non_zeros + 1
                column_index_temp(number_of_non_zeros) = i_basis_2
             end if

          end do

          my_index = (i_basis_1-1)/n_tasks + 1

          column_index(1,my_index)%size = number_of_non_zeros
          if(number_of_non_zeros>0) then
             allocate(column_index(1,my_index)%arr(number_of_non_zeros),stat=info)
             call check_allocation(info, 'column_index', func)

             column_index(1,my_index)%arr = column_index_temp(1:number_of_non_zeros)
          endif

          nnz_count = nnz_count + number_of_non_zeros

       end do

       deallocate(column_index_temp)

       if(use_local_index) then

          ! The index arrays can only be set after the grid is partitioned
          ! This is done in set_index_hamiltonian()
          ! We just get the size of the total hamiltonian (only for output)
          ! and do some safety settings here (although not necessary)

          call sync_integer(nnz_count)
          n_hamiltonian_matrix_size = nnz_count + 1
          index_hamiltonian(1, :, :) =  0
          index_hamiltonian(2, :, :) = -1

       else

          ! Gather the size of every column to be set in index_hamiltonian(2,.,.)

          index_hamiltonian = 0

          do i_basis_1 = 1,n_basis

             if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

             my_index = (i_basis_1-1)/n_tasks + 1
             index_hamiltonian(2,1,i_basis_1) = column_index(1,my_index)%size

          end do

          call sync_int_vector(index_hamiltonian, 2*1*n_basis)

          ! Set index_hamiltonian

          number_of_non_zeros = 0
          do i_basis_1 = 1,n_basis
             if(index_hamiltonian(2,1,i_basis_1)>0) then
                index_hamiltonian(1,1,i_basis_1) = number_of_non_zeros+1
                number_of_non_zeros = number_of_non_zeros + index_hamiltonian(2,1,i_basis_1)
                index_hamiltonian(2,1,i_basis_1) = number_of_non_zeros
             else
                ! empty column
                index_hamiltonian(1,1,i_basis_1) = 0
                index_hamiltonian(2,1,i_basis_1) =-1 
             endif
          end do

          ! Set column_index_hamiltonian

          call aims_allocate(column_index_hamiltonian, number_of_non_zeros+1, "column_index_hamiltonian" )
          column_index_hamiltonian = 0

          do i_basis_1 = 1,n_basis

             if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

             if(index_hamiltonian(1,1,i_basis_1)>0) then
                my_index = (i_basis_1-1)/n_tasks + 1
                column_index_hamiltonian(index_hamiltonian(1,1,i_basis_1):    &
                     index_hamiltonian(2,1,i_basis_1)) =  &
                     column_index(1,my_index)%arr
             endif
          end do

          ! number_of_non_zeros could be large, use packet-wise allreduces in sync_integer_vector
          call sync_integer_vector(column_index_hamiltonian,number_of_non_zeros)

          n_trash_start =  number_of_non_zeros + 1

          n_hamiltonian_matrix_size = number_of_non_zeros + 1 ! RJ: why "+ 1" ???

       endif ! use_local_index

    case default

       call localorb_info('Invalid packing type')
       call aims_stop

    end select ! packed_matrix_format-----------------------------------------------------------------------


    allocate(position_in_hamiltonian(1,1))
    position_in_hamiltonian = 1

    write(info_str,'(2X,A)') 'Hamiltonian matrix size:'
    call localorb_info(info_str,use_unit,'(A)')

    write(info_str,*) ' | Size of matrix non-packed:',&
       int(n_basis,8)*int((n_basis+1),8) / int(2,8)

    call localorb_info(info_str,use_unit,'(A)')


    select case(packed_matrix_format)

    case(PM_index)

       write(info_str,*) ' | Size of matrix packed:    ', n_hamiltonian_matrix_size
       call localorb_info(info_str,use_unit,'(A)')

       sparse_matrix_needs_reduce =.true.

    end select


  end subroutine initialize_packed_matrix_cluster
  !******
  !-----------------------------------------------------------------------------------------------------
  !****s* pbc_lists/set_index_hamiltonian
  !  NAME
  !     set_index_hamiltonian
  !  SYNOPSIS

  subroutine set_index_hamiltonian

    !  PURPOSE
    !     After the partitioning is known this can called to resize the overlap and hamiltonian
    use dimensions
    use aims_memory_tracking, only: aims_allocate
    use localorb_io, only: localorb_allinfo, use_unit, OL_low
    use mpi_tasks
    use synchronize_mpi_basic, only: sync_int_vector, sync_find_max
    use grids, only: batches
    implicit none
    !  ARGUMENTS

    !  none

    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE





    integer :: i_pass, i_cell_1, i_cell_2, i_cell, i, i_basis_1, idx, my_index, i_my_batch, info
    integer*8 :: number_of_non_zeros
    integer :: nnz_col, len_col, n_compute_a, root, mpierr
    integer, allocatable :: column_index_size(:,:), column_index_temp(:)
    integer, allocatable :: i_cell_pair(:,:)
    integer :: n_cell_pairs, n
    logical, allocatable :: need_basis(:,:), need_basis_ham(:)
    character*200 :: info_str
    character(*), parameter :: func = 'set_index_hamiltonian'

    ! Set up need_basis
  
    if(n_periodic>0) then
       allocate(need_basis(n_basis,n_cells),stat=info)
       call check_allocation(info, 'need_basis', func)

       allocate(i_cell_pair(2,n_cells*n_cells))
    else
       allocate(need_basis(n_basis,1),stat=info)
       call check_allocation(info, 'need_basis', func)

       allocate(i_cell_pair(2,1))
       n_cell_pairs = 1
       i_cell_pair(1,1) = 1
       i_cell_pair(2,1) = 1
    endif


    need_basis(:,:) = .false.

    do i_my_batch = 1, n_my_batches, 1
       n_compute_a = batches(i_my_batch)%batch_n_compute            
       do i = 1, n_compute_a
          i_basis_1 = batches(i_my_batch)%batch_i_basis(i)
          need_basis(Cbasis_to_basis(i_basis_1),center_to_cell(Cbasis_to_center(i_basis_1))) = .true.
       enddo
    enddo


    ! Distribute length of global column index

    allocate(column_index_size(n_cells_in_hamiltonian, n_basis),stat=info)
    call check_allocation(info, 'column_index_size', func)

    column_index_size = 0

    do i_cell = 1,n_cells_in_hamiltonian-1
       do i_basis_1 = 1,n_basis

          if(mod(i_basis_1-1,n_tasks) /= myid) cycle ! care only about my part

          my_index = (i_basis_1-1)/n_tasks + 1
          column_index_size(i_cell,i_basis_1) = column_index(i_cell,my_index)%size

       end do
    end do

    call sync_int_vector(column_index_size, n_cells_in_hamiltonian*n_basis)

    index_hamiltonian(1,:,:) = 0
    index_hamiltonian(2,:,:) =-1

    ! auxiliary array to store the column indices
    allocate(column_index_temp(n_basis),stat=info)
    call check_allocation(info, 'column_index_temp', func)

    allocate(need_basis_ham(n_basis),stat=info)
    call check_allocation(info, 'need_basis_ham                ')


    ! Broadcast column_index part by part, every processor builds up own_column_index

    do i_pass = 1,2 ! i_pass==1: count non zeros, i_pass==2: set column_index_hamiltonian

       number_of_non_zeros = 0
       do i_cell = 1,n_cells_in_hamiltonian-1

          if(n_periodic > 0) then
             ! Set up a list of cell pairs included in i_cell
             n_cell_pairs = 0
             do i_cell_1 = 1, n_cells
                do i_cell_2 = 1, n_cells
                   if(position_in_hamiltonian(i_cell_1, i_cell_2) == i_cell) then
                      n_cell_pairs = n_cell_pairs+1
                      i_cell_pair(1,n_cell_pairs) = i_cell_1
                      i_cell_pair(2,n_cell_pairs) = i_cell_2
                   endif
                enddo
             enddo
          endif

          do i_basis_1 = 1,n_basis

             len_col = column_index_size(i_cell,i_basis_1)
             if(len_col == 0) cycle

             root = mod(i_basis_1-1,n_tasks)
             if(myid==root) then
                my_index = (i_basis_1-1)/n_tasks + 1
                column_index_temp(1:len_col) = column_index(i_cell,my_index)%arr(1:len_col)
             endif

             call MPI_Bcast(column_index_temp,len_col,MPI_INTEGER,root,mpi_comm_global,mpierr)

             if(.not. any(need_basis(i_basis_1,:))) cycle

             ! Check which basis functions are actually needed for current cell

             need_basis_ham(:) = .false.
             do n = 1, n_cell_pairs
                if(need_basis(i_basis_1,i_cell_pair(1,n))) &
                   need_basis_ham(:) = need_basis_ham(:) .or. need_basis(:,i_cell_pair(2,n))
             enddo

             nnz_col = 0
             do idx=1,len_col
                if(need_basis_ham(column_index_temp(idx))) then
                   nnz_col = nnz_col+1
                   if(i_pass==2) &
                        column_index_hamiltonian(number_of_non_zeros+nnz_col) = column_index_temp(idx)
                endif
             enddo

             if(nnz_col>0 .and. i_pass==2) then
                index_hamiltonian(1,i_cell,i_basis_1) = number_of_non_zeros+1
                index_hamiltonian(2,i_cell,i_basis_1) = number_of_non_zeros+nnz_col
             endif
             number_of_non_zeros = number_of_non_zeros + nnz_col
          enddo
       enddo

       if(i_pass==1) then
          if ( number_of_non_zeros .le. huge(n_hamiltonian_matrix_size)) then
             n_hamiltonian_matrix_size = int(number_of_non_zeros)
          else
             call aims_stop('Integer overflow detected for n_hamiltonian_matrix_size. Stop here.', func)
          end if
          call sync_find_max(n_hamiltonian_matrix_size, max_matrix_size)
          if(myid==0) write(use_unit,*) ' use_local_index: Max. local matrix size = ',max_matrix_size
          call aims_allocate(column_index_hamiltonian, n_hamiltonian_matrix_size, "column_index_hamiltonian" )
          column_index_hamiltonian = 0
       endif
    enddo

    n_hamiltonian_matrix_size = number_of_non_zeros

    write(info_str, "(a,i5,a,i10)") &
    & '  Task: ',myid,' matrix size: ',n_hamiltonian_matrix_size
    call localorb_allinfo(info_str, use_unit, '(A)', OL_low)
    deallocate(column_index_size)
    deallocate(column_index_temp)
    deallocate(need_basis)
    deallocate(need_basis_ham)
    deallocate(i_cell_pair)

  end subroutine set_index_hamiltonian
    
  !******
  !-----------------------------------------------------------------------------------------------------
  !****s* pbc_lists/remove_small_numbers_in_hamiltonian_and_ovlp
  !  NAME
  !    remove_small_numbers_in_hamiltonian_and_ovlp
  !  SYNOPSIS

  subroutine remove_small_numbers_in_hamiltonian_and_ovlp(hamiltonian, overlap_matrix, info)

    !  PURPOSE
    !    If ovlp and hamiltonian matrix have values smaller than  packed_matrix_threshold
    !    then memory slot is removed from packed matrix.
    !  USES

    use dimensions
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use analyze_arrays
    use mpi_tasks, only: check_allocation
    use runtime_choices
    use localorb_io, only: localorb_info, use_unit, OL_norm
    implicit none

    !  ARGUMENTS

    real*8 :: hamiltonian(n_hamiltonian_matrix_size,n_spin)
    real*8 :: overlap_matrix(n_hamiltonian_matrix_size)
    logical :: info

    !  INPUTS
    ! o   hamiltonian    -- hamiltonian matrix (packed matrix form)
    ! o   overlap_matrix -- overlap matrix     (packed matrix form)
    !  OUTPUT
    ! o   hamiltonian    -- hamiltonian matrix (packed matrix form)
    ! o   overlap_matrix -- overlap matrix     (packed matrix form)
    ! o   info -- did anything change?
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE



    integer,dimension(:),allocatable:: temp_i
    integer:: n_size_of_hamiltonian_matrix_new, i_cell, i_index_new, i_index_old
    integer:: i_basis_2, i_begin_new, i_size, i_size_begin, i_size_end,iinfo
    character*100 :: info_str
    character(*), parameter :: func = 'remove_small_numbers_in_hamiltonian_and_ovlp'

    info = .false.

    if(sparse_matrix_needs_reduce)then

       write(info_str,'(2X,A)') 'Decreasing sparse matrix size:'
       call localorb_info(info_str,use_unit,'(A)',OL_norm)

       write(info_str,'(2X,A,E12.4)') '| Tolerance:', packed_matrix_threshold
       call localorb_info(info_str,use_unit,'(A)',OL_norm)

       write(info_str,'(2X,A)') 'Hamiltonian matrix'
       call localorb_info(info_str,use_unit,'(A)',OL_norm)
       call analyze_1d_array( hamiltonian, n_hamiltonian_matrix_size, packed_matrix_threshold, .true., i_size )

       write(info_str,'(2X,A)') 'Overlap matrix'
       call localorb_info(info_str,use_unit,'(A)',OL_norm)
       call analyze_1d_array( overlap_matrix, n_hamiltonian_matrix_size, packed_matrix_threshold, .true., i_size )


       i_index_new = 0

       do i_cell = 1,n_cells_in_hamiltonian-1

          do i_basis_2 = 1, n_basis

             if( index_hamiltonian(1,i_cell, i_basis_2) > 0 )then

                i_index_old = index_hamiltonian(1,i_cell, i_basis_2)-1

                i_begin_new = i_index_new + 1

                i_size_begin =  index_hamiltonian(1,i_cell, i_basis_2)
                i_size_end   =  index_hamiltonian(2,i_cell, i_basis_2)


                do i_size =i_size_begin,i_size_end, 1

                   i_index_old = i_index_old + 1

                   if(  maxval(abs(hamiltonian( i_index_old,1:n_spin))) >  packed_matrix_threshold  .or.  &
                        abs(overlap_matrix( i_index_old)) > packed_matrix_threshold .or. &
                        ( i_cell==1 .and. & 
                          (Cbasis_to_center(i_basis_2).eq.Cbasis_to_center(column_index_hamiltonian(i_index_old)))) ) then

                      i_index_new = i_index_new + 1

                      hamiltonian(i_index_new,1:n_spin) = hamiltonian( i_index_old,1:n_spin)
                      overlap_matrix( i_index_new) = overlap_matrix( i_index_old)
                      column_index_hamiltonian(i_index_new) =  column_index_hamiltonian(i_index_old)

                   end if
                end do

                if(i_begin_new > i_index_new)then
                   index_hamiltonian(1,i_cell, i_basis_2) = 0
                   index_hamiltonian(2,i_cell, i_basis_2) = -1
                else
                   index_hamiltonian(1,i_cell, i_basis_2) = i_begin_new
                   index_hamiltonian(2,i_cell, i_basis_2) = i_index_new
                end if


             else
                index_hamiltonian(1,i_cell, i_basis_2) = 0
                index_hamiltonian(2,i_cell, i_basis_2) = -1
             end if
          end do
       end do


       n_trash_start =  i_index_new + 1

       n_size_of_hamiltonian_matrix_new = i_index_new + 1

       allocate(temp_i(n_size_of_hamiltonian_matrix_new),stat=iinfo)
       call check_allocation(iinfo, 'temp_i                        ')

       temp_i = column_index_hamiltonian(1:n_size_of_hamiltonian_matrix_new)
       call aims_deallocate( column_index_hamiltonian,                                 "column_index_hamiltonian" )
       call aims_allocate( column_index_hamiltonian, n_size_of_hamiltonian_matrix_new, "column_index_hamiltonian" )

       column_index_hamiltonian = temp_i
       deallocate(temp_i)

       if(use_scalapack)then
          info = .true.
       elseif(n_size_of_hamiltonian_matrix_new < n_hamiltonian_matrix_size)then
          info = .true.
       else
          info = .false.
       end if



       !          allocate(temp_r(n_size_of_hamiltonian_matrix_new))
       !          temp_r = hamiltonian(1:n_size_of_hamiltonian_matrix_new)
       !          deallocate(hamiltonian)
       !          allocate(hamiltonian(n_size_of_hamiltonian_matrix_new))
       !          hamiltonian = temp_r

       !          temp_r = overlap_matrix(1:n_size_of_hamiltonian_matrix_new)
       !          deallocate(overlap_matrix)
       !          allocate(overlap_matrix(n_size_of_hamiltonian_matrix_new))
       !          overlap_matrix = temp_r



       !          temp_i = column_index_hamiltonian
       !          deallocate(column_index_hamiltonian)
       !          allocate(column_index_hamiltonian(n_size_of_hamiltonian_matrix_new))
       !          column_index_hamiltonian = temp_i(1:n_size_of_hamiltonian_matrix_new)

       n_hamiltonian_matrix_size = i_index_new + 1

       write(info_str,'(2X,A,I12)') 'New size of hamiltonian matrix:',n_hamiltonian_matrix_size 
       call localorb_info(info_str,use_unit,'(A)',OL_norm)

       !         deallocate(temp_i)
       !         deallocate(temp_r)

       sparse_matrix_needs_reduce = .false.


    end if


  end subroutine remove_small_numbers_in_hamiltonian_and_ovlp
    !******
  !-----------------------------------------------------------------------------------------------------

  !****f* pbc_lists/find_position_in_hamiltonian
  !  NAME
  !     find_position_in_hamiltonian
  !  SYNOPSIS

  function find_position_in_hamiltonian( i_basis_1, i_basis_2) result( position)

    !  PURPOSE
    !    Finds the position of the  i_basis_1, i_basis_2 in the packed matrix format.
    use dimensions
    implicit none
    !  ARGUMENTS

    integer :: i_basis_1, i_basis_2
    integer :: position

    !  INPUTS
    !   o i_basis_1, i_basis_2 -- index pair in non-packed matrix
    !  OUTPUT
    !   o  position -- position of the variable in a packed matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE



    integer :: i_up, i_down, i_center
    integer :: i_ite, i_cell_index
    integer :: i_bc_1, i_bc_2
    character(*), parameter :: func = 'find_position_in_hamiltonian'

    ! The Bi-sectioning method for finding data points in packed matrix format index.
    ! This routine is time critical. It is not also so much optimized for being fast.



    position = n_trash_start

    i_bc_1 = Cbasis_to_basis(i_basis_1)
    i_bc_2 = Cbasis_to_basis(i_basis_2)

    print *, "Basis 1: ", i_basis_1, " CenterTo: ", i_bc_1
    print *, "Basis 2: ", i_basis_2, " CenterTo: ", i_bc_2


    if(Cbasis_to_basis(i_basis_2) <= Cbasis_to_basis(i_basis_1))then

       !            write(use_unit,*) 'etsi'


       i_cell_index = position_in_hamiltonian( &
            center_to_cell(Cbasis_to_center(i_basis_1)) , center_to_cell(Cbasis_to_center(i_basis_2))) 
       print *, "CBasisCenter 1: ", Cbasis_to_center(i_basis_1), " CTC Basis 1:", center_to_cell(Cbasis_to_center(i_basis_1))
       print *, "CBasisCenter 2: ", Cbasis_to_center(i_basis_2), " CTC Basis 2:", center_to_cell(Cbasis_to_center(i_basis_2))
       print *, "Position in Hamiltonian: ", i_cell_index


       i_down =  index_hamiltonian(1,i_cell_index, i_bc_1)
       i_up   =  index_hamiltonian(2,i_cell_index, i_bc_1)
       i_center = (i_up+i_down)/2

       print *, "iDown: ", i_down
       print *, "iUp: ", i_up
       print *, "iCenter: ", i_center

       if(i_down >0)then




          loop: do i_ite = 1,n_basis+1

             !              write(use_unit,*)  column_index_hamiltonian(i_center) , i_bc_2
             !              write(use_unit,*) 'p', i_up, i_center, i_down


             print *, "coulumnIndexHamiltonian: ", &
                      column_index_hamiltonian(i_center), &
                     "ibc2: ", i_bc_2 
             if( column_index_hamiltonian(i_center) == i_bc_2) then

                position = i_center
                !                write(use_unit,*) 'loytyi', position
                exit loop

             else if(  column_index_hamiltonian(i_center) > i_bc_2) then

                i_up = i_center
                i_center = (i_up+i_down)/2

             else if(  column_index_hamiltonian(i_center) < i_bc_2)then

                i_down = i_center
                i_center = (i_up+i_down)/2

             else


                ! do i_center = index_hamiltonian(1,i_cell_index, i_bc_1), index_hamiltonian(2,i_cell_index, i_bc_1)
                !    if(column_index_hamiltonian(i_center) == i_bc_2)then
                !       write(use_unit,*) 'ERROR loytyi sittenkin'
                !       position = i_center
                !    end if
                ! end do



                exit loop

             end if


             if(abs(i_up -i_down)==1)then

                if( column_index_hamiltonian(i_down) == i_bc_2) then

                   position = i_down
                   exit loop

                else if( column_index_hamiltonian(i_up) == i_bc_2) then

                   position = i_up
                   exit loop

                else
                   exit loop
                end if
             end if
          end do loop
       end if
    end if


    print *, "Resulting Index: ", position

    !         write(use_unit,*) '-', position


  end function find_position_in_hamiltonian

  !******
  !-----------------------------------------------------------------------
  !****s* pbc_lists/deallocate_pbc_lists
  !  NAME
  !    deallocate_pbc_lists
  !  SYNOPSIS

  subroutine deallocate_pbc_lists()

    !  PURPOSE
    !    Deallocates pbc list variables
    use aims_memory_tracking, only : aims_deallocate
    use dimensions, only: n_states_k, ik2irred_map
    implicit none
    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE

    integer i_basis_1, i_cell_1
    character(*), parameter :: func = 'deallocate_pbc_lists'

    sparse_matrix_needs_reduce = .false.

    if ( allocated (coords_center) ) then
       deallocate ( coords_center )
    end if
    if ( allocated (species_center) ) then
       deallocate ( species_center )
    end if
    if ( allocated(center_to_atom))then
       deallocate(center_to_atom)
    end if
    if(allocated(center_to_cell))then
       deallocate(center_to_cell)
    end if
    if(allocated(center_to_cell_pairs))then !SVL
       deallocate(center_to_cell_pairs)
    end if
    if(allocated(basis_cell_to_all))then !SVL
       deallocate(basis_cell_to_all)
    end if
    if (allocated(k_weights)) then
       deallocate(k_weights)
    end if
    if (allocated(k_phase_base)) then
       deallocate(k_phase_base)
    end if
    if (allocated(k_phase)) then
       deallocate(k_phase)
    end if
    if (allocated(k_phase_exx)) then
       deallocate(k_phase_exx)
    end if
    if (allocated(k_point_list)) then
       deallocate(k_point_list)
    end if

    if (allocated(n_states_k)) then
       deallocate(n_states_k)
    end if

    if (allocated(ik2irred_map)) then
       deallocate(ik2irred_map)
    end if

    if (allocated(k_minus_q_point_list)) then
       deallocate(k_minus_q_point_list)
    end if

    if (allocated(irk_point_list)) then
       deallocate(irk_point_list)
    end if

    if (allocated(irk_weight)) then
       deallocate(irk_weight)
    end if

    if (allocated(irk_point_mapping)) then
       deallocate(irk_point_mapping)
    end if

    if (allocated(inv_irk_point_mapping)) then
       deallocate(inv_irk_point_mapping)
    end if

    if (allocated(irk_point_included)) then
       deallocate(irk_point_included)
    end if

    if(allocated( Cbasis_to_basis))then
       deallocate(Cbasis_to_basis)
    end if

    if(allocated( Cbasis_to_basis_s))then
       deallocate(Cbasis_to_basis_s)
    end if

    if(allocated( Cbasis_to_center))then
       deallocate(Cbasis_to_center)
    end if

    if(allocated( Cbasis_to_center_s))then
       deallocate(Cbasis_to_center_s)
    end if

    if(allocated(cell_index) ) then
       deallocate(cell_index)
    end if

    if(allocated( Cbasis_to_atom  ))then
       deallocate(Cbasis_to_atom)
    end if

    if(allocated( Cbasis_to_atom_s))then
       deallocate(Cbasis_to_atom_s)
    end if

    if(allocated( centers_hartree_potential  ))then
       deallocate(centers_hartree_potential  )
    end if

    if(allocated(centers_hartree_multipole   ))then
       deallocate(centers_hartree_multipole  )
    end if

    if(allocated(  centers_ele_summation ))then
       deallocate( centers_ele_summation )
    end if

    if(allocated( new_centers_basis_integrals  ))then
       deallocate( new_centers_basis_integrals )
    end if

    if(allocated( centers_basis_integrals  ))then
       deallocate( centers_basis_integrals )
    end if

    if(allocated( inv_centers_basis_integrals  ))then
       deallocate( inv_centers_basis_integrals )
    end if

    if(allocated(atom_radius_hartree))then
       deallocate(atom_radius_hartree)
    end if
    if(allocated(atom_radius_hartree_sq))then
       deallocate(atom_radius_hartree_sq)
    end if

    if(allocated(position_in_hamiltonian))then
       deallocate(position_in_hamiltonian)
    end if

 
   if( allocated(index_hamiltonian))       deallocate(index_hamiltonian)
    if(allocated(cell_index))               deallocate(cell_index)
    if(allocated(cell_index_pairs)) deallocate(cell_index_pairs) !SVL
    if(allocated(cell_index_bvk)) deallocate(cell_index_bvk) !SVL
    if(allocated(ovlp3fn_cells)) deallocate(ovlp3fn_cells) !SVL
    if(allocated(column_index_hamiltonian)) call aims_deallocate( column_index_hamiltonian, "column_index_hamiltonian" ) 

    if(allocated(column_index)) then
       do i_basis_1 = 1,ubound(column_index,2)
          do i_cell_1 = 1,ubound(column_index,1)
             if(column_index(i_cell_1,i_basis_1)%size > 0) &
                  deallocate(column_index(i_cell_1,i_basis_1)%arr)
          enddo
       enddo
       deallocate(column_index)
    endif

    ! BL: This deallocation statements were forgotton.
    ! Without the the code is not passing the regression tests due
    ! to double allocating the variables.
    ! @ shanghui: Please check that this is correct here and delete
    ! comment then.
    !---------shanghui add for DFPT_phonon------------------ 
    if(allocated(centers_in_hamiltonian)) then
       deallocate(centers_in_hamiltonian)
    end if
    
    if(allocated(inv_centers_in_hamiltonian)) then
       deallocate(inv_centers_in_hamiltonian)
    end if
    if(allocated(cell_and_atom_to_center)) then
       deallocate(cell_and_atom_to_center)
    end if
    if(allocated(position_in_hamiltonian_PBC))then
       deallocate(position_in_hamiltonian_PBC)
    end if
    if(allocated(position_in_hamiltonian_PBC_add))then
       deallocate(position_in_hamiltonian_PBC_add)
    end if
    if(allocated(cell_and_basis_to_basis_supercell)) then
       deallocate(cell_and_basis_to_basis_supercell)
    end if
    
    !----------begin for sc_DFPT-------------
    if(allocated(center_in_sc_DFPT_to_center)) then
       deallocate(center_in_sc_DFPT_to_center)
    end if
    
    if(allocated(center_to_center_in_sc_DFPT)) then
       deallocate(center_to_center_in_sc_DFPT)
    end if
    if(allocated(center_in_sc_DFPT_to_cell_in_sc_DFPT)) then
       deallocate(center_in_sc_DFPT_to_cell_in_sc_DFPT)
    end if
    if(allocated(center_in_sc_DFPT_to_atom)) then        
       deallocate(center_in_sc_DFPT_to_atom)
     endif

    if(allocated(cell_diff_sc_DFPT))then
       deallocate(cell_diff_sc_DFPT)
    end if
    if(allocated(cell_add_sc_DFPT))then
       deallocate(cell_add_sc_DFPT)
    end if
    if(allocated(cell_and_basis_to_basis_sc_DFPT)) then
       deallocate(cell_and_basis_to_basis_sc_DFPT)
    end if
    if(allocated(cell_and_atom_to_center_sc_DFPT)) then
       deallocate(cell_and_atom_to_center_sc_DFPT)
    end if
    if(allocated(cell_index_sc_DFPT)) then 
       deallocate(cell_index_sc_DFPT)
    endif
    if(allocated(cell_in_hamiltonian_to_cell_in_sc_DFPT)) then 
       deallocate(cell_in_hamiltonian_to_cell_in_sc_DFPT) 
    endif
    if(allocated(k_point_loc)) then
       deallocate(k_point_loc)
    end if
    !----------end for sc_DFPT-------------


   !---------shanghui end add for DFPT_phonon--------


  end subroutine deallocate_pbc_lists
  !******
  !------------------------------------------------------------------------------
  !****s* pbc_lists/kweight_occs
  !  NAME
  !    kweight_occs
  !  SYNOPSIS

  subroutine kweight_occs(caller, occ_numbers)

    !  PURPOSE
    !    Multiply the occupation numbers by the kweights.
    !  USES

    use dimensions, only: n_states, n_spin, n_k_points
    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: caller
    real*8, intent(INOUT)  :: occ_numbers(n_states, n_spin, n_k_points)

    !  INPUTS
    !    o caller -- calling procedure for error messages
    !    o occ_numbers -- unweighted occ_numbers
    !  OUTPUTS
    !    o occ_numbers -- weighted occ_numbers
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_k
    character(*), parameter :: func = 'kweight_occs'

    do i_k = 1, n_k_points
       occ_numbers(:,:, i_k) = occ_numbers(:,:,i_k) * k_weights(i_k)
    enddo
    call check_occs(caller, occ_numbers, .true.)

  end subroutine kweight_occs
  !******

  !------------------------------------------------------------------------------
  !****s* pbc_lists/de_kweight_occs
  !  NAME
  !    de_kweight_occs
  !  SYNOPSIS

  subroutine de_kweight_occs(caller, occ_numbers)

    !  PURPOSE
    !    Remove the kweights factor from the occupation numbers.
    !  USES

    use dimensions, only: n_states, n_spin, n_k_points
    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT)  :: occ_numbers(n_states, n_spin, n_k_points)
    character(*), intent(IN) :: caller

    !  INPUTS
    !    o caller -- calling procedure for error messages
    !    o occ_numbers -- weighted occ_numbers
    !  OUTPUTS
    !    o occ_numbers -- unweighted occ_numbers
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_k
    character(*), parameter :: func = 'de_kweight_occs'

    call check_occs(caller, occ_numbers, .true.)
    do i_k = 1, n_k_points
       occ_numbers(:,:, i_k) = occ_numbers(:,:,i_k) / k_weights(i_k)
    enddo

  end subroutine de_kweight_occs

  !******
  !------------------------------------------------------------------------------
  !******
 end module pbc_lists
