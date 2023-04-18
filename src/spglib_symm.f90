!****s* FHI-aims/spglib_symm.f90
!  NAME
!    spglib_symm.f90
!  SYNOPSIS


module spglib_symmetry

!  PURPOSE
  !
  !    Interfaces to spglib to obtain symmetry informations
  !
  !  Subroutines:
  !  o get_symmetry
  !  o write_symm_info
  !  o symmetrize_forces
  !  o ...
  !
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).

  !  USES

!  use geometry,        only: coords, lattice_vector, species, cart2frac
  use spglib    ! WPH: This is a wrapper module around spglib
  implicit none

  ! Variables

  integer, dimension(:),  allocatable :: spg_atom_type
  real*8, dimension(:,:), allocatable :: spg_lat_vec
  real*8, dimension(:,:), allocatable :: spg_coords
  integer, dimension(:,:,:), allocatable :: spg_rotations
  real*8, dimension(:,:),  allocatable :: spg_shift
  integer :: num_symmetry
  integer, dimension(:),  allocatable :: map
  integer, dimension(:),  allocatable :: map_sym
  integer, dimension(:),  allocatable :: map_inv
  type(SpglibDataset) :: dset
  integer, dimension(:), allocatable :: spg_partition_tab

  ! From spglib:
  !typedef struct {
  !int spacegroup_number;
  !int hall_number;
  !char international_symbol[11];
  !char hall_symbol[17];
  !char setting[6];
  !double transformation_matrix[3][3];
  !double origin_shift[3];
  !int n_operations;
  !int (*rotations)[3][3];
  !double (*translations)[3];
  !int n_atoms;
  !int *wyckoffs;
  !int *equivalent_atoms;
  !double brv_lattice[3][3];
  !int *brv_types;
  !double (*brv_positions)[3];
  !} SpglibDataset;

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

contains

subroutine check_spglib(flag_spg_lib)

  !  PURPOSE
  !
  !    Check if spglib has been compiled
  !

  implicit none

  !  ARGUMENTS
  !    flag_spg_lib : boolean here set to true since we are in the real module 
  !  INPUTS
     logical, intent(INOUT) :: flag_spg_lib
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  flag_spg_lib = .true.

end subroutine check_spglib


subroutine get_symmetry(lattice_vector, n_atoms, species, coords)

  !  PURPOSE
  !
  !    Interfaces to spglib to obtain symmetry informations
  !
  !  USES
  use applicable_citations, only: cite_reference
  use dimensions, only: sym_precision
  use geometry, only: cart2frac
  use mpi_tasks, only: check_allocation
  implicit none
  !  ARGUMENTS
  real*8 , intent(in) :: lattice_vector(3,3)
  integer, intent(in) :: n_atoms
  integer, intent(in) :: species(n_atoms)
  real*8,  intent(in) :: coords(3, n_atoms)
  !  INPUTS
  !    o lattice_vector - the lattice vectors
  !    o n_atoms - number of atoms
  !    o species - species for atoms
  !    o coords - Cartesian coordinates for atoms
  !  OUTPUTS
  !    o None (modifies module variables)
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables
  integer :: i_atom, i_type
  integer :: n_types
  integer :: info
  integer, dimension(n_atoms) :: spg_species_type
  logical :: found

  if (allocated(spg_coords)) deallocate(spg_coords)
  allocate (spg_coords(3,n_atoms) ,stat=info)
  call check_allocation(info, 'spg_coords                        ')

  if (allocated(spg_lat_vec)) deallocate(spg_lat_vec)
  allocate (spg_lat_vec(3,3) ,stat=info)
  call check_allocation(info, 'spg_lat_vec                       ')

  if (allocated(spg_atom_type)) deallocate(spg_atom_type)
  allocate (spg_atom_type(n_atoms) ,stat=info)
  call check_allocation(info, 'spg_atom_type                     ')


  call cite_reference("SPGlib")
  call cart2frac(lattice_vector, coords, spg_coords)
  spg_lat_vec = transpose(lattice_vector)

  spg_species_type = -1
  spg_atom_type = -1
  n_types = 0
  do i_atom = 1, n_atoms, 1
    found = .false.
    do i_type = 1, n_types, 1
      if ((species(i_atom) .eq. spg_species_type(i_type))) then
        found = .true.
        spg_atom_type(i_atom) = i_type
      end if
    end do
    if (.not.found) then
      n_types = n_types + 1
      spg_atom_type(i_atom) = n_types
      spg_species_type(n_types) = species(i_atom)
    end if
  end do

  ! The allocatable components of dset get deallocated on scope exit
  dset = spg_get_dataset(spg_lat_vec, &
                         spg_coords, spg_atom_type, n_atoms, sym_precision)

end subroutine get_symmetry
!******


subroutine write_symm_info()

  !  PURPOSE
  !
  !    Print symmetry informations
  !
  !  USES

  !  ARGUMENTS
  !    none
  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE
  use, intrinsic :: iso_c_binding, only: c_char
  use dimensions, only: n_atoms, sym_precision
  use geometry, only: lattice_vector, species, coords
  use localorb_io, only: localorb_info, use_unit
  implicit none

  ! Local Variables

  integer :: space_group
  character(len=120)   :: info_str
  character(len=10,kind=c_char)   :: schoenflies

  call get_symmetry(lattice_vector, n_atoms, species, coords)

  write(info_str,'(2X,A)') ''
  call localorb_info(info_str,use_unit,'(A)')  
  write(info_str,'(2X,A)') "Symmetry information by spglib:"
  call localorb_info(info_str,use_unit,'(A)')
  write(info_str,'(2X,A,E8.1)') "| Precision set to ", sym_precision
  call localorb_info(info_str,use_unit,'(A)')   

  ! The following might be interesting to debug
  ! Currently atoms are only distuinguished based on their species
  ! We might want to add spin constraints at some point
  ! the easy fix is to add different species information for now
  !write(info_str,'(2X,A)') "| Unit cell: "
  !call localorb_info(info_str,use_unit,'(A)')        
  !write(info_str,'(2X,A1,5X,3(2X,A12))') &
  !   "|","a","b","c"
  !call localorb_info(info_str,use_unit,'(A)')        
  !do i_periodic = 1, 3, 1
  !  write(info_str,'(2X,A1,5X,3(2X,F12.6))') "|", &
  !    (spg_lat_vec(i_coord,i_periodic)*bohr, i_coord=1,3,1)
  !  call localorb_info(info_str,use_unit,'(A)')        
  !enddo

  !write(info_str,'(2X,A)') "| Atomic structure: "
  !call localorb_info(info_str,use_unit,'(A)')        
  !write(info_str,'(2X,A1,A25,3(2X,A12))') &
  !   "|","fractional coords","x","y","z"
  !call localorb_info(info_str,use_unit,'(A)')      
  !do i_atom = 1, n_atoms, 1
  !  write(info_str,'(2X,A1,I10,A10,I5,3(2X,F12.6))') &
  !    "|",i_atom, ": Species ", &
  !    spg_atom_type(i_atom), (spg_coords(i_coord,i_atom), i_coord=1,3,1)
  !  call localorb_info(info_str,use_unit,'(A)')       
  !enddo

  if (dset % spacegroup_number /= 0) then
     schoenflies = ' '
     space_group = spg_get_schoenflies( schoenflies, spg_lat_vec, &
           spg_coords, spg_atom_type, n_atoms, sym_precision )
     call clean_c_string(schoenflies)
     call clean_c_string(dset % international_symbol)

     write(info_str,'(2X,A,I0)') "| Number of Operations  : ", &
       dset%n_operations
     call localorb_info(info_str,use_unit,'(A)') 
     write(info_str,'(2X,A,I0)') "| Space group           : ", &
       dset%spacegroup_number
     call localorb_info(info_str,use_unit,'(A)')  
     write(info_str,'(2X,2A)')   "| International         : ", &
       trim((dset%international_symbol))
     call localorb_info(info_str,use_unit,'(A)')  
     write(info_str,'(2X,2A)')   "| Schoenflies           : ", &
       trim(schoenflies)
     call localorb_info(info_str,use_unit,'(A)')  
  else
     print '("Space group could not be found,")'
  end if

end subroutine write_symm_info
!******

!****f* spglib_symm/get_symmetry_for_lattice
!*  NAME
!*    get_symmetry_for_lattice
!*  SYNOPSIS
subroutine get_symmetry_for_lattice &
           (lattice_vector, n_atoms, species, coords, symprec, &
            spacegroup_found, n_operations_sym, spacegroup_number, &
            international_symbol, schoenflies)
!*  PURPOSE
!*    Returns the symmetry information obtained from spglib, dissociated from
!*    the SpglibDataset derived type to allow for usage in other parts of the
!*    code.
!*  USES
  use, intrinsic :: iso_c_binding, only: c_char
  use dimensions, only: sym_precision
  implicit none
!*  ARGUMENTS
  real*8 , intent(in)  :: lattice_vector(3,3)
  integer, intent(in)  :: n_atoms
  integer, intent(in)  :: species(n_atoms)
  real*8,  intent(in)  :: coords(3, n_atoms)
  real*8,  intent(out) :: symprec
  logical, intent(out) :: spacegroup_found
  integer, intent(out) :: n_operations_sym
  integer, intent(out) :: spacegroup_number
  character(len=14, kind=c_char), intent(out) :: international_symbol
  character(len=10, kind=c_char), intent(out) :: schoenflies

!*  INPUTS
!*    o lattice_vector - the lattice vectors
!*    o n_atoms - number of atoms
!*    o species - species for atoms
!*    o coords - Cartesian coordinates for atoms
!*  OUTPUTS
!*    o symprec - Precision used for determining symmetry operations
!*    o spacegroup_found - Whether a space group was found
!*    o n_operations_sym - The number of symmetry operations
!*    o international_symbol - The international symbol for the space group
!*    o schoenflies - The Schoenflies symbol for the space group
!*  AUTHORS
!*    William Huhn (Duke University)
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject
!*    the terms and conditions of the respective license agreement."
!*  SOURCE
  character(*), parameter :: func = 'get_symmetry_for_lattice'

  integer :: space_group

  symprec = sym_precision
  call get_symmetry(lattice_vector, n_atoms, species, coords)

  if (dset % spacegroup_number /= 0) then
     schoenflies = ' '
     space_group = spg_get_schoenflies( schoenflies, spg_lat_vec, &
           spg_coords, spg_atom_type, n_atoms, sym_precision )
     call clean_c_string(schoenflies)
     call clean_c_string(dset % international_symbol)

     spacegroup_found = .true.
     n_operations_sym = dset%n_operations
     spacegroup_number = dset%spacegroup_number
     international_symbol = trim(adjustl(dset%international_symbol))
     schoenflies = trim(adjustl(schoenflies))
  else
     spacegroup_found = .false.
     n_operations_sym = -1
     spacegroup_number = -1
     international_symbol = ""
     schoenflies = ""
  end if

end subroutine get_symmetry_for_lattice
!******

subroutine out_symm_mats()

  !  PURPOSE
  !
  !    Get rotations and translations from spg_lib
  !
  !  USES

  !  ARGUMENTS
  !    none
  !  INPUTS
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE
  ! Local Variables

  use mpi_tasks, only: check_allocation
  implicit none

  integer :: n, i, j, info
  character(len=120)  :: info_str
  integer :: temp(3,3)
  real*8  :: temp_shift(3)



  num_symmetry = dset % n_operations

  if (.not.allocated(spg_rotations))then
    allocate ( spg_rotations(3,3,num_symmetry) ,stat=info)
      call check_allocation(info, 'spg_rotations                       ')
  endif
  if (.not.allocated(spg_shift))then
    allocate ( spg_shift(3,num_symmetry) ,stat=info)
      call check_allocation(info, 'spg_shift                      ')
  endif
  do n = 1, num_symmetry
     spg_shift(1:3 , n) = dset % translations(1:3, n)
     do i = 1, 3
       do j = 1, 3
        spg_rotations(i ,j , n) = dset % rotations(j,i, n)
       end do
     enddo
  end do
  ! Make sure Identity is 1 and Inversion 2
  do n = 1, num_symmetry
    if((spg_rotations(1, 1 , n).eq.1).and.(spg_rotations(2, 2 , n).eq.1).and.&
      &(spg_rotations(3, 3 , n).eq.1).and.(spg_rotations(1, 2 , n).eq.0).and.&
      &(spg_rotations(1, 3 , n).eq.0).and.(spg_rotations(2, 1 , n).eq.0).and.&
      &(spg_rotations(2, 3 , n).eq.0).and.(spg_rotations(3, 1 , n).eq.0).and.&
      &(spg_rotations(3, 2 , n).eq.0))then
      if(n.eq.1)then
        exit
      endif
      temp(:,:)=spg_rotations(:, : , 1)
      spg_rotations(:, : , 1)=spg_rotations(:, : , n)
      spg_rotations(:, : , n)=temp(:,:)
      temp_shift(:) = spg_shift(:,1)
      spg_shift(:,1) = spg_shift(:,n)
      spg_shift(:,n) = temp_shift(:)
      exit
    endif
  end do
  do n = 1, num_symmetry
    if((spg_rotations(1, 1 , n).eq.-1).and.(spg_rotations(2, 2 , n).eq.-1).and.&
      &(spg_rotations(3, 3 , n).eq.-1).and.(spg_rotations(1, 2 , n).eq.0).and.&
      &(spg_rotations(1, 3 , n).eq.0).and.(spg_rotations(2, 1 , n).eq.0).and.&
      &(spg_rotations(2, 3 , n).eq.0).and.(spg_rotations(3, 1 , n).eq.0).and.&
      &(spg_rotations(3, 2 , n).eq.0))then
      if(n.eq.2)then
        exit
      endif
      temp(:,:)=spg_rotations(:, : , 2)
      spg_rotations(:, : , 2)=spg_rotations(:, : , n)
      spg_rotations(:, : , n)=temp(:,:)
      temp_shift(:) = spg_shift(:,2)
      spg_shift(:,2) = spg_shift(:,n)
      spg_shift(:,n) = temp_shift(:)
      exit
    endif
  end do
end subroutine out_symm_mats

!******

subroutine map_kpoints_sym(n_k_xyz,k_off)

  !  PURPOSE
  !
  !    Use rotations from spglib to find irreducible set of k-points
  !
  !  USES
  use dimensions, only: n_k_points_nosym
  use localorb_io, only: localorb_info, use_unit
  use mpi_tasks, only: check_allocation, n_tasks
  use runtime_choices, only: n_k_points_xyz_nosym
  implicit none
  !  ARGUMENTS
  integer, intent(IN) :: n_k_xyz(3)
  real*8, intent(IN) :: k_off(3)
  !    none
  !  INPUTS
  !    o n_k_xyz -- number of original k-points in the three directions
  !    o k_off -- k-point offsets (parameter k_offset in control.in)
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE
  ! Local Variables


  integer :: info, num_ir_grid
  character(len=120)  :: info_str

  ! copy original dimensions of k-point grid
  n_k_points_xyz_nosym(1:3) =  n_k_xyz(1:3)
  n_k_points_nosym = product(n_k_xyz)

  
  ! Allocations
  ! Map between k-point sets
  if (.not.allocated(map))then
    allocate ( map(n_k_points_nosym) ,stat=info)
      call check_allocation(info, 'map                       ')
  endif
  ! index of operation in spg_rotations
  if (.not.allocated(map_sym))then
    allocate ( map_sym(n_k_points_nosym) ,stat=info)
      call check_allocation(info, 'map_sym                   ')
  endif

  ! Map for inversion only
  if (.not.allocated(map_inv))then
    allocate ( map_inv(n_k_points_nosym) ,stat=info)
      call check_allocation(info, 'map_inv                   ')
  endif
  call get_kpoint_map_sym(k_off, n_k_points_nosym, n_k_xyz, map,&
                          map_sym, map_inv, num_ir_grid)
  write(info_str, '(2X,A,I8,A,I8)') &
  & '| k-points reduced from: ', product(n_k_xyz), ' to ', num_ir_grid
  call localorb_info(info_str,use_unit,'(A)')
  if (num_ir_grid.lt.n_tasks) then
      write(info_str, '(2X,A,I8,A,I8,A)') &
    & 'The number of remaining irreducible k-points ', num_ir_grid, ' is less than the number of tasks ', n_tasks, '.'
    call localorb_info(info_str,use_unit,'(A)')
        write(info_str, '(2X,A)') &
    & 'Symmetry based k-point reduction failed. Increase number of k-points of reduce tasks'
    call localorb_info(info_str,use_unit,'(A)')
    !call aims_stop()
  endif


end subroutine map_kpoints_sym
!******
subroutine get_kpoint_map_sym(k_points_offset, n_k_points_nosym, &
                              n_k_points_xyz_nosym, map, map_sym, map_inv, n_kred)

!  PURPOSE
!
!
!  Apply rotation matrix of crystal to all k-points and find the irr. 
!  k-points and their mapping to the original set.
!  Box for mapping back to first BZ is set up.
!
!  USES
!
      implicit none

!  ARGUMENTS
      real*8, intent(in) ::  k_points_offset(3)
      integer, intent(in) :: n_k_points_nosym
      integer, intent(in) :: n_k_points_xyz_nosym(3)
      integer, intent(out) :: map(n_k_points_nosym)
      integer, intent(out) :: map_sym(n_k_points_nosym)
      integer, intent(out) :: map_inv(n_k_points_nosym)
      integer, intent(out) :: n_kred
! INPUTS
!
! o k_points_offset -- k-point offset
! o n_k_points_nosym -- original number of k-points
! o n_k_points_xyz_nosym -- number of k-points in x, y, z
!
! OUTPUTS
!
! o map -- mapping between new and old k-points
! o map -- mapping between new and old k-points - inversion only
! o map_sym -- index of rotation for transforming old in red. k-point
! o n_kred -- number of irreducible k-points
!
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
!

    Real*8 :: kbox (3, 4)

    kbox (:, 1) = k_points_offset / dble  (n_k_points_xyz_nosym)
    kbox (:, 2) = kbox (:, 1)
    kbox (:, 3) = kbox (:, 1)
    kbox (:, 4) = kbox (:, 1)
    kbox (1, 2) = kbox (1, 2) + 1.d0
    kbox (2, 3) = kbox (2, 3) + 1.d0
    kbox (3, 4) = kbox (3, 4) + 1.d0
    ! Map inverion symmetry only
    if((spg_rotations(1, 1 , 2).eq.-1).and.(spg_rotations(2, 2 , 2).eq.-1).and.&
      &(spg_rotations(3, 3 , 2).eq.-1).and.(spg_rotations(1, 2 , 2).eq.0).and.&
      &(spg_rotations(1, 3 , 2).eq.0).and.(spg_rotations(2, 1 , 2).eq.0).and.&
      &(spg_rotations(2, 3 , 2).eq.0).and.(spg_rotations(3, 1 , 2).eq.0).and.&
      &(spg_rotations(3, 2 , 2).eq.0))then
      call genppts_sym (n_k_points_xyz_nosym, kbox, n_kred, map_inv, map_sym, 2)
    else
      call genppts_sym (n_k_points_xyz_nosym, kbox, n_kred, map_inv, map_sym, 1)
    end if
    
    ! Full Mapping
    call genppts_sym (n_k_points_xyz_nosym, kbox, n_kred, map, map_sym, num_symmetry)

endsubroutine get_kpoint_map_sym

subroutine genppts_sym (n_k_points_xyz, kbox, nk_red, map, map_sym,sym_max)

      !  PURPOSE
      !
      !  Apply rotation matrix of crystal to all k-points and find the irr. 
      !  k-points and their mapping to the original set
      !
      !  USES
      Implicit None
      ! ARGUMETNS
      integer, intent (In) :: n_k_points_xyz (3)
      integer, intent (In) :: sym_max
      real*8, intent (In) :: kbox (3, 4)
      integer, intent (Out) :: nk_red
      integer,  intent (Out) :: map (n_k_points_xyz(1)*n_k_points_xyz(2)*n_k_points_xyz(3)) 
      integer,  intent (Out) :: map_sym (n_k_points_xyz(1)*n_k_points_xyz(2)*n_k_points_xyz(3)) 
      !  INPUTS
      !    o n_k_points_xyz -- number of k-points in x, y, z
      !    o kbox -- k-point grid box
      !  OUTPUT
      !    o nk_red -- number of reduced points
      !    o map -- mapping between new and old k-points
      !    o map_sym -- index of rotation for transforming old in red. k-point
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
      !    Release version, FHI-aims (2016).
      !  SOURCE
      ! local variables
      integer :: i1, i2, i3, ik, jk, nk
      integer :: isym
      real*8 :: v1 (3), v2 (3), v3 (3)
      real*8 :: b (3, 3), s (3, 3), diff
      real*8 :: lim = 1.e-6
      integer :: map_trans(n_k_points_xyz(1)*n_k_points_xyz(2)*n_k_points_xyz(3))
      real*8 :: klist (3, n_k_points_xyz(1)*n_k_points_xyz(2)*n_k_points_xyz(3))
      !nk=0
      !do i1 = 1, n_k_points_xyz (1)
      !do i2 = 1, n_k_points_xyz (2)
      !do i3 = 1, n_k_points_xyz (3)
      !nk=nk+1
      !knum (i1,i2,i3) = nk
      !enddo
      !enddo
      !enddo
      
      ! k-point box vector matrix
      b (:, 1) = kbox (:, 2) - kbox (:, 1)
      b (:, 2) = kbox (:, 3) - kbox (:, 1)
      b (:, 3) = kbox (:, 4) - kbox (:, 1)
      ik = 0
      nk = 0
      map = 0
      map_sym = 1
      diff = 1.d0
      do i1 = 0, n_k_points_xyz (1) - 1
         v1 (1) = dble (i1) / dble (n_k_points_xyz(1))                       
         do i2 = 0, n_k_points_xyz (2) - 1
            v1 (2) = dble (i2) / dble (n_k_points_xyz(2))
            do i3 = 0, n_k_points_xyz (3) - 1
                  v1 (3) = dble (i3) / dble (n_k_points_xyz(3)) 
                  nk=nk+1
                  !nk = knum (i1+1,i2+1,i3+1)
                  v2 = matmul(b, v1)
                  v2 (:) = v2 (:) + kbox (:, 1)
                  call map_to_center_cell_frac (v2)
                  ! find equivalent in list klist
                  do isym = 1, sym_max
                     s (:, :) = dble (spg_rotations(:, :, isym))
                     v3 = matmul(transpose(s), v2)
                     call map_to_center_cell_frac (v3)
                     do jk = 1, ik
                        diff = Abs (klist(1, jk)-v3(1)) + Abs (klist(2, jk)-v3(2)) &
                           + Abs (klist(3, jk)-v3(3))
                        if (diff .lt. lim) then
                           ! equivalent k-point found, update map and sym-op.
                           map (nk) = map_trans(jk)
                           map_sym(nk) = isym
                           ! debug write(use_unit,'(A)') '----------------------'
                           ! debug write(use_unit,'(A,I4,A,I4)') 'k-point: ', nk, ',  sym-op: ',isym
                           ! debug write(use_unit,'(A,3E11.4)') 'k-vector: ',v2
                           ! debug write(use_unit,'(A)') 'Rot. matrix:'
                           ! debug write(use_unit,'(3I4)') spg_rotations(1,1:3, isym)
                           ! debug  write(use_unit,'(3I4)') spg_rotations(2,1:3, isym)
                           ! debug write(use_unit,'(3I4)') spg_rotations(3,1:3, isym)
                           ! debug write(use_unit,'(A,3E11.4)') 'rot. k-vector:', matmul(transpose(s), v2)
                           ! debug write(use_unit,'(A,3E11.4)') 'rot. k-vector, mapped to BZ: ',v3
                           exit
                        endif
                     enddo
                     if (diff .lt. lim) exit
                  enddo
               ! add new k-point to list
               if (map (nk).eq.0) then
               ik = ik + 1
               map_trans(ik) = nk
               map(nk) = nk
               klist (:, ik) = v2 (:)
               ! debug write(use_unit,'(A)') '----------------------'
               ! debug write(use_unit,'(A,I4,A,I4)') 'k-point: ', nk, ',  sym-op: ',isym
               ! debug write(use_unit,'(A,3E11.4)') 'k-vector: ',v2
               endif
            enddo
         enddo
      enddo
      nk_red = ik

endsubroutine

subroutine map_to_center_cell_frac (v)

      !  PURPOSE
      !
      !  Map vector in fractional coordinates back to the interval [0,1)
      !
      !  USES
      Implicit None
      ! ARGUMETNS
      real*8, intent (INOUT) :: v (3)
      !  INPUTS
      !    o v -- input vector
      !  OUTPUT
      !    o v -- output vector
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
      ! local variables
      Integer :: i
      Real*8 :: lim = 1.0e-6
      Integer :: iv (3)
      do i = 1, 3
         iv (i) = Int (v(i))
         v (i) = v (i) - dble (iv(i))
         if (v(i) .lt. 0.d0) then
            v (i) = v (i) + 1.d0
            iv (i) = iv (i) - 1
         endif
         if (1.d0-v(i) .lt. lim) then
            v (i) = 0.d0
            iv (i) = iv (i) + 1
         endif
         if (v(i) .lt. lim) then
            v (i) = 0.d0
         endif
      enddo

End Subroutine map_to_center_cell_frac

subroutine k_point_symmetry_check_spg(n_k_xyz, k_off, k_number)

  !  PURPOSE
  !
  !    Idea: Figure out system symmetries and reduce the k-points accordingly.
  !
  !    Reduction is done solely done by putting positive values into k_number
  !    entries of irreducible k-points and zero out entries of redundant
  !    k-points.
  !
  !    JW: With no symmetry operations on basis functions (and real-valued
  !    spherical harmonics) at hand, making use of system symmetry is not
  !    trivial.  But time inversion is a different issue, so use only that.
  !
  !  USES

  use localorb_io, only: localorb_info, use_unit
  use mpi_tasks, only: myid
  use runtime_choices
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_k_xyz(3)
  real*8, intent(IN) :: k_off(3)
  integer, intent(OUT) :: k_number(n_k_xyz(1), n_k_xyz(2), n_k_xyz(3))
  !  INPUTS
  !    o n_k_xyz -- number of original k-points in the three directions
  !    o k_off -- k-point offsets (parameter k_offset in control.in)
  !  OUTPUT
  !    o  k_number(i_k_xyz(1), i_k_xyz(2), i_k_xyz(3))
  !       -- contains the number of k-points that have been mapped on to the
  !       original entries, this is required to determine the proper k-weights
  !       in the end
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

  integer :: i_kx,i_ky,i_kz, i
  character*140 :: info_str
  character(*), parameter :: func = 'k_point_symmetry_check'
  write(info_str,"(10X,A)") "Using symmetry for reducing the k-points"
  call localorb_info(info_str)
  write(info_str,"(10X,A)") "Only first 100 mappings shown."
  call localorb_info(info_str)


  call destroy_symmetry()
  call destroy_symmats()
  call destroy_sym_maps()
  call write_symm_info()
  call out_symm_mats()
  call map_kpoints_sym(n_k_xyz,k_off)
  i=0
  k_number = 0
  do i_kx=1, n_k_xyz(1)
    do i_ky=1, n_k_xyz(2)
      do i_kz=1, n_k_xyz(3)
        i=i+1
        k_number(i_kx,i_ky,i_kz) = count(map == i)
        if (myid.eq.0.and.i.lt.101)then
          write(info_str,'(A,I4,A,I4,A,I4,A,I4)') '  | # k-point: ', i, &
                          ', multiplicity: ',k_number(i_kx,i_ky,i_kz), &
                          ', Mapped to # k-points: ', map(i),', # Sym. op.: ',&
                          map_sym(i)
          call localorb_info(info_str,use_unit,'(A)')
        endif
      end do
    end do
  end do

  ! This is for debbuging!
  ! call remap(k_number,n_k_xyz,map,map_sym,8)

  ! Deallocate everything except sym mats
  call destroy_symmetry()
end subroutine k_point_symmetry_check_spg
!******

subroutine clean_c_string (string_in)
  use, intrinsic :: iso_c_binding, only: c_char, c_null_char

  ! Local Variables
  implicit none

  character(len=*,kind=c_char), intent(inout) :: string_in
  character(len=:), allocatable               :: test
  integer                                     :: i

  allocate (character(len=len(string_in)) :: test) 
  test = string_in
  do i = 1, len(test)
    select case (test(i:i))
    case ('A':'Z', 'a' : 'z', '0' : '9','_','/','^','-')
      string_in(i:i) = test(i:i)
    case default
      string_in(i:i) = ' '
    end select
  end do

  deallocate(test)
   
end subroutine clean_c_string

subroutine symmetrize_forces(forces)

  !  PURPOSE
  !
  !    Interfaces to spglib to obtain symmetry informations
  !
  !  USES

  use dimensions, only: n_atoms
  use localorb_io, only: localorb_info, use_unit
  implicit none

  !  ARGUMENTS
  !    o forces ({x,y,z},n_atoms) forces of the system
  !  INPUTS
  real*8, dimension(3,n_atoms), intent(inout) :: forces
  !  OUTPUTS
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables

  real*8, dimension(3) :: work
  integer              :: i_atom, i_symm
  character(len=120)   :: info_str

  do i_atom = 1, n_atoms
     write(info_str,'(2X,A,I3,A,I3)') "| Atom  : ", i_atom, &
       " maps to atom ", dset % equivalent_atoms(i_atom)
     call localorb_info(info_str,use_unit,'(A)')
     do i_symm = 1, dset % n_operations
        work(1:3) = dset % translations(1:3,i_symm) &
          + matmul(spg_coords(1:3,i_atom),dset % rotations(1:3,1:3,i_symm))
     !   write(info_str,'(2X,A6,I3,A10,I3,A,3(2X,F12.6))') &
     !      "|SYMM ",i_symm, &
     !      " maps atom ", i_atom, &
     !      " to ", (work(i_coord), i_coord = 1,3,1)
     !   call localorb_info(info_str,use_unit,'(A)')
         work = map_to_cell(work)
         write(info_str,'(2X,A6,I3,A10,I3,A,I3)') &
           "|SYMM ",i_symm, &
           " maps atom ", i_atom, &
           " to ", get_atom(work)
    call localorb_info(info_str,use_unit,'(A)')

     end do
  end do

  !do i_symm = 1, dset % n_operations
  !  write(info_str,'(2X,A6,I10,A10,3(2X,I5))') &
  !    "|ROT  ",i_symm, ": ", (dset % rotations(1,i_coord,i_symm), i_coord=1,3,1)
  !  call localorb_info(info_str,use_unit,'(A)')
  !  write(info_str,'(2X,A6,I10,A10,3(2X,I5))') &
  !    "|ROT  ",i_symm, ": ", (dset % rotations(2,i_coord,i_symm), i_coord=1,3,1)
  !  call localorb_info(info_str,use_unit,'(A)')
  !  write(info_str,'(2X,A6,I10,A10,3(2X,I5))') &
  !    "|ROT  ",i_symm, ": ", (dset % rotations(3,i_coord,i_symm), i_coord=1,3,1)
  !  call localorb_info(info_str,use_unit,'(A)')
  !  write(info_str,'(2X,A6,I10,A10,3(2X,F12.6))') &
  !    "|TRANS",i_symm, ": ", (dset % translations(i_coord,i_symm), i_coord=1,3,1)
  !  call localorb_info(info_str,use_unit,'(A)')
  !end do  
  

end subroutine symmetrize_forces


function map_to_cell(frac_coords) result(map_coords)

   implicit none
   
   real*8, dimension(3) :: frac_coords
   real*8, dimension(3) :: map_coords

   map_coords = frac_coords - 1d0*floor(frac_coords)

   return


end function map_to_cell

function get_atom(frac_coords) result(atom)

   use dimensions, only: n_atoms
   use localorb_io, only: localorb_info, use_unit
   use mpi_tasks, only: aims_stop
   implicit none
   
   real*8, dimension(3) :: frac_coords
   integer              :: atom
   integer              :: i_atom, i_coord
   real*8               :: delta
   character(len=120)   :: info_str

   atom = -1

   do i_atom = 1, n_atoms
      delta = sqrt((frac_coords(1) - spg_coords(1,i_atom))**2 &
               +   (frac_coords(2) - spg_coords(2,i_atom))**2 &
               +   (frac_coords(3) - spg_coords(3,i_atom))**2)
      if (delta .lt. 1e-2) then
        atom = i_atom
        exit
      end if       
   end do

   if (atom .eq. -1) then
    write(info_str,'(2X,A8,3(2X,F12.6))') &
      "|Error  ", (frac_coords(i_coord), i_coord=1,3,1)
    call localorb_info(info_str,use_unit,'(A)')
    call aims_stop()
   end if

   return


end function get_atom


subroutine destroy_symmetry()

  !  PURPOSE
  !
  !  Free Memory   
  !
  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables
  implicit none

  if (allocated(spg_coords))    deallocate ( spg_coords )
  if (allocated(spg_lat_vec))   deallocate ( spg_lat_vec )
  if (allocated(spg_atom_type)) deallocate ( spg_atom_type )

  if (allocated(dset % rotations))        deallocate( dset % rotations )
  if (allocated(dset % translations))     deallocate( dset % translations )
  if (allocated(dset % wyckoffs))         deallocate( dset % wyckoffs )
  if (allocated(dset % equivalent_atoms)) deallocate( dset % equivalent_atoms)
  if (allocated(dset % brv_types))        deallocate( dset % brv_types )
  if (allocated(dset % brv_positions))    deallocate( dset % brv_positions )


end subroutine destroy_symmetry

subroutine destroy_symmats()

  !  PURPOSE
  !
  !  Free Memory   
  !
  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables
  implicit none


  if (allocated(spg_rotations)) deallocate ( spg_rotations )
  if (allocated(spg_shift)) deallocate ( spg_shift )


end subroutine destroy_symmats

subroutine destroy_sym_maps()

  !  PURPOSE
  !
  !  Free Memory   
  !
  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables
  implicit none


  if (allocated(map)) deallocate ( map )
  if (allocated(map_sym)) deallocate ( map_sym )
  if (allocated(map_inv)) deallocate ( map_inv )

end subroutine destroy_sym_maps

subroutine remap(k_number,n_k_xyz,map,map_sym, max_sym)

  use localorb_io, only: use_unit
  implicit none
  !  PURPOSE
  !
  !  Only map k-points with rotation matrix not larger than #max_sym
  !  Routine is for debugging
  !
  !  ARGUMENTS
  !  o n_k_xyz - # k-points in x, y, z
  !  o max_sym - Maximum # rotation matrix
  !  o k_number - muliplicity of k-points
  !  o map - map between full and reduced grid
  !  o map_sym - rotations to transform irreduceble set to full set
  !  INPUTS
  integer, intent(IN) :: n_k_xyz(3)
  integer, intent(IN) :: max_sym
  !  OUTPUTS
  integer, intent(INOUT) :: k_number(n_k_xyz(1), n_k_xyz(2), n_k_xyz(3))
  integer, intent(INOUT) :: map(product(n_k_xyz))
  integer, intent(INOUT) :: map_sym(product(n_k_xyz))
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  ! Local Variables
  integer :: i_kx,i_ky,i_kz, i
  integer :: k_map(product(n_k_xyz),3)


  i=0
  do i_kx=1, n_k_xyz(1)
  do i_ky=1, n_k_xyz(2)
  do i_kz=1, n_k_xyz(3)
  i=i+1
  k_map(i,1) = i_kx
  k_map(i,2) = i_ky
  k_map(i,3) = i_kz
  end do
  end do
  end do

  i=0
  write(use_unit,'(A)') ''
  write(use_unit,'(A)') '  | Full k-point mapping'
  do i_kx=1, n_k_xyz(1)
  do i_ky=1, n_k_xyz(2)
  do i_kz=1, n_k_xyz(3)
    i=i+1
    write(use_unit,'(A,I4,A,I4,A,I4,A,I4)') '  | # k-point: ', i,', multiplicity: ',&
              k_number(i_kx,i_ky,i_kz), ', Mapped to # k-points: ', map(i),&
              ', # Sym. op.: ',map_sym(i)
    if (map_sym(i).gt.max_sym-1)then
      k_number(i_kx,i_ky,i_kz) = k_number(i_kx,i_ky,i_kz) + 1
      k_number(k_map(map(i),1),k_map(map(i),2),k_map(map(i),3))=&
                     k_number(k_map(map(i),1),k_map(map(i),2),k_map(map(i),3))-1
      map(i)=i
      map_sym(i)=1
    endif

  end do
  end do
  end do
  write(use_unit,'(A)') ''
  write(use_unit,'(A)') '  | Realised k-point mapping'
  i=0
  do i_kx=1, n_k_xyz(1)
  do i_ky=1, n_k_xyz(2)
  do i_kz=1, n_k_xyz(3)
    i=i+1
    write(use_unit,'(A,I4,A,I4,A,I4,A,I4)') '  | # k-point: ', i,', multiplicity: ',&
             k_number(i_kx,i_ky,i_kz), ', Mapped to # k-points: ', &
             map(i),', # Sym. op.: ',map_sym(i)
  end do
  end do
  end do
  write(use_unit,'(A)') ''

endsubroutine remap


end module spglib_symmetry
