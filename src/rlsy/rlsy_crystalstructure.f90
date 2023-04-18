module rlsy_crystalstructure
!!
!! Container for a crystal structure. Stripped down version of my real
!! container to make sure it is somewhat f2003-compliant.
!!
use rlsy_constants, only: r8,rl_huge,rl_hugeint
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_mpi_helper
use rlsy_helpers, only: rl_invert3x3matrix,rl_chop,rl_clean_fractional_coordinates,rl_determ

implicit none
private
public :: rl_crystalstructure

! Crystal structure class. A tiny fraction of what I actually have.
! Repeated here because my symmetry detection routines rely intimiately
! on this kind of object.
type rl_crystalstructure
    !> number of atoms in the cell
    integer :: n_atom=-rl_hugeint
    !> number of species in the cell
    integer :: n_species=-rl_hugeint
    !> basis vectors
    real(r8), dimension(3,3) :: latticevectors=rl_huge
    !> inverse of basis vectors
    real(r8), dimension(3,3) :: inv_latticevectors=rl_huge
    !> reciprocal lattice vectors
    real(r8), dimension(3,3) :: reciprocal_latticevectors=rl_huge
    !> inverse reciprocal lattice vectors
    real(r8), dimension(3,3) :: inv_reciprocal_latticevectors=rl_huge
    !> Atomic number for each atom (n_atom)
    integer, dimension(:), allocatable :: atomic_number
    !> What species is atom i? (n_atom)
    integer, dimension(:), allocatable :: species
    !> Positions, fractional coordinates (3,n_atom)
    real(r8), dimension(:,:), allocatable :: fractional_coordinate
    !> Positions, cartesian coordinates (3,n_atom)
    real(r8), dimension(:,:), allocatable :: cartesian_coordinate
    !> volume of cell
    real(r8) :: volume=-rl_huge
    contains
        !> create the structure
        procedure :: generate
        !> measure size in memory, approximately
        procedure :: size_in_mem=>structure_size_in_mem
        !> destroy the object
        procedure :: destroy
end type rl_crystalstructure

contains

!> Transform lists of stuff to a proper object. How a crystal structure should be initialized. This is but a tiny fraction of
subroutine generate(p,latticevectors,positions,atomic_number,species,mw,mem)
    !> crystal structure
    class(rl_crystalstructure), intent(out) :: p
    !> basis
    real(r8), dimension(3,3), intent(in) :: latticevectors
    !> positions of atoms, in fractional coordinates
    real(r8), dimension(:,:), intent(in) :: positions
    !> atomic numbers of each atom
    integer, dimension(:), intent(in) :: atomic_number
    !> which species is each atom
    integer, dimension(:), intent(in) :: species
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    !@TODO move sanity tests from real code here

    ! Man I miss block structures. Anyway, start filling out things in
    ! the order they are needed. This type is really barebones at the moment.
    p%latticevectors=latticevectors
    p%inv_latticevectors=rl_invert3x3matrix(p%latticevectors)
    p%reciprocal_latticevectors=transpose(rl_invert3x3matrix(p%latticevectors))
    p%inv_reciprocal_latticevectors=rl_invert3x3matrix(p%reciprocal_latticevectors)
    p%n_atom=size(positions,2)

    ! Store coordinates
    !call mem%allocate(p%cartesian_coordinate,[3,p%n_atom],persistent=.true.,scalable=.false.)
    !call mem%allocate(p%fractional_coordinate,[3,p%n_atom],persistent=.true.,scalable=.false.)
    allocate(p%cartesian_coordinate(3,p%n_atom))
    allocate(p%fractional_coordinate(3,p%n_atom))
    p%fractional_coordinate=rl_clean_fractional_coordinates(positions)
    p%fractional_coordinate=rl_chop(p%fractional_coordinate,1E-12_r8)
    p%fractional_coordinate=rl_clean_fractional_coordinates(positions)
    p%cartesian_coordinate=matmul(p%latticevectors,p%fractional_coordinate)

    ! Check that we agree
    call mw%check_and_sync(p%fractional_coordinate,0)
    call mw%check_and_sync(p%cartesian_coordinate,0)

    ! Store some metadata
    !call mem%allocate(p%atomic_number,p%n_atom,persistent=.true.,scalable=.false.)
    !call mem%allocate(p%species,p%n_atom,persistent=.true.,scalable=.false.)
    allocate(p%atomic_number(p%n_atom))
    allocate(p%species(p%n_atom))
    p%atomic_number=atomic_number
    p%species=species
    p%n_species=maxval(species)
    p%volume=abs(rl_determ(p%latticevectors))
end subroutine generate

!> measure size in memory, in bytes, roughly
function structure_size_in_mem(p) result(mem)
    !> dispersions
    class(rl_crystalstructure), intent(in) :: p
    !> memory in bytes
    integer :: mem

    mem=0
    mem=mem+storage_size(p)
    if ( allocated(p%atomic_number) )         mem=mem+storage_size(p%atomic_number)*size(p%atomic_number)
    if ( allocated(p%species) )               mem=mem+storage_size(p%species)*size(p%species)
    if ( allocated(p%fractional_coordinate) ) mem=mem+storage_size(p%fractional_coordinate)*size(p%fractional_coordinate)
    if ( allocated(p%cartesian_coordinate) )  mem=mem+storage_size(p%cartesian_coordinate)*size(p%cartesian_coordinate)
    mem=mem/8
end function

!> destroy the object
subroutine destroy(p,mem)
    !> crystal structure
    class(rl_crystalstructure), intent(inout) :: p
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    if ( allocated(p%atomic_number) )         call mem%deallocate(p%cartesian_coordinate,persistent=.true.,scalable=.false.)
    if ( allocated(p%species) )               call mem%deallocate(p%fractional_coordinate,persistent=.true.,scalable=.false.)
    if ( allocated(p%fractional_coordinate) ) call mem%deallocate(p%atomic_number,persistent=.true.,scalable=.false.)
    if ( allocated(p%cartesian_coordinate) )  call mem%deallocate(p%species,persistent=.true.,scalable=.false.)
end subroutine

end module rlsy_crystalstructure
