module generate_index_map
!!
!! When using symmetry you have to be able to collect and symmetrize the Hamiltonian,
!! overlap, densitymatrix and perhaps some other matrices.
!!
use rlsy_interface, only: rlsy_handle
use grids, only: batch_of_points

implicit none
private

public :: rl_generate_index_map
public :: rl_collect_overlap
public :: rl_collect_hamiltonian
public :: rl_collect_densitymatrix
public :: rl_solve_KS_equations
public :: rl_calculate_density
public :: rl_check_if_reinit_needed
public :: rl_inject_density

contains

!> check if a re-initialization of the symmetry handle is needed. Will have to add more things as we progress.
subroutine rl_check_if_reinit_needed(rlsy_h,frac_coords,lattice_vector,need_reinit)
    type(rlsy_handle), intent(inout) :: rlsy_h
    real*8, dimension(:,:), intent(in) :: frac_coords
    real*8, dimension(3,3), intent(in) :: lattice_vector
    logical, intent(out) :: need_reinit

    ! Do nothing
end subroutine

!> This routine collects the Overlap matrix from AIMS. It assumes that the index is already generated.
subroutine rl_collect_overlap(rlsy_h,overlap_matrix)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real*8, dimension(:), intent(in) :: overlap_matrix

    ! Does nothing
end subroutine

!> This routine collects the Overlap matrix from AIMS. It assumes that the index is already generated.
subroutine rl_collect_densitymatrix(rlsy_h,densitymatrix)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real*8, dimension(:), intent(in) :: densitymatrix

    ! Does nothing
end subroutine

!> This routine collects the Overlap matrix from AIMS. It assumes that the index is already generated.
subroutine rl_collect_hamiltonian(rlsy_h,hamiltonian)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real*8, dimension(:,:), intent(in) :: hamiltonian

    ! Does nothing
end subroutine

!> This routine injects the Hamiltonian back into AIMS. It assumes that the index is already generated.
subroutine rl_inject_hamiltonian(rlsy_h,hamiltonian)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> overlap matrix, in whatever format it may be
    real*8, dimension(:,:), intent(inout) :: hamiltonian

    ! Does nothing
end subroutine

!> This is very annoying, but I have to re-index the Hamiltonian into something useful.
subroutine rl_generate_index_map(rlsy_h,&
    coords_center, center_to_cell, center_to_atom, &
    index_hamiltonian, column_index_hamiltonian, cbasis_to_basis, &
    cbasis_to_center, centers_basis_integrals, frac_coords, lattice_vector,&
    position_in_hamiltonian)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> things from AIMS such that I can figure out the indexing:
    !> AIMS coords_centers, Where is center i
    real*8, intent(in), dimension(:,:) :: coords_center
    !> AIMS center_to_cell. What vector moves center i back to the unit cell
    integer, intent(in), dimension(:) :: center_to_cell
    !> AIMS center_to_atom, what index in the unit cell is center i?
    integer, intent(in), dimension(:) :: center_to_atom
    !> AIMS index Hamiltonian
    integer, intent(in), dimension(:,:,:) :: index_hamiltonian
    !> AIMS column index Hamiltonian
    integer, intent(in), dimension(:) :: column_index_hamiltonian
    !> AIMS indices that map the all the basis functions to something I can understand.
    integer, intent(in), dimension(:) :: cbasis_to_basis
    integer, intent(in), dimension(:) :: cbasis_to_center
    integer, intent(in), dimension(:) :: centers_basis_integrals
    !> AIMS coordinates of atoms
    real*8, intent(in), dimension(:,:) :: frac_coords
    real*8, intent(in), dimension(3,3) :: lattice_vector
    !> AIMS position_in_hamiltonian
    integer, dimension(:,:), intent(in) :: position_in_hamiltonian

    ! Does nothing
end subroutine

!> This routine injects the Hamiltonian back into AIMS. It assumes that the index is already generated.
subroutine rl_solve_KS_equations(rlsy_h)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h

    ! Does nothing
end subroutine

!> evaluates the density
subroutine rl_calculate_density(rlsy_h)
    !> MPI helper
    type(rlsy_handle), intent(inout) :: rlsy_h

    ! Does nothing
end subroutine

!> replace the AIMS density and gradient of density with the one I have calculated.
subroutine rl_inject_density(rlsy_h,&
    rho,rho_gradient,partition_tab,&
    n_aims_batches,aims_batches)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> Densities in different ways from AIMS
    real*8, dimension(:,:), intent(inout) :: rho
    real*8, dimension(:,:,:), intent(inout) :: rho_gradient
    !> grid things from AIMS
    real*8, dimension(:), intent(in) :: partition_tab
    !> Normal batches from AIMS
    integer, intent(in) :: n_aims_batches
    type(batch_of_points), dimension(:), intent(in) :: aims_batches

    ! Does nothing
end subroutine

end module
