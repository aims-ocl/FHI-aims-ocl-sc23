module test_symmetry_of_grid_and_things_on_grid
!!
!! Collection of routines to measury how symmetric, or not, properties in AIMS are.
!! This is intended for debugging/diagnostics. For production runs doing these tests
!! serves no purpose whatsoever, so don't do that.
!!
use rlsy_interface, only: rlsy_handle

! Normal AIMS modules
use grids, only: batch_of_points

implicit none
private

public :: rl_test_symmetry_of_grid_and_things_on_grid
public :: rl_track_progress

contains

subroutine rl_test_symmetry_of_grid_and_things_on_grid(rlsy_h,&
    rho,delta_rho,rho_free_superpos,rho_gradient,delta_rho_gradient,&
    weight_tab,partition_tab,&
    r_radial,r_angular,&
    n_aims_batches,aims_batches)
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> Densities in different ways from AIMS
    real*8, dimension(:,:), intent(in) :: rho
    real*8, dimension(:,:), intent(in) :: delta_rho
    real*8, dimension(:), intent(in) :: rho_free_superpos
    real*8, dimension(:,:,:), intent(in) :: rho_gradient
    real*8, dimension(:,:,:), intent(in) :: delta_rho_gradient
    !> grid things from AIMS
    real*8, dimension(:), intent(in) :: weight_tab
    real*8, dimension(:), intent(in) :: partition_tab
    !> point things from AIMS
    real*8, dimension(:,:), intent(in) :: r_radial
    real*8, dimension(:,:,:,:), intent(in) :: r_angular
    !> Normal batches from AIMS
    integer, intent(in) :: n_aims_batches
    type(batch_of_points), dimension(:), intent(in) :: aims_batches

    ! Does nothing
end subroutine

subroutine rl_track_progress(rlsy_h,&
    number_of_loops,total_energy,previous_total_energy,&
    rho_change,ev_sum,previous_ev_sum,diff_forces,diff_stress,&
    AS_stress_on,forces_on)

    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> things grabbed from inside scf_solver
    integer, intent(in) :: number_of_loops
    real*8, intent(in) :: total_energy
    real*8, intent(in) :: previous_total_energy
    real*8, dimension(:), intent(in) :: rho_change
    real*8, intent(in) :: ev_sum
    real*8, intent(in) :: previous_ev_sum
    real*8, intent(in) :: diff_forces
    real*8, intent(in) :: diff_stress
    logical, intent(in) :: AS_stress_on
    logical, intent(in) :: forces_on

    ! Does nothing
end subroutine

end module
