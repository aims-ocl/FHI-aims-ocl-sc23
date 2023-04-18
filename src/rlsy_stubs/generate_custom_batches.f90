module generate_custom_batches
!!
!! When using symmetry you have to use different sets of batches in diffent situations,
!! this routine handles that.
!!
use rlsy_interface, only: rlsy_handle

! Normal AIMS modules
use grids, only: batch_of_points

implicit none
private

public :: rl_generate_batches_for_density_update

contains

!> progress tracker, for debugging/understanding. More things will be added as needed.
subroutine rl_generate_batches_for_density_update(rlsy_h,n_batches,sym_batches,n_full_points,partition_tab) !,&
    !> symmetry handle
    type(rlsy_handle), intent(inout) :: rlsy_h
    !> how many batches
    integer, intent(out) :: n_batches
    !> actual batches
    type(batch_of_points), dimension(:), allocatable, intent(out) :: sym_batches
    !> how many full points
    integer, intent(out) :: n_full_points
    !> partition_tab
    real*8, dimension(:), allocatable, intent(out) :: partition_tab

    ! Does nothing
end subroutine

end module
