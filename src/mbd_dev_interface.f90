module mbd_interface

use synchronize_mpi_basic, only: &
    sync_vector, sync_real_number, sync_vector_complex, bcast_real_vector, &
    bcast_complex_vector
use localorb_io, only: OL_norm, OL_high, localorb_info

implicit none

real(8) :: pi
integer, parameter :: legendre_precision = 8

interface sync_sum
    module procedure sync_sum_dble_
    module procedure sync_sum_vector_dble_
    module procedure sync_sum_matrix_dble_
    module procedure sync_sum_3d_dble_
    module procedure sync_sum_4d_dble_
    module procedure sync_sum_vector_cmplx_
    module procedure sync_sum_matrix_cmplx_
    module procedure sync_sum_3d_cmplx_
    module procedure sync_sum_4d_cmplx_
end interface

interface broadcast
    module procedure broadcast_vector_dble_
    module procedure broadcast_matrix_dble_
    module procedure broadcast_3d_dble_
    module procedure broadcast_4d_dble_
    module procedure broadcast_vector_cmplx_
    module procedure broadcast_matrix_cmplx_
end interface

contains

subroutine sync_sum_dble_(x)
    real(8), intent(inout) :: x

    call sync_real_number(x)
end subroutine

subroutine sync_sum_vector_dble_(x)
    real(8), intent(inout) :: x(:)

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_matrix_dble_(x)
    real(8), intent(inout) :: x(:, :)

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_3d_dble_(x)
    real(8), intent(inout) :: x(:, :, :)

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_4d_dble_(x)
    real(8), intent(inout) :: x(:, :, :, :)

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_vector_cmplx_(x)
    complex(kind=8), intent(inout) :: x(:)

    call sync_vector_complex(x, size(x))
end subroutine

subroutine sync_sum_matrix_cmplx_(x)
    complex(kind=8), intent(inout) :: x(:, :)

    call sync_vector_complex(x, size(x))
end subroutine

subroutine sync_sum_3d_cmplx_(x)
    complex(kind=8), intent(inout) :: x(:, :, :)

    call sync_vector_complex(x, size(x))
end subroutine

subroutine sync_sum_4d_cmplx_(x)
    complex(kind=8), intent(inout) :: x(:, :, :, :)

    call sync_vector_complex(x, size(x))
end subroutine

subroutine broadcast_vector_dble_(x)
    real(8), intent(inout) :: x(:)

    call bcast_real_vector(x, size(x), 0)
end subroutine

subroutine broadcast_matrix_dble_(x)
    real(8), intent(inout) :: x(:, :)

    call bcast_real_vector(x, size(x), 0)
end subroutine

subroutine broadcast_3d_dble_(x)
    real(8), intent(inout) :: x(:, :, :)

    call bcast_real_vector(x, size(x), 0)
end subroutine

subroutine broadcast_4d_dble_(x)
    real(8), intent(inout) :: x(:, :, :, :)

    call bcast_real_vector(x, size(x), 0)
end subroutine

subroutine broadcast_vector_cmplx_(x)
    complex(8), intent(inout) :: x(:)

    call bcast_complex_vector(x, size(x), 0)
end subroutine

subroutine broadcast_matrix_cmplx_(x)
    complex(8), intent(inout) :: x(:, :)

    call bcast_complex_vector(x, size(x), 0)
end subroutine

subroutine print_log (str, mute)
    character(len=*), intent(in) :: str
    logical, optional :: mute

    if (present(mute)) then
        if (mute) return
    end if
    call localorb_info("  | " // str, priority=OL_norm)
end subroutine print_log

subroutine print_warning (str)
    character(len=*), intent(in) :: str
    call localorb_info("*** Warning: " // str, priority=OL_high)
end subroutine print_warning

subroutine print_error (str)
    character(len=*), intent(in) :: str
    call localorb_info("*** Error: " // str, priority=OL_high)
end subroutine print_error

end module mbd_interface
