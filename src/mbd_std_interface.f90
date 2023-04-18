module mbdvdw_interface_module

use constants, only: pi
use dimensions, only: nat => n_occ_atoms, mbd_vdw_forces => use_forces
use geometry, only: ityp => species
use species_data, only: atm => species_name
use mpi_tasks, only: &
    nproc_image => n_tasks, me_image => myid, &
    intra_image_comm => mpi_comm_global, stdout
use mpi_tasks
use synchronize_mpi_basic, only: &
    sync_real_number, sync_vector, sync_vector_integer
use hirshfeld, only: hirshfeld_wfforce => contract_with_dVdn

implicit none

integer, parameter :: dp = 8

type :: Fft_base_t
    integer :: nnr
end type

type(Fft_base_t) :: dfftp

integer :: root_image = 0

logical :: mbd_vdw_isolated = .true.
real(8) :: mbd_vdw_supercell = 0.d0  ! supercell radius
logical :: vdw_self_consistent = .false.
logical :: use_hirshfeld_deriv = .false.
real(8) :: mbd_vdw_beta = 0.d0  ! beta parameter for damping
real(8) :: mbd_vdw_econv_thr = 1.d-6
integer :: mbd_vdw_n_quad_pts = 15
logical :: use_old_hirshfeld = .false.
logical :: mbd_vdw_ewald = .true.
logical :: mbd_vdw_recip = .false.
logical :: vdw_debug = .false.
logical :: mbd_vdw_low_dim = .false.
integer :: mbd_vdw_kgrid(3) = (/ 0, 0, 0 /)
real(8) :: mbd_vdw_kgrid_shift = 0.5d0
logical :: mbd_vdw_vacuum(3) = .false.
integer :: mbd_vdw_verbosity = -1
logical :: timing = .false.
integer :: mbd_max_negative_eigvals = 3

type :: Timestamp_t
    character(len=30) :: label = ''
    integer :: cnt = 0
end type

integer, parameter :: n_timestamps = 30
type(Timestamp_t) :: timestamps(n_timestamps)

real(8), allocatable :: vefftsvdw(:)  ! (n_atoms, n_atoms, 3)
real(8), allocatable :: vfree(:)  ! (n_atoms, n_atoms, 3)
real(8), allocatable :: dveffdr(:, :, :)  ! (n_atoms, 3, 3)
real(8), allocatable :: dveffdh(:, :, :)  ! (n_atoms, 3, 3)

! QE-related flags
logical :: vdw_lstep0, vdw_lscreen
integer :: iverbosity

interface mp_sum
    module procedure sync_sum_dble_
    module procedure sync_sum_vector_dble_
    module procedure sync_sum_vector_int_
    module procedure sync_sum_matrix_dble_
    module procedure sync_sum_3d_dble_
    module procedure sync_sum_4d_dble_
end interface

external :: hirshfeld_analysis

contains

subroutine GetVdwParam(element, C6, alpha, R_vdw)
    use vdw_correction, only: get_vdw_param_aims => get_vdw_param, &
        vdw_hirshfeld_data_external, vdw_hirshfeld_C6, vdw_hirshfeld_alpha, &
        vdw_hirshfeld_R0

    character(len=2), intent(in) :: element
    real(8), intent(out) :: C6, alpha, R_vdw

    real(8) :: dummy
    integer :: i_specie

    do i_specie = 1, size(ityp)
        if (atm(i_specie) == element) exit
    end do
    if (vdw_hirshfeld_data_external(i_specie)) then
        C6 = vdw_hirshfeld_C6(i_specie)
        alpha = vdw_hirshfeld_alpha(i_specie)
        R_vdw = vdw_hirshfeld_R0(i_specie)
    else
        call get_vdw_param_aims( &
            element, dummy, C6, alpha, R_vdw &
        )
    end if
end subroutine

subroutine hirshfeld_initialize()
end subroutine

subroutine hirshfeld_calculate(dummy_a, rho, dummy_c, dummy_d, dummy_e)
    use vdw_correction, only: hirshfeldw, freeintegral
    use hirshfeld, only: run_hirshfeld, evaluate_free_atom_quantities
    use geometry, only: species
    use dimensions, only: n_periodic

    real(8) :: rho(:, :)
    real(8) :: dummy_a(:, :), dummy_c(:, :), dummy_d(:, :)
    logical :: dummy_e

    integer :: i_atom

    if (use_old_hirshfeld) then
        call hirshfeld_analysis()
        vefftsvdw = hirshfeldw
        do i_atom = 1, nat
            vfree(species(i_atom)) = freeintegral(i_atom)
        end do
        if (mbd_vdw_forces) then
            dveffdr(:, :, :) = 0.d0
            dveffdh(:, :, :) = 0.d0
        end if
    else
        if (mbd_vdw_forces) then
            if (use_hirshfeld_deriv) then
                if (n_periodic > 0) then
                    call run_hirshfeld(rho, vefftsvdw, dveffdr, dveffdh)
                else
                    call run_hirshfeld(rho, vefftsvdw, dveffdr)
                end if
            else
                call run_hirshfeld(rho, vefftsvdw)
                dveffdr(:, :, :) = 0.d0
                dveffdh(:, :, :) = 0.d0
            end if
        else
            call run_hirshfeld(rho, vefftsvdw)
        end if
        call evaluate_free_atom_quantities(vfree)
    end if
end subroutine

subroutine hirshfeld_cleanup()
end subroutine

integer function get_stamp(label) result(i_stamp)
    use mpi_tasks, only: aims_stop

    character(len=*), intent(in) :: label

    do i_stamp = 1, n_timestamps
        if (timestamps(i_stamp)%label == label) then
            exit
        else if (timestamps(i_stamp)%label == '') then
            timestamps(i_stamp)%label = label
            exit
        end if
    end do
    if (i_stamp > n_timestamps) call aims_stop('MBD: No free timestamps')
end function

subroutine start_clock(label)
    character(len=*), intent(in) :: label

    integer :: ts_cnt, ts_rate, ts_cnt_max, i_stamp

    if (.not. timing) return
    i_stamp = get_stamp(label)
    call system_clock(ts_cnt, ts_rate, ts_cnt_max)
    timestamps(i_stamp)%cnt = timestamps(i_stamp)%cnt - ts_cnt
end subroutine

subroutine stop_clock(label)
    character(len=*), intent(in) :: label

    integer :: ts_cnt, ts_rate, ts_cnt_max, i_stamp

    if (.not. timing) return
    i_stamp = get_stamp(label)
    call system_clock(ts_cnt, ts_rate, ts_cnt_max)
    timestamps(i_stamp)%cnt = timestamps(i_Stamp)%cnt + ts_cnt
end subroutine

subroutine errore(label, str, dummy_i)
    use localorb_io, only: OL_high, localorb_info

    character(len=*), intent(in) :: label, str
    integer, intent(in) :: dummy_i

    call localorb_info("*** Error: "//label//": "//str, priority=OL_high)
end subroutine

subroutine sync_sum_dble_(x, dummy)
    real(8), intent(inout) :: x
    integer, intent(in) :: dummy

    call sync_real_number(x)
end subroutine

subroutine sync_sum_vector_dble_(x, dummy)
    real(8), intent(inout) :: x(:)
    integer, intent(in) :: dummy

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_vector_int_(x, comm)
    integer, intent(inout) :: x(:)
    integer, intent(in) :: comm

    call sync_vector_integer(x, size(x))
end subroutine

subroutine sync_sum_matrix_dble_(x, dummy)
    real(8), intent(inout) :: x(:, :)
    integer, intent(in) :: dummy

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_3d_dble_(x, dummy)
    real(8), intent(inout) :: x(:, :, :)
    integer, intent(in) :: dummy

    call sync_vector(x, size(x))
end subroutine

subroutine sync_sum_4d_dble_(x, dummy)
    real(8), intent(inout) :: x(:, :, :, :)
    integer, intent(in) :: dummy

    call sync_vector(x, size(x))
end subroutine

end module
