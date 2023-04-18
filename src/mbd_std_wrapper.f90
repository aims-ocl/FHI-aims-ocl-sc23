module mbd_std_wrapper

use constants, only: pi
use localorb_io, only: OL_norm, OL_high, localorb_info
use dimensions, only: n_atoms, n_occ_atoms, n_species, n_full_points, n_periodic
use mbdvdw_module, only: &
    mbd_std_cleanup => mbdvdw_cleanup, mbd_first_step, &
    mbd_scf_converged => mbd_conv_elec
use mbdvdw_interface_module, only: &
    mbd_self_consistent => vdw_self_consistent

implicit none

private

public :: &
    mbd_std_parse, mbd_std_initialize, &
    mbd_std_calculate, mbd_std_potential, &
    mbd_std_cleanup, mbd_std_finalize, &
    mbd_self_consistent, mbd_scf_converged, mbd_first_step

integer, parameter :: n_flags = 10
character(len=40), public :: mbd_std_flags(n_flags)
logical, public :: run_numerical_pulay_forces = .false.

contains

integer function mbd_std_parse() result(retcode)
    use hirshfeld, only: grid_out, use_pulay
    use runtime_choices, only: flag_xc, n_k_points_xyz
    use mbdvdw_interface_module, only: &
        mbd_vdw_n_quad_pts, mbd_vdw_supercell, &
        use_old_hirshfeld, vdw_self_consistent, &
        use_hirshfeld_deriv, mbd_vdw_beta, timing, &
        mbd_vdw_recip, mbd_vdw_ewald, mbd_vdw_low_dim, &
        mbd_vdw_vacuum, mbd_vdw_kgrid, mbd_vdw_verbosity, &
        mbd_vdw_isolated, mbd_max_negative_eigvals

    integer :: stat, idx, i_flag
    character(len=len(mbd_std_flags(1))) :: flag
    character(len=20) :: key, val

    retcode = 1
    do i_flag = 1, n_flags
        flag = mbd_std_flags(i_flag)
        if (len(trim(flag)) == 0) exit
        if (flag(1:1) == '#') exit
        idx = index(flag, '=')
        if (idx == 0) then
            select case (lower(flag))
                case ('hirshfeld')
                case ('numerical_hirsh_deriv')
                case ('numerical_hirsh_pulay_deriv')
                    run_numerical_pulay_forces = .true.
                case default
                    call print_error("MBD: Unknown flag: " // trim(flag))
                    return
            end select
        else
            key = flag(1:idx-1)
            val = replace(flag(idx+1:), ':', ' ')
            stat = 0
            select case (lower(key))
                case ('supercell')
                    read (val, *, iostat=stat) mbd_vdw_supercell
                case ('k_grid')
                    read (val, *, iostat=stat) mbd_vdw_kgrid
                    mbd_vdw_recip = .true.
                case ('self_consistent')
                    read (val, *, iostat=stat) vdw_self_consistent
                case ('timing')
                    read (val, *, iostat=stat) timing
                case ('verbosity')
                    read (val, *, iostat=stat) mbd_vdw_verbosity
                case ('pulay!')
                    read (val, *, iostat=stat) use_pulay
                case ('hirshfeld_deriv!')
                    read (val, *, iostat=stat) use_hirshfeld_deriv
                case ('grid_out')
                    read (val, *, iostat=stat) grid_out
                case ('old_hirshfeld')
                    read (val, *, iostat=stat) use_old_hirshfeld
                case ('beta')
                    read (val, *, iostat=stat) mbd_vdw_beta
                case ('max_negative_eigvals')
                    read (val, *, iostat=stat) mbd_max_negative_eigvals
                case ('omega_grid')
                    read (val, *, iostat=stat) mbd_vdw_n_quad_pts
                case ('ewald')
                    read (val, *, iostat=stat) mbd_vdw_ewald
                case ('vacuum')
                    read (val, *, iostat=stat) mbd_vdw_vacuum
                    mbd_vdw_low_dim = any(mbd_vdw_vacuum)
                case default
                    call print_error("MBD: Unknown key: " // trim(key))
                    return
            end select
            if (stat /= 0) then
                call print_error('MBD: Invalid value for keyword ' // trim(key))
                return
            end if
        end if ! idx == 0
    end do ! i_flag
    if (n_periodic == 0) then
        if (mbd_vdw_recip .or. mbd_vdw_supercell > 0.d0) then
            call print_error( &
                'MBD: Non-periodic system, but periodic treatment specified' &
            )
            return
        end if
    end if
    if (mbd_vdw_beta <= 0.d0) then
        select case (flag_xc)
            case (1)  ! pbe0
                mbd_vdw_beta = 0.85d0
            case (6)  ! pbe
                mbd_vdw_beta = 0.83d0
            case (7)  ! hse
                mbd_vdw_beta = 0.85d0
            case(11)  ! rpbe HACK -- assuming the PBE holds
                mbd_vdw_beta = 0.83d0
            case default
                call print_error('MBD: Unknown XC functional and no beta specified')
                return
        end select
    end if
    if (n_periodic > 0) then
        if (.not. (mbd_vdw_recip .or. mbd_vdw_supercell > 0.d0)) then
            call print_log('MBD: Using the same k-grid sampling as DFT')
            mbd_vdw_kgrid = n_k_points_xyz
            mbd_vdw_recip = .true.
        end if
        mbd_vdw_isolated = .false.
    end if
    retcode = 0
end function

subroutine print_mbd_timing()
    use mbdvdw_interface_module, only: timing, n_timestamps, timestamps

    integer :: i_stamp
    character(len=200) :: strbuf
    integer :: ts_cnt, ts_rate, ts_cnt_max

    if (.not. timing) return
    call system_clock(ts_cnt, ts_rate, ts_cnt_max)
    do i_stamp = 1, n_timestamps
        if (timestamps(i_stamp)%label == '') exit
        write (strbuf, *) &
            timestamps(i_stamp)%label, dble(timestamps(i_stamp)%cnt)/ts_rate
        call localorb_info(strbuf)
    end do
end subroutine

subroutine mbd_std_calculate(rho, en_vdw, vdw_forces, vdw_stress)
    use geometry, only: coords, lattice_vector
    use xml_write, only: xml_open, xml_close, xml_elem, tostr
    use timing, only: start_timer, stop_timer, output_timer, output_timeheader
    use hirshfeld, only: numerical_hirshfeld_deriv, run_hirshfeld
    use mbdvdw_module, only: mbdvdw_calculate, Embdvdw, Fmbdvdw, Hmbdvdw
    use mbdvdw_interface_module, only: vefftsvdw
    use mpi_tasks, only: aims_stop_coll

    real(8), intent(in) :: rho(:, :)
    real(8), intent(out), optional :: en_vdw
    real(8), intent(out), optional :: vdw_forces(3, n_atoms)
    real(8), intent(out), optional :: vdw_stress(3, 3)

    integer :: i_flag
    real(8) :: timestamps(4)
    character(len=len(mbd_std_flags(1))) :: flag

    call print_log("Many-Body Dispersion (MBD@rsSCS) energy")
    call xml_open('mbd')
    call mbdvdw_calculate(coords, rho, lattice_vector)
    if (Embdvdw >= 0.d0) then
        call aims_stop_coll('MBD energy is zero', 'mbd_std_calculate')
    end if
    if (present(en_vdw)) en_vdw = Embdvdw
    if (present(vdw_forces)) then
        vdw_forces(:, :n_occ_atoms) = transpose(Fmbdvdw)
        vdw_forces(:, n_occ_atoms+1:) = 0.d0
    end if
    if (present(vdw_stress)) then
        vdw_stress = matmul(lattice_vector, transpose(Hmbdvdw))
    end if
    call xml_elem("energy", en_vdw, name="MBD@rsSCS", unit="Ha")
    do i_flag = 1, n_flags
        flag = mbd_std_flags(i_flag)
        if (len(trim(flag)) == 0) exit
        if (index(flag, '=') > 0) cycle
        call print_log('Evaluating ' // trim(flag) // '...')
        call start_timer(timestamps)
        select case (lower(flag))
            case ('hirshfeld')
                call run_hirshfeld(rho, vefftsvdw)
            case ('numerical_hirsh_deriv')
                call numerical_hirshfeld_deriv(rho)
            case ('numerical_hirsh_pulay_deriv')
                call print_log('Already preevaluated.')
            case default
                call print_error("MBD: Unknown flag: " // trim(flag))
                cycle
        end select
        call stop_timer(timestamps)
        call output_timer(trim(flag) // " evaluated", timestamps)
        call print_mbd_timing()
    end do
    call xml_close()
end subroutine

subroutine mbd_std_potential(rho, en_vdw, vdw_potential)
    use geometry, only: coords, lattice_vector
    use xml_write, only: xml_open, xml_close, xml_elem, tostr
    use timing, only: start_timer, stop_timer, output_timer, output_timeheader
    use mbdvdw_module, only: mbdvdw_calculate, Embdvdw, Umbdvdw

    real(8), intent(in) :: rho(:, :)
    real(8), intent(out) :: en_vdw
    real(8), intent(out) :: vdw_potential(n_full_points)

    call print_log('Evaluating MBD@rsSCS xc potential...')
    call mbdvdw_calculate(coords, rho, lattice_vector)
    en_vdw = Embdvdw
    vdw_potential = Umbdvdw
end subroutine

subroutine mbd_std_initialize()
    use geometry, only: lattice_vector
    use dimensions, only: n_full_points
    use mbdvdw_module, only: mbdvdw_initialize
    use mbdvdw_interface_module, only: &
        vefftsvdw, vfree, dveffdr, dveffdh, mbd_vdw_forces, dfftp

    allocate (vefftsvdw(n_atoms))
    allocate (vfree(n_species))
    if (mbd_vdw_forces) then
        allocate (dveffdr(n_atoms, n_atoms, 3))
        allocate (dveffdh(n_atoms, 3, 3))
    end if
    dfftp%nnr = n_full_points
    call mbdvdw_initialize(lattice_vector)
end subroutine

subroutine mbd_std_finalize()
    use mbdvdw_module, only: mbdvdw_finalize, mbdvdw_initialize
    use mbdvdw_interface_module, only: &
        vefftsvdw, vfree, dveffdr, dveffdh, mbd_vdw_forces

    call mbdvdw_finalize()
    if (allocated(vefftsvdw)) deallocate (vefftsvdw)
    if (allocated(vfree)) deallocate (vfree)
    if (mbd_vdw_forces) then
        if (allocated(dveffdr)) deallocate (dveffdr)
        if (allocated(dveffdh)) deallocate (dveffdh)
    end if
end subroutine

function lower(str)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower

    integer :: i

    do i = 1, len(str)
        select case (str(i:i))
            case ('A':'Z')
                lower(i:i) = achar(iachar(str(i:i))+32)
            case default
                lower(i:i) = str(i:i)
        end select
    end do
end function

function inverted3x3(A) result(A_inv)
    use mbdvdw_module, only: inv3x3_mat

    real(8), intent(in) :: A(:, :)
    real(8) :: A_inv(size(A, 1), size(A, 2))

    call inv3x3_mat(A, A_inv)
end function

function replace(str, from, to)
    character(len=*), intent(inout) :: str
    character(len=1), intent(in) :: from, to
    character(len=len(str)) :: replace

    integer :: i

    replace = str
    do i = 1, len(str)
        if (replace(i:i) == from) replace(i:i) = to
    end do
end function

subroutine print_log (str, mute)
    character(len=*), intent(in) :: str
    logical, optional :: mute

    if (present(mute)) then
        if (mute) return
    end if
    call localorb_info("  | " // str, priority=OL_norm)
end subroutine

subroutine print_warning (str)
    character(len=*), intent(in) :: str

    call localorb_info("*** Warning: " // str, priority=OL_high)
end subroutine

subroutine print_error (str)
    character(len=*), intent(in) :: str

    call localorb_info("*** Error: " // str, priority=OL_high)
end subroutine

end module
