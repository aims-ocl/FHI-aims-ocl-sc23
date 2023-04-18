module mbd_dev_wrapper

use localorb_io, only: OL_norm, OL_high, localorb_info
use dimensions, only: n_atoms => n_occ_atoms, n_full_points, n_periodic
use mbd_interface, only: print_error, print_warning, print_log

implicit none

private

integer, parameter :: n_flags = 10
character(len=40), public :: mbd_dev_flags(n_flags)

real(8) :: &
    ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
    mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta
real(8) :: beta_param = 0.d0
integer :: n_omega_grid = 20
real(8) :: delta = 1.d-4
integer :: k_grid(3), supercell(3)
character(len=3) :: mode
logical :: &
    do_supercell = .false., do_reciprocal = .false., &
    output_modes = .false., do_mbdrsscs_forces = .false., &
    do_ts_forces = .false., do_mbdrsscs = .false., &
    do_mbd = .false., do_ts = .false., use_scalapack = .false.

public :: mbd_dev_parse, mbd_dev_verify, mbd_dev_calculate

contains

integer function mbd_dev_parse() result(info)
    use hirshfeld, only: grid_out, grid_atoms
    use mbd_dev, only: &
        param_vacuum_axis, param_ts_energy_accuracy, &
        param_dipole_low_dim_cutoff, bohr, &
        param_ts_cutoff_radius, &
        param_zero_negative_eigs

    integer :: iostat, idx, i_flag
    character(len=20) :: key, val
    character(len=len(mbd_dev_flags(1))) :: flag

    info = 1
    do i_flag = 1, n_flags
        flag = mbd_dev_flags(i_flag)
        if (len(trim(flag)) == 0) exit
        idx = index(flag, '=')
        if (idx == 0) then
            select case (lower(flag))
                case ('mbd@rsscs')
                    do_mbdrsscs = .true.
                    do_mbd = .true.
                case ('mbd@ts')
                    do_mbd = .true.
                case ('ts')
                    do_ts = .true.
                case ('mbd@rsscs_forces')
                    do_mbdrsscs_forces = .true.
                case ('ts_forces')
                    do_ts_forces = .true.
                case ('rsscs', 'scs', 'ts@scs')
                case default
                    call print_error("MBD: Unknown flag: "//trim(flag))
                    return
            end select
        else
            key = flag(1:idx-1)
            val = replace(flag(idx+1:), ':', ' ')
            iostat = 0
            select case (lower(key))
                case ('supercell')
                    read (val, *, iostat=iostat) supercell
                    do_supercell = .true.
                case ('ts_ene_acc')
                    read (val, *, iostat=iostat) param_ts_energy_accuracy
                case ('ts_cutoff')
                    read (val, *, iostat=iostat) param_ts_cutoff_radius
                    param_ts_cutoff_radius = param_ts_cutoff_radius/bohr
                case ('vacuum')
                    read (val, *, iostat=iostat) param_vacuum_axis
                case ('scalapack')
                    read (val, *, iostat=iostat) use_scalapack
                    if (use_scalapack) call print_log('Will use ScaLAPACK for MBD')
                case ('delta')
                    read (val, *, iostat=iostat) delta
                case ('lowdim_cutoff')
                    read (val, *, iostat=iostat) param_dipole_low_dim_cutoff
                    param_dipole_low_dim_cutoff = param_dipole_low_dim_cutoff/bohr
                case ('zero_negative')
                    read (val, *, iostat=iostat) param_zero_negative_eigs
                case ('k_grid')
                    read (val, *, iostat=iostat) k_grid
                    do_reciprocal = .true.
                case ('grid_out')
                    read (val, *, iostat=iostat) grid_out
                case ('grid_atoms')
                    read (val, *, iostat=iostat) grid_atoms
                    if (iostat == -1) iostat = 0
                case ('modes') 
                    read (val, *, iostat=iostat) output_modes
                case ('beta') 
                    read (val, *, iostat=iostat) beta_param
                case ('omega_grid') 
                    read (val, *, iostat=iostat) n_omega_grid
                    if (n_omega_grid > 20) then
                        call print_error('MBD: omega_grid cannot be larger than 20')
                        return
                    end if
                case default
                    call print_error("MBD: Unknown key: " // trim(key))
                    return
            end select
            if (iostat /= 0) then
                call print_error('MBD: Invalid value for keyword '//trim(key))
                return
            end if
        end if ! idx == 0
    end do ! i_flag
    info = 0
end function

integer function mbd_dev_verify() result(info)
    use dimensions, only: use_forces
    use mpi_tasks, only: n_tasks
    use runtime_choices, only: flag_xc
    use mbd_dev, only: get_damping_parameters

    character(len=15) :: xc

    info = 1
    mode = ''
    select case (flag_xc)
        case (1)
            xc = "pbe0"
        case (6)
            xc = "pbe"
        case (7)
            xc = "hse"
        case default
            xc = ''
            if (beta_param == 0.d0 .or. do_ts .or. do_ts_forces) then
                call print_error('MBD: Unknown XC functional')
                return
            end if
    endselect
    call get_damping_parameters( &
        xc, ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta)
    if (beta_param /= 0.d0) mbd_rsscs_beta = beta_param
    if (n_tasks > 0) mode(1:1) = 'P'
    if (do_reciprocal) then
        mode(3:3) = 'R'
    end if
    if (n_periodic > 0) then
        mode(2:2) = 'C'
        if (do_mbd .and. .not. (do_reciprocal .or. do_supercell)) then
            call print_error( &
                'MBD: Periodic system, but no periodic method specified')
            return
        end if
    else
        if (do_reciprocal .or. do_supercell) then
            call print_error( &
                'MBD: Non-periodic system, but periodic method specified')
            return
        end if
    end if
    if ((do_mbdrsscs_forces .or. do_ts_forces) .and. .not. use_forces) then
        call print_error("MBD: Turn on forces in AIMS as well")
        return
    end if
    if (use_forces .and. .not. (do_mbdrsscs_forces .or. do_ts_forces)) then
        call print_error( &
            'MBD: Specify calculation of forces in MBD as well')
        return
    end if
    if (do_mbdrsscs_forces .and. do_ts_forces) then
        call print_error( &
            'MBD: Cannot do both MBD@rsSCS and TS forces')
        return
    end if
    info = 0
end function

subroutine mbd_dev_calculate(volume_ratio, en_vdw, vdw_forces, vdw_stress)
    use constants, only: pi
    use vdw_correction, only: get_vdw_param, vdw_hirshfeld_data_external, &
        vdw_hirshfeld_C6, vdw_hirshfeld_alpha, vdw_hirshfeld_R0
    use geometry, only: species, coords, lattice_vector
    use dimensions, only: n_atoms => n_occ_atoms, n_periodic
    use species_data, only: species_element
    use mpi_tasks, only: n_tasks, myid
    use xml_write, only: xml_open, xml_close, xml_elem, tostr
    use timing, only: start_timer, stop_timer, output_timer, output_timeheader
    use mbd_dev, only: &
        init_grid, destroy_grid, mbd_myid => my_task, mbd_n_tasks => n_tasks, &
        omega_eff, is_in, blanked, get_ts_energy, get_mbd_energy, &
        make_k_grid, make_g_grid, param_zero_negative_eigs
    use mbd_interface, only: mbd_pi => pi

    implicit none

    real(8), intent(in) :: volume_ratio(n_atoms)
    real(8), intent(out), optional :: en_vdw
    real(8), intent(out), optional :: vdw_forces(3, n_atoms)
    real(8), intent(out), optional :: vdw_stress(3, 3)

    integer :: i_flag, i_atom, i_specie
    real(8), dimension(n_atoms) :: alpha_0_free, C6_free, R_vdw_free
    real(8), dimension(n_atoms) :: alpha_0_scaled, C6_scaled, R_vdw_scaled
    real(8), dimension(n_atoms) :: alpha_0, C6, R_vdw
    real(8), allocatable :: k_points(:, :)
    real(8) :: timestamps(4)
    character(len=len(mbd_dev_flags(1))) :: flag
    real(8) :: ene
    real(8) :: forces(n_atoms, 3)
    real(8) :: stress(3, 3)
        
    character*300 :: info_str

    call print_log("Experimental MBD module with reciprocal-space implementation")
    call output_timeheader("2X", "| Timing") 
    call xml_open('mbd')

    mbd_myid = myid
    mbd_n_tasks = n_tasks
    mbd_pi = pi
    call init_grid(n_omega_grid)

    do i_atom = 1, n_atoms
        i_specie = species(i_atom)
        if (vdw_hirshfeld_data_external(i_specie)) then
            C6_free(i_atom) = vdw_hirshfeld_C6(i_specie)
            alpha_0_free(i_atom) = vdw_hirshfeld_alpha(i_specie)
            R_vdw_free(i_atom) = vdw_hirshfeld_R0(i_specie)
        else
            call get_vdw_param( &
                species_element(species(i_atom)), 0.d0, &
                C6_free(i_atom), alpha_0_free(i_atom), R_vdw_free(i_atom))
        end if
    enddo
    alpha_0_scaled = volume_ratio*alpha_0_free
    C6_scaled = (volume_ratio**2)*C6_free
    R_vdw_scaled = (volume_ratio**(1.d0/3))*R_vdw_free
    call xml_open('free')
    call xml_elem("alpha_0", alpha_0_free)
    call xml_elem("C6", C6_free)
    call xml_elem("omega", omega_eff(C6_free, alpha_0_free))
    call xml_elem("R_vdw", R_vdw_free)
    call xml_close()
    call xml_elem("volume_ratio", volume_ratio)
    call xml_open('TS')
    call xml_elem("alpha_0", alpha_0_scaled)
    call xml_elem("C6", C6_scaled)
    call xml_elem("omega", omega_eff(C6_scaled, alpha_0_scaled))
    call xml_elem("R_vdw", R_vdw_scaled)
    call xml_close()

    if (present(en_vdw)) en_vdw = 0.d0
    if (present(vdw_forces)) vdw_forces(:, :) = 0.d0
    if (present(vdw_stress)) vdw_stress(:, :) = 0.d0
    do i_flag = 1, n_flags
        flag = mbd_dev_flags(i_flag)
        if (len(trim(flag)) == 0) exit
        if (index(flag, '=') > 0) cycle
        call print_log('Evaluating '//trim(flag)//'...')
        ene = 0.d0
        call start_timer(timestamps)
        select case (lower(flag))
            case ('ts')
                ene = get_ts_energy( &
                    mode, 'fermi', transpose(coords(:, :n_atoms)), C6_scaled, &
                    alpha_0_scaled, R_vdw_scaled, s_R=ts_s_r, d=ts_d, &
                    unit_cell=transpose(lattice_vector))
                if (present(en_vdw)) en_vdw = ene
            case ('ts_forces')
                call ts_forces(forces, stress, &
                    transpose(coords(:, :n_atoms)), transpose(lattice_vector), &
                    alpha_0_scaled, C6_scaled, &
                    R_vdw_scaled, s_R=ts_s_r, d=ts_d, delta=delta)
                if (present(vdw_forces)) vdw_forces = transpose(forces)
                if (present(vdw_stress)) vdw_stress = stress
            case ('ts@scs')
                ene = ts_scs( &
                    transpose(coords(:, :n_atoms)), transpose(lattice_vector), &
                    alpha_0_scaled, C6_scaled, &
                    R_vdw_scaled, s_r=ts_s_r, d=ts_d)
            case ('mbd@rsscs')
                ene = mbd_rsscs( &
                    transpose(coords(:, :n_atoms)), transpose(lattice_vector), &
                    alpha_0_scaled, C6_scaled, &
                    R_vdw_scaled, beta=mbd_rsscs_beta, a=mbd_rsscs_a)
                if (present(en_vdw)) en_vdw = ene
            case ('mbd@ts')
                if (is_in('R', mode)) then
                    allocate (k_points(product(k_grid), 3))
                    k_points = make_k_grid( &
                        make_g_grid(k_grid(1), k_grid(2), k_grid(3)), &
                        transpose(lattice_vector)) 
                end if
                ene = get_mbd_energy( &
                    mode, 'fermi,dip', transpose(coords(:, :n_atoms)), alpha_0_scaled, &
                    omega_eff(C6_scaled, alpha_0_scaled), &
                    supercell, k_points, transpose(lattice_vector), &
                    R_vdw=R_vdw_scaled, beta=mbd_ts_fermi_beta, a=mbd_ts_a)
                if (allocated(k_points)) deallocate (k_points)
            case ('mbd@rsscs_forces')
                call mbd_rsscs_forces(forces, stress, &
                    transpose(coords(:, :n_atoms)), transpose(lattice_vector), &
                    alpha_0_scaled, C6_scaled, &
                    R_vdw_scaled, beta=mbd_rsscs_beta, a=mbd_rsscs_a, &
                    delta=delta)
                if (present(vdw_forces)) vdw_forces = transpose(forces)
                if (present(vdw_stress)) vdw_stress = stress
            case ('rsscs')
                call rsscs( &
                    transpose(coords(:, :n_atoms)), transpose(lattice_vector), &
                    alpha_0_scaled, C6_scaled, &
                    alpha_0, C6, R_vdw, &
                    R_vdw_scaled, beta=mbd_rsscs_beta, a=mbd_rsscs_a)
                call xml_open('rsSCS')
                call xml_elem("alpha_0", alpha_0)
                call xml_elem("C6", C6)
                call xml_elem("omega", omega_eff(C6, alpha_0))
                call xml_elem("R_vdw", R_vdw)
                call xml_close()
            case ('scs')
                call scs( &
                    transpose(coords(:, :n_atoms)), transpose(lattice_vector), &
                    alpha_0_scaled, C6_scaled, &
                    alpha_0, C6, R_vdw)
                call xml_open('SCS')
                call xml_elem("alpha_0", alpha_0)
                call xml_elem("C6", C6)
                call xml_elem("omega", omega_eff(C6, alpha_0))
                call xml_elem("R_vdw", R_vdw)
                call xml_close()
            case default
                call print_error("MBD: Unknown flag: "//trim(flag))
                cycle
        end select
        call stop_timer(timestamps)
        call output_timer(trim(flag)//" evaluated", timestamps)
        if (ene > 0.d0 .or. ene < 0.d0) then
            ! WPH:  For whatever reason, the PGI compilers on Titan did not like
            !       the tostr(ene) call in this one location and kept returning a
            !       "Ambiguous interfaces for generic procedure tostr" error message.  I
            !       can't for the life of me see what the ambiguity is, so I've
            !       rewritten the print_log statement to replicate the original behavior.
            write(info_str, '(A,A,g50.17e3)') trim(flag), " energy: ", ene
            call print_log(info_str)
            call xml_elem("energy", ene, name=trim(flag), unit="Ha")
        end if
    end do

    call destroy_grid()

    call xml_close()

contains

    real(8) function ts_scs(xyz, uc, alpha_0, C6, R_vdw, s_R, d) result(ene)
        use mbd_dev, only: &
            run_scs, alpha_dynamic_ts_all, get_C6_from_alpha, &
            get_single_mbd_energy, omega_eff, make_k_grid, make_g_grid, &
            get_reciprocal_mbd_energy, get_supercell_mbd_energy

        real(8), intent(in) :: xyz(n_atoms, 3)
        real(8), intent(in) :: uc(3, 3)
        real(8), dimension(n_atoms), intent(in) :: alpha_0, C6, R_vdw
        real(8), intent(in) :: s_r, d

        real(8), dimension(0:n_omega_grid, n_atoms) :: &
            alpha_dyn, alpha_dyn_scs
        real(8), dimension(n_atoms) :: C6_scs, R_vdw_scs

        alpha_dyn = alpha_dynamic_ts_all('C', n_omega_grid, alpha_0, C6=C6)
        alpha_dyn_scs = run_scs( &
            blanked('R', mode), 'dip,gg', xyz, alpha_dyn, unit_cell=uc)
        C6_scs = get_C6_from_alpha(alpha_dyn_scs)
        R_vdw_scs = R_vdw*(alpha_dyn_scs(0, :)/alpha_dyn(0, :))**(1.d0/3)
        ene = get_ts_energy( &
            mode, 'fermi', xyz, C6_scs, &
            alpha_dyn_scs(0, :), R_vdw_scs, s_R=s_R, d=d, &
            unit_cell=uc)
    end function

    real(8) function mbd_rsscs(xyz, uc, alpha_0, C6, R_vdw, beta, a) result(ene)
        use mbd_dev, only: &
            run_scs, alpha_dynamic_ts_all, get_C6_from_alpha, &
            get_single_mbd_energy, omega_eff, make_k_grid, make_g_grid, &
            get_reciprocal_mbd_energy, get_supercell_mbd_energy, &
            run_scs_s, get_single_mbd_energy_s, &
            get_reciprocal_mbd_energy_s, get_supercell_mbd_energy_s

        real(8), intent(in) :: xyz(n_atoms, 3)
        real(8), intent(in) :: uc(3, 3)
        real(8), dimension(n_atoms), intent(in) :: alpha_0, C6, R_vdw
        real(8), intent(in) :: beta, a

        real(8), dimension(0:n_omega_grid, n_atoms) :: &
            alpha_dyn, alpha_dyn_rsscs
        real(8), dimension(n_atoms) :: C6_rsscs, R_vdw_rsscs

        real(8), allocatable :: bands(:, :)
        complex(8), allocatable :: modes(:, :, :)

        alpha_dyn = alpha_dynamic_ts_all('C', n_omega_grid, alpha_0, C6=C6)
        if (use_scalapack) then
            alpha_dyn_rsscs = run_scs_s( &
                blanked('R', mode), 'fermi,dip,gg', xyz, alpha_dyn, &
                R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc)
        else
            alpha_dyn_rsscs = run_scs( &
                blanked('R', mode), 'fermi,dip,gg', xyz, alpha_dyn, &
                R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc)
        end if
        C6_rsscs = get_C6_from_alpha(alpha_dyn_rsscs)
        R_vdw_rsscs = R_vdw*(alpha_dyn_rsscs(0, :)/alpha_dyn(0, :))**(1.d0/3)
        if (n_periodic == 0) then
            if (use_scalapack) then
                ene = get_single_mbd_energy_s( &
                    mode, 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                    omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                    R_vdw=R_vdw_rsscs, beta=beta, a=a)
            else
                ene = get_single_mbd_energy( &
                    mode, 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                    omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                    R_vdw=R_vdw_rsscs, beta=beta, a=a)
            end if
        else
            if (do_reciprocal) then
                if (use_scalapack) then
                    ene = get_reciprocal_mbd_energy_s( &
                        mode//'R', 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                        omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                        make_k_grid( &
                            make_g_grid(k_grid(1), k_grid(2), k_grid(3)), uc), &
                        uc, &
                        R_vdw=R_vdw_rsscs, beta=beta, a=a)
                else if (output_modes) then
                    allocate (bands(product(k_grid), 3*n_atoms))
                    allocate (modes(product(k_grid), 3*n_atoms, 3*n_atoms))
                    ene = get_reciprocal_mbd_energy( &
                        mode//'REV', 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                        omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                        make_k_grid( &
                            make_g_grid(k_grid(1), k_grid(2), k_grid(3)), uc), &
                        uc, &
                        R_vdw=R_vdw_rsscs, beta=beta, a=a, &
                        mode_enes=bands, &
                        modes=modes)
                    call xml_elem('bands', bands)
                    deallocate (bands)
                    deallocate (modes)
                else
                    ene = get_reciprocal_mbd_energy( &
                        mode//'R', 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                        omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                        make_k_grid( &
                            make_g_grid(k_grid(1), k_grid(2), k_grid(3)), uc), &
                        uc, &
                        R_vdw=R_vdw_rsscs, beta=beta, a=a)
                end if
            else if (do_supercell) then
                if (use_scalapack) then
                    ene = get_supercell_mbd_energy_s( &
                        mode, 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                        omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), uc, &
                        supercell, &
                        R_vdw=R_vdw_rsscs, beta=beta, a=a)
                else
                    ene = get_supercell_mbd_energy( &
                        mode, 'fermi,dip', xyz, alpha_dyn_rsscs(0, :), &
                        omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), uc, &
                        supercell, &
                        R_vdw=R_vdw_rsscs, beta=beta, a=a)
                end if
            end if
        end if
    end function

    subroutine ts_forces(forces, stress, &
            xyz, uc, alpha_0, C6, R_vdw, delta, s_R, d)
        use mbd_dev, only: get_ts_energy
        use mbd_interface, only: sync_sum
        use runtime_choices, only: use_analytical_stress

        real(8), intent(in) :: xyz(n_atoms, 3)
        real(8), intent(in) :: uc(3, 3)
        real(8), dimension(n_atoms), intent(in) :: alpha_0, C6, R_vdw
        real(8), intent(in) :: s_R, d, delta
        real(8), intent(out) :: forces(n_atoms, 3)
        real(8), intent(out), optional :: stress(3, 3)

        real(8) :: &
            xyz_diff(n_atoms, 3), ene(-1:1), uc_diff(3, 3), strain(3, 3)
        integer :: i_atom, i_xyz, i_step, i, j_xyz

        forces(:, :) = 0.d0
        do i_atom = 1, n_atoms
            do i_xyz = 1, 3
                if (myid /= modulo(3*(i_atom-1)+i_xyz-1, n_tasks)) cycle
                do i_step = -1, 1, 2
                    xyz_diff = xyz
                    xyz_diff(i_atom, i_xyz) = xyz_diff(i_atom, i_xyz)+i_step*delta
                    ene(i_step) = get_ts_energy( &
                        blanked('P', mode)//'M', &
                        'fermi', xyz_diff, C6, &
                        alpha_0, R_vdw, s_R=s_R, d=d, &
                        unit_cell=transpose(lattice_vector))
                end do
                forces(i_atom, i_xyz) = -(ene(1)-ene(-1))/(2*delta)
            end do
        end do
        call sync_sum(forces)
        if (use_analytical_stress .and. present(stress)) then
            stress(:, :) = 0.d0
            do i_xyz = 1, 3
                do j_xyz = 1, 3
                    do i_step = -1, 1, 2
                        strain(:, :) = 0.d0
                        forall (i = 1:3) strain(i, i) = 1.d0
                        strain(i_xyz, j_xyz) = strain(i_xyz, j_xyz) &
                            +i_step*delta
                        xyz_diff = matmul(xyz, transpose(strain))
                        uc_diff = matmul(uc, transpose(strain))
                        ene(i_step) = get_ts_energy( &
                            mode//'M', &
                            'fermi', xyz_diff, C6, &
                            alpha_0, R_vdw, s_R=s_R, d=d, &
                            unit_cell=uc_diff)
                    end do
                    stress(i_xyz, j_xyz) = (ene(1)-ene(-1))/(2*delta)
                end do
            end do
        end if
    end subroutine

    subroutine mbd_rsscs_forces(forces, stress, &
            xyz, uc, alpha_0, C6, R_vdw, delta, beta, a)
        use mbd_dev, only: &
            run_scs, alpha_dynamic_ts_all, get_C6_from_alpha, &
            get_mbd_energy, omega_eff, make_k_grid, make_g_grid, &
            run_scs_s, get_mbd_energy_s
        use mbd_interface, only: sync_sum
        use runtime_choices, only: use_analytical_stress

        real(8), intent(in) :: xyz(n_atoms, 3)
        real(8), intent(in) :: uc(3, 3)
        real(8), dimension(n_atoms), intent(in) :: alpha_0, C6, R_vdw
        real(8), intent(in) :: beta, a, delta
        real(8), intent(out) :: forces(n_atoms, 3)
        real(8), intent(out), optional :: stress(3, 3)

        real(8), dimension(0:n_omega_grid, n_atoms) :: &
            alpha_dyn, alpha_dyn_rsscs
        real(8), dimension(n_atoms) :: C6_rsscs, R_vdw_rsscs

        real(8) :: &
            xyz_diff(n_atoms, 3), ene(-1:1), uc_diff(3, 3), strain(3, 3)
        integer :: i_atom, i_xyz, i_step, i, j_xyz
        real(8), allocatable :: k_points(:, :)

        alpha_dyn = alpha_dynamic_ts_all('C', n_omega_grid, alpha_0, C6=C6)
        if (is_in('R', mode)) then
            allocate (k_points(product(k_grid), 3))
            k_points = make_k_grid(make_g_grid(k_grid(1), k_grid(2), k_grid(3)), uc)
        end if

        forces(:, :) = 0.d0
        do i_atom = 1, n_atoms
            do i_xyz = 1, 3
                if (.not. use_scalapack .and. myid /= modulo(3*(i_atom-1)+i_xyz-1, n_tasks)) cycle
                do i_step = -1, 1, 2
                    xyz_diff = xyz
                    xyz_diff(i_atom, i_xyz) = xyz_diff(i_atom, i_xyz)+i_step*delta
                    if (use_scalapack) then
                        alpha_dyn_rsscs = run_scs_s( &
                            blanked('PR', mode)//'M', &
                            'fermi,dip,gg', xyz_diff, alpha_dyn, &
                            R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc)
                    else
                        alpha_dyn_rsscs = run_scs( &
                            blanked('PR', mode)//'M', &
                            'fermi,dip,gg', xyz_diff, alpha_dyn, &
                            R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc)
                    end if
                    C6_rsscs = get_C6_from_alpha(alpha_dyn_rsscs)
                    R_vdw_rsscs = R_vdw*(alpha_dyn_rsscs(0, :) &
                        /alpha_dyn(0, :))**(1.d0/3)
                    if (use_scalapack) then
                        ene(i_step) = get_mbd_energy_s( &
                            blanked('P', mode)//'M', &
                            'fermi,dip', xyz_diff, alpha_dyn_rsscs(0, :), &
                            omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                            supercell, &
                            k_points, &
                            uc, &
                            R_vdw=R_vdw_rsscs, beta=beta, a=a)
                    else
                        ene(i_step) = get_mbd_energy( &
                            blanked('P', mode)//'M', &
                            'fermi,dip', xyz_diff, alpha_dyn_rsscs(0, :), &
                            omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                            supercell, &
                            k_points, &
                            uc, &
                            R_vdw=R_vdw_rsscs, beta=beta, a=a)
                    end if
                end do
                forces(i_atom, i_xyz) = -(ene(1)-ene(-1))/(2*delta)
            end do
        end do
        if (.not. use_scalapack) then
            call sync_sum(forces)
        end if
        if (use_analytical_stress .and. present(stress)) then
            stress(:, :) = 0.d0
            do i_xyz = 1, 3
                do j_xyz = 1, 3
                    do i_step = -1, 1, 2
                        strain(:, :) = 0.d0
                        forall (i = 1:3) strain(i, i) = 1.d0
                        strain(i_xyz, j_xyz) = strain(i_xyz, j_xyz) &
                            +i_step*delta
                        xyz_diff = matmul(xyz, transpose(strain))
                        uc_diff = matmul(uc, transpose(strain))
                        if (use_scalapack) then
                            alpha_dyn_rsscs = run_scs_s( &
                                blanked('R', mode)//'M', &
                                'fermi,dip,gg', xyz_diff, alpha_dyn, &
                                R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc_diff)
                        else
                            alpha_dyn_rsscs = run_scs( &
                                blanked('R', mode)//'M', &
                                'fermi,dip,gg', xyz_diff, alpha_dyn, &
                                R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc_diff)
                        end if
                        C6_rsscs = get_C6_from_alpha(alpha_dyn_rsscs)
                        R_vdw_rsscs = R_vdw*(alpha_dyn_rsscs(0, :) &
                            /alpha_dyn(0, :))**(1.d0/3)
                        if (is_in('R', mode)) then
                            k_points = make_k_grid( &
                                make_g_grid(k_grid(1), k_grid(2), k_grid(3)), &
                                uc_diff &
                            )
                        end if
                        if (use_scalapack) then
                            ene(i_step) = get_mbd_energy_s( &
                                mode//'M', &
                                'fermi,dip', xyz_diff, alpha_dyn_rsscs(0, :), &
                                omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                                supercell, &
                                k_points, &
                                uc_diff, &
                                R_vdw=R_vdw_rsscs, beta=beta, a=a)
                        else
                            ene(i_step) = get_mbd_energy( &
                                mode//'M', &
                                'fermi,dip', xyz_diff, alpha_dyn_rsscs(0, :), &
                                omega_eff(C6_rsscs, alpha_dyn_rsscs(0, :)), &
                                supercell, &
                                k_points, &
                                uc_diff, &
                                R_vdw=R_vdw_rsscs, beta=beta, a=a)
                        end if
                    end do
                    stress(i_xyz, j_xyz) = (ene(1)-ene(-1))/(2*delta)
                end do
            end do
        end if

        if (allocated(k_points)) deallocate (k_points)
    end subroutine

    subroutine rsscs( &
            xyz, uc, alpha_0, C6, &
            alpha_0_rsscs, C6_rsscs, R_vdw_rsscs, &
            R_vdw, beta, a)
        use mbd_dev, only: &
            run_scs, alpha_dynamic_ts_all, get_C6_from_alpha

        real(8), intent(in) :: xyz(n_atoms, 3)
        real(8), intent(in) :: uc(3, 3)
        real(8), dimension(n_atoms), intent(in) :: alpha_0, C6, R_vdw
        real(8), dimension(n_atoms), intent(out) :: &
            alpha_0_rsscs, C6_rsscs, R_vdw_rsscs
        real(8), intent(in) :: beta, a

        real(8), dimension(0:n_omega_grid, n_atoms) :: alpha_dyn, alpha_dyn_rsscs

        alpha_dyn = alpha_dynamic_ts_all('C', n_omega_grid, alpha_0, C6=C6)
        alpha_dyn_rsscs = run_scs( &
            blanked('R', mode), 'fermi,dip,gg', xyz, alpha_dyn, &
            R_vdw=R_vdw, beta=beta, a=a, unit_cell=uc)
        alpha_0_rsscs = alpha_dyn_rsscs(0, :)
        C6_rsscs = get_C6_from_alpha(alpha_dyn_rsscs)
        R_vdw_rsscs = R_vdw*(alpha_0_rsscs/alpha_dyn(0, :))**(1.d0/3)
    end subroutine

    subroutine scs( &
            xyz, uc, alpha_0, C6, &
            alpha_0_scs, C6_scs, R_vdw_scs)
        use mbd_dev, only: &
            run_scs, alpha_dynamic_ts_all, get_C6_from_alpha

        real(8), intent(in) :: xyz(n_atoms, 3)
        real(8), intent(in) :: uc(3, 3)
        real(8), dimension(n_atoms), intent(in) :: alpha_0, C6
        real(8), dimension(n_atoms), intent(out) :: &
            alpha_0_scs, C6_scs, R_vdw_scs

        real(8), dimension(0:n_omega_grid, n_atoms) :: alpha_dyn, alpha_dyn_scs

        alpha_dyn = alpha_dynamic_ts_all('C', n_omega_grid, alpha_0, C6=C6)
        alpha_dyn_scs = run_scs( &
            blanked('R', mode), 'dip,gg', xyz, alpha_dyn, unit_cell=uc)
        alpha_0_scs = alpha_dyn_scs(0, :)
        C6_scs = get_C6_from_alpha(alpha_dyn_scs)
        R_vdw_scs = R_vdw*(alpha_0_scs/alpha_dyn(0, :))**(1.d0/3)
    end subroutine

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

end module
