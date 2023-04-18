module mbd_dev

! modes:
! P: parallel (MPI)
! C: crystal
! R: reciprocal
! Q: do RPA
! E: get eigenvalues
! V: get eigenvectors
! O: get RPA orders
! L: scale dipole with lambda
! M: mute

use mbd_interface, only: &
    sync_sum, broadcast, print_error, print_warning, print_log, pi
use mbd_helper, only: &
    is_in, blanked
use elsi_wrapper, only: aims_elsi_stdevp
use mpi_tasks, only: &
    mpi_comm_global, myid
use tools, only: norm2

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!  This is an implementation of the Many-Body Dispersion (MBD)       !!
!!  dispersion model [Tkatchenko, et al. PRL 108, 236402 (2012)]      !!
!!  Including vdW(TS) [Tkatchenko, Scheffler PRL 102, 073005 (2009)]  !!
!!  and some advanced additional capabilities for development         !!
!!  based on code by J. Hermann (FHI) ammended by M. Stoehr (UL)      !!
!!                                                                    !!
!!  The code is structured as follows:                                !!
!!  . vdW(TS) energy                                                  !!
!!                                                                    !!
!!  . LAPACK/MPI-based routines to build dipole interaction tensor    !!
!!  . LAPACK/MPI-based routines for self-consistent screening         !!
!!  . LAPACK/MPI-based routines for MBD/RPA energy                    !!
!!                                                                    !!
!!  . additional tools (development)                                  !!
!!  . helper functions [LAPACK/MPI framework]                         !!
!!                                                                    !!
!!  . vdW(TS) energy (low memory version, less customizability)       !!
!!                                                                    !!
!!  . ScaLAPACK-based routines to build dipole interaction tensor     !!
!!  . ScaLAPACK-based routines for self-consistent screening          !!
!!  . ScaLAPACK-based routines for MBD energy (TODO?: RPA)            !!
!!  . helper functions [ScaLAPACK framework]                          !!
!!                                                                    !!
!!  . dipole-dipole coupling and damping                              !!
!!  . general operations/routines                                     !!
!!  . timing                                                          !!
!!                                                                    !!
!!  Available eigensolvers for CFDM -- in both frameworks             !!
!!  (controlled via module variable 'eigensolver' = ):                !!
!!  . 'qr   ': standard (P)DSYEV/(P)ZHEEV QR algorithm [default]      !!
!!  . 'dandc': Divide & Conquer algorithm (P)DSYEVD/(P)ZHEEVD         !!
!!  . 'mrrr ': MRRR algorithm (P)DSYEVR/(P)ZHEEVR                     !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8), parameter :: bohr = 0.529177249d0

real(8) :: &
    param_ts_energy_accuracy = 1.d-10, &
    param_ts_cutoff_radius = 50.d0/bohr, &
    param_dipole_low_dim_cutoff = 100.d0/bohr, &
    param_dipole_cutoff = 400.d0/bohr, &  ! used only when Ewald is off
    param_mayer_scaling = 1.d0, &
    param_ewald_real_cutoff_scaling = 1.d0, &
    param_ewald_rec_cutoff_scaling = 1.d0, &
    param_k_grid_shift = 0.5d0
logical :: &
    param_ewald_on = .true., &
    param_zero_negative_eigs = .false.
integer :: &
    param_mbd_nbody_max = 3, &
    param_rpa_order_max = 10
logical :: &
    param_vacuum_axis(3) = (/ .false., .false., .false. /)

integer :: n_grid_omega
real(8), allocatable :: omega_grid(:), omega_grid_w(:)

integer, parameter :: n_timestamps = 100
logical :: measure_time = .true.
integer :: timestamps(n_timestamps), ts_counts(n_timestamps)
integer :: ts_cnt, ts_rate, ts_cnt_max, ts_aid

integer :: my_task, n_tasks

integer, parameter :: minFileID = 20
integer, parameter :: maxFileID = 65535

! VY: Now always use ELPA solver in parallel cases
character(len=5) :: eigensolver = 'qr   '

interface operator(.cprod.)
    module procedure cart_prod_
end interface

interface diag
    module procedure get_diag_
    module procedure get_diag_cmplx_
    module procedure make_diag_
end interface

interface invert
    module procedure invert_ge_dble_
    module procedure invert_ge_cmplx_
end interface

interface diagonalize
    module procedure diagonalize_ge_dble_
    module procedure diagonalize_ge_cmplx_
end interface

interface sdiagonalize
    module procedure diagonalize_sym_dble_
    module procedure diagonalize_he_cmplx_
end interface

interface diagonalized
    module procedure diagonalized_ge_dble_
end interface

interface sdiagonalized
    module procedure diagonalized_sym_dble_
end interface

interface tostr
    module procedure tostr_int_
    module procedure tostr_dble_
end interface

! external :: ZHEEV, DGEEV, DSYEV, DGETRF, DGETRI, DGESV, ZGETRF, &
!     ZGETRI, ZGEEV, BLACS_BARRIER, BLACS_GRIDINFO, &
!     INFOG2L, DGERV2D, DGESD2D, PDGETRI, PDSYEV, &
!     PZHEEV, BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT, &
!     descinit, PDGETRF, DGSUM2D, BLACS_GRIDEXIT, pdelget, &
!     PDSYEVD, PZHEEVD, PDSYEVR, PZHEEVR, DSYEVD, ZHEEVD, &
!     DSYEVR, ZHEEVR, DGEMM, blacs_exit, ZGESD2D, ZGERV2D, &
!     ZGEMM

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     vdW(TS) ENERGY     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_ts_energy( mode, version, xyz, C6, alpha_0, R_vdw, s_R, &
        d, overlap, damping_custom, unit_cell) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        C6(size(xyz, 1)), &
        alpha_0(size(xyz, 1))
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        s_R, &
        d, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        unit_cell(3, 3)
    real(8) :: ene

    real(8) :: C6_ij, r(3), r_norm, f_damp, R_vdw_ij, overlap_ij, &
        ene_shell, ene_pair, R_cell(3)
    integer :: i_shell, i_cell, i_atom, j_atom, range_cell(3), idx_cell(3)
    real(8), parameter :: shell_thickness = 10.d0
    logical :: is_crystal, is_parallel

    is_crystal = is_in('C', mode)
    is_parallel = is_in('P', mode)

    ene = 0.d0
    i_shell = 0
    do
        i_shell = i_shell+1
        ene_shell = 0.d0
        if (is_crystal) then
            range_cell = supercell_circum(unit_cell, i_shell*shell_thickness)
        else
            range_cell = (/ 0, 0, 0 /)
        end if
        idx_cell = (/ 0, 0, -1 /)
        do i_cell = 1, product(1+2*range_cell)
            call shift_cell(idx_cell, -range_cell, range_cell)
            ! MPI code begin
            if (is_parallel .and. is_crystal) then
                if (my_task /= modulo(i_cell, n_tasks)) cycle
            end if
            ! MPI code end
            if (is_crystal) then
                R_cell = matmul(idx_cell, unit_cell)
            else
                R_cell = (/ 0.d0, 0.d0, 0.d0 /)
            end if
            do i_atom = 1, size(xyz, 1)
                ! MPI code begin
                if (is_parallel .and. .not. is_crystal) then
                    if (my_task /= modulo(i_atom, n_tasks)) cycle
                end if
                ! MPI code end
                do j_atom = 1, i_atom
                    if (i_cell == 1) then
                        if (i_atom == j_atom) cycle
                    end if
                    r = xyz(i_atom, :)-xyz(j_atom, :)-R_cell
                    r_norm = sqrt(sum(r**2))
                    if (r_norm > param_ts_cutoff_radius) cycle
                    if (r_norm >= i_shell*shell_thickness &
                        .or. r_norm < (i_shell-1)*shell_thickness) then
                        cycle
                    end if
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha_0(i_atom), alpha_0(j_atom))
                    if (present(R_vdw)) then
                        R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
                    end if
                    if (present(overlap)) then
                        overlap_ij = overlap(i_atom, j_atom)
                    end if
                    select case (version)
                        case ("fermi")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)
                        case ("fermi2")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)**2
                        case ("erf")
                            f_damp = damping_erf(r_norm, s_R*R_vdw_ij, d)
                        case ("1mexp")
                            f_damp = damping_1mexp(r_norm, s_R*R_vdw_ij, d)
                        case ("overlap")
                            f_damp = damping_overlap( &
                                r_norm, overlap_ij, C6_ij, s_R, d)
                        case ("custom")
                            f_damp = damping_custom(i_atom, j_atom)
                    end select
                    ene_pair = -C6_ij*f_damp/r_norm**6
                    if (i_atom == j_atom) then
                        ene_shell = ene_shell+ene_pair/2
                    else
                        ene_shell = ene_shell+ene_pair
                    endif
                end do ! j_atom
            end do ! i_atom
        end do ! i_cell
        ! MPI code begin
        if (is_parallel) then
            call sync_sum(ene_shell)
        end if
        ! MPI code end
        ene = ene+ene_shell
        if (.not. is_crystal) exit
        if (i_shell > 1 .and. abs(ene_shell) < param_ts_energy_accuracy) then
            call print_log("Periodic TS converged in " &
                //trim(tostr(i_shell))//" shells, " &
                //trim(tostr(i_shell*shell_thickness*bohr))//" angstroms")
            exit
        endif
    end do ! i_shell
end function get_ts_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     DIPOLE INTERACTION TENSOR (LAPACK/MPI)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine add_dipole_matrix(mode, version, xyz, alpha, R_vdw, beta, a, &
        overlap, C6, damping_custom, potential_custom, unit_cell, k_point, &
        relay, relay_c)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :)
    real(8), intent(in), optional :: &
        alpha(size(xyz, 1)), &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3), &
        unit_cell(3, 3), &
        k_point(3)
    real(8), intent(inout), optional :: relay(3*size(xyz, 1), 3*size(xyz, 1))
    complex(8), intent(inout), optional :: &
        relay_c(3*size(xyz, 1), 3*size(xyz, 1))

    real(8) :: Tpp(3, 3), R_cell(3), r(3), r_norm, R_vdw_ij, C6_ij, &
        overlap_ij, sigma_ij, volume, ewald_alpha, real_space_cutoff
    complex(8) :: Tpp_c(3, 3)
    character(len=1) :: parallel_mode
    integer :: i_atom, j_atom, i_cell, idx_cell(3), range_cell(3), i, j
    logical :: is_crystal, is_parallel, is_reciprocal, is_low_dim, mute, &
        do_ewald

    is_parallel = is_in('P', mode)
    is_reciprocal = is_in('R', mode)
    is_crystal = is_in('C', mode) .or. is_reciprocal
    is_low_dim = any(param_vacuum_axis)
    do_ewald = .false.
    mute = is_in('M', mode)
    if (is_parallel) then
        parallel_mode = 'A' ! atoms
        if (is_crystal .and. size(xyz, 1) < n_tasks) then
            parallel_mode = 'C' ! cells
        end if
    else
        parallel_mode = ''
    end if

    ! MPI code begin
    if (is_parallel) then
        ! will be restored by syncing at the end
        if (is_reciprocal) then
            relay_c = relay_c/n_tasks
        else
            relay = relay/n_tasks
        end if
    end if
    ! MPI code end
    if (is_crystal) then
        if (is_low_dim) then
            real_space_cutoff = param_dipole_low_dim_cutoff
        else if (param_ewald_on) then
            do_ewald = .true.
            volume = max(abs(dble(product(diagonalized(unit_cell)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1.d0/3)
            real_space_cutoff = 6.d0/ewald_alpha*param_ewald_real_cutoff_scaling
            call print_log('Ewald: using alpha = '//trim(tostr(ewald_alpha)) &
                //', real cutoff = '//trim(tostr(real_space_cutoff)), mute)
        else
            real_space_cutoff = param_dipole_cutoff
        end if
        range_cell = supercell_circum(unit_cell, real_space_cutoff)
    else
        range_cell = (/ 0, 0, 0 /)
    end if
    if (is_crystal) then
        call print_log('Ewald: summing real part in cell vector range of ' &
            //trim(tostr(1+2*range_cell(1)))//'x' &
            //trim(tostr(1+2*range_cell(2)))//'x' &
            //trim(tostr(1+2*range_cell(3))), mute)
    end if
    call ts(11)
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        ! MPI code begin
        if (parallel_mode == 'C') then
            if (my_task /= modulo(i_cell, n_tasks)) cycle
        end if
        ! MPI code end
        if (is_crystal) then
            R_cell = matmul(idx_cell, unit_cell)
        else
            R_cell(:) = 0.d0
        end if
        do i_atom = 1, size(xyz, 1)
            ! MPI code begin
            if (parallel_mode == 'A') then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            ! MPI code end
            !$omp parallel do private(r, r_norm, R_vdw_ij, sigma_ij, overlap_ij, C6_ij, &
            !$omp    Tpp, i, j, Tpp_c)
            do j_atom = 1, i_atom
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = xyz(i_atom, :)-xyz(j_atom, :)-R_cell
                r_norm = sqrt(sum(r**2))
                if (is_crystal .and. r_norm > real_space_cutoff) cycle
                if (present(R_vdw)) then
                    R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
                end if
                if (present(alpha)) then
                    sigma_ij = sqrt(sum(get_sigma_selfint( &
                        alpha((/ i_atom , j_atom /)))**2))
                end if
                if (present(overlap)) then
                    overlap_ij = overlap(i_atom, j_atom)
                end if
                if (present(C6)) then
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha(i_atom), alpha(j_atom))
                end if
                select case (version)
                    case ("bare")
                        Tpp = T_bare(r)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,erf")
                        Tpp = T_erf_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,fermi")
                        Tpp = T_fermi_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,overlap")
                        Tpp = T_overlap_coulomb(r, overlap_ij, C6_ij, beta, a)
                    case ("1mexp,dip")
                        Tpp = damping_1mexp(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("erf,dip")
                        Tpp = damping_erf(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi^2,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)**2*T_bare(r)
                    case ("overlap,dip")
                        Tpp = damping_overlap( &
                            r_norm, overlap_ij, C6_ij, beta, a)*T_bare(r)
                    case ("custom,dip")
                        Tpp = damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp = potential_custom(i_atom, j_atom, :, :)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("1mexp,dip,gg")
                        Tpp = (1.d0-damping_1mexp(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("erf,dip,gg")
                        Tpp = (1.d0-damping_erf(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("fermi,dip,gg")
                        Tpp = (1.d0-damping_fermi(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("custom,dip,gg")
                        Tpp = (1.d0-damping_custom(i_atom, j_atom)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                end select
                if (do_ewald) then
                    Tpp = Tpp + T_erfc(r, ewald_alpha) - T_bare(r)
                end if
                if (is_reciprocal) then
                    Tpp_c = Tpp*exp(-cmplx(0.d0, 1.d0, 8)*( &
                        dot_product(k_point, r)))
                end if
                i = 3*i_atom-2
                j = 3*j_atom-2
                if (is_reciprocal) then
                    relay_c(i:i+2, j:j+2) = relay_c(i:i+2, j:j+2) &
                        + Tpp_c
                    if (i_atom /= j_atom) then
                        relay_c(j:j+2, i:i+2) = relay_c(j:j+2, i:i+2) &
                            + conjg(transpose(Tpp_c))
                    end if
                else
                    relay(i:i+2, j:j+2) = relay(i:i+2, j:j+2) &
                        + Tpp
                    if (i_atom /= j_atom) then
                        relay(j:j+2, i:i+2) = relay(j:j+2, i:i+2) &
                            + transpose(Tpp)
                    end if
                end if
            end do ! j_atom
            !$omp end parallel do
        end do ! i_atom
    end do ! i_cell
    call ts(-11)
    ! MPI code begin
    if (is_parallel) then
        if (is_reciprocal) then
            call sync_sum(relay_c)
        else
            call sync_sum(relay)
        end if
    end if
    ! MPI code end
    if (do_ewald) then
        call add_ewald_dipole_parts( &
            mode, xyz, unit_cell, ewald_alpha, k_point, &
            relay, relay_c)
    end if
end subroutine add_dipole_matrix


subroutine add_ewald_dipole_parts(mode, xyz, unit_cell, alpha, &
        k_point, relay, relay_c)
    character(len=*), intent(in) :: mode
    real(8), intent(in) :: &
        xyz(:, :), &
        unit_cell(3, 3), &
        alpha
    real(8), intent(in), optional :: k_point(3)
    real(8), intent(inout), optional :: relay(3*size(xyz, 1), 3*size(xyz, 1))
    complex(8), intent(inout), optional :: &
        relay_c(3*size(xyz, 1), 3*size(xyz, 1))

    logical :: is_parallel, is_reciprocal, mute, do_surface
    real(8) :: rec_unit_cell(3, 3), volume, G_vector(3), r(3), k_total(3), &
        k_sq, rec_space_cutoff, Tpp(3, 3), k_prefactor(3, 3), elem
    complex(8) :: Tpp_c(3, 3)
    integer :: &
        i_atom, j_atom, i, j, i_xyz, j_xyz, idx_G_vector(3), i_G_vector, &
        range_G_vector(3)
    character(len=1) :: parallel_mode

    is_parallel = is_in('P', mode)
    is_reciprocal = is_in('R', mode)
    mute = is_in('M', mode)
    if (is_parallel) then
        parallel_mode = 'A' ! atoms
        if (size(xyz, 1) < n_tasks) then
            parallel_mode = 'G' ! G vectors
        end if
    else
        parallel_mode = ''
    end if

    ! MPI code begin
    if (is_parallel) then
        ! will be restored by syncing at the end
        if (is_reciprocal) then
            relay_c = relay_c/n_tasks
        else
            relay = relay/n_tasks
        end if
    end if
    ! MPI code end
    rec_unit_cell = 2*pi*inverted(transpose(unit_cell))
    volume = abs(dble(product(diagonalized(unit_cell))))
    rec_space_cutoff = 10.d0*alpha*param_ewald_rec_cutoff_scaling
    range_G_vector = supercell_circum(rec_unit_cell, rec_space_cutoff)
    call print_log('Ewald: using reciprocal cutoff = ' &
        //trim(tostr(rec_space_cutoff)), mute)
    call print_log('Ewald: summing reciprocal part in G vector range of ' &
        //trim(tostr(1+2*range_G_vector(1)))//'x' &
        //trim(tostr(1+2*range_G_vector(2)))//'x' &
        //trim(tostr(1+2*range_G_vector(3))), mute)
    call ts(12)
    idx_G_vector = (/ 0, 0, -1 /)
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        ! MPI code begin
        if (parallel_mode == 'G') then
            if (my_task /= modulo(i_G_vector, n_tasks)) cycle
        end if
        ! MPI code end
        G_vector = matmul(idx_G_vector, rec_unit_cell)
        if (is_reciprocal) then
            k_total = k_point+G_vector
        else
            k_total = G_vector
        end if
        k_sq = sum(k_total**2)
        if (sqrt(k_sq) > rec_space_cutoff) cycle
        k_prefactor(:, :) = 4*pi/volume*exp(-k_sq/(4*alpha**2))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz) &
                *k_total(i_xyz)*k_total(j_xyz)/k_sq
        do i_atom = 1, size(xyz, 1)
            ! MPI code begin
            if (parallel_mode == 'A') then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            ! MPI code end
            !$omp parallel do private(r, Tpp, i, j, Tpp_c)
            do j_atom = 1, i_atom
                r = xyz(i_atom, :)-xyz(j_atom, :)
                if (is_reciprocal) then
                    Tpp_c = k_prefactor*exp(cmplx(0.d0, 1.d0, 8) &
                        *dot_product(G_vector, r))
                else
                    Tpp = k_prefactor*cos(dot_product(G_vector, r))
                end if
                i = 3*i_atom-2
                j = 3*j_atom-2
                if (is_reciprocal) then
                    relay_c(i:i+2, j:j+2) = relay_c(i:i+2, j:j+2) + Tpp_c
                    if (i_atom /= j_atom) then
                        relay_c(j:j+2, i:i+2) = relay_c(j:j+2, i:i+2) &
                            + conjg(transpose(Tpp_c))
                    end if
                else
                    relay(i:i+2, j:j+2) = relay(i:i+2, j:j+2) + Tpp
                    if (i_atom /= j_atom) then
                        relay(j:j+2, i:i+2) = relay(j:j+2, i:i+2) &
                            + transpose(Tpp)
                    end if
                end if
            end do ! j_atom
            !$omp end parallel do
        end do ! i_atom
    end do ! i_G_vector
    ! MPI code begin
    if (is_parallel) then
        if (is_reciprocal) then
            call sync_sum(relay_c)
        else
            call sync_sum(relay)
        end if
    end if
    ! MPI code end
    do i_atom = 1, size(xyz, 1) ! self energy
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            if (is_reciprocal) then
                relay_c(i, i) = relay_c(i, i)-4*alpha**3/(3*sqrt(pi))
            else
                relay(i, i) = relay(i, i)-4*alpha**3/(3*sqrt(pi))
            end if
        end do
    end do
    do_surface = .true.
    if (present(k_point)) then
        k_sq = sum(k_point**2)
        if (sqrt(k_sq) > 1.d-15) then
            do_surface = .false.
            do i_atom = 1, size(xyz, 1)
            do j_atom = 1, i_atom
                do i_xyz = 1, 3
                do j_xyz = 1, 3
                    i = 3*(i_atom-1)+i_xyz
                    j = 3*(j_atom-1)+j_xyz
                    elem = 4*pi/volume*k_point(i_xyz)*k_point(j_xyz)/k_sq &
                        *exp(-k_sq/(4*alpha**2))
                    if (is_reciprocal) then
                        relay_c(i, j) = relay_c(i, j) + elem
                        if (i_atom /= j_atom) then
                            relay_c(j, i) = relay_c(j, i) + elem
                        end if
                    else
                        relay(i, j) = relay(i, j) + elem
                        if (i_atom /= j_atom) then
                            relay(j, i) = relay(j, i) + elem
                        end if
                    end if ! is_reciprocal
                end do ! j_xyz
                end do ! i_xyz
            end do ! j_atom
            end do ! i_atom
        end if ! k_sq >
    end if ! k_point present
    if (do_surface) then ! surface energy
        do i_atom = 1, size(xyz, 1)
        do j_atom = 1, i_atom
            do i_xyz = 1, 3
                i = 3*(i_atom-1)+i_xyz
                j = 3*(j_atom-1)+i_xyz
                if (is_reciprocal) then
                    relay_c(i, j) = relay_c(i, j)+4*pi/(3*volume)
                    relay_c(j, i) = relay_c(i, j)
                else
                    relay(i, j) = relay(i, j)+4*pi/(3*volume)
                    relay(j, i) = relay(i, j)
                end if
            end do ! i_xyz
        end do ! j_atom
        end do ! i_atom
    end if
    call ts(-12)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     SELF-CONSISTENT SCREENING (LAPACK/MPI)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_grid(n)
    integer, intent(in) :: n

    n_grid_omega = n
    if (allocated(omega_grid)) deallocate(omega_grid)
    if (allocated(omega_grid_w)) deallocate(omega_grid_w)
    allocate (omega_grid(0:n))
    allocate (omega_grid_w(0:n))
    omega_grid(0) = 0.d0
    omega_grid_w(0) = 0.d0
    call get_omega_grid(n, 0.6d0, omega_grid(1:n), omega_grid_w(1:n))
    call print_log( &
        "Initialized a radial integration grid of "//trim(tostr(n))//" points." &
    )
    call print_log( &
        "Relative quadrature error in C6 of carbon atom: "// &
        trim(tostr(test_frequency_grid())) &
    )
end subroutine


real(8) function test_frequency_grid() result(error)
    real(8) :: alpha(0:n_grid_omega, 1)

    alpha = alpha_dynamic_ts_all('C', n_grid_omega, (/ 21.d0 /), C6=(/ 99.5d0 /))
    error = abs(get_total_C6_from_alpha(alpha)/99.5d0-1.d0)
end function


subroutine get_omega_grid(n, L, x, w)
    integer, intent(in) :: n
    real(8), intent(in) :: L
    real(8), intent(out) :: x(n), w(n)

    call gauss_legendre(n, x, w)
    w = 2*L/(1-x)**2*w
    x = L*(1+x)/(1-x)
    w = w(n:1:-1)
    x = x(n:1:-1)
end subroutine get_omega_grid


subroutine gauss_legendre(n, r, w)
    use mbd_interface, only: legendre_precision

    integer, intent(in) :: n
    real(8), intent(out) :: r(n), w(n)

    integer, parameter :: q = legendre_precision
    integer, parameter :: n_iter = 1000
    real(q) :: x, f, df, dx
    integer :: k, iter, i
    real(q) :: Pk(0:n), Pk1(0:n-1), Pk2(0:n-2)

    if (n == 1) then
        r(1) = 0.d0
        w(1) = 2.d0
        return
    end if
    Pk2(0) = 1._q  ! k = 0
    Pk1(0:1) = (/ 0._q, 1._q /)  ! k = 1
    do k = 2, n
        Pk(0:k) = ((2*k-1)*(/ 0.0_q, Pk1(0:k-1) /)-(k-1)*(/ Pk2(0:k-2), 0._q, 0._q /))/k
        if (k < n) then
            Pk2(0:k-1) = Pk1(0:k-1)
            Pk1(0:k) = Pk(0:k)
        end if
    end do
    ! now Pk contains k-th Legendre polynomial
    do i = 1, n
        x = cos(pi*(i-0.25_q)/(n+0.5_q))
        do iter = 1, n_iter
            df = 0._q
            f = Pk(n)
            do k = n-1, 0, -1
                df = f + x*df
                f = Pk(k) + x*f
            end do
            dx = f/df
            x = x-dx
            if (abs(dx) < 10*epsilon(dx)) exit
        end do
        r(i) = dble(x)
        w(i) = dble(2/((1-x**2)*df**2))
    end do
end subroutine


subroutine destroy_grid()
    deallocate (omega_grid)
    deallocate (omega_grid_w)
end subroutine


subroutine init_eqi_grid(n, a, b)
    integer, intent(in) :: n
    real(8), intent(in) :: a, b

    real(8) :: delta
    integer :: i

    n_grid_omega = n
    if (allocated(omega_grid)) deallocate(omega_grid)
    if (allocated(omega_grid_w)) deallocate(omega_grid_w)
    allocate (omega_grid(0:n))
    allocate (omega_grid_w(0:n))
    omega_grid(0) = 0.d0
    omega_grid_w(0) = 0.d0
    delta = (b-a)/n
    omega_grid(1:n) = (/ (a+delta/2+i*delta, i = 0, n-1) /)
    omega_grid_w(1:n) = delta
end subroutine


function run_scs(mode, version, xyz, alpha, R_vdw, beta, a, unit_cell) &
        result(alpha_scs)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha(:, :)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        unit_cell(3, 3)
    real(8) :: alpha_scs(size(alpha, 1), size(alpha, 2))

    real(8) :: alpha_full(3*size(xyz, 1), 3*size(xyz, 1))
    integer :: i_grid_omega
    logical :: is_parallel
    character(len=1) :: mute

    is_parallel = is_in('P', mode)
    if (is_in('M', mode)) then
        mute = 'M'
    else
        mute = ''
    end if

    alpha_scs(:, :) = 0.d0
    do i_grid_omega = 0, n_grid_omega
        ! MPI code begin
        if (is_parallel) then
            if (my_task /= modulo(i_grid_omega, n_tasks)) cycle
        end if
        ! MPI code end
        alpha_full = do_scs( &
            blanked('P', mode)//mute, &
            version, &
            xyz, &
            alpha=alpha(i_grid_omega+1, :), &
            R_vdw=R_vdw, &
            beta=beta, &
            a=a, &
            unit_cell=unit_cell)
        alpha_scs(i_grid_omega+1, :) = contract_polarizability(alpha_full)
        mute = 'M'
    end do
    ! MPI code begin
    if (is_parallel) then
        call sync_sum(alpha_scs)
    end if
    ! MPI code end
end function run_scs


function do_scs(mode, version, xyz, alpha, R_vdw, beta, a, lam, unit_cell) &
        result(alpha_full)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha(:)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        unit_cell(3, 3), &
        lam
    real(8) :: alpha_full(3*size(xyz, 1), 3*size(xyz, 1))

    logical :: scale_lambda
    integer :: i_atom, i_xyz, i

    scale_lambda = is_in('L', mode)

    alpha_full(:, :) = 0.d0
    call add_dipole_matrix( &
        blanked('R', mode), &
        version, &
        xyz=xyz, &
        alpha=alpha, &
        R_vdw=R_vdw, &
        beta=beta, &
        a=a, &
        unit_cell=unit_cell, &
        relay=alpha_full)
    if (scale_lambda) then
        alpha_full = lam*alpha_full
    end if
    do i_atom = 1, size(xyz, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            alpha_full(i, i) = alpha_full(i, i)+1.d0/alpha(i_atom)
        end do
    end do
    call ts(31)
    call invert(alpha_full)
    call ts(-31)
end function do_scs


function do_scs_k_point(mode, version, xyz, alpha, k_point, R_vdw, &
        beta, a, lam, unit_cell) result(alpha_full)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha(:), &
        k_point(3)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        unit_cell(3, 3), &
        lam
    complex(8) :: alpha_full(3*size(xyz, 1), 3*size(xyz, 1))

    logical :: scale_lambda
    integer :: i_atom, i_xyz, i

    scale_lambda = is_in('L', mode)

    alpha_full(:, :) = cmplx(0.d0, 0.d0, 8)
    call add_dipole_matrix( &
        mode//'R', &
        version, &
        xyz=xyz, &
        alpha=alpha, &
        R_vdw=R_vdw, &
        beta=beta, &
        a=a, &
        unit_cell=unit_cell, &
        k_point=k_point, &
        relay_c=alpha_full)
    if (scale_lambda) then
        alpha_full = lam*alpha_full
    end if
    do i_atom = 1, size(xyz, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            alpha_full(i, i) = alpha_full(i, i)+1.d0/alpha(i_atom)
        end do
    end do
    call ts(32)
    call invert(alpha_full)
    call ts(-32)
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     MBD ENERGY (LAPACK/MPI)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_mbd_energy(mode, version, xyz, alpha_0, omega, &
        supercell, k_grid, unit_cell, R_vdw, beta, a, overlap, C6, damping_custom, &
        potential_custom) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha_0(size(xyz, 1)), &
        omega(size(xyz, 1))
    integer, intent(in), optional :: supercell(3)
    real(8), intent(in), optional :: &
        k_grid(:, :), &
        unit_cell(3, 3)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3)
    real(8) :: ene

    logical :: is_parallel, do_rpa, is_reciprocal, is_crystal
    real(8), allocatable :: alpha(:, :)

    is_parallel = is_in('P', mode)
    is_crystal = is_in('C', mode)
    do_rpa = is_in('Q', mode)
    is_reciprocal = is_in('R', mode)
    if (.not. is_crystal) then
        if (.not. do_rpa) then
            ene = get_single_mbd_energy(mode, version, xyz, alpha_0, omega, R_vdw, &
                beta, a, overlap, C6, damping_custom, potential_custom, unit_cell)
        else
            allocate (alpha(0:n_grid_omega, size(alpha_0)))
            alpha = alpha_dynamic_ts_all(mode, n_grid_omega, alpha_0, C6, omega)
            ene = get_single_rpa_energy(mode, version, xyz, alpha, R_vdw, beta, &
                a, overlap, C6, damping_custom, potential_custom, unit_cell)
            deallocate (alpha)
        end if
    else
        if (is_reciprocal) then
            ene = get_reciprocal_mbd_energy(mode, version, xyz, alpha_0, omega, &
                k_grid, unit_cell, R_vdw, beta, a, overlap, C6, damping_custom, &
                potential_custom)
        else
            ene = get_supercell_mbd_energy(mode, version, xyz, alpha_0, omega, &
                unit_cell, supercell, R_vdw, beta, a, C6)
        end if
    end if
end function get_mbd_energy


function get_supercell_mbd_energy(mode, version, xyz, alpha_0, omega, &
        unit_cell, supercell, R_vdw, beta, a, C6, rpa_orders) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha_0(size(xyz, 1)), &
        omega(size(xyz, 1)), &
        unit_cell(3, 3)
    integer, intent(in) :: supercell(3)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        C6(size(xyz, 1))
    real(8), intent(out), optional :: rpa_orders(20)
    real(8) :: ene

    logical :: do_rpa, get_orders
    real(8) :: R_cell(3)
    integer :: i_atom, i
    integer :: i_cell
    integer :: idx_cell(3), n_cells

    real(8), allocatable :: &
        xyz_super(:, :), alpha_0_super(:), omega_super(:), &
        R_vdw_super(:), C6_super(:), alpha_ts_super(:, :)
    real(8) :: unit_cell_super(3, 3)

    do_rpa = is_in('Q', mode)
    get_orders = is_in('O', mode)

    n_cells = product(supercell)
    do i = 1, 3
        unit_cell_super(i, :) = unit_cell(i, :)*supercell(i)
    end do
    allocate (xyz_super(n_cells*size(xyz, 1), 3))
    allocate (alpha_0_super(n_cells*size(alpha_0)))
    allocate (alpha_ts_super(0:n_grid_omega, n_cells*size(alpha_0)))
    allocate (omega_super(n_cells*size(omega)))
    allocate (R_vdw_super(n_cells*size(R_vdw)))
    allocate (C6_super(n_cells*size(C6)))
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, n_cells
        call shift_cell(idx_cell, (/ 0, 0, 0 /), supercell-1)
        R_cell = matmul(idx_cell, unit_cell)
        do i_atom = 1, size(xyz, 1)
            i = (i_cell-1)*size(xyz, 1)+i_atom
            xyz_super(i, :) = xyz(i_atom, :)+R_cell
            alpha_0_super(i) = alpha_0(i_atom)
            omega_super(i) = omega(i_atom)
            if (present(R_vdw)) then
                R_vdw_super(i) = R_vdw(i_atom)
            end if
            if (present(C6)) then
                C6_super(i) = C6(i_atom)
            end if
        end do
    end do
    if (do_rpa) then
        alpha_ts_super = alpha_dynamic_ts_all( &
            'O', n_grid_omega, alpha_0_super, omega=omega_super)
        ene = get_single_rpa_energy( &
            mode, &
            version, &
            xyz_super, &
            alpha_ts_super, &
            R_vdw=R_vdw_super, &
            beta=beta, &
            a=a, &
            C6=C6_super, &
            unit_cell=unit_cell_super, &
            rpa_orders=rpa_orders)
    else
        ene = get_single_mbd_energy( &
            mode, &
            version, &
            xyz_super, &
            alpha_0_super, &
            omega_super, &
            R_vdw=R_vdw_super, &
            beta=beta, &
            a=a, &
            C6=C6_super, &
            unit_cell=unit_cell_super)
    end if
    deallocate (xyz_super)
    deallocate (alpha_0_super)
    deallocate (alpha_ts_super)
    deallocate (omega_super)
    deallocate (R_vdw_super)
    deallocate (C6_super)
    ene = ene/n_cells
    if (get_orders) then
        rpa_orders = rpa_orders/n_cells
    end if
end function get_supercell_mbd_energy


function get_single_mbd_energy(mode, version, xyz, alpha_0, omega, R_vdw, &
        beta, a, overlap, C6, damping_custom, potential_custom, unit_cell, &
        mode_enes, modes) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha_0(size(xyz, 1)), &
        omega(size(xyz, 1))
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3), &
        unit_cell(3, 3)
    real(8), intent(out), optional :: &
        mode_enes(3*size(xyz, 1)), &
        modes(3*size(xyz, 1), 3*size(xyz, 1))
    real(8) :: ene

    real(8) :: relay(3*size(xyz, 1), 3*size(xyz, 1))
    real(8) :: eigs(3*size(xyz, 1))
    integer :: i_atom, j_atom, i_xyz, i, j
    integer :: n_negative_eigs
    logical :: get_eigenvalues, get_eigenvectors, is_parallel

    get_eigenvalues = is_in('E', mode)
    get_eigenvectors = is_in('V', mode)
    is_parallel = is_in('P', mode)

    relay(:, :) = 0.d0
    call add_dipole_matrix( & ! relay = T
        mode, &
        version, &
        xyz, &
        alpha=alpha_0, &
        R_vdw=R_vdw, &
        beta=beta, &
        a=a, &
        overlap=overlap, &
        C6=C6, &
        damping_custom=damping_custom, &
        potential_custom=potential_custom, &
        unit_cell=unit_cell, &
        relay=relay)
    do i_atom = 1, size(xyz, 1)
        do j_atom = 1, size(xyz, 1)
            i = 3*i_atom-2
            j = 3*j_atom-2
            relay(i:i+2, j:j+2) = & ! relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                *sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                relay(i:i+2, j:j+2)
        end do
    end do
    do i_atom = 1, size(xyz, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            relay(i, i) = relay(i, i)+omega(i_atom)**2
            ! relay = w^2+sqrt(a*a)*w*w*T
        end do
    end do
    call ts(21)
    if (.not. is_parallel .or. my_task == 0) then
        if (get_eigenvectors) then
            call sdiagonalize('V', relay, eigs)
            modes = relay
        else
            call sdiagonalize('N', relay, eigs)
        end if
    end if
    ! MPI code begin
    if (is_parallel) then
        if (get_eigenvectors) call broadcast(modes)
        call broadcast(eigs)
    end if
    ! MPI code end
    call ts(-21)
    if (get_eigenvalues) then
        mode_enes(:) = 0.d0
        where (eigs > 0) mode_enes = sqrt(eigs)
    end if
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call print_warning( &
            "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues" &
        )
        if (param_zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    end if
    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)
end function get_single_mbd_energy


function get_reciprocal_mbd_energy(mode, version, xyz, alpha_0, omega, &
        k_grid, unit_cell, R_vdw, beta, a, overlap, C6, damping_custom, &
        potential_custom, mode_enes, modes, rpa_orders) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha_0(size(xyz, 1)), &
        omega(size(xyz, 1)), &
        k_grid(:, :), &
        unit_cell(3, 3)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3)
    real(8), intent(out), optional :: &
        mode_enes(size(k_grid, 1), 3*size(xyz, 1)), &
        rpa_orders(size(k_grid, 1), 20)
    complex(8), intent(out), optional :: &
        modes(size(k_grid, 1), 3*size(xyz, 1), 3*size(xyz, 1))
    real(8) :: ene

    logical :: &
        is_parallel, do_rpa, get_orders, get_eigenvalues, get_eigenvectors
    integer :: i_kpt
    real(8) :: k_point(3), alpha_ts(0:n_grid_omega, size(xyz, 1))
    character(len=1) :: mute

    is_parallel = is_in('P', mode)
    do_rpa = is_in('Q', mode)
    get_eigenvalues= is_in('E', mode)
    get_eigenvectors= is_in('V', mode)
    get_orders = is_in('O', mode)
    if (is_in('M', mode)) then
        mute = 'M'
    else
        mute = ''
    end if
    if (do_rpa) then
        alpha_ts = alpha_dynamic_ts_all('O', n_grid_omega, alpha_0, omega=omega)
    endif
    ene = 0.d0
    if (get_eigenvalues) mode_enes(:, :) = 0.d0
    if (get_eigenvectors) modes(:, :, :) = 0.d0
    if (get_orders) rpa_orders(:, :) = 0.d0
    do i_kpt = 1, size(k_grid, 1)
        ! MPI code begin
        if (is_parallel) then
            if (my_task /= modulo(i_kpt, n_tasks)) cycle
        end if
        ! MPI code end
        k_point = k_grid(i_kpt, :)
        if (do_rpa) then
            if (get_orders) then
                ene = ene+get_single_reciprocal_rpa_ene( &
                    blanked('P', mode)//mute, &
                    version, &
                    xyz, &
                    alpha_ts, &
                    k_point, &
                    unit_cell, &
                    R_vdw=R_vdw, &
                    beta=beta, &
                    a=a, &
                    overlap=overlap, &
                    C6=C6, &
                    damping_custom=damping_custom, &
                    potential_custom=potential_custom, &
                    rpa_orders=rpa_orders(i_kpt, :))
            else
                ene = ene+get_single_reciprocal_rpa_ene( &
                    blanked('P', mode)//mute, &
                    version, &
                    xyz, &
                    alpha_ts, &
                    k_point, &
                    unit_cell, &
                    R_vdw=R_vdw, &
                    beta=beta, &
                    a=a, &
                    overlap=overlap, &
                    C6=C6, &
                    damping_custom=damping_custom, &
                    potential_custom=potential_custom)
            end if
        else
            if (get_eigenvalues .and. get_eigenvectors) then
                ene = ene+get_single_reciprocal_mbd_ene( &
                    blanked('P', mode)//mute, &
                    version, &
                    xyz, &
                    alpha_0, &
                    omega, &
                    k_point, &
                    unit_cell, &
                    R_vdw=R_vdw, &
                    beta=beta, &
                    a=a, &
                    overlap=overlap, &
                    C6=C6, &
                    damping_custom=damping_custom, &
                    potential_custom=potential_custom, &
                    mode_enes=mode_enes(i_kpt, :), &
                    modes=modes(i_kpt, :, :))
            elseif (get_eigenvalues) then
                ene = ene+get_single_reciprocal_mbd_ene( &
                    blanked('P', mode)//mute, &
                    version, &
                    xyz, &
                    alpha_0, &
                    omega, &
                    k_point, &
                    unit_cell, &
                    R_vdw=R_vdw, &
                    beta=beta, &
                    a=a, &
                    overlap=overlap, &
                    C6=C6, &
                    damping_custom=damping_custom, &
                    potential_custom=potential_custom, &
                    mode_enes=mode_enes(i_kpt, :))
            elseif (get_eigenvectors) then
                ene = ene+get_single_reciprocal_mbd_ene( &
                    blanked('P', mode)//mute, &
                    version, &
                    xyz, &
                    alpha_0, &
                    omega, &
                    k_point, &
                    unit_cell, &
                    R_vdw=R_vdw, &
                    beta=beta, &
                    a=a, &
                    overlap=overlap, &
                    C6=C6, &
                    damping_custom=damping_custom, &
                    potential_custom=potential_custom, &
                    modes=modes(i_kpt, :, :))
            else
                ene = ene+get_single_reciprocal_mbd_ene( &
                    blanked('P', mode)//mute, &
                    version, &
                    xyz, &
                    alpha_0, &
                    omega, &
                    k_point, &
                    unit_cell, &
                    R_vdw=R_vdw, &
                    beta=beta, &
                    a=a, &
                    overlap=overlap, &
                    C6=C6, &
                    damping_custom=damping_custom, &
                    potential_custom=potential_custom)
            end if
        end if
        mute = 'M'
    end do ! k_point loop
    ! MPI code begin
    if (is_parallel) then
        call sync_sum(ene)
        if (get_eigenvalues) call sync_sum(mode_enes)
        if (get_eigenvectors) call sync_sum(modes)
        if (get_orders) call sync_sum(rpa_orders)
    end if
    ! MPI code end
    ene = ene/size(k_grid, 1)
    if (get_orders) rpa_orders = rpa_orders/size(k_grid, 1)
end function get_reciprocal_mbd_energy


function get_single_reciprocal_mbd_ene(mode, version, xyz, alpha_0, omega, &
        k_point, unit_cell, R_vdw, beta, a, overlap, C6, damping_custom, &
        potential_custom, mode_enes, modes) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha_0(size(xyz, 1)), &
        omega(size(xyz, 1)), &
        k_point(3), &
        unit_cell(3, 3)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3)
    real(8), intent(out), optional :: mode_enes(3*size(xyz, 1))
    complex(8), intent(out), optional :: modes(3*size(xyz, 1), 3*size(xyz, 1))
    real(8) :: ene

    complex(8) :: relay(3*size(xyz, 1), 3*size(xyz, 1))
    real(8) :: eigs(3*size(xyz, 1))
    integer :: i_atom, j_atom, i_xyz, i, j
    integer :: n_negative_eigs
    logical :: get_eigenvalues, get_eigenvectors, is_parallel

    get_eigenvalues = is_in('E', mode)
    get_eigenvectors = is_in('V', mode)
    is_parallel = is_in('P', mode)

    relay(:, :) = cmplx(0.d0, 0.d0, 8)
    call add_dipole_matrix( & ! relay = T
        mode, &
        version, &
        xyz, &
        alpha=alpha_0, &
        R_vdw=R_vdw, &
        beta=beta, &
        a=a, &
        overlap=overlap, &
        C6=C6, &
        damping_custom=damping_custom, &
        potential_custom=potential_custom, &
        unit_cell=unit_cell, &
        k_point=k_point, &
        relay_c=relay)
    do i_atom = 1, size(xyz, 1)
        do j_atom = 1, size(xyz, 1)
            i = 3*i_atom-2
            j = 3*j_atom-2
            relay(i:i+2, j:j+2) = & ! relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                *sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                relay(i:i+2, j:j+2)
        end do
    end do
    do i_atom = 1, size(xyz, 1)
        do i_xyz = 1, 3
            i = 3*(i_atom-1)+i_xyz
            relay(i, i) = relay(i, i)+omega(i_atom)**2
            ! relay = w^2+sqrt(a*a)*w*w*T
        end do
    end do
    call ts(22)
    if (.not. is_parallel .or. my_task == 0) then
        if (get_eigenvectors) then
            call sdiagonalize('V', relay, eigs)
            modes = relay
        else
            call sdiagonalize('N', relay, eigs)
        end if
    end if
    ! MPI code begin
    if (is_parallel) then
        if (get_eigenvectors) call broadcast(modes)
        call broadcast(eigs)
    end if
    ! MPI code end
    call ts(-22)
    if (get_eigenvalues) then
        mode_enes(:) = 0.d0
        where (eigs > 0) mode_enes = sqrt(eigs)
    end if
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call print_warning( &
            "CDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues" &
        )
        if (param_zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    end if
    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)
end function get_single_reciprocal_mbd_ene


function get_single_rpa_energy(mode, version, xyz, alpha, R_vdw, beta, &
        a, overlap, C6, damping_custom, potential_custom, unit_cell, &
        rpa_orders) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha(:, :)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3), &
        unit_cell(3, 3)
    real(8), intent(out), optional :: rpa_orders(20)
    real(8) :: ene

    real(8), dimension(3*size(xyz, 1), 3*size(xyz, 1)) :: relay, AT
    complex(8) :: eigs(3*size(xyz, 1))
    integer :: i_atom, i_xyz, i_grid_omega, i
    integer :: n_order, n_negative_eigs
    logical :: is_parallel, get_orders
    character(len=1) :: mute

    is_parallel = is_in('P', mode)
    get_orders = is_in('O', mode)
    if (is_in('M', mode)) then
        mute = 'M'
    else
        mute = ''
    end if

    ene = 0.d0
    do i_grid_omega = 0, n_grid_omega
        ! MPI code begin
        if (is_parallel) then
            if (my_task /= modulo(i_grid_omega, n_tasks)) cycle
        end if
        ! MPI code end
        relay(:, :) = 0.d0
        call add_dipole_matrix( & ! relay = T
            blanked('P', mode)//mute, &
            version, &
            xyz, &
            alpha=alpha(i_grid_omega+1, :), &
            R_vdw=R_vdw, &
            beta=beta, &
            a=a, &
            overlap=overlap, &
            C6=C6, &
            damping_custom=damping_custom, &
            potential_custom=potential_custom, &
            unit_cell=unit_cell, &
            relay=relay)
        do i_atom = 1, size(xyz, 1)
            do i_xyz = 1, 3
                i = (i_atom-1)*3+i_xyz
                relay(i, :) = alpha(i_grid_omega+1, i_atom)*relay(i, :)
                ! relay = alpha*T
            end do
        end do
        AT = relay
        do i = 1, 3*size(xyz, 1)
            relay(i, i) = 1.d0+relay(i, i) ! relay = 1+alpha*T
        end do
        call ts(23)
        call diagonalize('N', relay, eigs)
        call ts(-23)
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            call print_warning("1+AT matrix has " &
                //trim(tostr(n_negative_eigs))//" negative eigenvalues")
        end if
        ene = ene+1.d0/(2*pi)*sum(log(dble(eigs)))*omega_grid_w(i_grid_omega)
        if (get_orders) then
            call ts(24)
            call diagonalize('N', AT, eigs)
            call ts(-24)
            do n_order = 2, param_rpa_order_max
                rpa_orders(n_order) = rpa_orders(n_order) &
                    +(-1.d0/(2*pi)*(-1)**n_order &
                    *sum(dble(eigs)**n_order)/n_order) &
                    *omega_grid_w(i_grid_omega)
            end do
        end if
        mute = 'M'
    end do
    if (is_parallel) then
        call sync_sum(ene)
        if (get_orders) then
            call sync_sum(rpa_orders)
        end if
    end if
end function get_single_rpa_energy


function get_single_reciprocal_rpa_ene(mode, version, xyz, alpha, k_point, &
        unit_cell, R_vdw, beta, a, overlap, C6, damping_custom, &
        potential_custom, rpa_orders) result(ene)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha(:, :), &
        k_point(3), &
        unit_cell(3, 3)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        overlap(size(xyz, 1), size(xyz, 1)), &
        C6(size(xyz, 1)), &
        damping_custom(size(xyz, 1), size(xyz, 1)), &
        potential_custom(size(xyz, 1), size(xyz, 1), 3, 3)
    real(8), intent(out), optional :: rpa_orders(20)
    real(8) :: ene

    complex(8), dimension(3*size(xyz, 1), 3*size(xyz, 1)) :: relay, AT
    complex(8) :: eigs(3*size(xyz, 1))
    integer :: i_atom, i_xyz, i_grid_omega, i
    integer :: n_order, n_negative_eigs
    logical :: is_parallel, get_orders
    character(len=1) :: mute

    is_parallel = is_in('P', mode)
    get_orders = is_in('O', mode)
    if (is_in('M', mode)) then
        mute = 'M'
    else
        mute = ''
    end if

    ene = 0.d0
    if (get_orders) rpa_orders(:) = 0.d0
    do i_grid_omega = 0, n_grid_omega
        ! MPI code begin
        if (is_parallel) then
            if (my_task /= modulo(i_grid_omega, n_tasks)) cycle
        end if
        ! MPI code end
        relay(:, :) = cmplx(0.d0, 0.d0, 8)
        call add_dipole_matrix( & ! relay = T
            blanked('P', mode)//mute, &
            version, &
            xyz, &
            alpha=alpha(i_grid_omega+1, :), &
            R_vdw=R_vdw, &
            beta=beta, &
            a=a, &
            overlap=overlap, &
            C6=C6, &
            damping_custom=damping_custom, &
            potential_custom=potential_custom, &
            k_point=k_point, &
            unit_cell=unit_cell, &
            relay_c=relay)
        do i_atom = 1, size(xyz, 1)
            do i_xyz = 1, 3
                i = (i_atom-1)*3+i_xyz
                relay(i, :) = alpha(i_grid_omega+1, i_atom)*relay(i, :)
                ! relay = alpha*T
            end do
        end do
        AT = relay
        do i = 1, 3*size(xyz, 1)
            relay(i, i) = 1.d0+relay(i, i) ! relay = 1+alpha*T
        end do
        call ts(25)
        call diagonalize('N', relay, eigs)
        call ts(-25)
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            call print_warning("1+AT matrix has " &
                //trim(tostr(n_negative_eigs))//" negative eigenvalues")
        end if
        ene = ene+1.d0/(2*pi)*dble(sum(log(eigs)))*omega_grid_w(i_grid_omega)
        if (get_orders) then
            call ts(26)
            call diagonalize('N', AT, eigs)
            call ts(-26)
            do n_order = 2, param_rpa_order_max
                rpa_orders(n_order) = rpa_orders(n_order) &
                    +(-1.d0)/(2*pi)*(-1)**n_order &
                    *dble(sum(eigs**n_order))/n_order &
                    *omega_grid_w(i_grid_omega)
            end do
        end if
        mute = 'M'
    end do
    if (is_parallel) then
        call sync_sum(ene)
        if (get_orders) then
            call sync_sum(rpa_orders)
        end if
    end if
end function get_single_reciprocal_rpa_ene


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     ADDITIONAL ANALYSIS TOOLS     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function mbd_nbody( &
!         xyz, &
!         alpha_0, &
!         omega, &
!         version, &
!         R_vdw, beta, a, &
!         my_task, n_tasks) &
!         result(ene_orders)
!     real(8), intent(in) :: &
!         xyz(:, :), &
!         alpha_0(size(xyz, 1)), &
!         omega(size(xyz, 1)), &
!         R_vdw(size(xyz, 1)), &
!         beta, a
!     character(len=*), intent(in) :: version
!     integer, intent(in), optional :: my_task, n_tasks
!     real(8) :: ene_orders(20)
!
!     integer :: &
!         multi_index(param_mbd_nbody_max), i_body, j_body, i_tuple, &
!         i_atom_ind, j_atom_ind, i_index
!     real(8) :: ene
!     logical :: is_parallel
!     
!     is_parallel = .false.
!     if (present(n_tasks)) then
!         if (n_tasks > 0) then
!             is_parallel = .true.
!         end if
!     end if
!     ene_orders(:) = 0.d0
!     do i_body = 2, param_mbd_nbody_max
!         i_tuple = 0
!         multi_index(1:i_body-1) = 1
!         multi_index(i_body:param_mbd_nbody_max) = 0
!         do
!             multi_index(i_body) = multi_index(i_body)+1
!             do i_index = i_body, 2, -1
!                 if (multi_index(i_index) > size(xyz, 1)) then
!                     multi_index(i_index) = 1
!                     multi_index(i_index-1) = multi_index(i_index-1)+1
!                 end if
!             end do
!             if (multi_index(1) > size(xyz, 1)) exit
!             if (any(multi_index(1:i_body-1)-multi_index(2:i_body) >= 0)) cycle
!             i_tuple = i_tuple+1
!             if (is_parallel) then
!                 if (my_task /= modulo(i_tuple, n_tasks)) cycle
!             end if
!             ene = get_mbd_energy( &
!                 xyz(multi_index(1:i_body), :), &
!                 alpha_0(multi_index(1:i_body)), &
!                 omega(multi_index(1:i_body)), &
!                 version, &
!                 R_vdw(multi_index(1:i_body)), &
!                 beta, a)
!             ene_orders(i_body) = ene_orders(i_body) &
!                 +ene+3.d0/2*sum(omega(multi_index(1:i_body)))
!         end do ! i_tuple
!     end do ! i_body
!     if (is_parallel) then
!         call sync_sum(ene_orders, size(ene_orders))
!     end if
!     ene_orders(1) = 3.d0/2*sum(omega)
!     do i_body = 2, min(param_mbd_nbody_max, size(xyz, 1))
!         do j_body = 1, i_body-1
!             ene_orders(i_body) = ene_orders(i_body) &
!                 -nbody_coeffs(j_body, i_body, size(xyz, 1))*ene_orders(j_body)
!         end do
!     end do
!     ene_orders(1) = sum(ene_orders(2:param_mbd_nbody_max))
! end function mbd_nbody


function eval_mbd_nonint_density(pts, xyz, charges, masses, omegas) result(rho)
    real(8), intent(in) :: &
        pts(:, :), &
        xyz(:, :), &
        charges(size(xyz, 1)), &
        masses(size(xyz, 1)), &
        omegas(size(xyz, 1))
    real(8) :: rho(size(pts, 1))

    integer :: i_pt, i_atom, n_atoms
    real(8), dimension(size(xyz, 1)) :: pre, kernel, rsq

    pre = charges*(masses*omegas/pi)**(3.d0/2)
    kernel = masses*omegas
    n_atoms = size(xyz, 1)
    rho(:) = 0.d0
    do i_pt = 1, size(pts, 1)
        if (my_task /= modulo(i_pt, n_tasks)) cycle
        forall (i_atom = 1:n_atoms)
            rsq(i_atom) = sum((pts(i_pt, :)-xyz(i_atom, :))**2)
        end forall
        rho(i_pt) = sum(pre*exp(-kernel*rsq))
    end do
    call sync_sum(rho)
end function


function eval_mbd_int_density(pts, xyz, charges, masses, omegas, modes) result(rho)
    real(8), intent(in) :: &
        pts(:, :), &
        xyz(:, :), &
        charges(size(xyz, 1)), &
        masses(size(xyz, 1)), &
        omegas(3*size(xyz, 1)), &
        modes(3*size(xyz, 1), 3*size(xyz, 1))
    real(8) :: rho(size(pts, 1))

    integer :: i_pt, i_atom, n_atoms, i, i_xyz, j_xyz
    integer :: self(3), other(3*(size(xyz, 1)-1))
    real(8) :: &
        pre(size(xyz, 1)), &
        factor(size(xyz, 1)), &
        rdiffsq(3, 3), &
        omegas_p(3*size(xyz, 1), 3*size(xyz, 1)), &
        kernel(3, 3, size(xyz, 1)), &
        rdiff(3)

    omegas_p = matmul(matmul(modes, diag(omegas)), transpose(modes))
    n_atoms = size(xyz, 1)
    kernel(:, :, :) = 0.d0
    pre(:) = 0.d0
    do i_atom = 1, n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        self(:) = (/ (3*(i_atom-1)+i, i = 1, 3) /)
        other(:) = (/ (i, i = 1, 3*(i_atom-1)),  (i, i = 3*i_atom+1, 3*n_atoms) /)
        kernel(:, :, i_atom) = masses(i_atom) &
            *(omegas_p(self, self) &
                -matmul(matmul(omegas_p(self, other), inverted(omegas_p(other, other))), &
                    omegas_p(other, self)))
        pre(i_atom) = charges(i_atom)*(masses(i_atom)/pi)**(3.d0/2) &
            *sqrt(product(omegas)/product(sdiagonalized(omegas_p(other, other))))
    end do
    call sync_sum(kernel)
    call sync_sum(pre)
    rho(:) = 0.d0
    do i_pt = 1, size(pts, 1)
        if (my_task /= modulo(i_pt, n_tasks)) cycle
        do i_atom = 1, n_atoms
            rdiff(:) = pts(i_pt, :)-xyz(i_atom, :)
            forall (i_xyz = 1:3, j_xyz = 1:3)
                rdiffsq(i_xyz, j_xyz) = rdiff(i_xyz)*rdiff(j_xyz)
            end forall
            factor(i_atom) = sum(kernel(:, :, i_atom)*rdiffsq(:, :))
        end do
        rho(i_pt) = sum(pre*exp(-factor))
    end do
    call sync_sum(rho)
end function


function eval_mbd_int_density_io(gridpts, xyz, charges, masses, &
                                 omega_int, fname_modes) &
                                 result(rho)
    real(8),          intent(in) :: gridpts(:, :), &
                                     xyz(:, :), &
                                     charges(size(xyz, 1)), &
                                     masses(size(xyz, 1)), &
                                     omega_int(3*size(xyz, 1))
    character(len=*), intent(in) :: fname_modes

    real(8) :: rho(size(gridpts, 1))

    integer :: i_pt, i_atom, n_atoms, i, i_xyz, j_xyz
    integer :: self(3), other(3*(size(xyz, 1)-1))
    real(8) :: pre(size(xyz, 1)), &
                factor(size(xyz, 1)), &
                rdiffsq(3,3), &
                kernel(3, 3, size(xyz, 1)), &
                rdiff(3), om_3_3(3,3), &
                evals_oo(3*(size(xyz, 1)-1)), &
                om_ss(3,3)

    real(8), allocatable :: modes(:,:), mm_helper(:,:), omegas_p(:,:), &
                             om_oo(:,:), om_oo_inv(:,:), om_helper(:,:),&
                             om_so(:,:)


    n_atoms = size(xyz, 1)
    allocate( modes(3*n_atoms, 3*n_atoms) )
    call read_mbd_modes(modes, fname_modes)
    allocate( mm_helper(3*n_atoms, 3*n_atoms) )
    mm_helper(:,:) = 0.d0
    do i_atom=1, 3*n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        mm_helper(:,i_atom) = sqrt(omega_int(i_atom))*modes(:,i_atom)
    enddo
    deallocate( modes )
    call sync_sum(mm_helper)
    allocate( omegas_p(3*n_atoms, 3*n_atoms) )
    if ( (.not. n_tasks>1) .or. (my_task == 0) ) then
        call DGEMM('N', 'T', 3*n_atoms, 3*n_atoms, 3*n_atoms, 1.d0, &
                   mm_helper, 3*n_atoms, mm_helper, 3*n_atoms, 0.d0, &
                   omegas_p, 3*n_atoms)
    endif
    deallocate( mm_helper )
    if (n_tasks > 1) call broadcast(omegas_p)
    kernel(:, :, :) = 0.d0
    pre(:) = 0.d0
    allocate( om_oo_inv(3*(n_atoms-1), 3*(n_atoms-1)) )
    allocate( om_helper(3, 3*(n_atoms-1)) )
    allocate( om_oo(3*(n_atoms-1), 3*(n_atoms-1)) )
    allocate( om_so(3, 3*(n_atoms-1)) )
    do i_atom=1, n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        self(:) = (/ (3*(i_atom-1)+i, i = 1, 3) /)
        other(:) = (/ (i, i = 1, 3*(i_atom-1)), &
                      (i, i = 3*i_atom+1, 3*n_atoms) /)
        om_oo = omegas_p(other,other)
        om_so = omegas_p(self, other)
        om_ss = omegas_p(self, self)
        om_oo_inv = inverted(om_oo)
        call DGEMM('N', 'N', 3, 3*(n_atoms-1), 3*(n_atoms-1), 1.d0, &
                  om_so, 3, om_oo_inv, 3*(n_atoms-1), 0.d0, om_helper, 3)
        om_so = omegas_p(other, self)
        call DGEMM('N', 'N', 3, 3, 3*(n_atoms-1), 1.d0, om_helper, 3, &
                   om_so, 3*(n_atoms-1), 0.d0, om_3_3, 3)

        kernel(:, :, i_atom) = masses(i_atom)*( om_ss - om_3_3 )

        evals_oo = sdiagonalized(om_oo)
        pre(i_atom) = charges(i_atom)*(masses(i_atom)/pi)**(3.d0/2)*&
                      sqrt(product(omega_int(:3*n_atoms-3)/evals_oo)*&
                           product(omega_int(3*n_atoms-2:)))
    enddo
    deallocate( omegas_p )
    deallocate( om_oo_inv )
    deallocate( om_helper )
    deallocate( om_oo )
    deallocate( om_so )
    call sync_sum(kernel)
    call sync_sum(pre)

    rho(:) = 0.d0
    do i_pt=1, size(gridpts, 1)
        if (my_task /= modulo(i_pt, n_tasks)) cycle
        do i_atom=1, n_atoms
            rdiff(:) = gridpts(i_pt, :)-xyz(i_atom, :)
            forall (i_xyz = 1:3, j_xyz = 1:3)
                rdiffsq(i_xyz, j_xyz) = rdiff(i_xyz)*rdiff(j_xyz)
            end forall
            factor(i_atom) = sum(kernel(:, :, i_atom)*rdiffsq(:, :))
        enddo
        rho(i_pt) = sum(pre*exp(-factor))
    enddo
    call sync_sum(rho)

end function


function eval_mbd_drho_int_nonint_io(gridpts, xyz, charges, masses, &
                                     omega_int, omegas0, fname_modes) &
                                     result(drho)
    real(8),          intent(in) :: gridpts(:, :), &
                                     xyz(:, :), &
                                     charges(size(xyz, 1)), &
                                     masses(size(xyz, 1)), &
                                     omega_int(3*size(xyz, 1)), &
                                     omegas0(size(xyz, 1))
    character(len=*), intent(in) :: fname_modes

    real(8) :: drho(size(gridpts, 1))

    integer :: i_pt, i_atom, n_atoms, i, i_xyz, j_xyz
    integer :: self(3), other(3*(size(xyz, 1)-1))
    real(8) :: pre(size(xyz, 1)), &
                factor(size(xyz, 1)), &
                rdiffsq(3, 3), &
                kernel(3, 3, size(xyz, 1)), &
                rdiff(3), om_3_3(3,3), rho_ni_ipt, &
                evals_oo(3*(size(xyz, 1)-1)), &
                om_ss(3, 3)

    real(8), allocatable :: modes(:,:), mm_helper(:,:), omegas_p(:,:), &
                             om_oo_inv(:,:), om_helper(:,:), pre_ni(:), &
                             kernel_ni(:), rsq(:), om_oo(:,:), &
                             om_so(:,:), om_os(:,:)


    n_atoms = size(xyz, 1)
    allocate( modes(3*n_atoms, 3*n_atoms) )
    call read_mbd_modes(modes, fname_modes)
    allocate( mm_helper(3*n_atoms, 3*n_atoms) )
    mm_helper(:,:) = 0.d0
    do i_atom=1, 3*n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        mm_helper(:,i_atom) = sqrt(omega_int(i_atom))*modes(:,i_atom)
    enddo
    deallocate( modes )
    call sync_sum(mm_helper)
    allocate( omegas_p(3*n_atoms, 3*n_atoms) )
    if ( (.not. n_tasks>1) .or. (my_task == 0) ) then
        call DGEMM('N', 'T', 3*n_atoms, 3*n_atoms, 3*n_atoms, 1.d0, &
                   mm_helper, 3*n_atoms, mm_helper, 3*n_atoms, 0.d0, &
                   omegas_p, 3*n_atoms)
    endif
    deallocate( mm_helper )
    if (n_tasks > 1) call broadcast(omegas_p)
    kernel(:, :, :) = 0.d0
    pre(:) = 0.d0
    allocate( om_oo_inv(3*(n_atoms-1), 3*(n_atoms-1)) )
    allocate( om_helper(3, 3*(n_atoms-1)) )
    allocate( om_oo(3*(n_atoms-1), 3*(n_atoms-1)) )
    allocate( om_so(3, 3*(n_atoms-1)) )
    allocate( om_os(3*(n_atoms-1), 3) )
    do i_atom=1, n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        self(:) = (/ (3*(i_atom-1)+i, i = 1, 3) /)
        other(:) = (/ (i, i = 1, 3*(i_atom-1)), &
                      (i, i = 3*i_atom+1, 3*n_atoms) /)
        om_oo = omegas_p(other,other)
        om_so = omegas_p(self, other)
        om_os = omegas_p(other, self)
        om_ss = omegas_p(self, self)
        om_oo_inv = inverted(om_oo)
        call DGEMM('N', 'N', 3, 3*(n_atoms-1), 3*(n_atoms-1), 1.d0, &
                  om_so, 3, om_oo_inv, 3*(n_atoms-1), 0.d0, om_helper, 3)
        call DGEMM('N', 'N', 3, 3, 3*(n_atoms-1), 1.d0, om_helper, 3, &
                   om_so, 3*(n_atoms-1), 0.d0, om_3_3, 3)

        kernel(:, :, i_atom) = masses(i_atom)*( om_ss - om_3_3 )

        evals_oo = sdiagonalized(om_oo)
        pre(i_atom) = charges(i_atom)*(masses(i_atom)/pi)**(3.d0/2)*&
                      sqrt(product(omega_int(:3*n_atoms-3)/evals_oo)*&
                           product(omega_int(3*n_atoms-2:)))
    enddo
    deallocate( omegas_p )
    deallocate( om_oo_inv )
    deallocate( om_helper )
    deallocate( om_oo )
    deallocate( om_so )
    deallocate( om_os )
    call sync_sum(kernel)
    call sync_sum(pre)

    allocate( pre_ni(n_atoms) )
    allocate( kernel_ni(n_atoms) )
    allocate( rsq(n_atoms) )
    pre_ni = charges*(masses*omegas0/pi)**(3.d0/2)
    kernel_ni = masses*omegas0

    drho(:) = 0.d0
    do i_pt=1, size(gridpts, 1)
        if (my_task /= modulo(i_pt, n_tasks)) cycle
        do i_atom=1, n_atoms
            rdiff(:) = gridpts(i_pt, :)-xyz(i_atom, :)
            forall (i_xyz = 1:3, j_xyz = 1:3)
                rdiffsq(i_xyz, j_xyz) = rdiff(i_xyz)*rdiff(j_xyz)
            end forall
            factor(i_atom) = sum(kernel(:, :, i_atom)*rdiffsq(:, :))
            rsq(i_atom) = sum(rdiff*rdiff)
        enddo
        rho_ni_ipt = sum(pre_ni*exp(-kernel_ni*rsq))
        drho(i_pt) = sum(pre*exp(-factor)) - rho_ni_ipt
    enddo
    call sync_sum(drho)
    deallocate( rsq )
    deallocate( kernel_ni )
    deallocate( pre_ni )

end function


function nbody_coeffs(k, m, N) result(a)
    integer, intent(in) :: k, m, N
    integer :: a

    integer :: i

    a = 1
    do i = N-m+1, N-k
        a = a*i
    end do
    do i = 1, m-k
        a = a/i
    end do
end function nbody_coeffs


subroutine get_interaction_modes(mode, version, xyz, alpha_0, omega, &
                                 R_vdw, beta, a, C6, indices_subsystem1, &
                                 indices_subsystem2, ene_full, unit_cell, &
                                 k_grid, prefix_mbd_modes)
    !! calculates contribution to full MBD mode, which describes
    !! interaction between subsytems
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:,:), alpha_0(size(xyz, 1)), &
                                      omega(size(xyz,1)), &
                                      R_vdw(size(xyz, 1)), &
                                      beta, a, &
                                      C6(size(xyz, 1))
    integer, intent(in) :: indices_subsystem1(:), &
                                      indices_subsystem2(:)
    real(8), intent(in) :: ene_full(3*size(xyz,1))

    real(8), intent(in), optional :: unit_cell(3,3), k_grid(:,:)
    character(len=*), intent(in), optional :: prefix_mbd_modes

    real(8), allocatable :: modes_int(:,:)
    complex(8), allocatable :: modes_int_c(:,:)
    character(len=50) :: prefix_in
    integer :: ik
    logical :: is_crystal
    real(8) :: k_point(3)


    if (present(prefix_mbd_modes)) then
        prefix_in = prefix_mbd_modes
    else
        prefix_in = "mbd_eigenmodes"
    endif
    is_crystal = is_in('C', mode)
    if (.not. is_crystal) then
        allocate( modes_int(3*size(xyz, 1), 3*size(xyz, 1)) )
        call get_single_int_modes(mode, version, xyz, alpha_0, omega, &
                         R_vdw, beta, a, C6, indices_subsystem1, &
                         indices_subsystem2, trim(prefix_in)//".out", &
                         ene_full, modes_int)

        call DWRITEEIGVEC(modes_int, "mbd_interactionmodes.out")
        deallocate(modes_int)
    else
        allocate( modes_int_c(3*size(xyz, 1), 3*size(xyz, 1)) )
        do ik = 1, size(k_grid, 1)
            modes_int_c(:,:) = cmplx(0.d0, 0.d0, 8)
            k_point = k_grid(ik,:)
            call get_single_reciprocal_int_modes(mode, version, xyz, &
                   alpha_0, omega, R_vdw, beta, a, C6, unit_cell, &
                   k_point, indices_subsystem1, indices_subsystem2, &
                   trim(prefix_in)//"_kpt"//trim(tostr(ik))//".out", &
                   ene_full, modes_int_c)

            call ZWRITEEIGVEC(modes_int_c, &
                 "mbd_interactionmodes_kpt"//trim(tostr(ik))//".out",&
                 kpoint=k_point)
        enddo
        deallocate(modes_int_c)
    endif

end subroutine get_interaction_modes


subroutine get_single_int_modes(mode, version, xyz, alpha_0, omega, &
                                R_vdw, beta, a, C6, indices_subsystem1, &
                                indices_subsystem2, fname_modes, &
                                ene_full, modes_int)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:,:), alpha_0(size(xyz, 1)), &
                                     omega(size(xyz,1)), &
                                     R_vdw(size(xyz, 1)), &
                                     beta, a, &
                                     C6(size(xyz, 1))
    integer, intent(in) :: indices_subsystem1(:), &
                                     indices_subsystem2(:)
    character(len=*), intent(in) :: fname_modes
    real(8), intent(in) :: ene_full(3*size(xyz, 1))

    real(8), intent(inout) :: &
                               modes_int(3*size(xyz, 1), 3*size(xyz, 1))

    real(8), allocatable :: modes(:,:), h_int(:,:)
    integer :: i, nModes


    nModes = 3*size(xyz,1)
    allocate( h_int(nModes,nModes) )
    call build_h_int(mode, version, xyz, alpha_0, omega, R_vdW, beta, &
                     a, C6, indices_subsystem1, indices_subsystem2, &
                     h_int=h_int)

    allocate( modes(nModes, nModes) )
    call read_mbd_modes(modes, fname_modes)
    call DGEMM('N','N', nModes, nModes, nModes, 1.d0, h_int, nModes, &
               modes, nModes, 0.d0, modes_int, nModes)

    deallocate(modes)
    do i=1,nModes
        modes_int(i,:) = modes_int(i,:)/ene_full(i)
    enddo

end subroutine get_single_int_modes


subroutine get_single_reciprocal_int_modes(mode, version, xyz, alpha_0, &
                         omega, R_vdw, beta, a, C6, unit_cell, k_point, &
                         indices_subsystem1, indices_subsystem2, &
                         fname_modes, ene_full, modes_int)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:,:), alpha_0(size(xyz, 1)), &
                                     omega(size(xyz,1)), &
                                     R_vdw(size(xyz, 1)), &
                                     beta, a, &
                                     C6(size(xyz, 1)), &
                                     unit_cell(3,3), &
                                     k_point(3)
    integer, intent(in) :: indices_subsystem1(:), &
                                     indices_subsystem2(:)
    character(len=*), intent(in) :: fname_modes
    real(8), intent(in) :: ene_full(3*size(xyz, 1))

    complex(8), intent(inout) :: modes_int(3*size(xyz, 1), 3*size(xyz, 1))

    complex(8), allocatable :: modes(:,:), h_int(:,:)
    integer :: i, nModes


    nModes = 3*size(xyz,1)
    allocate( h_int(nModes,nModes) )
    call build_h_int(mode, version, xyz, alpha_0, omega, R_vdW, beta, &
                     a, C6, indices_subsystem1, indices_subsystem2, &
                     unit_cell=unit_cell, k_point=k_point, h_int_c=h_int)

    allocate( modes(nModes,nModes) )
    call read_mbd_modes_reciprocal(modes, fname_modes)
    call ZGEMM('N','N', nModes, nModes, nModes, 1.d0, h_int, nModes, &
               modes, nModes, 0.d0, modes_int, nModes)

    deallocate(modes)
    do i=1,nModes
        modes_int(i, :) = modes_int(i,:)/ene_full(i)
    enddo

end subroutine get_single_reciprocal_int_modes


subroutine build_h_int(mode, version, xyz, alpha, omega, R_vdw, beta, &
                       a, C6, indices_subsystem1, indices_subsystem2, &
                       unit_cell, k_point, h_int, h_int_c)
    !! returns block off-diagonal part of Hamiltonian describing
    !! interaction between subsystems:
    !!              /   0     V_{21} \
    !!   H_{int} = |                  |
    !!              \ V_{12}    0    /

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha(size(xyz, 1)), &
                                      omega(size(xyz,1)), &
                                      R_vdw(size(xyz, 1)), &
                                      beta, a, &
                                      C6(size(xyz, 1))
    integer, intent(in) :: indices_subsystem1(:), &
                                      indices_subsystem2(:)

    real(8), intent(in), optional :: unit_cell(3, 3), k_point(3)

    real(8), intent(inout), optional :: &
                                 h_int(3*size(xyz, 1), 3*size(xyz, 1))
    complex(8), intent(inout), optional :: &
                                 h_int_c(3*size(xyz, 1), 3*size(xyz, 1))

    real(8) :: Tpp(3, 3), R_cell(3), r(3), r_norm, R_vdw_ij, &
                         C6_ij, sigma_ij, volume, ewald_alpha, &
                         real_space_cutoff, wwsqrtaa
    complex(8) :: Tpp_c(3, 3)
    character(len=1) :: parallel_mode
    integer :: i_atom, j_atom, i_cell, idx_cell(3), &
                         range_cell(3), i, j, ii, jj
    logical :: is_crystal, is_parallel, is_reciprocal, &
                         is_low_dim, mute, do_ewald


    is_parallel = is_in('P', mode)
    is_reciprocal = is_in('R', mode)
    is_crystal = is_in('C', mode) .or. is_reciprocal
    is_low_dim = any(param_vacuum_axis)
    do_ewald = .false.
    mute = is_in('M', mode)
    if (is_parallel) then
        parallel_mode = 'A' ! atoms
        if (is_crystal .and. size(xyz, 1) < n_tasks) then
            parallel_mode = 'C' ! cells
        end if
    else
        parallel_mode = ''
    end if

    ! MPI code begin
    if (is_parallel) then
        ! will be restored by syncing at the end
        if (is_reciprocal) then
            h_int_c = h_int_c/n_tasks
        else
            h_int = h_int/n_tasks
        end if
    end if
    ! MPI code end
    if (is_crystal) then
        if (is_low_dim) then
            real_space_cutoff = param_dipole_low_dim_cutoff
        else if (param_ewald_on) then
            do_ewald = .true.
            volume = max(abs(dble(product(diagonalized(unit_cell)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1.d0/3)
            real_space_cutoff = 6.d0/ewald_alpha*param_ewald_real_cutoff_scaling
            call print_log('Ewald: using alpha = '//trim(tostr(ewald_alpha)) &
                //', real cutoff = '//trim(tostr(real_space_cutoff)), mute)
        else
            real_space_cutoff = param_dipole_cutoff
        end if
        range_cell = supercell_circum(unit_cell, real_space_cutoff)
    else
        range_cell = (/ 0, 0, 0 /)
    end if
    if (is_crystal) then
        call print_log('Ewald: summing real part in cell vector range of '&
            //trim(tostr(1+2*range_cell(1)))//'x' &
            //trim(tostr(1+2*range_cell(2)))//'x' &
            //trim(tostr(1+2*range_cell(3))), mute)
    end if
    call ts(11)
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        ! MPI code begin
        if (parallel_mode == 'C') then
            if (my_task /= modulo(i_cell, n_tasks)) cycle
        end if
        ! MPI code end
        if (is_crystal) then
            R_cell = matmul(idx_cell, unit_cell)
        else
            R_cell(:) = 0.d0
        end if
        do ii=1, size(indices_subsystem1)
            i_atom = indices_subsystem1(ii)
            ! MPI code begin
            if (parallel_mode == 'A') then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            ! MPI code end
            do jj=1, size(indices_subsystem2)
                j_atom = indices_subsystem2(jj)
                r = xyz(i_atom, :)-xyz(j_atom, :)-R_cell
                r_norm = sqrt(sum(r*r))
                if (is_crystal .and. r_norm > real_space_cutoff) cycle
                R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
                sigma_ij = sqrt(sum(get_sigma_selfint( &
                        alpha((/ i_atom , j_atom /)))**2))
                C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha(i_atom), alpha(j_atom))
                select case (version)
                    case ("bare")
                        Tpp = T_bare(r)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,erf")
                        Tpp = T_erf_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,fermi")
                        Tpp = T_fermi_coulomb(r, beta*R_vdw_ij, a)
                    case ("1mexp,dip")
                        Tpp = damping_1mexp(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("erf,dip")
                        Tpp = damping_erf(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi^2,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)**2*T_bare(r)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("1mexp,dip,gg")
                        Tpp = (1.d0-damping_1mexp(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("erf,dip,gg")
                        Tpp = (1.d0-damping_erf(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("fermi,dip,gg")
                        Tpp = (1.d0-damping_fermi(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                end select
                if (do_ewald) then
                    Tpp = Tpp + T_erfc(r, ewald_alpha) - T_bare(r)
                end if
                if (is_reciprocal) then
                    Tpp_c = Tpp*exp(-cmplx(0.d0, 1.d0, 8)*( &
                                    dot_product(k_point, r)))
                end if
                i = 3*i_atom-2
                j = 3*j_atom-2
                if (is_reciprocal) then
                    h_int_c(i:i+2, j:j+2) = h_int_c(i:i+2, j:j+2) + Tpp_c
                    h_int_c(j:j+2, i:i+2) = h_int_c(j:j+2, i:i+2) +&
                                            transpose(Tpp_c)
                else
                    h_int(i:i+2, j:j+2) = h_int(i:i+2, j:j+2) + Tpp
                    h_int(j:j+2, i:i+2) = h_int(j:j+2, i:i+2) +&
                                          transpose(Tpp)
                end if
            end do ! jj (j_atom)
        end do ! ii (i_atom)
    end do ! i_cell
    call ts(-11)
    ! MPI code begin
    if (is_parallel) then
        if (is_reciprocal) then
            call sync_sum(h_int_c)
        else
            call sync_sum(h_int)
        end if
    end if
    ! MPI code end
    if (do_ewald) then
        call add_ewald_h_int_parts(mode, xyz, unit_cell, ewald_alpha, &
                              indices_subsystem1, indices_subsystem2, &
                              k_point, h_int, h_int_c)
    end if

    !! update off-diagonal blocks such that H_int = w*w*sqrt(a*a)*T
    do ii=1, size(indices_subsystem1)
        i_atom = indices_subsystem1(ii)
        i = 3*i_atom-2
        do jj=1, size(indices_subsystem2)
            j_atom = indices_subsystem2(jj)
            wwsqrtaa = omega(i_atom)*omega(j_atom)*&
                       sqrt(alpha(i_atom)*alpha(j_atom))
            j = 3*j_atom-2
            if (is_reciprocal) then
                h_int_c(i:i+2, j:j+2) = wwsqrtaa*h_int_c(i:i+2, j:j+2)
                h_int_c(j:j+2, i:i+2) = wwsqrtaa*h_int_c(j:j+2, i:i+2)
            else
                h_int(i:i+2, j:j+2) = wwsqrtaa*h_int(i:i+2, j:j+2)
                h_int(j:j+2, i:i+2) = wwsqrtaa*h_int(j:j+2, i:i+2)
            endif
        enddo
    enddo

end subroutine build_h_int


subroutine add_ewald_h_int_parts(mode, xyz, unit_cell, alpha, &
                                 indices_subsystem1, indices_subsystem2,&
                                 k_point, h_int, h_int_c)
    character(len=*), intent(in) :: mode
    real(8), intent(in) :: xyz(:, :), unit_cell(3, 3), alpha
    integer, intent(in) :: indices_subsystem1(:), &
                                      indices_subsystem2(:)

    real(8), intent(in), optional :: k_point(3)
    real(8), intent(inout), optional :: &
                                 h_int(3*size(xyz, 1), 3*size(xyz, 1))
    complex(8), intent(inout), optional :: &
                                 h_int_c(3*size(xyz, 1), 3*size(xyz, 1))

    logical :: is_parallel, is_reciprocal, mute, do_surface
    real(8) :: rec_unit_cell(3, 3), volume, G_vector(3), r(3), &
                        k_total(3), k_sq, rec_space_cutoff, Tpp(3, 3), &
                        k_prefactor(3, 3), elem
    complex(8) :: Tpp_c(3, 3)
    integer :: i_atom, j_atom, i, j, i_xyz, j_xyz, ii, jj, &
                        idx_G_vector(3), i_G_vector, range_G_vector(3)
    character(len=1) :: parallel_mode


    is_parallel = is_in('P', mode)
    is_reciprocal = is_in('R', mode)
    mute = is_in('M', mode)
    if (is_parallel) then
        parallel_mode = 'A' ! atoms
        if (size(xyz, 1) < n_tasks) then
            parallel_mode = 'G' ! G vectors
        end if
    else
        parallel_mode = ''
    end if

    ! MPI code begin
    if (is_parallel) then
        ! will be restored by syncing at the end
        if (is_reciprocal) then
            h_int_c = h_int_c/n_tasks
        else
            h_int = h_int/n_tasks
        end if
    end if
    ! MPI code end
    rec_unit_cell = 2*pi*inverted(transpose(unit_cell))
    volume = abs(dble(product(diagonalized(unit_cell))))
    rec_space_cutoff = 10.d0*alpha*param_ewald_rec_cutoff_scaling
    range_G_vector = supercell_circum(rec_unit_cell, rec_space_cutoff)
    call print_log('Ewald: using reciprocal cutoff = ' &
        //trim(tostr(rec_space_cutoff)), mute)
    call print_log('Ewald: summing reciprocal part in G vector range of ' &
        //trim(tostr(1+2*range_G_vector(1)))//'x' &
        //trim(tostr(1+2*range_G_vector(2)))//'x' &
        //trim(tostr(1+2*range_G_vector(3))), mute)
    call ts(12)
    idx_G_vector = (/ 0, 0, -1 /)
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        ! MPI code begin
        if (parallel_mode == 'G') then
            if (my_task /= modulo(i_G_vector, n_tasks)) cycle
        end if
        ! MPI code end
        G_vector = matmul(idx_G_vector, rec_unit_cell)
        if (is_reciprocal) then
            k_total = k_point+G_vector
        else
            k_total = G_vector
        end if
        k_sq = sum(k_total*k_total)
        if (sqrt(k_sq) > rec_space_cutoff) cycle
        k_prefactor(:, :) = 4*pi/volume*exp(-k_sq/(4*alpha*alpha))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz)*&
                k_total(i_xyz)*k_total(j_xyz)/k_sq
        do ii=1, size(indices_subsystem1)
            i_atom = indices_subsystem1(ii)
            ! MPI code begin
            if (parallel_mode == 'A') then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            ! MPI code end
            do jj=1, size(indices_subsystem2)
                j_atom = indices_subsystem2(jj)
                r = xyz(i_atom, :)-xyz(j_atom, :)
                if (is_reciprocal) then
                    Tpp_c = k_prefactor*exp(cmplx(0.d0, 1.d0, 8)*&
                            dot_product(G_vector, r))
                else
                    Tpp = k_prefactor*cos(dot_product(G_vector, r))
                end if
                i = 3*i_atom-2
                j = 3*j_atom-2
                if (is_reciprocal) then
                    h_int_c(i:i+2, j:j+2) = h_int_c(i:i+2, j:j+2) + Tpp_c
                    h_int_c(j:j+2, i:i+2) = h_int_c(j:j+2, i:i+2) +&
                                            transpose(Tpp_c)
                else
                    h_int(i:i+2, j:j+2) = h_int(i:i+2, j:j+2) + Tpp
                    h_int(j:j+2, i:i+2) = h_int(j:j+2, i:i+2) +&
                                          transpose(Tpp)
                end if
            end do ! jj (j_atom)
        end do ! ii (i_atom)
    end do ! i_G_vector
    ! MPI code begin
    if (is_parallel) then
        if (is_reciprocal) then
            call sync_sum(h_int_c)
        else
            call sync_sum(h_int)
        end if
    end if
    ! MPI code end
    do_surface = .true.
    if (present(k_point)) then
        k_sq = sum(k_point*k_point)
        if (sqrt(k_sq) > 1.d-15) then
            do_surface = .false.
            do ii=1, size(indices_subsystem1)
                i_atom = indices_subsystem1(ii)
                do jj=1, size(indices_subsystem2)
                    j_atom = indices_subsystem2(jj)
                    do i_xyz=1, 3
                        do j_xyz=1, 3
                            i = 3*(i_atom-1)+i_xyz
                            j = 3*(j_atom-1)+j_xyz
                            elem = 4*pi/volume*k_point(i_xyz)*&
                                   k_point(j_xyz)/k_sq*&
                                   exp(-k_sq/(4*alpha*alpha))
                            if (is_reciprocal) then
                                h_int_c(i, j) = h_int_c(i, j) + elem
                                h_int_c(j, i) = h_int_c(j, i) + elem
                            else
                                h_int(i, j) = h_int(i, j) + elem
                                h_int(j, i) = h_int(j, i) + elem
                            end if ! is_reciprocal
                        end do ! j_xyz
                    end do ! i_xyz
                end do ! jj (j_atom)
            end do ! ii (i_atom)
        end if ! k_sq >
    end if ! k_point present
    if (do_surface) then ! surface energy
        do ii=1, size(indices_subsystem1)
            do jj=1, size(indices_subsystem2)
                do i_xyz = 1, 3
                    i = 3*(indices_subsystem1(ii)-1)+i_xyz
                    j = 3*(indices_subsystem2(jj)-1)+i_xyz
                    if (is_reciprocal) then
                        h_int_c(i, j) = h_int_c(i, j)+4*pi/(3*volume)
                        h_int_c(j, i) = h_int_c(i, j)
                    else
                        h_int(i, j) = h_int(i, j)+4*pi/(3*volume)
                        h_int(j, i) = h_int(i, j)
                    end if
                end do ! i_xyz
            end do ! jj (j_atom)
        end do ! ii (i_atom)
    end if
    call ts(-12)
end subroutine add_ewald_h_int_parts


subroutine get_onsite_modes(mode, version, xyz, alpha_0, omega, R_vdw, &
                            beta, a, C6, indices_subsystem1, &
                            indices_subsystem2, ene_full, unit_cell, &
                            k_grid, prefix_mbd_modes)
    !! calculates contribution to full MBD mode, which describes
    !! interaction between subsytems
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:,:), alpha_0(size(xyz, 1)), &
                                      omega(size(xyz,1)), &
                                      R_vdw(size(xyz, 1)), &
                                      beta, a, &
                                      C6(size(xyz, 1))
    integer, intent(in) :: indices_subsystem1(:), &
                                      indices_subsystem2(:)
    real(8), intent(in) :: ene_full(3*size(xyz, 1))

    real(8), intent(in), optional :: unit_cell(3,3), k_grid(:,:)
    character(len=*), intent(in), optional :: prefix_mbd_modes

    real(8), allocatable :: modes_on(:,:)
    complex(8), allocatable :: modes_on_c(:,:)
    character(len=100) :: prefix_in
    logical :: is_crystal
    integer :: ik
    real(8) :: k_point(3)



    if (present(prefix_mbd_modes)) then
        prefix_in = prefix_mbd_modes
    else
        prefix_in = "mbd_eigenmodes"
    endif
    is_crystal = is_in('C', mode)
    if (.not. is_crystal) then
        allocate( modes_on(3*size(xyz, 1), 3*size(xyz, 1)) )
        call get_single_on_modes(mode, version, xyz, alpha_0, omega, &
                 R_vdw, beta, a, C6, ene_full, indices_subsystem1, &
                 indices_subsystem2, trim(prefix_in)//".out", &
                 modes_on)

        call DWRITEEIGVEC(modes_on, "mbd_onsitemodes.out")
        deallocate(modes_on)
    else
        allocate( modes_on_c(3*size(xyz, 1), 3*size(xyz, 1)) )
        do ik = 1, size(k_grid, 1)
            k_point = k_grid(ik,:)
            modes_on_c(:,:) = cmplx(0.d0, 0.d0, 8)
            call get_single_reciprocal_on_modes(mode, version, xyz, &
                     alpha_0, omega, R_vdw, beta, a, C6, &
                     ene_full, unit_cell, k_point, &
                     indices_subsystem1, indices_subsystem2, &
                     trim(prefix_in)//"_kpt"//trim(tostr(ik))//".out", &
                     modes_on_c)

            call ZWRITEEIGVEC(modes_on_c, &
                         "mbd_onsitemodes_kpt"//trim(tostr(ik))//".out",&
                         kpoint=k_point)
        enddo
        deallocate(modes_on_c)
    endif

end subroutine get_onsite_modes


subroutine get_single_on_modes(mode, version, xyz, alpha_0, omega, &
                               R_vdw, beta, a, C6, ene_full, &
                               indices_subsystem1, indices_subsystem2, &
                               fname_modes, modes_on)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:,:), alpha_0(size(xyz, 1)), &
                                     omega(size(xyz,1)), &
                                     R_vdw(size(xyz, 1)), &
                                     beta, a, &
                                     C6(size(xyz, 1)), &
                                     ene_full(3*size(xyz,1))
    integer, intent(in) :: indices_subsystem1(:), &
                                     indices_subsystem2(:)
    character(len=*), intent(in) :: fname_modes

    real(8), intent(inout) :: modes_on(3*size(xyz,1),3*size(xyz,1))

    real(8), allocatable :: modes(:,:)
    integer :: i, nModes


    nModes = 3*size(xyz,1)
    call get_single_int_modes(mode, version, xyz, alpha_0, omega, &
                              R_vdw, beta, a, C6, indices_subsystem1, &
                              indices_subsystem2, fname_modes, &
                              ene_full, modes_on)

    allocate( modes(nModes,nModes) )
    call read_mbd_modes(modes, fname_modes)
    do i=1, nModes
        modes_on(i,:) = modes(i,:)-modes_on(i,:)
    enddo
    deallocate(modes)

end subroutine get_single_on_modes


subroutine get_single_reciprocal_on_modes(mode, version, xyz, alpha_0, &
                     omega, R_vdw, beta, a, C6, ene_full, unit_cell, &
                     k_point, indices_subsystem1, indices_subsystem2, &
                     fname_modes, modes_on)
    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:,:), alpha_0(size(xyz, 1)), &
                                     omega(size(xyz,1)), &
                                     R_vdw(size(xyz, 1)), &
                                     beta, a, &
                                     C6(size(xyz, 1)), &
                                     ene_full(3*size(xyz,1)), &
                                     unit_cell(3,3), &
                                     k_point(3)
    integer, intent(in) :: indices_subsystem1(:), &
                                     indices_subsystem2(:)
    character(len=*), intent(in) :: fname_modes

    complex(8), intent(inout) :: modes_on(3*size(xyz,1),3*size(xyz,1))

    complex(8), allocatable :: modes(:,:)
    integer :: i, nModes


    nModes = 3*size(xyz,1)
    call get_single_reciprocal_int_modes(mode, version, xyz, alpha_0, &
                       omega, R_vdw, beta, a, C6, unit_cell, k_point, &
                       indices_subsystem1, indices_subsystem2, &
                       fname_modes, ene_full, modes_on)

    allocate( modes(nModes,nModes) )
    call read_mbd_modes_reciprocal(modes, fname_modes)
    do i=1, nModes
        modes_on(i,:) = modes(i,:)-modes_on(i,:)
    enddo
    deallocate(modes)

end subroutine get_single_reciprocal_on_modes


subroutine read_mbd_modes(modes, filenam)
    real(8), intent(inout) :: modes(:,:)
    character(len=*), intent(in) :: filenam
    integer :: fID_V, nModes, buffer_n, i, iMode


    fID_V = get_FileID()
    open( file=trim(filenam), unit=fID_V, status='old', &
          form='unformatted' )
    read(fID_V) nModes, buffer_n
    if (nModes /= size(modes,1)) then
        call print_error(" Size of eigenmode input file does not match system.")
    endif
    !! read modes [NOTE: in output format mode i = modes(i,:)]
    do i=1, nModes
        do iMode=1, nModes
            read(fID_V) modes(iMode, i)
        enddo
    enddo
    close(fID_V)

end subroutine read_mbd_modes


subroutine read_mbd_modes_reciprocal(modes, filenam)
    complex(8), intent(inout) :: modes(:,:)
    character(len=*), intent(in) :: filenam
    real(8) :: buffer_k(3)
    integer :: fID_V, nModes, buffer_n, i, iMode


    fID_V = get_FileID()
    open( file=trim(filenam), unit=fID_V, status='old', &
          form='unformatted' )
    read(fID_V) nModes, buffer_n
    read(fID_V) buffer_k

    if (nModes /= size(modes,1)) then
        call print_error(" Size of eigenmode input file does not match system.")
    endif
    !! read modes [NOTE: in output format mode i = modes(i,:)]
    do i=1, nModes
        do iMode=1, nModes
            read(fID_V) modes(iMode, i)
        enddo
    enddo
    close(fID_V)

end subroutine read_mbd_modes_reciprocal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    SCS/MBD (LAPACK/MPI) HELPER FUNCTIONS     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function contract_polarizability(alpha_3n_3n) result(alpha_n)
    real(8), intent(in) :: alpha_3n_3n(:, :)
    real(8) :: alpha_n(size(alpha_3n_3n, 1)/3)

    integer :: i_atom, i_xyz, j_xyz, dim_3n
    real(8) :: alpha_3_3(3, 3), alpha_diag(3)

    dim_3n = size(alpha_3n_3n, 1)
    do i_atom = 1, size(alpha_n)
        forall (i_xyz = 1:3, j_xyz = 1:3) alpha_3_3(i_xyz, j_xyz) &
                = sum(alpha_3n_3n(i_xyz:dim_3n:3, 3*(i_atom-1)+j_xyz))
        alpha_diag = sdiagonalized(alpha_3_3)
        alpha_n(i_atom) = sum(alpha_diag)/3
    end do
end function contract_polarizability


subroutine invert_ge_dble_(A)
    real(8), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    real(8), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    real(8) :: n_work_arr_optim
    integer :: error_flag

    n = size(A, 1)
    if (n == 0) return
    call DGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    call DGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


subroutine invert_ge_cmplx_(A)
    complex(8), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    complex(8), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    complex(8) :: n_work_arr_optim
    integer :: error_flag

    n = size(A, 1)
    call ZGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    call ZGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(dble(n_work_arr_optim))
    allocate (work_arr(n_work_arr))
    call ZGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_error( &
            "Matrix inversion failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


function inverted(A) result(A_inv)
    real(8), intent(in) :: A(:, :)
    real(8) :: A_inv(size(A, 1), size(A, 2))

    A_inv = A
    call invert(A_inv)
end function


subroutine diagonalize_sym_dble_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    real(8), intent(inout) :: A(:, :)
    real(8), intent(out) :: eigs(size(A, 1))

    real(8), allocatable :: work_arr(:), vecs(:,:)
    integer :: n, vl, vu, il, iu, nvals
    integer, allocatable :: iwork_arr(:), isupp(:)
    integer :: error_flag, lwork, liwork

    n = size(A, 1)
    if ( allocated(work_arr) ) deallocate(work_arr)
    allocate( work_arr(1) )
    selectcase(trim(eigensolver))
        case("qr")        !! QR algorithm: medium fast, minimal memory requirement
            call DSYEV(mode, "U", n, A, n, eigs, work_arr, -1, error_flag)
            lwork = nint(work_arr(1))
            deallocate(work_arr)
            allocate ( work_arr(lwork) )
            call DSYEV(mode, "U", n, A, n, eigs, work_arr, lwork, error_flag)

        case("dandc")     !! Divide & Conquer: fast, higher memory requirement
            if ( allocated(iwork_arr) ) deallocate(iwork_arr)
            allocate( iwork_arr(1) )
            call DSYEVD(mode, "U", n, A, n, eigs, work_arr, -1, &
                        iwork_arr, -1, error_flag)
            lwork = nint(work_arr(1))
            liwork = iwork_arr(1)
            deallocate(work_arr)
            deallocate(iwork_arr)
            allocate (work_arr(lwork))
            allocate (iwork_arr(liwork))
            call DSYEVD(mode, "U", n, A, n, eigs, work_arr, lwork, &
                        iwork_arr, liwork, error_flag)
            if ( allocated(iwork_arr) ) deallocate(iwork_arr)

        case("mrrr")      !! MRRR algorithm: medium memory requirement
            if ( allocated(isupp) ) deallocate(isupp)
            allocate( isupp(2*n), stat=error_flag)
            if ( allocated(vecs) ) deallocate(vecs)
            allocate( vecs(n,n) )
            if ( allocated(iwork_arr) ) deallocate(iwork_arr)
            allocate( iwork_arr(1) )
            call DSYEVR(mode, "A", "U", n, A, n, vl, vu, il, iu, 1.d-10, &
                        nvals, eigs, vecs, n, isupp, work_arr, -1, &
                        iwork_arr, -1, error_flag)
            lwork = nint(work_arr(1))
            liwork = iwork_arr(1)
            deallocate(work_arr)
            deallocate(iwork_arr)
            allocate (work_arr(lwork))
            allocate (iwork_arr(liwork))
            call DSYEVR(mode, "A", "U", n, A, n, vl, vu, il, iu, 1.d-10, &
                        nvals, eigs, vecs, n, isupp, work_arr, lwork, &
                        iwork_arr, liwork, error_flag)
            A = vecs
            deallocate(iwork_arr)
            deallocate(isupp)
            deallocate(vecs)

    endselect
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_log( &
            "DSYEVx failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


function diagonalized_sym_dble_(A, eigvecs) result(eigs)
    real(8), intent(in) :: A(:, :)
    real(8), intent(out), optional, target :: eigvecs(size(A, 1), size(A, 2))
    real(8) :: eigs(size(A, 1))

    real(8), pointer :: eigvecs_p(:, :)
    character(len=1) :: mode

    if (present(eigvecs)) then
        mode = 'V'
        eigvecs_p => eigvecs
    else
        mode = 'N'
        allocate (eigvecs_p(size(A, 1), size(A, 2)))
    end if
    eigvecs_p = A
    call sdiagonalize(mode, eigvecs_p, eigs)
    if (.not. present(eigvecs)) then
        deallocate (eigvecs_p)
    end if
end function


subroutine diagonalize_ge_dble_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    real(8), intent(inout) :: A(:, :)
    complex(8), intent(out) :: eigs(size(A, 1))

    real(8), allocatable :: work_arr(:)
    integer :: n
    real(8) :: n_work_arr
    integer :: error_flag
    real(8) :: eigs_r(size(A, 1)), eigs_i(size(A, 1))
    real(8) :: dummy
    real(8) :: vectors(size(A, 1), size(A, 2))

    n = size(A, 1)
    call DGEEV('N', mode, n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DGEEV('N', mode, n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, work_arr, size(work_arr), error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_log( &
            "DGEEV failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    eigs = cmplx(eigs_r, eigs_i, 8)
    A = vectors
end subroutine


function diagonalized_ge_dble_(A, eigvecs) result(eigs)
    real(8), intent(in) :: A(:, :)
    real(8), intent(out), optional, target :: eigvecs(size(A, 1), size(A, 2))
    complex(8) :: eigs(size(A, 1))

    real(8), pointer :: eigvecs_p(:, :)
    character(len=1) :: mode

    if (present(eigvecs)) then
        mode = 'V'
        eigvecs_p => eigvecs
    else
        mode = 'N'
        allocate (eigvecs_p(size(A, 1), size(A, 2)))
    end if
    eigvecs_p = A
    call diagonalize(mode, eigvecs_p, eigs)
    if (.not. present(eigvecs)) then
        deallocate (eigvecs_p)
    end if
end function


subroutine diagonalize_he_cmplx_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    complex(8), intent(inout) :: A(:, :)
    real(8), intent(out) :: eigs(size(A, 1))

    complex(8), allocatable :: work(:), vecs(:,:)
    real(8), allocatable :: rwork(:)
    integer, allocatable :: iwork(:), isupp(:)
    integer :: vl, vu, il, iu, nvals
    integer :: n, lwork, lrwork, liwork
    integer :: error_flag


    n = size(A, 1)
    if ( allocated(rwork) ) deallocate(rwork)
    if ( allocated(work) ) deallocate(work)
    selectcase(trim(eigensolver))
        case("qr")        !! QR algorithm: less fast, minimal memory consumption
            allocate( rwork(max(1, 3*n-2)) )
            allocate( work(1) )
            call ZHEEV(mode, "U", n, A, n, eigs, work, -1, rwork, &
                       error_flag)
            lwork = nint( dble(work(1)) )
            deallocate(work)
            allocate( work(lwork) )
            call ZHEEV(mode, "U", n, A, n, eigs, work, lwork, rwork, &
                       error_flag)

        case("dandc")     !! Divide & Conquer: fast, high memory consumption
            allocate(work(1))
            allocate(rwork(1))
            if ( allocated(iwork) ) deallocate(iwork)
            allocate(iwork(1))
            call ZHEEVD(mode, "U", n, A, n, eigs, work, -1, rwork, -1, &
                        iwork, -1, error_flag)
            lwork = nint( dble( work(1) ) )
            lrwork = nint( rwork(1) )
            liwork = iwork(1)
            deallocate(work)
            deallocate(rwork)
            deallocate(iwork)
            allocate( work(lwork) )
            allocate( rwork(lrwork) )
            allocate( iwork(liwork) )
            call ZHEEVD(mode, "U", n, A, n, eigs, work, lwork, rwork, &
                        lrwork, iwork, liwork, error_flag)
            deallocate(iwork)

        case("mrrr")      !! MRRR algorithm: fast, medium memory consumption
            allocate(work(1))
            allocate(rwork(1))
            if ( allocated(iwork) ) deallocate(iwork)
            allocate(iwork(1))
            if ( allocated(vecs) ) deallocate(vecs)
            allocate( vecs(n,n) )
            if ( allocated(isupp) ) deallocate(isupp)
            allocate( isupp(2*n) )
            call ZHEEVR(mode, "A", "U", n, A, n, vl, vu, il, iu, 1.d-10, &
                        nvals, eigs, vecs, n, isupp, work, -1, rwork, &
                        -1, iwork, -1, error_flag)
            lwork = nint( dble( work(1) ) )
            lrwork = nint( rwork(1) )
            liwork = iwork(1)
            deallocate(work)
            deallocate(rwork)
            deallocate(iwork)
            allocate (work(lwork))
            allocate (rwork(lrwork))
            allocate (iwork(liwork))
            call ZHEEVR(mode, "A", "U", n, A, n, vl, vu, il, iu, 1.d-10, &
                        nvals, eigs, vecs, n, isupp, work, lwork, rwork, &
                        lrwork, iwork, liwork, error_flag)
            A = vecs
            deallocate(iwork)
            deallocate(isupp)
            deallocate(vecs)

    endselect
    deallocate(rwork)
    deallocate(work)
    if (error_flag /= 0) then
        call print_error( &
            "ZHEEVx failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
end subroutine


subroutine diagonalize_ge_cmplx_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    complex(8), intent(inout) :: A(:, :)
    complex(8), intent(out) :: eigs(size(A, 1))

    complex(8), allocatable :: work(:)
    real(8) :: rwork(2*size(A, 1))
    integer :: n, lwork
    complex(8) :: lwork_arr
    integer :: error_flag
    complex(8) :: dummy
    complex(8) :: vectors(size(A, 1), size(A, 2))

    n = size(A, 1)
    call ZGEEV('N', mode, n, A, n, eigs, dummy, 1, &
        vectors, n, lwork_arr, -1, rwork, error_flag)
    lwork = nint(dble(lwork_arr))
    allocate (work(lwork))
    call ZGEEV('N', mode, n, A, n, eigs, dummy, 1, &
        vectors, n, work, lwork, rwork, error_flag)
    deallocate (work)
    if (error_flag /= 0) then
        call print_log( &
            "ZGEEV failed in module mbd with error code " &
            //trim(tostr(error_flag)))
    endif
    A = vectors
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!  !!!!!!  !!!!!!!  !!      !!!!!!!  !!!!!!  !!!!!!!  !!!!!!  !!   !!
!!      !!      !!   !!  !!      !!   !!  !!  !!  !!   !!  !!      !!  !!
!!!!!!  !!      !!!!!!!  !!      !!!!!!!  !!!!!!  !!!!!!!  !!      !!!!!
    !!  !!      !!   !!  !!      !!   !!  !!      !!   !!  !!      !!  !!
!!!!!!  !!!!!!  !!   !!  !!!!!!  !!   !!  !!      !!   !!  !!!!!!  !!   !!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     vdW(TS) ENERGY     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_ts_energy_lowmem(mode, version, xyz, C6, alpha_0, R_vdw, &
                              s_R, d, unit_cell) result(ene)
!! REDUCED CUSTOMIZABILITY FOR THE SAKE OF MEMORY REQUIREMENTS
    !!!!! INPUT !!!!!
    !! mode:      runtime option 'C': use periodic boundary conditions
    !! version:   which damping function to use
    !! xyz:       atomic positions [a.u.]
    !! C6:        C6 interaction coefficients to use
    !! alpha_0:   static, isotropic atomic polarizabilities
    !! R_vdw:     van der Waals radii [a.u.]
    !! s_R, d:    damping parameters (for range and steepness)
    !! unit_cell: well, unit cell [a.u.]
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT !!!!!
    !! ene:       vdW(TS) dispersion energy [a.u.]
    !!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! C6_ij:           combined C6 interaction coefficient
    !! r/r_norm:        interatomic distance and corresponding vector
    !! f_damp:          damping value for given atoms and distance
    !! R_vdw_ij:        combined van der Waals radius
    !! ene_shell:       contribution to energy per shell of periodic cells
    !! ene_pair:        contribution of current pair of atoms
    !! R_cell:          distance to current periodic image
    !! i/j_X:           loop indices
    !! range_cell:      range of neighboring cells to consider
    !! idx_cell:        index of current neighboring cell
    !! shell_thickness: range of neighboring cells to be considered
    !! is_X:            flags whether to do/use X
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), C6(size(xyz, 1)), &
                                      alpha_0(size(xyz, 1))
    real(8), intent(in), optional :: R_vdw(size(xyz, 1)), s_R, d, &
                                      unit_cell(3, 3)

    real(8) :: ene

    real(8) :: C6_ij, r(3), r_norm, f_damp, R_vdw_ij, ene_shell, &
                ene_pair, R_cell(3)
    integer :: i_shell, i_cell, i_atom, j_atom, range_cell(3), &
                minr_cell(3), idx_cell(3)
    real(8), parameter :: shell_thickness = 10.d0
    logical :: is_crystal, is_parallel


    is_crystal = is_in('C', mode)
    is_parallel = is_in('P', mode)

    ene = 0.d0
    i_shell = 0
    do          ! until energy cutoff criterion is reached
        i_shell = i_shell+1
        ene_shell = 0.d0
        if (is_crystal) then
            range_cell = supercell_circum(unit_cell, i_shell*shell_thickness)
        else
            range_cell = (/ 0, 0, 0 /)
        endif
        minr_cell = -1*range_cell
        idx_cell = (/ 0, 0, -1 /)
        do i_cell = 1, product(1+2*range_cell)
            call shift_cell(idx_cell, minr_cell, range_cell)
            ! MPI code begin
            if (is_parallel .and. is_crystal) then
                if (my_task /= modulo(i_cell, n_tasks)) cycle
            endif
            ! MPI code end
            if (is_crystal) then
                R_cell = matmul(idx_cell, unit_cell)
            else
                R_cell = (/ 0.d0, 0.d0, 0.d0 /)
            endif
            do i_atom = 1, size(xyz, 1)
                ! MPI code begin
                if (is_parallel .and. .not. is_crystal) then
                    if (my_task /= modulo(i_atom, n_tasks)) cycle
                endif
                ! MPI code end
                !! loop over pairwise interaction terms
                do j_atom = 1, i_atom
                    if ( (i_cell == 1) .and. (i_atom == j_atom) ) cycle
                    r = xyz(i_atom, :)-xyz(j_atom, :)-R_cell
                    r_norm = sqrt(sum(r**2))
                    if (r_norm > param_ts_cutoff_radius) cycle
                    if (r_norm >= i_shell*shell_thickness &
                        .or. r_norm < (i_shell-1)*shell_thickness) then
                        cycle
                    endif
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha_0(i_atom), alpha_0(j_atom))
!                    if (present(R_vdw)) then
                    R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
!                    endif
                    select case (version)
                        case ("fermi")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)
                        case ("fermi2")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)
                            f_damp = f_damp*f_damp
                        case ("erf")
                            f_damp = damping_erf(r_norm, s_R*R_vdw_ij, d)
                        case ("1mexp")
                            f_damp = damping_1mexp(r_norm, s_R*R_vdw_ij, d)
                    endselect
                    ene_pair = -C6_ij*f_damp/r_norm**6
                    if (i_atom == j_atom) then
                        ene_shell = ene_shell+ene_pair/2
                    else
                        ene_shell = ene_shell+ene_pair
                    endif
                enddo ! j_atom
            enddo ! i_atom
        enddo ! i_cell
        ! MPI code begin
        if (is_parallel) then
            call sync_sum(ene_shell)
        endif
        ! MPI code end
        ene = ene+ene_shell
        if (.not. is_crystal) exit
        if (i_shell > 1 .and. abs(ene_shell) < param_ts_energy_accuracy) then
            call print_log("Periodic TS converged in " &
                //trim(tostr(i_shell))//" shells, " &
                //trim(tostr(i_shell*shell_thickness*bohr))//" angstroms")
            exit
        endif
    enddo ! i_shell
end function get_ts_energy_lowmem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     DIPOLE INTERACTION TENSOR (ScaLAPACK)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine add_dipole_matrix_s(mode, version, xyz, row2glob, col2glob, &
                               do_cmplx, do_real, alpha, R_vdw, beta,a,&
                               C6, unit_cell, k_point, relay, relay_c)
!! REDUCED CUSTOMIZABILITY FOR THE SAKE OF MEMORY REQUIREMENTS
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                C: periodic boundary conditions
    !!                R: use reciprocal space formalism
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! row2glob:  mapping local row index to global row index (atoms)
    !! col2glob:  mapping local column index to global column (atoms)
    !! do_cmplx:  controls whether complex result array should be used
    !! do_real:   controls whether real result array should be used
    !! alpha:     (rescaled) isotropic atomic polarizability
    !! R_vdw:     (rescaled) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole interaction
    !! a:         damping factor for dipole-dipole coupling
    !! C6:        (rescaled) C6 coefficients
    !! unit_cell: well, unit cell [a.u.]
    !! k_point:   current k-point
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT(S) !!!!!
    !! relay:   real dipole-dipole tensor (real-space formalism)
    !! relay_c: complex dipole-dipole tensor (reciprocal-space formalism)
    !!!!!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! Tpp:               atom-pairwise (3x3) interaction tensor
    !! R_cell:            radius of cell
    !! r:                 position iAtom - positions jAtom
    !! r_norm:            interatomic distance
    !! R_vdw_ij:          combined van der Waals radius
    !! C6_ij:             combined C6 coefficient
    !! sigma_ij:          combined width of gaussian dipole density
    !! volume:            volume of unit cell
    !! ewald_alpha:       real/reciprocal parameter (ewald summation)
    !! real_space_cutoff: real-space cutoff radius [a.u.]
    !! Tpp_c:             complex atom-pairwise (3x3) interaction tensor
    !! i/j_local:         loop index for local matrix entries
    !! i/j_atom:          (most often corresponding) index for atom
    !! i_cell:            loop index over periodic images
    !! idx_cell:          index of current periodic image
    !! range_cell:        range of periodic images included
    !! i, j, i_xyz:       additional loop indices
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :)
    integer, intent(in) :: row2glob(:), col2glob(:), do_cmplx, do_real
    real(8), intent(in), optional :: alpha(size(xyz, 1))
    real(8), intent(in), optional :: R_vdw(size(xyz, 1)), beta, a
    real(8), intent(in), optional :: C6(size(xyz, 1)), unit_cell(3, 3)
    real(8), intent(in), optional :: k_point(3)
    real(8), intent(inout), optional :: &
          relay(do_real*3*size(row2glob), do_real*3*size(col2glob))
    complex(8), intent(inout), optional :: &
          relay_c(do_cmplx*3*size(row2glob), do_cmplx*3*size(col2glob))

    real(8) :: Tpp(3, 3), R_cell(3), r(3), r_norm, R_vdw_ij, C6_ij
    real(8) :: sigma_ij, volume, ewald_alpha, real_space_cutoff
    complex(8) :: Tpp_c(3, 3)
    integer :: i_local, j_local, i_atom, j_atom, i_cell, idx_cell(3),&
                range_cell(3), i, j
    logical :: is_crystal, is_reciprocal, is_low_dim, mute, do_ewald


    is_reciprocal = is_in('R', mode)
    is_crystal = is_in('C', mode) .or. is_reciprocal
    is_low_dim = any(param_vacuum_axis)
    do_ewald = .false.
    mute = is_in('M', mode)

    if (is_crystal) then
        if (is_low_dim) then
            real_space_cutoff = param_dipole_low_dim_cutoff
        else if (param_ewald_on) then
            do_ewald = .true.
            volume = max(abs(dble(product(diagonalized(unit_cell)))), 0.2d0)
            ewald_alpha = 2.5d0/(volume)**(1.d0/3)
            real_space_cutoff = 6.d0/ewald_alpha*param_ewald_real_cutoff_scaling
            call print_log('Ewald: using alpha = '//trim(tostr(ewald_alpha)) &
                //', real cutoff = '//trim(tostr(real_space_cutoff)), mute)
        else
            real_space_cutoff = param_dipole_cutoff
        endif
        range_cell = supercell_circum(unit_cell, real_space_cutoff)
        call print_log('Ewald: summing real part in cell vector range of ' &
            //trim(tostr(1+2*range_cell(1)))//'x' &
            //trim(tostr(1+2*range_cell(2)))//'x' &
            //trim(tostr(1+2*range_cell(3))), mute)
    else
        range_cell = (/ 0, 0, 0 /)
    endif
    call ts(11)
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_cell)
        call shift_cell(idx_cell, -range_cell, range_cell)
        if (is_crystal) then
            R_cell = matmul(idx_cell, unit_cell)
        else
            R_cell(:) = 0.d0
        endif
        do i_local = 1, size(row2glob)
            i_atom = row2glob(i_local)
            do j_local = 1, size(col2glob)
                j_atom = col2glob(j_local)
                if ( (i_cell==1) .and. (i_atom==j_atom) ) cycle

                r = xyz(i_atom, :)-xyz(j_atom, :)-R_cell
                r_norm = norm2(r)
                if (is_crystal .and. r_norm > real_space_cutoff) cycle
!                if (present(R_vdw)) then
                R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
!                endif
!                if (present(alpha)) then
                sigma_ij = norm2( get_sigma_selfint( &
                                  alpha((/ i_atom, j_atom /)) ) )
!                endif
                if (present(C6)) then
                    C6_ij = combine_C6( C6(i_atom), C6(j_atom), &
                                        alpha(i_atom), alpha(j_atom) )
                endif
                select case (version)
                    case ("fermi,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi,dip,gg")
                        Tpp = (1.d0-damping_fermi(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("bare")
                        Tpp = T_bare(r)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,erf")
                        Tpp = T_erf_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,fermi")
                        Tpp = T_fermi_coulomb(r, beta*R_vdw_ij, a)
                    case ("1mexp,dip")
                        Tpp = damping_1mexp(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("erf,dip")
                        Tpp = damping_erf(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi^2,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)**2*T_bare(r)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("1mexp,dip,gg")
                        Tpp = (1.d0-damping_1mexp(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                    case ("erf,dip,gg")
                        Tpp = (1.d0-damping_erf(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                        do_ewald = .false.
                endselect
                if (do_ewald) then
                    Tpp = Tpp+T_erfc(r, ewald_alpha)-T_bare(r)
                endif
                if (is_reciprocal) then
                    Tpp_c = Tpp*exp(-cmplx(0.d0, 1.d0, 8)*( &
                        dot_product(k_point, r)))
                endif
                i = 3*i_local-2
                j = 3*j_local-2
                if (is_reciprocal) then
                    relay_c(i:i+2, j:j+2) = relay_c(i:i+2, j:j+2) &
                        + Tpp_c
                else
                    relay(i:i+2, j:j+2) = relay(i:i+2, j:j+2) &
                        + Tpp
                endif
            enddo ! j_local
        enddo ! i_local
    enddo ! i_cell
    call ts(-11)
    if (do_ewald) then
        call add_ewald_dipole_parts_s(mode, xyz, unit_cell,ewald_alpha,&
                                 row2glob, col2glob, do_cmplx, do_real,&
                                 k_point, relay, relay_c)
    endif

end subroutine add_dipole_matrix_s


subroutine add_ewald_dipole_parts_s(mode, xyz, unit_cell, alpha, &
                          row2glob, col2glob, do_cmplx, do_real, &
                          k_point, relay, relay_c)
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                C: periodic boundary conditions
    !!                R: use reciprocal space formalism
    !! xyz:       atomic positions [a.u.]
    !! unit_cell: well, unit cell [a.u.]
    !! alpha:     Ewald parameter
    !! row2glob:  mapping local row index to global row index (atoms)
    !! col2glob:  mapping local column index to global column (atoms)
    !! do_cmplx:  controls whether complex result array should be used
    !! do_real:   controls whether real result array should be used
    !! k_point:   current k-point
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT(S) !!!!!
    !! relay:   real dipole-dipole tensor (real-space formalism)
    !! relay_c: complex dipole-dipole tensor (reciprocal-space formalism)
    !!!!!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! rec_unit_cell:     unit_cell in reciprocal space
    !! volume:            volume of unit cell
    !! G_vector:          current lattice (G) vector
    !! r:                 position iAtom - positions jAtom
    !! real_space_cutoff: real-space cutoff radius [a.u.]
    !! Tpp:               real atom-pairwise (3x3) interaction tensor
    !! Tpp_c:             complex atom-pairwise (3x3) interaction tensor
    !! i/j_local:         loop index for local matrix entries
    !! i/j_atom:          (most often corresponding) index for atom
    !! idx_G_vector:      indices for current lattice (G) vector
    !! range_G_vector:    range of lattice (G) vectors included
    !! i, j, i/j_xyz:     additional loop indices
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode
    real(8), intent(in) :: xyz(:, :), unit_cell(3, 3), &
                                      alpha
    integer, intent(in) :: row2glob(:), col2glob(:), &
                                      do_cmplx, do_real
    real(8), intent(in), optional :: k_point(3)

    real(8), intent(inout), optional :: &
          relay(do_real*3*size(row2glob), do_real*3*size(col2glob))
    complex(8), intent(inout), optional :: &
          relay_c(do_cmplx*3*size(row2glob), do_cmplx*3*size(col2glob))

    logical :: is_reciprocal, mute, do_surface
    real(8) :: rec_unit_cell(3, 3), volume, G_vector(3), r(3), &
                   k_total(3), k_sq, rec_space_cutoff, Tpp(3, 3), &
                   k_prefactor(3, 3), elem
    complex(8) :: Tpp_c(3, 3)
    integer :: i_local, j_local, i_atom, j_atom, i, j, i_xyz, &
                   j_xyz, idx_G_vector(3), i_G_vector, range_G_vector(3)


    is_reciprocal = is_in('R', mode)
    mute = is_in('M', mode)

    rec_unit_cell = 2*pi*inverted(transpose(unit_cell))
    volume = abs(dble(product(diagonalized(unit_cell))))
    rec_space_cutoff = 10.d0*alpha*param_ewald_rec_cutoff_scaling
    range_G_vector = supercell_circum(rec_unit_cell, rec_space_cutoff)
    call print_log('Ewald: using reciprocal cutoff = ' &
        //trim(tostr(rec_space_cutoff)), mute)
    call print_log('Ewald: summing reciprocal part in G vector range of ' &
        //trim(tostr(1+2*range_G_vector(1)))//'x' &
        //trim(tostr(1+2*range_G_vector(2)))//'x' &
        //trim(tostr(1+2*range_G_vector(3))), mute)
    call ts(12)
    idx_G_vector = (/ 0, 0, -1 /)
    do i_G_vector = 1, product(1+2*range_G_vector)
        call shift_cell(idx_G_vector, -range_G_vector, range_G_vector)
        if (i_G_vector == 1) cycle
        G_vector = matmul(idx_G_vector, rec_unit_cell)
        if (is_reciprocal) then
            k_total = k_point+G_vector
        else
            k_total = G_vector
        endif
        k_sq = sum(k_total*k_total)
        if (sqrt(k_sq) > rec_space_cutoff) cycle
        k_prefactor(:, :) = 4*pi/volume*exp(-k_sq/(4*alpha*alpha))
        forall (i_xyz = 1:3, j_xyz = 1:3) &
                k_prefactor(i_xyz, j_xyz) = k_prefactor(i_xyz, j_xyz) &
                *k_total(i_xyz)*k_total(j_xyz)/k_sq
        do i_local = 1, size(row2glob)
            i_atom = row2glob(i_local)
            do j_local = 1, size(col2glob)
                j_atom = col2glob(j_local)
                r = xyz(i_atom, :)-xyz(j_atom, :)
                if (is_reciprocal) then
                    Tpp_c = k_prefactor*exp(cmplx(0.d0, 1.d0, 8) &
                        *dot_product(G_vector, r))
                else
                    Tpp = k_prefactor*cos(dot_product(G_vector, r))
                endif
                i = 3*i_local-2
                j = 3*j_local-2
                if (is_reciprocal) then
                    relay_c(i:i+2, j:j+2) = relay_c(i:i+2, j:j+2)+Tpp_c
                else
                    relay(i:i+2, j:j+2) = relay(i:i+2, j:j+2) + Tpp
                endif
            enddo ! j_atom
        enddo ! i_atom
    enddo ! i_G_vector
    do i_local = 1, size(row2glob)  ! self energy
        i_atom = row2glob(i_local)
        do j_local = 1, size(col2glob)
            j_atom = col2glob(j_local)
            if (i_atom == j_atom) then
                do i_xyz = 1, 3
                    i = 3*(i_local - 1) + i_xyz
                    j = 3*(j_local - 1) + i_xyz
                    if (is_reciprocal) then
                        relay_c(i, j) = relay_c(i, j) - &
                                        4*alpha**3/(3*sqrt(pi))
                    else
                        relay(i, j) = relay(i, j) - 4*alpha**3/(3*sqrt(pi))
                    endif
                enddo
                !! there is no (iAtom==jAtom) beyond this point
                cycle
            endif
        enddo
    enddo
    do_surface = .true.
    if (present(k_point)) then
        k_sq = sum(k_point*k_point)
        if (sqrt(k_sq) > 1.d-15) then
            do_surface = .false.
            do i_local = 1, size(row2glob)
                do j_local = 1, size(col2glob)
                    do i_xyz = 1, 3
                        do j_xyz = 1, 3
                            i = 3*(i_local-1)+i_xyz
                            j = 3*(j_local-1)+j_xyz
                            elem = 4*pi/volume*k_point(i_xyz)*&
                                   k_point(j_xyz)/k_sq &
                                   *exp(-k_sq/(4*alpha*alpha))
                            if (is_reciprocal) then
                                relay_c(i, j) = relay_c(i, j) + elem
                            else
                                relay(i, j) = relay(i, j) + elem
                            endif ! is_reciprocal
                        enddo ! j_xyz
                    enddo ! i_xyz
                enddo ! j_local
            enddo ! i_local
        endif ! k_sq > threshold
    endif ! k_point present
    if (do_surface) then ! surface energy
        do i_local = 1, size(row2glob)
            do j_local = 1, size(col2glob)
                do i_xyz = 1, 3
                    i = 3*(i_local-1)+i_xyz
                    j = 3*(j_local-1)+i_xyz
                    if (is_reciprocal) then
                        relay_c(i, j) = relay_c(i, j) + 4*pi/(3*volume)
                    else
                        relay(i, j) = relay(i, j) + 4*pi/(3*volume)
                    endif
                enddo ! i_xyz
            enddo ! j_local
        enddo ! i_local
    endif
    call ts(-12)
end subroutine add_ewald_dipole_parts_s





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     SELF-CONSISTENT SCREENING (ScaLAPACK)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function run_scs_s(mode, version, xyz, alpha, R_vdw, beta, a,unit_cell)&
        result(alpha_scs)
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options (here: useless?)
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! alpha:     (rescaled) frequency-dependent isotropic atomic polarizabilities
    !! R_vdw:     (rescaled) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole coupling
    !! a:         damping factor for dipole-dipole coupling?
    !! unit_cell: well, unit cell [a.u]
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT !!!!!
    !! alpha_scs: frequency-dependent, isotropic atomic polarizabilities
    !!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! i_grid_omega: loop index for frequencies
    !! mute:         print messages during calculation?
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: &
        xyz(:, :), &
        alpha(:, :)
    real(8), intent(in), optional :: &
        R_vdw(size(xyz, 1)), &
        beta, a, &
        unit_cell(3, 3)
    real(8), dimension(size(alpha, 1), size(alpha, 2)) :: alpha_scs

    integer :: i_grid_omega
    character(len=1) :: mute


    if (is_in('M', mode)) then
        mute = 'M'
    else
        mute = ''
    endif

    alpha_scs(:, :) = 0.d0
    !! for each frequency: run screening
    do i_grid_omega = 0, n_grid_omega
        alpha_scs(i_grid_omega+1, :) = &
            do_scs_s(mode//mute, &
                     version, &
                     xyz, &
                     alpha(i_grid_omega+1,:), &
                     R_vdw=R_vdw, &
                     beta=beta, &
                     a=a, &
                     unit_cell=unit_cell)

        mute = 'M'
    enddo

end function run_scs_s


function do_scs_s(mode, version, xyz, alpha, R_vdw, beta, a, &
                  lam, unit_cell) result(alpha_scs_i)
    !!!!! INPUT !!!!!
    !! mode:          controls output (useless in this subroutine?)
    !! version:       which damping to use for dipole-dipole coupling
    !! xyz:           atomic positions [a.u]
    !! alpha:         (rescaled) isotropic atomic polarizabilities
    !! R_vdw:         (rescaled) van der Waals radii
    !! beta:          range-separation parameter
    !! a:             damping factor for dipole-dipole interaction
    !! lam:           scaling factor of bare dipole-dipole tensor
    !!!!!!!!!!!!!!!!!!

    !!!!! RESULT !!!!!
    !! alpha_scs_i:   screened, isotropic atomic polarizabilities
    !!!!!!!!!!!!!!!!!!

    !!!!! PRIVATE VARIABLES !!!!!
    !! n_tasks_bak:   back-up for number of MPI-tasks (= number of cores)
    !! size1n:        number of atoms
    !! size3n:        dimension of interaction tensor (= 3*size1n)
    !! my_task_blacs: processor on which each process is
    !! n_tasks_blacs: number of available processors for ScaLAPACK
    !! scs_ctxt:      BLACS context
    !! nprows:        number of rows on processor grid
    !! npcols:        number of columns on processor grid
    !! my_prow:       row on processor grid on which current task is
    !! my_pcol:       column on processor grid on which current task is
    !! fb3n:          block size in which (3Nx3N) matrix is distributed
    !! fb1n:          block size in which (Nx1) matrix is distributed
    !! my_nrows3n:    number of rows of local subarray of the (3Nx3N) array
    !! my_nrows1n:    number of rows of local subarray of the (Nx1) array
    !! my_ncols3n:    number of columns of local subarray of the (2Nx3N) array
    !! my_ncols1n:    number of rows of local subarray of the (Nx1) array
    !! numroc:        external routine to determine size of local subarrays
    !! desc3n3n:      BLACS descriptor for (3Nx3N) matrix
    !! desc1n:        BLACS descriptor for (Nx1) matrix
    !! my_row2glob:   mapping local atom index to global atom index (rows)
    !! my_row2glob:   mapping local atom index to global atom index (columns)
    !! i/j_local:     loop index over local subarrays
    !! i/j_atom:      (most often the corresponding) atom index
    !! i/j_xyz:       loop index over (x,y,z) components
    !! i/j:           additional (auxiliary) loop indices
    !! ierr:          error status
    !! my_scs_lwork:  size of work array for inversion
    !! my_scs_liwork: size of (integer) work array for inversion
    !! my_alphafull:  local subarray of dipole-dipole tensor
    !! my_scs_work:   work array for inversion
    !! my_ipiv:       pivoting array of local subarray for L/U decomposition
    !! my_scs_iwork:  (integer) work array for inversion
    !! my_alphascs:   subarray of screened atomic polarizability
    !! alpha_3_3:     (3x3) screened anisotropic atomic polarizability
    !! alpha_diag:    eigenvalues of alpha_3_3 (mean = isotropic alpha)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha(:)
    real(8), intent(in), optional :: R_vdw(size(xyz, 1))
    real(8), intent(in), optional :: beta, a, unit_cell(3, 3), lam
    real(8) :: alpha_scs_i(size(xyz, 1))
    logical :: scale_lambda
    integer :: n_tasks_bak, size1n, size3n, my_task_blacs, n_tasks_blacs
    integer :: scs_ctxt, nprows, npcols, my_prow, my_pcol, fb3n, fb1n
    integer :: my_nrows3n, my_nrows1n, my_ncols3n, my_ncols1n, numroc
    integer, dimension(9) :: desc3n3n, desc1n
    integer, dimension(:), allocatable :: my_row2glob, my_col2glob
    integer :: i_local, i_atom, i, i_xyz, j_local, j_atom, j, j_xyz
    integer :: ierr, my_scs_lwork, my_scs_liwork
    real(8), dimension(:,:), allocatable :: my_alphafull
    real(8), dimension(:), allocatable :: my_scs_work
    integer, dimension(:), allocatable :: my_ipiv, my_scs_iwork
    real(8), dimension(:,:), allocatable :: my_alphascs
    real(8) :: alpha_3_3(3,3), alpha_diag(3)


    !! should I scale bare dipole-dipole tensor?
    scale_lambda = is_in('L', mode)

    !!!!!  BACKUP NTASKS, DEFINE GLOBAL ARRAY SIZES  !!!!!
    n_tasks_bak = n_tasks
    n_tasks = 1
    size1n = size(xyz,1)
    size3n = 3*size1n


    !!!!!  BLACS INITIALIZATION  !!!!!
    !! get current processor and total number of processors
    call BLACS_PINFO(my_task_blacs, n_tasks_blacs)

    !! build square processor grid from available processors
    do nprows = int(  sqrt( dble(n_tasks_blacs) )  ), 1, -1
        if (mod(n_tasks_blacs, nprows) == 0) exit
    enddo
    npcols = n_tasks_blacs/nprows

!    call BLACS_GET(0, 0, scs_ctxt)
    scs_ctxt = mpi_comm_global
    call BLACS_GRIDINIT(scs_ctxt, 'R', nprows, npcols)
    call BLACS_GRIDINFO(scs_ctxt, nprows, npcols, my_prow, my_pcol)

    !! block sizes, local shapes for (3n by 3n) and (n by n) matrices
    !! fb3n can be adjusted for performance (assert: mod(fb3n, 3) = 0)
    fb3n = 3
    fb1n = fb3n/3
    my_nrows3n = max(1, numroc(size3n, fb3n, my_prow, 0, nprows))
    my_nrows1n = max(1, numroc(size1n, fb1n, my_prow, 0, nprows))
    my_ncols3n = max(1, numroc(size3n, fb3n, my_pcol, 0, npcols))
    my_ncols1n = max(1, numroc(size1n, fb1n, my_pcol, 0, npcols))

    !! BLACS descriptors
    call descinit(desc3n3n, size3n, size3n, fb3n, fb3n, 0,0,scs_ctxt, &
                  my_nrows3n, ierr)
    call descinit(desc1n, size1n, 1, fb1n, 1, 0, 0, scs_ctxt, &
                  my_nrows1n, ierr)


    !!!!!  GET ATOM MAPPINGS  !!!!!
    if ( allocated(my_row2glob) ) deallocate(my_row2glob)
    allocate( my_row2glob(my_nrows1n), stat=ierr)
    if ( allocated(my_col2glob) ) deallocate(my_col2glob)
    allocate( my_col2glob(my_ncols1n), stat=ierr)
    !! build my_row(col)2glob and my_diagelems (TODO)
    do i_local = 1, my_nrows1n
         my_row2glob(i_local) = map_local2global(i_local, nprows, &
                                                 fb1n, my_prow, 0)
    enddo
    do i_local = 1, my_ncols1n
         my_col2glob(i_local) = map_local2global(i_local, npcols, &
                                                 fb1n, my_pcol, 0)
    enddo


    !!!!!  BUILD (LOCAL) DIPOLE-DIPOLE (SUB)TENSOR  !!!!!
    if ( allocated(my_alphafull) ) deallocate(my_alphafull)
    allocate( my_alphafull(my_nrows3n, my_ncols3n), stat=ierr)
    my_alphafull(:, :) = 0.d0

    call add_dipole_matrix_s(blanked('R', mode), &
                             version, &
                             xyz=xyz, &
                             row2glob=my_row2glob, &
                             col2glob=my_col2glob, &
                             do_cmplx=0, do_real=1, &
                             alpha=alpha, &
                             R_vdw=R_vdw, &
                             beta=beta, &
                             a=a, &
                             unit_cell=unit_cell, &
                             relay=my_alphafull)

    if (scale_lambda) then
        my_alphafull = lam*my_alphafull
    endif
    !! update diagonal terms (this will change to loop over my_diagelems)
    do i_local = 1, my_nrows1n
        i_atom = my_row2glob(i_local)
        do j_local = 1, my_ncols1n
            j_atom = my_col2glob(j_local)
            if (i_atom == j_atom) then
                do i_xyz = 1, 3
                    i = 3*(i_local-1) + i_xyz
                    j = 3*(j_local-1) + i_xyz
                    my_alphafull(i, j) = my_alphafull(i, j) + &
                                         1.d0/alpha(i_atom)
                enddo ! i_xyz
                !! there is no (i_atom==j_atom) beyond this point
                cycle
            endif
        enddo ! j_local
    enddo ! j_local


    !!!!!  INVERT (LOCAL) DIPOLE-DIPOLE (SUB)TENSOR  !!!!!
    call ts(31)

    !! step 1: L/U decomposition
    if ( allocated(my_ipiv) ) deallocate(my_ipiv)
    allocate( my_ipiv(size3n) , stat=ierr)
    call PDGETRF(size3n, size3n, my_alphafull, 1, 1, desc3n3n, &
                 my_ipiv, ierr)

    !! step 2a: work space query for inversion
    if ( allocated(my_scs_work) ) deallocate(my_scs_work)
    allocate(my_scs_work(1), stat=ierr)
    if ( allocated(my_scs_iwork) ) deallocate(my_scs_iwork)
    allocate(my_scs_iwork(1), stat=ierr)
    call PDGETRI(size3n, my_alphafull, 1, 1, desc3n3n, my_ipiv, &
                 my_scs_work, -1, my_scs_iwork, -1, ierr)

    !! step 2b: actual inversion
    my_scs_lwork = nint(my_scs_work(1))
    my_scs_liwork = int(my_scs_iwork(1))
    if ( allocated(my_scs_work) ) deallocate(my_scs_work)
    allocate(my_scs_work(my_scs_lwork), stat=ierr)
    if ( allocated(my_scs_iwork) ) deallocate(my_scs_iwork)
    allocate(my_scs_iwork(my_scs_liwork), stat=ierr)
    call PDGETRI(size3n, my_alphafull, 1, 1, desc3n3n, my_ipiv, &
                my_scs_work, my_scs_lwork, my_scs_iwork, my_scs_liwork,&
                ierr)


    !!!!!  CONTRACT ATOMIC POLARIZABILITIES  !!!!!
    !!
    !! (diagonalized) SCS-tensor is composed of NxN atomic 3x3 submatrices.
    !! Thus, in terms of SCS-tensor elements ij(AtomA,AtomB):
    !!
    !!  / xx(1,1) xy(1,1) xz(1,1)   xx(1,2) xy(1,2) xz(1,2)   xx(1,3) ... \
    !!  | yx(1,1) yy(1,1) yz(1,1)   yx(1,2) yy(1,2) yz(1,2)   yx(1,3) ... |
    !!  | zx(1,1) zy(1,1) zz(1,1)   zx(1,2) zy(1,2) zz(1,2)   zx(1,3) ... |
    !!  |                                                                 |
    !!  | xx(2,1) xy(2,1) xz(2,1)   xx(2,2) xy(2,2) xz(2,2)   xx(2,3) ... |
    !!  |   ...     ...     ...       ...     ...     ...       ...   ... |
    !!  \ xx(N,1) xy(N,1) xz(N,1)   xx(N,2) xy(N,2) xz(N,2)   xx(N,3) ... /
    !!
    !! step 1: ij element of the (3x3) screened anisotropic polarizability
    !! for atom iAtom is given by sum of all elements ij(iAtom,jAtom) ,
    !! where ij \in {xx, xy, xz, yx, yy, yz, zx, zy, zz}
    !!
    !! step 2: contracted isotropic polarizability is given by the sum of
    !! eigenvalues of the (3x3) anisotropic atomic polarizability tensor.
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( allocated(my_alphascs) ) deallocate(my_alphascs)
    allocate( my_alphascs(my_nrows1n,1), stat=ierr)
    my_alphascs(:,1) = 0.d0

    do i_local = 1, my_nrows1n
        alpha_3_3(:,:) = 0.d0
        !! sum over contributions of all atoms (columns) in LOCAL (sub)array
        do i_xyz = 1, 3
            do j_xyz = 1, 3
                i = 3*(i_local-1)+i_xyz
                alpha_3_3(i_xyz, j_xyz) = &
                             sum( my_alphafull(i, j_xyz:my_ncols3n:3) )
            enddo
        enddo

        !! -> each process has one (3x3) matrix, which corresponds to the
        !! anisotropic polarizability of atom iAtom including screening
        !! by all atoms contained in LOCAL submatrix
        !! next: sum all anisotropic alpha_3_3 and communicate
        !! VY: Could be rewritten with MPI_Allreduce
        call DGSUM2D(scs_ctxt, 'R', ' ', 3, 3, alpha_3_3, 3, my_prow, -1)

        !! now on all columns of the processor grid:
        !! anisotropic alpha of atom i_atom including screening by ALL atoms
        !! next: diagonalize anisotropic alpha to get isotropic alpha
        call diagonalize_sym_dble_('N', alpha_3_3, alpha_diag)
        my_alphascs(i_local,1) = sum(alpha_diag)/3.d0
    enddo


    !! gather isotropic atomic polarizabilities into results array
    do i_atom=1, size1n
        call pdelget('A',' ', alpha_scs_i(i_atom), my_alphascs, &
                     i_atom, 1, desc1n)
    enddo

    !!!!!  DEALLOCATIONS  !!!!!
    if ( allocated(my_alphafull) ) deallocate(my_alphafull)
    if ( allocated(my_row2glob) ) deallocate(my_row2glob)
    if ( allocated(my_col2glob) ) deallocate(my_col2glob)
    if ( allocated(my_scs_work) ) deallocate(my_scs_work)
    if ( allocated(my_scs_iwork) ) deallocate(my_scs_iwork)
    if ( allocated(my_ipiv) ) deallocate(my_ipiv)
    if ( allocated(my_alphascs) ) deallocate(my_alphascs)

    !!!!!  RESTORE NTASKS AND FINISH UP  !!!!!
    n_tasks = n_tasks_bak
    call ts(-31)
    call BLACS_GRIDEXIT(scs_ctxt)

end function do_scs_s





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     MBD ENERGY (ScaLAPACK)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_mbd_energy_s(mode, version, xyz, alpha_0, omega, supercell, &
                          k_grid, unit_cell, R_vdw, beta, a, C6)&
                          result(ene)
!! REDUCED CUSTOMIZABILITY FOR THE SAKE OF MEMORY REQUIREMENTS
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                C: use periodic boundary conditions
    !!                R: use reciprocal space formulation
    !!                E: output eigenvalues to file
    !!                V: output eigenvectors to file
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! alpha_0:   (rescaled and screened) isotropic atomic polarizabilities
    !! omega:     effective frequencies/force constants of QHOs
    !! supercell: size of supercell to use (if applicable)
    !! k_grid:    k-point grid to use (if applicable)
    !! unit_cell: well, unit cell [a.u.]
    !! R_vdw:     (rescaled and screened) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole interaction
    !! a:         damping factor for dipole-dipole interaction
    !! C6:        (rescaled and screened) C6 coefficients
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT !!!!!
    !! ene : MBD energy [a.u.]
    !!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! is_X:    flags wheter to do/use X (see mode)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha_0(size(xyz, 1)),&
                                      omega(size(xyz, 1))
    integer, intent(in), optional :: supercell(3)
    real(8), intent(in), optional :: k_grid(:, :), unit_cell(3, 3), &
                                      R_vdw(size(xyz, 1)), beta, a, &
                                      C6(size(xyz, 1))

    real(8) :: ene

    logical :: is_reciprocal, is_crystal, is_parallel


    is_parallel = is_in('P', mode)
    is_crystal = is_in('C', mode)
    is_reciprocal = is_in('R', mode)
!    if (.not. is_in('M', mode)) then
!        call print_log('---------------------------------------------------&
!                    &------------')
!        call print_log('Framework:   ScaLAPACK')
!        selectcase(trim(eigensolver))
!            case('qr')
!                call print_log('eigensolver: QR algorithm')
!            case('dandc')
!                call print_log('eigensolver: Divide and Conquer')
!            case('mrrr')
!                call print_log('eigensolver: Multiple Relatively Robust &
!                            &Representations (MRRR)')
!        endselect
!        call print_log('---------------------------------------------------&
!                    &------------')
!    end if

    if (.not. is_crystal) then
        ene = get_single_mbd_energy_s(blanked('R',mode), version, xyz, &
                               alpha_0, omega, R_vdw=R_vdw, beta=beta, &
                               a=a, C6=C6, unit_cell=unit_cell)
    else
        if (is_reciprocal) then
            ene = get_reciprocal_mbd_energy_s(mode, version, xyz, alpha_0,&
                                   omega, k_grid, unit_cell, R_vdw=R_vdw, &
                                   beta=beta, a=a, C6=C6)
        else
            ene = get_supercell_mbd_energy_s(mode, version, xyz, alpha_0,&
                                             omega, unit_cell, supercell,&
                                             R_vdw=R_vdw, beta=beta, a=a,&
                                             C6=C6)
        endif
    endif
end function get_mbd_energy_s


function get_supercell_mbd_energy_s(mode, version, xyz, alpha_0, omega, &
                  unit_cell, supercell, R_vdw, beta, a, C6) result(ene)
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                C: use periodic boundary conditions
    !!                E: output eigenvalues to file
    !!                   (TODO: reduce output to one unit cell)
    !!                V: output eigenvectors to file
    !!                   (TODO: reduce output to one unit cell)
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! alpha_0:   (rescaled and screened) isotropic atomic polarizabilities
    !! omega:     effective frequencies/force constants of QHOs
    !! unit_cell: well, unit cell [a.u.]
    !! supercell: size of supercell to use (if applicable)
    !! R_vdw:     (rescaled and screened) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole interaction
    !! a:         damping factor for dipole-dipole interaction
    !! C6:        (rescaled and screened) C6 coefficients
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT !!!!!
    !! ene : MBD energy [a.u.]
    !!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! R_cell:  radii of current supercell [a.u.]
    !! i_atom:  loop index over atoms
    !! i:       loop index
    !! is_X:    flags wheter to do/use X (see mode)
    !! alpha:   frequenct-dependent isotropic atomic polarizabilities (for RPA)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha_0(size(xyz, 1)),&
                                      omega(size(xyz, 1)), unit_cell(3, 3)
    integer, intent(in) :: supercell(3)
    real(8), intent(in), optional :: R_vdw(size(xyz, 1)), beta, a, &
                                      C6(size(xyz, 1))

    real(8) :: ene

    real(8) :: R_cell(3)
    integer :: i_atom, i, i_cell, idx_cell(3), n_cells

    real(8), allocatable :: xyz_super(:, :), alpha_0_super(:), &
                             omega_super(:), R_vdw_super(:), C6_super(:)
    real(8) :: unit_cell_super(3, 3)


    n_cells = product(supercell)
    !! build supercell
    do i = 1, 3
        unit_cell_super(i, :) = unit_cell(i, :)*supercell(i)
    enddo
    allocate (xyz_super(n_cells*size(xyz, 1), 3))
    allocate (alpha_0_super(n_cells*size(alpha_0)))
    allocate (omega_super(n_cells*size(omega)))
    allocate (R_vdw_super(n_cells*size(R_vdw)))
    allocate (C6_super(n_cells*size(C6)))
    idx_cell = (/ 0, 0, -1 /)
    do i_cell = 1, n_cells
        call shift_cell(idx_cell, (/ 0, 0, 0 /), supercell-1)
        R_cell = matmul(idx_cell, unit_cell)
        do i_atom = 1, size(xyz, 1)
            i = (i_cell-1)*size(xyz, 1)+i_atom
            xyz_super(i, :) = xyz(i_atom, :)+R_cell
            alpha_0_super(i) = alpha_0(i_atom)
            omega_super(i) = omega(i_atom)
            if (present(R_vdw)) then
                R_vdw_super(i) = R_vdw(i_atom)
            endif
            if (present(C6)) then
                C6_super(i) = C6(i_atom)
            endif
        enddo
    enddo
    !! get MBD energy for supercell
    ene = get_single_mbd_energy_s(mode, &
                                  version, &
                                  xyz_super, &
                                  alpha_0_super, &
                                  omega_super, &
                                  R_vdw=R_vdw_super, &
                                  beta=beta, &
                                  a=a, &
                                  C6=C6_super, &
                                  unit_cell=unit_cell_super)

    deallocate (xyz_super)
    deallocate (alpha_0_super)
    deallocate (omega_super)
    deallocate (R_vdw_super)
    deallocate (C6_super)

    !! correct energy down to one unit cell
    ene = ene/n_cells

end function get_supercell_mbd_energy_s


function get_single_mbd_energy_s(mode, version, xyz, alpha_0, omega, &
                            R_vdw, beta, a, C6, unit_cell) result(ene)
!! REDUCED CUSTOMIZABILITY FOR THE SAKE OF MEMORY REQUIREMENTS
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                C: periodic boundary conditions
    !!                R: use reciprocal space formalism
    !!                E: output eigenvalues to file
    !!                V: output eigenvectors to file
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! alpha_0:   (rescaled and screened) isotropic atomic polarizabilities
    !! omega:     effective frequencies/force constants of QHOs
    !! R_vdw:     (rescaled) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole interaction
    !! a:         damping factor for dipole-dipole coupling
    !! C6:        (rescaled) C6 coefficients
    !! unit_cell: well, unit cell [a.u.]
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT(S) !!!!!
    !! ene: MBD energy [a.u.]
    !!!!!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! my_relay:       local subarray of interaction tensor
    !! my_evecs:       local subarray of diagonalized tensor (eigenvectors)
    !! my_cfdm_work:   work array for diagonalization
    !! my_row2glob:    mapping local row index to global index (row)
    !! my_col2glob:    mapping local column index to global index (column)
    !! n_tasks_bak:    back-up for number of MPI-tasks (= number of cores)
    !! size1n:         number of atoms
    !! size3n:         dimension of interaction tensor (= 3*size1n)
    !! my_task_blacs:  processor on which each process is
    !! n_tasks_blacs:  number of available processors for ScaLAPACK
    !! cfdm_ctxt:      BLACS context
    !! nprows:         number of rows on processor grid
    !! npcols:         number of columns on processor grid
    !! my_prow:        row on processor grid on which current task is
    !! my_pcol:        column on processor grid on which current task is
    !! fb3n:           block size in which (3Nx3N) matrix is distributed
    !! fb1n:           block size in which (Nx1) matrix is distributed
    !! my_nrows3n:     number of rows of local subarray of the (3Nx3N) array
    !! my_nrows1n:     number of rows of local subarray of the (Nx1) array
    !! my_ncols3n:     number of columns of local subarray of the (2Nx3N) array
    !! my_ncols1n:     number of rows of local subarray of the (Nx1) array
    !! numroc:         external routine to determine size of local subarrays
    !! desc3n3n:       BLACS descriptor for (3Nx3N) matrix
    !! desc1n:         BLACS descriptor for (Nx1) matrix
    !! my_cfdm_lwork:  size of work array for inversion
    !! eigs:           eigenvectors of interaction tensor
    !! i/j_local:      loop index over local subarrays
    !! i/j_atom:       (most often the corresponding) atom index
    !! i_xyz:          loop index over (x,y,z) components
    !! i/j:            additional (auxiliary) loop indices
    !! ierr:           error status
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha_0(size(xyz, 1)),&
                                      omega(size(xyz, 1))
    real(8), intent(in), optional :: R_vdw(size(xyz, 1)), beta, a, &
                                      C6(size(xyz, 1)), unit_cell(3, 3)

    real(8) :: ene

    real(8), dimension(:,:), allocatable :: my_relay, my_evecs
    real(8), dimension(:), allocatable :: my_cfdm_work, my_write_work
    integer, dimension(:), allocatable :: my_row2glob, my_col2glob
    integer, dimension(:), allocatable :: my_cfdm_iwork
    integer :: vl, vu, il, iu, nvals, nvecs
    integer :: n_tasks_bak, size1n, size3n, my_task_blacs
    integer :: cfdm_ctxt, nprows, npcols, my_prow, my_pcol, fb3n, fb1n
    integer :: my_nrows3n, my_nrows1n, my_ncols3n, my_ncols1n, numroc
    integer, dimension(9) :: desc3n3n, desc1n
    integer :: my_cfdm_lwork, my_cfdm_liwork, n_tasks_blacs
    integer :: nev

    real(8) :: eigs(3*size(xyz, 1))
    integer :: i_local, j_local, i_atom, j_atom, i_xyz, i, j
    integer :: n_negative_eigs, ierr
    logical :: get_eigenvalues, get_eigenvectors, io_proc

    get_eigenvalues = is_in('E', mode)
    get_eigenvectors = is_in('V', mode)

    n_tasks_bak = n_tasks
    n_tasks = 1

    size1n = size(xyz,1)
    size3n = 3*size1n

    !! blacs context, grid info, ...
    call BLACS_PINFO(my_task_blacs, n_tasks_blacs)
    do nprows = int(  sqrt( dble(n_tasks_blacs) )  ), 1, -1
        if (mod(n_tasks_blacs, nprows) == 0) exit
    enddo
    npcols = n_tasks_blacs/nprows
!    call BLACS_GET(0, 0, cfdm_ctxt)
    cfdm_ctxt = mpi_comm_global
    call BLACS_GRIDINIT(cfdm_ctxt, 'R', nprows, npcols)
    call BLACS_GRIDINFO(cfdm_ctxt, nprows, npcols, my_prow, my_pcol)
    io_proc = ( (my_prow==0) .and. (my_pcol==0) )

    !! block sizes, local shapes for (3n by 3n) and (n by n) matrices
    !! fb3n can be adjusted for performance (assert: mod(fb3n, 3) = 0)
    fb3n = 3

    fb1n = fb3n/3
    my_nrows3n = max(1, numroc(size3n, fb3n, my_prow, 0, nprows))
    my_nrows1n = max(1, numroc(size1n, fb1n, my_prow, 0, nprows))
    my_ncols3n = max(1, numroc(size3n, fb3n, my_pcol, 0, npcols))
    my_ncols1n = max(1, numroc(size1n, fb1n, my_pcol, 0, npcols))

    !! blacs descriptors
    call descinit(desc3n3n, size3n, size3n, fb3n, fb3n,0,0,cfdm_ctxt,&
                  my_nrows3n, ierr)

    call descinit(desc1n, size1n, 1, fb1n, 1, 0, 0, cfdm_ctxt, &
                  my_nrows1n, ierr)

    if ( allocated(my_row2glob) ) deallocate(my_row2glob)
    allocate( my_row2glob(my_nrows1n), stat=ierr)
    if ( allocated(my_col2glob) ) deallocate(my_col2glob)
    allocate( my_col2glob(my_ncols1n), stat=ierr)

    !! build my_row(col)2glob (and my_diagidx TODO)
    do i_atom=1, my_nrows1n
         my_row2glob(i_atom) = map_local2global(i_atom, nprows, fb1n, &
                                                my_prow, 0)
    enddo
    do i_atom=1, my_ncols1n
         my_col2glob(i_atom) = map_local2global(i_atom, npcols, fb1n, &
                                                my_pcol, 0)
    enddo
    if ( allocated(my_relay) ) deallocate(my_relay)
    allocate( my_relay(my_nrows3n, my_ncols3n), stat=ierr)

    my_relay(:, :) = 0.d0
    call add_dipole_matrix_s(mode, &
                             version, &
                             xyz, &
                             my_row2glob, &
                             my_col2glob, &
                             do_cmplx=0, do_real=1, &
                             alpha=alpha_0, &
                             R_vdw=R_vdw, &
                             beta=beta, &
                             a=a, &
                             C6=C6, &
                             unit_cell=unit_cell, &
                             relay=my_relay)

    do i_local = 1, size(my_row2glob)
        i_atom = my_row2glob(i_local)
        do j_local = 1, size(my_col2glob)
            j_atom = my_col2glob(j_local)
            i = 3*i_local-2
            j = 3*j_local-2
            my_relay(i:i+2, j:j+2) = & ! my_relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                *sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                my_relay(i:i+2, j:j+2)

            if (i_atom == j_atom) then
                do i_xyz = 1, 3
                    i = 3*i_local-3 + i_xyz
                    j = 3*j_local-3 + i_xyz
                    my_relay(i, j) = my_relay(i, j) + &
                                     omega(i_atom)*omega(i_atom)
                    ! V(i,i) = w^2+sqrt(a*a)*w*w*T
                enddo
            endif
        enddo
    enddo
    call ts(21)
    if ( allocated(my_cfdm_work) ) deallocate(my_cfdm_work)
    allocate(my_cfdm_work(1))

    if(allocated(my_evecs)) deallocate(my_evecs)
    allocate(my_evecs(my_nrows3n,my_ncols3n),stat=ierr)

    if(get_eigenvectors) then
       nev = size3n
    else
       nev = 1
    endif

    call aims_elsi_stdevp(size3n, my_nrows3n, my_ncols3n, nev, my_relay, &
        eigs, my_evecs, cfdm_ctxt, fb3n, mpi_comm_global)

    !! deallocations
    deallocate(my_evecs)
    deallocate(my_cfdm_work, stat=ierr)
    deallocate(my_row2glob)
    deallocate(my_col2glob)
    deallocate(my_relay, stat=ierr)
    call ts(-21)

    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call print_warning( &
            "CFDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues" &
        )
        if (param_zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    endif

    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)

    if (get_eigenvalues) then
        where (eigs < 0) eigs = 0.d0
        eigs = sqrt(eigs)
        if (io_proc) call write_eigenenergies(eigs, "mbd_eigenenergies.out")
    endif

    ! destroy blacs grid, restore n_tasks
    call BLACS_GRIDEXIT(cfdm_ctxt)
    n_tasks = n_tasks_bak

end function get_single_mbd_energy_s


function get_reciprocal_mbd_energy_s(mode, version, xyz,alpha_0, &
                   omega, k_grid, unit_cell, R_vdw, beta, a, C6) &
                   result(ene)
!! REDUCED CUSTOMIZABILITY FOR THE SAKE OF MEMORY REQUIREMENTS
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                E: output eigenvalues to file
    !!                V: output eigenvectors to file (TODO: DEBUG)
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! alpha_0:   (rescaled and screened) isotropic atomic polarizabilities
    !! omega:     effective frequencies/force constants of QHOs
    !! k_grid:    k-point grid to use
    !! unit_cell: well, unit cell [a.u.]
    !! R_vdw:     (rescaled and screened) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole interaction
    !! a:         damping factor for dipole-dipole coupling
    !! C6:        (rescaled and screened) C6 coefficients
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT !!!!!
    !! ene:       MBD energy [a.u.]
    !!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! i_kpt:     loop index for k-point grid
    !! k_point:   current k-point
    !! mute:      print messages during calculation?
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha_0(size(xyz, 1))
    real(8), intent(in) :: omega(size(xyz, 1)), &
                                      k_grid(:, :), unit_cell(3, 3)
    real(8), intent(in), optional :: R_vdw(size(xyz, 1)), beta, a
    real(8), intent(in), optional :: C6(size(xyz, 1))

    real(8) :: ene

    integer :: i_kpt
    real(8) :: k_point(3)
    character(len=1) :: mute


    if (is_in('M', mode)) then
        mute = 'M'
    else
        mute = ''
    endif

    ene = 0.d0
    do i_kpt = 1, size(k_grid, 1)
        k_point = k_grid(i_kpt, :)
        ene = ene+get_single_reciprocal_mbd_ene_s( &
                         blanked('P', mode)//mute, &
                         version, &
                         xyz, &
                         alpha_0, &
                         omega, &
                         k_point, &
                         unit_cell, &
                         i_kpt, &
                         R_vdw=R_vdw, &
                         beta=beta, &
                         a=a, &
                         C6=C6)

        mute = 'M'
    enddo ! k_point loop
    ene = ene/size(k_grid, 1)

end function get_reciprocal_mbd_energy_s


function get_single_reciprocal_mbd_ene_s(mode, version, xyz, alpha_0, &
                omega, k_point, unit_cell, i_kpt, R_vdw, beta, a, C6) &
                result(ene)
!! REDUCED CUSTOMIZABILITY FOR THE SAKE OF MEMORY REQUIREMENTS
    !!!!! INPUT !!!!!
    !! mode:      controls output and runtime options
    !!                C: periodic boundary conditions
    !!                R: use reciprocal space formalism
    !!                E: output eigenvalues to file
    !!                V: return eigenvectors to file (TODO: check output)
    !! version:   which damping to use for dipole-dipole coupling
    !! xyz:       atomic positions [a.u.]
    !! alpha_0:   (rescaled and screened) isotropic atomic polarizabilities
    !! omega:     effective frequencies/force constants of QHOs
    !! k_point:   current k-point
    !! unit_cell: well, unit cell [a.u.]
    !! i_kpt:     index of current k-point
    !! R_vdw:     (rescaled) van der Waals radii [a.u.]
    !! beta:      range-separation parameter for dipole-dipole interaction
    !! a:         damping factor for dipole-dipole coupling
    !! C6:        (rescaled) C6 coefficients
    !!!!!!!!!!!!!!!!!

    !!!!! RESULT(S) !!!!!
    !! ene: MBD energy [a.u.]
    !!!!!!!!!!!!!!!!!!!!!

    !!!!! LOCAL VARIABLES !!!!!
    !! my_relay:       local subarray of interaction tensor
    !! my_evecs:       local subarray of diagonalized tensor (eigenvectors)
    !! my_cfdm_work:   work array for diagonalization
    !! my_cfdm_rwork:  real work array for diagonalization
    !! my_row2glob:    mapping local row index to global index (row)
    !! my_col2glob:    mapping local column index to global index (column)
    !! n_tasks_bak:    back-up for number of MPI-tasks (= number of cores)
    !! size1n:         number of atoms
    !! size3n:         dimension of interaction tensor (= 3*size1n)
    !! my_task_blacs:  processor on which each process is
    !! n_tasks_blacs:  number of available processors for ScaLAPACK
    !! cfdm_ctxt:      BLACS context
    !! nprows:         number of rows on processor grid
    !! npcols:         number of columns on processor grid
    !! my_prow:        row on processor grid on which current task is
    !! my_pcol:        column on processor grid on which current task is
    !! fb3n:           block size in which (3Nx3N) matrix is distributed
    !! fb1n:           block size in which (Nx1) matrix is distributed
    !! my_nrows3n:     number of rows of local subarray of the (3Nx3N) array
    !! my_nrows1n:     number of rows of local subarray of the (Nx1) array
    !! my_ncols3n:     number of columns of local subarray of the (2Nx3N) array
    !! my_ncols1n:     number of rows of local subarray of the (Nx1) array
    !! numroc:         external routine to determine size of local subarrays
    !! desc3n3n:       BLACS descriptor for (3Nx3N) matrix
    !! desc1n:         BLACS descriptor for (Nx1) matrix
    !! my_cfdm_lwork:  size of work array for inversion
    !! my_cfdm_lrwork: size of (real) work array for inversion
    !! eigs:           eigenvectors of interaction tensor
    !! i/j_local:      loop index over local subarrays
    !! i/j_atom:       (most often the corresponding) atom index
    !! i_xyz:          loop index over (x,y,z) components
    !! i/j:            additional (auxiliary) loop indices
    !! ierr:           error status
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=*), intent(in) :: mode, version
    real(8), intent(in) :: xyz(:, :), alpha_0(size(xyz, 1)),&
                                      omega(size(xyz, 1)), k_point(3), &
                                      unit_cell(3, 3)
    integer, intent(in), optional :: i_kpt
    real(8), intent(in), optional :: R_vdw(size(xyz, 1)), beta, a, &
                                      C6(size(xyz, 1))

    real(8) :: ene

    complex(8), dimension(:,:), allocatable :: my_relay, my_evecs
    complex(8), dimension(:), allocatable :: my_cfdm_work
    complex(8), dimension(:), allocatable :: my_write_work
    real(8), dimension(:), allocatable :: my_cfdm_rwork
    integer, dimension(:), allocatable :: my_row2glob, my_col2glob
    integer, dimension(:), allocatable :: my_cfdm_iwork
    integer :: vl, vu, il, iu, nvals, nvecs
    integer :: n_tasks_bak, size1n, size3n, my_task_blacs, n_tasks_blacs
    integer :: cfdm_ctxt, nprows, npcols, my_prow, my_pcol, fb3n, fb1n
    integer :: my_nrows3n, my_nrows1n, my_ncols3n, my_ncols1n, numroc
    integer, dimension(9) :: desc3n3n, desc1n
    integer :: nev

    integer :: my_cfdm_lwork, my_cfdm_lrwork, my_cfdm_liwork

    real(8) :: eigs(3*size(xyz, 1))
    integer :: i_local, j_local, i_atom, j_atom, i_xyz, i, j
    integer :: n_negative_eigs, ierr, kptID, nrows_max, ncols_max
    logical :: get_eigenvalues, get_eigenvectors, io_proc, &
                output_data, fexists

    get_eigenvalues = is_in('E', mode)
    get_eigenvectors = is_in('V', mode)
    output_data = (get_eigenvalues .or. get_eigenvectors)
    if ( (output_data) .and. (.not. present(i_kpt)) ) then
        call print_warning("You chose to output eigenvectors but did"//&
                          &" not specify the current k-point index."//&
                          &"Using random index for output file!")

        kptID = get_FileID()
    else
        kptID = i_kpt
    endif

    n_tasks_bak = n_tasks
    n_tasks = 1

    size1n = size(xyz,1)
    size3n = 3*size1n

    !! blacs context, grid info, ...
    call BLACS_PINFO(my_task_blacs, n_tasks_blacs)
    do nprows = int(  sqrt( dble(n_tasks_blacs) )  ), 1, -1
        if (mod(n_tasks_blacs, nprows) == 0) exit
    enddo
    npcols = n_tasks_blacs/nprows
!    call BLACS_GET(0, 0, cfdm_ctxt)
    cfdm_ctxt = mpi_comm_global
    call BLACS_GRIDINIT(cfdm_ctxt, 'R', nprows, npcols)
    call BLACS_GRIDINFO(cfdm_ctxt, nprows, npcols, my_prow, my_pcol)
    io_proc = ( (my_prow==0) .and. (my_pcol==0) )

    !! block sizes, local shapes for (3n by 3n) and (n by n) matrices
    !! fb3n can be adjusted for performance (assert: mod(fb3n, 3) = 0)
    fb3n = 3

    fb1n = fb3n/3
    my_nrows3n = max(1, numroc(size3n, fb3n, my_prow, 0, nprows))
    my_nrows1n = max(1, numroc(size1n, fb1n, my_prow, 0, nprows))
    my_ncols3n = max(1, numroc(size3n, fb3n, my_pcol, 0, npcols))
    my_ncols1n = max(1, numroc(size1n, fb1n, my_pcol, 0, npcols))
    nrows_max = max(1, numroc(size3n, fb3n, 0, 0, nprows))
    ncols_max = max(1, numroc(size3n, fb3n, 0, 0, npcols))

    !! blacs descriptors
    call descinit(desc3n3n, size3n, size3n, fb3n, fb3n,0,0,cfdm_ctxt,&
                  my_nrows3n, ierr)

    call descinit(desc1n, size1n, 1, fb1n, 1, 0, 0, cfdm_ctxt, &
                  my_nrows1n, ierr)

    if ( allocated(my_row2glob) ) deallocate(my_row2glob)
    allocate( my_row2glob(my_nrows1n), stat=ierr)
    if ( allocated(my_col2glob) ) deallocate(my_col2glob)
    allocate( my_col2glob(my_ncols1n), stat=ierr)

    !! build my_row(col)2glob (and my_diagidx TODO)
    do i_atom=1, my_nrows1n
         my_row2glob(i_atom) = map_local2global(i_atom, nprows, fb1n, &
                                                my_prow, 0)
    enddo
    do i_atom=1, my_ncols1n
         my_col2glob(i_atom) = map_local2global(i_atom, npcols, fb1n, &
                                                my_pcol, 0)
    enddo
    if ( allocated(my_relay) ) deallocate(my_relay)
    allocate( my_relay(my_nrows3n, my_ncols3n), stat=ierr)

    my_relay(:, :) = cmplx(0.d0, 0.d0, 8)
    call add_dipole_matrix_s(mode, &
                             version, &
                             xyz, &
                             my_row2glob, &
                             my_col2glob, &
                             do_cmplx=1, do_real=0, &
                             alpha=alpha_0, &
                             R_vdw=R_vdw, &
                             beta=beta, &
                             a=a, &
                             C6=C6, &
                             unit_cell=unit_cell, &
                             k_point=k_point, &
                             relay_c=my_relay)

    do i_local = 1, size(my_row2glob)
        i_atom = my_row2glob(i_local)
        do j_local = 1, size(my_col2glob)
            j_atom = my_col2glob(j_local)
            i = 3*i_local-2
            j = 3*j_local-2
            my_relay(i:i+2, j:j+2) = & ! my_relay = sqrt(a*a)*w*w*T
                omega(i_atom)*omega(j_atom) &
                *sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                my_relay(i:i+2, j:j+2)

            if (i_atom == j_atom) then
                do i_xyz = 1, 3
                    i = 3*i_local-3 + i_xyz
                    j = 3*j_local-3 + i_xyz
                    my_relay(i, j) = my_relay(i, j) + &
                                     omega(i_atom)*omega(i_atom)
                    ! V(i,i) = w^2+sqrt(a*a)*w*w*T
                enddo
            endif
        enddo
    enddo
    call ts(22)
    if ( allocated(my_cfdm_work) ) deallocate(my_cfdm_work)
    allocate(my_cfdm_work(1))
    if ( allocated(my_cfdm_rwork) ) deallocate(my_cfdm_rwork)
    allocate(my_cfdm_rwork(1))

    if(allocated(my_evecs)) deallocate(my_evecs)
    allocate(my_evecs(my_nrows3n,my_ncols3n),stat=ierr)

    if(get_eigenvectors) then
       nev = size3n
    else
       nev = 1
    endif

    call aims_elsi_stdevp(size3n, my_nrows3n, my_ncols3n, nev, my_relay, &
        eigs, my_evecs, cfdm_ctxt, fb3n, mpi_comm_global)

    !! deallocations
    deallocate(my_evecs)
    deallocate(my_cfdm_work, stat=ierr)
    deallocate(my_cfdm_rwork, stat=ierr)
    deallocate(my_row2glob)
    deallocate(my_col2glob)
    deallocate(my_relay, stat=ierr)
    call ts(-22)

    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        call print_warning( &
            "CFDM Hamiltonian has " // trim(tostr(n_negative_eigs)) // &
            " negative eigenvalues" &
        )
        if (param_zero_negative_eigs) where (eigs < 0) eigs = 0.d0
    endif

    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)

    !! write eigenvalues to output file
    if (get_eigenvalues) then
        where (eigs < 0) eigs = 0.d0
        eigs = sqrt(eigs)
        if (io_proc) then
            !! remove old outputfile of eigenenergies
            if (kptID == 1) then
                inquire(file='mbd_eigenenergies_reciprocal.out', exist=fexists)
                if ( fexists ) then
                    i = get_FileID()
                    open(file="mbd_eigenenergies_reciprocal.out", unit=i, &
                         status='old')
                    close(i, status='delete')
                endif
            endif

            call write_eigenenergies(eigs, &
                                   "mbd_eigenenergies_reciprocal.out", &
                                   k_point=k_point)
        endif
    endif

    ! destroy blacs grid, restore n_tasks
    call BLACS_GRIDEXIT(cfdm_ctxt)
    n_tasks = n_tasks_bak

end function get_single_reciprocal_mbd_ene_s





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     SCS/MBD (ScaLAPACK) HELPER FUNCTIONS     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function map_local2global(idx_local, npX, blocksize, my_pX, pX_start) &
         result(idx_global)
        !!!!! INPUT !!!!!
        !! idx_local:  local row/column index
        !! npX:        number of processor rows/columns
        !! blocksize:  blocksize used to distribute rows/columns of matrix
        !! my_pX:      processor row/column of current task
        !! pX_start:   processor row/column where 1st row/column of matrix
        !!             is distributed
        !!!!!!!!!!!!!!!!!

        !!!!! RESULT !!!!!
        !! idx_global: corresponding row/column index of global matrix
        !!!!!!!!!!!!!!!!!!

        integer, intent(in) :: idx_local,npX,blocksize,my_pX,pX_start
        integer :: idx_global


        idx_global = npX*floor(real(idx_local-1)/blocksize)
        idx_global = idx_global + mod(my_pX-pX_start, npX)
        idx_global = blocksize*idx_global + mod(idx_local-1,blocksize)+1

end function map_local2global


subroutine write_eigenenergies(eigenes, evals_fname, k_point)
    real(8), intent(in) :: eigenes(:)
    character(len=*), intent(in) :: evals_fname

    real(8), intent(in), optional :: k_point(3)
    integer :: fid_evals, i_eval


    fid_evals = get_FileID()
    if ( present(k_point) )  then
        open(file=evals_fname, unit=fid_evals, position='append', &
             status='unknown', form='formatted')
    else
        open(file=evals_fname, unit=fid_evals, status='replace', &
             form='formatted')
    endif

    if (present(k_point)) write(fid_evals, "(A10,3F10.6)") "k point = ",&
                                      &k_point(1), k_point(2), k_point(3)
    do i_eval=1, size(eigenes)
        write(fid_evals, *) eigenes(i_eval)
    enddo
    if ( present(k_point) ) write(fid_evals, *) ""
    close(fid_evals)

end subroutine write_eigenenergies


function get_FileID()
    integer :: get_FileID
    integer, save :: currentID = minFileID

    if (currentID > maxFileID) then
        call print_log("error('get_FileID: No more free file identification numbers')")
    end if
    get_FileID = currentID
    currentID = currentID + 1

end function get_FileID


subroutine ZWRITEEIGVEC(eigenvectors, filenam, kpoint)
    complex(8), intent(in) :: eigenvectors(:,:)
    character(len=*),intent(in) :: filenam
    real(8), intent(in), optional :: kpoint(3)

    integer :: fID, n, i, j


    fID = get_FileId()
    n = size(eigenvectors,1)
    open( fID, file=filenam, status='unknown', form='unformatted' )
    write( fID ) n, n
    if ( present(kpoint) ) write( fID ) kpoint
    do i=1, n
        do j=1, n
            write( fID ) eigenvectors(j, i)
        enddo
    enddo
    close( fID )

end subroutine ZWRITEEIGVEC


subroutine DWRITEEIGVEC(eigenvectors, filenam)
    real(8), intent(in) :: eigenvectors(:,:)
    character(len=*), intent(in) :: filenam

    integer :: fID, n, i, j


    fID = get_FileId()
    n = size(eigenvectors,1)
    open( fID, file=filenam, status='unknown', form='unformatted' )
    write( fID ) n, n
    do i=1, n
        do j=1,n
            write( fID ) eigenvectors(j, i)
        enddo
    enddo
    close( fID )

end subroutine DWRITEEIGVEC


subroutine exit_blacs_and_finalize(continuation)
    integer, intent(in), optional :: continuation
    integer :: blacs_continue


    if ( present(continuation) ) then
        blacs_continue = continuation
    else
        blacs_continue = 0
    endif

    call blacs_exit(blacs_continue)

end subroutine exit_blacs_and_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     DIPOLE-DIPOLE COUPLING AND DAMPING     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function T_bare(rxyz) result(T)
    real(8), intent(in) :: rxyz(3)
    real(8) :: T(3, 3)

    integer :: i, j
    real(8) :: r_sq, r_5

    r_sq = sum(rxyz(:)**2)
    r_5 = sqrt(r_sq)**5
    do i = 1, 3
        T(i, i) = (3.d0*rxyz(i)**2-r_sq)/r_5
        do j = i+1, 3
            T(i, j) = 3.d0*rxyz(i)*rxyz(j)/r_5
            T(j, i) = T(i, j)
        end do
    end do
    T = -T
end function


real(8) function B_erfc(r, a) result(B)
    real(8), intent(in) :: r, a

    B = (erfc(a*r)+(2*a*r/sqrt(pi))*exp(-(a*r)**2))/r**3
end function


real(8) elemental function C_erfc(r, a) result(C)
    real(8), intent(in) :: r, a

    C = (3*erfc(a*r)+(2*a*r/sqrt(pi))*(3.d0+2*(a*r)**2)*exp(-(a*r)**2))/r**5
end function


function T_erfc(rxyz, alpha) result(T)
    real(8), intent(in) :: rxyz(3), alpha
    real(8) :: T(3, 3)

    integer :: i, j
    real(8) :: r, B, C

    r = sqrt(sum(rxyz(:)**2))
    B = B_erfc(r, alpha)
    C = C_erfc(r, alpha)
    do i = 1, 3
        do j = i, 3
            T(i, j) = -C*rxyz(i)*rxyz(j)
            if (i /= j) T(j, i) = T(i, j)
        end do
        T(i, i) = T(i, i)+B
    end do
end function


function damping_fermi(r, sigma, a) result(f)
    real(8), intent(in) :: r, sigma, a
    real(8) :: f

    f = 1.d0/(1+exp(-a*(r/sigma-1)))
end function


function damping_erf(r, sigma, a) result(f)
    real(8), intent(in) :: r, sigma, a
    real(8) :: f

    f = erf((r/sigma)**a)
end function


function damping_1mexp(r, sigma, a) result(f)
    real(8), intent(in) :: r, sigma, a
    real(8) :: f

    f = 1-exp(-(r/sigma)**a)
end function


function damping_overlap(r, overlap, C6, beta, a) result(f)
    real(8), intent(in) :: r, overlap, C6, beta, a
    real(8) :: f

    f = 1.d0-terf(-overlap/(erf(r/6)**6*C6/r**6), beta, a)
end function


function T_overlap_coulomb(rxyz, overlap, C6, beta, a) result(T)
    real(8), intent(in) :: rxyz(3), overlap, C6, beta, a
    real(8) :: T(3, 3)

    real(8) :: zeta_1, zeta_2
    real(8) :: r, erff, exp36, qene, qenep, qenepp

    r = sqrt(sum(rxyz**2))
    erff = erf(r/6)
    exp36 = exp(r**2/36)
    qene = overlap*r**6/(C6*erff**6)
    qenep = 2.d0*overlap*r**5*(-(1.d0/exp36)*r/sqrt(pi)+3.d0*erff)/(C6*erff**7)
    qenepp = (1.d0/exp36**2)*overlap*r**4/(9*C6*pi*erff**8) &
        *(42*r**2+exp36*sqrt(pi)*r*(-216+r**2)*erff &
        +270.d0*exp36**2*pi*erff**2)
    zeta_1 = 1.d0/2*(2.d0-erf(a*(beta-qene))+erf(a*(beta+qene)) &
        +2*a*r*qenep/sqrt(pi)*(-exp(-a**2*(beta-qene)**2) &
        -exp(-a**2*(beta+qene)**2)))
    zeta_2 = 1.d0/sqrt(pi)*a*exp(-a**2*(beta+qene)**2)*r**2 &
        *(2*a**2*qenep**2*(beta*(-1.d0+exp(4*a**2*beta*qene)) &
        -qene*(1.d0+exp(4*a**2*beta*qene))) &
        +qenepp*(1.d0+exp(4*a**2*beta*qene)))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function


function T_fermi_coulomb(rxyz, sigma, a) result(T)
    real(8), intent(in) :: rxyz(3), sigma, a
    real(8) :: T(3, 3)

    real(8) :: r_sigma, d_r_sigma_m_1, zeta_1, zeta_2

    r_sigma = sqrt(sum(rxyz**2))/sigma
    d_r_sigma_m_1 = a*(r_sigma-1)
    zeta_1 = 1.d0/(1.d0+exp(-d_r_sigma_m_1)) &
        -a/2.d0*r_sigma/(1.d0+cosh(-d_r_sigma_m_1))
    zeta_2 = 2.d0*a**2*r_sigma**2/sinh(-d_r_sigma_m_1)**3 &
        *sinh(-d_r_sigma_m_1/2.d0)**4
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function


function T_erf_coulomb(rxyz, sigma, a) result(T)
    real(8), intent(in) :: rxyz(3), sigma, a
    real(8) :: T(3, 3)

    real(8) :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = erf(r_sigma)-2.d0/sqrt(pi)*a*r_sigma*exp(-r_sigma**2)
    zeta_2 = -2.d0/sqrt(pi)*a*r_sigma*exp(-r_sigma**2) &
        *(1.d0+a*(-1.d0+2.d0*r_sigma**2))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function


function T_1mexp_coulomb(rxyz, sigma, a) result(T)
    real(8), intent(in) :: rxyz(3), sigma, a
    real(8) :: T(3, 3)

    real(8) :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = 1.d0-exp(-r_sigma)-a*r_sigma*exp(-r_sigma)
    zeta_2 = -r_sigma*a*exp(-r_sigma)*(1+a*(-1+r_sigma))
    T = zeta_1*T_bare(rxyz)-zeta_2*(rxyz .cprod. rxyz)/sqrt(sum(rxyz**2))**5
end function


subroutine get_damping_parameters(xc, ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, &
        mbd_ts_erf_beta, mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta)
    character(len=*), intent(in) :: xc
    real(8), intent(out) :: &
        ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta

    ts_d = 20.d0
    ts_s_r = 1.d0
    mbd_scs_a = 2.d0
    mbd_ts_a = 6.d0
    mbd_ts_erf_beta = 1.d0
    mbd_ts_fermi_beta = 1.d0
    mbd_rsscs_a = 6.d0
    mbd_rsscs_beta = 1.d0
    select case (xc)
        case ("pbe")
            ts_s_r = 0.94d0
            mbd_scs_a = 2.56d0
            mbd_ts_erf_beta = 1.07d0
            mbd_ts_fermi_beta = 0.81d0
            mbd_rsscs_beta = 0.83d0
        case ("pbe0")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("hse")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("blyp")
            ts_s_r = 0.62d0
        case ("b3lyp")
            ts_s_r = 0.84d0
        case ("revpbe")
            ts_s_r = 0.60d0
        case ("am05")
            ts_s_r = 0.84d0
    endselect
end subroutine get_damping_parameters


elemental function terf(r, r0, a)
    real(8), intent(in) :: r, r0, a
    real(8) :: terf

    terf = 0.5d0*(erf(a*(r+r0))+erf(a*(r-r0)))
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     GENERAL OPERATIONS/ROUTINES     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function alpha_dynamic_ts_all(mode, n, alpha_0, C6, omega) result(alpha_dyn)
    character(len=1), intent(in) :: mode
    integer, intent(in) :: n
    real(8), intent(in) :: alpha_0(:)
    real(8), intent(in), optional :: C6(size(alpha_0)), omega(size(alpha_0))
    real(8) :: alpha_dyn(0:n, size(alpha_0))

    integer :: i_grid_omega

    do i_grid_omega = 0, n_grid_omega
        alpha_dyn(i_grid_omega, :) = alpha_dynamic_ts( &
            mode, alpha_0, omega_grid(i_grid_omega), C6, omega)
    end do
end function alpha_dynamic_ts_all


function alpha_dynamic_ts(mode, alpha_0, u, C6, omega) result(alpha)
    character(len=1), intent(in) :: mode
    real(8), intent(in) :: alpha_0(:), u
    real(8), intent(in), optional :: C6(size(alpha_0)), omega(size(alpha_0))
    real(8) :: alpha(size(alpha_0))

    select case (mode)
        case  ('O')
            alpha(:) = alpha_osc(alpha_0, omega, u)
        case ('C')
            alpha(:) = alpha_osc(alpha_0, omega_eff(C6, alpha_0), u)
    end select
end function


elemental function alpha_osc(alpha_0, omega, u) result(alpha)
    real(8), intent(in) :: alpha_0, omega, u
    real(8) :: alpha

    alpha = alpha_0/(1+(u/omega)**2)
end function


elemental function combine_C6 (C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6_ij)
    real(8), intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    real(8) :: C6_ij

    C6_ij = 2*C6_i*C6_j/(alpha_0_j/alpha_0_i*C6_i+alpha_0_i/alpha_0_j*C6_j)
end function


elemental function V_to_R(V) result(R)
    real(8), intent(in) :: V
    real(8) :: R

    R = (3.d0*V/(4.d0*pi))**(1.d0/3)
end function


function omega_eff(C6, alpha) result(omega)
    real(8), intent(in) :: C6(:), alpha(size(C6))
    real(8) :: omega(size(C6))

    omega = 4.d0/3*C6/alpha**2
end function


elemental function get_sigma_selfint(alpha) result(sigma)
    real(8), intent(in) :: alpha
    real(8) :: sigma

    sigma = param_mayer_scaling*(sqrt(2.d0/pi)*alpha/3.d0)**(1.d0/3)
end function


function get_C6_from_alpha(alpha) result(C6)
    real(8), intent(in) :: alpha(:, :)
    real(8) :: C6(size(alpha, 2))
    integer :: i_atom

    do i_atom = 1, size(alpha, 2)
        C6(i_atom) = 3.d0/pi*sum((alpha(:, i_atom)**2)*omega_grid_w(:))
    end do
end function


function get_total_C6_from_alpha(alpha) result(C6)
    real(8), intent(in) :: alpha(:, :)
    real(8) :: C6

    C6 = 3.d0/pi*sum((sum(alpha, 2)**2)*omega_grid_w(:))
end function


function supercell_circum(uc, radius) result(sc)
    real(8), intent(in) :: uc(3, 3), radius
    integer :: sc(3)

    real(8) :: ruc(3, 3), layer_sep(3)
    integer :: i

    ruc = 2*pi*inverted(transpose(uc))
    forall (i = 1:3) layer_sep(i) = sum(uc(i, :)*ruc(i, :)/sqrt(sum(ruc(i, :)**2)))
    sc = ceiling(radius/layer_sep+0.5d0)
    where (param_vacuum_axis) sc = 0
end function


subroutine shift_cell(ijk, first_cell, last_cell)
    integer, intent(inout) :: ijk(3)
    integer, intent(in) :: first_cell(3), last_cell(3)

    integer :: i_dim, i

    do i_dim = 3, 1, -1
        i = ijk(i_dim)+1
        if (i <= last_cell(i_dim)) then
            ijk(i_dim) = i
            return
        else
            ijk(i_dim) = first_cell(i_dim)
        end if
    end do
end subroutine


function make_g_grid(n1, n2, n3) result(g_grid)
    integer, intent(in) :: n1, n2, n3
    real(8) :: g_grid(n1*n2*n3, 3)

    integer :: g_kpt(3), i_kpt, kpt_range(3)
    real(8) :: g_kpt_shifted(3)

    g_kpt = (/ 0, 0, -1 /)
    kpt_range = (/ n1, n2, n3 /)
    do i_kpt = 1, n1*n2*n3
        call shift_cell(g_kpt, (/ 0, 0, 0 /), kpt_range-1)
        g_kpt_shifted = dble(g_kpt)+param_k_grid_shift
        where (2*g_kpt_shifted > kpt_range)
            g_kpt_shifted = g_kpt_shifted-dble(kpt_range)
        end where
        g_grid(i_kpt, :) = g_kpt_shifted/kpt_range
    end do
end function make_g_grid


function make_k_grid(g_grid, uc) result(k_grid)
    real(8), intent(in) :: g_grid(:, :), uc(3, 3)
    real(8) :: k_grid(size(g_grid, 1), 3)

    integer :: i_kpt
    real(8) :: ruc(3, 3)

    ruc = 2*pi*inverted(transpose(uc))
    do i_kpt = 1, size(g_grid, 1)
        k_grid(i_kpt, :) = matmul(g_grid(i_kpt, :), ruc)
    end do
end function make_k_grid


function cart_prod_(a, b) result(c)
    real(8), intent(in) :: a(:), b(:)
    real(8) :: c(size(a), size(b))

    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i)*b(j)
        end do
    end do
end function


function get_diag_(A) result(d)
    real(8), intent(in) :: A(:, :)
    real(8) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function


function get_diag_cmplx_(A) result(d)
    complex(8), intent(in) :: A(:, :)
    complex(8) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function


function make_diag_(d) result(A)
    real(8), intent(in) :: d(:)
    real(8) :: A(size(d), size(d))

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:size(d)) A(i, i) = d(i)
end function


function solve_lin_sys(A, b) result(x)
    real(8), intent(in) :: A(:, :), b(size(A, 1))
    real(8) :: x(size(b))

    real(8) :: A_(size(b), size(b))
    integer :: i_pivot(size(b))
    integer :: n
    integer :: error_flag

    A_ = A
    x = b
    n = size(b)
    call DGESV(n, 1, A_, n, i_pivot, x, n, error_flag)
end function


function eye(n) result(A)
    integer, intent(in) :: n
    real(8) :: A(n, n)

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:n) A(i, i) = 1.d0
end function


character(len=50) elemental function tostr_int_(k, format)
    implicit none

    integer, intent(in) :: k
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_int_, format) k
    else
        write (tostr_int_, "(i20)") k
    end if
    tostr_int_ = adjustl(tostr_int_)
end function tostr_int_


character(len=50) elemental function tostr_dble_(x, format)
    implicit none

    double precision, intent(in) :: x
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_dble_, format) x
    else
        write (tostr_dble_, "(g50.17e3)") x
    end if
    tostr_dble_ = adjustl(tostr_dble_)
end function tostr_dble_

!!!!!!!!!!!!!!!!!!!!
!!     TIMING     !!
!!!!!!!!!!!!!!!!!!!!

subroutine ts(id, always)
    integer, intent(in) :: id
    logical, intent(in), optional :: always

    if (measure_time .or. present(always)) then
        call system_clock(ts_cnt, ts_rate, ts_cnt_max)
        if (id > 0) then
            timestamps(id) = timestamps(id)-ts_cnt
        else
            ts_aid = abs(id)
            timestamps(ts_aid) = timestamps(ts_aid)+ts_cnt
            ts_counts(ts_aid) = ts_counts(ts_aid)+1
        end if
    end if
end subroutine ts


function clock_rate() result(rate)
    integer :: cnt, rate, cnt_max

    call system_clock(cnt, rate, cnt_max)
end function clock_rate


end module mbd_dev
