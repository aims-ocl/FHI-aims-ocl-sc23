module hirshfeld
! This module should contain all extensions of the older hirshfeld_analysis
! subroutine. --JH, 2016-01-20

use geometry, only: coords, lattice_vector, species, empty
use species_data, only: species_name
use dimensions, only: &
    n_atoms, n_full_points, n_max_radial, n_species, n_max_spline, n_periodic
use constants, only: pi, hartree
use mpi_tasks, only: myid, aims_stop
use synchronize_mpi_basic, only: sync_vector
use dimensions, only: n_my_batches, n_spin
use localorb_io, only: localorb_info
use grids, only: &
    batches, grid_point, r_grid_min, r_grid_inc, invert_log_grid, &
    n_grid, r_radial, w_radial, w_angular, n_radial, scale_radial, &
    invert_radial_grid, get_radial_weight
use spline, only: val_spline
use free_atoms, only: free_rho_spl, free_drho_dr_spl
use pbc_lists, only: &
    centers_basis_integrals, coords_center, species_center, center_to_atom
use xml_write, only: xml_open_file, Xml_file_t, xml_elem, xml_open, &
    xml_close, xml_close_file, tostr
use sym_base, only: evaluate_densmat_sym
use mbd_geom, only: get_freq_grid

implicit none

private

logical, public :: grid_out = .false.
logical, public :: use_pulay = .false.
integer, public :: grid_atoms(100) = -1
real(8), parameter, public :: C_vv = 0.0093d0  ! see evaluate_free_atom_quantities()

integer, parameter :: n_freq = 15

type(Xml_file_t) :: xml_file_grid

public :: &
    contract_with_dVdn, run_hirshfeld, &
    evaluate_free_atom_quantities, numerical_hirshfeld_deriv, &
    eval_hirsh_vol_pulay_deriv_dm, numerical_hirsh_pulay_forces, &
    hirshfeld_initialize, hirshfeld_finalize

real(8), allocatable :: &
    hirshfeld_pulay_forces(:, :, :), hirshfeld_pulay_stress(:, :, :)

contains

subroutine hirshfeld_initialize()
    if (use_pulay) then
        allocate (hirshfeld_pulay_forces(n_atoms, n_atoms, 3))
        if (n_periodic > 0) allocate (hirshfeld_pulay_stress(n_atoms, 3, 3))
    end if
end subroutine

subroutine hirshfeld_finalize()
    if (allocated(hirshfeld_pulay_forces)) deallocate (hirshfeld_pulay_forces)
    if (allocated(hirshfeld_pulay_stress)) deallocate (hirshfeld_pulay_stress)
end subroutine

subroutine run_hirshfeld(rho_grid, volumes, dvoldr, dvoldh, shift, strain, &
        rho_grad_grid, kinetic_density_grid, alpha_0_vv_nm, C6_vv_nm)
    real(8), intent(in) :: rho_grid(n_spin, n_full_points)
    real(8), intent(out) :: volumes(n_atoms)
    real(8), intent(out), optional :: dvoldr(n_atoms, n_atoms, 3)
    real(8), intent(out), optional :: dvoldh(n_atoms, 3, 3)
    real(8), intent(in), optional :: shift(n_atoms, 3)
    real(8), intent(in), optional :: strain(3, 3)
    real(8), intent(in), optional :: rho_grad_grid(3, n_spin, n_full_points)
    real(8), intent(in), optional :: kinetic_density_grid(:, :)
    real(8), intent(out), optional :: alpha_0_vv_nm(:)
    real(8), intent(out), optional :: C6_vv_nm(:)

    real(8) :: rho, rho_grad_norm, part_weight, weight, dweightdr(n_atoms, 3), &
        dweightdh(3, 3), h_inv(3, 3), kin_dens, ion_pot, alpha
    integer :: i_full_point, i_batch, i_batch_point, i_atom, i, i_freq
    type(grid_point) :: point
    real(8), allocatable :: vv_pols(:, :), vv_pol(:), vv_pols_nm(:, :)
    logical, allocatable :: print_grid_atom(:)
    real(8), allocatable :: omega_grid(:), omega_grid_w(:)

    if (grid_out) then
        call xml_open_file("grid-"//trim(tostr(myid))//".xml", "grid", xml_file_grid)
        allocate (print_grid_atom(n_atoms), source=.false.)
        if (grid_atoms(1) == -1) print_grid_atom(:) = .true.
        do i = 1, size(grid_atoms)
            if (grid_atoms(i) == -1) exit
            print_grid_atom(grid_atoms(i)) = .true.
        end do
    end if
    volumes(:) = 0.d0
    if (present(dvoldr)) dvoldr(:, :, :) = 0.d0
    if (present(dvoldh)) dvoldh(:, :, :) = 0.d0
    if (present(alpha_0_vv_nm) .or. present(C6_vv_nm)) then
        allocate (omega_grid(0:n_freq), source=0.d0)
        allocate (omega_grid_w(0:n_freq), source=0.d0)
        call get_freq_grid(n_freq, omega_grid(1:n_freq), omega_grid_w(1:n_freq))
        allocate (vv_pol(0:n_freq))
        allocate (vv_pols(0:n_freq, n_atoms), source=0d0)
        allocate (vv_pols_nm(0:n_freq, n_atoms), source=0d0)
    end if
    if (n_periodic > 0) h_inv = inverted3x3(lattice_vector)
    i_full_point = 0
    do i_batch = 1, n_my_batches
        do i_batch_point = 1, batches(i_batch)%size
            i_full_point = i_full_point+1
            point = batches(i_batch)%points(i_batch_point)
            i_atom = point%index_atom
            if (empty(i_atom)) cycle
            if (grid_out) xml_file_grid%muted = .not. print_grid_atom(i_atom)
            call evaluate_hirshfeld_factor( &
                'W', point, weight, part_weight, dweightdr, dweightdh, shift, strain &
            )
            ! also opens the xml elem if weight > 0
            if (.not. weight > 0.d0) cycle
            rho = sum(rho_grid(:, i_full_point))
            call xml_elem('i_atom', i_atom, file=xml_file_grid)
            call xml_elem('rho', rho, file=xml_file_grid)
            volumes(i_atom) = volumes(i_atom)+rho*weight
            if (present(alpha_0_vv_nm) .or. present(C6_vv_nm)) then
                rho_grad_norm = sqrt(sum(sum(rho_grad_grid(:, :, i_full_point), 2)**2))
                call xml_elem('rho_grad_norm', rho_grad_norm, file=xml_file_grid)
                kin_dens = sum(kinetic_density_grid(:, i_full_point))
                call xml_elem('kin_dens', kin_dens, file=xml_file_grid)
                ion_pot = (rho_grad_norm/rho)**2/8
                alpha = (kin_dens/2 - rho_grad_norm**2/(8*rho)) / &
                    (3d0/10*(3*pi**2)**(2d0/3)*rho**(5d0/3))
                vv_pol = vv_pol_func(rho, rho_grad_norm, omega_grid, C_vv)
                call xml_elem('vv_pol', vv_pol(0), file=xml_file_grid)
                if (rho > 0) then
                    vv_pols(:, i_atom) = vv_pols(:, i_atom) + vv_pol*part_weight
                    vv_pols_nm(:, i_atom) = vv_pols_nm(:, i_atom) + &
                        vv_pol*nm_cutoff(ion_pot, alpha)*part_weight
                end if
            end if
            if (present(dvoldr)) then
                dvoldr(i_atom, :, :) = dvoldr(i_atom, :, :)+rho*dweightdr
                call xml_elem('rho_dweightdr', rho*dweightdr, file=xml_file_grid)
                if (present(dvoldh)) then
                    dvoldh(i_atom, :, :) = dvoldh(i_atom, :, :)+rho*dweightdh
                    call xml_elem('rho_dweightdh', rho*dweightdh, file=xml_file_grid)
                end if
            end if
            call xml_close(file=xml_file_grid)
        end do
    end do
    xml_file_grid%muted = .false.
    call sync_vector(volumes, size(volumes))
    if (.not. (present(shift) .or. present(strain))) then
        call xml_elem('volumes', volumes)
    end if
    if (present(alpha_0_vv_nm) .or. present(C6_vv_nm)) then
        call sync_vector(vv_pols, size(vv_pols))
        call xml_elem('vv_pols', vv_pols)
        call sync_vector(vv_pols_nm, size(vv_pols_nm))
        call xml_elem('vv_pols_nm', vv_pols_nm)
        call xml_elem('omega_grid', omega_grid)
        call xml_elem('omega_grid_w', omega_grid_w)
        if (present(alpha_0_vv_nm)) then
            alpha_0_vv_nm = vv_pols_nm(0, :)
            call xml_elem('alpha_0_vv_nm', alpha_0_vv_nm)
        end if
        if (present(C6_vv_nm)) then
            C6_vv_nm(:) = 0.d0
            do i_freq = 1, n_freq
                C6_vv_nm = C6_vv_nm + 3d0/pi*vv_pols_nm(i_freq, :)**2*omega_grid_w(i_freq)
            end do
            call xml_elem('C6_vv_nm', C6_vv_nm)
        end if
    end if
    if (present(dvoldr)) then
        call sync_vector(dvoldr, size(dvoldr))
        call xml_elem('dvoldr', dvoldr)
        if (use_pulay) then
            dvoldr = dvoldr + hirshfeld_pulay_forces
            call xml_elem('dvoldr_total', dvoldr)
        end if
        if (present(dvoldh)) then
            call sync_vector(dvoldh, size(dvoldh))
            h_inv = inverted3x3(lattice_vector)
            forall (i_atom = 1:n_atoms)
                dvoldh(i_atom, :, :) = matmul(h_inv, dvoldh(i_atom, :, :))
            end forall
            call xml_elem('dvoldh', dvoldh)
            if (use_pulay) then
                dvoldh = dvoldh + hirshfeld_pulay_stress
                call xml_elem('dvoldh_total', dvoldh)
            end if
        end if
    end if

    if (grid_out) call xml_close_file(xml_file_grid)
end subroutine run_hirshfeld

elemental real(8) function logistic(x, w, x0) result(y)
    real(8), intent(in) :: x, w, x0

    y = 1d0/(1d0+exp(-4/w*(x-x0)))
end function

elemental real(8) function scanintp(x, c) result(y)
    real(8), intent(in) :: x, c

    if (x < 0) then
        y = 1
    else if (x > 1) then
        y = 0
    else
        y = exp(-c*x/(1-x))
    end if
end function

elemental real(8) function nm_cutoff(ion_pot, alpha) result(cutoff)
    real(8), intent(in) :: ion_pot, alpha

    cutoff = 1d0 - ( &
        (1d0 - logistic(ion_pot, w=1d0/hartree, x0=5d0/hartree)) &
        * (1d0 - scanintp(alpha-3*sqrt(ion_pot), c=0.1d0)) &
    )
end function

elemental function vv_pol_func(rho, rho_grad, omega, C) result(alpha)
    real(8), intent(in) :: rho, rho_grad, omega, C
    real(8) :: alpha

    alpha = rho/(4*pi/3*rho+C*(rho_grad/rho)**4+omega**2)
end function

elemental function terf(r, r0, a)
    real(8), intent(in) :: r, r0, a
    real(8) :: terf

    if (r > 0) then
        terf = 0.5d0*(erf(a*(r+r0))+erf(a*(r-r0)))
    else
        terf = 0.d0
    end if
end function

subroutine contract_with_dVdn(x, y)
    ! Calculates $y(\mathbf r)=\sum_{i\in\text{atoms}}x_jw^\text{H}_j(r)|\mathbf
    ! r-\mathbf R_j|^3$.
    real(8), intent(in) :: x(n_atoms)
    real(8), intent(out) :: y(n_full_points)

    real(8) :: r(3), weight
    integer :: i_full_point, i_batch, i_batch_point, i_atom
    type(grid_point) :: point

    i_full_point = 0
    do i_batch = 1, n_my_batches
        do i_batch_point = 1, batches(i_batch)%size
            i_full_point = i_full_point+1
            point = batches(i_batch)%points(i_batch_point)
            i_atom = point%index_atom
            if (empty(i_atom)) cycle
            call evaluate_hirshfeld_factor('N', point, weight)
            r(:) = point%coords(:)
            if (.not. weight >= 0.d0) cycle
            y(i_full_point) = x(i_atom)*weight
        end do
    end do
end subroutine

subroutine evaluate_hirshfeld_factor( &
        mode, point, weight, part_weight, dweightdr, dweightdh, shift, strain &
    )
    ! Calculates $w^\text{H}_i(\mathbf r)|\mathbf{r}-\mathbf{R}_i|^3$ if `mode
    ! == 'N'` or multiplied with $w_\text{int}(\mathbf r)$ if `mode == 'W'`.
    ! Optionally calcualte (dw/dR*r^3+w*dr^3/dR)*w_\text{int}.

    use dimensions, only: n_centers_basis_integrals

    type(grid_point), intent(in) :: point
    character(len=1), intent(in) :: mode
    real(8), intent(out) :: weight
    real(8), intent(out), optional :: part_weight
    real(8), intent(out), optional :: dweightdr(n_atoms, 3)
    real(8), intent(out), optional :: dweightdh(3, 3)
    real(8), intent(in), optional :: shift(n_atoms, 3)
    real(8), intent(in), optional :: strain(3, 3)

    real(8) :: &
        r(3), r_log, r_diff(3), x(3), &
        dist_to_atoms(n_centers_basis_integrals), &
        rho_0_i(n_centers_basis_integrals), &
        drho_0_i(n_centers_basis_integrals), &
        weight_hirsh, &
        rho_0_total, weight_int, &
        dweight_hirsh_dR(n_atoms, 3), dr3dr(3), &
        dweight_hirsh_dh(3, 3)
    integer :: &
        specie_i, i_atom, i_radial, i_angular, j_atom, j_center_integrals, &
        specie_j, j_center, i_vec

    weight = 0.d0
    r = point%coords
    do j_center_integrals = 1, n_centers_basis_integrals
        j_center = centers_basis_integrals(j_center_integrals)
        r_diff(:) = r(:)-coords_center(:, j_center)
        if (present(shift)) then
            r_diff(:) = r_diff(:)-shift(center_to_atom(j_center), :)
        else if (present(strain)) then
            r_diff(:) = matmul(strain, r_diff)
        end if
        dist_to_atoms(j_center_integrals) = sqrt(sum(r_diff**2))
    enddo
    if ( &
        any(dist_to_atoms < r_grid_min(species_center(centers_basis_integrals))) &
    ) return
    do j_center_integrals = 1, n_centers_basis_integrals
        j_center = centers_basis_integrals(j_center_integrals)
        if (empty(center_to_atom(j_center))) then
            rho_0_i(j_center_integrals) = 0.d0
            drho_0_i(j_center_integrals) = 0.d0
            cycle
        end if
        specie_j = species_center(j_center)
        r_log = invert_log_grid( &
            dist_to_atoms(j_center_integrals), &
            r_grid_min(specie_j), &
            r_grid_inc(specie_j) &
        )
        rho_0_i(j_center_integrals) = 1.d0/(4*pi)*val_spline( &
            r_log, free_rho_spl(:, :, specie_j), n_grid(specie_j) &
        )
        drho_0_i(j_center_integrals) = -1.d0/(4*pi)*val_spline( &
            r_log, free_drho_dr_spl(:, :, specie_j), n_grid(specie_j) &
        )
    enddo
    where (rho_0_i < 0.d0) rho_0_i = 0.d0
    if (all(rho_0_i <= 0.d0)) return
    i_atom = point%index_atom
    i_radial = point%index_radial
    i_angular = point%index_angular
    specie_i = species(i_atom)
    rho_0_total = sum(rho_0_i)
    weight_hirsh = rho_0_i(i_atom)/rho_0_total
    part_weight = weight_hirsh
    weight = weight_hirsh*dist_to_atoms(i_atom)**3
    if (.not. weight > 0.d0) return
    call xml_open('point', file=xml_file_grid)
    call xml_elem('x', r(1), file=xml_file_grid)
    call xml_elem('y', r(2), file=xml_file_grid)
    call xml_elem('z', r(3), file=xml_file_grid)
    call xml_elem('i_radial', i_radial, file=xml_file_grid)
    call xml_elem('i_angular', i_angular, file=xml_file_grid)
    call xml_elem('rho_0', rho_0_i(i_atom), file=xml_file_grid)
    call xml_elem('rho_0_total', rho_0_total, file=xml_file_grid)
    call xml_elem('drho_0', drho_0_i(i_atom), file=xml_file_grid)
    call xml_elem('weight_hirsh', weight_hirsh, file=xml_file_grid)
    select case (mode)
    case ('N')
    case ('W')
        weight_int = 4*pi*r_radial(i_radial, specie_i)**2 &
            * w_radial(i_radial, specie_i) &
            * w_angular(i_angular, i_radial, specie_i)
        call xml_elem('weight_int', weight_int, file=xml_file_grid)
        weight = weight*weight_int
        part_weight = part_weight*weight_int
    end select
    call xml_elem('weight', weight, file=xml_file_grid)
    call xml_elem('part_weight', part_weight, file=xml_file_grid)
    if (present(dweightdr)) then
        r_diff(:) = r(:)-coords(:, i_atom)
        dr3dr(:) = -3*dist_to_atoms(i_atom)*r_diff(:)
        x(:) = drho_0_i(i_atom)/rho_0_total*r_diff(:)/dist_to_atoms(i_atom)
        dweight_hirsh_dR(:, :) = 0.d0
        dweight_hirsh_dR(i_atom, :) = dweight_hirsh_dR(i_atom, :) + x(:)
        if (n_periodic > 0) then
            forall (i_vec = 1:3) dweight_hirsh_dh(i_vec, :) = -x(:)*r_diff(i_vec)
        end if
        do j_center_integrals = 1, n_centers_basis_integrals
            j_center = centers_basis_integrals(j_center_integrals)
            r_diff(:) = r(:)-coords_center(:, j_center)
            j_atom = center_to_atom(j_center)
            x(:) = -weight_hirsh*drho_0_i(j_center_integrals)/rho_0_total * &
                r_diff(:)/dist_to_atoms(j_center_integrals)
            dweight_hirsh_dR(j_atom, :) = dweight_hirsh_dR(j_atom, :) + x(:)
            if (n_periodic > 0) then
                forall (i_vec = 1:3)
                    dweight_hirsh_dh(i_vec, :) = dweight_hirsh_dh(i_vec, :) - &
                        x(:)*r_diff(i_vec)
                end forall
            end if
        end do
        x(:) = weight_hirsh*dr3dr(:)
        dweightdr(:, :) = dweight_hirsh_dR(:, :)*dist_to_atoms(i_atom)**3
        dweightdr(i_atom, :) = dweightdr(i_atom, :) + x(:)
        call xml_elem('dweightdr', dweightdr, file=xml_file_grid)
        if (n_periodic > 0) then
            r_diff = r(:)-coords(:, i_atom)
            dweightdh(:, :) = dweight_hirsh_dh(:, :)*dist_to_atoms(i_atom)**3
            forall (i_vec = 1:3)
                dweightdh(i_vec, :) = dweightdh(i_vec, :) - x(:)*r_diff(i_vec)
            end forall
            call xml_elem('dweightdh', dweightdr, file=xml_file_grid)
        end if
        select case (mode)
        case ('N')
        case ('W')
            dweightdr = dweightdr*weight_int
            if (n_periodic > 0) dweightdh = dweightdh*weight_int
        end select
    end if
end subroutine

subroutine evaluate_free_atom_quantities(volumes, alpha_0_vv, C6_vv)
    real(8), intent(out) :: volumes(n_species)
    real(8), intent(out), optional :: alpha_0_vv(n_species)
    real(8), intent(out), optional :: C6_vv(n_species)

    integer :: i_radial, i_specie, i_freq
    real(8) :: rho, r, r_log, weight_int, rho_grad
    real(8), allocatable :: vv_pols(:, :)
    real(8), allocatable :: omega_grid(:), omega_grid_w(:)
    type(Xml_file_t) :: xml_file_grid

    if (grid_out .and. myid == 0) then
        call xml_open_file("grid_free.xml", "grid", xml_file_grid)
    end if
    volumes(:) = 0.d0
    if (present(alpha_0_vv) .or. present(C6_vv)) then
        allocate (omega_grid(0:n_freq), source=0.d0)
        allocate (omega_grid_w(0:n_freq), source=0.d0)
        call get_freq_grid(n_freq, omega_grid(1:n_freq), omega_grid_w(1:n_freq))
        allocate (vv_pols(0:n_freq, n_species), source=0.d0)
    end if
    do i_specie = 1, n_species
        call xml_open('species', file=xml_file_grid)
        call xml_elem('name', trim(species_name(i_specie)), file=xml_file_grid)
        do i_radial = 1, n_radial(i_specie)
            r = r_radial(i_radial, i_specie)
            r_log = invert_log_grid( &
                r, r_grid_min(i_specie), r_grid_inc(i_specie) &
            )
            rho = val_spline( &
                r_log, free_rho_spl(:, :, i_specie), n_grid(i_specie) &
            )/(4d0*pi)
            if (rho <= 0.d0) cycle
            call xml_open('point', file=xml_file_grid)
            call xml_elem('r', r, file=xml_file_grid)
            call xml_elem('rho', rho, file=xml_file_grid)
            weight_int = 4d0*pi*r_radial(i_radial, i_specie)**2 &
                * w_radial(i_radial, i_specie)
            call xml_elem('weight_int', weight_int, file=xml_file_grid)
            rho_grad = val_spline(&
                r_log, free_drho_dr_spl(:, :, i_specie), n_grid(i_specie) &
            )/(4d0*pi)
            call xml_elem('rho_grad', rho_grad, file=xml_file_grid)
            call xml_close(file=xml_file_grid)
            volumes(i_specie) = volumes(i_specie) &
                + rho*r_radial(i_radial, i_specie)**3*weight_int
            if (present(alpha_0_vv) .or. present(C6_vv)) then
                vv_pols(:, i_specie) = vv_pols(:, i_specie) &
                    + vv_pol_func(rho, rho_grad, omega_grid, C_vv)*weight_int
            end if
        end do
        call xml_close(file=xml_file_grid)
    end do
    call xml_open('free_atoms')
    call xml_elem('volumes', volumes)
    call xml_elem('species', species)
    if (present(alpha_0_vv) .or. present(C6_vv)) then
        call xml_elem('vv_pols', vv_pols)
        if (present(alpha_0_vv)) alpha_0_vv = vv_pols(0, :)
        if (present(C6_vv)) then
            C6_vv(:) = 0.d0
            do i_freq = 1, n_freq
                C6_vv = C6_vv + 3d0/pi*vv_pols(i_freq, :)**2*omega_grid_w(i_freq)
            end do
        end if
    end if
    if (present(alpha_0_vv) .and. present(C6_vv)) then
        ! TODO use nonspherical atom densities calculated by atom_sphere for lighter atoms
        ! for lighter elements, the error in VV polarizability from using a
        ! spherical free atom densities is too large, so we use hardcoded values
        ! precalculated by atom_sphere atomic solver with the PBE functional and
        ! C = 0.0093 (so this has to remain fixed for now). for heavier elements,
        ! this doesn't seem to be an issue, but relativity starts to play a role,
        ! so sratom should be used for those
        do i_specie = 1, n_species
            associate (a0 => alpha_0_vv(i_specie), C6 => C6_vv(i_specie))
                select case (species_name(i_specie))
                case ('H'); a0 = 4.601300; C6 = 6.916624
                case ('He'); a0 = 1.307572; C6 = 1.397397
                case ('Li'); a0 = 79.382990; C6 = 575.475693
                case ('Be'); a0 = 31.560229; C6 = 179.386328
                case ('B'); a0 = 18.938659; C6 = 98.951839
                case ('C'); a0 = 11.205754; C6 = 50.583489
                case ('N'); a0 = 7.116014; C6 = 27.934184
                case ('O'); a0 = 5.269107; C6 = 18.739490
                case ('F'); a0 = 3.834770; C6 = 12.341105
                case ('Ne'); a0 = 2.853014; C6 = 8.315834
                case ('Na'); a0 = 85.572388; C6 = 650.263509
                case ('Mg'); a0 = 52.277221; C6 = 393.328069
                case ('Al'); a0 = 48.471309; C6 = 405.658261
                case ('Si'); a0 = 33.849222; C6 = 276.749048
                case ('P'); a0 = 23.790392; C6 = 184.279342
                case ('S'); a0 = 18.445013; C6 = 135.418497
                case ('Cl'); a0 = 14.102330; C6 = 96.937274
                case ('Ar'); a0 = 10.928427; C6 = 70.413260
                end select
            end associate
        end do
    end if
    if (present(alpha_0_vv)) call xml_elem('alpha_0_vv', alpha_0_vv)
    if (present(C6_vv)) call xml_elem('C6_vv', C6_vv)
    call xml_close()
    if (grid_out .and. myid == 0) call xml_close_file(xml_file_grid)
end subroutine

subroutine numerical_hirshfeld_deriv(rho)
    real(8), intent(in) :: rho(:, :)

    real(8) :: dvoldr(n_atoms, n_atoms, 3)
    real(8) :: dvoldh(n_atoms, 3, 3)
    real(8) :: shift(n_atoms, 3), strain(3, 3)
    real(8) :: delta = 1.d-5
    real(8) :: dvolumes(-3:3, n_atoms)
    real(8) :: h_inv(3, 3)

    integer :: j_atom, i_xyz, i_step, i_atom, j_xyz, i

    do j_atom = 1, n_atoms
        do i_xyz = 1, 3
            shift(:, :) = 0.d0
            do i_step = -3, 3
                shift(j_atom, i_xyz) = i_step*delta
                call run_hirshfeld(rho, dvolumes(i_step, :), shift=shift)
            end do
            do i_atom = 1, n_atoms
                dvoldr(i_atom, j_atom, i_xyz) = diff7(dvolumes(:, i_atom), delta)
            end do
        end do
    end do
    call xml_elem('dvoldr_numeric', dvoldr)
    if (n_periodic > 0) then
        do i_xyz = 1, 3
            do j_xyz = 1, 3
                do i_step = -3, 3
                    strain(:, :) = 0.d0
                    forall (i = 1:3) strain(i, i) = 1.d0
                    strain(i_xyz, j_xyz) = strain(i_xyz, j_xyz) + i_step*delta
                    call run_hirshfeld(rho, dvolumes(i_step, :), strain=strain)
                end do
                do i_atom = 1, n_atoms
                    dvoldh(i_atom, i_xyz, j_xyz) = diff7(dvolumes(:, i_atom), delta)
                end do
            end do
        end do
        h_inv = inverted3x3(lattice_vector)
        forall (i_atom = 1:n_atoms)
            dvoldh(i_atom, :, :) = matmul(h_inv, dvoldh(i_atom, :, :))
        end forall
        call xml_elem('dvoldh_numeric', dvoldh)
    end if
end subroutine

subroutine eval_hirsh_vol_pulay_deriv_dm( &
        KS_eigenvector, KS_eigenvector_cmplx, occ_numbers, &
        densmat_sparse_work, shift, strain &
    )
    ! Calculate the Pulay part of Hirshfeld volume derivatives $dV_i/dR_j$,
    ! which is $\int(dn/dR_j)w_ir_i^3=2\sum_{\mu\nu}P_{\mu\nu}\langle
    ! \phi_\mu|w_ir_i^3|\nabla\phi_\nu\rangle$. Since the Hirshfeld volume of an
    ! atom is integrated only over points of that atom, the same has to happen
    ! also for the forces, and hence at each point only one Hirshfeld factor
    ! $w_ir_i^3$ needs to be calculated. This is a huge time saving. A python
    ! pseudocode for the routine is
    !
    ! ``` python
    ! for batch in batches:
    !   pruned_basis = basis.pruned(batch)
    !   pt_to_at = [point.atom for point in batch]
    !   construct(bra[:n_basis, :n_points])
    !   construct(ket[:n_basis, :n_points, :3])
    !   for atom in unique(pt_to_at):
    !     for dim in range(3):
    !       braket[:n_basis, :n_basis] = dgemm(
    ! 		bra[:n_basis, pt_to_at == atom],
    ! 		ket[:n_basis, pt_to_at == atom])
    !       for mu in pruned_basis:
    !         for nu in pruned_basis:
    !           hirshfeld_forces[atom, nu.atom, dim] \
    !           	+= dens_mat[mu.full_idx, nu.full_idx]*psi_grad_psi[mu, nu]
    ! ```
    use density_matrix_evaluation, only: evaluate_densmat
    use species_data, only: l_shell_max
    use dimensions, only: &
        n_basis, n_states, n_spin, n_k_points_task, n_k_points, &
        n_hamiltonian_matrix_size, n_centers_integrals, &
        n_centers_basis_I, n_my_batches, &
        n_centers, n_max_batch_size, n_basis_fns, l_wave_max
    use basis, only: basis_deriv_ordered, basis_wave_ordered
    use mpi_tasks, only: check_allocation
    use pbc_lists, only: &
        inv_centers_basis_integrals, kweight_occs, de_kweight_occs, &
        cbasis_to_center
    use localorb_io, only: OL_norm, use_unit
    use runtime_choices, only: use_symmetry_reduced_spg

    implicit none

    real(8), intent(in) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex(8), intent(in) :: KS_eigenvector_cmplx(n_basis, n_states, n_spin, n_k_points_task)
    real(8), intent(inout) :: occ_numbers(n_states, n_spin, n_k_points)
    real(8) :: densmat_sparse_work(n_hamiltonian_matrix_size, n_spin)
    real(8), intent(in), optional :: shift(n_atoms, 3)
    real(8), intent(in), optional :: strain(3, 3)

    external :: map_to_center_cell, tab_atom_centered_coords_p0, prune_basis_p0, &
        prune_radial_basis_p0, tab_local_geometry_p0, tab_trigonom_p0, &
        tab_wave_ylm_p0, evaluate_radial_functions_p0, evaluate_waves_p0, &
        tab_gradient_ylm_p0, evaluate_wave_gradient_p0

    integer :: &
        i_batch, i_spin, i_index, i_l, i_m, i_point, i_dim, &
        i_atom, i_center, i_center_integrals, nu, i_xyz, j_xyz, j_center
    integer :: &
        n_points, n_pruned, n_compute_fns, n_compute_atoms, n_dummy, &
        n_points_atom, n_dims, n_max_compute_atoms, n_max_compute_fns_dens, &
        n_max_compute_dens
    integer :: &
        i_basis_fns(n_basis_fns*n_centers_integrals), &
        i_basis_fns_inv(n_basis_fns, n_centers), &
        i_atom_fns(n_basis_fns*n_centers_integrals), &
        atom_index(n_centers_integrals), &
        atom_index_inv(n_centers), &
        spline_array_start(n_centers_integrals), &
        spline_array_end(n_centers_integrals), &
        point_to_atom(n_max_batch_size), &
        points_atom(n_max_batch_size)
    real(8) :: &
        r(3), r_diff(3), &
        dummy_matrix(1, 1), &
        dist_tab(n_centers_integrals), &
        dist_tab_sq(n_centers_integrals, n_max_batch_size), &
        dir_tab(3, n_centers_integrals, n_max_batch_size), &
        hirshfeld_factor(n_max_batch_size), hirsh_stress(1), h_inv(3, 3), &
        r_grid_min_sq(n_species)
    logical :: atom_in_batch(n_atoms)
    type(grid_point) :: point
    logical :: finite_diff

    integer, allocatable :: index_lm(:, :), i_basis(:), i_basis_fn(:)
    real(8), allocatable :: &
        ylm_tab(:, :), &
        dylm_dtheta_tab(:, :), &
        scaled_dylm_dphi_tab(:, :), &
        bra(:, :),  &
        ket(:, :, :), &
        radial_wave(:), &
        radial_wave_deriv(:), &
        braket(:, :), &
        gradient_basis_wave(:, :, :), &
        trigonom_tab(:, :), &
        i_r(:)

    call localorb_info( &
        "Evaluating Hirshfeld Pulay forces using the density matrix.", &
        use_unit, '(2X, A)', OL_norm)

    r_grid_min_sq(:) = r_grid_min(:)**2

    finite_diff = present(shift) .or. present(strain)
    hirshfeld_pulay_forces(:, :, :) = 0.d0
    hirshfeld_pulay_stress(:, :, :) = 0.d0
    if (finite_diff) then
        n_dims = 1  ! only |nu> needs to be computed
    else if (n_periodic == 0) then
        n_dims = 3  ! grad_i|nu> needs to be computed
    else
        n_dims = 12  ! grad_i|nu> & grad_i|nu>*(r_j-R_j) are computed
        ! 1-3: grad_i|nu>, 4-6: grad_i|nu>*(r_1-R_1),...
    end if

    call get_max_dimensions( &
        n_max_compute_atoms, n_max_compute_dens, n_max_compute_fns_dens, &
        shift, strain &
    )

    allocate (i_basis(n_centers_basis_I))
    allocate (i_basis_fn(n_centers_basis_I))
    allocate (bra(n_max_compute_dens, n_max_batch_size))
    allocate (ket(n_max_compute_dens, n_max_batch_size, n_dims))
    allocate (braket(n_max_compute_dens, n_max_compute_dens))
    allocate (radial_wave(n_max_compute_fns_dens))
    allocate (index_lm(-l_wave_max:l_wave_max, 0:l_wave_max))
    allocate (ylm_tab((l_wave_max+1)**2, n_max_compute_atoms))
    allocate (trigonom_tab(4, n_max_compute_atoms))
    allocate (i_r(n_max_compute_atoms))
    if (.not. finite_diff) then
        allocate (radial_wave_deriv(n_max_compute_fns_dens))
        allocate (gradient_basis_wave(n_max_compute_dens, 3, n_max_batch_size))
        allocate (dylm_dtheta_tab((l_wave_max+1)**2, n_max_compute_atoms))
        allocate (scaled_dylm_dphi_tab((l_wave_max+1)**2, n_max_compute_atoms))
    end if

    densmat_sparse_work(:, :) = 0.d0  ! TODO is this necessary?
    call kweight_occs('eval_hirsh_vol_pulay_deriv_dm', occ_numbers)
    do i_spin = 1, n_spin
        if (use_symmetry_reduced_spg) then
            call evaluate_densmat_sym( &
                KS_eigenvector, KS_eigenvector_cmplx, &
                occ_numbers, dummy_matrix, densmat_sparse_work(:, i_spin), i_spin, .true. &
            )
        else
            call evaluate_densmat( &
                KS_eigenvector, KS_eigenvector_cmplx, &
                occ_numbers, dummy_matrix, densmat_sparse_work(:, i_spin), i_spin, .true. &
            )
        endif
    end do
    call de_kweight_occs('eval_hirsh_vol_pulay_deriv_dm', occ_numbers)

    i_index = 0
    do i_l = 0, l_wave_max
        do i_m = -i_l, i_l
            i_index = i_index+1
            index_lm(i_m, i_l) = i_index
        end do
    end do

    batch_loop: do i_batch = 1, n_my_batches
        n_points = batches(i_batch)%size
        n_pruned = 0
        i_basis(:) = 0
        atom_in_batch(:) = .false.
        do i_point = 1, batches(i_batch)%size
            point = batches(i_batch)%points(i_point)
            point_to_atom(i_point) = point%index_atom
            atom_in_batch(point%index_atom) = .true.
            r = point%coords
            if (n_periodic > 0) call map_to_center_cell(r)
            do i_center_integrals = 1, n_centers_integrals
                i_center = centers_basis_integrals(i_center_integrals)
                r_diff(:) = r(:)-coords_center(:, i_center)
                if (present(shift)) then
                    r_diff(:) = r_diff(:)-shift(center_to_atom(i_center), :)
                else if (present(strain)) then
                    r_diff(:) = matmul(strain, r_diff)
                end if
                dir_tab(:, i_center_integrals, i_point) = r_diff(:)
                dist_tab_sq(i_center_integrals, i_point) = sum(r_diff**2)
            end do
            call prune_basis_p0( &
                dist_tab_sq(:, i_point), &
                n_dummy, n_pruned, i_basis, &
                n_centers_basis_I, n_centers_integrals, &
                inv_centers_basis_integrals &
            )
        end do
        if (n_pruned == 0) cycle
        do i_point = 1, batches(i_batch)%size
            if (any(dist_tab_sq(:, i_point) < r_grid_min_sq( &
                species_center(centers_basis_integrals(1:n_centers_integrals)) &
            ))) then
                hirshfeld_factor(i_point) = 0.d0
                cycle
            end if
            point = batches(i_batch)%points(i_point)
            n_compute_atoms = 0
            n_compute_fns = 0
            i_basis_fns_inv(:, :) = 0
            atom_index_inv(:) = 0
            call prune_radial_basis_p0( &
                dist_tab_sq(:, i_point), dist_tab, dir_tab(:, :, i_point), &
                n_compute_atoms, atom_index, atom_index_inv, n_compute_fns, &
                i_basis_fns, i_basis_fns_inv, i_atom_fns, spline_array_start, &
                spline_array_end, n_centers_integrals, centers_basis_integrals &
            )
            call tab_local_geometry_p0( &
                dist_tab_sq(:, i_point), n_compute_atoms, atom_index, &
                dir_tab(:, :, i_point), dist_tab, i_r &
            )
            call tab_trigonom_p0( &
                n_compute_atoms, dir_tab(:, :, i_point), trigonom_tab &
            )
            call tab_wave_ylm_p0( &
                n_compute_atoms, atom_index, trigonom_tab, l_shell_max, &
                l_wave_max, ylm_tab &
            )
            call evaluate_radial_functions_p0( &
                spline_array_start, spline_array_end, n_compute_atoms, &
                n_compute_fns, dist_tab, i_r, atom_index, i_basis_fns_inv, &
                basis_wave_ordered, radial_wave, .false., n_pruned, &
                n_max_compute_fns_dens &
            )
            call evaluate_waves_p0( &
                l_wave_max, ylm_tab, dist_tab, index_lm, n_pruned, i_basis, &
                radial_wave, bra(:, i_point), n_compute_atoms, atom_index_inv, &
                n_compute_fns, i_basis_fns_inv, n_max_compute_fns_dens &
            )
            if (finite_diff) then
                ket(:n_pruned, i_point, 1) = bra(:n_pruned, i_point)
            else
                call evaluate_radial_functions_p0( &
                    spline_array_start, spline_array_end, n_compute_atoms, &
                    n_compute_fns, dist_tab, i_r, atom_index, i_basis_fns_inv, &
                    basis_deriv_ordered, radial_wave_deriv, .true., n_pruned, &
                    n_max_compute_fns_dens &
                )
                call tab_gradient_ylm_p0( &
                    trigonom_tab, l_shell_max, l_wave_max, n_compute_atoms, &
                    atom_index, ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab &
                )
                call evaluate_wave_gradient_p0( &
                    dist_tab, dir_tab(:, :, i_point), trigonom_tab, l_wave_max, &
                    ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab, index_lm, &
                    n_pruned, i_basis, radial_wave, radial_wave_deriv, &
                    gradient_basis_wave(:, :, i_point), n_compute_atoms, &
                    atom_index_inv, n_compute_fns, i_basis_fns_inv, &
                    n_max_compute_fns_dens &
                )
                ! For efficient array slicing, we calculated ket(n_basis, n_dim, n_points),
                ! with the first two dimensions packed, but now we need
                ! ket(n_basis, n_points, n_dim) (and hence unpacked).
                forall (i_xyz = 1:3)
                    ket(:n_pruned, i_point, i_xyz) = permute_gradient( &
                        n_pruned, gradient_basis_wave(:, :, i_point), i_xyz &
                    )
                end forall
                if (n_periodic > 0) then
                    do nu = 1, n_pruned
                        j_center = atom_index_inv(cbasis_to_center(i_basis(nu)))
                        if (j_center == 0) then
                            ket(nu, i_point, 4:) = 0.d0
                            cycle
                        end if
                        forall (j_xyz = 1:3)
                            ket(nu, i_point, 3*j_xyz+1:3*j_xyz+3) = &
                                -ket(nu, i_point, 1:3) &
                                *dir_tab(j_xyz, j_center, i_point) &
                                *dist_tab(j_center)
                        end forall
                    end do
                end if
            end if
            call evaluate_hirshfeld_factor('W', point, hirshfeld_factor(i_point))
        end do
        if (.not. finite_diff) then
            hirshfeld_factor = 2*hirshfeld_factor ! d(mu*nu)=2*d(mu)*nu
            ket = -ket ! d(mu)=-grad(mu)
        end if
        do i_atom = 1, n_atoms
            if (.not. atom_in_batch(i_atom)) cycle
            call find(point_to_atom(:n_points) == i_atom, n_points_atom, points_atom)
            do i_dim  = 1, n_dims
                call integrate_braket( &
                    n_pruned, n_points_atom, &
                    bra(:n_pruned, points_atom(:n_points_atom)), &
                    ket(:n_pruned, points_atom(:n_points_atom), i_dim), &
                    hirshfeld_factor(points_atom(:n_points_atom)), braket &
                )
                if (i_dim <= 3) then
                    call update_hirshfeld_pulay_forces( &
                        n_pruned, sum(densmat_sparse_work, 2), braket, i_basis, &
                        hirshfeld_pulay_forces(:, i_atom, i_dim), forces=.true. &
                    )
                else
                    i_xyz = (i_dim-1)/3
                    j_xyz = i_dim-i_xyz*3
                    hirsh_stress(:) = 0.d0
                    call update_hirshfeld_pulay_forces( &
                        n_pruned, sum(densmat_sparse_work, 2), braket, i_basis, &
                        hirsh_stress, forces=.false. &
                    )
                    hirshfeld_pulay_stress(i_atom, i_xyz, j_xyz) = &
                        hirshfeld_pulay_stress(i_atom, i_xyz, j_xyz) &
                        + hirsh_stress(1)
                end if
            end do
        end do
    end do batch_loop

    call sync_vector(hirshfeld_pulay_forces, size(hirshfeld_pulay_forces))
    ! For efficient array slicing, we calculated dV_i/dR_ja as
    ! hirshfeld_forces(j, i, a), but we need hirshfeld_forces(i, j, a).
    hirshfeld_pulay_forces = reshape( &
        hirshfeld_pulay_forces, (/ n_atoms, n_atoms, 3 /), order=(/ 2, 1, 3 /) &
    )
    if (n_periodic > 0) then
        call sync_vector(hirshfeld_pulay_stress, size(hirshfeld_pulay_stress))
        h_inv = inverted3x3(lattice_vector)
        forall (i_atom = 1:n_atoms)
            hirshfeld_pulay_stress(i_atom, :, :) = &
                matmul(h_inv, hirshfeld_pulay_stress(i_atom, :, :))
        end forall
    end if

    if (.not. finite_diff) then
        call xml_elem('hirshfeld_pulay_forces', hirshfeld_pulay_forces)
        if (n_periodic > 0) then
            call xml_elem('hirshfeld_pulay_stress', hirshfeld_pulay_stress)
        end if
    end if

    deallocate (i_basis)
    deallocate (bra)
    deallocate (ket)
    deallocate (braket)
    deallocate (radial_wave)
    deallocate (index_lm)
    deallocate (ylm_tab)
    deallocate (trigonom_tab)
    deallocate (i_r)
    if (.not. finite_diff) then
        deallocate (radial_wave_deriv)
        deallocate (gradient_basis_wave)
        deallocate (dylm_dtheta_tab)
        deallocate (scaled_dylm_dphi_tab)
    end if

end subroutine

subroutine find(mask, n_found, i_found)
    logical, intent(in) :: mask(:)
    integer, intent(out) :: n_found
    integer, intent(out) :: i_found(:)

    integer :: i

    n_found = 0
    do i = 1, size(mask)
        if (mask(i)) then
            n_found = n_found+1
            i_found(n_found) = i
        end if
    end do
end subroutine

pure function permute_gradient(n_pruned, arr_basis_dim, i_dim) result(arr_basis)
    integer, intent(in) :: n_pruned
    real(8), intent(in) :: arr_basis_dim(n_pruned, 3)
    integer, intent(in) :: i_dim
    real(8) :: arr_basis(n_pruned)

    arr_basis(:) = arr_basis_dim(:, i_dim)
end function

subroutine integrate_braket(n_pruned, n_points, bra, ket, weight_int, braket)
    ! Calculate $\langle\hat O_1\phi_\mu|\hat
    ! O_2\phi_\nu\rangle_\text{batch}=\sum_p \hat O_1\phi_\mu(\mathbf
    ! r_p)\hat O_2\phi_\nu(\mathbf r_p)w_p$ where $\hat O_1\phi_\mu$ is
    ! supplied as `bra` and $\hat O_2\phi_\nu$ as `ket` and $w_p$ is some
    ! quadrature weight. $\hat O_1$, $\hat O_2$ might not be Hermitian,
    ! hence the result may not be symmetric.
    integer, intent(in) :: n_pruned, n_points
    real(8), intent(in) :: bra(n_pruned, n_points)
    real(8), intent(in) :: ket(n_pruned, n_points)
    real(8), intent(in) :: weight_int(n_points)
    real(8), intent(out) :: braket(n_pruned, n_pruned)

    external :: DGEMM

    integer :: i_point
    real(8) :: ket_times_w(n_pruned, n_points)

    forall (i_point = 1:n_points)
        ket_times_w(:, i_point) = ket(:, i_point)*weight_int(i_point)
    end forall
    call DGEMM( &
        'N', 'T', n_pruned, n_pruned, n_points, 1.0d0, &
        bra, n_pruned, ket_times_w, n_pruned, &
        0.0d0, braket, n_pruned &
    )
    ! braket = matmul(bra, transpose(ket_times_w))
end subroutine

subroutine update_hirshfeld_pulay_forces( &
        n_pruned, density_matrix, braket, i_basis, pulay_force, forces &
    )
    ! Add a Pulay contribution to $\partial V_i^P/\partial
    ! R_{j\alpha}=\int(\partial n/\partial R_{j\alpha})w_ir_i^3$, from the
    ! current batch. This is calculated as
    !
    ! $$ \frac{\partial V_i^P}{\partial R_{j\alpha}}=-2\sum_{\mu,\nu\in
    ! j}P_{\mu\nu}\langle\phi_\mu|w_ir_i^3|\nabla_\alpha\phi_\nu\rangle $$
    use runtime_choices, only: PM_none, PM_index, packed_matrix_format
    use pbc_lists, only: n_cells, cbasis_to_atom, cbasis_to_basis, &
        cbasis_to_center, center_to_cell, position_in_hamiltonian, &
        index_hamiltonian, column_index_hamiltonian
    use dimensions, only: n_hamiltonian_matrix_size

    integer, intent(in) :: n_pruned
    real(8), intent(in) :: density_matrix(n_hamiltonian_matrix_size)
    real(8), intent(in) :: braket(n_pruned, n_pruned)
    integer, intent(in) :: i_basis(n_pruned)  ! pruned to full basis mapping
    real(8), intent(inout) :: pulay_force(:)  ! single dimension
    logical, intent(in) :: forces

    integer :: &
        mu, nu, i_offset, mu_nu, i_start, &
        i_end, mu_full, nu_full, i_cell_mu, i_cell_nu, i_cell, &
        mu_nu_current, idx
    integer :: i_starts(n_cells), i_ends(n_cells)

    idx = 1  ! doesn't change when evaluating stress
    select case (packed_matrix_format)
    case (PM_none)
        mu_nu = 0
        do mu = 1, n_pruned
            mu_full = i_basis(mu)
            i_offset = (mu_full-1)*mu_full/2
            do nu = 1, mu
                nu_full = i_basis(nu)
                mu_nu = i_offset + nu_full
                if (forces) idx = cbasis_to_atom(nu_full)
                pulay_force(idx) = pulay_force(idx) &
                    + density_matrix(mu_nu)*braket(mu, nu)
                if (forces .and. mu_full /= nu_full) then
                    idx = cbasis_to_atom(mu_full)
                    pulay_force(idx) = pulay_force(idx) &
                        + density_matrix(mu_nu)*braket(nu, mu)
                end if
            end do
        end do
    case (PM_index)
        if (n_periodic == 0) then
            do mu = 1, n_pruned
                mu_full = i_basis(mu)
                i_start = index_hamiltonian(1, 1, mu_full)
                i_end = index_hamiltonian(2, 1, mu_full)
                do nu = 1, n_pruned
                    nu_full = i_basis(nu)
                    do mu_nu = i_start, i_end
                        if (column_index_hamiltonian(mu_nu) == nu_full) then
                            if (forces) idx = cbasis_to_atom(nu_full)
                            pulay_force(idx) = pulay_force(idx) &
                                + density_matrix(mu_nu)*braket(mu, nu)
                            if (mu_full /= nu_full) then
                                idx = cbasis_to_atom(mu_full)
                                pulay_force(idx) = pulay_force(idx) &
                                    + density_matrix(mu_nu)*braket(nu, mu)
                            end if
                            mu_nu_current = mu_nu
                            exit
                        else if (column_index_hamiltonian(mu_nu) > nu_full) then
                            mu_nu_current = mu_nu
                            exit
                        end if
                    end do
                    i_start = mu_nu_current
                end do
            end do
        else
            do mu = 1, n_pruned
                mu_full = cbasis_to_basis(i_basis(mu))
                i_cell_mu = center_to_cell(cbasis_to_center(i_basis(mu)))
                i_starts(:) = -1
                i_ends(:) = -1
                do i_cell_nu = 1, n_cells
                    i_cell = position_in_hamiltonian(i_cell_mu, i_cell_nu)
                    i_starts(i_cell_nu) = index_hamiltonian(1, i_cell, mu_full)
                    i_ends(i_cell_nu) = index_hamiltonian(2, i_cell, mu_full)
                end do
                do nu = 1, n_pruned
                    nu_full = cbasis_to_basis(i_basis(nu))
                    if (nu_full > mu_full) cycle
                    i_cell_nu = center_to_cell(cbasis_to_center(i_basis(nu)))
                    do mu_nu = i_starts(i_cell_nu), i_ends(i_cell_nu)
                        if (column_index_hamiltonian(mu_nu) == nu_full) then
                            if (forces) idx = cbasis_to_atom(nu_full)
                            pulay_force(idx) = pulay_force(idx) &
                                + density_matrix(mu_nu)*braket(mu, nu)
                            if (mu_full /= nu_full) then
                                if (forces) idx = cbasis_to_atom(mu_full)
                                pulay_force(idx) = pulay_force(idx) &
                                    + density_matrix(mu_nu)*braket(nu, mu)
                            end if
                            exit
                        else if (column_index_hamiltonian(mu_nu) > nu_full) then
                            exit
                        end if
                    end do
                    i_starts(i_cell_nu) = mu_nu
                end do
            end do
        end if
    end select
end subroutine

subroutine get_max_dimensions( &
    n_max_compute_atoms, n_max_compute_dens, n_max_compute_fns_dens, &
    shift, strain &
)
    use dimensions, only: &
        n_centers_integrals, n_centers_basis_I, n_my_batches, &
        n_centers, n_max_batch_size, n_basis_fns
    use mpi_tasks, only: check_allocation
    use pbc_lists, only: &
        inv_centers_basis_integrals, kweight_occs, de_kweight_occs

    integer, intent(out) :: n_max_compute_atoms
    integer, intent(out) :: n_max_compute_dens
    integer, intent(out) :: n_max_compute_fns_dens
    real(8), intent(in), optional :: shift(n_atoms, 3)
    real(8), intent(in), optional :: strain(3, 3)

    integer :: i_batch, i_point, i_center, i_center_integrals
    integer :: n_pruned, n_dummy, n_compute_fns, n_compute_atoms
    integer :: i_basis_fns(n_basis_fns*n_centers_integrals), &
        i_basis_fns_inv(n_basis_fns, n_centers), &
        i_atom_fns(n_basis_fns*n_centers_integrals), &
        atom_index(n_centers_integrals), &
        atom_index_inv(n_centers), &
        spline_array_start(n_centers_integrals), &
        spline_array_end(n_centers_integrals), &
        i_basis(n_centers_basis_I)
    real(8) :: r(3), r_diff(3), &
        dist_tab(n_centers_integrals), &
        dist_tab_sq(n_centers_integrals, n_max_batch_size), &
        dir_tab(3, n_centers_integrals, n_max_batch_size), &
        r_grid_min_sq(n_species)
    type(grid_point) :: point

    external :: map_to_center_cell, tab_atom_centered_coords_p0, prune_basis_p0, &
        prune_radial_basis_p0, tab_local_geometry_p0, tab_trigonom_p0, &
        tab_wave_ylm_p0, evaluate_radial_functions_p0, evaluate_waves_p0, &
        tab_gradient_ylm_p0, evaluate_wave_gradient_p0

    n_max_compute_dens = 0
    n_max_compute_fns_dens = 0
    n_max_compute_atoms = 0

    r_grid_min_sq(:) = r_grid_min(:)**2

    do i_batch = 1, n_my_batches
        n_pruned = 0
        i_basis(:) = 0
        do i_point = 1, batches(i_batch)%size
            point = batches(i_batch)%points(i_point)
            r = point%coords
            if (n_periodic > 0) call map_to_center_cell(r)
            do i_center_integrals = 1, n_centers_integrals
                i_center = centers_basis_integrals(i_center_integrals)
                r_diff(:) = r(:)-coords_center(:, i_center)
                if (present(shift)) then
                    r_diff(:) = r_diff(:)-shift(center_to_atom(i_center), :)
                else if (present(strain)) then
                    r_diff(:) = matmul(strain, r_diff)
                end if
                dir_tab(:, i_center_integrals, i_point) = r_diff(:)
                dist_tab_sq(i_center_integrals, i_point) = sum(r_diff**2)
            end do
            call prune_basis_p0( &
                dist_tab_sq(:, i_point), &
                n_dummy, n_pruned, i_basis, &
                n_centers_basis_I, n_centers_integrals, &
                inv_centers_basis_integrals &
            )
        end do
        n_max_compute_dens = max(n_pruned, n_max_compute_dens)
        if (n_pruned == 0) cycle
        do i_point = 1, batches(i_batch)%size
            n_compute_atoms = 0
            n_compute_fns = 0
            call prune_radial_basis_p0( &
                dist_tab_sq(:, i_point), dist_tab, dir_tab(:, :, i_point), &
                n_compute_atoms, atom_index, atom_index_inv, &
                n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                i_atom_fns, spline_array_start, spline_array_end, &
                n_centers_integrals, centers_basis_integrals)
            n_max_compute_fns_dens = max(n_compute_fns, n_max_compute_fns_dens)
            n_max_compute_atoms = max(n_compute_atoms, n_max_compute_atoms)
        end do
    end do
end subroutine

subroutine numerical_hirsh_pulay_forces( &
        KS_eigenvector, KS_eigenvector_cmplx, occ_numbers, density_matrix &
    )
    use dimensions, only: &
        n_basis, n_states, n_spin, n_k_points_task, n_k_points, &
        n_hamiltonian_matrix_size

    real(8), intent(in) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex(8), intent(in) :: KS_eigenvector_cmplx(n_basis, n_states, n_spin, n_k_points_task)
    real(8), intent(inout) :: occ_numbers(n_states, n_spin, n_k_points)
    real(8), intent(out) :: density_matrix(n_hamiltonian_matrix_size, n_spin)

    real(8) :: pulay_forces(n_atoms, n_atoms, 3)
    real(8) :: pulay_stress(n_atoms, 3, 3)
    real(8) :: pulay_work(n_atoms, -2:2)
    real(8) :: shift(n_atoms, 3), strain(3, 3)
    real(8) :: delta = 1.d-4
    real(8) :: h_inv(3, 3)

    integer :: j_atom, i_xyz, i_step, i_atom, j_xyz, i

    do j_atom = 1, n_atoms
        do i_xyz = 1, 3
            shift(:, :) = 0.d0
            do i_step = -2, 2
                shift(j_atom, i_xyz) = i_step*delta
                call eval_hirsh_vol_pulay_deriv_dm( &
                    KS_eigenvector, KS_eigenvector_cmplx, occ_numbers, &
                    density_matrix, shift=shift &
                )
                pulay_work(:, i_step) = sum(hirshfeld_pulay_forces(:, :, 1), 2)
            end do
            do i_atom = 1, n_atoms
                pulay_forces(i_atom, j_atom, i_xyz) = &
                    diff5(pulay_work(i_atom, :), delta)
            end do
        end do
    end do
    call xml_elem('numerical_pulay_forces', pulay_forces)
    if (n_periodic > 0) then
        do i_xyz = 1, 3
            do j_xyz = 1, 3
                do i_step = -2, 2
                    strain(:, :) = 0.d0
                    forall (i = 1:3) strain(i, i) = 1.d0
                    strain(i_xyz, j_xyz) = strain(i_xyz, j_xyz) + i_step*delta
                    call eval_hirsh_vol_pulay_deriv_dm( &
                        KS_eigenvector, KS_eigenvector_cmplx, occ_numbers, &
                        density_matrix, strain=strain &
                    )
                    pulay_work(:, i_step) = sum(hirshfeld_pulay_forces(:, :, 1), 2)
                end do
                do i_atom = 1, n_atoms
                    pulay_stress(i_atom, i_xyz, j_xyz) = &
                        diff5(pulay_work(i_atom, :), delta)
                end do
            end do
        end do
        h_inv = inverted3x3(lattice_vector)
        forall (i_atom = 1:n_atoms)
            pulay_stress(i_atom, :, :) = matmul(h_inv, pulay_stress(i_atom, :, :))
        end forall
        call xml_elem('numerical_pulay_stress', pulay_stress)
    end if
end subroutine

real(8) function diff5(y, delta)
    real(8), intent(in) :: y(5)
    real(8), intent(in) :: delta

    real(8) :: coeffs(5) = (/ 1d0/12, -2d0/3, 0d0, 2d0/3, -1d0/12 /)

    diff5 = sum(coeffs*y)/delta
end function

real(8) function diff7(y, delta)
    real(8), intent(in) :: y(7)
    real(8), intent(in) :: delta

    real(8) :: coeffs(7) = (/ -1d0/60, 3d0/20, -3d0/4, 0d0, 3d0/4, -3d0/20, 1d0/60 /)

    diff7 = sum(coeffs*y)/delta
end function

function inverted3x3(A) result(B)
    real(8), intent(in) :: A(3, 3)
    real(8), dimension(3,3) :: B

    real(8) :: D

    D = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)- &
        A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
        A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
    B(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/D
    B(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/D
    B(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/D
    B(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/D
    B(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/D
    B(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/D
    B(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/D
    B(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/D
    B(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/D
end function inverted3x3

end module
