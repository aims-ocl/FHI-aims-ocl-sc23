

! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.


module mbd
!! High-level Fortran API.

use mbd_constants
use mbd_damping, only: damping_t
use mbd_formulas, only: scale_with_ratio
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_request_t
use mbd_methods, only: get_mbd_energy, get_mbd_scs_energy
use mbd_ts, only: ts_energy
use mbd_utils, only: result_t, exception_t, printer
use mbd_vdw_param, only: ts_vdw_params, tssurf_vdw_params, species_index

implicit none

private

type, public :: mbd_input_t
    !! Contains user input to an MBD calculation.
    character(len=30) :: method = 'mbd-rsscs'
        !! VdW method to use to calculate energy and gradients.
        !!
        !! - `mbd-rsscs`: The MBD@rsSCS method.
        !! - `mbd-nl`: The MBD-NL method.
        !! - `ts`: The TS method.
        !! - `mbd`: Generic MBD method (without any screening).
    integer :: comm = -1
        !! MPI communicator.
        !!
        !! Only used when compiled with MPI. Leave as is to use the
        !! MPI_COMM_WORLD communicator.
    logical :: debug = .false.
        !! Whether debugging info should be printer
    logical :: calculate_forces = .true.
        !! Whether to calculate forces.
    logical :: calculate_spectrum = .false.
        !! Whether to keep MBD eigenvalues.
    logical :: do_rpa = .false.
        !! Whether to evalulate the MBD energy as an RPA integral over frequency.
    logical :: rpa_orders = .false.
        !! Whether to calculate individual RPA orders
    logical :: rpa_rescale_eigs = .false.
        !! Whether to rescale RPA eigenvalues as in 10.1021/acs.jctc.6b00925.
    real(dp) :: ts_ene_acc = 1d-6
        !! Required accuracy of the TS energy.
    real(dp) :: ts_f_acc = 1d-7
        !! Required accuracy of the TS gradients.
    integer :: n_omega_grid = 15
        !! Number of imaginary frequency grid points.
    real(dp) :: k_grid_shift = 0.5d0
        !! Off-\(\Gamma\) shift of the \(k\)-point grid in units of
        !! inter-\(k\)-point distance.
    logical :: zero_negative_eigvals = .false.
        !! Whether to zero out negative eigenvalues.
    character(len=20) :: xc = ''
        !! XC functional for automatic setting of damping parameters.
    real(dp) :: ts_d = 20d0
        !! TS damping parameter \(d\).
    real(dp) :: ts_sr = -1
        !! Custom TS damping parameter \(s_R\).
        !!
        !! Leave as is to use a value based on the XC functional.
    real(dp) :: mbd_a = 6d0
        !! MBD damping parameter \(a\).
    real(dp) :: mbd_beta = -1
        !! Custom MBD damping parameter \(\beta\).
        !!
        !! Leave as is to use a value based on the XC functional.
    character(len=10) :: vdw_params_kind = 'ts'
        !! Which free-atom reference vdW parameters to use for scaling.
        !!
        !! - `ts`: Values from original TS method.
        !! - `tssurf`: Values from the TS\(^\text{surf}\) approach.
    character(len=3), allocatable :: atom_types(:)
        !! (\(N\)) Atom types used for picking free-atom reference values.
    real(dp), allocatable :: free_values(:, :)
        !! (\(N\times3\), a.u.) Custom free-atom vdW paramters to use for
        !! scaling.
        !!
        !! Columns contain static polarizabilities, C6 coefficients, and vdW
        !! radii.
    real(dp), allocatable :: coords(:, :)
        !! (\(3\times N\), a.u.) Atomic coordinates.
    real(dp), allocatable :: lattice_vectors(:, :)
        !! (\(3\times 3\), a.u.) Lattice vectors in columns, unallocated if not
        !! periodic.
    integer :: k_grid(3) = [-1, -1, -1]
        !! Number of \(k\)-points along reciprocal axes.
    character(len=10) :: parallel_mode = 'auto'
        !! Parallelization scheme.
        !!
        !! - `auto`: Pick based on system system size and number of \(k\)-points.
        !! - `kpoints`: Parallelize over \(k\)-points.
        !! - `atoms`: Parallelize over atom pairs.
end type

type, public :: mbd_calc_t
    !! Represents an MBD calculation.
    private
    type(geom_t) :: geom
    type(damping_t) :: damp
    real(dp), allocatable :: alpha_0(:)
    real(dp), allocatable :: C6(:)
    character(len=30) :: method
    type(result_t) :: results
    logical :: calculate_gradients
    real(dp), allocatable :: free_values(:, :)
    logical :: debug
contains
    procedure :: init => mbd_calc_init
    procedure :: destroy => mbd_calc_destroy
    procedure :: update_coords => mbd_calc_update_coords
    procedure :: update_lattice_vectors => mbd_calc_update_lattice_vectors
    procedure :: update_vdw_params_custom => mbd_calc_update_vdw_params_custom
    procedure :: update_vdw_params_from_ratios => mbd_calc_update_vdw_params_from_ratios
    procedure :: update_vdw_params_nl => mbd_calc_update_vdw_params_nl
    procedure :: evaluate_vdw_method => mbd_calc_evaluate_vdw_method
    procedure :: get_gradients => mbd_calc_get_gradients
    procedure :: get_lattice_derivs => mbd_calc_get_lattice_derivs
    procedure :: get_lattice_stress => mbd_calc_get_lattice_stress
    procedure :: get_spectrum_modes => mbd_calc_get_spectrum_modes
    procedure :: get_rpa_orders => mbd_calc_get_rpa_orders
    procedure :: get_exception => mbd_calc_get_exception
end type

contains

subroutine mbd_calc_init(this, input)
    !! Initialize an MBD calculation from an MBD input.
    class(mbd_calc_t), target, intent(inout) :: this
    type(mbd_input_t), intent(in) :: input
        !! MBD input.


    if (input%comm /= -1) this%geom%mpi_comm = input%comm

    this%method = input%method
    this%calculate_gradients = input%calculate_forces
    this%geom%get_eigs = input%calculate_spectrum
    this%geom%get_modes = input%calculate_spectrum
    this%geom%do_rpa = input%do_rpa
    this%geom%get_rpa_orders = input%rpa_orders
    this%geom%param%rpa_rescale_eigs = input%rpa_rescale_eigs
    this%geom%param%ts_energy_accuracy = input%ts_ene_acc
    ! TODO ... = input%ts_f_acc
    this%geom%param%n_freq = input%n_omega_grid
    this%geom%param%k_grid_shift = input%k_grid_shift
    this%geom%param%zero_negative_eigvals = input%zero_negative_eigvals
    if (.not. all(input%k_grid == -1)) this%geom%k_grid = input%k_grid
    this%geom%coords = input%coords
    if (allocated(input%lattice_vectors)) then
        if (.not. allocated(this%geom%k_grid)) then
            this%geom%exc = exception_t( &
                MBD_EXC_INPUT, &
                'calc%init()', &
                'Lattice vectors present but no k-grid specified' &
            )
            return
        end if
        this%geom%lattice = input%lattice_vectors
    end if
    this%geom%parallel_mode = input%parallel_mode
    call this%geom%init()
    if (allocated(input%free_values)) then
        this%free_values = input%free_values
    else
        select case (input%vdw_params_kind)
        case ('ts')
            this%free_values = ts_vdw_params(:, species_index(input%atom_types))
        case ('tssurf')
            this%free_values = tssurf_vdw_params(:, species_index(input%atom_types))
        end select
    end if
    if (input%xc == '') then
        this%damp%beta = input%mbd_beta
        this%damp%a = input%mbd_a
        this%damp%ts_d = input%ts_d
        this%damp%ts_sr = input%ts_sr
        select case (input%method)
        case ('ts')
            if (input%ts_sr < 0) then
                this%geom%exc%code = MBD_EXC_DAMPING
                this%geom%exc%msg = 'Damping parameter S_r for TS not specified'
            end if
        case default
            if (input%mbd_beta < 0) then
                this%geom%exc%code = MBD_EXC_DAMPING
                this%geom%exc%msg = 'Damping parameter beta for MBD not specified'
            end if
        end select
    else
        this%geom%exc = this%damp%set_params_from_xc(input%xc, input%method)
    end if
    if (this%geom%has_exc()) return
    this%debug = input%debug
end subroutine

subroutine mbd_calc_destroy(this)
    !! Finalize an MBD calculation.
    class(mbd_calc_t), target, intent(inout) :: this

    call this%geom%destroy()
end subroutine

subroutine mbd_calc_update_coords(this, coords)
    !! Update atomic coordinates.
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(in) :: coords(:, :)
        !! (\(3\times N\), a.u.) New atomic coordinates.

    this%geom%coords = coords
end subroutine

subroutine mbd_calc_update_lattice_vectors(this, latt_vecs)
    !! Update unit-cell lattice vectors.
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(in) :: latt_vecs(:, :)
        !! (\(3\times 3\), a.u.) New lattice vectors in columns.

    this%geom%lattice = latt_vecs
end subroutine

subroutine mbd_calc_update_vdw_params_custom(this, alpha_0, C6, r_vdw)
    !! Update vdW parameters in a custom way.
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(in) :: alpha_0(:)
        !! (a.u.) New atomic static polarizabilities.
    real(dp), intent(in) :: C6(:)
        !! (a.u.) New atomic \(C_6\) coefficients.
    real(dp), intent(in) :: r_vdw(:)
        !! (a.u.) New atomic vdW radii.

    this%alpha_0 = alpha_0
    this%C6 = C6
    this%damp%r_vdw = r_vdw
end subroutine

subroutine mbd_calc_update_vdw_params_from_ratios(this, ratios)
    !! Update vdW parameters based on scaling of free-atom values.
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(in) :: ratios(:)
        !! Ratios of atomic volumes in the system and in vacuum.

    real(dp), allocatable :: ones(:)

    allocate (ones(size(ratios)), source=1d0)
    this%alpha_0 = scale_with_ratio(this%free_values(1, :), ratios, ones, 1d0)
    this%C6 = scale_with_ratio(this%free_values(2, :), ratios, ones, 2d0)
    this%damp%r_vdw = scale_with_ratio(this%free_values(3, :), ratios, ones, 1d0/3)
end subroutine

subroutine mbd_calc_update_vdw_params_nl(this, alpha_0_ratios, C6_ratios)
    !! Update vdW parameters for the MBD-NL method.
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(in) :: alpha_0_ratios(:)
        !! Ratios of free-atom exact static polarizabilities and those from the
        !! VV functional.
    real(dp), intent(in) :: C6_ratios(:)
        !! Ratios of free-atom exact \(C_6\) coefficients and those from the VV
        !! functional.

    this%alpha_0 = this%free_values(1, :)*alpha_0_ratios
    this%C6 = this%free_values(2, :)*C6_ratios
    this%damp%r_vdw = 2.5d0*this%free_values(1, :)**(1d0/7)*alpha_0_ratios**(1d0/3)
end subroutine

subroutine mbd_calc_evaluate_vdw_method(this, energy)
    !! Evaluate a given vdW method for a given system and vdW parameters,
    !! retrieve energy.
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(out) :: energy
        !! (a.u.) VdW energy.

    type(grad_request_t) :: grad

    if (this%calculate_gradients) then
        grad%dcoords = .true.
        if (allocated(this%geom%lattice)) grad%dlattice = .true.
    end if
    select case (this%method)
    case ('mbd', 'mbd-nl')
        this%damp%version = 'fermi,dip'
        this%results = get_mbd_energy( &
            this%geom, this%alpha_0, this%C6, this%damp, grad &
        )
        energy = this%results%energy
    case ('mbd-rsscs')
        this%results = get_mbd_scs_energy( &
            this%geom, 'rsscs', this%alpha_0, this%C6, this%damp, grad &
        )
        energy = this%results%energy
    case ('ts')
        energy = ts_energy(this%geom, this%alpha_0, this%C6, this%damp)
    end select
    if (this%debug) call this%geom%timer%print()
end subroutine

subroutine mbd_calc_get_gradients(this, gradients)  ! 3 by N  dE/dR
    !! Retrieve nuclear energy gradients if they were requested in the MBD
    !! input.
    !!
    !! The gradients are calculated together with the energy, so a call to this
    !! method must be preceeded by a call to
    !! [[mbd_calc_t:evaluate_vdw_method]].  For the same reason, the
    !! gradients must be requested prior to this called via
    !! [[mbd_input_t:calculate_forces]].
    class(mbd_calc_t), intent(in) :: this
    real(dp), intent(out) :: gradients(:, :)
        !! (\(3\times N\), a.u.) Energy gradients, \(\mathrm dE/\mathrm d\mathbf
        !! R_i\), index \(i\) runs over columns.

    gradients = transpose(this%results%dE%dcoords)
end subroutine

subroutine mbd_calc_get_lattice_derivs(this, latt_derivs)
    !! Provide lattice-vector energy gradients if they were requested in the MBD
    !! input.
    !!
    !! The gradients are actually calculated together with the energy, so a call
    !! to this method must be preceeded by a call to
    !! [[mbd_calc_t:evaluate_vdw_method]].  For the same reason, the
    !! gradients must be requested prior to this called via
    !! [[mbd_input_t:calculate_forces]].
    class(mbd_calc_t), intent(in) :: this
    real(dp), intent(out) :: latt_derivs(:, :)
        !! (\(3\times 3\), a.u.) Energy gradients, \(\mathrm dE/\mathrm d\mathbf
        !! a_i\), index \(i\) runs over columns.

    latt_derivs = transpose(this%results%dE%dlattice)
end subroutine

subroutine mbd_calc_get_lattice_stress(this, stress)
    !! Provide stress tensor of the lattice.
    !!
    !! This is a utility function wrapping [[mbd_calc_t:get_lattice_derivs]].
    !! The lattice vector gradients are coverted to the stress tensor.
    class(mbd_calc_t), intent(in) :: this
    real(dp), intent(out) :: stress(:, :)
        !! (\(3\times 3\), a.u.) Stress tensor.

    stress = ( &
        matmul(this%geom%lattice, this%results%dE%dlattice) &
        + matmul(this%geom%coords, this%results%dE%dcoords) &
    )
end subroutine

subroutine mbd_calc_get_spectrum_modes(this, spectrum, modes)
    !! Provide MBD spectrum if it was requested in the MBD input.
    !!
    !! The spectrum is actually calculated together with the energy, so a call
    !! to this method must be preceeded by a call to
    !! [[mbd_calc_t:evaluate_vdw_method]].  For the same reason, the
    !! spectrum must be requested prior to this call via
    !! [[mbd_input_t:calculate_spectrum]].
    class(mbd_calc_t), intent(inout) :: this
    real(dp), intent(out) :: spectrum(:)
        !! (\(3N\), a.u.) Energies (frequencies) of coupled MBD modues,
        !! \(\omega_i\).
    real(dp), intent(out), allocatable, optional :: modes(:, :)
        !! (\(3N\times 3N\)) Coupled-mode wave functions (MBD eigenstates),
        !! \(\psi_j\), in the basis of uncoupled states,
        !! \(C_{ij}=\langle\phi_i|\psi_j\rangle\), index \(j\) runs over
        !! columns.
        !!
        !! To save memory, the argument must be allocatable, and the method
        !! transfers allocation from the internal state to the argument. For
        !! this reason, the method can be called only once wih this optional
        !! argument per calculation.

    spectrum = this%results%mode_eigs
    if (present(modes)) call move_alloc(this%results%modes, modes)
end subroutine

subroutine mbd_calc_get_rpa_orders(this, rpa_orders)
    !! Provide RPA orders if they were requested in the MBD input.
    !!
    !! The orders are actually calculated together with the energy, so a call
    !! to this method must be preceeded by a call to
    !! [[mbd_calc_t:evaluate_vdw_method]]. For the same reason, the
    !! spectrum must be requested prior to this call via
    !! [[mbd_input_t:do_rpa]] and [[mbd_input_t:rpa_orders]].
    class(mbd_calc_t), intent(inout) :: this
    real(dp), allocatable, intent(out) :: rpa_orders(:)
        !! (a.u.) MBD energy decomposed to RPA orders.

    rpa_orders = this%results%rpa_orders
end subroutine

subroutine mbd_calc_get_exception(this, code, origin, msg)
    !! Retrieve an exception in the MBD calculation if it occured.
    class(mbd_calc_t), intent(inout) :: this
    integer, intent(out) :: code
        !! Exception code, values defined in [[mbd_constants]].
    character(*), intent(out) :: origin
        !! Exception origin.
    character(*), intent(out) :: msg
        !! Exception message.

    code = this%geom%exc%code
    if (code == 0) return
    origin = this%geom%exc%origin
    msg = this%geom%exc%msg
    this%geom%exc%code = 0
    this%geom%exc%origin = ''
    this%geom%exc%msg = ''
end subroutine

end module
