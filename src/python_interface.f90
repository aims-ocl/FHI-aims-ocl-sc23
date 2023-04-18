!****h* FHI-aims/python_interface
! PURPOSE
!   Implements support for calling the embedded Python.
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2016-04-24
! NOTES
!   Requires compiler support for the iso_c_binding module (Fortran 2003).
!   `register_python_hook()` is called by `read_control`. A hook is invoked by
!   calling `run_python_hook()`. See cffi/python_interface.py for the Python
!   side of the interface.
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note
!   that any use of the "FHI-aims-Software" is subject to the terms and
!   conditions of the respective license agreement.
!******
module python_interface

use iso_c_binding
use localorb_io, only: localorb_info, default_unit

implicit none

private

public :: register_python_hook, run_python_hook

integer, parameter :: hookname_len = 50

type :: Hook_t
    private
    character(len=hookname_len) :: filename
    character(len=hookname_len) :: label
    logical :: parallel = .false.
    logical, public :: registered = .false.
end type

type :: HookRegister_t
    type(Hook_t) :: &
        post_scf, post_hirshfeld
end type

type(HookRegister_t), public :: python_hooks

! Wraps the Fortran `grid_point` type to C
type, bind(c) :: GridPoint_t
    type(c_ptr) :: coords
    integer(c_int) :: index_atom, index_radial, index_angular
end type

! Wraps the Fortran `batch_of_points` type to C
type, bind(c) :: BatchOfPoints_t
    integer(c_int) :: size
    type(c_ptr) :: points
end type

! C struct that is passed to the Python part of the interface
type, bind(c) :: AimsContext_t
    logical(c_bool) :: real_eigenvectors
    type(c_ptr) :: elements
    integer(c_int) :: &
        n_atoms, n_full_points, n_spin, n_my_batches, n_k_points, &
        n_k_points_task, n_basis, n_states, n_hamiltonian_matrix_size
    type(c_ptr) :: &
        species
    type(c_ptr) :: &
        coords, rho, rho_gradient, kinetic_density, &
        partition_tab, hirshfeld_volume, &
        hirshfeldw, freeintegral, &
        KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, &
        hamiltonian
    type(c_ptr) :: batches
end type

! Fortran interface for the C signature of the corresponding Python function
interface
    integer(c_int) function call_python_cffi(ctx, filename, event, my_task) &
            result(retcode) bind(c)
        import AimsContext_t, c_char, c_int
        type(AimsContext_t), intent(inout) :: ctx
        character(kind=c_char), intent(in) :: filename(*)
        character(kind=c_char), intent(in) :: event(*)
        integer(c_int), intent(in), value :: my_task
    end function
end interface

contains

!****f* python_interface/register_python_hook
! PURPOSE
!   Register a Python hook to a specified location.
! INPUTS
!   * label (in) string -- Hook location.
!   * filename (in) string -- Path to Python source file.
!   * attrib (in) string -- Hook attribute (may be empty string).
!******
integer function register_python_hook(label, filename, attrib) result(retcode)
    use mpi_tasks

    character(len=*), intent(in) :: label, filename, attrib

    type(AimsContext_t) :: ctx
    character(len=hookname_len+1) :: c_filename
    type(Hook_t) :: hook

    retcode = 0
    hook%registered = .true.
    if (filename == 'REPL' .and. n_tasks > 1) then
        call localorb_info(' *** Error: Cannot execute REPL hook in an MPI run')
        return
    end if
    hook%filename = filename
    hook%label = label
    select case (attrib)
    case ('')
    case ('parallel')
        hook%parallel = .true.
    case default
        call localorb_info(' *** Error: Unknown hook attribute: ' // attrib)
        return
    end select
    select case (label)
    case ('post_scf')
        python_hooks%post_scf = hook
    case ('post_hirshfeld')
        python_hooks%post_hirshfeld = hook
    case default
        call localorb_info(' *** Error: Unknown hook label: ' // label)
        return
    end select
    c_filename = trim(hook%filename) // c_null_char
    ! In the parse call, the context object is not used and can be uninitialized.
    retcode = call_python_cffi(ctx, c_filename, 'parse' // c_null_char, myid)
    if (retcode /= -1) then
        call aims_stop()
    end if
    call localorb_info('  Registered python hook ' // trim(label))
end function

!****f* python_interface/run_python_hook
! PURPOSE
!   Executes a Python hook.
! INPUTS
!   * hook (in) PythonHook_t -- Python hook to be executed.
!******
subroutine run_python_hook(hook)
    use mpi_tasks
    use xml_write, only: tostr

    type(Hook_t), intent(in) :: hook

    external :: MPI_Bcast

    character(len=hookname_len+1) :: c_filename
    type(AimsContext_t) :: ctx
    integer :: retcode, mpicode, rootcode

    if (.not. hook%registered) then
        call localorb_info(' *** Error: Python hook not registered')
        call aims_stop()
    end if
    call localorb_info( &
        '  Executing Python hook ' // trim(hook%label) // ' (' &
        // trim(hook%filename) // ')...' &
    )
    if (hook%parallel .or. myid == 0) then
        ctx = get_aims_context()
        c_filename = trim(hook%filename) // c_null_char
        retcode = call_python_cffi(ctx, c_filename, 'run' // c_null_char, myid)
        call destroy_aims_context(ctx)
    else
        retcode = -1
    end if
    rootcode = retcode
    if (use_mpi) then  ! send root return code to all tasks
        call MPI_Bcast(rootcode, 1, MPI_INTEGER, 0, mpi_comm_global, mpicode)
    end if
    if (rootcode /= -1) then  ! no need to print anything from Fortran when root fails
        call aims_stop()
    else if (retcode /= -1) then
        write (default_unit, *) &
            ' *** Error: Python hook return code on process ' &
            // trim(tostr(myid)) // ': ' // trim(tostr(retcode))
        call aims_stop()
    end if
    call localorb_info('  Successfully executed Python hook ' // trim(hook%label))
end subroutine

! Wrap Fortran objects as C objects.
type(AimsContext_t) function get_aims_context() result(ctx)
    use runtime_choices, only: real_eigenvectors
    use species_data, only: species_element
    use dimensions, only: n_atoms, n_full_points, n_spin, n_my_batches, &
        n_k_points, n_k_points_task, n_basis, n_states, n_hamiltonian_matrix_size
    use physics, only: &
        rho, rho_gradient, kinetic_density, partition_tab, &
        KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, &
        hamiltonian
    use geometry, only: coords, species
    use vdw_correction, only: hirshfeld_volume, hirshfeldw, freeintegral
    use grids, only: batches

    integer :: i_batch, i_point
    type(BatchOfPoints_t), pointer :: c_batches(:)
    type(GridPoint_t), pointer :: c_points(:)

    ! The fortran arrays need to have `target` attribute to obtain C pointers,
    ! we do this by wrapping them as dummy arguments in a subroutine.
    ctx = target_wrapper( &
        coords, species, rho, rho_gradient, kinetic_density, &
        partition_tab, hirshfeld_volume, species_element, &
        KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, &
        hamiltonian, hirshfeldw, freeintegral &
    )
    ctx%n_atoms = n_atoms
    ctx%n_hamiltonian_matrix_size = n_hamiltonian_matrix_size
    ctx%n_full_points = n_full_points
    ctx%n_spin = n_spin
    ctx%n_my_batches = n_my_batches
    ctx%n_k_points = n_k_points
    ctx%n_k_points_task = n_k_points_task
    ctx%n_basis = n_basis
    ctx%n_states = n_states
    ctx%real_eigenvectors = real_eigenvectors

    ! Allocate memory for the bathches->points tree. This needs to be
    ! deallocated manually when the python hook finishes.
    allocate (c_batches(n_my_batches))
    ctx%batches = c_loc(c_batches)
    do i_batch = 1, n_my_batches
        allocate (c_points(batches(i_batch)%size))
        c_batches(i_batch)%size = batches(i_batch)%size
        c_batches(i_batch)%points = c_loc(c_points)
        do i_point = 1, batches(i_batch)%size
            c_points(i_point)%coords = c_loc(batches(i_batch)%points(i_point)%coords)
            c_points(i_point)%index_atom = batches(i_batch)%points(i_point)%index_atom
            c_points(i_point)%index_radial = &
                batches(i_batch)%points(i_point)%index_radial
            c_points(i_point)%index_angular = &
                batches(i_batch)%points(i_point)%index_angular
        end do
    end do

    contains

    type(AimsContext_t) function target_wrapper( &
        coords, species, rho, rho_gradient, kinetic_density, &
        partition_tab, hirshfeld_volume, elements, &
        KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, &
        hamiltonian, hirshfeldw, freeintegral &
    ) result(ctx)
        real(c_double), target :: &
            coords(:, :), rho(:, :), rho_gradient(:, :, :), &
            kinetic_density(:, :), partition_tab(:), hirshfeld_volume(:), &
            KS_eigenvector(:, :, :, :), KS_eigenvalue(:, :, :), &
            occ_numbers(:, :, :), hamiltonian(:, :), hirshfeldw(:), freeintegral(:)
        integer(c_int), target :: &
            species(:)
        complex(c_double_complex), target :: &
            KS_eigenvector_complex(:, :, :, :)
        character(c_char), target :: &
            elements(2, n_atoms)
        ctx%coords = c_loc(coords)
        ctx%species = c_loc(species)
        ctx%rho = c_loc(rho)
        ctx%rho_gradient = c_loc(rho_gradient)
        ctx%kinetic_density = c_loc(kinetic_density)
        ctx%partition_tab = c_loc(partition_tab)
        ctx%hirshfeld_volume = c_loc(hirshfeld_volume)
        ctx%hirshfeldw = c_loc(hirshfeldw)
        ctx%freeintegral = c_loc(freeintegral)
        ctx%KS_eigenvector = c_loc(KS_eigenvector)
        ctx%KS_eigenvector_complex = c_loc(KS_eigenvector_complex)
        ctx%KS_eigenvalue = c_loc(KS_eigenvalue)
        ctx%occ_numbers = c_loc(occ_numbers)
        ctx%hamiltonian = c_loc(hamiltonian)
        ctx%elements = c_loc(elements)
    end function
end function

! Deallocate all allocated pointers.
subroutine destroy_aims_context(ctx)
    type(AimsContext_t), intent(in) :: ctx

    type(BatchOfPoints_t), pointer :: c_batches(:)
    type(GridPoint_t), pointer :: c_points(:)
    integer :: i_batch

    call c_f_pointer(ctx%batches, c_batches, [ctx%n_my_batches])
    do i_batch = 1, ctx%n_my_batches
        call c_f_pointer( &
            c_batches(i_batch)%points, &
            c_points, &
            [c_batches(i_batch)%size] &
        )
        deallocate (c_points)
    end do
    deallocate (c_batches)
end subroutine

end module
