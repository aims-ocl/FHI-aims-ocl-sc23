!****h* FHI-aims/json_output
!*  NAME
!*    json_output
!*  SYNOPSIS
module json_output
!*  PURPOSE
!*    This module is responsible for outputting JSON texts using FortJSON.
!*  USES
  use FortJSON, only: fjson_handle
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University)
!*  HISTORY
!*    January 2018 - Created.
!*    April 2018 - Converted to use FortJSON.
!*    May 2018 - Flattened to make more amenable to the data framed data structure.
!*    June 2018 - Refactoring to eliminate hidden dependencies introduced by
!*      changes to ELSI uild system.
!*  TODO
!*    Known missing:
!*    - Forces and stress tensor
!*    - Band structure record object
!*    Known issues:
!*    - When SOC is enabled, the Mulliken object will be output twice, once for
!*      the SR version and once for the SOC version.  This is desirable
!*      behavior, but there's currently no way to distinguish the two.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is
!*    subject to the terms and conditions of the respective license
!*    agreement.
!*  SOURCE

  private

  type(fjson_handle), private :: json_log_handle ! IO handle for JSON log

  integer, parameter          :: JSON_LOG_UNIT = 67
  character(len=*), parameter :: JSON_LOG_NAME = "aims.json"

  ! File Management
  public :: open_json_log_file
  public :: close_json_log_file
  public :: get_json_handle
  ! Helper subroutines to write blocks of data to output
  public :: write_core_aims_vars_to_json_handle
  public :: write_versioning_to_json_handle
  public :: write_gpu_info_to_json_handle
  public :: write_symmetry_to_json_handle
  public :: write_lattice_to_json_handle
  public :: write_atoms_to_json_handle
  public :: write_elec_struct_to_json_handle
  public :: write_stats_to_json_handle
  ! Writing (well-formed) record objects
  public :: write_geometry_to_json_log
  public :: write_scf_iteration_to_json_log
  public :: write_mulliken_to_json_log
  public :: write_final_output_to_json_log

  interface write_stats_to_json_handle
    module procedure write_stats_int_to_json_handle
  end interface

contains

  pure type(fjson_handle) function get_json_handle() result(y)
    y = json_log_handle
  end function get_json_handle

  !****f* json_output/open_json_log_file
  !*  NAME
  !*    open_json_log_file
  !*  SYNOPSIS
  subroutine open_json_log_file()
  !*  PURPOSE
  !*    Wrapper around fjson_open_file to open the file where the JSON log will
  !*    be written
  !*  USES
    use FortJSON, only: fjson_open_file, fjson_start_array
    use mpi_tasks, only: myid
    use runtime_choices, only: out_aims_json_log
    implicit none
  !*  ARGUMENTS
  !*    o None (modifies module variables)
  !*  INPUTS
  !*    o None
  !*  OUTPUTS
  !*    o None (modifies module variables)
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    if (myid.ne.0 .or. .not.out_aims_json_log) return

    call fjson_open_file(json_log_handle, JSON_LOG_UNIT, JSON_LOG_NAME)
    call fjson_start_array(json_log_handle)
  end subroutine open_json_log_file
  !******

  !****f* json_output/close_json_log_file
  !*  NAME
  !*    close_json_log_file
  !*  SYNOPSIS
  subroutine close_json_log_file()
  !*  PURPOSE
  !*    Wrapper around fjson_close_file to close the file where the JSON log
  !*    is being written
  !*  USES
    use FortJSON, only: fjson_close_file, fjson_finish_array
    use mpi_tasks, only: myid
    use runtime_choices, only: out_aims_json_log
    implicit none
  !*  ARGUMENTS
  !*    o None (modifies module variables)
  !*  INPUTS
  !*    o None
  !*  OUTPUTS
  !*    o None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    if (myid.ne.0 .or. .not.out_aims_json_log) return

    call fjson_finish_array(json_log_handle)
    call fjson_close_file(json_log_handle)
  end subroutine close_json_log_file
  !******

  subroutine synchronize_batch_statistics( )
    use batch_statistics, only: n_max_compute_ham_stats, ld_ham_integ_stats, &
        ld_ham_hpot_stats, ld_ham_density_stats, batch_n_compute_stats
    use dimensions, only: n_grid_batches, n_max_compute_ham, n_my_batches
    use grids, only: batches
    use load_balancing, only: batch_perm, batch_permutation, n_bp_integ, &
        n_bp_hpot, n_bp_density
    use mpi_tasks, only: myid, n_tasks
    use mpi_utilities, only: batch_task_list
    use statistics, only: calc_stats_array
    use synchronize_mpi_basic, only: sync_int_vector
    implicit none
    integer :: i_grid_batch, i_my_batch, task_for_batch, offset = 0
    integer, allocatable :: n_batches_on_task(:), array(:)

    allocate( array(n_tasks) )

    ! Calculate distribution for various quantities spread across MPI tasks
    array = 0
    array(myid+1) = n_max_compute_ham
    call sync_int_vector(array, n_tasks)
    call calc_stats_array(array, n_max_compute_ham_stats, n_fractiles = 10, &
                          n_bins = 10)

    ! Local matrix size for integration batch permutation
    if (batch_perm(n_bp_integ)%initialized) then
      array = 0
      array(myid+1) = batch_perm(n_bp_integ)%n_local_matrix_size
      call sync_int_vector(array, n_tasks)
      call calc_stats_array(array, ld_ham_integ_stats, n_fractiles = 10, &
                            n_bins = 10)
    end if

    ! Local matrix size for Hartree potential batch permutation
    if (batch_perm(n_bp_hpot)%initialized) then
      array = 0
      array(myid+1) = batch_perm(n_bp_hpot)%n_local_matrix_size
      call sync_int_vector(array, n_tasks)
      call calc_stats_array(array, ld_ham_hpot_stats, n_fractiles = 10, &
                            n_bins = 10)
    end if

    ! Local matrix size for density batch permutation
    if (batch_perm(n_bp_density)%initialized) then
      array = 0
      array(myid+1) = batch_perm(n_bp_density)%n_local_matrix_size
      call sync_int_vector(array, n_tasks)
      call calc_stats_array(array, ld_ham_density_stats, n_fractiles = 10, &
                            n_bins = 10)
    end if

    deallocate( array )

    ! The batch size statistics are computed in partition_grid(), so we don't
    ! need to recompute them here.

    ! Here's where we compute the number of basis functions touching a batch.
    ! This is a bit more difficult than previous statistics, because as far as
    ! I know, there's no good place in the code where we calculate this quantity
    ! and still have the global index for a batch, like we did for batch sizes.
    ! So we have to go through a somewhat awkward construction of extracting the
    ! number of batches per MPI task from batch_task_list, using that knowledge
    ! to determine how many batches all MPI tasks with lower ranks than the
    ! current one have, which then allows us to create a list of batches across
    ! all MPI ranks, ordered by MPI rank.
    allocate( array(n_grid_batches) )

    ! Create a list of how many batches each MPI task has, essentially the
    ! inverse array of batch_task_list
    allocate( n_batches_on_task(0:n_tasks-1) )
    n_batches_on_task = 0
    do i_grid_batch = 1, n_grid_batches
      task_for_batch = batch_task_list(i_grid_batch)
      n_batches_on_task(task_for_batch) = n_batches_on_task(task_for_batch) + 1
    end do
    ! Now use the number of batches all previous MPI tasks have to determine
    ! the current task's offset in the global list of batches
    offset = 0
    if (myid > 0) offset = sum( n_batches_on_task(0:(myid-1)) )
    deallocate( n_batches_on_task )

    ! Now that we know where the batches on this MPI task reside in the global
    ! array of batches, generate statistics on batch_n_compute
    array = 0
    do i_my_batch = 1, n_my_batches
      array(offset + i_my_batch) = batches(i_my_batch)%batch_n_compute
    end do
    call sync_int_vector(array, n_grid_batches)
    call calc_stats_array(array, batch_n_compute_stats, n_fractiles = 10, &
                          n_bins = 10)

    deallocate( array )
  end subroutine synchronize_batch_statistics

  !****f* json_output/write_core_aims_vars_to_json_handle
  !*  NAME
  !*    write_core_aims_vars_to_json_handle
  !*  SYNOPSIS
  subroutine write_core_aims_vars_to_json_handle(fj_h)
  !*  PURPOSE
  !*    Writes out a set of aims variables that are considered to be fundamental
  !*    to the inner workings of aims as a sequence of name/value piars.  These
  !*    variables can be broadly classified as common array dimensions, physical
  !*    levels of theory used, and parallelization/data structure choices.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use dimensions
    use FortJSON, only: fjson_write_name_value, fjson_start_name_object, &
        fjson_finish_object, fjson_start_name_array, fjson_write_value, &
        fjson_finish_array
    use geometry, only: cell_volume, species
    use mpi_tasks, only: aims_stop, n_tasks
    use physics
    use runtime_choices
    use species_data, only: atoms_in_structure, species_name
    use xc_library
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
  !*  INPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    character(LEN=200) :: info_str
    integer :: i

    character(*), parameter :: func = 'write_core_aims_vars_to_json_handle'

    ! Physics
    call fjson_write_name_value(fj_h, "n_atoms", n_atoms)
    call fjson_write_name_value(fj_h, "n_species", n_species)
    ! Should add the ability to print arrays of strings to FortJSON...
    call fjson_start_name_array(fj_h, "species_name")
    do i = 1, n_species
      call fjson_write_value(fj_h, trim(adjustl(species_name(i))))
    end do
    call fjson_finish_array(fj_h)
    call fjson_start_name_object(fj_h, "atoms_in_structure")
    do i = 1, n_species
      call fjson_write_name_value(fj_h, trim(adjustl(species_name(i))), atoms_in_structure(i))
    end do
    call fjson_finish_object(fj_h)
    call fjson_write_name_value(fj_h, "n_basis", n_basis)
    call fjson_write_name_value(fj_h, "n_states", n_states)
    call fjson_write_name_value(fj_h, "n_spin", n_spin)
    call fjson_write_name_value(fj_h, "n_k_points", n_k_points)
    call fjson_write_name_value(fj_h, "n_electrons", n_electrons)
    call fjson_write_name_value(fj_h, "n_k_points_task", n_k_points_task)
    select case (flag_rel)
    case (REL_none)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_none")
    case (REL_zora)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_zora")
    case (REL_atomic_zora)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_atomic_zora")
    case (REL_own)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_own")
    case (REL_kolning_harmon)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_kolning_harmon")
    case (REL_zora_spinor)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_zora_spinor")
    case (REL_at_zora_spinor)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_at_zora_spinor")
    case (REL_x2c)
      call fjson_write_name_value(fj_h, "flag_rel", "REL_x2c")
    case default
      call aims_stop("Unsupported relativistic handling", func)
    end select
    call fjson_write_name_value(fj_h, "xc", get_xc_name(flag_xc))
    if (is_xc_lda(flag_xc)) then
      call fjson_write_name_value(fj_h, "xc_type", "lda")
    else if (is_xc_gga(flag_xc)) then
      call fjson_write_name_value(fj_h, "xc_type", "gga")
    else if (is_xc_meta_gga(flag_xc)) then
      call fjson_write_name_value(fj_h, "xc_type", "meta_gga")
    else if (is_xc_hybrid(flag_xc)) then
      call fjson_write_name_value(fj_h, "xc_type", "hybrid")
    else if (is_xc_double_hybrid(flag_xc)) then
      call fjson_write_name_value(fj_h, "xc_type", "double_hybrid")
    else if (is_xc_nonlocal(flag_xc)) then
      call fjson_write_name_value(fj_h, "xc_type", "nonlocal")
    else
      call aims_stop("Unsupported xc functional", func)
    end if
    if (n_periodic.gt.0) then
      call fjson_write_name_value(fj_h, "is_periodic", .true.)
      ! k_grid is the name of the control.in keyword, which is more familiar
      call fjson_write_name_value(fj_h, "k_grid", n_k_points_xyz, 3)
      call fjson_write_name_value(fj_h, "cell_volume", cell_volume*bohr*bohr*bohr)
    else
      call fjson_write_name_value(fj_h, "is_periodic", .false.)
    end if

    ! Runtime Settings
    call fjson_write_name_value(fj_h, "real_eigenvectors", real_eigenvectors)
    call fjson_write_name_value(fj_h, "collect_eigenvectors", collect_eigenvectors)
    call fjson_write_name_value(fj_h, "use_scalapack", use_scalapack)
    call fjson_write_name_value(fj_h, "n_tasks", n_tasks)
    select case (packed_matrix_format)
    case (PM_none)
      call fjson_write_name_value(fj_h, "packed_matrix_format", "PM_none")
    case (PM_index)
      call fjson_write_name_value(fj_h, "packed_matrix_format", "PM_index")
    case default
      call aims_stop("Unsupported matrix format", func)
    end select
    call fjson_write_name_value(fj_h, "use_local_index", use_local_index)
    call fjson_write_name_value(fj_h, "use_load_balancing", use_load_balancing)
  end subroutine write_core_aims_vars_to_json_handle
  !******

  !****f* json_output/write_versioning_to_json_handle
  !*  NAME
  !*    write_versioning_to_json_handle
  !*  SYNOPSIS
  subroutine write_versioning_to_json_handle(fj_h)
  !*  PURPOSE
  !*    Writes out a set of aims variables that are useful for versioning as a
  !*    sequence of name/value pairs.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use FortJSON, only: fjson_write_name_value, fjson_get_datetime_rfc3339, &
                        DATETIME_LEN
    use generate_aims_uuid, only: write_bare_aims_uuid
    use mpi_tasks, only: myid, MPI_MAX_PROCESSOR_NAME
    use mpi_utilities, only: get_my_processor
    use runtime_choices, only: calculation_set, calculation_subset
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
  !*  INPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    character(LEN=37)  :: uuid_bare_str
    character(LEN=MPI_MAX_PROCESSOR_NAME) :: proc_name
    character(len=DATETIME_LEN) :: datetime_of_record
    integer :: i
    logical :: spglib_was_compiled

    character(*), parameter :: func = 'write_versioning_to_json_handle'

    call get_my_processor(proc_name, i) ! i passed as dummy
    call write_bare_aims_uuid(uuid_bare_str)
    call fjson_get_datetime_rfc3339(datetime_of_record)

    ! Versioning
    call fjson_write_name_value(fj_h, "data_source", "FHI-aims")
    call fjson_write_name_value(fj_h, "data_source_version", "UNKNOWN")
    call fjson_write_name_value(fj_h, "git_commit", "UNKNOWN")
    call fjson_write_name_value(fj_h, "git_commit_modified", "UNKNOWN")
    call fjson_write_name_value(fj_h, "git_message_abbrev", "UNKNOWN")
    call fjson_write_name_value(fj_h, "executable_created_on_hostname", "UNKNOWN")
    call fjson_write_name_value(fj_h, "executable_created_on_datetime", "UNKNOWN")
    call fjson_write_name_value(fj_h, "uuid", trim(uuid_bare_str))
    call fjson_write_name_value(fj_h, "reporting_task_rank", myid)
    call fjson_write_name_value(fj_h, "reporting_task_name", trim(adjustl(proc_name)))
    if (allocated(calculation_set)) then
      call fjson_write_name_value(fj_h, "calculation_set", trim(calculation_set))
    end if
    if (allocated(calculation_subset)) then
      call fjson_write_name_value(fj_h, "calculation_subset", trim(calculation_subset))
    end if
    call fjson_write_name_value(fj_h, "datetime_of_record", trim(datetime_of_record))
    call fjson_write_name_value(fj_h, "spglib_was_compiled", spglib_was_compiled)
  end subroutine write_versioning_to_json_handle
  !******

  !****f* json_output/write_gpu_info_to_json_handle
  !*  NAME
  !*    write_gpu_info_to_json_handle
  !*  SYNOPSIS
  subroutine write_gpu_info_to_json_handle(fj_h)
  !*  PURPOSE
  !*    Writes out information about the GPU (including whether it's being used)
  !*    as a sequence of name/value pairs.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use aims_gpu
    use FortJSON, only: fjson_write_name_value
    use physics, only: forces_on
    use runtime_choices, only: use_gpu
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
  !*  INPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    character(*), parameter :: func = 'write_gpu_info_to_json_handle'

    call fjson_write_name_value(json_log_handle, "use_gpu", use_gpu)
    if (use_gpu) then
      call fjson_write_name_value(json_log_handle, "gpu_density_used", gpu_density_used)
      call fjson_write_name_value(json_log_handle, "gpu_hamiltonian_used", gpu_hamiltonian_used)
      if (forces_on) then
        call fjson_write_name_value(json_log_handle, "gpu_forces_used", gpu_forces_used)
      end if
      call fjson_write_name_value(json_log_handle, "gpu_name", gpu_name)
      call fjson_write_name_value(json_log_handle, "n_gpus_avail", n_gpus_avail)
      call fjson_write_name_value(json_log_handle, "cuda_major", cuda_major)
      call fjson_write_name_value(json_log_handle, "cuda_minor", cuda_minor)
      call fjson_write_name_value(json_log_handle, "total_global_gpu_mem", total_global_gpu_mem)
      call fjson_write_name_value(json_log_handle, "total_constant_gpu_mem", total_constant_gpu_mem)
      call fjson_write_name_value(json_log_handle, "shared_gpu_mem_per_block", shared_gpu_mem_per_block)
      call fjson_write_name_value(json_log_handle, "gpu_clock_rate", gpu_clock_rate)
      call fjson_write_name_value(json_log_handle, "gpu_multiprocessor_count", gpu_multiprocessor_count)
      call fjson_write_name_value(json_log_handle, "gpu_max_threads_per_multiprocessor", gpu_max_threads_per_multiprocessor)
    endif
  end subroutine write_gpu_info_to_json_handle
  !******

  !****f* spglib_symm/write_symmetry_to_json_handle
  !*  NAME
  !*    write_symmetry_to_json_handle
  !*  SYNOPSIS
  subroutine write_symmetry_to_json_handle &
             (fj_h, lattice_name, lattice_vector, n_atoms, species, coords)
  !*  PURPOSE
  !*    Writes out symmetry information about the geometry as a sequence of
  !*    name/value pairs.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use, intrinsic :: iso_c_binding, only: c_char
    use dimensions, only: sym_precision
    use FortJSON, only: fjson_write_name_value, fjson_handle
    use spglib_symmetry, only: check_spglib, get_symmetry_for_lattice
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
    character(len=*), intent(in) :: lattice_name
    real*8 , intent(in) :: lattice_vector(3,3)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: species(n_atoms)
    real*8,  intent(in) :: coords(3, n_atoms)
  !*  INPUTS
  !*    o fj_h - FortJSON handle to write to
  !*    o lattice_name - the name of the lattice, to be output as a prefix
  !*    o lattice_vector - the lattice vectors
  !*    o n_atoms - number of atoms
  !*    o species - species for atoms
  !*    o coords - Cartesian coordinates for atoms
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    character(*), parameter :: func = 'write_symmetry_to_json_handle'

    logical :: spglib_was_compiled
    real*8  :: symprec
    logical :: spacegroup_found
    integer :: n_operations_sym
    integer :: spacegroup_number
    character(len=14, kind=c_char) :: international_symbol
    character(len=10, kind=c_char) :: schoenflies

    ! When spglib is not compiled in, do nothing
    call check_spglib(spglib_was_compiled)
    if (.not.spglib_was_compiled) return

    call get_symmetry_for_lattice &
         (lattice_vector, n_atoms, species, coords, symprec, spacegroup_found, &
          n_operations_sym, spacegroup_number, international_symbol, schoenflies)

    call fjson_write_name_value(fj_h, lattice_name // "_sym_precision", symprec)

    if (spacegroup_number /= 0) then
       call fjson_write_name_value(fj_h, lattice_name // "_spacegroup_found", .true.)
       call fjson_write_name_value(fj_h, lattice_name // "_n_operations_sym", n_operations_sym)
       call fjson_write_name_value(fj_h, lattice_name // "_spacegroup_number", spacegroup_number)
       call fjson_write_name_value(fj_h, lattice_name // "_international_symbol", trim(adjustl(international_symbol)))
       call fjson_write_name_value(fj_h, lattice_name // "_schoenflies", trim(adjustl(schoenflies)))
    else
       call fjson_write_name_value(fj_h, lattice_name // "_spacegroup_found", .false.)
    end if

  end subroutine write_symmetry_to_json_handle
  !******

  !****f* json_output/write_lattice_to_json_handle
  !*  NAME
  !*    write_lattice_to_json_handle
  !*  SYNOPSIS
  subroutine write_lattice_to_json_handle(fj_h, lattice_name, lat_vec, unit_conv, output_symm)
  !*  PURPOSE
  !*    Writes out the lattice and associated information as a sequence of
  !*    name/value pairs.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use constants, only: pi
    use dimensions, only: n_atoms
    use FortJSON, only: fjson_write_value, fjson_start_name_array, &
        fjson_finish_array, fjson_write_name_value
    use geometry, only: species, coords
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
    character(len=*), intent(in) :: lattice_name
    real*8, intent(in) :: lat_vec(3,3)
    real*8, intent(in) :: unit_conv
    logical, intent(in) :: output_symm
  !*  INPUTS
  !*    o fj_h  - FortJSON handle to write to
  !*    o lattice_name - the name of the lattice, to be output as a prefix
  !*    o lat_vec - lattice vectors to extract information from
  !*    o unit_conv - Factor for unit conversion to apply to lattice
  !*    o output_symm - whether symmetry should be output
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    real*8 :: a, b, c
    real*8 :: alpha, beta, gamma_ang ! gamma is a built-in function it turns out
    real*8 :: volume
    real*8 :: lat_vec_conv(3,3) ! Lattice vectors with units converted
    integer :: i

    character(*), parameter :: func = 'write_lattice_to_json_handle'

    lat_vec_conv = lat_vec * unit_conv

    a = sqrt(dot_product(lat_vec_conv(:,1), lat_vec_conv(:,1)))
    b = sqrt(dot_product(lat_vec_conv(:,2), lat_vec_conv(:,2)))
    c = sqrt(dot_product(lat_vec_conv(:,3), lat_vec_conv(:,3)))

    if (b .gt. 0 .and. c .gt. 0) then
      alpha = acos(dot_product(lat_vec_conv(:,2), lat_vec_conv(:,3)) / b / c) * 180. / pi
    else
      alpha = 0.
    end if
    if (a .gt. 0 .and. c .gt. 0) then
      beta  = acos(dot_product(lat_vec_conv(:,1), lat_vec_conv(:,3)) / a / c) * 180. / pi
    else
      beta  = 0.
    end if
    if (a .gt. 0 .and. b .gt. 0) then
      gamma_ang = acos(dot_product(lat_vec_conv(:,1), lat_vec_conv(:,2)) / a / b) * 180. / pi
    else
      gamma_ang = 0.
    end if

    ! Inlined from get_cell_volume() in the bravais module (didn't want to
    ! trigger aims_stop)
    volume = &
       abs(   lat_vec_conv(1,1) * lat_vec_conv(2,2) * lat_vec_conv(3,3)  &
            + lat_vec_conv(1,2) * lat_vec_conv(2,3) * lat_vec_conv(3,1)  &
            + lat_vec_conv(1,3) * lat_vec_conv(2,1) * lat_vec_conv(3,2)  &
            - lat_vec_conv(1,1) * lat_vec_conv(2,3) * lat_vec_conv(3,2)  &
            - lat_vec_conv(1,2) * lat_vec_conv(2,1) * lat_vec_conv(3,3)  &
            - lat_vec_conv(1,3) * lat_vec_conv(2,2) * lat_vec_conv(3,1)  &
          )

    call fjson_start_name_array(fj_h, lattice_name // "_vectors")
    do i = 1, 3
      call fjson_write_value(fj_h, lat_vec_conv(i,:), 3)
    end do
    call fjson_finish_array(fj_h)

    call fjson_write_name_value(fj_h, lattice_name // "_a", a)
    call fjson_write_name_value(fj_h, lattice_name // "_b", b)
    call fjson_write_name_value(fj_h, lattice_name // "_c", c)
    call fjson_write_name_value(fj_h, lattice_name // "_alpha", alpha)
    call fjson_write_name_value(fj_h, lattice_name // "_beta", beta)
    call fjson_write_name_value(fj_h, lattice_name // "_gamma", gamma_ang)
    call fjson_write_name_value(fj_h, lattice_name // "_volume", volume)
    if (output_symm) then
      ! Nothing in symmetry output depends on units, so use original quantities
      call write_symmetry_to_json_handle(fj_h, lattice_name, lat_vec, n_atoms, species, coords)
    end if
  end subroutine write_lattice_to_json_handle
  !******

  !****f* json_output/write_lattice_to_json_handle
  !*  NAME
  !*    write_lattice_to_json_handle
  !*  SYNOPSIS
  subroutine write_stats_int_to_json_handle(fj_h, stats_name, stats)
  !*  PURPOSE
  !*    Writes out the the information contained in a stats_int object.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use FortJSON, only: fjson_start_name_object, fjson_finish_object, &
        fjson_write_name_value
    use statistics, only: stats_int
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
    character(len=*), intent(in) :: stats_name
    type(stats_int), intent(in) :: stats
  !*  INPUTS
  !*    o fj_h  - FortJSON handle to write to
  !*    o stats_name - the name of the stats object, to be output as a prefix
  !*    o stats - the stats object to be written out
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    call fjson_write_name_value(fj_h, stats_name // "_is_valid", stats%is_valid)
    if (stats%is_valid) then
      call fjson_write_name_value(fj_h, stats_name // "_count", stats%num_values)
      call fjson_write_name_value(fj_h, stats_name // "_mean", stats%mean_value)
      call fjson_write_name_value(fj_h, stats_name // "_std", stats%std_dev)
      call fjson_write_name_value(fj_h, stats_name // "_min", stats%min_value)
      call fjson_write_name_value(fj_h, stats_name // "_max", stats%max_value)
      call fjson_write_name_value(fj_h, stats_name // "_median", stats%median_value)
      call fjson_write_name_value(fj_h, stats_name // "_mode_value", stats%mode_value)
      call fjson_write_name_value(fj_h, stats_name // "_mode_count", stats%mode_count)
      call fjson_write_name_value(fj_h, stats_name // "_sum", stats%sum_values)

      ! Fractiles information
      call fjson_write_name_value(fj_h, stats_name //"_fractiles_count", stats%n_fractiles)
      if (stats%n_fractiles .gt. 0) then
        call fjson_write_name_value(fj_h, stats_name // "_fractiles", stats%fractiles, stats%n_fractiles+1)
      end if

      ! Histogram information
      call fjson_write_name_value(fj_h, stats_name // "_hist_bins", stats%n_bins)
      if (stats%n_bins .gt. 0) then
        call fjson_write_name_value(fj_h, stats_name // "_hist_bounds", stats%histogram_bounds, stats%n_bins+1)
        call fjson_write_name_value(fj_h, stats_name // "_hist", stats%histogram, stats%n_bins)
      end if
    end if
  end subroutine write_stats_int_to_json_handle
  !******

  !****f* json_output/write_lattice_to_json_handle
  !*  NAME
  !*    write_lattice_to_json_handle
  !*  SYNOPSIS
  subroutine write_atoms_to_json_handle(fj_h, unit_conv)
  !*  PURPOSE
  !*    Writes out the position of atoms and associated information as a
  !*    sequence of name/value pairs.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use dimensions, only: n_atoms, n_periodic
    use FortJSON, only: fjson_write_name_value, fjson_start_name_array, &
        fjson_finish_array, fjson_write_value
    use geometry, only: species, coords, frac_coords
    use species_data, only: species_name, species_m, species_z
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
    real*8, intent(in) :: unit_conv
  !*  INPUTS
  !*    o fj_h  - FortJSON handle to write to
  !*    o unit_conv - Factor for unit conversion to apply to coordinates
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer :: i_atom, i_species
    real*8 :: coords_conv(3) ! Cartesian coordinates with units converted

    character(*), parameter :: func = 'write_atoms_to_json_handle'

    ! Thanksfully, species is an array with n_atoms elements, so we don't need
    ! to process it
    call fjson_write_name_value(fj_h, "atoms_species", species, n_atoms)

    ! Most other atom-specific quantities are species-dependent, so we need to
    ! index them using the species array
    call fjson_start_name_array(fj_h, "atoms_name")
    do i_atom = 1, n_atoms
      i_species = species(i_atom)
      call fjson_write_value(fj_h, trim(adjustl(species_name(i_species))))
    end do
    call fjson_finish_array(fj_h)

    call fjson_start_name_array(fj_h, "atoms_atomic_number")
    do i_atom = 1, n_atoms
      i_species = species(i_atom)
      call fjson_write_value(fj_h, species_z(i_species))
    end do
    call fjson_finish_array(fj_h)

    call fjson_start_name_array(fj_h, "atoms_mass")
    do i_atom = 1, n_atoms
      i_species = species(i_atom)
      call fjson_write_value(fj_h, species_m(i_species))
    end do
    call fjson_finish_array(fj_h)

    ! Coordinates
    call fjson_start_name_array(fj_h, "atoms_cart_coords")
    do i_atom = 1, n_atoms
      coords_conv(:) = coords(:, i_atom) * unit_conv
      call fjson_write_value(fj_h, coords_conv, 3)
    end do
    call fjson_finish_array(fj_h)

    if (n_periodic .gt. 0) then
      call fjson_start_name_array(fj_h, "atoms_frac_coords")
      do i_atom = 1, n_atoms
        call fjson_write_value(fj_h, frac_coords(:, i_atom), 3)
      end do
      call fjson_finish_array(fj_h)
    end if

  end subroutine write_atoms_to_json_handle
  !******

  !****f* json_output/write_elec_struct_to_json_handle
  !*  NAME
  !*    write_elec_struct_to_json_handle
  !*  SYNOPSIS
  subroutine write_elec_struct_to_json_handle(fj_h, KS_eigenvalue, occ_numbers)
  !*  PURPOSE
  !*    Writes out details about the electronic structure as a sequence of
  !*    name/value pairs.
  !*
  !*    This subroutine does *not* return a well-formed JSON text and should not
  !*    be called standalone; it is intended for use when writing out a record
  !*    object.
  !*  USES
    use constants, only: hartree
    use dimensions, only: n_states, n_spin, n_k_points, n_periodic
    use FortJSON, only: fjson_write_name_value
    use pbc_lists, only: k_point_list
    implicit none
  !*  ARGUMENTS
    type(fjson_handle), intent(inout) :: fj_h
    real*8, intent(in) :: KS_eigenvalue(n_states, n_spin, n_k_points)
    real*8, intent(in) :: occ_numbers(n_states, n_spin, n_k_points)
  !*  INPUTS
  !*    o fj_h  - FortJSON handle to write to
  !*  OUTPUTS
  !*    o fj_h - FortJSON handle to write to
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer :: i_spin, i_state, i_k_point
    real*8, dimension(n_spin) :: n_elec_per_spin_channel
    real*8 :: overall_spin_moment

    real*8  :: homo_level
    real*8  :: lumo_level
    real*8  :: homo_occ
    real*8  :: lumo_occ
    integer :: i_kpt_homo
    integer :: i_kpt_lumo
    integer :: i_spin_homo
    integer :: i_spin_lumo
    logical :: found_min_direct_gap
    real*8  :: min_direct_gap
    integer :: i_kpt_min_direct_gap
    integer :: i_spin_min_direct_homo
    integer :: i_spin_min_direct_lumo

    character(*), parameter :: func = 'write_elec_struct_to_json_handle'

    ! The following is essentially a rewrite of output_real_eigenfunctions
    ! to use JSON instead of stdio

    n_elec_per_spin_channel = 0.d0

    ! Calculate the # electrons per spin channel
    if (n_spin .eq. 2) then
      n_elec_per_spin_channel = 0.d0
      do i_k_point = 1, n_k_points
        do i_spin = 1, n_spin, 1
          do i_state = 1, n_states, 1
            n_elec_per_spin_channel(i_spin) = n_elec_per_spin_channel(i_spin) + &
                 occ_numbers(i_state, i_spin, i_k_point)
          end do
        end do
      end do
      if (n_k_points .gt. 1) then
        n_elec_per_spin_channel(:) = n_elec_per_spin_channel(:) / (n_k_points * 1.0d0)
      end if
      overall_spin_moment = 0.5d0 * abs(n_elec_per_spin_channel(2) - n_elec_per_spin_channel(1))
    end if

    call find_homo_lumo_gap &
      ( KS_eigenvalue, occ_numbers, homo_level, lumo_level, homo_occ, lumo_occ,&
        i_kpt_homo, i_kpt_lumo, i_spin_homo, i_spin_lumo, found_min_direct_gap,&
        min_direct_gap, i_kpt_min_direct_gap, i_spin_min_direct_homo, &
        i_spin_min_direct_lumo )

    ! Overall spin polarization information
    if (n_spin .eq. 2) then
      call fjson_write_name_value(fj_h, "n_elec_per_spin_channel", n_elec_per_spin_channel, n_spin)
      call fjson_write_name_value(fj_h, "overall_spin_moment", overall_spin_moment)
    end if

    ! HOMO/LUMO information
    call fjson_write_name_value(fj_h, "homo_level", homo_level*hartree)
    call fjson_write_name_value(fj_h, "lumo_level", lumo_level*hartree)
    call fjson_write_name_value(fj_h, "homo_occ", homo_occ)
    call fjson_write_name_value(fj_h, "lumo_occ", lumo_occ)
    if (n_periodic.gt.0) then
      call fjson_write_name_value(fj_h, "i_kpt_homo", i_kpt_homo)
      call fjson_write_name_value(fj_h, "i_kpt_lumo", i_kpt_lumo)
      call fjson_write_name_value(fj_h, "kpt_homo", k_point_list(i_kpt_homo,:), 3)
      call fjson_write_name_value(fj_h, "kpt_lumo", k_point_list(i_kpt_lumo,:), 3)
    end if
    if (n_spin.gt.1) then
      call fjson_write_name_value(fj_h, "i_spin_homo", i_spin_homo)
      call fjson_write_name_value(fj_h, "i_spin_lumo", i_spin_lumo)
    end if

    ! Gap information
    if ((i_kpt_homo.gt.0) .and. (i_kpt_lumo.gt.0)) then
      call fjson_write_name_value(fj_h, "homo_lumo_gap", (lumo_level-homo_level)*hartree)
      if (i_kpt_homo .ne. i_kpt_lumo ) then
        call fjson_write_name_value(fj_h, "gap_type", "indirect")
        call fjson_write_name_value(fj_h, "found_min_direct_gap", found_min_direct_gap)
        if (found_min_direct_gap) then
          call fjson_write_name_value(fj_h, "min_direct_gap", min_direct_gap*hartree)
          call fjson_write_name_value(fj_h, "i_kpt_min_direct_gap", i_kpt_min_direct_gap)
          call fjson_write_name_value(fj_h, "kpt_min_direct_gap", k_point_list(i_kpt_min_direct_gap,:), 3)
          if (i_spin.gt.1) then
            call fjson_write_name_value(fj_h, "i_spin_min_direct_homo", i_spin_min_direct_homo)
            call fjson_write_name_value(fj_h, "i_spin_min_direct_lumo", i_spin_min_direct_lumo)
          end if
        end if
      else
        call fjson_write_name_value(fj_h, "gap_type", "direct")
      end if
    end if
  end subroutine write_elec_struct_to_json_handle
  !******

  !****f* json_output/write_geometry_to_json_log
  !*  NAME
  !*    write_geometry_to_json_log
  !*  SYNOPSIS
  subroutine write_geometry_to_json_log(geometry_type)
  !*  PURPOSE
  !*    Writes out information about the current geometry as a record object in
  !*    the JSON log.
  !*  USES
    use constants, only: bohr
    use dimensions, only: n_periodic
    use FortJSON, only: fjson_start_object, fjson_write_name_value, &
        fjson_start_name_object, fjson_finish_object
    use geometry, only: lattice_vector, recip_lattice_vector
    use mpi_tasks, only: myid
    use runtime_choices, only: out_aims_json_log
    implicit none
  !*  ARGUMENTS
    character(len=*), intent(in) :: geometry_type
  !*  INPUTS
  !*    o geometry_type - description of what kind of geometry this is ("initial", "final", etc.)
  !*  OUTPUTS
  !*    o None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    character(*), parameter :: func = 'write_geometry_to_json_log'

    if (myid.ne.0 .or. .not.out_aims_json_log) return

    call fjson_start_object(json_log_handle)

    call fjson_write_name_value(json_log_handle, "record_type", "geometry")
    call fjson_write_name_value(json_log_handle, "geometry_type", geometry_type)

    call write_versioning_to_json_handle(json_log_handle)

    call write_atoms_to_json_handle(json_log_handle, bohr)
    if (n_periodic .gt. 0) then
      call write_lattice_to_json_handle(json_log_handle, "lattice", lattice_vector, bohr, .true.)
      call write_lattice_to_json_handle(json_log_handle, "recip_lattice", recip_lattice_vector, 1.0d0/bohr, .false.)
    end if

    call fjson_finish_object(json_log_handle)

  end subroutine write_geometry_to_json_log
  !******

  !****f* json_output/write_scf_iteration_to_json_log
  !*  NAME
  !*    write_scf_iteration_to_json_log
  !*  SYNOPSIS
  subroutine write_scf_iteration_to_json_log(AS_counter)
  !*  PURPOSE
  !*    Writes out information about the SCF iteration just completed as a
  !*    record object in the JSON log.
  !*  USES
    use batch_statistics, only: batch_sizes_stats, n_max_compute_ham_stats, &
        ld_ham_integ_stats, ld_ham_hpot_stats, ld_ham_density_stats, &
        batch_n_compute_stats
    use dimensions
    use FortJSON, only: fjson_write_name_value, fjson_start_object, &
        fjson_finish_object
    use mpi_tasks, only: aims_stop, myid
    use physics
    use precondition, only: kerker_preconditioner_on
    use runtime_choices
    use timing
    implicit none
  !*  ARGUMENTS
    integer, intent(in) :: AS_counter
  !*  INPUTS
  !*    o AS_counter - Current step in analytical stress tensor calculation
  !*  OUTPUTS
  !*    o None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer :: current_total_number_loops
    integer :: i
    character(LEN=200) :: info_str

    character(*), parameter :: func = 'write_scf_iteration_to_json_log'

    if (myid.ne.0 .or. .not.out_aims_json_log) return

    ! Generally I use the variable name in aims, but in cases where it's unclear
    ! or deceptive, I've renamed the output

    ! total_number_of_loops isn't updated until the end of an SCF cycle
    current_total_number_loops = total_number_of_loops + number_of_loops

    call fjson_start_object(json_log_handle)

    call fjson_write_name_value(json_log_handle, "record_type", "scf_iteration")

    ! Iterations
    call fjson_write_name_value(json_log_handle, "iteration", current_total_number_loops)
    call fjson_write_name_value(json_log_handle, "scf_cycle_number", number_of_scf_inits)
    call fjson_write_name_value(json_log_handle, "iteration_this_cycle", number_of_loops)

    ! Write out standard blocks of informations
    call write_versioning_to_json_handle(json_log_handle)
    call write_core_aims_vars_to_json_handle(json_log_handle)
    call write_gpu_info_to_json_handle(json_log_handle)
    call write_elec_struct_to_json_handle(json_log_handle, KS_eigenvalue, occ_numbers)

    ! Runtime Settings
    call fjson_write_name_value(json_log_handle, "forces_on", forces_on)
    call fjson_write_name_value(json_log_handle, "AS_stress_on", AS_stress_on)
    if (AS_stress_on) then
      call fjson_write_name_value(json_log_handle, "AS_counter", AS_counter)
    end if

    ! Statistics related to batch integration
    call fjson_write_name_value(json_log_handle, "n_hamiltonian_matrix_size", n_hamiltonian_matrix_size)
    call write_stats_to_json_handle(json_log_handle, "batch_sizes", batch_sizes_stats)
    call write_stats_to_json_handle(json_log_handle, "batch_n_compute", batch_n_compute_stats)
    call write_stats_to_json_handle(json_log_handle, "n_max_compute_ham", n_max_compute_ham_stats)
    call write_stats_to_json_handle(json_log_handle, "ld_ham_integ", ld_ham_integ_stats)
    call write_stats_to_json_handle(json_log_handle, "ld_ham_hpot", ld_ham_hpot_stats)
    call write_stats_to_json_handle(json_log_handle, "ld_ham_density", ld_ham_density_stats)

    ! Technical settings (mostly preconditioning and mixing)
    call fjson_write_name_value(json_log_handle, "ewald_radius", Ewald_radius)
    call fjson_write_name_value(json_log_handle, "basis_threshold", basis_threshold)
    if (use_kerker_preconditioner.and.kerker_preconditioner_on) then
      call fjson_write_name_value(json_log_handle, "kerker_preconditioner_used", .true.)
    else
      call fjson_write_name_value(json_log_handle, "kerker_preconditioner_used", .false.)
    end if
    if (use_density_matrix) then
      if (split_updates) then
        call fjson_write_name_value(json_log_handle, "density_update_method", "split_update_methods")
      else
        call fjson_write_name_value(json_log_handle, "density_update_method", "density_matrix")
      end if
    else
      call fjson_write_name_value(json_log_handle, "density_update_method", "orbital")
    end if
    ! This variable is defined by hard-coded numbers throughout aims.
    select case (occupation_type)
    case (0)
      call fjson_write_name_value(json_log_handle, "occupation_type", "gaussian")
    case (1)
      call fjson_write_name_value(json_log_handle, "occupation_type", "fermi")
    case (2)
      call fjson_write_name_value(json_log_handle, "occupation_type", "methfessel-paxton")
    case (3)
      call fjson_write_name_value(json_log_handle, "occupation_type", "integer")
    case (4)
      call fjson_write_name_value(json_log_handle, "occupation_type", "cubic")
    case default
      call aims_stop("Unsupported occupation type", func)
    end select
    call fjson_write_name_value(json_log_handle, "occupation_width", occupation_width*hartree)
    select case (mixer)
    case (MIX_LINEAR)
      call fjson_write_name_value(json_log_handle, "mixer", "MIX_LINEAR")
      call fjson_write_name_value(json_log_handle, "linear_mix_param", linear_mix_param(0))
      if (n_spin.eq.2) then
        call fjson_write_name_value(json_log_handle, "linear_mix_param_spin", linear_mix_param(1))
      end if
    case (MIX_PULAY)
      call fjson_write_name_value(json_log_handle, "mixer", "MIX_PULAY")
      ! WPH: From reading the code, the Pulay mixer uses charge_mix_param during
      ! its linear phase as well as its Pulay phase
      call fjson_write_name_value(json_log_handle, "ini_linear_mixing", ini_linear_mixing)
      call fjson_write_name_value(json_log_handle, "charge_mix_param", charge_mix_param(1))
      if (n_spin.eq.2) then
        call fjson_write_name_value(json_log_handle, "change_mix_param_spin", charge_mix_param(2))
      end if
    case (MIX_BROYDEN)
      call fjson_write_name_value(json_log_handle, "mixer", "MIX_BROYDEN")
      call fjson_write_name_value(json_log_handle, "charge_mix_param", charge_mix_param(1))
      if (n_spin.eq.2) then
        call fjson_write_name_value(json_log_handle, "change_mix_param_spin", charge_mix_param(2))
      end if
    case default
      call aims_stop("Unsupported mixer type", func)
    end select

    ! Energy and SCF convergence parameters
    call fjson_write_name_value(json_log_handle, "total_energy", total_energy*hartree)
    call fjson_write_name_value(json_log_handle, "sc_accuracy_rho", sc_accuracy_rho)
    call fjson_write_name_value(json_log_handle, "sc_accuracy_eev", sc_accuracy_eev)
    call fjson_write_name_value(json_log_handle, "sc_accuracy_etot", sc_accuracy_etot)
    if (forces_on) then
      call fjson_write_name_value(json_log_handle, "sc_accuracy_forces", sc_accuracy_forces)
    end if
    if ( AS_stress_on.and. AS_counter.ge.3) then
      call fjson_write_name_value(json_log_handle, "sc_accuracy_stress", sc_accuracy_stress)
    end if
    call fjson_write_name_value(json_log_handle, "change_charge_density", rho_change(1))
    if (n_spin.eq. 2) then
      call fjson_write_name_value(json_log_handle, "change_spin_density", rho_change(2))
    end if
    call fjson_write_name_value(json_log_handle, "change_sum_eigenvalues", (ev_sum - previous_ev_sum)*hartree)
    call fjson_write_name_value(json_log_handle, "change_total_energy", (total_energy - previous_total_energy)*hartree)
    if (forces_on) then
      call fjson_write_name_value(json_log_handle, "change_forces", diff_forces)
    end if
    if ( AS_stress_on.and. AS_counter.ge.3) then
      call fjson_write_name_value(json_log_handle, "change_anal_stress", diff_stress)
    end if

    ! Time for this iteration
    call fjson_write_name_value(json_log_handle, "time_sc_loop", time_sc_loop)
    call fjson_write_name_value(json_log_handle, "clock_time_sc_loop", clock_time_sc_loop)

    ! Charge density & force component update
    ! TODO: Change name based on forces/stress
    call fjson_write_name_value(json_log_handle, "time_density", time_density)
    call fjson_write_name_value(json_log_handle, "clock_time_density", clock_time_density)

    ! Density mixing & preconditioning
    if (use_kerker_preconditioner.and.kerker_preconditioner_on) then
      call fjson_write_name_value(json_log_handle, "time_mixing_preconditioner", time_mixing)
      call fjson_write_name_value(json_log_handle, "clock_time_mixing_perconditioner", clock_time_mixing)
    else
      call fjson_write_name_value(json_log_handle, "time_mixing", time_mixing)
      call fjson_write_name_value(json_log_handle, "clock_time_mixing", clock_time_mixing)
    endif

    ! Hartree multipole update
    call fjson_write_name_value(json_log_handle, "time_hartree_multi", time_hartree_multi)
    call fjson_write_name_value(json_log_handle, "clock_time_hartree_multi", clock_time_hartree_multi)

    ! Hartree multipole summation
    call fjson_write_name_value(json_log_handle, "time_hartree_sum", time_hartree_sum)
    call fjson_write_name_value(json_log_handle, "clock_time_hartree_sum", clock_time_hartree_sum)

    ! Integration
    call fjson_write_name_value(json_log_handle, "time_intg", time_intg)
    call fjson_write_name_value(json_log_handle, "clock_time_intg", clock_time_intg)

    ! Fock matrix evaluation
    if (use_hartree_fock) then
      call fjson_write_name_value(json_log_handle, "time_fock", time_fock)
      call fjson_write_name_value(json_log_handle, "clock_time_fock", clock_time_fock)
    end if

    ! Solution of K.-S. eqns.
    call fjson_write_name_value(json_log_handle, "time_diag", time_diag)
    call fjson_write_name_value(json_log_handle, "clock_time_diag", clock_time_diag)

    ! Total energy evaluation
    call fjson_write_name_value(json_log_handle, "time_etot", time_etot)
    call fjson_write_name_value(json_log_handle, "clock_time_etot", clock_time_etot)

    call fjson_finish_object(json_log_handle)
  end subroutine write_scf_iteration_to_json_log
  !******

  !****f* json_output/write_mulliken_to_json_log
  !*  NAME
  !*    write_mulliken_to_json_log
  !*  SYNOPSIS
  subroutine write_mulliken_to_json_log &
       (n_spin_proj, at_projected_charge, l_projected_charge, spin_per_atom, &
        total_charge, charge_difference)
  !*  PURPOSE
  !*    Writes out (per-atom) Mulliken information as a record object in the
  !*    JSON log.  This is essentially an extended version of
  !*    write_atoms_to_json_handle.
  !*  USES
    use constants, only: bohr
    use dimensions, only: n_atoms, n_periodic, l_wave_max
    use FortJSON, only: fjson_start_object, fjson_write_name_value, &
        fjson_start_name_array, fjson_finish_array, fjson_write_value, &
        fjson_start_object, fjson_finish_object, fjson_start_array
    use geometry, only: lattice_vector, recip_lattice_vector, species, coords, &
        frac_coords
    use mpi_tasks, only: myid
    use runtime_choices, only: out_aims_json_log, spin_treatment
    use species_data, only: l_shell_max, species_name, species_m, species_z, &
        species_pseudoized
    implicit none
  !*  ARGUMENTS
    integer, intent(in) :: n_spin_proj
    real*8,  intent(in) :: at_projected_charge(n_atoms, n_spin_proj)
    real*8 , intent(in) :: l_projected_charge(0:l_wave_max, n_atoms, n_spin_proj)
    real*8,  intent(in) :: spin_per_atom(n_atoms)
    real*8,  intent(in) :: total_charge(n_spin_proj)
    real*8,  intent(in) :: charge_difference
  !*  INPUTS
  !*    o n_spin_proj - Number of (projected) spin channels
  !*    o at_project_charge - Charge projected on each atom
  !*    o l_projected_charge - Charge projected on each atom and l channel
  !*    o spin_per_atom- Spin projected on each atom
  !*    o at_project_charge - Total charge (per spin channel)
  !*    o charge_deviation - Deviation from ideality
  !*  OUTPUTS
  !*    o None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer :: i_atom, i_species, i_l
    real*8 :: electrons
    real*8 :: coords_conv(3) ! Coordinates converted from bohr to Angstrom

    character(*), parameter :: func = 'write_mulliken_to_json_log'

    if (myid.ne.0 .or. .not.out_aims_json_log) return

    call fjson_start_object(json_log_handle)

    call fjson_write_name_value(json_log_handle, "record_type", "mulliken")

    call write_versioning_to_json_handle(json_log_handle)
    call write_core_aims_vars_to_json_handle(json_log_handle)
    call write_atoms_to_json_handle(json_log_handle, bohr)

    ! Projected # of electrons
    call fjson_start_name_array(json_log_handle, "atoms_proj_electrons")
    do i_atom = 1, n_atoms
      i_species = species(i_atom)
      if(species_pseudoized(i_species)) cycle

      electrons = sum(at_projected_charge(i_atom, 1:n_spin_proj))
      call fjson_write_value(json_log_handle, electrons)
    end do
    call fjson_finish_array(json_log_handle)

    ! Projected charge
    call fjson_start_name_array(json_log_handle, "atoms_proj_charge")
    do i_atom = 1, n_atoms
      i_species = species(i_atom)
      if(species_pseudoized(i_species)) cycle

      electrons = sum(at_projected_charge(i_atom, 1:n_spin_proj))
      call fjson_write_value(json_log_handle, electrons - species_z(i_species))
    end do
    call fjson_finish_array(json_log_handle)

    ! Projected electrons per l-channel
    call fjson_start_name_array(json_log_handle, "atoms_proj_electrons_l")
    do i_atom = 1, n_atoms
      i_species = species(i_atom)
      if(species_pseudoized(i_species)) cycle

      call fjson_start_array(json_log_handle)
      do i_l = 0, l_shell_max(species(i_atom))
        call fjson_write_value(json_log_handle, sum(l_projected_charge(i_l,i_atom,1:n_spin_proj)))
      end do
      call fjson_finish_array(json_log_handle)
    end do
    call fjson_finish_array(json_log_handle)

    if (spin_treatment .eq. 1) then
      ! Projected overall spin
      call fjson_write_name_value(json_log_handle, "atoms_proj_spin", spin_per_atom, n_atoms)

      ! Projected spin per l-channel
      call fjson_start_name_array(json_log_handle, "spin_proj_projected_l")
      do i_atom = 1, n_atoms
        call fjson_start_array(json_log_handle)
        do i_l = 0, l_shell_max(species(i_atom))
          call fjson_write_value(json_log_handle, l_projected_charge(i_l,i_atom,1) - l_projected_charge(i_l,i_atom,2))
        end do
        call fjson_finish_array(json_log_handle)
      end do
      call fjson_finish_array(json_log_handle)
    end if

    call fjson_write_name_value(json_log_handle, "total_electrons", sum(total_charge(1:n_spin_proj)))
    call fjson_write_name_value(json_log_handle, "charge_difference", -charge_difference)
    call fjson_write_name_value(json_log_handle, "total_spin", sum(spin_per_atom(1:n_atoms)))

    call fjson_finish_object(json_log_handle)
  end subroutine write_mulliken_to_json_log
  !******

  !****f* json_output/write_final_output_to_json_log
  !*  NAME
  !*    write_final_output_to_json_log
  !*  SYNOPSIS
  subroutine write_final_output_to_json_log(is_a_nice_day)
  !*  PURPOSE
  !*    Writes out information about the final energies and final timings
  !*    together as an record object in the JSON log.
  !*  USES
    use dimensions
    use FortJSON, only: fjson_write_name_value, fjson_start_object, &
        fjson_finish_object
    use mpi_tasks, only: aims_stop, myid
    use physics
    use runtime_choices
    use timing
    implicit none
  !*  ARGUMENTS
    logical, intent(in) :: is_a_nice_day
  !*  INPUTS
  !*    o is_a_nice_day - Did the SCF cycle converge?
  !*  OUTPUTS
  !*    o None
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    integer :: current_total_number_loops
    integer :: i
    real*8  :: time_total_json, clock_time_total_json
    character(LEN=200) :: info_str

    character(*), parameter :: func = 'write_final_output_to_json_log'

    ! WPH: We do not output time_total and clock_time_total directly for
    !      architectural reasons.  The final_output record must be output before
    !      final_deallocations() occurs, as there are certain arrays that will
    !      be undefined otherwise.  However, final_deallocations() must occur
    !      before final_timings_and_memory(), because otherwise the tracked
    !      residual memory will not be zero, which is needed for the aims memory
    !      tracking to give useful information.
    !      Therefore, I have opted to output the total time that has elapsed up
    !      to this record being written out.  This will cause a divergence
    !      between the stdio and JSON output, but this divergence should be
    !      small unless there's serious IO issues.
    call get_total_time_elapsed(time_total_json, clock_time_total_json)

    if (myid.ne.0 .or. .not.out_aims_json_log) return

    ! Start the record
    call fjson_start_object(json_log_handle)

    call fjson_write_name_value(json_log_handle, "record_type", "final_output")

    ! Write out standard blocks of informations
    call write_versioning_to_json_handle(json_log_handle)
    call write_core_aims_vars_to_json_handle(json_log_handle)
    call write_gpu_info_to_json_handle(json_log_handle)

    ! Final energy output
    if (.not.calculate_atom_bsse) then
      call fjson_write_name_value(json_log_handle, "total_energy", total_energy*hartree)
      if (en_vdw.ne.0.d0) then
        call fjson_write_name_value(json_log_handle, "en_vdw", en_vdw*hartree)
      end if
      call fjson_write_name_value(json_log_handle, "final_zero_broad_corr_energy", &
           (total_energy + 2 * (dble(n_methfessel_paxton + 1) / &
            dble(n_methfessel_paxton + 2)) * entropy_correction) &
            * hartree)
      if (post_scf_total_energy.ne.0.d0) then
        call fjson_write_name_value(json_log_handle, "post_scf_total_energy", post_scf_total_energy*hartree)
      end if
      if (SOC_non_sc_total_energy.ne.0.0d0) then
        call fjson_write_name_value(json_log_handle, "SOC_non_sc_total_energy_DO_NOT_USE", SOC_non_sc_total_energy*hartree)
      end if
      call fjson_write_name_value(json_log_handle, "hartree", hartree)
    end if

    ! Iterations
    call fjson_write_name_value(json_log_handle, "total_number_of_loops", total_number_of_loops)
    call fjson_write_name_value(json_log_handle, "number_of_scf_inits", number_of_scf_inits)
    if (use_geo_relaxation) then
      call fjson_write_name_value(json_log_handle, "relaxation_step_number", relaxation_step_number)
    end if
    if (calculate_atom_bsse) then
      call fjson_write_name_value(json_log_handle, "current_atom_num_bsse", current_atom_num_bsse)
    end if
    if (use_molecular_dynamics) then
      call fjson_write_name_value(json_log_handle, "MD_stepcount", MD_stepcount)
      call fjson_write_name_value(json_log_handle, "MD_force_evaluations", MD_force_evaluations)
    end if

    ! Timings
    call fjson_write_name_value(json_log_handle, "time_total", time_total_json)
    call fjson_write_name_value(json_log_handle, "clock_time_total", clock_time_total_json)
    if (use_gpu) then
      call fjson_write_name_value(json_log_handle, "tot_time_gpu_init", tot_time_gpu_init)
      call fjson_write_name_value(json_log_handle, "tot_time_clock_gpu_init", tot_clock_time_gpu_init)
    end if
    call fjson_write_name_value(json_log_handle, "time_prep", time_prep)
    call fjson_write_name_value(json_log_handle, "clock_time_prep", clock_time_prep)
    call fjson_write_name_value(json_log_handle, "tot_time_bc_init", tot_time_bc_init)
    call fjson_write_name_value(json_log_handle, "tot_clock_time_bc_init", tot_clock_time_bc_init)
    call fjson_write_name_value(json_log_handle, "tot_time_partition", tot_time_partition)
    call fjson_write_name_value(json_log_handle, "tot_clock_time_partition", tot_clock_time_partition)
    call fjson_write_name_value(json_log_handle, "tot_time_free_offset", tot_time_free_offset)
    call fjson_write_name_value(json_log_handle, "tot_clock_time_free_offset", tot_clock_time_free_offset)
    call fjson_write_name_value(json_log_handle, "tot_time_superpos", tot_time_superpos)
    call fjson_write_name_value(json_log_handle, "tot_clock_time_superpos", tot_clock_time_superpos)
    call fjson_write_name_value(json_log_handle, "tot_time_intg", tot_time_intg)
    call fjson_write_name_value(json_log_handle, "tot_clock_time_intg", tot_clock_time_intg)
    if(use_prodbas) then
       call fjson_write_name_value(json_log_handle, "tot_time_prodbas_total", tot_time_prodbas_total)
       call fjson_write_name_value(json_log_handle, "tot_clock_time__prodbas_total", tot_clock_time_prodbas_total)

       if (.not.use_hf_realspace) then
          if (.not. use_lvl .or. .not. use_logsbt_lvltriples) then
             call fjson_write_name_value(json_log_handle, "tot_time_ovlp3fn", tot_time_ovlp3fn)
             call fjson_write_name_value(json_log_handle, "tot_clock_time_ovlp3fn", tot_clock_time_ovlp3fn)
          end if
          call fjson_write_name_value(json_log_handle, "tot_time_coulomb_matr", tot_time_coulomb_matr)
          call fjson_write_name_value(json_log_handle, "tot_clock_time_coulomb_matr", tot_clock_time_coulomb_matr)
          if (RI_type /= RI_LVL) then
             call fjson_write_name_value(json_log_handle, "tot_time_inv_coulomb_matr", tot_time_inv_coulomb_matr)
             call fjson_write_name_value(json_log_handle, "tot_clock_time_inv_coulomb_matr", tot_clock_time_inv_coulomb_matr)
          end if
          if (use_lvl .and. use_logsbt_lvltriples) then
             call fjson_write_name_value(json_log_handle, "tot_time_ovlp_multi", tot_time_ovlp_multi)
             call fjson_write_name_value(json_log_handle, "tot_clock_time_ovlp_multi", tot_clock_time_ovlp_multi)
          else
             call fjson_write_name_value(json_log_handle, "tot_time_ovlp_multi", tot_time_ovlp_multi)
             call fjson_write_name_value(json_log_handle, "tot_clock_time_ovlp_multi", tot_clock_time_ovlp_multi)
          end if
       end if

       if ((n_periodic .gt. 0) .and. (use_hf_kspace .or. use_rpa_ene .or. use_gw)) then
            call fjson_write_name_value(json_log_handle, "tot_time_lvl_triple_recip", tot_time_lvl_triple_recip)
            call fjson_write_name_value(json_log_handle, "tot_clock_time_lvl_triple_recip", tot_clock_time_lvl_triple_recip)
       endif
       if(use_hartree_fock) then
          call fjson_write_name_value(json_log_handle, "tot_time_fock", tot_time_fock)
          call fjson_write_name_value(json_log_handle, "tot_clock_time_fock", tot_clock_time_fock)
       endif
       if(use_ovlp_swap) then
          call fjson_write_name_value(json_log_handle, "time_ovlp_swap", time_ovlp_swap)
          call fjson_write_name_value(json_log_handle, "clock_time_ovlp_swap", clock_time_ovlp_swap)
       endif
    endif
    if(use_full_spectrum) then
       call fjson_write_name_value(json_log_handle, "time_ovlp3fn_to_ovlp3KS", time_ovlp3fn_to_ovlp3KS)
       call fjson_write_name_value(json_log_handle, "clock_time_ovlp3fn_to_ovlp3KS", clock_time_ovlp3fn_to_ovlp3KS)
    endif
    if(use_gw) then
       call fjson_write_name_value(json_log_handle, "time_total_gw", time_total_gw)
       call fjson_write_name_value(json_log_handle, "clock_time_total_gw", clock_time_total_gw)
       call fjson_write_name_value(json_log_handle, "time_self_energy", time_self_energy)
       call fjson_write_name_value(json_log_handle, "clock_time_self_energy", clock_time_self_energy)
       call fjson_write_name_value(json_log_handle, "time_Wmn_imagfreq", time_Wmn_imagfreq)
       call fjson_write_name_value(json_log_handle, "clock_time_Wmn_imagfreq", clock_time_Wmn_imagfreq)
       call fjson_write_name_value(json_log_handle, "time_polar_imagfreq", time_polar_imagfreq)
       ! The lack of clock time isn't an oversight here, the original timings are similar
       if(out_band) then
         call fjson_write_name_value(json_log_handle, "time_band_self_energy", time_band_self_energy)
         call fjson_write_name_value(json_log_handle, "clock_time_band_self_energy", clock_time_band_self_energy)
       endif
       if(use_gw2ox) then
         call fjson_write_name_value(json_log_handle, "time_2ox_selfenergy", time_2ox_selfenergy)
         call fjson_write_name_value(json_log_handle, "clock_time_2ox_selfenergy", clock_time_2ox_selfenergy)
       elseif(use_gws2ox) then
         call fjson_write_name_value(json_log_handle, "time_s2ox_selfenergy", time_s2ox_selfenergy)
         call fjson_write_name_value(json_log_handle, "clock_time_s2ox_selfenergy", clock_time_s2ox_selfenergy)
       elseif(use_gwsoxw) then
         call fjson_write_name_value(json_log_handle, "time_soxw_selfenergy", time_soxw_selfenergy)
         call fjson_write_name_value(json_log_handle, "clock_time_soxw_selfenergy", clock_time_soxw_selfenergy)
       endif
       if(use_contour_def_gw) then
         call fjson_write_name_value(json_log_handle, "time_qp_energy_cd", time_qp_energy_cd)
         call fjson_write_name_value(json_log_handle, "clock_time_qp_energy_cd", clock_time_qp_energy_cd)
         call fjson_write_name_value(json_log_handle, "time_self_energy_cd", time_self_energy_cd)
         call fjson_write_name_value(json_log_handle, "clock_time_self_energy_cd", clock_time_self_energy_cd)
         call fjson_write_name_value(json_log_handle, "time_WPQ_realfreq_cd", time_WPQ_realfreq_cd)
         call fjson_write_name_value(json_log_handle, "clock_time_WPQ_realfreq_cd", clock_time_WPQ_realfreq_cd)
         if(flag_print_self_energy) then
           call fjson_write_name_value(json_log_handle, "time_print_sigma_cd", time_print_sigma_cd)
           call fjson_write_name_value(json_log_handle, "clock_time_print_sigma_cd", clock_time_print_sigma_cd)
         endif
         if(flag_calc_spectral_func) then
           call fjson_write_name_value(json_log_handle, "time_spec_func_cd", time_spec_func_cd)
           call fjson_write_name_value(json_log_handle, "clock_time_spec_func_cd", clock_time_spec_func_cd)
         endif
       endif
    endif
    if(use_rpa_ene) then
       call fjson_write_name_value(json_log_handle, "time_polar", time_polar)
       call fjson_write_name_value(json_log_handle, "clock_time_polar", clock_time_polar)
       call fjson_write_name_value(json_log_handle, "time_rpa_corr", time_rpa_corr)
       call fjson_write_name_value(json_log_handle, "clock_time_rpa_corr", clock_time_rpa_corr)
       if(n_periodic .gt. 0) then
         call fjson_write_name_value(json_log_handle, "time_rse_corr", time_rse_corr)
         call fjson_write_name_value(json_log_handle, "clock_time_rse_corr", clock_time_rse_corr)
       endif
    endif
    if(use_rpa_plus_2ox) then
       call fjson_write_name_value(json_log_handle, "time_2oex_corr", time_2oex_corr)
       call fjson_write_name_value(json_log_handle, "clock_time_2oex_corr", clock_time_2oex_corr)
    endif
    if(use_rpa_plus_sosex) then
       call fjson_write_name_value(json_log_handle, "time_sosex_corr", time_sosex_corr)
       call fjson_write_name_value(json_log_handle, "clock_time_sosex_corr", clock_time_sosex_corr)
    endif
    if(use_C6_coef) then
       call fjson_write_name_value(json_log_handle, "time_C6_coef", time_C6_coef)
       call fjson_write_name_value(json_log_handle, "clock_time_C6_coef", clock_time_C6_coef)
    endif
    call fjson_write_name_value(json_log_handle, "tot_time_diag", tot_time_diag)
    call fjson_write_name_value(json_log_handle, "tot_clock_time_diag", tot_clock_time_diag)
    if (orthonormalize_evs .and. &
        (use_geo_relaxation.or.use_molecular_dynamics)) then
       call fjson_write_name_value(json_log_handle, "tot_time_on_evs", tot_time_on_evs)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_on_evs", tot_clock_time_on_evs)
    end if
    if (force_potential.eq.0) then
       if (.not.use_forces) then
          call fjson_write_name_value(json_log_handle, "tot_time_density", tot_time_density)
          call fjson_write_name_value(json_log_handle, "tot_clock_time_density", tot_clock_time_density)
       else
          call fjson_write_name_value(json_log_handle, "tot_time_density", tot_time_density)
          call fjson_write_name_value(json_log_handle, "tot_clock_time_density", tot_clock_time_density)
       end if
       if (use_kerker_preconditioner) then
          call fjson_write_name_value(json_log_handle, "tot_time_mixing", tot_time_mixing)
          call fjson_write_name_value(json_log_handle, "tot_clock_time_mixing", tot_clock_time_mixing)
       else
          call fjson_write_name_value(json_log_handle, "tot_time_mixing", tot_time_mixing)
          call fjson_write_name_value(json_log_handle, "tot_clock_time_mixing", tot_clock_time_mixing)
       end if
       call fjson_write_name_value(json_log_handle, "tot_time_hartree_multi", tot_time_hartree_multi)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_hartree_multi", tot_clock_time_hartree_multi)
       call fjson_write_name_value(json_log_handle, "tot_time_hartree_sum", tot_time_hartree_sum)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_hartree_sum", tot_clock_time_hartree_sum)
       call fjson_write_name_value(json_log_handle, "tot_time_etot", tot_time_etot)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_etot", tot_clock_time_etot)
    end if
    if (flag_rel.gt.0) then
       call fjson_write_name_value(json_log_handle, "tot_time_scaled_zora", tot_time_scaled_zora)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_scaled_zora", tot_clock_time_scaled_zora)
    endif
    if (use_vdw_correction_hirshfeld .or. use_mbd_old &
        .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd) then
       call fjson_write_name_value(json_log_handle, "tot_time_vdw", tot_time_vdw)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_vdw", tot_clock_time_vdw)
    endif
    if (use_wf_extrapolation) then
       call fjson_write_name_value(json_log_handle, "tot_time_wf_extra_in", tot_time_wf_extra_in)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_wf_extra_in", tot_clock_time_wf_extra_in)
       call fjson_write_name_value(json_log_handle, "tot_time_wf_extra_out", tot_time_wf_extra_out)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_wf_extra_out", tot_clock_time_wf_extra_out)
    endif

    if (use_embedding_pp.or.use_embedding_potential) then
       call fjson_write_name_value(json_log_handle, "tot_time_embedding", tot_time_embedding)
       call fjson_write_name_value(json_log_handle, "tot_clock_time_embedding", tot_clock_time_embedding)
    endif

    if (magnetic_response) then
       call fjson_write_name_value(json_log_handle, "time_magnetic_response", time_magnetic_response)
       call fjson_write_name_value(json_log_handle, "clock_time_magnetic_response", clock_time_magnetic_response)
    end if

    if(flag_neutral_excitation) then
      call fjson_write_name_value(json_log_handle, "time_total_neutral_excitation", time_total_neutral_excitation)
      call fjson_write_name_value(json_log_handle, "clock_time_total_neutral_excitation", clock_time_total_neutral_excitation)
      call fjson_write_name_value(json_log_handle, "time_neutral_excitation_dipole", time_neutral_excitation_dipole)
      call fjson_write_name_value(json_log_handle, "clock_time_neutral_excitation_dipole", clock_time_neutral_excitation_dipole)
      call fjson_write_name_value(json_log_handle, "time_neutral_excitation_orbitals", time_neutral_excitation_orbitals)
      call fjson_write_name_value(json_log_handle, "clock_time_neutral_excitation_orbitals", clock_time_neutral_excitation_orbitals)
      call fjson_write_name_value(json_log_handle, "time_fxc_integr", time_fxc_integr)
      call fjson_write_name_value(json_log_handle, "clock_time_fxc_integr", clock_time_fxc_integr)
      call fjson_write_name_value(json_log_handle, "time_coeffs_fxc3fn", time_coeffs_fxc3fn)
      call fjson_write_name_value(json_log_handle, "clock_time_coeffs_fxc3fn", clock_time_coeffs_fxc3fn)
      call fjson_write_name_value(json_log_handle, "time_ovlp3fn_to_ovlp3KS", time_ovlp3fn_to_ovlp3KS)
      call fjson_write_name_value(json_log_handle, "clock_time_ovlp3fn_to_ovlp3KS", clock_time_ovlp3fn_to_ovlp3KS)
      call fjson_write_name_value(json_log_handle, "time_matvec_cont", time_matvec_cont)
      call fjson_write_name_value(json_log_handle, "clock_time_matvec_cont", clock_time_matvec_cont)
      call fjson_write_name_value(json_log_handle, "time_neutral_excitation_casida", time_neutral_excitation_casida)
      call fjson_write_name_value(json_log_handle, "clock_time_neutral_excitation_casida", clock_time_neutral_excitation_casida)
      call fjson_write_name_value(json_log_handle, "time_omegabuild", time_omegabuild)
      call fjson_write_name_value(json_log_handle, "clock_time_omegabuild", clock_time_omegabuild)
      call fjson_write_name_value(json_log_handle, "time_3ks_comm", time_3ks_comm)
      call fjson_write_name_value(json_log_handle, "clock_time_3ks_comm", clock_time_3ks_comm)
      call fjson_write_name_value(json_log_handle, "time_casida_matrix_calculation", time_casida_matrix_calculation)
      call fjson_write_name_value(json_log_handle, "clock_time_casida_matrix_calculation", clock_time_casida_matrix_calculation)
      call fjson_write_name_value(json_log_handle, "time_omegaresdist", time_omegaresdist)
      call fjson_write_name_value(json_log_handle, "clock_time_omegaresdist", clock_time_omegaresdist)
      call fjson_write_name_value(json_log_handle, "time_omega_solve", time_omega_solve)
      call fjson_write_name_value(json_log_handle, "clock_time_omega_solve", clock_time_omega_solve)
    endif

    if(calculate_perturbative_soc) then
      call fjson_write_name_value(json_log_handle, "tot_time_soc", tot_time_soc)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_soc", tot_clock_time_soc)
    end if

    if (out_band.or.out_dos) then
      call fjson_write_name_value(json_log_handle, "tot_time_band_dos", tot_time_band_dos)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_band_dos", tot_clock_time_band_dos)
    end if

    if (flag_run_mulliken) then
      call fjson_write_name_value(json_log_handle, "tot_time_mulliken", tot_time_mulliken)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_mulliken", tot_clock_time_mulliken)
    end if

    if (out_cube) then
      call fjson_write_name_value(json_log_handle, "tot_time_cube_output", tot_time_cube_output)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_cube_output", tot_clock_time_cube_output)
    end if

    if(flag_out_dielectric .or. flag_out_absorption) then
      call fjson_write_name_value(json_log_handle, "tot_time_dielectric", tot_time_dielectric)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_dielectric", tot_clock_time_dielectric)
    end if

    if(use_DFPT_reduce_memory) then
      call fjson_write_name_value(json_log_handle, "tot_time_first_order_density", tot_time_first_order_density)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_first_order_density", tot_clock_time_first_order_density)
      call fjson_write_name_value(json_log_handle, "tot_time_first_order_potential", tot_time_first_order_potential)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_first_order_potential", tot_clock_time_first_order_potential)
      call fjson_write_name_value(json_log_handle, "tot_time_first_order_H", tot_time_first_order_H)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_first_order_H", tot_clock_time_first_order_H)
      call fjson_write_name_value(json_log_handle, "tot_time_Sternhemier", tot_time_Sternheimer)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_Sternhemier", tot_clock_time_Sternheimer)
      call fjson_write_name_value(json_log_handle, "tot_time_first_order_S", tot_time_first_order_S)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_first_order_S", tot_clock_time_first_order_S)
      call fjson_write_name_value(json_log_handle, "tot_time_Hessian", tot_time_Hessian)
    endif

    if(fo_finalise) then
      call fjson_write_name_value(json_log_handle, "tot_time_fodft", tot_time_fodft)
      call fjson_write_name_value(json_log_handle, "tot_clock_time_fodft", tot_clock_time_fodft)
    endif

    if (use_gpu) then
      if (use_gpu_density) then
        call fjson_write_name_value(json_log_handle, "gpu_density_always_used", gpu_density_always_used)
      end if
      if (use_gpu_hamiltonian) then
        call fjson_write_name_value(json_log_handle, "gpu_hamiltonian_always_used", gpu_hamiltonian_always_used)
      end if
      if (use_forces.and.use_gpu_forces) then
        call fjson_write_name_value(json_log_handle, "gpu_forces_always_used", gpu_forces_always_used)
      end if
    end if

    ! Did we have a nice day?
    call fjson_write_name_value(json_log_handle, "is_a_nice_day", is_a_nice_day)

    call fjson_finish_object(json_log_handle)
  end subroutine write_final_output_to_json_log
  !******
end module
!******
