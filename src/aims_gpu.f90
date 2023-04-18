!****h* FHI-aims/aims_gpu
!*  NAME
!*    aims_gpu
!*  SYNOPSIS
module aims_gpu
!*  PURPOSE
!*    This module handles the low-level details of the GPU implementation in
!*    FHI-aims from the Fortran side.  The actual GPU acceleration is done via
!*    various subroutines written in C CUDA, which (as of 2018-01-25) may be
!*    found in the cuda subdirectory of the main aims source directory.
!*  USES
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University)
!*  HISTORY
!*    January 2018 - Created.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is
!*    subject to the terms and conditions of the respective license
!*    agreement.
!*  SOURCE

  private

  integer, public :: n_gpus_avail = -1
  integer, public :: gpu_number = -1

  ! Hardware details, only used for output and light error checking
  integer, public :: cuda_major = -1
  integer, public :: cuda_minor = -1
  character(len=:), allocatable, public :: gpu_name
  integer, public :: total_global_gpu_mem = -1
  integer, public :: total_constant_gpu_mem = -1
  integer, public :: shared_gpu_mem_per_block = -1
  integer, public :: gpu_clock_rate = -1
  integer, public :: gpu_multiprocessor_count = -1
  integer, public :: gpu_max_threads_per_multiprocessor = -1

  ! The flags specifying whether GPUs should be used for a calculation are
  ! found in runtime_choices

  ! Were GPUs were actually used for the specified operation?
  logical, public :: gpu_density_used = .false.
  logical, public :: gpu_hamiltonian_used = .false.
  logical, public :: gpu_forces_used = .false.

  ! Flags specifying whether GPUs were *always* used for the specified operation
  ! have been placed in timing, since they are used as part of the final timing
  ! output

  ! Length of name field in cudaDeviceProp struct
  integer, parameter :: GPU_NAME_LEN = 256

  public :: initialize_gpu
  public :: output_mem_array_gpu
  public :: finalize_gpu

contains

!****s* FHI-aims/initialize_gpu
!  NAME
!    initialize_gpu
!  SYNOPSIS
  subroutine initialize_gpu()
!  PURPOSE
!    Distribute available GPUs to CPUs and initialize them.
!  USES
    use c_helper, only: c_char_array_to_fort_str
    use localorb_io, only: localorb_info, localorb_allinfo
    use mpi_tasks, only: myid, aims_stop
    use, intrinsic :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
    use runtime_choices, only: use_gpu
    use timing, only: get_times, get_timestamps, time_gpu_init, &
        clock_time_gpu_init, tot_time_gpu_init, tot_clock_time_gpu_init
    implicit none
!  ARGUMENTS
!    o None
!  INPUTS
!    o None
!  OUTPUT
!    o None
!  AUTHOR
!    Written originally by Bjoern Lange in 20XX, extensively modified by William
!    Huhn.
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    January 2018 - created from cuda_management
!  NOTES
!    Various subroutines called in this subroutine are found in the C CUDA
!    code.
!  SOURCE
    integer :: i_char, i_gpu
    character(kind=c_char, len=1), dimension(GPU_NAME_LEN) :: gpu_name_c
    character(len=GPU_NAME_LEN) :: gpu_name_fort
    character*400 :: info_str

    character(*), parameter :: func = 'initialize_gpu'

    call get_timestamps( time_gpu_init, clock_time_gpu_init )

    write(info_str,'(2X,A)') "Initializing GPUs..."
    call localorb_info( info_str )
    write(info_str,'(2X,A)') ""
    call localorb_info( info_str )

    ! Find out how many GPUs are accessible by the current MPI task
    call get_num_gpus (n_gpus_avail)
    write(info_str,'(2X,A,I9,A,I3)') "Number of GPUs available to MPI task ", &
         myid, " : ", n_gpus_avail
    call localorb_allinfo(info_str)

    ! If there's more than one GPU accessible by the current MPI task, then
    ! decide which one to use via a round-robin-esque distribution
    ! This implicitly assumes that nodes containing GPUs are homogeneous and
    ! that the ranks of MPI tasks are sequentially ordered, i.e. the first X
    ! ranks reside on node 1, the next X ranks reside on node 2, and so on
    if (n_gpus_avail > 0) then
      use_gpu = .true.
      gpu_number = mod(myid,n_gpus_avail)
    else
      use_gpu = .false.
      gpu_number = -1
    end if

    ! Activate CUDA support
    do i_gpu = 0, n_gpus_avail - 1
      if (gpu_number == i_gpu) then
        call initialize_cuda_and_cublas(i_gpu)
      end if
    end do
    if (gpu_number >= 0) call set_gpu(gpu_number)

    ! Obtain relevant versioning information for the GPU attached to this
    ! MPI task
    if (use_gpu) then
      call get_gpu_specs(gpu_number, cuda_major, cuda_minor, &
           total_global_gpu_mem, total_constant_gpu_mem, &
           shared_gpu_mem_per_block, gpu_clock_rate, gpu_multiprocessor_count, &
           gpu_max_threads_per_multiprocessor, gpu_name_c)

      gpu_name = c_char_array_to_fort_str(GPU_NAME_LEN, gpu_name_c)

      if (cuda_major < 2) then
        write(info_str,'(2X,2A)') "*** Error: GPU does not", &
             " support double arithmetic!"
        call aims_stop (info_str, func)
      endif
    else
      gpu_name = "N/A"
    end if

    ! Output versioning information and exit
    write(info_str,'(2X,A)')
    call localorb_info( info_str )
    write (info_str,'(A,I9,A,L1,A,A)') "  | MPI task ", myid, &
         ", GPU activated = ", use_gpu, &
         ", Device Name = ", trim(gpu_name)
    call localorb_allinfo(info_str)

    call get_times( time_gpu_init,     clock_time_gpu_init, &
                    tot_time_gpu_init, tot_clock_time_gpu_init, .true. )
  end subroutine initialize_gpu
!******

  subroutine output_mem_array_gpu(name_array, n_elements, bytes_data_type)
    use localorb_io, only: localorb_info, OL_norm, use_unit
    implicit none

    character(len=*), intent(in) :: name_array
    integer, intent(in) :: n_elements
    integer, intent(in) :: bytes_data_type

    character*400 :: info_str

    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ', dble(n_elements) * dble(bytes_data_type) / 1.0d6, &
         ' MB on GPU for ', trim(adjustl(name_array))
    call localorb_info ( info_str, use_unit,'(2X,A)', OL_norm  )
  end subroutine output_mem_array_gpu

!****s* FHI-aims/finalize_gpu
!  NAME
!    finalize_gpu
!  SYNOPSIS
  subroutine finalize_gpu()
!  PURPOSE
!    Tear down the GPU and uninitialize relevant variables.
!  USES
    use mpi_tasks, only: aims_stop
    use runtime_choices, only: use_gpu
    implicit none
!  ARGUMENTS
!    o None
!  INPUTS
!    o None
!  OUTPUT
!    o None
!  AUTHOR
!    William Huhn (Duke University)
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    January 2018 - Created.
!  NOTES
!    The cuda_* subroutines called in this subroutine are found in the C CUDA
!    code.
!  SOURCE
    character(*), parameter :: func = 'finalize_gpu'

    if (.not.use_gpu) then
      call aims_stop("GPU is not in use, can't deinitialize them.",func)
    end if

    ! Turn off the GPUs
    call finalize_cuda_and_cublas()

    ! And reset most relevant variables to their initial values.
    gpu_density_used = .false.
    gpu_hamiltonian_used = .false.
    gpu_forces_used = .false.

    n_gpus_avail = -1
    cuda_major = -1
    cuda_minor = -1
    if (allocated(gpu_name)) deallocate(gpu_name)
    total_global_gpu_mem = -1
    total_constant_gpu_mem = -1
    shared_gpu_mem_per_block = -1
    gpu_clock_rate = -1
    gpu_multiprocessor_count = -1
    gpu_max_threads_per_multiprocessor = -1
  end subroutine finalize_gpu
!******

end module
!******
