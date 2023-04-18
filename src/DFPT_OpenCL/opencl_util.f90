module opencl_util
    implicit none
    logical :: use_c_version = .true.
    logical :: use_opencl_version = .true.
    integer :: opencl_util_debug = 1

    integer, dimension(:), allocatable :: batches_size_sumup              ! (n_my_batches)
    real*8, dimension(:,:,:), allocatable :: batches_points_coords_sumup  ! (3, n_max_batch_size, n_my_batches)
    integer :: n_my_batches_work_sumup
    integer :: n_full_points_work_sumup
    logical :: opencl_sumup_fortran_init = .false.
    
    integer :: max_n_batch_centers = 10

    integer :: n_my_batches_work_rho
    integer :: n_full_points_work_rho
    integer, dimension(:), allocatable :: batches_size_rho              ! (n_my_batches)
    integer, dimension(:), allocatable :: batches_batch_n_compute_rho   ! (n_my_batches)
    integer, dimension(:,:), allocatable :: batches_batch_i_basis_rho   ! (n_centers_basis_I, n_my_batches)
    real*8, dimension(:,:,:), allocatable :: batches_points_coords_rho  ! (3, n_max_batch_size, n_my_batches)
    logical :: opencl_rho_fortran_init = .false.

    integer :: n_my_batches_work_h
    integer :: n_full_points_work_h
    integer, dimension(:), allocatable :: batches_size_h                ! (n_my_batches)
    integer, dimension(:), allocatable :: batches_batch_n_compute_h     ! (n_my_batches)
    integer, dimension(:,:), allocatable :: batches_batch_i_basis_h     ! (n_centers_basis_I, n_my_batches)
    real*8, dimension(:,:,:), allocatable :: batches_points_coords_h    ! (3, n_max_batch_size, n_my_batches)
    logical :: opencl_h_fortran_init = .false.
    
    logical :: opencl_util_init = .false.
    logical :: load_balance_finished = .false.

    logical :: use_rho_c_cl_version = .true.
    logical :: use_sumup_c_cl_version = .true.
    logical :: use_sumup_pre_c_cl_version = .false.
    logical :: use_h_c_cl_version = .true.

    ! true 
    ! false

    real*8, dimension(:,:,:), allocatable:: wave_batches_h
    real*8, dimension(:,:,:), allocatable:: wave_batches_rho

    ! integer :: mpi_platform_num = 1
    ! integer :: mpi_platform_tasks = 32
    ! integer :: mpi_platform_device_num = 4
    integer :: mpi_platform_relative_id = 0     ! will be set in main.f90
    integer :: mpi_per_node = 32
    integer :: mpi_task_per_gpu = 8

    integer :: useless_var

    character*128  evalue
    integer :: tmpint = 0

contains

    subroutine read_opencl_settings()
        ! implicit none
        ! open (2, file='opencl_settings.config', status='old')
        ! read (2, *) mpi_platform_num
        ! print*, "mpi_platform_num=", mpi_platform_num
        ! close(2)

        ! if(use_rho_c_cl_version .or. use_sumup_c_cl_version .or. use_h_c_cl_version) then
        !     call getenv('MPI_PER_NODE', evalue)
        !     if(myid .eq. 0) print*, "MPI_PER_NODE: '", evalue, "'"
        !     read(evalue,"i0") tmpint
        !     if(tmpint .le. 0) then
        !         if(myid .eq. 0) print*, "Envionment variable MPI_PER_NODE should be set!"
        !         stop
        !     endif
        !     mpi_per_node = tmpint


        !     call getenv('MPI_TASK_PER_GPU', evalue)
        !     if(myid .eq. 0) print*, "MPI_TASK_PER_GPU: '", evalue, "'"
        !     read(evalue,"i0") tmpint
        !     if(tmpint .le. 0) then
        !         if(myid .eq. 0) print*, "Envionment variable MPI_TASK_PER_GPU should be set!"
        !         stop
        !     endif
        !     mpi_task_per_gpu = tmpint
        ! endif
    end subroutine read_opencl_settings

end module opencl_util