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

    logical :: use_rho_c_cl_version = .false.
    logical :: use_sumup_c_cl_version = .false.
    logical :: use_sumup_pre_c_cl_version = .false.
    logical :: use_h_c_cl_version = .false.

    logical :: use_opencl = .false.

    integer :: debug_io = 0

    ! true 
    ! false

    real*8, dimension(:,:,:), allocatable:: wave_batches_h
    real*8, dimension(:,:,:), allocatable:: wave_batches_rho

    ! integer :: mpi_platform_num = 1
    ! integer :: mpi_platform_tasks = 32
    ! integer :: mpi_platform_device_num = 4

    integer :: useless_var

    character*128  evalue
    integer :: tmpint = 0

contains

    subroutine read_opencl_settings()
        use mpi_tasks, only : myid, mpi_platform_relative_id, mpi_per_node, mpi_task_per_gpu, mpi_gpu_per_node
        use hartree_potential_storage, only : use_rho_multipole_shmem, rho_multipole_tile_size
        implicit none

        character(len=100) :: var_value ! 声明字符型变量，用于存储环境变量的值
        integer :: ierr                ! 存储GET_ENVIRONMENT_VARIABLE的返回值
        logical :: is_empty            ! 用于判断环境变量是否为空

        call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE', var_value)
        is_empty = TRIM(var_value) == ''
        if (.not. is_empty) then
            read(var_value, *) rho_multipole_tile_size
            if(myid .eq. 0) write(*, *) 'The environment variable FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE is:', rho_multipole_tile_size
        else
            if(myid .eq. 0) write(*, *) 'The environment variable FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE is empty, default to ', rho_multipole_tile_size
        endif

        call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM', var_value)
        if (TRIM(var_value) == "ON") then
            use_rho_multipole_shmem = .true.
            if(myid .eq. 0) write(*, *) 'Enable FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM.'

            call GET_ENVIRONMENT_VARIABLE('MPI_PER_NODE', var_value)
            is_empty = TRIM(var_value) == ''
            if (is_empty) then
                if(myid .eq. 0) write(*, *) 'The environment variable MPI_PER_NODE is empty, it should be set to the value of tasks-per-node or 1 if needed, to enable shmem or opencl.'
                stop
            else
                read(var_value, *) mpi_per_node
                if(myid .eq. 0) write(*, *) 'The environment variable MPI_PER_NODE is:', mpi_per_node
            endif

        else
            if(myid .eq. 0) write(*, *) 'The environment variable FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM is not ON, set to OFF.'
        endif

        call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL', var_value)
        if (TRIM(var_value) == "ON") then
            use_opencl = .true.

            use_rho_c_cl_version = .true.
            use_sumup_c_cl_version = .true.
            use_h_c_cl_version = .true.
            if(myid .eq. 0) write(*, *) 'Enable FHI_AIMS_OCL'

            call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL_DFPT_POLAR_RHO', var_value)
            if (TRIM(var_value) == "OFF") then
                use_rho_c_cl_version = .false.
                if(myid .eq. 0) write(*, *) '  Disable FHI_AIMS_OCL_DFPT_POLAR_RHO'
            else
                if(myid .eq. 0) write(*, *) '  Enable FHI_AIMS_OCL_DFPT_POLAR_RHO'
            endif
            call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL_DFPT_DIEL_SUMUP', var_value)
            if (TRIM(var_value) == "OFF") then
                use_sumup_c_cl_version = .false.
                if(myid .eq. 0) write(*, *) '  Disable FHI_AIMS_OCL_DFPT_DIEL_SUMUP'
            else
                if(myid .eq. 0) write(*, *) '  Enable FHI_AIMS_OCL_DFPT_DIEL_SUMUP'
            endif
            call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL_DFPT_POLAR_H', var_value)
            if (TRIM(var_value) == "OFF") then
                use_h_c_cl_version = .false.
                if(myid .eq. 0) write(*, *) '  Disable FHI_AIMS_OCL_DFPT_POLAR_H'
            else
                if(myid .eq. 0) write(*, *) '  Enable FHI_AIMS_OCL_DFPT_POLAR_H'
            endif


            call GET_ENVIRONMENT_VARIABLE('MPI_PER_NODE', var_value)
            is_empty = TRIM(var_value) == ''
            if (is_empty) then
                if(myid .eq. 0) write(*, *) 'The environment variable MPI_PER_NODE is empty, it should be set to the value of tasks-per-node.'
                stop
            else
                read(var_value, *) mpi_per_node
                if(myid .eq. 0) write(*, *) 'The environment variable MPI_PER_NODE is:', mpi_per_node
            endif

            call GET_ENVIRONMENT_VARIABLE('GPU_PER_NODE', var_value)
            is_empty = TRIM(var_value) == ''
            if (is_empty) then
                if(myid .eq. 0) write(*, *) 'The environment variable GPU_PER_NODE is empty, it should be set to the number of gpus used on each node.'
                stop
            else
                read(var_value, *) mpi_gpu_per_node
                if(myid .eq. 0) write(*, *) 'The environment variable GPU_PER_NODE is:', mpi_gpu_per_node
                mpi_task_per_gpu = mpi_per_node / mpi_gpu_per_node
                if(myid .eq. 0) write(*, *) 'mpi_task_per_gpu is:', mpi_task_per_gpu
                if(myid .eq. 0 .and. (mpi_task_per_gpu .lt. 0 .or. mpi_task_per_gpu .gt. 16)) then
                    write(*, *) 'mpi_task_per_gpu must <= 16'
                    stop
                endif
            endif

            call GET_ENVIRONMENT_VARIABLE('FHI_AIMS_OCL_DEBUG_IO', var_value)
            if (TRIM(var_value) == "ON") then
                debug_io = 1
                if(myid .eq. 0) write(*, *) 'Enable FHI_AIMS_OCL_DEBUG_IO.'
            else
                if(myid .eq. 0) write(*, *) 'The environment variable FHI_AIMS_OCL_DEBUG_IO is not ON, set to OFF.'
            endif

        else
            if(myid .eq. 0) write(*, *) 'The environment variable FHI_AIMS_OCL is not ON, will not use DFPT-OpenCL accelerate.'
        endif



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