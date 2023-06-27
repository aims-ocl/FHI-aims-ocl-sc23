!****s* FHI-aims/integrate_first_order_rho_polar_reduce_memory
!  NAME
!   integrate_first_order_rho_polar_reduce_memory
!  SYNOPSIS

subroutine integrate_first_order_rho_polar_reduce_memory &
     (partition_tab_std, basis_l_max, first_order_density_matrix, &
      first_order_rho, first_order_density_matrix_size)

!  PURPOSE
!  
!  calculate the sum(u,v){DM(1)*phi_u*phu_v},
!  using a fixed basis set. 

!  called by a SCF subroutine

!  shanghui 2012.06.28
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data, only: species_name
  use load_balancing
  use pbc_lists
  use scalapack_wrapper
  use mpi_tasks
  use opencl_util
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: partition_tab_std
  integer, dimension(n_species), intent(in)          :: basis_l_max 
  !shanghui------------------------------------------------------------------------
  !real*8,  dimension(n_hamiltonian_matrix_size), intent(in) :: first_order_density_matrix
  real*8,  intent(in) :: first_order_density_matrix(*)
  real*8,  dimension(n_full_points), intent(out) ::  first_order_rho 
  integer :: first_order_density_matrix_size
  !shanghui------------------------------------------------------------------------

!  INPUTS
!  o partition_tab_std -- values of partition functions
!  o basis_l_max -- maximum l of basis functions.
!  o first_order_DM -- first-order density matrix
!
!  OUTPUT
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




  !  local variables

  real*8, dimension(n_spin) :: local_potential_parts

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab


  real*8 coord_current(3)

!  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
!  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)

  real*8,dimension(:,:),allocatable:: dist_tab
  real*8,dimension(:,:),allocatable:: dist_tab_sq

  real*8 i_r(n_max_compute_atoms)

!  real*8 dir_tab(3,n_centers_integrals, n_max_batch_size)
  real*8, dimension(:,:,:),allocatable:: dir_tab


  real*8 trigonom_tab(4,n_max_compute_atoms)

  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: kinetic_wave
  real*8,dimension(:,:)  ,allocatable:: wave

 !------------shanghui add for first_order_rho-----------------
  real*8,dimension(:), allocatable :: local_first_order_rho
  real*8, dimension(:,:),  allocatable :: first_order_density_matrix_con 
  ! wyj add for prune manually
  real*8, dimension(:),  allocatable :: my_first_order_density_matrix_con
  ! wyj: add tmp array for first_order_rho
  real*8, pointer :: first_order_rho_tmp(:)
 !------------shanghui end add for first_order_rho-------------


  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)




  !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

  integer :: n_compute_c, n_compute_a
!  integer :: i_basis(n_centers_basis_I)
  integer,dimension(:),allocatable :: i_basis

  integer :: n_compute_fns

!  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
!  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
!  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer,dimension(:),  allocatable :: i_basis_fns
  integer,dimension(:,:),allocatable :: i_basis_fns_inv
  integer,dimension(:),  allocatable :: i_atom_fns

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms)

  ! indices for basis functions that are nonzero at current point

  integer :: rad_index(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_ham)
  integer :: l_index(n_max_compute_fns_ham)
  integer :: l_count(n_max_compute_fns_ham)
  integer :: fn_atom(n_max_compute_fns_ham)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute
  integer :: zero_index_point(n_max_compute_ham)

  ! active atoms in current batch
  integer :: n_batch_centers
  integer :: batch_center(n_centers_integrals)

  !     for splitting of angular shells into "octants"

  integer division_low
  integer division_high

  !  counters

  integer i_basis_1
  integer i_basis_2
  integer i_atom, i_atom_2
  integer i_grid
  integer i_index, i_l, i_m
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_DM_rho
  integer :: i_full_points_2

  integer :: i_spin
  character*200 :: info_str

  integer :: i_my_batch
  integer :: i_my_batch2

  integer :: i_radial, i_angular, info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

!  integer ld_hamiltonian  ! leading dimension of hamiltonian in calling routine

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i_off, i, j, n_bp
  integer, allocatable :: ins_idx(:)
  integer, allocatable :: ins_idx_all_batches(:,:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all
  real*8 time_f(12)
  real*8 time2, time3

  integer, dimension(:), allocatable :: n_batch_centers_all_batches
  integer, dimension(:,:), allocatable :: batch_center_all_batches
  integer, dimension(:,:), allocatable :: batch_center_all_batches_for_copy
  integer, dimension(:), allocatable :: n_points_all_batches
  integer, dimension(:,:), allocatable :: batch_point_to_i_full_point

  integer :: mpi_buffer
  integer :: tag, count
  integer :: ierr, status(MPI_STATUS_SIZE)
  real*8 :: gemm_flop
  real*8 :: time_evaluate, time1

  integer :: n_basis_local_actual = 0

  ! begin work

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating first_order_rho: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating first_order_rho: batch-based integration."
  endif
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

  ! begin with general allocations

  allocate(dist_tab(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab                      ')

  allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab_sq                   ')

  allocate(dir_tab(3,n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dir_tab                       ')

  allocate(i_basis_fns(n_basis_fns*n_centers_integrals), stat=info)
  call check_allocation(info, 'i_basis_fns                   ')

  allocate(i_basis_fns_inv(n_basis_fns,n_centers), stat=info)
  call check_allocation(info, 'i_basis_fns_inv               ')

  allocate(i_atom_fns(n_basis_fns*n_centers_integrals),stat=info)
  call check_allocation(info, 'i_atom_fns                    ')


  l_ylm_max = l_wave_max


  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info )
  call check_allocation(info, 'index_lm                      ')


  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(kinetic_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'kinetic_wave                  ')

  allocate(wave(n_max_compute_ham, n_max_batch_size), STAT=info )
  call check_allocation(info, 'wave                          ')

 !------------shanghui add for first_order_rho-----------------
  
  allocate(local_first_order_rho(n_max_batch_size), STAT=info ) 
  call check_allocation(info, 'local_first_order_rho   ') 

 ! if(.not. allocated(first_order_density_matrix_con))then
 ! allocate(first_order_density_matrix_con(n_max_compute_dens,n_max_compute_dens),stat=info)
 ! call check_allocation(info, 'first_order_density_matrix_con            ')
 ! end if
 ! wyj: allocate memory 
 !if (use_local_index .and. use_load_balancing) then
 if (n_bp > 0) then
     if(.not. allocated(my_first_order_density_matrix_con))then
         allocate(my_first_order_density_matrix_con(n_max_compute_dens*n_max_compute_dens),stat=info)
         call check_allocation(info, 'my_first_order_density_matrix_con            ')
     end if
 else
     if(.not. allocated(first_order_density_matrix_con))then
         allocate(first_order_density_matrix_con(n_max_compute_dens,n_max_compute_dens),stat=info)
         !allocate(first_order_density_matrix_con_check(n_max_compute_dens,n_max_compute_dens),stat=info)
         call check_allocation(info, 'first_order_density_matrix_con            ')
     end if
 endif
 !------------shanghui end add for first_order_rho------------- 


  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')


  first_order_rho = 0.0d0



  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab


    n_basis_local_actual = batch_perm(n_bp)%n_basis_local
    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))
    allocate(ins_idx_all_batches(batch_perm(n_bp)%n_basis_local, n_my_batches_work))

    ! wyj:  tmp array 
    allocate(first_order_rho_tmp(batch_perm(n_bp)%n_full_points))
    first_order_rho_tmp = 0.0d0

    !call init_comm_full_local_matrix_scalapack(&
    !    batch_perm(n_bp)%n_basis_local, &
    !    batch_perm(n_bp)%i_basis_local )

    !ld_hamiltonian = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

    ! wyj: add for rho_tmp
    allocate(first_order_rho_tmp(n_full_points))
    first_order_rho_tmp = 0.0d0
    !ld_hamiltonian = n_hamiltonian_matrix_size

  endif

  position_in_hamiltonian_dim1 = size(position_in_hamiltonian,dim=1)
  position_in_hamiltonian_dim2 = size(position_in_hamiltonian,dim=2)
  column_index_hamiltonian_size = size(column_index_hamiltonian)

  if(allocated(batches_size_rho) .or. (n_my_batches_work .ne. n_my_batches_work_rho)) then
    if(allocated(batches_size_rho) .and. (n_my_batches_work .ne. n_my_batches_work_rho)) then
      deallocate(batches_size_rho)
      deallocate(batches_points_coords_rho)
      deallocate(batches_batch_n_compute_rho)
      deallocate(batches_batch_i_basis_rho)
    endif
    if(.not. allocated(batches_size_rho)) then
      allocate(batches_size_rho(n_my_batches_work), stat=info)
      call check_allocation(info, 'batches_size_rho')
      allocate(batches_points_coords_rho(3, n_max_batch_size, n_my_batches_work), stat=info)
      call check_allocation(info, 'batches_points_coords_rho')
      allocate(batches_batch_n_compute_rho(n_my_batches_work), stat=info)
      call check_allocation(info, 'batches_batch_n_compute_rho')
      ! allocate(batches_batch_i_basis_rho(n_centers_basis_I, n_my_batches_work), stat=info)
      allocate(batches_batch_i_basis_rho(n_max_compute_dens, n_my_batches_work), stat=info)
      call check_allocation(info, 'batches_batch_i_basis_rho')
    endif
    if(n_my_batches_work .ne. n_my_batches_work_rho) then
      do i_my_batch = 1, n_my_batches_work, 1
        batches_size_rho(i_my_batch) = batches_work(i_my_batch)%size
      end do
    endif
  endif

  n_my_batches_work_rho = n_my_batches_work
  n_full_points_work_rho = n_full_points
  ! if(load_balance_finished .and. .not. allocated(wave_batches_rho)) then
  !    allocate(wave_batches_rho(n_max_compute_ham, n_max_batch_size, n_my_batches_work_h), STAT=info )
  !    wave_batches_rho = 0
  !    call check_allocation(info, 'wave_batches_rho')
  ! endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  allocate(n_points_all_batches(n_my_batches_work),stat=info)
  call check_allocation(info, 'n_points_all_batches                 ')
  allocate(n_batch_centers_all_batches(n_my_batches_work),stat=info)
  call check_allocation(info, 'n_batch_centers_all_batches          ')
  allocate(batch_center_all_batches(max_n_batch_centers, n_my_batches_work),stat=info)
  call check_allocation(info, 'batch_center_all_batches             ')
  allocate(batch_center_all_batches_for_copy(max_n_batch_centers, n_my_batches_work),stat=info)
  call check_allocation(info, 'batch_center_all_batches_for_copy             ')
  allocate(batch_point_to_i_full_point(n_max_batch_size, n_my_batches_work),stat=info)
  call check_allocation(info, 'batch_point_to_i_full_point          ')

  gemm_flop = 0
  i_full_points_2 = 0
  do i_my_batch = 1, n_my_batches_work, 1
   n_compute_c = 0
   n_compute_a = 0
   i_basis = 0
   i_point = 0
   do i_index = 1, batches_work(i_my_batch)%size, 1
      i_full_points_2 = i_full_points_2 + 1
      if (partition_tab(i_full_points_2).gt.0.d0) then
         i_point = i_point+1
         batch_point_to_i_full_point(i_point, i_my_batch) = i_full_points_2
        !  if(i_my_batch .le. 3) print*, i_my_batch, i_point, batch_point_to_i_full_point(i_point, i_my_batch), i_full_points_2
   !       points_match_all_batches(i_point, i_my_batch) = i_full_points_2
   !       ! get current integration point coordinate
         ! TODO 注意这里 ！！！！！ 第二维可能是 i_point 也可能是 i_index, 要与 .c/.cl 适配 ！！！！
         batches_points_coords_rho(:,i_point,i_my_batch) = batches_work(i_my_batch) % points(i_index) % coords(:)
         if(n_periodic > 0)then
            call map_to_center_cell(batches_points_coords_rho(:,i_point, i_my_batch) )
         end if
      end if
   enddo
   n_points_all_batches(i_my_batch) = i_point
   if (prune_basis_once) then
     n_compute_c = batches_work(i_my_batch)%batch_n_compute
     i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
     if(n_bp > 0) then
      do i=1,n_compute_c
        ins_idx_all_batches(i, i_my_batch) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
      enddo
     endif
   end if
   gemm_flop = gemm_flop + 2 * i_point * n_compute_c * n_compute_c
   batches_batch_n_compute_rho(i_my_batch) = n_compute_c
   batches_batch_i_basis_rho(1:n_compute_c, i_my_batch) = i_basis(1:n_compute_c)
   call collect_batch_centers_p2 &
   ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, &
   n_batch_centers_all_batches(i_my_batch), batch_center &
   )

   ! if the n_batch_centers is lager than max_n_batch_centers, the resize the arrays and copy
   if (n_batch_centers_all_batches(i_my_batch) .ge. max_n_batch_centers) then
    ! new, with factor = 1.2
    max_n_batch_centers = (n_batch_centers_all_batches(i_my_batch)+1) * 1.2d0 + 4

    do i_my_batch2 = 1, i_my_batch-1, 1
      batch_center_all_batches_for_copy(:, i_my_batch2) = batch_center_all_batches(:, i_my_batch2)
    enddo

    deallocate(batch_center_all_batches)
    allocate(batch_center_all_batches(max_n_batch_centers, n_my_batches_work),stat=info)
    call check_allocation(info, 'batch_center_all_batches             ')

    do i_my_batch2 = 1, i_my_batch-1, 1
      batch_center_all_batches(:, i_my_batch2) = batch_center_all_batches_for_copy(:, i_my_batch2)
    enddo

    deallocate(batch_center_all_batches_for_copy)
    allocate(batch_center_all_batches_for_copy(max_n_batch_centers, n_my_batches_work),stat=info)
    call check_allocation(info, 'batch_center_all_batches_for_copy             ')
   endif

   ! only copy the batch_center that will be used latter
   batch_center_all_batches(1:n_batch_centers_all_batches(i_my_batch), i_my_batch) = batch_center(1:n_batch_centers_all_batches(i_my_batch))


  end do

  !-----------------------------------------------------------------------------


  ! initialize


  i_basis_fns_inv = 0

  ! initialize index_lm

  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo


  i_full_points = 0
  i_full_points_DM_rho = 0
  i_full_points_2 = 0

  ! perform partitioned integration, batch by batch of integration point.
  ! This will be the outermost loop, to save evaluations of the potential.
  ! and the Y_lm functions

  call m_save_load_not_count()

  !print *, 'cpscf_rho=', myid
  if(myid .eq. 0) print*, "n_bp= ", n_bp
  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  time_evaluate = 0



! if(myid .eq. 0) then
  call rho_pass_vars( &
    l_ylm_max, batch_perm(n_bp)%n_local_matrix_size, n_basis_local_actual,&
    batch_perm(n_bp)%n_full_points, first_order_density_matrix_size,&
    basis_l_max, n_points_all_batches,&
    n_batch_centers_all_batches, batch_center_all_batches,&
    batch_point_to_i_full_point,&
    ins_idx_all_batches,&
    first_order_rho_tmp, first_order_density_matrix,&
    partition_tab(1))
! endif

if((use_c_version .or. use_opencl_version) .and. opencl_rho_fortran_init .and. use_rho_c_cl_version .and. opencl_util_init) then
  tag = 10203
  count = 1
  ! mpi_platform_relative_id = myid / mpi_platform_num
  ! mpi_platform_relative_id = mod(myid, 32)

  if(mod(mpi_platform_relative_id, mpi_task_per_gpu) .ne. 0) then
    call MPI_Recv(mpi_buffer, count, MPI_Integer, myid-1, tag, MPI_COMM_WORLD, status, ierr)
    ! print*, "myid=", myid, "recv '", mpi_buffer, "' from id=", myid-1
  endif

  call integrate_first_order_rho_sub_tmp2( &
          l_ylm_max, batch_perm(n_bp)%n_local_matrix_size, batch_perm(n_bp)%n_basis_local,&
          batch_perm(n_bp)%n_full_points, first_order_density_matrix_size,&
          basis_l_max, n_points_all_batches,&
          n_batch_centers_all_batches, batch_center_all_batches,&
          batch_point_to_i_full_point,&
          ins_idx_all_batches,&
          first_order_rho_tmp, first_order_density_matrix,&
          partition_tab(1))

  if(mod(mpi_platform_relative_id, mpi_task_per_gpu) .ne. (mpi_task_per_gpu-1) .and. myid .ne. (n_tasks - 1)) then
    mpi_buffer = 2000 + myid
    call MPI_Send(mpi_buffer, count, MPI_Integer, myid+1, tag, MPI_COMM_WORLD, ierr)
    ! print*, "myid=", myid, "send '", mpi_buffer, "' to id=", myid+1
  endif

  ! call integrate_first_order_rho_sub_tmp2( &
  !         l_ylm_max, -1, -1, n_hamiltonian_matrix_size,&
  !         basis_l_max, n_points_all_batches,&
  !         n_batch_centers_all_batches, batch_center_all_batches,&
  !         batch_point_to_i_full_point,&
  !         -1,&
  !         first_order_rho_tmp, first_order_density_matrix,&
  !         partition_tab(1))
  ! 注意是 first_order_rho_tmp 而非 first_order_rho !!!
  ! 拿 n_hamiltonian_matrix_size 作为 first_order_density_matrix 的大小，尽管没标

  ! call m_save_check_rho(first_order_rho)
else

  time_f(:) = 0

  do i_my_batch = 1, n_my_batches_work, 1
     if(get_batch_weights) time_start = mpi_wtime()

     n_compute_c = 0
     n_compute_a = 0
     i_basis = 0

     i_point = 0

     time2 = mpi_wtime()

     ! loop over one batch
     do i_index = 1, batches_work(i_my_batch)%size, 1

        i_full_points_2 = i_full_points_2 + 1

        if (partition_tab(i_full_points_2).gt.0.d0) then

           i_point = i_point+1

           ! get current integration point coordinate
           coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

           if(n_periodic > 0)then
              call map_to_center_cell(coord_current(1:3) )
           end if

           ! compute atom-centered coordinates of current integration point,
           ! as viewed from all atoms
           call tab_atom_centered_coords_p0 &
                ( coord_current,  &
                dist_tab_sq(1,i_point),  &
                dir_tab(1,1,i_point), &
                n_centers_integrals, centers_basis_integrals )

           ! determine which basis functions are relevant at current integration point,
           ! and tabulate their indices

           ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
           if (.not.prune_basis_once) then
              call prune_basis_p2 &
                   ( dist_tab_sq(1,i_point), &
                   n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
           endif

        end if
     enddo

     time_f(1) = time_f(1) + (mpi_wtime() - time2)
     time2 = mpi_wtime()
     if (prune_basis_once) then
        n_compute_c = batches_work(i_my_batch)%batch_n_compute
        i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
        ! wyj: use prune manually
        if(n_bp > 0) then
            do i=1,n_compute_c
                ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
            enddo
            my_first_order_density_matrix_con(1:n_compute_c*n_compute_c) = 0
            !first_order_density_matrix_con(1:n_compute_c, 1:n_compute_c) = 0
            do i=1,n_compute_c
                i_off = (ins_idx(i)*(ins_idx(i)-1))/2
                do j=1,i
                    if (ins_idx(j) + i_off > batch_perm(n_bp)%n_local_matrix_size) stop
                    my_first_order_density_matrix_con(j+(i-1)*n_compute_c) = first_order_density_matrix(ins_idx(j)+i_off)
                    my_first_order_density_matrix_con(i+(j-1)*n_compute_c) = first_order_density_matrix(ins_idx(j)+i_off)
                    !first_order_density_matrix_con(j, (i-1)) = first_order_density_matrix_sparse(ins_idx(j)+i_off)
                enddo
            enddo
        endif
     end if
     !write(use_unit,*) prune_basis_once,i_my_batch,n_compute_c,i_basis(1:n_compute_c)
     !-----------shanghui begin test prune_matrix here-------------------------
     !if (.not. (use_local_index .and. use_load_balancing)) then
     if (n_bp <= 0) then
         call  prune_density_matrix_sparse_polar_reduce_memory(first_order_density_matrix, & 
             first_order_density_matrix_con, &
             n_compute_c, i_basis)  
     endif

     time_f(2) = time_f(2) + (mpi_wtime() - time2)
     time2 = mpi_wtime()
     !-----------shanghui end test prune_matrix here-------------------------

     ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
     call collect_batch_centers_p2 &
     ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, &
       n_batch_centers, batch_center &
     )

     time_f(3) = time_f(3) + (mpi_wtime() - time2)

     n_points = i_point
     !------shanghui make all wave to 0.0d0---------------
     !wave=0.0d0
     !gradient_basis_wave_npoints=0.0d0
     !------shanghui end make all wave to 0.0d0----------


     ! Perform actual integration if more than 0 basis functions
     ! are actually relevant on the present angular shell ...
     if (n_compute_c.gt.0) then

        n_rel_points = 0
        i_point = 0

        ! loop over one batch of integration points
        do i_index = 1, batches_work(i_my_batch)%size, 1

           ! Increment the (global) counter for the grid, to access storage arrays
           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1

              coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)!SAG


              ! for all integrations
              partition(i_point) = partition_tab(i_full_points)


              n_compute_atoms = 0
              n_compute_fns = 0

              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
              ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
              ! without any copying and without doing any unnecessary operations.
              ! The price is that the interface is no longer explicit in terms of physical
              ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
              time2 = mpi_wtime()
              call prune_radial_basis_p2 &
                   ( n_max_compute_atoms, n_max_compute_fns_ham, &
                     dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
                     n_compute_atoms, atom_index, atom_index_inv, &
                     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                     i_atom_fns, spline_array_start, spline_array_end, &
                     n_centers_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                     fn_atom, n_zero_compute, zero_index_point &
                    )

              time_f(4) = time_f(4) + (mpi_wtime() - time2)
              time2 = mpi_wtime()  
              ! Tabulate distances, unit vectors, and inverse logarithmic grid units
              ! for all atoms which are actually relevant
              call tab_local_geometry_p2 &
                   ( n_compute_atoms, atom_index, &
                     dist_tab(1,i_point), i_r )
              time_f(5) = time_f(5) + (mpi_wtime() - time2)
              time2 = mpi_wtime()
              ! compute trigonometric functions of spherical coordinate angles
              ! of current integration point, viewed from all atoms
              call tab_trigonom_p0 &
                   ( n_compute_atoms, dir_tab(1,1,i_point), trigonom_tab )
              time_f(6) = time_f(6) + (mpi_wtime() - time2)
              time2 = mpi_wtime()

                call tab_wave_ylm_p0 &
                   ( n_compute_atoms, atom_index,  &
                   trigonom_tab, basis_l_max,  &
                   l_ylm_max, ylm_tab )
              time_f(7) = time_f(7) + (mpi_wtime() - time2)
              time2 = mpi_wtime()
              ! Now evaluate radial functions
              ! from the previously stored compressed spline arrays
              call evaluate_radial_functions_p0  &
                   (   spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_wave_ordered, radial_wave,  &
                   .false. , n_compute_c, n_max_compute_fns_ham )
              time_f(8) = time_f(8) + (mpi_wtime() - time2)
              time2 = mpi_wtime()
              ! tabulate total wave function value for each basis function
              call evaluate_waves_p2  &
                   ( n_compute_c, n_compute_atoms, n_compute_fns, &
                     l_ylm_max, ylm_tab, one_over_dist_tab,   &
                     radial_wave, wave(1,i_point), &
                     rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point &
                   )
              time_f(9) = time_f(9) + (mpi_wtime() - time2)
              time2 = mpi_wtime()



              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0


           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch

        time1 = mpi_wtime()

        ! Now add all contributions to the full first_order_rho

        ! wyj: transfer correspond arugment
        !if (use_local_index .and. use_load_balancing) then
        if (n_bp > 0) then
            !wyj:TODO
            !my_first_order_density_matrix_con(:) = 1.0d0
            call  evaluate_first_order_rho_polar_reduce_memory( &
                n_points,     &
                n_compute_c,i_basis,& 
                wave, &
                my_first_order_density_matrix_con(1:n_compute_c*n_compute_c), & 
                local_first_order_rho)
        else
            !first_order_density_matrix_con(:,:) = 1.0d0
            call  evaluate_first_order_rho_polar_reduce_memory( &
                n_points,     &
                n_compute_c,i_basis,& 
                wave, &
                first_order_density_matrix_con, & 
                local_first_order_rho)

        endif

        time_evaluate = time_evaluate + (mpi_wtime() - time1)
        ! wave_batches_rho(:, 1:n_points, i_my_batch) = wave(:, 1:n_points)

       i_point = 0
       ! loop over one batch of integration points
       do i_index = 1, batches_work(i_my_batch)%size, 1
          ! Increment the (global) counter for the grid, to access storage arrays
          i_full_points_DM_rho = i_full_points_DM_rho + 1
          if (partition_tab(i_full_points_DM_rho).gt.0.d0) then
                 i_point = i_point+1

                 !first_order_rho(i_full_points_DM_rho)= & 
                 first_order_rho_tmp(i_full_points_DM_rho)= & 
                 local_first_order_rho(i_point) 


         endif  

       enddo
!----------------------------shanghui-------------------------------------------
        !
     else

       i_full_points = i_full_points + batches_work(i_my_batch)%size

        !first_order_rho(i_full_points_DM_rho+1:& 
        first_order_rho_tmp(i_full_points_DM_rho+1:& 
       i_full_points_DM_rho + batches_work(i_my_batch)%size) =0.0d0


       i_full_points_DM_rho= i_full_points_DM_rho + batches_work(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

     if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  end do ! end loop over batches
  if(myid .eq. 0) call m_save_check_rho(first_order_rho_tmp)
endif

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  ! if(mod(myid, (n_tasks / 8)) .eq. 0) print*, "myid=", myid, " time_work=", time_work
  if(myid .le. 10 .or. myid .eq. (n_tasks-1)) print*, "myid=", myid, " time_work=", time_work
  ! print*, "myid=", myid, " gemm_flop_dp=", gemm_flop, " time_evaluate=", time_evaluate
  ! print*, "myid=", myid, " time_f=", time_f(1:9)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  if(time_all>time_work*1.3 .and. .not.use_load_balancing) &
    info_str = trim(info_str) // ' => Consider using load balancing!'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)




  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( wave                 )) deallocate( wave                 )
  if(allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if(allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )

  if(allocated( n_points_all_batches         )) deallocate( n_points_all_batches         )
  if(allocated( n_batch_centers_all_batches  )) deallocate( n_batch_centers_all_batches  )
  if(allocated( batch_center_all_batches     )) deallocate( batch_center_all_batches     )
  if(allocated( batch_center_all_batches_for_copy     )) deallocate( batch_center_all_batches_for_copy )
  if(allocated( batch_point_to_i_full_point  )) deallocate( batch_point_to_i_full_point  )

  if(allocated( first_order_density_matrix_con   )) deallocate( first_order_density_matrix_con   )
  ! wyj: deallocate intime
  if (allocated( my_first_order_density_matrix_con   )) deallocate( my_first_order_density_matrix_con   )
  if(allocated( local_first_order_rho)) deallocate( local_first_order_rho )

  !if(use_batch_permutation > 0) then
  !  deallocate(ins_idx)
  !endif

  !if(get_batch_weights) deallocate(batch_times)

  ! wyj: merge with deallocate
  !if(get_batch_weights) call set_batch_weights(n_bp, batch_times)
  if(use_batch_permutation > 0) then
      call permute_point_array_back(n_bp, 1, first_order_rho_tmp, first_order_rho)
      deallocate(first_order_rho_tmp)
      deallocate(ins_idx)
      deallocate(ins_idx_all_batches)
  else
      first_order_rho(1:n_full_points) = first_order_rho_tmp(1:n_full_points)
      deallocate(first_order_rho_tmp)
  endif
  ! wyj: set new times
  if(get_batch_weights) then
      if (allocated(batch_times)) then
          !print *, myid, 'first_order_rho begin set_batch_weights '
          call set_batch_weights(n_bp, batch_times)
          deallocate(batch_times)
          !print *, myid, 'first_order_rho end set_batch_weights '
      endif
  endif

  if(.not. opencl_rho_fortran_init) opencl_rho_fortran_init = .true.

end subroutine integrate_first_order_rho_polar_reduce_memory
!******
