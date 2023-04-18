!****h* FHI-aims/load_balancing
!  NAME
!    load_balancing -- routines for load balancing the integrations
!  SYNOPSIS
module load_balancing
!  PURPOSE
!    This module provides routines for load balancing the integrations
!  USES
  use grids, only: batch_of_points
  implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
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

  private ! By default all variables below are private

  ! Number of batch permutation to use, 0 = no permutation
  integer, public :: use_batch_permutation
  ! Flag if batch weights should be measured
  logical, public :: get_batch_weights

  integer, parameter :: max_batch_permutations = 3

  ! load balancing pattern for integrations
  integer, parameter, public :: n_bp_integ = 1
  ! load balancing pattern for Hartree potential
  integer, parameter, public :: n_bp_hpot  = 2
  ! load balancing pattern for density update
  integer, parameter, public :: n_bp_density = 3
  ! weights of original batches for distribution
  real*8, allocatable :: batch_weight(:)
  ! Offset of my original batches within global batches
  integer :: my_batch_off

  public batch_permutation
  type batch_permutation

    logical :: initialized

    integer :: n_my_batches  ! My number of batches in batch permutation
    integer :: n_full_points ! My number of points in batch permutation

    type (batch_of_points), pointer :: batches(:) ! Permuted batches

    ! Owners in permutation for original batches
    integer, pointer :: perm_batch_owner(:)

    ! Arrays for permuting the points
    integer, pointer :: point_send_cnt(:), point_send_off(:)
    integer, pointer :: point_recv_cnt(:), point_recv_off(:)

    ! Related to basis functions of permuted batches
    integer :: n_basis_local
    integer :: n_local_matrix_size
    integer, pointer :: i_basis_local(:)
    integer, pointer :: i_basis_glb_to_loc(:)

    ! Since partition_tab doesn't change and is needed almost everywhere,
    ! we keep a copy of the permuted partition_tab here

    real*8, pointer :: partition_tab(:)

  end type

  type(batch_permutation), public, target :: &
     batch_perm(max_batch_permutations)

  ! public routines

  public :: compute_balanced_batch_distribution
  public :: set_batch_weights
  public :: reset_batch_permutation
  public :: reset_load_balancing
  public :: permute_point_array
  public :: permute_point_array_back

  public :: set_full_local_ovlp
  public :: set_full_local_ham
  public :: get_full_local_matrix
  public :: init_comm_full_local_matrix

  public :: print_sparse_to_dense_local_index_cpscf ! wyj add for CPSCF
  public :: print_matrix_real8_2d ! wyj add for debug

   interface permute_point_array
      module procedure permute_point_array_real8_2d
      module procedure permute_point_array_real8_1d
      module procedure permute_point_array_real8_1d_legacy
      module procedure permute_point_array_real8_3d_legacy
      module procedure permute_point_array_int_2d
      module procedure permute_point_array_int_1d
   end interface
   interface permute_point_array_back
      module procedure permute_point_array_back_real8_2d
      module procedure permute_point_array_back_real8_1d
      module procedure permute_point_array_back_real8_1d_legacy
      module procedure permute_point_array_back_real8_3d_legacy
      module procedure permute_point_array_back_int_2d
      module procedure permute_point_array_back_int_1d
   end interface

contains

  !-----------------------------------------------------------------------------
  !****s* load_balancing/compute_balanced_batch_distribution
  !  NAME
  !    compute_balanced_batch_distribution
  !  SYNOPSIS

  subroutine compute_balanced_batch_distribution(n_bp)

    !  PURPOSE
    !    Computes a new batch distribution for given batch weights
    !  USES
    !    none
    use dimensions, only: n_centers_basis_I, n_centers_basis_T, n_full_points, &
        n_grid_batches, n_my_batches, n_periodic
    use runtime_choices, only: output_level, prune_basis_once
    use physics, only: partition_tab
    use synchronize_mpi_basic, only: sync_vector, sync_integer_vector
    use mpi_utilities, only: distribute_batches_by_location
    use mpi_tasks
    use grids, only: batches
    use localorb_io, only: use_unit
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !  OUTPUTS
    !    o none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i, j, n, nt, my_off, my_len, n_compute, mpierr, mpi_status(mpi_status_size)
    integer num_batches(0:n_tasks-1)
    integer, allocatable :: batch_sizes(:), batch_owner(:), i_basis(:), send_req(:), orig_batch_owner(:)
    real*8 center(3), coord_current(3)
    real*8, allocatable :: batch_desc(:,:), tmp(:,:), tmp_perm(:,:)
    logical, allocatable :: have_basis_local(:)

    type (batch_of_points), pointer :: p_batches(:)
    real*8 t_start

    character(*), parameter :: func = 'compute_balanced_batch_distribution'

    call mpi_barrier(mpi_comm_global, mpierr) ! Just for timing
    t_start = mpi_wtime()

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
       write(use_unit,"(2X,A)") "Start compute_balanced_batch_distribution"

    ! Check n_bp

    if(n_bp<=0 .or. n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if(batch_perm(n_bp)%initialized) &
      call aims_stop('value for n_bp already in use', func)

    if(.not.allocated(batch_weight)) &
      call aims_stop('batch_weight not set', func)

    !-------------------------------------------------------------------
    ! Calculate new batch distribution
    !-------------------------------------------------------------------

    ! Get the number of batches on all MPI tasks

    call mpi_allgather(n_my_batches, 1, MPI_INTEGER,num_batches, &
                       1, MPI_INTEGER, mpi_comm_global, mpierr)
    if(sum(num_batches) /= n_grid_batches) &
       ! Must never happen
       call aims_stop ('n_grid_batches does not match number of batches', func)
    ! Get my offset within global batch list

    my_batch_off = 0
    do n=0,myid-1
      my_batch_off = my_batch_off + num_batches(n)
    enddo

    ! Set up batch description with my contribution

    allocate(batch_desc(5,n_grid_batches))
    batch_desc(:,:) = 0

    do i=1,n_my_batches

      ! Get center of batch
      center(1:3) = 0.
      do j=1,batches(i)%size
        coord_current(:) = batches(i)%points(j)%coords(1:3)
        if(n_periodic > 0) call map_to_center_cell(coord_current)
        center(1:3) = center(1:3) + coord_current(1:3)
      enddo
      center(1:3) = center(1:3) / batches(i)%size

      batch_desc(1,i+my_batch_off) = center(1)
      batch_desc(2,i+my_batch_off) = center(2)
      batch_desc(3,i+my_batch_off) = center(3)
      batch_desc(4,i+my_batch_off) = batch_weight(i)
      batch_desc(5,i+my_batch_off) = batches(i)%size ! for getting the global sizes only
    enddo

    ! Get complete batch_desc

    call sync_vector(batch_desc, 5*n_grid_batches)

    ! Get batch sizes of all batches in global batch list
    allocate(batch_sizes(n_grid_batches))
    batch_sizes(:) = batch_desc(5,:)

    ! Reset batch_desc(5,:) for the call of distribute_batches_by_location

    do i=1,n_grid_batches
      batch_desc(5,i) = i
    enddo

    ! Distribute batches according to batch location
    my_off = 0
    my_len = n_grid_batches
    call distribute_batches_by_location(batch_desc, 0, n_tasks-1, my_off, my_len)

    ! Number of permuted batches for my task

    batch_perm(n_bp)%n_my_batches = my_len

    ! Set new batch owners for all the original batches

    allocate(batch_owner(n_grid_batches))
    batch_owner(:) = 0

    do i = my_off+1, my_off+my_len
      batch_owner(nint(batch_desc(5,i))) = myid
    end do

    ! batch_desc is not needed any more; since it is rather big, deallocate immediatly
    deallocate(batch_desc)

    call sync_integer_vector(batch_owner, n_grid_batches)

    !-------------------------------------------------------------------
    ! Initialize communication patterns old/new distribution
    !-------------------------------------------------------------------

    ! Set the new owners of my batches in the batch permutation

    allocate(batch_perm(n_bp)%perm_batch_owner(n_my_batches))
    do i=1,n_my_batches
      batch_perm(n_bp)%perm_batch_owner(i) = batch_owner(i+my_batch_off)
    enddo

    ! Get send counts/offsets for permuting the point arrays

    allocate(batch_perm(n_bp)%point_send_cnt(0:n_tasks-1))
    allocate(batch_perm(n_bp)%point_send_off(0:n_tasks-1))

    batch_perm(n_bp)%point_send_cnt(:) = 0
    do i=1,n_my_batches
      n = batch_perm(n_bp)%perm_batch_owner(i)
      batch_perm(n_bp)%point_send_cnt(n) = batch_perm(n_bp)%point_send_cnt(n) + batches(i)%size
    enddo
    batch_perm(n_bp)%point_send_off(0) = 0
    do n=1,n_tasks-1
      batch_perm(n_bp)%point_send_off(n) = batch_perm(n_bp)%point_send_off(n-1) + batch_perm(n_bp)%point_send_cnt(n-1)
    enddo

    ! Get recv counts/offsets for permuting the point arrays

    allocate(batch_perm(n_bp)%point_recv_cnt(0:n_tasks-1))
    allocate(batch_perm(n_bp)%point_recv_off(0:n_tasks-1))

    call mpi_alltoall(batch_perm(n_bp)%point_send_cnt, 1, MPI_INTEGER, batch_perm(n_bp)%point_recv_cnt, 1, MPI_INTEGER, &
                      mpi_comm_global, mpierr)
    batch_perm(n_bp)%point_recv_off(0) = 0
    do n=1,n_tasks-1
      batch_perm(n_bp)%point_recv_off(n) = batch_perm(n_bp)%point_recv_off(n-1) + batch_perm(n_bp)%point_recv_cnt(n-1)
    enddo

    !-------------------------------------------------------------------
    ! allocate and initialize permuted batches in batch_perm(n_bp)%batches
    !-------------------------------------------------------------------

    allocate(batch_perm(n_bp)%batches(batch_perm(n_bp)%n_my_batches))
    allocate(orig_batch_owner(batch_perm(n_bp)%n_my_batches))

    ! Please note: The following ponter assignment simplifies the code AND
    ! is necessary for the BG/P compiler which produces completely broken code otherways:

    p_batches => batch_perm(n_bp)%batches

    i = 0 ! counts all batches
    n = 0 ! counts my batches
    batch_perm(n_bp)%n_full_points = 0
    do nt = 0, n_tasks-1
      do j=1,num_batches(nt)
        i = i+1
        if(batch_owner(i) == myid) then
          n = n+1
          p_batches(n)%size = batch_sizes(i)
          allocate(p_batches(n)%points(batch_sizes(i)))
          batch_perm(n_bp)%n_full_points = batch_perm(n_bp)%n_full_points + batch_sizes(i)
          orig_batch_owner(n) = nt
        endif
      enddo
    enddo

    ! batch_sizes/batch_owner are not needed any more; since rather big, deallocate immediatly

    deallocate(batch_sizes)
    deallocate(batch_owner)

    !-------------------------------------------------------------------
    ! Set the points structure of the permuted batches from originals
    !-------------------------------------------------------------------

    batch_perm(n_bp)%initialized = .true. ! must be set before calling permute_point_array

    allocate(tmp(3,n_full_points),tmp_perm(3,batch_perm(n_bp)%n_full_points))

    n = 0
    do i=1,n_my_batches
      do j=1,batches(i)%size
        n = n+1
        tmp(:,n) = batches(i)%points(j)%coords(:)
      enddo
    enddo

    call permute_point_array(n_bp,3,tmp,tmp_perm)

    n = 0
    do i=1,batch_perm(n_bp)%n_my_batches
      do j=1,p_batches(i)%size
        n = n+1
        p_batches(i)%points(j)%coords(:) = tmp_perm(:,n)
      enddo
    enddo

    ! The following might not be necessary

    n = 0
    do i=1,n_my_batches
      do j=1,batches(i)%size
        n = n+1
        tmp(1,n) = batches(i)%points(j)%index_atom
        tmp(2,n) = batches(i)%points(j)%index_radial
        tmp(3,n) = batches(i)%points(j)%index_angular
      enddo
    enddo

    call permute_point_array(n_bp,3,tmp,tmp_perm)

    n = 0
    do i=1,batch_perm(n_bp)%n_my_batches
      do j=1,p_batches(i)%size
        n = n+1
        p_batches(i)%points(j)%index_atom    = tmp_perm(1,n)
        p_batches(i)%points(j)%index_radial  = tmp_perm(2,n)
        p_batches(i)%points(j)%index_angular = tmp_perm(3,n)
      enddo
    enddo

    deallocate(tmp, tmp_perm)

    !-------------------------------------------------------------------
    ! Set the basis functions of the permuted batches from originals
    !-------------------------------------------------------------------

    if(.not.prune_basis_once) &
       call aims_stop ('load balancing needs prune_basis_once', func)

    allocate(send_req(n_my_batches))
    allocate(i_basis(n_centers_basis_T))

    do i=1,n_my_batches
      call mpi_isend(batches(i)%batch_i_basis, batches(i)%batch_n_compute, MPI_INTEGER, &
                     batch_perm(n_bp)%perm_batch_owner(i), 1, mpi_comm_global, send_req(i), mpierr)
    enddo

    do i=1,batch_perm(n_bp)%n_my_batches
      call mpi_recv(i_basis, n_centers_basis_I, MPI_INTEGER, orig_batch_owner(i), 1, mpi_comm_global, mpi_status, mpierr)
      call mpi_get_count(mpi_status, MPI_INTEGER, n_compute, mpierr)
      p_batches(i)%batch_n_compute = n_compute
      allocate(p_batches(i)%batch_i_basis(n_compute))
      p_batches(i)%batch_i_basis(1:n_compute) = i_basis(1:n_compute)
    enddo

!    VB: Due to Mac OSX/ifort/OpenMPI bug, replace:
!    call mpi_waitall(n_my_batches, send_req, MPI_STATUSES_IGNORE, mpierr)
!    by:
    do i=1,n_my_batches
      call mpi_wait(send_req(i), mpi_status, mpierr)
    enddo

    deallocate(send_req)
    deallocate(i_basis)
    deallocate(orig_batch_owner)

    !-------------------------------------------------------------------
    ! Set local basis related variables
    !-------------------------------------------------------------------

    allocate(have_basis_local(n_centers_basis_T))

    have_basis_local(:) = .false.
    do i = 1, batch_perm(n_bp)%n_my_batches, 1
      n_compute = p_batches(i)%batch_n_compute
      have_basis_local(p_batches(i)%batch_i_basis(1:n_compute)) = .true.
    enddo

    n = count(have_basis_local)
    batch_perm(n_bp)%n_basis_local = n
    batch_perm(n_bp)%n_local_matrix_size =  n*(n+1)/2

    allocate(batch_perm(n_bp)%i_basis_local(batch_perm(n_bp)%n_basis_local))
    allocate(batch_perm(n_bp)%i_basis_glb_to_loc(n_centers_basis_T))
    batch_perm(n_bp)%i_basis_glb_to_loc(:) = 0

    n = 0
    do i=1,n_centers_basis_T
      if(have_basis_local(i))then
        n = n+1
        batch_perm(n_bp)%i_basis_glb_to_loc(i) = n
        batch_perm(n_bp)%i_basis_local(n) = i
      endif
    enddo

    deallocate(have_basis_local)

    !-------------------------------------------------------------------
    ! Set permuted partition_tab
    !-------------------------------------------------------------------

    allocate(batch_perm(n_bp)%partition_tab(batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,1,partition_tab,batch_perm(n_bp)%partition_tab)

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
       write(use_unit,"(2X,A,F10.3,A)") "Done load balancing, time needed: ",mpi_wtime()-t_start," s"

  end subroutine compute_balanced_batch_distribution
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/set_batch_weights
  !  NAME
  !    set_batch_weights
  !  SYNOPSIS

  subroutine set_batch_weights(n_bp, weights)

    !  PURPOSE
    !    Sets the batch weights for later evaluation in compute_balanced_batch_distribution
    !  USES
    use dimensions, only: n_my_batches
    use mpi_tasks
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp
    real*8, intent(in) :: weights(:)

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !    o weights -- weights to be set
    !  OUTPUTS
    !    o none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i, n, mpierr
    integer, allocatable :: batch_send_cnt(:), batch_send_off(:)
    integer, allocatable :: batch_recv_cnt(:), batch_recv_off(:)
    real*8, allocatable :: w_tmp(:)

    character(*), parameter :: func = 'set_batch_weights'

    ! Safety checks
    if(n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if(allocated(batch_weight)) deallocate(batch_weight)
    allocate(batch_weight(n_my_batches))

    if(n_bp>0) then

      ! weights is for permuted batches, we need it for the original ones

      if(.not. batch_perm(n_bp)%initialized) &
        call aims_stop('value for n_bp not in use', func)

      ! Build up cnt/off arrays for communicating batches back and forth
      ! Currently this is local in this routine, it could go to the
      ! batch_permutation structure if needed also somewhere else.

      ! Get send counts/offsets

      allocate(batch_send_cnt(0:n_tasks-1))
      allocate(batch_send_off(0:n_tasks-1))

      batch_send_cnt(:) = 0
      do i=1,n_my_batches
        n = batch_perm(n_bp)%perm_batch_owner(i)
        batch_send_cnt(n) = batch_send_cnt(n) + 1
      enddo
      batch_send_off(0) = 0
      do n=1,n_tasks-1
        batch_send_off(n) = batch_send_off(n-1) + batch_send_cnt(n-1)
      enddo

      ! Get recv counts/offsets

      allocate(batch_recv_cnt(0:n_tasks-1))
      allocate(batch_recv_off(0:n_tasks-1))

      call mpi_alltoall(batch_send_cnt, 1, MPI_INTEGER, batch_recv_cnt, 1, MPI_INTEGER, mpi_comm_global, mpierr)
      batch_recv_off(0) = 0
      do n=1,n_tasks-1
        batch_recv_off(n) = batch_recv_off(n-1) + batch_recv_cnt(n-1)
      enddo

      ! communicate weights back(!) to original distribution,
      ! thus send/recv are exchanged

      allocate(w_tmp(n_my_batches))

      call mpi_alltoallv(weights, batch_recv_cnt, batch_recv_off, MPI_REAL8,    &
                         w_tmp,   batch_send_cnt, batch_send_off, MPI_REAL8,    &
                         mpi_comm_global, mpierr)

      ! Sort into batch_weight. Attention: batch_send_off is destroyed here

      batch_weight(:) = 1.d300 ! Safety only
      do i=1,n_my_batches
        n = batch_perm(n_bp)%perm_batch_owner(i)
        batch_send_off(n) = batch_send_off(n) + 1
        batch_weight(i) = w_tmp(batch_send_off(n))
      enddo
      if(any(batch_weight(:)==1.d300)) &
         call aims_stop ('batch_weight not set correctly', func)

      deallocate(batch_send_cnt)
      deallocate(batch_send_off)
      deallocate(batch_recv_cnt)
      deallocate(batch_recv_off)
      deallocate(w_tmp)

    else

      do i=1,n_my_batches
        batch_weight(i) = weights(i)
      enddo

    endif

  end subroutine set_batch_weights
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/reset_batch_permutation
  !  NAME
  !    reset_batch_permutation
  !  SYNOPSIS

  subroutine reset_batch_permutation(n_bp)

    !  PURPOSE
    !    Resets a batch permutation
    !  USES
    use mpi_tasks, only: aims_stop
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !  OUTPUTS
    !    o none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i

    character(*), parameter :: func = 'reset_batch_permutation'

    ! Safety checks
    if(n_bp<=0 .or. n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if (associated(batch_perm(n_bp)%batches)) then
      do i=1,batch_perm(n_bp)%n_my_batches
        if (associated(batch_perm(n_bp)%batches(i)%points)) then
             deallocate(batch_perm(n_bp)%batches(i)%points)
             nullify(batch_perm(n_bp)%batches(i)%points)
        end if
        if (associated(batch_perm(n_bp)%batches(i)%batch_i_basis)) then
             deallocate(batch_perm(n_bp)%batches(i)%batch_i_basis)
             nullify(batch_perm(n_bp)%batches(i)%batch_i_basis)
        end if

        batch_perm(n_bp) % batches(i) % size = 0
        batch_perm(n_bp) % batches(i) % batch_n_compute = 0
      enddo

      deallocate(batch_perm(n_bp)%batches)
      nullify(batch_perm(n_bp)%batches)
    end if

    if (associated(batch_perm(n_bp)%perm_batch_owner)) then
         deallocate(batch_perm(n_bp)%perm_batch_owner)
         nullify(batch_perm(n_bp)%perm_batch_owner)
    end if
    if (associated(batch_perm(n_bp)%point_send_cnt)) then
         deallocate(batch_perm(n_bp)%point_send_cnt)
         nullify(batch_perm(n_bp)%point_send_cnt)
    end if
    if (associated(batch_perm(n_bp)%point_send_off)) then
         deallocate(batch_perm(n_bp)%point_send_off)
         nullify(batch_perm(n_bp)%point_send_off)
    end if
    if (associated(batch_perm(n_bp)%point_recv_cnt)) then
         deallocate(batch_perm(n_bp)%point_recv_cnt)
         nullify(batch_perm(n_bp)%point_recv_cnt)
    end if
    if (associated(batch_perm(n_bp)%point_recv_off)) then
         deallocate(batch_perm(n_bp)%point_recv_off)
         nullify(batch_perm(n_bp)%point_recv_off)
    end if
    if (associated(batch_perm(n_bp)%i_basis_local)) then
         deallocate(batch_perm(n_bp)%i_basis_local)
         nullify(batch_perm(n_bp)%i_basis_local)
    end if
    if (associated(batch_perm(n_bp)%i_basis_glb_to_loc)) then
         deallocate(batch_perm(n_bp)%i_basis_glb_to_loc)
         nullify(batch_perm(n_bp)%i_basis_glb_to_loc)
    end if
    if (associated(batch_perm(n_bp)%partition_tab)) then
         deallocate(batch_perm(n_bp)%partition_tab)
         nullify(batch_perm(n_bp)%partition_tab)
    end if

    ! for safety only
    batch_perm(n_bp)%n_my_batches  = 0
    batch_perm(n_bp)%n_full_points = 0
    batch_perm(n_bp)%n_basis_local = 0
    batch_perm(n_bp)%n_local_matrix_size = 0

    batch_perm(n_bp)%initialized = .false.

  end subroutine reset_batch_permutation
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/reset_load_balancing
  !  NAME
  !    reset_load_balancing
  !  SYNOPSIS

  subroutine reset_load_balancing

    !  PURPOSE
    !    Deletes all batch permutations and anything else stored for load balanciing
    !  USES
    !    none
    !  ARGUMENTS
    !    o none
    !  INPUTS
    !    o none
    !  OUTPUTS
    !    o none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer n

    character(*), parameter :: func = 'reset_load_balancing'

    do n = 1, max_batch_permutations
      call reset_batch_permutation(n)
    enddo

    use_batch_permutation = 0
    get_batch_weights = .false.
    my_batch_off = 0
    if (allocated(batch_weight)) deallocate(batch_weight)

  end subroutine reset_load_balancing
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_real8_2d
  !  NAME
  !    permute_point_array_real8_2d
  !  SYNOPSIS

  subroutine permute_point_array_real8_2d(n_bp, ndim, arr, arr_perm)

    !  PURPOSE
    !    Permutes an array of integration points from original to permuted distribution
    !  USES
    use dimensions, only: n_full_points, n_my_batches
    use grids, only: batches
    use mpi_tasks
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    real*8, intent(in)  :: arr(ndim,n_full_points)
    real*8, intent(out) :: arr_perm(ndim,batch_perm(n_bp)%n_full_points)

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !    o ndim -- leading dimension of arr/arr_perm
    !    o arr -- array to be send
    !  OUTPUTS
    !    o arr_perm -- array to be received
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i, j, n, mpierr
    real*8, allocatable :: tmp(:,:)
    integer my_point_send_cnt(0:n_tasks-1), my_point_send_off(0:n_tasks-1)
    integer my_point_recv_cnt(0:n_tasks-1), my_point_recv_off(0:n_tasks-1)
    integer my_off(0:n_tasks-1), n_off

    character(*), parameter :: func = 'permute_point_array_real8_2d'

    ! Safety checks
    if(n_bp<=0 .or. n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if(.not. batch_perm(n_bp)%initialized) &
      call aims_stop('value for n_bp not in use', func)

    allocate(tmp(ndim,n_full_points))

    n_off = 0
    my_off(:) = batch_perm(n_bp)%point_send_off(:)
    do i=1,n_my_batches
      n = batch_perm(n_bp)%perm_batch_owner(i)
      do j=1, batches(i)%size
        tmp(:,my_off(n)+j) = arr(:,n_off+j)
      enddo
      n_off = n_off+batches(i)%size
      my_off(n) = my_off(n)+batches(i)%size
    enddo

    my_point_send_cnt(:) = batch_perm(n_bp)%point_send_cnt(:)*ndim
    my_point_send_off(:) = batch_perm(n_bp)%point_send_off(:)*ndim
    my_point_recv_cnt(:) = batch_perm(n_bp)%point_recv_cnt(:)*ndim
    my_point_recv_off(:) = batch_perm(n_bp)%point_recv_off(:)*ndim


    call mpi_alltoallv(tmp,      my_point_send_cnt, my_point_send_off, MPI_REAL8,    &
                       arr_perm, my_point_recv_cnt, my_point_recv_off, MPI_REAL8,    &
                       mpi_comm_global, mpierr)

    deallocate(tmp)

  end subroutine permute_point_array_real8_2d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_real8_1d
  !  NAME
  !    permute_point_array_real8_1d
  !  SYNOPSIS

  subroutine permute_point_array_real8_1d(n_bp, arr, arr_perm)

    !  PURPOSE
    !   see permute_point_array_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    real*8, intent(in)  :: arr(n_full_points)
    real*8, intent(out) :: arr_perm(batch_perm(n_bp)%n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_real8_1d'

    call permute_point_array_real8_2d( &
         n_bp, 1, arr, arr_perm)
  end subroutine permute_point_array_real8_1d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_real8_1d_legacy
  !  NAME
  !    permute_point_array_real8_1d_legacy
  !  SYNOPSIS

  subroutine permute_point_array_real8_1d_legacy(n_bp, ndim, arr, arr_perm)

    !  PURPOSE
    !   see permute_point_array_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    real*8, intent(in)  :: arr(n_full_points)
    real*8, intent(out) :: arr_perm(batch_perm(n_bp)%n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_real8_1d_legacy'

    call permute_point_array_real8_2d( &
         n_bp, ndim, arr, arr_perm)
  end subroutine permute_point_array_real8_1d_legacy
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_real8_3d_legacy
  !  NAME
  !    permute_point_array_real8_3d_legacy
  !  SYNOPSIS

  subroutine permute_point_array_real8_3d_legacy(n_bp, ndim, arr, arr_perm)

    !  PURPOSE
    !   see permute_point_array_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    real*8, intent(in)  :: arr(ndim,1,n_full_points)
    real*8, intent(out) :: arr_perm(ndim,1,batch_perm(n_bp)%n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_real8_3d_legacy'

    call permute_point_array_real8_2d( &
         n_bp, ndim, arr, arr_perm)
  end subroutine permute_point_array_real8_3d_legacy
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_int_2d
  !  NAME
  !    permute_point_array_int_2d
  !  SYNOPSIS

  subroutine permute_point_array_int_2d(n_bp, ndim, arr, arr_perm)

    !  PURPOSE
    !    Permutes an array of integration points from original to permuted distribution
    !  USES
    use dimensions, only: n_full_points, n_my_batches
    use grids, only: batches
    use mpi_tasks
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    integer, intent(in)  :: arr(ndim,n_full_points)
    integer, intent(out) :: arr_perm(ndim,batch_perm(n_bp)%n_full_points)

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !    o ndim -- leading dimension of arr/arr_perm
    !    o arr -- array to be send
    !  OUTPUTS
    !    o arr_perm -- array to be received
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i, j, n, mpierr
    integer, allocatable :: tmp(:,:)
    integer my_point_send_cnt(0:n_tasks-1), my_point_send_off(0:n_tasks-1)
    integer my_point_recv_cnt(0:n_tasks-1), my_point_recv_off(0:n_tasks-1)
    integer my_off(0:n_tasks-1), n_off

    character(*), parameter :: func = 'permute_point_array_int_2d'

    ! Safety checks
    if(n_bp<=0 .or. n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if(.not. batch_perm(n_bp)%initialized) &
      call aims_stop('value for n_bp not in use', func)

    allocate(tmp(ndim,n_full_points))

    n_off = 0
    my_off(:) = batch_perm(n_bp)%point_send_off(:)
    do i=1,n_my_batches
      n = batch_perm(n_bp)%perm_batch_owner(i)
      do j=1, batches(i)%size
        tmp(:,my_off(n)+j) = arr(:,n_off+j)
      enddo
      n_off = n_off+batches(i)%size
      my_off(n) = my_off(n)+batches(i)%size
    enddo

    my_point_send_cnt(:) = batch_perm(n_bp)%point_send_cnt(:)*ndim
    my_point_send_off(:) = batch_perm(n_bp)%point_send_off(:)*ndim
    my_point_recv_cnt(:) = batch_perm(n_bp)%point_recv_cnt(:)*ndim
    my_point_recv_off(:) = batch_perm(n_bp)%point_recv_off(:)*ndim


    call mpi_alltoallv(tmp,      my_point_send_cnt, my_point_send_off, MPI_INTEGER,    &
                       arr_perm, my_point_recv_cnt, my_point_recv_off, MPI_INTEGER,    &
                       mpi_comm_global, mpierr)

    deallocate(tmp)

  end subroutine permute_point_array_int_2d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_int_1d
  !  NAME
  !    permute_point_array_int_1d
  !  SYNOPSIS

  subroutine permute_point_array_int_1d(n_bp, arr, arr_perm)

    !  PURPOSE
    !   see permute_point_array_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    integer, intent(in)  :: arr(n_full_points)
    integer, intent(out) :: arr_perm(batch_perm(n_bp)%n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_int_1d'

    call permute_point_array_int_2d( &
         n_bp, 1, arr, arr_perm)
  end subroutine permute_point_array_int_1d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_back_real8_2d
  !  NAME
  !    permute_point_array_back_real8_2d
  !  SYNOPSIS

  subroutine permute_point_array_back_real8_2d(n_bp, ndim, arr_perm, arr)

    !  PURPOSE
    !    Permutes an array of integration points from permuted to original distribution
    !  USES
    use dimensions, only: n_full_points, n_my_batches
    use grids, only: batches
    use mpi_tasks
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    real*8, intent(in)  :: arr_perm(ndim,batch_perm(n_bp)%n_full_points)
    real*8, intent(out) :: arr(ndim,n_full_points)

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !    o ndim -- leading dimension of arr/arr_perm
    !    o arr_perm -- array to be send
    !  OUTPUTS
    !    o arr -- array to be received
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i, j, n, mpierr
    real*8, allocatable :: tmp(:,:)
    integer my_point_send_cnt(0:n_tasks-1), my_point_send_off(0:n_tasks-1)
    integer my_point_recv_cnt(0:n_tasks-1), my_point_recv_off(0:n_tasks-1)
    integer my_off(0:n_tasks-1), n_off

    character(*), parameter :: func = 'permute_point_array_back_real8_2d'

    ! Safety checks
    if(n_bp<=0 .or. n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if(.not. batch_perm(n_bp)%initialized) &
      call aims_stop('value for n_bp not in use', func)

    allocate(tmp(ndim,n_full_points))

    my_point_send_cnt(:) = batch_perm(n_bp)%point_recv_cnt(:)*ndim
    my_point_send_off(:) = batch_perm(n_bp)%point_recv_off(:)*ndim
    my_point_recv_cnt(:) = batch_perm(n_bp)%point_send_cnt(:)*ndim
    my_point_recv_off(:) = batch_perm(n_bp)%point_send_off(:)*ndim

    call mpi_alltoallv(arr_perm, my_point_send_cnt, my_point_send_off, MPI_REAL8,    &
                       tmp,      my_point_recv_cnt, my_point_recv_off, MPI_REAL8,    &
                       mpi_comm_global, mpierr)

    n_off = 0
    my_off(:) = batch_perm(n_bp)%point_send_off(:)
    do i=1,n_my_batches
      n = batch_perm(n_bp)%perm_batch_owner(i)
      do j=1, batches(i)%size
        arr(:,n_off+j) = tmp(:,my_off(n)+j)
      enddo
      n_off = n_off+batches(i)%size
      my_off(n) = my_off(n)+batches(i)%size
    enddo

    deallocate(tmp)

  end subroutine permute_point_array_back_real8_2d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_back_real8_1d
  !  NAME
  !    permute_point_array_back_real8_1d
  !  SYNOPSIS

  subroutine permute_point_array_back_real8_1d(n_bp, arr_perm, arr)

    !  PURPOSE
    !   see permute_point_array_back_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    real*8, intent(in) :: arr_perm(batch_perm(n_bp)%n_full_points)
    real*8, intent(out)  :: arr(n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_back_real8_1d'

    call permute_point_array_back_real8_2d( &
         n_bp, 1, arr_perm, arr)
  end subroutine permute_point_array_back_real8_1d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_back_real8_1d_legacy
  !  NAME
  !    permute_point_array_back_real8_1d_legacy
  !  SYNOPSIS

  subroutine permute_point_array_back_real8_1d_legacy(n_bp, ndim, arr_perm, arr)

    !  PURPOSE
    !   see permute_point_array_back_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    real*8, intent(in) :: arr_perm(batch_perm(n_bp)%n_full_points)
    real*8, intent(out)  :: arr(n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_back_real8_1d_legacy'

    call permute_point_array_back_real8_2d( &
         n_bp, ndim, arr_perm, arr)
  end subroutine permute_point_array_back_real8_1d_legacy
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_back_real8_3d_legacy
  !  NAME
  !    permute_point_array_back_real8_3d_legacy
  !  SYNOPSIS

  subroutine permute_point_array_back_real8_3d_legacy(n_bp, ndim, arr_perm, arr)

    !  PURPOSE
    !   see permute_point_array_back_real8_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    real*8, intent(in) :: arr_perm(ndim,1,batch_perm(n_bp)%n_full_points)
    real*8, intent(out)  :: arr(ndim,1,n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_back_real8_3d_legacy'

    call permute_point_array_back_real8_2d( &
         n_bp, ndim, arr_perm, arr)
  end subroutine permute_point_array_back_real8_3d_legacy
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_back_int_2d
  !  NAME
  !    permute_point_array_back_int_2d
  !  SYNOPSIS

  subroutine permute_point_array_back_int_2d(n_bp, ndim, arr_perm, arr)

    !  PURPOSE
    !    Permutes an array of integration points from permuted to original distribution
    !  USES
    use dimensions, only: n_full_points, n_my_batches
    use grids, only: batches
    use mpi_tasks
    implicit none
    !  ARGUMENTS

    integer, intent(in) :: n_bp
    integer, intent(in) :: ndim
    integer, intent(in)  :: arr_perm(ndim,batch_perm(n_bp)%n_full_points)
    integer, intent(out) :: arr(ndim,n_full_points)

    !  INPUTS
    !    o n_bp -- number of batch distribution to use
    !    o ndim -- leading dimension of arr/arr_perm
    !    o arr_perm -- array to be send
    !  OUTPUTS
    !    o arr -- array to be received
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer i, j, n, mpierr
    integer, allocatable :: tmp(:,:)
    integer my_point_send_cnt(0:n_tasks-1), my_point_send_off(0:n_tasks-1)
    integer my_point_recv_cnt(0:n_tasks-1), my_point_recv_off(0:n_tasks-1)
    integer my_off(0:n_tasks-1), n_off

    character(*), parameter :: func = 'permute_point_array_back_int_2d'

    ! Safety checks
    if(n_bp<=0 .or. n_bp>max_batch_permutations) &
      call aims_stop('illegal value for n_bp', func)

    if(.not. batch_perm(n_bp)%initialized) &
      call aims_stop('value for n_bp not in use', func)

    allocate(tmp(ndim,n_full_points))

    my_point_send_cnt(:) = batch_perm(n_bp)%point_recv_cnt(:)*ndim
    my_point_send_off(:) = batch_perm(n_bp)%point_recv_off(:)*ndim
    my_point_recv_cnt(:) = batch_perm(n_bp)%point_send_cnt(:)*ndim
    my_point_recv_off(:) = batch_perm(n_bp)%point_send_off(:)*ndim

    call mpi_alltoallv(arr_perm, my_point_send_cnt, my_point_send_off, MPI_INTEGER,    &
                       tmp,      my_point_recv_cnt, my_point_recv_off, MPI_INTEGER,    &
                       mpi_comm_global, mpierr)

    n_off = 0
    my_off(:) = batch_perm(n_bp)%point_send_off(:)
    do i=1,n_my_batches
      n = batch_perm(n_bp)%perm_batch_owner(i)
      do j=1, batches(i)%size
        arr(:,n_off+j) = tmp(:,my_off(n)+j)
      enddo
      n_off = n_off+batches(i)%size
      my_off(n) = my_off(n)+batches(i)%size
    enddo

    deallocate(tmp)

  end subroutine permute_point_array_back_int_2d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/permute_point_array_back_int_1d
  !  NAME
  !    permute_point_array_back_int_1d
  !  SYNOPSIS

  subroutine permute_point_array_back_int_1d(n_bp, arr_perm, arr)
    !  PURPOSE
    !   see permute_point_array_back_int_2d
    !  USES
    use dimensions, only: n_full_points
    implicit none
    !  ARGUMENTS
    integer, intent(in) :: n_bp
    integer, intent(in) :: arr_perm(batch_perm(n_bp)%n_full_points)
    integer, intent(out)  :: arr(n_full_points)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Development version, FHI-aims (2017).
    !  SOURCE
    character(*), parameter :: func = 'permute_point_array_back_int_1d'

   call permute_point_array_back_int_2d( &
         n_bp, 1, arr_perm, arr)
  end subroutine permute_point_array_back_int_1d
  !******
  !------------------------------------------------------------------------------
  !****s* load_balancing/set_full_local_ovlp
  !  NAME
  !    set_full_local_ovlp
  !  SYNOPSIS
  subroutine set_full_local_ovlp(matrix, overlap_matrix, overlap_matrix_complex, k_point_global)
  !  PURPOSE
  !    Create the overlap matrix (real or complex) that will passed to the
  !    eigensolver from the real-space overlap matrix stored in the full local
  !    matrix format
  !  USES
    use mpi_tasks, only: aims_stop
    use runtime_choices, only: use_scalapack, use_wf_extrapolation, &
        real_eigenvectors
    use scalapack_wrapper, only: get_set_full_local_matrix_scalapack, &
        n_local_matrix_size, mxld, mxcol, ovlp, ovlp_complex
    use full_local_mat_lapack, only: get_set_full_local_matrix_lapack
    implicit none
  !  ARGUMENTS
    real*8, intent(in) :: matrix(n_local_matrix_size)
    real*8, intent(out), optional :: overlap_matrix(:)
    complex*16, intent(out), optional :: overlap_matrix_complex(:)
    integer, intent(in), optional :: k_point_global
  !  INPUTS
  !    o matrix -- real-space overlap matrix in full local format
  !  OUTPUT
  !    o hamiltonian -- real overlap matrix in eigensolver format (only for LAPACK)
  !    o hamiltonian_complex -- complex overlap matrix in eigensolver format (only for LAPACK)
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
    character(*), parameter :: func = 'set_full_local_ovlp'

    real*8 :: dummy_hamiltonian(1)
    complex*16 :: dummy_hamiltonian_complex(1)

    if (use_scalapack) then
       ! ScaLAPACK branch updates ovlp directly in scalapack_wrapper, no need to
       ! pass in
       call get_set_full_local_matrix_scalapack(matrix, 0, 1)

       if (use_wf_extrapolation) then
          if (real_eigenvectors) then
             call wf_save_overlap(mxld, mxcol, ovlp)
          else
             call wf_save_overlap_cmplx(mxld, mxcol, ovlp_complex)
          end if
       end if
    else
       if (.not.present(overlap_matrix) .or. &
           .not.present(overlap_matrix_complex)) then
          call aims_stop("Missing mandatory arguments for load balancing &
                         &when using LAPACK, exiting.", func)
       end if
       call get_set_full_local_matrix_lapack &
            (matrix, dummy_hamiltonian, overlap_matrix, &
             dummy_hamiltonian_complex, overlap_matrix_complex, 0, 1, &
             k_point_global)
    end if

  end subroutine set_full_local_ovlp
  !-----------------------------------------------------------------------------------
  !****s* load_balancing/set_full_local_ham
  !  NAME
  !    set_full_local_ham
  !  SYNOPSIS
  subroutine set_full_local_ham(matrix, hamiltonian, hamiltonian_complex, k_point_global)
  !  PURPOSE
  !    Create the Hamiltonian (real or complex) that will passed to the eigensolver
  !    from the real-space Hamiltonian stored in the full local matrix format
  !  USES
    use dimensions, only: n_spin
    use mpi_tasks, only: aims_stop
    use runtime_choices, only: use_scalapack
    use scalapack_wrapper, only: get_set_full_local_matrix_scalapack, &
        n_local_matrix_size
    use full_local_mat_lapack, only: get_set_full_local_matrix_lapack
    implicit none
  !  ARGUMENTS
    real*8, intent(in) :: matrix(n_local_matrix_size, n_spin)
    real*8, intent(out), optional :: hamiltonian(:,:)
    complex*16, intent(out), optional :: hamiltonian_complex(:,:)
    integer, intent(in), optional :: k_point_global
  !  INPUTS
  !    o matrix -- real-space Hamiltonian in full local format
  !    o k_point_global -- global index for k-point (i.e. for k_phase)
  !  OUTPUT
  !    o hamiltonian -- real Hamiltonian in eigensolver format (only for LAPACK)
  !    o hamiltonian_complex -- complex Hamiltonian in eigensolver format (only for LAPACK)
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

    real*8 :: dummy_overlap(1)
    complex*16 :: dummy_overlap_complex(1)
    integer :: i_spin

    character(*), parameter :: func = 'set_full_local_ham'

    if (use_scalapack) then
       ! ScaLAPACK branch updates ham directly in scalapack_wrapper, no need to
       ! pass in
       do i_spin = 1, n_spin
          call get_set_full_local_matrix_scalapack(matrix(1,i_spin), 1, i_spin)
       end do
    else
       if (.not.present(hamiltonian) .or. &
           .not.present(hamiltonian_complex) .or. &
           .not.present(k_point_global)) then
          call aims_stop("Missing mandatory arguments for load balancing &
                         &when using LAPACK, exiting.", func)
       end if
       do i_spin = 1, n_spin
          call get_set_full_local_matrix_lapack &
               (matrix, hamiltonian, dummy_overlap, &
                hamiltonian_complex, dummy_overlap_complex, 1, i_spin, &
                k_point_global)
       end do
    end if

    ! VB - comment on the matrix dimension for the case of more than one spin, which
    !      is very tricky. (with load_balancing.f90)
    !
    ! This routine is used by the 'load_balancing' infrastructure.
    ! It would appear that n_local_matrix_size is the correct variable to use for the dimension.
    ! However, notice that this 'matrix' is filled with values by subroutine
    !
    ! call update_full_matrix_p0X(  &
    !               n_compute_c, n_compute_c, i_basis(1), hamiltonian_shell, &
    !               hamiltonian((i_spin-1)*ld_hamiltonian+1) )
    !
    ! in integrate_hamiltonian_matrix_p2.f90
    ! i.e., the formal handling of the i_spin index (where does spin number 2 start in the array) is
    ! different.(!)
    ! where ld_hamiltonian = batch_perm(n_bp)%n_local_matrix_size
    !
    ! Here, instead, we use a dimension called n_local_matrix_size directly.
    !
    ! Checking in load_balancing.f90 seems to indicate that both definitions are
    ! consistent as of this writing, Oct 18 2014. However, the use of two
    ! different variable names makes me uncomfortable, so I record the definitions here. - VB
  end subroutine set_full_local_ham
  !-----------------------------------------------------------------------------------
  !****s* load_balancing/get_full_local_matrix
  !  NAME
  !    get_full_local_matrix
  !  SYNOPSIS
  subroutine get_full_local_matrix(matrix_local, i_spin, matrix_eigen)
  !  PURPOSE
  !    Converts a matrix stored in the eigensolver's format (BLACS for
  !    ScaLAPACK, upper triangular for LAPACK) into the full local matrix format
  !  USES
    use mpi_tasks, only: aims_stop
    use runtime_choices, only: use_scalapack
    use scalapack_wrapper, only: get_set_full_local_matrix_scalapack, &
        n_local_matrix_size
    use full_local_mat_lapack, only: get_set_full_local_matrix_lapack
    implicit none
  !  ARGUMENTS
    real*8, intent(out) :: matrix_local(n_local_matrix_size)
    integer, intent(in) :: i_spin
    real*8, intent(in), optional :: matrix_eigen(:)
  !  INPUTS
  !    o matrix_eigen -- real matrix in eigensolver's format (only for LAPACK)
  !    o matrix_eigen_complex -- complex matrix in eigensolver's format (only for LAPACK)
  !    o i_spin -- spin channel to use in matrix_eigen/matrix_eigen_complex
  !  OUTPUT
  !    o matrix_local -- matrix in full local format
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

    real*8 :: dummy_overlap(1)
    complex*16 :: dummy_matrix_eigen_complex(1)
    complex*16 :: dummy_overlap_complex(1)

    character(*), parameter :: func = 'get_full_local_matrix'

    if (use_scalapack) then
       ! For the BLACS format, the matrix to be converted is assumed to be
       ! pre-emptively stored into the ham/ham_complex module variable
       call get_set_full_local_matrix_scalapack(matrix_local, 2, i_spin)
    else
       if (.not.present(matrix_eigen)) then
          call aims_stop("Missing mandatory arguments for load balancing &
                         &when using LAPACK, exiting.", func)
       end if
       call get_set_full_local_matrix_lapack &
            (matrix_local, matrix_eigen, dummy_overlap, &
             dummy_matrix_eigen_complex, dummy_overlap_complex, 2, i_spin, 1)
    end if

  end subroutine get_full_local_matrix
  !-----------------------------------------------------------------------------------
  !****s* load_balancing/init_comm_full_local_matrix
  !  NAME
  !    init_comm_full_local_matrix
  !  SYNOPSIS
  subroutine init_comm_full_local_matrix(n_basis_local_full, i_basis_local_full)
  !  PURPOSE
  !    Initializes the communication for the conversion of full local matrices
  !    into formats suitable for eigensolvers
  !  USES
    use runtime_choices, only: use_scalapack
    use scalapack_wrapper, only: init_comm_full_local_matrix_scalapack
    use full_local_mat_lapack, only: init_comm_full_local_matrix_lapack
    implicit none
  !  ARGUMENTS

    integer, intent(in) :: n_basis_local_full
    integer, intent(in) :: i_basis_local_full(n_basis_local_full)
  !  INPUTS
  !    o n_basis_local_full -- number of local basis functions for current task
  !    o i_basis_local_full -- array with local basis functions for current task
  !  OUTPUT
  !    all index arrays for local_index communication are set
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
    character(*), parameter :: func = 'init_comm_full_local_matrix'

    if (use_scalapack) then
       call init_comm_full_local_matrix_scalapack(n_basis_local_full, i_basis_local_full)
    else
       call init_comm_full_local_matrix_lapack(n_basis_local_full, i_basis_local_full)
    end if
  end subroutine init_comm_full_local_matrix
  !******

  subroutine print_sparse_to_dense_local_index_cpscf(n_bp, first_order_H, which)

      !  PURPOSE
      !  print local dense matrix into global dense matrix format

      !  USES
      use dimensions, only: n_basis, n_centers_basis_I 
      use synchronize_mpi_basic, only: sync_vector
      use mpi_tasks, only: aims_stop, myid

      implicit none
      ! ARGUMENTS
      integer :: n_bp
      real*8 :: first_order_H(*)
      integer, optional :: which

      ! local variables
      integer n, i, j, i_global, j_global
      real*8, allocatable :: tmp_matrix(:,:)
      integer cnt
      integer :: id_H=1,id_DM=2

      ! in H2, local index, this must equal(just for H2 debuging)
      if (n_basis .ne. n_centers_basis_I) call aims_stop("cpscf local_index error")
      !allocate(tmp_matrix(n_basis,n_basis))
      cnt = 20
      if (cnt > n_basis) cnt = n_basis
      allocate(tmp_matrix(cnt,1))
      tmp_matrix = 0.0d0

      n = batch_perm(n_bp)%n_basis_local
      do i=1,n
        do j=1,i
            i_global = batch_perm(n_bp)%i_basis_local(i)
            j_global = batch_perm(n_bp)%i_basis_local(j)
            
            if (i_global >= 1 .and. i_global <= cnt .and. j_global == 1) then
                tmp_matrix(i_global, j_global) = first_order_H((i*(i-1))/2 + j)
                if (which == id_DM) print *, myid, 'local_DM(', i_global, j_global, ')', first_order_H((i*(i-1))/2 + j)
                !if (i_global .ne. j_global) then
                !    !tmp_matrix(j_global, i_global) = first_order_H((i*(i-1))/2 + j)
                !    print *, myid, 'H/DM(', j_global, i_global, ')', first_order_H((i*(i-1))/2 + j)
                !endif
            endif
        enddo
      enddo

      if (which == id_H) then
          !call sync_vector(tmp_matrix, n_basis * n_basis)
          call sync_vector(tmp_matrix, cnt * 1)
      endif
      if (myid == 0) call print_matrix_real8_2d(tmp_matrix, cnt, 1)
      if (allocated(tmp_matrix)) deallocate(tmp_matrix)

  end subroutine print_sparse_to_dense_local_index_cpscf
  
  subroutine print_matrix_real8_2d(matrix, len_d1, len_d2)
      ! @wuyangjun
      ! PURPOSE: print `real*8 matrix(len_d1,len_d2)` 

      ! ARGUMENTS
      real*8 :: matrix(len_d1, len_d2)
      integer :: len_d1, len_d2
      
      ! LOCAL VARIABLES
      integer :: i, j

      do i = 1, len_d1
        do j = 1, len_d2
        print *, 'ALL_REDUCE_H/DM:(', i, j, ')', matrix(i,j)
        enddo
      enddo

  end subroutine print_matrix_real8_2d
end module load_balancing
!******
