!****h* FHI-aims/sparse_tensor
!  NAME
!    sparse_tensor
!  SYNOPSIS

module sparse_tensor

  !  PURPOSE
  !
  !    Provide facilities to store and manipulate sparse and distributed
  !    tensors of general rank.  Best to think of the rank to be 3.
  !
  !    The idea is to keep rectangular blocks of values distributed over all
  !    nodes.  The distribution scheme can be anything; each node only knows
  !    about the blocks it owns.
  !
  !  USES

  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !
  !  CONVENTIONS
  !
  !    ten -- The sparse tensor derived type instance
  !    rank -- Number of indices needed to choose one tensor entry
  !    nzb -- non-zero block
  !           specified by top(1:rank) and shp(1:rank) in the sense that the
  !           block spans  top <= ind < top + shp
  !    nzv -- non-zero value
  !           n_nzv == sum(product(shp(:,:), 1)) == size(val)
  !    top -- Indices of top (top, left, front, ...) element.
  !           E.g. (/1, 1, 1/)
  !    shp -- Sizes in each dimension.
  !    bot -- Indices /below/ bottom (beyond right, behind back, ...) element.
  !           bot := top + shp
  !           E.g. (/1+shp(1), 1+shp(2), 1+shp(3)/)
  !
  !  PUBLIC SUBROUTINES
  !
  !     * Constructor/destructor, ...
  !     alloc_sp_ten -- Prepare sparse tensor data structure
  !     dealloc_sp_ten -- Free allocated memory and reset data structure
  !     copy_sp_ten -- Make a clone
  !     join_blocks -- Simplify local blocking
  !     simplified_tensor -- Simplify local blocking of a full tensor
  !     upgrade_blocking -- Pure blocking -> descriptor/tensor conversion
  !
  !     * Communication
  !     redistribute_sp_ten -- Collectively redistribute ten according to wants
  !
  !     * Wants/descriptors
  !     get_scalapack_blocking -- Blocking for block-cyclic distribution
  !     get_scalapack_descriptor -- Get descriptor for 2D scalapack array
  !
  !     * Extractors
  !     extract_matricization -- tensor -> 2D-matrix
  !     update_matricization -- 2D-matrix -> tensor
  !
  !     * Checker and output
  !     check_global_aliasing -- Check for aliasing in global matrix
  !     check_local_aliasing -- Check for aliasing in local storage
  !     debug_output_blocking -- Output blocking
  !     debug_local_sparse2python_dok -- Output local parts as Python dict
  !     debug_global_sparse2python_dok -- Output global matrix as Python dict
  !
  !  SOURCE

  ! What kind of data does ten contain?
  integer, parameter :: SPTEN_INVALID = 0
  ! SPTEN_INVALID: Not valid at all
  integer, parameter :: SPTEN_BLOCKING = 1
  ! SPTEN_BLOCKING: Contains blocking (i.e. local part of nonzero pattern) only
  integer, parameter :: SPTEN_DESCRIPTOR = 2
  ! SPTEN_DESCRIPTOR: Additionally contains offsets (and maybe non-default
  !                   strides) to allow access into an external array
  integer, parameter :: SPTEN_FULL = 3
  ! SPTEN_FULL: Additionally contains the actual values in ten%val.
  integer, parameter :: SPTEN_N_MODE = 3

  logical, parameter :: SPTEN_HAS_OFF(SPTEN_N_MODE) = &
  & (/.false., .true., .true./)
  logical, parameter :: SPTEN_HAS_VAL(SPTEN_N_MODE) = &
  & (/.false., .false., .true./)

  type sp_ten
     character*20 :: name
     ! General (global) properties
     integer :: rank     ! Rank (typically 3)
     integer, allocatable :: glb_shp(:) ! (rank): Shape of global tensor
     integer :: mode

     ! Local properties
     integer :: n_nzb    ! Number of locally stored nonzero blocks
     integer :: n_nzv    ! Size of internal %val(:) or external var(:) array
     ! Position of blocks in global array:
     !   top(:, i_nzb) : (top(:, i_nzb) + shp(:, i_nzb) - 1)
     integer, allocatable :: top(:,:) ! (rank, n_nzb)
     integer, allocatable :: shp(:,:) ! (rank, n_nzb)
     ! The strides array may or may not be allocated.
     ! If it is not, Fortran ordering is assumed.
     ! High-level use should always be done through get_strides().
     integer, allocatable :: str(:,:) ! (rank, n_nzb)
     ! Offset of blocks in local array (only if SPTEN_HAS_OFF(mode)):
     integer, allocatable :: off(:)   ! (n_nzb)
     ! Local storage (only if SPTEN_HAS_VAL(mode)):
     real*8, allocatable :: val(:)    ! (n_nzv)
  end type sp_ten

contains

  !----------------------------------------------------------------------------
  !****s* sparse_tensor/alloc_sp_ten
  !  NAME
  !    alloc_sp_ten
  !  SYNOPSIS

  subroutine alloc_sp_ten(ten, rank, glb_shp, mode, n_nzb, n_nzv, &
  &                       name, need_str)

    !  PURPOSE
    !    Allocate tensor ten
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: ten
    integer, intent(IN) :: rank
    integer, intent(IN) :: glb_shp(rank)
    integer, intent(IN) :: mode
    integer, intent(IN) :: n_nzb
    integer, intent(IN) :: n_nzv
    character(*), intent(IN) :: name
    logical, intent(IN), optional :: need_str

    !  INPUTS
    !    o rank -- Rank of the tensor (usually 3)
    !    o n_nzb -- Number of nonzero blocks
    !    o n_nzv -- Number of nonzero values
    !    o need_str -- If present and true, also allocate the %str field
    !  OUTPUTS
    !    o ten -- [mode] Allocated sparse tensor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: info
    character(*), parameter :: func = 'alloc_sp_ten'

    call dealloc_sp_ten(ten)

    if (mode < 1 .or. mode > SPTEN_N_MODE) then
       call aims_stop('invalid mode', func)
    end if

    ten%rank = rank
    ten%n_nzb = n_nzb
    ten%mode = mode
    allocate(ten%glb_shp(rank), stat=info)
    call check_allocation(info, 'alloc_sp_ten: glb_shp')
    ten%glb_shp = glb_shp
    allocate(ten%top(rank, n_nzb), stat=info)
    call check_allocation(info, 'alloc_sp_ten: top')
    ten%top = 0   ! Invalid
    allocate(ten%shp(rank, n_nzb), stat=info)
    call check_allocation(info, 'alloc_sp_ten: shp')
    ten%shp = -1  ! Invalid

    if (SPTEN_HAS_OFF(mode)) then
       ten%n_nzv = n_nzv
       allocate(ten%off(n_nzb), stat=info)
       call check_allocation(info, 'alloc_sp_ten: off')
       ten%off = -1
    else
       if (n_nzv /= 0) call aims_stop('n_nzv /= 0 for blocking', func)
       ten%n_nzv = 0
    end if
    if (SPTEN_HAS_VAL(mode)) then
       allocate(ten%val(n_nzv), stat=info)
       call check_allocation(info, 'alloc_sp_ten: val')
       ten%val = 0.d0
    end if
    ten%name = name
    if (present(need_str)) then
       if (need_str) then
          if (.not. SPTEN_HAS_OFF(mode)) call aims_stop('Strides w/o off',func)
          allocate(ten%str(rank, n_nzb), stat=info)
          call check_allocation(info, 'alloc_sp_ten: top')
       end if
    end if

  end subroutine alloc_sp_ten
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/dealloc_sp_ten
  !  NAME
  !    dealloc_sp_ten
  !  SYNOPSIS

  subroutine dealloc_sp_ten(ten)

    !  PURPOSE
    !    Deallocate tensor ten
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(INOUT) :: ten

    !  INPUTS
    !    o ten -- [any] Sparse tensor []
    !  OUTPUTS
    !    o ten -- [invalid] Deallocated sparse tensor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    ten%rank = 0
    ten%n_nzb = 0
    ten%n_nzv = 0
    ten%mode = SPTEN_INVALID
    if (allocated(ten%glb_shp)) deallocate(ten%glb_shp)
    if (allocated(ten%top)) deallocate(ten%top)
    if (allocated(ten%shp)) deallocate(ten%shp)
    if (allocated(ten%str)) deallocate(ten%str)
    if (allocated(ten%off)) deallocate(ten%off)
    if (allocated(ten%val)) deallocate(ten%val)

  end subroutine dealloc_sp_ten
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/copy_sp_ten
  !  NAME
  !    copy_sp_ten
  !  SYNOPSIS

  subroutine copy_sp_ten(ten_in, ten_out, name)

    !  PURPOSE
    !    Copy sparse tensor ten_in to ten_out
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_in
    type(sp_ten), intent(OUT) :: ten_out
    character(*), intent(IN) :: name

    !  INPUTS
    !    o ten_in -- [any] Sparse tensor
    !  OUTPUTS
    !    o ten_out -- [ten_in%mode] Copy of ten_in
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: mode

    mode = ten_in%mode

    call alloc_sp_ten(ten_out, ten_in%rank, ten_in%glb_shp, ten_in%mode, &
    &                 ten_in%n_nzb, ten_in%n_nzv, name)

    ten_out%top = ten_in%top
    ten_out%shp = ten_in%shp
    if (SPTEN_HAS_OFF(mode)) ten_out%off = ten_in%off
    if (allocated(ten_in%str)) then
       allocate(ten_out%str(ten_in%rank, ten_in%n_nzb))
       ten_out%str = ten_in%str
    end if
    if (SPTEN_HAS_VAL(mode)) ten_out%val = ten_in%val

  end subroutine copy_sp_ten
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/join_blocks
  !  NAME
  !    join_blocks
  !  SYNOPSIS

  subroutine join_blocks(ten_in, ten_out, name)

    !  PURPOSE
    !
    !    Regroup and join the blocks of ten_in.  Therefore, ten_out contains
    !    the same nonzero pattern as ten_in, but less blocks.  Only outputs
    !    the blocking.
    !
    !    While this routine could in principle be made quite intellegent, it
    !    is extremely dumb at the moment and only will join undubious cases.
    !    Still, this might be enough for our purposes.
    !
    !  USES

    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_in
    type(sp_ten), intent(OUT) :: ten_out
    character(*), intent(IN) :: name

    !  INPUTS
    !    o ten_in -- [>=blocking] Input blocking
    !  OUTPUTS
    !   o ten_out -- [blocking] Packed blocking.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank, i, glb_shp(ten_in%rank)
    integer :: i_dim_sort
    integer :: i_dim, j_dim, i_nzb, j_nzb, n_nzb_uptonow, n_nzb
    integer, allocatable :: new2old(:)
    integer, allocatable :: topshp(:,:)
    integer, allocatable :: alltop(:,:), allshp(:,:)
    integer, allocatable :: sorttop(:,:), sortshp(:,:)
    integer :: info
    logical, allocatable :: done(:)
    integer :: this_top(ten_in%rank), this_bot(ten_in%rank)
    integer :: this_shp(ten_in%rank)
    character(*), parameter :: func = 'join_blocks'

    rank = ten_in%rank
    glb_shp = ten_in%glb_shp
    allocate(alltop(rank, ten_in%n_nzb), stat=info)
    call check_allocation(info, 'alltop')
    allocate(allshp(rank, ten_in%n_nzb), stat=info)
    call check_allocation(info, 'allshp')
    allocate(sorttop(rank, ten_in%n_nzb), stat=info)
    call check_allocation(info, 'sorttop')
    allocate(sortshp(rank, ten_in%n_nzb), stat=info)
    call check_allocation(info, 'sortshp')
    allocate(topshp(2*rank, ten_in%n_nzb), stat=info)
    call check_allocation(info, 'topshp')
    allocate(new2old(ten_in%n_nzb), stat=info)
    call check_allocation(info, 'new2old')
    allocate(done(ten_in%n_nzb), stat=info)
    call check_allocation(info, 'done')

    alltop(:,:) = ten_in%top
    allshp(:,:) = ten_in%shp
    n_nzb = ten_in%n_nzb

    ! One pass per index
    do i_dim = 1, rank
       ! Glue together subsequent quarders in dimension i_dim which agree in
       ! all other dimensions.
       !
       ! Gather all entries where the other dimensions agree by sorting:
       do i = 1, rank
          topshp(2*i-1, 1:n_nzb) = alltop(i, 1:n_nzb)
          topshp(2*i,   1:n_nzb) = allshp(i, 1:n_nzb)
       end do
       i_dim_sort = modulo(i_dim, rank) + 1
       call heapsort_index(topshp, 2*rank, n_nzb, &
       &                   2*i_dim_sort-1, new2old, .true.)
       do i_nzb = 1, n_nzb
          sorttop(:, i_nzb) = alltop(:, new2old(i_nzb))
          sortshp(:, i_nzb) = allshp(:, new2old(i_nzb))
       end do

       ! Now do the actual work
       done = .false.
       n_nzb_uptonow = 0
       do i_nzb = 1, n_nzb
          if (done(i_nzb)) cycle

          this_top = sorttop(:, i_nzb)
          this_shp = sortshp(:, i_nzb)
          this_bot = this_top + this_shp

          SEARCH_J: do j_nzb = i_nzb+1, n_nzb
             if (done(j_nzb)) cycle
             ! Check if they still match:
             do i = 1, rank-1
                j_dim = modulo(i_dim + i - 1, rank) + 1
                if (this_top(j_dim) /= sorttop(j_dim, j_nzb)) exit SEARCH_J
                if (this_shp(j_dim) /= sortshp(j_dim, j_nzb)) exit SEARCH_J
             end do
             if (this_bot(i_dim) == sorttop(i_dim, j_nzb)) then
                done(j_nzb) = .true.
                this_shp(i_dim) = this_shp(i_dim) + sortshp(i_dim, j_nzb)
                this_bot(i_dim) = this_bot(i_dim) + sortshp(i_dim, j_nzb)
             end if
          end do SEARCH_J
          ! Done; save.
          done(i_nzb) = .true.
          n_nzb_uptonow = n_nzb_uptonow + 1
          alltop(:, n_nzb_uptonow) = this_top
          allshp(:, n_nzb_uptonow) = this_shp
       end do
       n_nzb = n_nzb_uptonow
    end do

    call alloc_sp_ten(ten_out, rank, glb_shp, SPTEN_BLOCKING, n_nzb, 0, name)
    ten_out%top = alltop(:, 1:n_nzb)
    ten_out%shp = allshp(:, 1:n_nzb)
    deallocate(alltop, allshp, sorttop, sortshp, topshp, new2old, done)

  end subroutine join_blocks
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/simplified_tensor
  !  NAME
  !    simplified_tensor
  !  SYNOPSIS

  subroutine simplified_tensor(ten_in, val_in, ten_out, name)

    !  PURPOSE
    !
    !    Regroup and join the blocks of ten_in.  Therefore, ten_out contains
    !    the same data as ten_in, but less blocks, and may therefore be more
    !    efficient.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_in
    real*8, intent(IN) :: val_in(ten_in%n_nzv)
    type(sp_ten), intent(OUT) :: ten_out
    character(*), intent(IN) :: name

    !  INPUTS
    !    o ten_in -- [>=descriptor] Input Tensor
    !    o val_in -- Corresponding data (if in doubt, use ten_in%val)
    !  OUTPUTS
    !    o ten_out -- Packed tensor.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'simplified_tensor'

    call join_blocks(ten_in, ten_out, name)
    call extract_tensor(ten_in, val_in, ten_out, name)

  end subroutine simplified_tensor
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/upgrade_blocking
  !  NAME
  !    upgrade_blocking
  !  SYNOPSIS

  subroutine upgrade_blocking(ten, mode, name)

    !  PURPOSE
    !
    !    Upgrade a SPTEN_BLOCKING tensor to either SPTEN_DESCRIPTOR or
    !    SPTEN_FULL by allocating/constructing offsets and values.
    !
    !    For (mode == SPTEN_FULL), the values are initialized to zero.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(INOUT) :: ten
    integer, intent(IN) :: mode
    character(*), intent(IN) :: name

    !  INPUTS
    !    o ten -- [blocking] Tensor containing well-defined blocking info.
    !    o mode -- Either SPTEN_DESCRIPTOR or SPTEN_FULL
    !    o name -- new name for ten
    !  OUTPUTS
    !    o ten -- [mode] Upgraded tensor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: info
    character(*), parameter :: func = 'upgrade_blocking'

    ten%name = name

    if (ten%mode == SPTEN_BLOCKING) then
       if (mode /= SPTEN_DESCRIPTOR .and. mode /= SPTEN_FULL) then
          call aims_stop('Invalid mode', func)
       end if
       ten%mode = SPTEN_DESCRIPTOR
       allocate(ten%off(ten%n_nzb), stat=info)
       call check_allocation(info, 'off', func)
       call construct_off(ten%rank, ten%n_nzb, ten%shp, ten%off)
       call count_nzv(ten%rank, ten%n_nzb, ten%shp, ten%n_nzv)
    end if

    if (ten%mode == SPTEN_DESCRIPTOR .and. mode == SPTEN_FULL) then
       ten%mode = SPTEN_FULL
       allocate(ten%val(ten%n_nzv), stat=info)
       call check_allocation(info, 'val', func)
       ten%val = 0.d0
    end if

  end subroutine upgrade_blocking
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/redistribute_sp_ten
  !  NAME
  !    redistribute_sp_ten
  !  SYNOPSIS

  subroutine redistribute_sp_ten(ten, loc_want, name, val_in)

    !  PURPOSE
    !
    !    Redistribute ten.
    !
    !    Each processor specifies what he wants by a sequence of cuboids
    !    (blocks).  In general, the union of all blocks on all nodes will be
    !    the whole (dense) tensor.  Still, the resulting redistributed ten
    !    will have the same sparsity pattern as the input, though some blocks
    !    may need to be split.
    !
    !    What a particular processor owns is given in ten.
    !
    !  USES
    use mpi_tasks, only: MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUCCESS, &
        n_tasks, mpi_comm_global, use_mpi, aims_stop, check_allocation
    use localorb_io, only: localorb_info, OL_norm, use_unit
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(INOUT) :: ten
    type(sp_ten), intent(IN) :: loc_want
    character(*), intent(IN) :: name
    real*8, intent(IN), optional :: val_in(ten%n_nzv)

    !  INPUTS
    !    o ten -- [full/descriptor] Sparse tensor in old distribution
    !    o loc_want -- [blocking] locally wanted blocks
    !  OUTPUTS
    !    o ten -- [full] Redistributed sparse tensor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank, glb_shp(ten%rank)

    integer :: i_task, i_nzb
    integer :: n_nzb_new, n_nzv_new
    integer :: left, right
    character*50 :: glb_want_name

    ! Global wants
    type(sp_ten) :: glb_want
    integer, allocatable :: task2n_nzb_want(:), want2task(:)

    ! Intersection (its) of local holds and global wants
    type(sp_ten) :: its
    integer :: n_nzb_its, n_nzb_its_max
    integer, allocatable :: ptr_its(:,:)
    ! Send buffers
    integer :: n_size
    integer, allocatable :: recver2n_nzbv(:,:), sender2n_nzbv(:,:)
    integer, allocatable :: scnt_nzbi(:), scnt_nzv(:)
    integer, allocatable :: rcnt_nzbi(:), rcnt_nzv(:)
    integer, allocatable :: sdsp_nzbi(:), sdsp_nzv(:)
    integer, allocatable :: rdsp_nzbi(:), rdsp_nzv(:)

    ! General (error) thingys
    integer :: info, mpierr
    integer :: dummy_i2(2,1)
    character*150 :: info_str
    character(*), parameter :: func = 'redistribute_sp_ten'

    write(info_str, &
    &     "(2X,'| Redistributing ',A,' to ',A,' according to ',A,'.')") &
    & trim(ten%name), trim(name), trim(loc_want%name)
    call localorb_info(info_str, use_unit, '(A)', OL_norm)

    call check_global_aliasing(ten, 'ten', func)
    call check_local_aliasing(ten, 'ten', func)
    call check_global_aliasing(loc_want, 'loc_want', func)

    ! --- Initialize

    rank = ten%rank
    glb_shp = ten%glb_shp

    ! --- Synchronize wants

    allocate(task2n_nzb_want(0:n_tasks-1), stat=info)
    call check_allocation(info, 'task2n_nzb_want')

    write(glb_want_name, "('global ',A)") trim(loc_want%name)
    call synchronize_wants(loc_want, glb_want, task2n_nzb_want, glb_want_name)

    ! Map nzb of glb_want to task
    allocate(want2task(glb_want%n_nzb), stat=info)
    call check_allocation(info, 'want2task')
    right = 0
    do i_task = 0, n_tasks-1
       left = right+1
       right = right + task2n_nzb_want(i_task)
       want2task(left:right) = i_task
    end do
    if (right /= glb_want%n_nzb) call aims_stop('glob n_nzb mismatch', func)
    deallocate(task2n_nzb_want)


    ! --- Prepare send buffers

    ! - Find intersections
    ! Count them
    call intersect_sp_ten(glb_want, ten, -1, n_nzb_its, its, dummy_i2)
    n_nzb_its_max = n_nzb_its
    ! Get them
    allocate(ptr_its(2, n_nzb_its_max), stat=info)
    call check_allocation(info, 'ptr_its', func)
    call intersect_sp_ten(glb_want, ten, n_nzb_its_max, n_nzb_its, its,ptr_its)
    if (n_nzb_its /= n_nzb_its_max) call aims_stop('Too many n_nzb', func)

    ! - Serialize intersection contents
    if (present(val_in)) then
       call extract_from_sp_ten(ten, val_in, its, ptr_its(2,:), 'its-val')
    else
       if (.not. allocated(ten%val)) call aims_stop('need val', func)
       call extract_from_sp_ten(ten, ten%val, its, ptr_its(2,:), 'its-val')
    end if
    call dealloc_sp_ten(ten)    ! Done with original data.

    if (.not. use_mpi) then
       call copy_sp_ten(its, ten, name)
       deallocate(want2task)
       return
    end if

    ! - Get send displacements
    allocate(recver2n_nzbv(2, 0:n_tasks-1), stat=info)
    call check_allocation(info, 'recver2n_nzbv')
    recver2n_nzbv = 0
    do i_nzb = 1, its%n_nzb
       i_task = want2task(ptr_its(1, i_nzb))
       ! Update communication sizes
       n_size = product(its%shp(:, i_nzb))
       recver2n_nzbv(1, i_task) = recver2n_nzbv(1, i_task) + 1
       recver2n_nzbv(2, i_task) = recver2n_nzbv(2, i_task) + n_size
    end do
    deallocate(ptr_its)


    ! --- sync n_nzb, n_nzv

    allocate(sender2n_nzbv(2, 0:n_tasks-1), stat=info)
    call check_allocation(info, 'sender2n_nzbv')

    call MPI_Alltoall(recver2n_nzbv, 2, MPI_INTEGER, &
    &                 sender2n_nzbv, 2, MPI_INTEGER, mpi_comm_global, mpierr)
    if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_Alltoall failed')

    n_nzb_new = sum(sender2n_nzbv(1,:))
    n_nzv_new = sum(sender2n_nzbv(2,:))
    allocate(rcnt_nzbi(0:n_tasks-1), rcnt_nzv(0:n_tasks-1), stat=info)
    call check_allocation(info, 'rcnts', func)
    allocate(scnt_nzbi(0:n_tasks-1), scnt_nzv(0:n_tasks-1), stat=info)
    call check_allocation(info, 'scnts', func)
    allocate(rdsp_nzbi(0:n_tasks-1), rdsp_nzv(0:n_tasks-1), stat=info)
    call check_allocation(info, 'rdsps', func)
    allocate(sdsp_nzbi(0:n_tasks-1), sdsp_nzv(0:n_tasks-1), stat=info)
    call check_allocation(info, 'sdsps', func)

    ! sender2n_nzbv is about what I receive from the distinct senders
    rcnt_nzbi = rank * sender2n_nzbv(1,:)
    rcnt_nzv = sender2n_nzbv(2,:)
    ! recver2n_nzbv is about what I send to the distinct receivers
    scnt_nzbi = rank * recver2n_nzbv(1,:)
    scnt_nzv = recver2n_nzbv(2,:)

    rdsp_nzbi(0) = 0
    rdsp_nzv(0) = 0
    sdsp_nzbi(0) = 0
    sdsp_nzv(0) = 0
    do i_task = 1, n_tasks-1
       rdsp_nzbi(i_task) = rdsp_nzbi(i_task-1) + rcnt_nzbi(i_task-1)
       rdsp_nzv(i_task)  = rdsp_nzv(i_task-1)  + rcnt_nzv(i_task-1)
       sdsp_nzbi(i_task) = sdsp_nzbi(i_task-1) + scnt_nzbi(i_task-1)
       sdsp_nzv(i_task)  = sdsp_nzv(i_task-1)  + scnt_nzv(i_task-1)
    end do

!!$    write(0,*) myid, 'rcnt_nzbi:', rcnt_nzbi
!!$    write(0,*) myid, 'rdsp_nzbi:', rdsp_nzbi
!!$    write(0,*) myid, 'scnt_nzbi:', scnt_nzbi
!!$    write(0,*) myid, 'sdsp_nzbi:', sdsp_nzbi

    ! --- Allocate

    call alloc_sp_ten(ten, rank, glb_shp, SPTEN_FULL, n_nzb_new, n_nzv_new, &
    &                 name)

    ! --- Alltoall

    ! The integer communication can probably contracted to one.

    ! redistribute top
    call MPI_Alltoallv(its%top, scnt_nzbi, sdsp_nzbi, MPI_INTEGER, &
    &                  ten%top, rcnt_nzbi, rdsp_nzbi, MPI_INTEGER, &
    &                  mpi_comm_global, mpierr)
    if (mpierr /= MPI_SUCCESS) call aims_stop('Failed alltoallv on top')

    ! redistribute shp
    call MPI_Alltoallv(its%shp, scnt_nzbi, sdsp_nzbi, MPI_INTEGER, &
    &                  ten%shp, rcnt_nzbi, rdsp_nzbi, MPI_INTEGER, &
    &                  mpi_comm_global, mpierr)
    if (mpierr /= MPI_SUCCESS) call aims_stop('Failed alltoallv on shp')

    ! redistribute val
    call MPI_Alltoallv(its%val, scnt_nzv, sdsp_nzv, MPI_DOUBLE_PRECISION, &
    &                  ten%val,  rcnt_nzv, rdsp_nzv, MPI_DOUBLE_PRECISION, &
    &                  mpi_comm_global, mpierr)
    if (mpierr /= MPI_SUCCESS) call aims_stop('Failed alltoallv on val')
    call construct_off(rank, ten%n_nzb, ten%shp, ten%off)

    ! --- tidy up

    call dealloc_sp_ten(its)
    deallocate(rdsp_nzv, rdsp_nzbi, sdsp_nzv, sdsp_nzbi)
    deallocate(recver2n_nzbv, sender2n_nzbv)
    deallocate(want2task)

    call check_global_aliasing(ten, 'ten (out)', func)
    call check_local_aliasing(ten, 'ten (out)', func)

  end subroutine redistribute_sp_ten
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_scalapack_blocking
  !  NAME
  !    get_scalapack_blocking
  !  SYNOPSIS

  subroutine get_scalapack_blocking(want, rank, glb_shp, i_dim_M, i_dim_N, &
  &                                 MB, NB, nprow, npcol, myprow, mypcol, name)

    !  PURPOSE
    !
    !    Construct a blocking (as a type(sp_ten) with mode=SPTEN_BLOCKING)
    !    corresponding to a Scalapack block-cyclic distribution over the
    !    dimensions i_dim_M and i_dim_N.
    !
    !    Distribution is done only over the specified dimensions.  Please note
    !    that the resulting "want" could only be extended to a "descriptor"
    !    (by adding a valid %off field) for rank==2.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    use scalapack_utils, only: sclpck_num_loc_block, sclpck_top_block, &
        sclpck_shp_block
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: want
    integer, intent(IN) :: rank
    integer, intent(IN) :: glb_shp(rank)
    integer, intent(IN) :: i_dim_M, i_dim_N
    integer, intent(IN) :: MB, NB
    integer, intent(IN) :: nprow, npcol
    integer, intent(IN) :: myprow, mypcol
    character(*), intent(IN) :: name

    !  INPUTS
    !    o rank -- Rank of the tensor (e.g. 3)
    !    o glb_shp -- Global shape of the tensor
    !    o i_dim_M, i_dim_N -- Dimensions of the sparse tensor to distribute
    !    o MB, NB -- Blocking for rows and columns
    !    o nprow, npcol -- Size of 2D processor grid
    !    o myprow, mypcol -- Processor id (0 <= myp... < np...)
    !  OUTPUTS
    !    o want -- [blocking] Sparse tensor structure with mode=SPTEN_BLOCKING
    !              specifying the locally stored blocks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: M, N
    integer :: n_loc_row, n_loc_col, i_loc_row, i_loc_col
    integer :: n_nzb, n_nzb_uptonow
    integer :: top(rank), shp(rank)
    character(*), parameter :: func = 'get_scalapack_blocking'

    if (i_dim_M == i_dim_N) call aims_stop('i_dim_M == i_dim_N', func)
    M = glb_shp(i_dim_M)
    N = glb_shp(i_dim_N)
    n_loc_row = sclpck_num_loc_block(M, MB, nprow, myprow)
    n_loc_col = sclpck_num_loc_block(N, NB, npcol, mypcol)
    n_nzb = n_loc_row * n_loc_col
    call alloc_sp_ten(want, rank, glb_shp, SPTEN_BLOCKING, n_nzb, 0, name)
    shp = glb_shp
    top = 1

    n_nzb_uptonow = 0
    do i_loc_col = 1, n_loc_col
       top(i_dim_N) = sclpck_top_block(i_loc_col, N, NB, npcol, mypcol)
       shp(i_dim_N) = sclpck_shp_block(i_loc_col, N, NB, npcol, mypcol)
       do i_loc_row = 1, n_loc_row
          top(i_dim_M) = sclpck_top_block(i_loc_row, M, MB, nprow, myprow)
          shp(i_dim_M) = sclpck_shp_block(i_loc_row, M, MB, nprow, myprow)
          n_nzb_uptonow = n_nzb_uptonow + 1
          want%top(:, n_nzb_uptonow) = top
          want%shp(:, n_nzb_uptonow) = shp
       end do
    end do
    if (n_nzb /= n_nzb_uptonow) call aims_stop('n_nzb mismatch', func)

  end subroutine get_scalapack_blocking
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_scalapack_descriptor
  !  NAME
  !    get_scalapack_descriptor
  !  SYNOPSIS

  subroutine get_scalapack_descriptor(desc, M, N, LLD, &
  &                                   MB, NB, nprow, npcol, myprow, mypcol, &
  &                                   name)

    !  PURPOSE
    !
    !    Construct a descriptor for a scalapack-like block-cyclic distributed
    !    2D matrix.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    use scalapack_utils, only: sclpck_num_loc_block, sclpck_top_block
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(OUT) :: desc
    integer, intent(IN) :: M, N, LLD
    integer, intent(IN) :: MB, NB
    integer, intent(IN) :: nprow, npcol
    integer, intent(IN) :: myprow, mypcol
    character(*), intent(IN) :: name

    !  INPUTS
    !    o M, N -- Global shape of the matrix
    !    o LLD -- Local leading dimension [ A == A(LLD, *) ]
    !    o MB, NB -- Blocking for rows and columns
    !    o nprow, npcol -- Size of 2D processor grid
    !    o myprow, mypcol -- Processor id (0 <= myp... < np...)
    !  OUTPUTS
    !    o desc -- [descriptor] Sparse tensor structure with
    !              specifying the locally stored blocks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: glb_shp(2), top(2)
    integer :: n_loc_row, n_loc_col, i_loc_row, i_loc_col
    integer :: i_nzb
    integer :: info
    character(*), parameter :: func = 'get_scalapack_descriptor'

    glb_shp(1) = M
    glb_shp(2) = N
    call get_scalapack_blocking(desc, 2, glb_shp, 1, 2, &
    &                           MB, NB, nprow, npcol, myprow, mypcol, name)
    desc%mode = SPTEN_DESCRIPTOR
    allocate(desc%off(desc%n_nzb), desc%str(2, desc%n_nzb), stat=info)
    desc%str(1,:) = 1
    desc%str(2,:) = LLD
    desc%n_nzv = LLD*N

    n_loc_row = sclpck_num_loc_block(M, MB, nprow, myprow)
    n_loc_col = sclpck_num_loc_block(N, NB, npcol, mypcol)
    if (desc%n_nzb /= n_loc_row * n_loc_col) then
       call aims_stop('n_nzb mismatch', func)
    end if

    i_nzb = 0
    do i_loc_col = 1, n_loc_col
       top(2) = sclpck_top_block(i_loc_col, N, NB, npcol, mypcol)
       do i_loc_row = 1, n_loc_row
          top(1) = sclpck_top_block(i_loc_row, M, MB, nprow, myprow)
          i_nzb = i_nzb + 1
          if (any(top /= desc%top(:, i_nzb))) then
             call aims_stop('top mismatch', func)
          end if
          desc%off(i_nzb) = (i_loc_row-1)*MB + (i_loc_col-1)*NB*LLD
       end do
    end do
    if (i_nzb /= desc%n_nzb) call aims_stop('i_nzb mismatch', func)
    call check_local_aliasing(desc, 'desc', func)
    call check_global_aliasing(desc, 'desc', func)

  end subroutine get_scalapack_descriptor
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/extract_matricization
  !  NAME
  !    extract_matricization
  !  SYNOPSIS

  subroutine extract_matricization(ten, val, arr, &
  &                                n_dim_row, dim_row, n_ind_row, ind_row, &
  &                                n_dim_col, dim_col, n_ind_col, ind_col)

    !  PURPOSE
    !
    !    Get an orthogonal "view" of the "matricization" of the global tensor.
    !
    !    In a first (mental) step, the tensor is matricized.  The indices of
    !    the global array are split into "row" and "column" tuples of rank
    !    n_dim_row and n_dim_col, respectively.  The arrays dim_row and
    !    dim_col contain the actual indices (1 <= dim_... <= ten%rank).
    !
    !    In the (mental) step, we are interested only in a subset of n_ind_row
    !    rows and n_ind_col columns.  The actual subsets are given by ind_row
    !    and ind_col.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    real*8, intent(IN) :: val(ten%n_nzv)
    real*8, intent(OUT) :: arr(:,:)
    integer, intent(IN) :: n_dim_row, n_dim_col
    integer, intent(IN) :: dim_row(n_dim_row), dim_col(n_dim_col)
    integer, intent(IN) :: n_ind_row, n_ind_col
    integer, intent(IN) :: ind_row(n_dim_row, n_ind_row)
    integer, intent(IN) :: ind_col(n_dim_col, n_ind_col)

    !  INPUTS
    !    o ten (+ val) -- [full] Tensor containing the data
    !    o n_dim_row, n_dim_col -- Number of global indices for rows and cols.
    !                              n_dim_row + n_dim_col == ten%rank.
    !    o dim_row, dim_col -- Which indices are for rows and for cols
    !                          sort(dim_row // dim_col) == (/1, 2, ..., rank/)
    !    o n_ind_... -- Number of entries for rows/columns
    !    o ind_...... -- Actual entries for rows/columns
    !  OUTPUTS
    !    o arr -- Values of ten at orthogonal view
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank
    integer :: i, i_nzb, i_row, i_col
    logical :: have_dim(ten%rank)
    integer :: n_col_use, n_row_use
    integer, allocatable :: row_use(:), col_use(:), row_off(:), col_off(:)
    integer :: top(ten%rank), shp(ten%rank), str(ten%rank), off
    integer :: i_row_use, i_col_use, i_val
    integer :: info
    character(*), parameter :: func = 'extract_matricization'

    rank = ten%rank
    if (n_dim_row + n_dim_col /= rank) call aims_stop('n_dim mismatch', func)
    have_dim = .false.
    do i = 1, n_dim_row
       have_dim(dim_row(i)) = .true.
    end do
    do i = 1, n_dim_col
       have_dim(dim_col(i)) = .true.
    end do
    if (any(.not. have_dim)) call aims_stop('dim clash', func)

    allocate(row_use(n_ind_row), col_use(n_ind_col), stat=info)
    call check_allocation(info, 'row_use, col_use', func)
    allocate(row_off(n_ind_row), col_off(n_ind_col), stat=info)
    call check_allocation(info, 'row_off, col_off', func)

    arr = 0.d0
    do i_nzb = 1, ten%n_nzb
       top = ten%top(:, i_nzb)
       shp = ten%shp(:, i_nzb)
       call get_strides(ten, i_nzb, str)
       off = ten%off(i_nzb)

       ! Figure out which columns are needed and from where
       call get_rowcol_view_off(rank, n_dim_row, dim_row, n_ind_row, ind_row, &
       &                        top, shp, str, n_row_use, row_use, row_off)
       if (n_row_use == 0) cycle
       ! Figure out which rows are needed and from where
       call get_rowcol_view_off(rank, n_dim_col, dim_col, n_ind_col, ind_col, &
       &                        top, shp, str, n_col_use, col_use, col_off)
       if (n_col_use == 0) cycle

       ! Cycle over columns and rows and fill arr(:,:)
       do i_col_use = 1, n_col_use
          i_col = col_use(i_col_use)
          do i_row_use = 1, n_row_use
             i_row = row_use(i_row_use)
             i_val = off + row_off(i_row_use) + col_off(i_col_use) + 1
             arr(i_row, i_col) = val(i_val)
          end do
       end do
    end do

    deallocate(row_off, col_off, row_use, col_use)

  end subroutine extract_matricization
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/update_matricization
  !  NAME
  !    update_matricization
  !  SYNOPSIS

  subroutine update_matricization(ten, val, arr, &
  &                               n_dim_row, dim_row, n_ind_row, ind_row, &
  &                               n_dim_col, dim_col, n_ind_col, ind_col)

    !  PURPOSE
    !
    !    Use an orthogonal "view" of the "matricization" of the global tensor
    !    to update some of its values.
    !
    !    This is to some extent the reverse operation to
    !    extract_matricization().  But note that the blocking of the tensor is
    !    not updated and that entries not within the view are not touched.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    real*8, intent(INOUT) :: val(ten%n_nzv)
    real*8, intent(IN) :: arr(:,:)
    integer, intent(IN) :: n_dim_row, n_dim_col
    integer, intent(IN) :: dim_row(n_dim_row), dim_col(n_dim_col)
    integer, intent(IN) :: n_ind_row, n_ind_col
    integer, intent(IN) :: ind_row(n_dim_row, n_ind_row)
    integer, intent(IN) :: ind_col(n_dim_col, n_ind_col)

    !  INPUTS
    !    o ten (+ val) -- [full] Tensor containing the data
    !    o arr -- Matricized values to be used to update
    !    o n_dim_row, n_dim_col -- Number of global indices for rows and cols.
    !                              n_dim_row + n_dim_col == ten%rank.
    !    o dim_row, dim_col -- Which indices are for rows and for cols
    !                          sort(dim_row // dim_col) == (/1, 2, ..., rank/)
    !    o n_ind_... -- Number of entries for rows/columns
    !    o ind_...... -- Actual entries for rows/columns
    !  OUTPUTS
    !    o val -- Updated value field of ten
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank
    integer :: i, i_nzb, i_row, i_col
    logical :: have_dim(ten%rank)
    integer :: n_col_use, n_row_use
    integer, allocatable :: row_use(:), col_use(:), row_off(:), col_off(:)
    integer :: top(ten%rank), shp(ten%rank), str(ten%rank), off
    integer :: i_row_use, i_col_use, i_val
    integer :: info
    character(*), parameter :: func = 'update_matricization'

    rank = ten%rank
    if (n_dim_row + n_dim_col /= rank) call aims_stop('n_dim mismatch', func)
    have_dim = .false.
    do i = 1, n_dim_row
       have_dim(dim_row(i)) = .true.
    end do
    do i = 1, n_dim_col
       have_dim(dim_col(i)) = .true.
    end do
    if (any(.not. have_dim)) call aims_stop('dim clash', func)

    allocate(row_use(n_ind_row), col_use(n_ind_col), stat=info)
    call check_allocation(info, 'row_use, col_use', func)
    allocate(row_off(n_ind_row), col_off(n_ind_col), stat=info)
    call check_allocation(info, 'row_off, col_off', func)

    do i_nzb = 1, ten%n_nzb
       top = ten%top(:, i_nzb)
       shp = ten%shp(:, i_nzb)
       call get_strides(ten, i_nzb, str)
       off = ten%off(i_nzb)

       ! Figure out which columns are needed and from where
       call get_rowcol_view_off(rank, n_dim_row, dim_row, n_ind_row, ind_row, &
       &                        top, shp, str, n_row_use, row_use, row_off)
       if (n_row_use == 0) cycle
       ! Figure out which rows are needed and from where
       call get_rowcol_view_off(rank, n_dim_col, dim_col, n_ind_col, ind_col, &
       &                        top, shp, str, n_col_use, col_use, col_off)
       if (n_col_use == 0) cycle

       ! Cycle over columns and rows and fill arr(:,:)
       do i_col_use = 1, n_col_use
          i_col = col_use(i_col_use)
          do i_row_use = 1, n_row_use
             i_row = row_use(i_row_use)
             i_val = off + row_off(i_row_use) + col_off(i_col_use) + 1
             val(i_val) = arr(i_row, i_col)
          end do
       end do
    end do

    deallocate(row_off, col_off, row_use, col_use)

  end subroutine update_matricization
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/update_block
  !  NAME
  !    update_block
  !  SYNOPSIS

  subroutine update_block(ten, val, n_arr, arr, arr_top, arr_shp, arr_str)

    !  PURPOSE
    !
    !    Use a block of the global tensor to update some of its values.
    !
    !    This is to some extent the reverse operation to
    !    extract_block().  But note that the blocking of the tensor is
    !    not updated and that entries not within the view are not touched.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    real*8, intent(INOUT) :: val(ten%n_nzv)
    integer, intent(IN) :: n_arr
    real*8, intent(IN) :: arr(n_arr)
    integer, intent(IN) :: arr_top(ten%rank)
    integer, intent(IN) :: arr_shp(ten%rank)
    integer, intent(IN) :: arr_str(ten%rank)

    !  INPUTS
    !    o ten (+ val) -- [full] Tensor containing the data
    !    o arr -- Values to be used to update
    !    o arr_top -- First element of block (arr(1) -> ten[top])
    !    o arr_shp -- Shape of block to update
    !    o arr_str -- Strides in arr
    !  OUTPUTS
    !    o val -- Updated value field of ten
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank, i_nzb
    integer :: arr_bot(ten%rank)
    integer :: ten_top(ten%rank), ten_shp(ten%rank), ten_bot(ten%rank)
    integer :: ten_str(ten%rank), ten_off
    integer :: its_top(ten%rank), its_shp(ten%rank), its_bot(ten%rank)
    integer :: its_str(ten%rank)
    integer :: i_its, i_ten, i_arr, ind(ten%rank)
    integer :: info
    character(*), parameter :: func = 'update_block'

    rank = ten%rank
    arr_bot = arr_top + arr_shp

    do i_nzb = 1, ten%n_nzb
       ten_top = ten%top(:, i_nzb)
       ten_shp = ten%shp(:, i_nzb)
       ten_bot = ten_top + ten_shp
       its_top = max(ten_top, arr_top)
       its_bot = min(ten_bot, arr_bot)
       its_shp = its_bot - its_top

       if (all(its_shp > 0)) then
          call get_strides(ten, i_nzb, ten_str)
          ten_off = ten%off(i_nzb)
          call construct_strides(rank, its_shp, its_str)  ! Only for looping
          do i_its = 1, product(its_shp)
             call get_global_ind(rank, its_top, its_shp, its_str, ind, i_its)
             if (any(ind <= 0)) call aims_stop('Invalid ind', func)
             call get_local_ind(rank, arr_top, arr_shp, arr_str, ind, i_arr)
             call get_local_ind(rank, ten_top, ten_shp, ten_str, ind, i_ten)
             if (i_arr <= 0) call aims_stop('Invalid i_arr', func)
             if (i_ten <= 0) call aims_stop('Invalid i_ten', func)
             val(ten_off + i_ten) = arr(i_arr)
          end do
       end if
    end do

  end subroutine update_block
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/extract_block
  !  NAME
  !    extract_block
  !  SYNOPSIS

  subroutine extract_block(ten_in, val_in, top, shp, val_out, off, str)

    !  PURPOSE
    !    Extract a block of given position (top) and shape (shp) from
    !    a given sparse input tensor (ten_in).  Beware that the extraction
    !    is *added* to val_out.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_in
    real*8, intent(IN) :: val_in(ten_in%n_nzv)
    integer, intent(IN) :: top(ten_in%rank), shp(ten_in%rank)
    real*8, intent(INOUT) :: val_out(*)
    integer, intent(IN), optional :: off, str(ten_in%rank)

    !  INPUTS
    !    o ten_in -- [full/descriptor] Input sparse tensor
    !    o val_in -- Value array for ten_in (if in doubt, use ten_in%val)
    !    o top, shp -- Block to extract
    !    o val_out -- Contains initial values of block (do not forget to reset)
    !    o str, off -- Optional specification for val_out
    !  OUTPUTS
    !    o val_out -- Output block.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_nzb, i_its, i_in, i_out, ind(ten_in%rank)
    integer :: str_out(ten_in%rank), off_out
    integer :: rank, off_in
    integer :: top_in(ten_in%rank), shp_in(ten_in%rank), str_in(ten_in%rank)
    integer :: top_its(ten_in%rank), shp_its(ten_in%rank), str_its(ten_in%rank)
    integer :: bot_its(ten_in%rank)
    character(*), parameter :: func = 'extract_block'

    rank = ten_in%rank
    if (present(off)) then
       off_out = off
    else
       off_out = 0
    end if
    if (present(str)) then
       str_out = str
    else
       call construct_strides(rank, shp, str_out)
    end if

    do i_nzb = 1, ten_in%n_nzb
       top_in = ten_in%top(:, i_nzb)
       shp_in = ten_in%shp(:, i_nzb)

       top_its = max(top_in, top)
       bot_its = min(top_in+shp_in, top+shp)
       shp_its = bot_its - top_its
       if (any(shp_its <= 0)) cycle

       off_in = ten_in%off(i_nzb)
       call get_strides(ten_in, i_nzb, str_in)

       call construct_strides(rank, shp_its, str_its)  ! Only for looping
       do i_its = 1, product(shp_its)
          call get_global_ind(rank, top_its, shp_its, str_its, ind, i_its)
          if (any(ind <= 0)) call aims_stop('Invalid ind', func)
          call get_local_ind(rank, top_in, shp_in, str_in, ind, i_in)
          if (i_in <= 0) call aims_stop('Invalid i_in', func)
          call get_local_ind(rank, top, shp, str_out, ind, i_out)
          if (i_out <= 0) call aims_stop('Invalid i_out', func)
          val_out(off_out + i_out) = val_out(off_out + i_out) + &
          &                          val_in(off_in + i_in)
       end do
    end do

  end subroutine extract_block
  !******
  !----------------------------------------------------------------------------
  !----------------------------- Larger helpers -------------------------------
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/extract_tensor
  !  NAME
  !    extract_tensor
  !  SYNOPSIS

  subroutine extract_tensor(ten_in, val_in, ten_out, name)

    !  PURPOSE
    !    Extract the portion of ten_in which are wanted in ten_out.
    !    Use the blocking given in ten_out.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_in
    real*8, intent(IN) :: val_in(ten_in%n_nzv)
    type(sp_ten), intent(INOUT) :: ten_out
    character(*), intent(IN) :: name

    !  INPUTS
    !    o ten_in -- [full/descriptor] Input tensor
    !    o val_in -- Value array for ten_in (if in doubt, use ten_in%val)
    !    o ten_out -- [blocking] Wanted values
    !    o name -- new name for ten_out
    !  OUTPUTS
    !    o ten_out -- [full] Output tensor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_nzb, rank
    integer :: top(ten_in%rank), shp(ten_in%rank), off
    character(*), parameter :: func = 'extract_tensor'

    rank = ten_in%rank
    if (SPTEN_HAS_OFF(ten_out%mode)) then
       call aims_stop('ten_out should be blocking only', func)
    end if
    call upgrade_blocking(ten_out, SPTEN_FULL, name)

    do i_nzb = 1, ten_out%n_nzb
       top = ten_out%top(:, i_nzb)
       shp = ten_out%shp(:, i_nzb)
       off = ten_out%off(i_nzb)
       call extract_block(ten_in, val_in, top, shp, ten_out%val, off)
    end do

  end subroutine extract_tensor
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/synchronize_wants
  !  NAME
  !    synchronize_wants
  !  SYNOPSIS

  subroutine synchronize_wants(loc_want, glb_want, task2n_nzb_want, name)

    !  PURPOSE
    !
    !    Collectively broadcast all the local wants.  The resulting blocking
    !    is special in the sense that overlapping blocks may make perfect
    !    sense.
    !
    !  USES
    use mpi_tasks, only: n_tasks, myid, use_mpi
    use synchronize_mpi_basic, only: sync_int_vector
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: loc_want
    type(sp_ten), intent(OUT) :: glb_want
    integer, intent(OUT) :: task2n_nzb_want(0:n_tasks-1)
    character(*), intent(IN) :: name

    !  INPUTS
    !    o loc_want -- sp_ten struct containing all locally wanted blocks
    !  OUTPUTS
    !    o glb_want -- sp_ten struct containing all globally wanted blocks
    !    o task2n_nzb_want -- number of blocks each task wants
    !          as the wanted blocks are stored in task order
    !          this array can be used to reconstruct who wants a particular
    !          block.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: n_nzb_want_glob, myoff, rank
    character(*), parameter :: func = 'synchronize_wants'

    if (.not. use_mpi) then
       call copy_sp_ten(loc_want, glb_want, name)
       task2n_nzb_want(0) = loc_want%n_nzb
       return
    end if

    rank = loc_want%rank

    ! Number of blocks each task wants (for sizes)
    task2n_nzb_want = 0
    task2n_nzb_want(myid) = loc_want%n_nzb
    call sync_int_vector(task2n_nzb_want, n_tasks)
    n_nzb_want_glob = sum(task2n_nzb_want)
    myoff = sum(task2n_nzb_want(0:myid-1))  ! myid == 0  -> myoff=0

    call alloc_sp_ten(glb_want, loc_want%rank, loc_want%glb_shp, &
    &                 SPTEN_BLOCKING, n_nzb_want_glob, 0, name)

    ! Synchronize wants
    glb_want%top = 0
    glb_want%shp = 0
    glb_want%top(:, myoff+1:myoff+loc_want%n_nzb) = loc_want%top
    glb_want%shp(:, myoff+1:myoff+loc_want%n_nzb) = loc_want%shp
    call sync_int_vector(glb_want%top, rank*n_nzb_want_glob)
    call sync_int_vector(glb_want%shp, rank*n_nzb_want_glob)

  end subroutine synchronize_wants
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/intersect_sp_ten
  !  NAME
  !    intersect_sp_ten
  !  SYNOPSIS

  subroutine intersect_sp_ten(ten_1, ten_2, n_nzb_max, n_nzb, its, ptr)

    !  PURPOSE
    !
    !    Find intersection of the arrays of blocks 1 and 2.
    !
    !    Naive brute force O(ten_1%n_nzb*ten_2%n_nzb) implementation.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_1
    type(sp_ten), intent(IN) :: ten_2
    integer, intent(IN) :: n_nzb_max
    integer, intent(OUT) :: n_nzb
    type(sp_ten), intent(OUT) :: its
    integer, intent(OUT) :: ptr(2, n_nzb_max)

    !  INPUTS
    !    o ten_1 -- [blocking] First blocking to intersect
    !    o ten_2 -- [blocking] Second blocking to intersect
    !    o n_nzb_max -- Output dimension (if -1, only n_nzb query)
    !  OUTPUTS
    !    o n_nzb -- Number of intersection blocks (<= n_nzb_1*n_nzb_2)
    !    o its -- [blocking] intersections of ten_1 and ten_2
    !    o ptr -- i_nzb -> i_nzb_1, i_nzb_2
    !    + If (n_nzb > n_nzb_max), its and ptr will not be complete.
    !    + The intersection blocks are sorted by ptr(1,:).
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank, glb_shp(ten_1%rank)
    integer :: i_nzb_1, i_nzb_2, n_match_uptonow
    integer :: top_1(ten_1%rank), top_2(ten_1%rank)
    integer :: bot_1(ten_1%rank), bot_2(ten_1%rank)
    integer :: top_its(ten_1%rank), bot_its(ten_1%rank), shp_its(ten_1%rank)
    character(*), parameter :: func = 'intersect_sp_ten'

    ! In order to reduce the scaling with the number of blocks, one could put
    ! all blocks of one of the structures into a tree of hierarchical boxes,
    ! where only the smallest blocks get to the leave boxes...  One won't get
    ! n log n with this, but at least less than n^2 in many cases.

    rank = ten_1%rank
    glb_shp = ten_1%glb_shp
    if (n_nzb_max >= 0) then
       call alloc_sp_ten(its, rank, glb_shp, SPTEN_BLOCKING, n_nzb_max, 0, &
       &                 'its')
    end if

    n_match_uptonow = 0
    do i_nzb_1 = 1, ten_1%n_nzb
       top_1 = ten_1%top(:, i_nzb_1)
       bot_1 = top_1 + ten_1%shp(:, i_nzb_1)
       do i_nzb_2 = 1, ten_2%n_nzb
          top_2 = ten_2%top(:, i_nzb_2)
          bot_2 = top_2 + ten_2%shp(:, i_nzb_2)
          top_its = max(top_1, top_2)
          bot_its = min(bot_1, bot_2)
          shp_its = bot_its - top_its
          if (all(shp_its > 0)) then
             n_match_uptonow = n_match_uptonow + 1
             if (n_match_uptonow <= n_nzb_max) then
                its%top(:, n_match_uptonow) = top_its
                its%shp(:, n_match_uptonow) = shp_its
                ptr(1, n_match_uptonow) = i_nzb_1
                ptr(2, n_match_uptonow) = i_nzb_2
             end if
          end if
       end do
    end do
    n_nzb = n_match_uptonow
    if (n_nzb <= n_nzb_max) then
       its%n_nzb = n_nzb
    else if (n_nzb_max >= 0) then
       call dealloc_sp_ten(its)
    end if

  end subroutine intersect_sp_ten
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/extract_from_sp_ten
  !  NAME
  !    extract_from_sp_ten
  !  SYNOPSIS

  subroutine extract_from_sp_ten(ten_in, val_in, ten_out, out2in, name)

    !  PURPOSE
    !     Special cased version of extract_tensor() where we know in advance
    !     that only one (and which) block of ten_in is needed for each needed
    !     block of ten_out.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten_in
    real*8, intent(IN) :: val_in(ten_in%n_nzv)
    type(sp_ten), intent(INOUT) :: ten_out
    integer, intent(IN) :: out2in(ten_out%n_nzb)
    character(*), intent(IN) :: name

    !  INPUTS
    !    o ten_in -- [descriptor] Sparse input tensor
    !    o ten_out -- [blocking] structure of output tensor
    !    o out2in -- The i_out-th block of ten_out is a subset of the
    !                out2in(i_out)-th block of ten_in.
    !  OUTPUTS
    !    o ten_out -- [full] Sparse output tensor
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank
    integer :: i_nzb_in, i_nzb_out
    integer :: off_in, off_out
    integer :: top_in(ten_in%rank), shp_in(ten_in%rank), str_in(ten_in%rank)
    integer :: top_out(ten_in%rank), shp_out(ten_in%rank), str_out(ten_in%rank)
    integer :: top_its(ten_in%rank), shp_its(ten_in%rank), str_its(ten_in%rank)
    integer :: bot_its(ten_in%rank), ind(ten_in%rank)
    integer :: i_in, i_out, i_its
    integer :: info
    character(*), parameter :: func = 'extract_from_sp_ten'

    rank = ten_in%rank
    call upgrade_blocking(ten_out, SPTEN_FULL, name)

    do i_nzb_out = 1, ten_out%n_nzb
       off_out = ten_out%off(i_nzb_out)
       top_out = ten_out%top(:, i_nzb_out)
       shp_out = ten_out%shp(:, i_nzb_out)
       call get_strides(ten_out, i_nzb_out, str_out)

       i_nzb_in = out2in(i_nzb_out)
       off_in = ten_in%off(i_nzb_in)
       top_in = ten_in%top(:, i_nzb_in)
       shp_in = ten_in%shp(:, i_nzb_in)
       call get_strides(ten_in, i_nzb_in, str_in)

       top_its = max(top_in, top_out)
       bot_its = min(top_in+shp_in, top_out+shp_out)
       shp_its = bot_its - top_its
       call construct_strides(rank, shp_its, str_its)
       if (any(shp_its <= 0)) call aims_stop('Blocks do not overlap', func)

       do i_its = 1, product(shp_its)
          call get_global_ind(rank, top_its, shp_its, str_its, ind, i_its)
          if (any(ind <= 0)) call aims_stop('Invalid ind', func)
          call get_local_ind(rank, top_in, shp_in, str_in, ind, i_in)
          call get_local_ind(rank, top_out, shp_out, str_out, ind, i_out)
          if (i_in <= 0) cycle
          if (i_out <= 0) cycle
          ten_out%val(off_out + i_out) = val_in(off_in + i_in)
       end do
    end do

  end subroutine extract_from_sp_ten
  !******
  !----------------------------------------------------------------------------
  !------------------------------ small helpers -------------------------------
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/construct_off
  !  NAME
  !    construct_off
  !  SYNOPSIS

  subroutine construct_off(rank, n_nzb, shp, off)

    !  PURPOSE
    !     Generate default off from shp.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: rank, n_nzb
    integer, intent(IN) :: shp(rank, n_nzb)
    integer, intent(OUT) :: off(n_nzb)

    !  INPUTS
    !    o rank -- Rank of tensor
    !    o n_nzb -- Number of non-zero blocks
    !    o shp -- Shapes of non-zero blocks
    !  OUTPUTS
    !    o off -- Offset of non-zero blocks in storage array
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_nzb
    character(*), parameter :: func = 'construct_off'

    if (n_nzb <= 0) return
    off(1) = 0
    do i_nzb = 2, n_nzb
       off(i_nzb) = off(i_nzb-1) + product(shp(:, i_nzb-1))
    end do

  end subroutine construct_off
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/count_nzv
  !  NAME
  !    count_nzv
  !  SYNOPSIS

  subroutine count_nzv(rank, n_nzb, shp, n_nzv)

    !  PURPOSE
    !     Get the number of nonzero values (n_nzv) from the shape (shp).
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: rank
    integer, intent(IN) :: n_nzb
    integer, intent(IN) :: shp(rank, n_nzb)
    integer, intent(OUT) :: n_nzv

    !  INPUTS
    !    o rank -- Rank of tensor
    !    o n_nzb -- Number of nonzero blocks
    !    o shp -- Shape of blocks
    !  OUTPUTS
    !    o n_nzv -- Number of nonzero values
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_nzb, n_nzv_uptonow
    character(*), parameter :: func = 'count_nzv'

    n_nzv_uptonow = 0
    do i_nzb = 1, n_nzb
       n_nzv_uptonow = n_nzv_uptonow + product(shp(:, i_nzb))
    end do
    n_nzv = n_nzv_uptonow

  end subroutine count_nzv
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/construct_strides
  !  NAME
  !    construct_strides
  !  SYNOPSIS

  subroutine construct_strides(rank, shp, str)

    !  PURPOSE
    !    Construct the (Fortran-contiguous) strides from a given shape.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: rank
    integer, intent(IN) :: shp(rank)
    integer, intent(OUT) :: str(rank)

    !  INPUTS
    !    o rank -- Rank of (sub-)tensor
    !    o shp -- Shape of block
    !    o str -- Strides of block if saved contiguous (as it is)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i

    str(1) = 1
    do i = 2, rank
       str(i) = str(i-1) * shp(i-1)
    end do

  end subroutine construct_strides
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_strides
  !  NAME
  !    get_strides
  !  SYNOPSIS

  subroutine get_strides(ten, i_nzb, str)

    !  PURPOSE
    !    Construct the (Fortran-contiguous) strides from a given shape.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    integer, intent(IN) :: i_nzb
    integer, intent(OUT) :: str(ten%rank)

    !  INPUTS
    !    o ten -- Sparse tensor or sparse descriptor
    !    o i_nzb -- Block id
    !  OUTPUTS
    !    o str -- Strides of block
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'get_strides'

    if (.not. SPTEN_HAS_OFF(ten%mode)) then
       call aims_stop('No offset information present', func)
    end if
    if (allocated(ten%str)) then
       str = ten%str(:, i_nzb)
    else
       call construct_strides(ten%rank, ten%shp(:, i_nzb), str)
    end if

  end subroutine get_strides
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_local_ind
  !  NAME
  !    get_local_ind
  !  SYNOPSIS

  subroutine get_local_ind(rank, top, shp, str, ind, i_val)

    !  PURPOSE
    !    Get the (linearized) position in local val(:) array from
    !    global indices ind(:).
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: rank
    integer, intent(IN) :: top(rank), shp(rank), str(rank)
    integer, intent(IN) :: ind(rank)
    integer, intent(OUT) :: i_val

    !  INPUTS
    !    o rank -- Rank of (sub-)tensor
    !    o shp -- Shape of block
    !    o str -- Strides of block if saved contiguous (as it is)
    !    o ind -- Indices into global array
    !  OUTPUTS
    !    o i_val -- Position in local val(:) array
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_dim, off

    i_val = 0
    do i_dim = 1, rank
       off = ind(i_dim) - top(i_dim)   ! 0 ... shp(i_dim)-1
       if (off < 0 .or. off >= shp(i_dim)) then
          i_val = 0
          return
       end if
       i_val = i_val + off * str(i_dim)
    end do
    i_val = i_val + 1  ! Fortran starts with 1.

  end subroutine get_local_ind
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_global_ind
  !  NAME
  !    get_global_ind
  !  SYNOPSIS

  subroutine get_global_ind(rank, top, shp, str, ind, i_val)

    !  PURPOSE
    !    Get the global position from local topex in val(:) array.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: rank
    integer, intent(IN) :: top(rank), shp(rank), str(rank)
    integer, intent(OUT) :: ind(rank)
    integer, intent(IN) :: i_val

    !  INPUTS
    !    o rank -- Rank of (sub-)tensor
    !    o shp -- Shape of block
    !    o str -- Strides of block if saved contiguous (as it is)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_val_m1, off, i_dim

    if (i_val < 1 .or. i_val > product(shp)) then
       ind = 0
       return
    end if

    i_val_m1 = i_val - 1  ! Fortran starts with 1.
    do i_dim = 1, rank
       off = modulo(i_val_m1 / str(i_dim), shp(i_dim))
       ind(i_dim) = off + top(i_dim)
    end do

  end subroutine get_global_ind
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_rowcol_view_off
  !  NAME
  !    get_rowcol_view_off
  !  SYNOPSIS

  subroutine get_rowcol_view_off(rank, n_dim, dim2gdim, n_ind, ind, &
  &                              top, shp, str, n_use, use2ind, use2off)

    !  PURPOSE
    !
    !    This subroutine is used to merge a set of n_dim dimensions (named by
    !    dim2gdim(1:n_dim)) into a single one and for a given block retrieve a
    !    subset of used indices in this dimension and the corresponding
    !    offsets in the ten%val array.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: rank
    integer, intent(IN) :: n_dim
    integer, intent(IN) :: dim2gdim(n_dim)
    integer, intent(IN) :: n_ind
    integer, intent(IN) :: ind(n_dim, n_ind)
    integer, intent(IN) :: top(rank), shp(rank), str(rank)
    integer, intent(OUT) :: n_use
    integer, intent(OUT) :: use2ind(n_ind), use2off(n_ind)

    !  INPUTS
    !    o rank -- Rank of (global) tensor
    !    o n_dim -- [1<=n_dim<=rank] Number of dimensions for this row/col
    !    o dim2gdim -- Actual dimensions for this row/col
    !    o n_dim -- Number of rows/cols
    !    o ind -- Adresses in tensor
    !    o top, shp, str -- Properties of nonzero block in (global) tensor
    !  OUTPUTS
    !    o n_use -- Number of relevant rows/cols
    !    o use2ind -- i_use -> i_ind
    !    o use2off -- Corresponding offsets in %val
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_ind   ! Counter of rows/cols
    integer :: i_dim   ! Counter of dimensions in row/col
    integer :: i_gdim  ! Dimension in global tensor
    integer :: n_use_uptonow      ! Counter of used rows/cols
    integer :: this_off           ! Offset in %val in construction

    character(*), parameter :: func = 'get_rowcol_view_off'

    n_use_uptonow = 0
    use2ind = 0
    THIS_IND: do i_ind = 1, n_ind
       this_off = 0
       do i_dim = 1, n_dim
          i_gdim = dim2gdim(i_dim)
          if (ind(i_dim, i_ind) <  top(i_gdim) .or. &
          &   ind(i_dim, i_ind) >= top(i_gdim) + shp(i_gdim)) then
             cycle THIS_IND
          end if
          this_off = this_off + &
          &          (ind(i_dim, i_ind) - top(i_gdim)) * str(i_gdim)
       end do
       ! If we reach here, this index is actually used.
       n_use_uptonow = n_use_uptonow + 1
       use2ind(n_use_uptonow) = i_ind
       use2off(n_use_uptonow) = this_off
    end do THIS_IND
    n_use = n_use_uptonow

  end subroutine get_rowcol_view_off
  !******
  !----------------------------------------------------------------------------
  !-------------------------------- checkers ----------------------------------
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/check_global_aliasing
  !  NAME
  !    check_global_aliasing
  !  SYNOPSIS

  subroutine check_global_aliasing(ten, name, caller)

    !  PURPOSE
    !
    !     Check for aliasing i.e. entries in the global array which is
    !     referenced multiple times.
    !
    !     As the behavior on tensors with aliasing is currently more or less
    !     ill defined, error out for such a case.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    character(*), intent(IN) :: name
    character(*), intent(IN) :: caller

    !  INPUTS
    !    o ten -- [blocking+] Tensor to check
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: n_nzb_its
    integer, allocatable :: ptr(:,:)
    type(sp_ten) :: its
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'check_global_aliasing'

    allocate(ptr(2, ten%n_nzb), stat=info)
    call check_allocation(info, 'ptr', func)
    call intersect_sp_ten(ten, ten, ten%n_nzb, n_nzb_its, its, ptr)
    if (n_nzb_its > ten%n_nzb) then
       write(info_str, "(7A,A,I7,A,I7,A)") &
       & 'Tensor ', trim(ten%name), '/', trim(name), ' from ', trim(caller), &
       & ': ', &
       & 'Among the', ten%n_nzb, ' blocks there are', n_nzb_its - ten%n_nzb, &
       & ' spurious overlaps.'
       call aims_stop(info_str, func)
    end if

  end subroutine check_global_aliasing
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/check_local_aliasing
  !  NAME
  !    check_local_aliasing
  !  SYNOPSIS

  subroutine check_local_aliasing(ten, name, caller)

    !  PURPOSE
    !
    !     Check for local aliasing i.e. entries in the value array which are
    !     mapped to several global entries.
    !
    !     As the behavior on tensors with aliasing is currently more or less
    !     ill defined, error out for such a case.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    use localorb_io, only: localorb_info, OL_norm, use_unit
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    character(*), intent(IN) :: name
    character(*), intent(IN) :: caller

    !  INPUTS
    !    o ten -- [descriptor+] Tensor to check
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    logical, allocatable :: val_used(:)
    integer :: rank, i_nzb, i, off, i_val, n_aliasing
    integer :: top(ten%rank), shp(ten%rank), str(ten%rank), ind(ten%rank)
    integer :: van_str(ten%rank)
    character*150 :: info_str
    integer :: info, n_nzv_count
    character(*), parameter :: func = 'check_local_aliasing'

    rank = ten%rank
    call count_nzv(rank, ten%n_nzb, ten%shp, n_nzv_count)
    if (ten%mode == SPTEN_BLOCKING) call aims_stop('No pure blocking', func)
    if (ten%n_nzv <= 0 .and. n_nzv_count > 0) then
       call aims_stop('Need valid n_nzv', func)
    end if

    allocate(val_used(ten%n_nzv), stat=info)
    if (info /= 0) then
       call localorb_info('  * Not enough memory for local aliasing check.')
       return
    end if
    val_used = .false.
    n_aliasing = 0
    do i_nzb = 1, ten%n_nzb
       top = ten%top(:, i_nzb)
       shp = ten%shp(:, i_nzb)
       off = ten%off(i_nzb)
       call get_strides(ten, i_nzb, str)
       call construct_strides(rank, shp, van_str)
       do i = 1, product(shp)
          call get_global_ind(rank, top, shp, van_str, ind, i)
          call get_local_ind(rank, top, shp, str, ind, i_val)
          if (val_used(off + i_val)) n_aliasing = n_aliasing + 1
          val_used(i_val) = .true.
       end do
    end do
    if (n_aliasing > 0) then
       ! call debug_output_blocking(use_unit, ten)
       write(info_str, "(7A,A,I7,A,I7,A,I7,A)") &
       & 'Tensor ', trim(ten%name), '/', trim(name), ' from ', trim(caller), &
       & ': ', 'Among the', ten%n_nzv, ' values there are', n_aliasing, &
       & ' excess and ', ten%n_nzv - count(val_used), ' missing references.'
       call aims_stop(info_str, func)
    end if
    deallocate(val_used)

  end subroutine check_local_aliasing
  !******
  !----------------------------------------------------------------------------
  !---------------------------------- debug -----------------------------------
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/get_memory_footprint
  !  NAME
  !    get_memory_footprint
  !  SYNOPSIS

  subroutine get_memory_footprint(ten, n_max_bytes, &
  &                               n_nzv_sum, n_nzv_min, n_nzv_max, n_nzv_avg, &
  &                               n_nzb_sum, n_nzb_min, n_nzb_max, n_nzb_avg)

    !  PURPOSE
    !
    !    Calculate the memory footprint of a sparse tensor.  Use real*8
    !    instead of integer to avoid overflows in the case of 32 bit integers.
    !
    !  USES
    use mpi_tasks, only: n_tasks
    use synchronize_mpi_basic, only: get_max_double, sync_real_number, &
        get_min_double
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    real*8, intent(OUT) :: n_max_bytes
    real*8, intent(OUT), optional :: n_nzv_sum, n_nzv_min, n_nzv_max
    real*8, intent(OUT), optional :: n_nzv_avg
    real*8, intent(OUT), optional :: n_nzb_sum, n_nzb_min, n_nzb_max
    real*8, intent(OUT), optional :: n_nzb_avg
    
    !  INPUTS
    !    o ten -- Sparse tensor
    !  OUTPUTS
    !    o n_max_bytes -- Total storage for this tensor (max over procs)
    !    o n_nzv_sum -- Total number of non-zero entries
    !    o n_nzv_min, n_nzv_max -- minimum and maximum non-zeros on a proc
    !    o n_nzv_avg -- average non-zeros per proc
    !    o n_nzb_sum -- Total number of blocs entries
    !    o n_nzb_min, n_nzb_max -- minimum and maximum block count on a proc
    !    o n_nzb_avg -- average block count per proc
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: my_n_nzv_max, my_n_nzb_max
    real*8 :: my_n_nzv_sum, my_n_nzb_sum
    character(*), parameter :: func = 'get_memory_footprint'

    call get_max_double(my_n_nzv_max, dble(ten%n_nzv))
    call get_max_double(my_n_nzb_max, dble(ten%n_nzb))
    my_n_nzv_sum = ten%n_nzv
    call sync_real_number(my_n_nzv_sum)
    my_n_nzb_sum = ten%n_nzb
    call sync_real_number(my_n_nzb_sum)

    if (present(n_nzv_max)) n_nzv_max = my_n_nzv_max
    if (present(n_nzb_max)) n_nzb_max = my_n_nzb_max
    if (present(n_nzv_sum)) n_nzv_sum = my_n_nzv_sum
    if (present(n_nzb_sum)) n_nzb_sum = my_n_nzb_sum
    if (present(n_nzv_min)) call get_min_double(n_nzv_min, dble(ten%n_nzv))
    if (present(n_nzb_min)) call get_min_double(n_nzb_min, dble(ten%n_nzb))
    if (present(n_nzv_avg)) n_nzv_avg = my_n_nzv_sum / n_tasks
    if (present(n_nzb_avg)) n_nzb_avg = my_n_nzb_sum / n_tasks

    ! Assume 64bit ints.
    n_max_bytes = 8*my_n_nzv_max + 8*(3*ten%rank+1)*my_n_nzb_max

  end subroutine get_memory_footprint
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/debug_output_blocking
  !  NAME
  !    debug_output_blocking
  !  SYNOPSIS

  subroutine debug_output_blocking(unit, want)

    !  PURPOSE
    !    Output sparse tensor as python dictionary.
    !    Should by called ony by one node.
    !  USES
    use mpi_tasks, only: myid
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: unit
    type(sp_ten), intent(IN) :: want

    !  INPUTS
    !    o unit -- Fortran I/O unit to output to (0: stderr, 6: stdout)
    !    o want -- Blocking to output
    !  OUTPUTS
    !    none (writes to unit)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*50 :: buffer
    integer :: rank
    integer :: i_nzb, i_dim
    integer :: top(want%rank), shp(want%rank), str(want%rank)

    rank = want%rank
    write(unit, "('# Tensor ""',A,'""')") trim(want%name)
    do i_nzb = 1, want%n_nzb
       top = want%top(:, i_nzb)
       shp = want%shp(:, i_nzb)
       write(unit, "('myid=',I4)", advance='NO') myid
       do i_dim = 1, rank
          if (shp(i_dim) > 1) then
             write(buffer, "(I6,'..',I0)") top(i_dim), top(i_dim)+shp(i_dim)-1
          else if (shp(i_dim) == 1) then
             write(buffer, "(6X,'  ',I0)") top(i_dim)
          else if (shp(i_dim) == 0) then
             write(buffer, "(6X,' --')")
          else
             write(buffer, "(3X,8('*'))")
          end if
          write(unit, "(A13)", advance='NO') trim(buffer)
       end do
       if (allocated(want%off)) then
          call get_strides(want, i_nzb, str)
          write(unit, "('  ->',I6,'|')", advance='NO') want%off(i_nzb)
          do i_dim = 1, rank
             write(unit, "(I6)", advance='NO') str(i_dim)
          end do
       end if
       write(unit, "()")
    end do

  end subroutine debug_output_blocking
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/debug_local_sparse2python_dok
  !  NAME
  !    debug_local_sparse2python_dok
  !  SYNOPSIS

  subroutine debug_local_sparse2python_dok(filename, varname, ten, val_in)

    !  PURPOSE
    !    Output sparse tensor as python dictionary.
    !    Should by called ony by one node.
    !  USES

    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: filename
    character(*), intent(IN) :: varname
    type(sp_ten), intent(IN) :: ten
    real*8, intent(IN), optional :: val_in(ten%n_nzv)

    !  INPUTS
    !    o filename -- file to write to (gets replaced)
    !    o varname -- variable name   ->  <varname> = [[...]]
    !    o ten -- array to output
    !  OUTPUTS
    !    none (writes to disk)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: rank
    integer :: i_nzb, i_nzv, i_val, i_dim
    integer :: top(ten%rank), shp(ten%rank), ind(ten%rank)
    integer :: loop_str(ten%rank), val_str(ten%rank)

    rank = ten%rank

    open(20, FILE=filename, status='replace')
    write(20,"('# Tensor ""',A,'""')") trim(ten%name)
    write(20,"(A,A)") trim(varname), ' = {'
    do i_nzb = 1, ten%n_nzb
       top = ten%top(:, i_nzb)
       shp = ten%shp(:, i_nzb)
       call construct_strides(ten%rank, ten%shp(:, i_nzb), loop_str)
       if  (SPTEN_HAS_OFF(ten%mode)) call get_strides(ten, i_nzb, val_str)

       write(20,"('  # block no.',I5)") i_nzb
       do i_nzv = 1, product(ten%shp(:, i_nzb))
          call get_global_ind(rank, top, shp, loop_str, ind, i_nzv)
          write(20, "('    (')", advance='NO')
          do i_dim = 1, rank
             write(20, "(I6,',')", advance='NO') ind(i_dim)
          end do
          if (SPTEN_HAS_OFF(ten%mode)) then
             call get_local_ind(ten%rank, top, shp, val_str, ind, i_val)
             i_val = ten%off(i_nzb) + i_val
             if (present(val_in)) then
                write(20, "('): ',ES24.16,', #',I7)") val_in(i_val), i_val
             else if (SPTEN_HAS_VAL(ten%mode)) then
                write(20, "('): ',ES24.16,', #',I7)") ten%val(i_val), i_val
             else
                write(20, "('): None, # i_val:',I7)") i_val
             end if
          else
             write(20, "('): None,')")
          end if
       end do
    end do
    write(20,"('}')")
    close(20)

  end subroutine debug_local_sparse2python_dok
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/debug_global_sparse2python_dok
  !  NAME
  !    debug_global_sparse2python_dok
  !  SYNOPSIS

  subroutine debug_global_sparse2python_dok(filename, varname, ten, val_in)

    !  PURPOSE
    !    Output global distributed sparse tensor as python dictionary.
    !    Must be called collectively.
    !  USES
    use mpi_tasks, only: myid
    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: filename
    character(*), intent(IN) :: varname
    type(sp_ten), intent(IN) :: ten
    real*8, intent(IN), optional :: val_in(ten%n_nzv)

    !  INPUTS
    !    o filename -- file to write to (gets replaced)
    !    o varname -- variable name   ->  <varname> = [[...]]
    !    o ten -- array to output
    !  OUTPUTS
    !    none (writes to disk)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*50 :: glb_name
    type(sp_ten) :: myten, want
    integer :: rank, n_nzb_want

    rank = ten%rank
    write(glb_name, "('debug ',A)") trim(ten%name)
    call copy_sp_ten(ten, myten, glb_name)

    if (myid == 0) then
       n_nzb_want = 1
    else
       n_nzb_want = 0
    end if
    call alloc_sp_ten(want, rank, ten%glb_shp, SPTEN_BLOCKING, n_nzb_want, 0, &
    &                 'full')
    if (myid == 0) then
       want%top = 1
       want%shp(:, 1) = ten%glb_shp
    end if

    write(glb_name, "('global ',A)") trim(ten%name)
    call redistribute_sp_ten(myten, want, glb_name, val_in)

    if (myid == 0) then
       call debug_local_sparse2python_dok(filename, varname, myten)
    end if

  end subroutine debug_global_sparse2python_dok
  !******
  !----------------------------------------------------------------------------
  !****s* sparse_tensor/debug_global_size
  !  NAME
  !    debug_global_size
  !  SYNOPSIS

  subroutine debug_global_size(ten, n_max_bytes)

    !  PURPOSE
    !    In a collective operation, get the minimum, maximum, and average
    !    numbers of blocks and values for the given sparse tensor
    !  USES

    use localorb_io, only: localorb_info
    implicit none

    !  ARGUMENTS

    type(sp_ten), intent(IN) :: ten
    real*8, intent(OUT) :: n_max_bytes

    !  INPUTS
    !    o ten -- sparse tensor
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: n_nzv_min, n_nzv_max, n_nzv_sum
    real*8 :: n_nzb_min, n_nzb_max, n_nzb_sum
    real*8 :: n_nzv_avg, n_nzb_avg
    character*150 :: info_str
    character(*), parameter :: func = 'debug_global_size'

    call get_memory_footprint(ten, n_max_bytes, &
    &                         n_nzv_sum, n_nzv_min, n_nzv_max, n_nzv_avg, &
    &                         n_nzb_sum, n_nzb_min, n_nzb_max, n_nzb_avg)

    write(info_str, "(2X,'| ',3A,I7,A)") &
    & 'Tensor "', trim(ten%name), '" needs', &
    & ceiling(n_max_bytes / 2**20), ' MiB per node.'
    call localorb_info(info_str)

    write(info_str, "(2X,'|   ',A,ES10.2,A,A,ES10.2,A,ES10.2,A,ES10.2)") &
    & 'It stores ', n_nzv_sum, ' numbers. ', &
    & 'Per node: min:', n_nzv_min, '; avg:', n_nzv_avg, '; max:', n_nzv_max
    call localorb_info(info_str)

    write(info_str, "(2X,'|   ',A,ES10.2,A,A,ES10.2,A,ES10.2,A,ES10.2)") &
    & 'It stores ', n_nzb_sum, ' blocks.  ', &
    & 'Per node: min:', n_nzb_min, '; avg:', n_nzb_avg, '; max:', n_nzb_max
    call localorb_info(info_str)

    call localorb_info('')

  end subroutine debug_global_size
  !******
end module sparse_tensor
!******
