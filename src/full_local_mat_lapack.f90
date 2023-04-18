!****h* FHI-aims/full_local_mat_lapack
!  NAME
!    full_local_mat_lapack
!  SYNOPSIS
module full_local_mat_lapack
!  PURPOSE
!    This module provides routines for converting the full local matrix format
!    into dense matrix formats suitable for LAPACK, and vice versa.
!  USES
  implicit none
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
  private

  public :: get_set_full_local_matrix_lapack
  public :: init_comm_full_local_matrix_lapack
contains
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!****s* lapack_wrapper/get_set_full_local_matrix_lapack
!  NAME
!    get_set_full_local_matrix_lapack
!  SYNOPSIS
  subroutine get_set_full_local_matrix_lapack( matrix, hamiltonian_w, &
             overlap_matrix_w, hamiltonian_w_complex, overlap_matrix_w_complex,&
             which, i_spin, k_point_global )
!  PURPOSE
!    Gets or sets a LAPACK matrix in local-index mode (working with full local matrices)
!    which = 0: Set overlap_matrix_w's from local overlap matrix (in matrix)
!    which = 1: Set hamiltonian_w's from local hamiltonian (in matrix)
!    which = 2: Set matrix from hamiltonian_w's (hamiltonian's must be set to
!               something usefull before call)
!
!    WPH:  To understand this subroutine, it is recommend that you first read
!    and absorb the init_comm_full_local_matrix_lapack subroutine, which sets
!    up the communication for this subroutine.  It is in that subroutine where
!    most of the actual logic is performed to determine which chunks of data
!    goes where and the necessary data structures are defined and set up. This
!    subroutine then performs the heavy-and-stupid lifting of shuffling the data
!    and multiplying by k-point phases when needed.
!
!    TODO:  I do not believe which == 2 is currently working correctly.
!  USES
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use dimensions, only: n_basis, n_periodic, n_spin, n_k_points
    use localorb_io, only: localorb_info
    use pbc_lists, only: cbasis_to_basis, center_to_cell, cbasis_to_center, &
        position_in_hamiltonian, n_cells_in_hamiltonian, k_phase, k_point_loc
    use mpi_tasks, only: n_tasks, mpi_wtime, mpi_logical, mpi_land, &
        mpi_comm_global, mpi_real8, myid, mpi_status_ignore, aims_stop, &
        myid
    use runtime_choices, only: use_alltoall, output_level, real_eigenvectors
    use scalapack_wrapper, only: n_local_matrix_size, n_basis_local, &
        send_mat_displ, send_mat_count, recv_mat_count, recv_mat_displ, &
        send_mat_count_tot, recv_mat_count_tot, i_basis_local, basis_col_limit,&
        basis_col, basis_row_limit, basis_row
    implicit none
!  ARGUMENTS
    real*8 :: matrix(n_local_matrix_size)
    real*8 :: hamiltonian_w   ( n_basis*(n_basis+1)/2, n_spin)
    real*8, intent(out) :: overlap_matrix_w( n_basis*(n_basis+1)/2)
    complex*16 :: hamiltonian_w_complex   ( n_basis*(n_basis+1)/2, n_spin)
    complex*16, intent(out) :: overlap_matrix_w_complex( n_basis*(n_basis+1)/2)
    integer, intent(in) :: which, i_spin
    integer, intent(in) :: k_point_global
!  INPUTS
!    o matrix -- matrix in full local format (which != 2)
!    o hamiltonian_w -- arbitrary real triangular matrix (which == 2)
!    o hamiltonian_w_complex -- arbitrary complex triangular matrix (which == 2)
!    o which -- the operation to perform
!    o i_spin -- the spin channel
!    o k_point_global -- global index for k-point (i.e. for k_phase)
!  OUTPUT
!    o overlap_matrix_w -- the real Kohn-Sham overlap matrix (which == 0)
!    o overlap_matrix_w_complex -- the complex Kohn-Sham overlap matrix (which == 0)
!    o hamiltonian_w -- the real Kohn-Sham Hamiltonian matrix (which == 1)
!    o hamiltonian_w_complex -- the complex Kohn-Sham Hamiltonian matrix (which == 1)
!    o matrix -- arbitrary matrix in full local format (which == 2)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i, j, ip, i_k_point, i_cell_r, i_cell_c
    integer :: i_cell, i_row, i_col, i_diff, i_send, i_recv, i_cnt_mat
    integer :: jj, n_off, n_rows_local, j_row_local(n_basis_local)
    integer :: istat1, istat2, mpierr
    integer :: send_mat_off(0:n_tasks-1)
    real*8, allocatable :: matrix_send(:), matrix_recv(:)
    logical alloc_success, use_alltoall_really
    integer :: i_index_new, i_basis_1, i_basis_2
    integer, allocatable :: index_b_new(:,:)
    real*8 ttt0

    character(*), parameter :: func = "get_set_full_local_matrix_lapack"

    ttt0 = mpi_wtime()

    ! Reset the matrix to be set, check parameter which

    if(which==0) then
       if(real_eigenvectors)then
          overlap_matrix_w(:) = 0
       else
          overlap_matrix_w_complex(:) = 0
       end if
    else if(which==1) then
       if(real_eigenvectors)then
          hamiltonian_w(:,i_spin) = 0
       else
          hamiltonian_w_complex(:,i_spin) = 0
       end if
    else if(which==2) then
       matrix(:) = 0
    else
       call aims_stop("Illegal parameter for 'which'", func)
    endif

    call aims_allocate(index_b_new, n_basis, n_basis, "index_b_new")

    index_b_new = 0
    i_index_new = 0
    do i_basis_2 = 1,n_basis, 1
      do i_basis_1 = 1,i_basis_2,1
        i_index_new = i_index_new + 1
        index_b_new(i_basis_1, i_basis_2) = i_index_new
      end do
    end do

    use_alltoall_really = use_alltoall

    ! When use_alltoall is set we check if we can allocate the memory needed
    ! and if this doesn't succeed we fall back to sendrecv.

    ! The reason why we don't do that always this way but use a config switch is the fact that
    ! the system might be configured with swap space. In this case the allocate may succeed
    ! but the machine will start swapping - which is absolutely nothing we want!

    if(use_alltoall) then

      ! Try to allocate large matrices for using mpi_alltoallv ...
      allocate(matrix_send(send_mat_count_tot), stat=istat1)
      allocate(matrix_recv(recv_mat_count_tot), stat=istat2)

      ! ... check if allocation succeeded on ALL procs, otherways use mpi_sendrecv

      alloc_success = (istat1==0 .and. istat2==0)
      call mpi_allreduce(alloc_success, use_alltoall_really, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_global, mpierr)

      if(.not.use_alltoall_really) then
        ! fall back to sendrecv
        if(allocated(matrix_send)) deallocate(matrix_send)
        if(allocated(matrix_recv)) deallocate(matrix_recv)
        call localorb_info('  *** Not enough memory for using mpi_alltoall, falling back to mpi_sendrecv')
      endif

    endif

    if(.not.use_alltoall_really) then
       allocate(matrix_send(maxval(send_mat_count)))
       allocate(matrix_recv(maxval(recv_mat_count)))
    endif


    if(which<=1 .and. use_alltoall_really) then
      ! Send the complete local matrix to the receivers using 1 mpi_alltoallv call
      ! This is the most effective way but needs a large amount of memory

      ! Put matrix into matrix_send ...
      send_mat_off(:) = send_mat_displ(:)
      do i_k_point = 1, n_k_points
        if (n_periodic == 0) then
          ip = 0
        else
          ip = k_point_loc(1, i_k_point)
        end if

        do i=1,n_basis_local
          i_col = i_basis_local(i)
          do j=1,n_basis_local
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              send_mat_off(ip) = send_mat_off(ip)+1
              if(j<=i) then
                matrix_send(send_mat_off(ip)) = matrix(i*(i-1)/2+j)
              else
                matrix_send(send_mat_off(ip)) = matrix(j*(j-1)/2+i)
              endif
            endif
          enddo
        enddo
      enddo

      ! ... and send it away
      if (n_tasks .eq. 1) then
        matrix_recv(1:send_mat_count(0)) = matrix_send(1:send_mat_count(0))
      else
        call mpi_alltoallv(matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
                           matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
                           mpi_comm_global, mpierr)
      end if
    endif

    if(which>1 .and. use_alltoall_really) matrix_recv(:) = 0

    ! Insert data from local matrix into lapack matrix (which<=1)
    ! or gather data from lapack matrix for putting into local matrix (which>1)
    ! This is done for every remote task separatly for the case that
    ! mpi_alltoallv cannot be used

    do i_diff = 0, n_tasks-1

      i_send = mod(myid+i_diff,n_tasks)         ! Task to which we send data
      i_recv = mod(myid+n_tasks-i_diff,n_tasks) ! Task from which we get data

      if(.not.use_alltoall_really) then

        ! Gather all local rows going to proc i_send
        ! (for which > 1 this is needed below)
        n_rows_local = 0
        do j=1,n_basis_local
          i_row = i_basis_local(j)
          n_rows_local = n_rows_local+1
          j_row_local(n_rows_local) = j
        enddo

        if(which<=1) then

          n_off = 0
          do i=1,n_basis_local
            i_col = i_basis_local(i)
            do jj=1,n_rows_local
              j = j_row_local(jj)
              i_row = i_basis_local(j)
              if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
                n_off = n_off+1
                if(j<=i) then
                  matrix_send(n_off) = matrix(i*(i-1)/2+j)
                else
                  matrix_send(n_off) = matrix(j*(j-1)/2+i)
                endif
              endif
            enddo
          enddo

          ! Gather and send data for remote task i_send,
          ! receive corresponding data from task i_recv

          if (n_tasks .eq. 1) then
             matrix_recv(1:recv_mat_count(i_recv)) = matrix_send(1:send_mat_count(i_send))
          else
            call mpi_sendrecv(matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
                              matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
                              mpi_comm_global, mpi_status_ignore, mpierr)
          end if
        else

          matrix_recv(1:recv_mat_count(i_recv)) = 0

        endif

      endif

      if(use_alltoall_really) then
        i_cnt_mat = recv_mat_displ(i_recv) ! Counter in matrix_recv
      else
        i_cnt_mat = 0 ! Counter in matrix_recv, reset for every task
      endif

      if(n_periodic==0) then

        ! We could also use the code for the periodic case here,
        ! but the code below is much more simple (and hopefully faster!)
        do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

          i_col = basis_col(j)

          do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

            i_row = basis_row(i)
            if(i_row > i_col) exit ! done with this column

            i_cnt_mat = i_cnt_mat+1
            i_index_new = index_b_new(i_row, i_col)

            if(which==0) then
              overlap_matrix_w(i_index_new) = overlap_matrix_w(i_index_new) + matrix_recv(i_cnt_mat)
            else if(which==1) then
              hamiltonian_w(i_index_new, i_spin) = hamiltonian_w(i_index_new, i_spin) + matrix_recv(i_cnt_mat)
            else
              matrix_recv(i_cnt_mat) = hamiltonian_w(i_index_new, i_spin)
            endif
          enddo
        enddo

      else ! periodic case

        do j=basis_col_limit(i_recv)+1,basis_col_limit(i_recv+1)

          i_col = basis_col(j)
          i_cell_c = center_to_cell(Cbasis_to_center(i_col))
          i_basis_2 = Cbasis_to_basis(i_col)

          do i=basis_row_limit(i_recv)+1,basis_row_limit(i_recv+1)

            i_row = basis_row(i)
            i_basis_1 = Cbasis_to_basis(i_row)
            if(i_basis_1 > i_basis_2) cycle
            i_cell_r = center_to_cell(Cbasis_to_center(i_row))
            i_index_new = index_b_new(i_basis_1, i_basis_2)

            i_cell = position_in_hamiltonian(i_cell_c, i_cell_r) ! Attention: position_in_hamiltonian is not symmetric!

            i_cnt_mat = i_cnt_mat+1

            if(i_cell == n_cells_in_hamiltonian) cycle

            if(which==0) then
              if(real_eigenvectors)then
                overlap_matrix_w(i_index_new) = overlap_matrix_w(i_index_new) + dble(k_phase(i_cell,k_point_global)) * matrix_recv(i_cnt_mat)
              else
                overlap_matrix_w_complex(i_index_new) = overlap_matrix_w_complex(i_index_new) + k_phase(i_cell,k_point_global) * matrix_recv(i_cnt_mat)
              end if
            else if(which==1) then
              if(real_eigenvectors)then
                hamiltonian_w(i_index_new, i_spin) = hamiltonian_w(i_index_new, i_spin) + dble(k_phase(i_cell,k_point_global)) * matrix_recv(i_cnt_mat)
              else
                hamiltonian_w_complex(i_index_new, i_spin) = hamiltonian_w_complex(i_index_new, i_spin) + k_phase(i_cell,k_point_global) * matrix_recv(i_cnt_mat)
              end if
            else
              ! The density matrix is always real, so no complex case here
              matrix_recv(i_cnt_mat) = hamiltonian_w(i_index_new, i_spin)*dble(k_phase(i_cell,k_point_global))
            endif
          enddo
        enddo

      endif

      if(which>1 .and. .not.use_alltoall_really) then

        ! Send matrix_recv immediatly back to owner of local matrix
        if (n_tasks .eq. 1) then
          matrix_send(1:send_mat_count(i_send)) = matrix_recv(1:recv_mat_count(i_recv))
        else
          call mpi_sendrecv(matrix_recv, recv_mat_count(i_recv), MPI_REAL8, i_recv, 111, &
                            matrix_send, send_mat_count(i_send), MPI_REAL8, i_send, 111, &
                            mpi_comm_global, mpi_status_ignore, mpierr)
        end if

        n_off = 0
        do i=1,n_basis_local
          i_col = i_basis_local(i)
          do jj=1,n_rows_local
            j = j_row_local(jj)
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              n_off = n_off+1
              ! some of the elements we got are duplicates, we must be careful
              ! not to add them twice!!!!
              if(Cbasis_to_basis(i_row)==Cbasis_to_basis(i_col) .and. j>i) cycle
              if(j<=i) then
                matrix(i*(i-1)/2+j) = matrix(i*(i-1)/2+j) + matrix_send(n_off)
              else
                matrix(j*(j-1)/2+i) = matrix(j*(j-1)/2+i) + matrix_send(n_off)
              endif
            endif
          enddo
        enddo

      endif

    enddo

    if(which>1 .and. use_alltoall_really) then
      ! Send the matrix gathered in matrix_recv to owners of local matrix with mpi_alltoallv
      ! This is the most effective way but needs a large amount of memory
      if (n_tasks .eq. 1) then
        matrix_send(1:send_mat_count(i_send)) = matrix_recv(1:recv_mat_count(i_recv))
      else
        call mpi_alltoallv(matrix_recv, recv_mat_count, recv_mat_displ, MPI_REAL8, &
                           matrix_send, send_mat_count, send_mat_displ, MPI_REAL8, &
                           mpi_comm_global, mpierr)
      end if

      ! Insert received matrix
      send_mat_off(:) = send_mat_displ(:)
      do i_k_point = 1, n_k_points
        if (n_periodic == 0) then
          ip = 0
        else
          ip = k_point_loc(1, i_k_point)
        end if

        do i=1,n_basis_local
          i_col = i_basis_local(i)
          do j=1,n_basis_local
            i_row = i_basis_local(j)
            if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
              send_mat_off(ip) = send_mat_off(ip)+1
              ! some of the elements we got are duplicates, we must be careful
              ! not to add them twice!!!!
              if(Cbasis_to_basis(i_row)==Cbasis_to_basis(i_col) .and. j>i) cycle
              if(j<=i) then
                matrix(i*(i-1)/2+j) = matrix(i*(i-1)/2+j) + matrix_send(send_mat_off(ip))
              else
                matrix(j*(j-1)/2+i) = matrix(j*(j-1)/2+i) + matrix_send(send_mat_off(ip))
              endif
            endif
          enddo
        enddo

      enddo
    endif

    deallocate(matrix_send)
    deallocate(matrix_recv)

    if((myid==0).and.(output_level .ne. 'MD_light' )) &
       print '(a,f13.6,a)','  Time get_set_full_local_matrix_lapack:',mpi_wtime()-ttt0, ' s'

    call aims_deallocate(index_b_new, "index_b_new")

  end subroutine get_set_full_local_matrix_lapack
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!****s* lapack_wrapper/init_comm_full_local_matrix_lapack
!  NAME
!    init_comm_full_local_matrix_lapack
!  SYNOPSIS
  subroutine init_comm_full_local_matrix_lapack(n_basis_local_full, i_basis_local_full)
!  PURPOSE
!    Initializes the communication for get_set_full_local_matrix_lapack
!  USES
    use dimensions, only: n_centers_basis_T, n_k_points, n_periodic
    use mpi_tasks, only: n_tasks, myid, mpi_comm_global, mpi_integer, aims_stop
    use pbc_lists, only: Cbasis_to_basis, k_point_loc
    use synchronize_mpi_basic, only: sync_integer_vector
    use scalapack_wrapper, only: basis_col_limit, basis_row_limit, &
        i_basis_local, basis_row, basis_col, recv_mat_displ, recv_mat_count, &
        send_mat_displ, send_mat_count, n_basis_local, n_local_matrix_size, &
        recv_mat_count_tot, send_mat_count_tot
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

    ! Variables used for setting up receiving data
    integer :: irow, icol

    ! Variables used for setting up sending data
    integer :: i_k_point, i_col, i_row

    ! Generic multipurpose counters and dummy parameters
    integer :: i, j, ip, ncnt, mpierr

    ! The dimension of the local matrix for a given MPI task
    integer :: n_basis_all(0:n_tasks)

    ! The total number of copies of basis elements that the current MPI task
    ! will receive from other MPI tasks (and itself)
    integer :: ncnt_row, ncnt_col

    ! The number of MPI tasks whose local matrix contains a given basis element
    integer :: i_basis(n_centers_basis_T)

    character(*), parameter :: func = "init_comm_full_local_matrix_lapack"

    ! Set up the dimensions and basis functions contained in the local matrix
    ! for *this* MPI tasks (a.k.a. myid)
    n_basis_local = n_basis_local_full
    if(allocated(i_basis_local)) deallocate(i_basis_local)
    allocate(i_basis_local(n_basis_local))
    i_basis_local(:) = i_basis_local_full(:)
    n_local_matrix_size = n_basis_local*(n_basis_local+1)/2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               Receiving Data                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! WPH: We here set up the communication where myid obtains local matrix
    ! elements from other MPI tasks so that it may construct its dense matrix
    ! suitable for LAPACK.
    !
    ! This requires some communication from other MPI tasks to construct various
    ! indexing arrays as, by design, myid has no method to determine the details
    ! of the local matrices for other MPI tasks without asking them.

    ! Get total number of tasks having each basis function

    i_basis(:) = 0
    i_basis(i_basis_local(1:n_basis_local)) = 1
    call sync_integer_vector(i_basis,n_centers_basis_T)

    ! Count how much basis functions will go to my task

    ncnt_row = 0
    ncnt_col = 0
    do i=1,n_centers_basis_T
      ncnt_row = ncnt_row+i_basis(i)
      ncnt_col = ncnt_col+i_basis(i)
    enddo

    ! Allocate arrays for remote basis functions

    if(allocated(basis_row)) deallocate(basis_row)
    if(allocated(basis_col)) deallocate(basis_col)
    allocate(basis_row(ncnt_row))
    allocate(basis_col(ncnt_col))

    if(.not.allocated(basis_col_limit)) allocate(basis_col_limit(0:n_tasks))
    if(.not.allocated(basis_row_limit)) allocate(basis_row_limit(0:n_tasks))
    basis_col_limit(0) = 0
    basis_row_limit(0) = 0

    ! Get all remote basis functions

    ! First, tell every MPI task the dimensions of the local matrices on every
    ! other MPI task
    call mpi_allgather(n_basis_local,1,MPI_INTEGER,n_basis_all,1,MPI_INTEGER,mpi_comm_global,mpierr)

    ! Next, we loop through each MPI task, broadcasting the basis functions of
    ! its local matrices to every other MPI task
    do ip=0,n_tasks-1
      if(ip==myid) i_basis(1:n_basis_local) = i_basis_local(:)
      call mpi_bcast(i_basis,n_basis_all(ip),mpi_integer,ip,mpi_comm_global,mpierr)
      irow = 0
      icol = 0
      do i=1,n_basis_all(ip)
        irow = irow+1
        basis_row(basis_row_limit(ip)+irow) = i_basis(i)
        icol = icol+1
        basis_col(basis_col_limit(ip)+icol) = i_basis(i)
      enddo
      ! Update the offset for the next broadcasting MPI task based on the number
      ! of elements for the current MPI task pair
      basis_row_limit(ip+1) = basis_row_limit(ip) + irow
      basis_col_limit(ip+1) = basis_col_limit(ip) + icol
    enddo

    ! Plausibility checks
    if(basis_row_limit(n_tasks) /= ncnt_row) call aims_stop('INTERNAL ERROR: basis_row_limit', func)
    if(basis_col_limit(n_tasks) /= ncnt_col) call aims_stop('INTERNAL ERROR: basis_col_limit', func)

    ! set up count/displ arrays for matrix send/recv

    if(.not.allocated(send_mat_count)) allocate(send_mat_count(0:n_tasks-1))
    if(.not.allocated(send_mat_displ)) allocate(send_mat_displ(0:n_tasks-1))
    if(.not.allocated(recv_mat_count)) allocate(recv_mat_count(0:n_tasks-1))
    if(.not.allocated(recv_mat_displ)) allocate(recv_mat_displ(0:n_tasks-1))

    !---- what we recv:

    ! Knowing now the indices of the basis elements that will be communicated
    ! between the MPI task pairs, we'll now loop through broadcasting MPI tasks
    ! again to determine the number of matrix elements that will be communicated
    recv_mat_displ(0) = 0
    do ip=0,n_tasks-1
      ncnt = 0
      ! Loop over all pairs of basis elements involved in communication, only
      ! keeping those that lie in the upper triangular portion of the matrix
      do i=basis_col_limit(ip)+1,basis_col_limit(ip+1)
        i_col = basis_col(i)
        do j=basis_row_limit(ip)+1,basis_row_limit(ip+1)
           i_row = basis_row(j)
           if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) ncnt = ncnt+1
        enddo
      enddo
      ! Update the number of matrix elements that will be received by myid from
      ! the current MPI task pair as well as the offset for the next MPI task
      ! pair
      recv_mat_count(ip) = ncnt
      if(ip>0) recv_mat_displ(ip) = recv_mat_displ(ip-1) + recv_mat_count(ip-1)
    enddo
    recv_mat_count_tot = recv_mat_displ(n_tasks-1)+recv_mat_count(n_tasks-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                Sending Data                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! WPH: We here set up the communication where myid sends its local matrix
    ! elements to other MPI tasks so that they may construct their
    ! dense matrix suitable for LAPACK.
    !
    ! This section requires no communication with other MPI tasks due to the
    ! deterministic setup of the dense matrix format.
    !
    ! This section is independent of the "Receiving Data" section and could
    ! easily be moved to its own subroutine.

    !---- what we send

    send_mat_count(:) = 0
    do ip=0,n_tasks-1
       do i=1,n_basis_local
         i_col = i_basis_local(i)
         do j=1,n_basis_local
           i_row = i_basis_local(j)
           ! As always, we only consider the upper triangular matrix elements
           if(Cbasis_to_basis(i_row) <= Cbasis_to_basis(i_col)) then
             send_mat_count(ip) = send_mat_count(ip)+1
           endif
         enddo
       enddo
    enddo

    send_mat_displ(0) = 0
    do ip=1,n_tasks-1
       send_mat_displ(ip) = send_mat_displ(ip-1) + send_mat_count(ip-1)
    enddo
    send_mat_count_tot = send_mat_displ(n_tasks-1)+send_mat_count(n_tasks-1)

  end subroutine init_comm_full_local_matrix_lapack
!******
end module full_local_mat_lapack
