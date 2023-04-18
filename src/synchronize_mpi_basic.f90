!****h* FHI-aims/synchronize_mpi_basic
!  NAME
!    synchronize_mpi_basic
!  SYNOPSIS

module synchronize_mpi_basic

  !  PURPOSE
  !    This module performs the high-level synchronization of MPI-tasks.
  !    It is designed to no as little as possible about FHI-aims data
  !    structures.
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
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

contains

  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_vector
  !  NAME
  !    sync_vector
  !  SYNOPSIS
  subroutine sync_vector(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double vector over all tasks in a given
    !    MPI-communicator in manageable parts.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim
    real*8 :: vector(dim)
    integer, optional :: mpi_comm
    !  INPUTS
    !    o dim -- dimension of the array to be synchronized
    !    o vector -- the array itself
    !    o mpi_comm -- the communicator over which the synchronization
    !                  takes place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks
    !                in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    real*8, allocatable, dimension(:) :: temp_mpi

    integer:: comm, i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 1000000


    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    allocate(temp_mpi(min(dim,max_len)),stat=i)
    call check_allocation(i,'temp_mpi                      ')

    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_vector
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_vector
  !  NAME
  !    sync_vector
  !  SYNOPSIS
  subroutine sync_vector_no_limit(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double vector over all tasks in a given
    !    MPI-communicator in manageable parts.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim
    real*8 :: vector(dim)
    integer, optional :: mpi_comm
    !  INPUTS
    !    o dim -- dimension of the array to be synchronized
    !    o vector -- the array itself
    !    o mpi_comm -- the communicator over which the synchronization
    !                  takes place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks
    !                in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    real*8, allocatable, dimension(:) :: temp_mpi

    integer:: comm, i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 10000000


    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    allocate(temp_mpi(min(dim,max_len)),stat=i)
    call check_allocation(i,'temp_mpi                      ')

    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_vector_no_limit
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_vector_integer
  !  NAME
  !    sync_vector_integer
  !  SYNOPSIS
  subroutine sync_vector_integer(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double vector over all tasks in a given
    !    MPI-communicator in manageable parts.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim
    integer :: vector(dim)
    integer, optional :: mpi_comm
    !  INPUTS
    !    o dim -- dimension of the array to be synchronized
    !    o vector -- the array itself
    !    o mpi_comm -- the communicator over which the synchronization
    !                  takes place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks
    !                in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer, allocatable, dimension(:) :: temp_mpi

    integer:: comm, i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 1000000


    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    allocate(temp_mpi(min(dim,max_len)),stat=i)
    call check_allocation(i,'temp_mpi                      ')

    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_INTEGER, MPI_SUM, comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_vector_integer
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_vector_complex
  !  NAME
  !    sync_complex_vector
  !  SYNOPSIS
  subroutine sync_vector_complex(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double complex vector over all tasks in
    !    a given MPI-communicator in manageable parts.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim
    complex*16 :: vector(dim)
    integer, optional :: mpi_comm
    !  INPUTS
    !    o dim -- dimension of the array to be synchronized
    !    o vector -- the array itself
    !    o mpi_comm -- the communicator over which the synchronization takes place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    complex*16, allocatable, dimension(:) :: temp_mpi

    integer:: comm, i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 500000


    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    allocate(temp_mpi(min(dim,max_len)),stat=i)
    call check_allocation(i, 'temp_mpi                      ')

    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_DOUBLE_COMPLEX, MPI_SUM, comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_vector_complex
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi/collect_complex_vector
  !  NAME
  !    collect_complex_vector
  !  SYNOPSIS
      subroutine collect_vector_complex(vector, results, dim, targed, mpi_comm)
  !  PURPOSE
  !    Synchronize an arbitrary double complex vector over all tasks in a given
  !    MPI-communicator in manageable parts.
      use mpi_tasks
      implicit none
  !  ARGUMENTS
      integer:: dim, targed
      complex*16 :: vector(dim)
      complex*16 :: results(dim)
      integer, optional :: mpi_comm
  !  INPUTS
  !    o dim -- dimension of the array to be synchronized
  !    o vector -- the array itself
  !    o mpi_comm -- the communicator over which the synchronization takes place
  !  OUTPUT
  !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks in mpi_comm
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  !      complex*16, allocatable, dimension(:) :: temp_mpi

      integer:: comm, i, len
      integer :: mpierr

  !     adjust if needed
      integer, parameter :: max_len = 500000


      if (.not.use_mpi) return

      if(present(mpi_comm)) then
         comm = mpi_comm
      else
         comm = mpi_comm_global
      endif

  !      allocate(temp_mpi(max_len),stat=i)
  !      call check_allocation(i, 'temp_mpi                      ')

      if(myid==targed)then
         results = (0.d0,0.d0)
      end if

      do i=1,dim,max_len

         if(dim-i+1 < max_len) then
            len = dim-i+1
         else
            len = max_len
         endif

  !         temp_mpi = (0.0d0,0.0d0)
         call MPI_REDUCE(vector(i), results(i), len, &
              MPI_DOUBLE_COMPLEX, MPI_SUM, targed, comm, mpierr)
  !         vector(i:i-1+len) = temp_mpi(1:len)

      enddo

   !     deallocate(temp_mpi)

    end subroutine collect_vector_complex
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_logical_vector
  !  NAME
  !    sync_logical_vector
  !  SYNOPSIS
  subroutine sync_logical_vector(vector,size__, operation)
    !  PURPOSE
    !    Synchronize an arbitrary logical vector over all tasks
    !    in the global communicator.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer, intent(IN) :: size__
    logical, intent(INOUT) :: vector(size__)
    integer, intent(IN) :: operation

    !  INPUTS
    !    o size__ -- size__ of the array to be synchronized
    !    o vector -- the array itself
    !    o operation -- either SYNC_AND or SYNC_OR
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE with MPI_LAND over
    !                all tasks in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: mpierr,info
    logical,allocatable,dimension(:):: vector_temp
    character(*), parameter :: func = 'sync_logical_vector'


    if (.not.use_mpi)then

       return

    else

       if (operation /= SYNC_AND .and. operation /= SYNC_OR) then
          call aims_stop('invalid operation', func)
       end if
       allocate(vector_temp(size__),stat=info)
       call check_allocation(info, 'vector_temp                   ')

       vector_temp = .true.

       call MPI_ALLREDUCE( vector, vector_temp, &
       size__ , MPI_LOGICAL, &
       operation, mpi_comm_global, mpierr)

       vector = vector_temp
       deallocate(vector_temp)
    end if

  end subroutine  sync_logical_vector
  !******
  !----------------------------------------------------------------------

  !****s* synchronize_mpi_basic/sync_logical
  !  NAME
  !    sync_logical
  !  SYNOPSIS
  subroutine sync_logical(bool, operation)
    !  PURPOSE
    !    Synchronize an arbitrary logical variable over all tasks
    !    in the global communicator.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    logical, intent(INOUT) :: bool
    integer, intent(IN) :: operation

    !  INPUTS
    !    o vector -- the array itself
    !    o operation -- either SYNC_AND or SYNC_OR
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE with 'operation' over
    !                all tasks in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: mpierr
    logical :: bool_temp
    character(*), parameter :: func = 'sync_logical'

    if (.not.use_mpi)then

       return

    else

       if (operation /= SYNC_AND .and. operation /= SYNC_OR) then
          call aims_stop('invalid operation', func)
       end if


       bool_temp = .true.

       call MPI_ALLREDUCE( bool, bool_temp, &
       1, MPI_LOGICAL, &
       operation, mpi_comm_global , mpierr)

       bool = bool_temp

    end if

  end subroutine  sync_logical
  !******
  !----------------------------------------------------------------------



  !****s* synchronize_mpi_basic/sync_integer_vector
  !  NAME
  !    sync_integer_vector
  !  SYNOPSIS
  subroutine sync_integer_vector(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary integer vector over all tasks in a given
    !    MPI-communicator in manageable parts.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim
    integer :: vector(dim)
    integer, optional :: mpi_comm
    !  INPUTS
    !    o dim -- dimension of the array to be synchronized
    !    o vector -- the array itself
    !    o mpi_comm -- the communicator over which the synchronization takes
    !                  place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks
    !                in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer, allocatable, dimension(:) :: temp_mpi

    integer:: comm, i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 1000000


    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    allocate(temp_mpi(min(dim,max_len)),stat=i)
    call check_allocation(i, 'temp_mpi                      ')


    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_INTEGER, MPI_SUM, comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_integer_vector
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_integer
  !  NAME
  !    sync_integer
  !  SYNOPSIS
  subroutine sync_integer(vector, mpi_comm)
    !  PURPOSE
    !    Synchronize an integer over all tasks in a given
    !    MPI-communicator.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer :: vector
    integer, optional :: mpi_comm
    !  INPUTS
    !    o vector -- the integer
    !    o mpi_comm -- the communicator over which the synchronization takes place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer :: temp_mpi

    integer:: comm
    integer :: mpierr

    !     adjust if needed

    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    temp_mpi = 0
    call MPI_ALLREDUCE(vector, temp_mpi, 1, &
    MPI_INTEGER, MPI_SUM, comm, mpierr)
    vector  = temp_mpi


  end subroutine sync_integer
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_long_integer
  !  NAME
  !    sync_long_integer
  !  SYNOPSIS
  subroutine sync_long_integer(vector, mpi_comm)
    !  PURPOSE
    !    Synchronize an long integer (i.e kind=8) over all tasks in a given
    !    MPI-communicator.
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer(kind=8) :: vector
    integer, optional :: mpi_comm
    !  INPUTS
    !    o vector -- the integer
    !    o mpi_comm -- the communicator over which the synchronization takes place
    !  OUTPUT
    !    o vector -- is set to the result of MPI_ALLREDUCE over all tasks in mpi_comm
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer(kind=8) :: temp_mpi

    integer:: comm
    integer :: mpierr

    !     adjust if needed

    if (.not.use_mpi) return

    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = mpi_comm_global
    endif

    temp_mpi = 0
    call MPI_ALLREDUCE(vector, temp_mpi, 1, &
    MPI_INTEGER8, MPI_SUM, comm, mpierr)
    vector  = temp_mpi


  end subroutine sync_long_integer
  !******

  !------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_matrix
  !  NAME
  !    sync_matrix
  !  SYNOPSIS
  subroutine sync_matrix(matrix, dim1,dim2 )
    !  PURPOSE
    !    Synchronize a general double matrix.
    !  ARGUMENTS
    integer:: dim1,dim2
    real*8 :: matrix(dim1, dim2)
    !  INPUTS
    !    o dim1 -- leading dimension of the array
    !    o dim2 -- second dimension of the array
    !    o matrix -- the array itself
    !  OUTPUT
    !    the array is synchronized over tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    call sync_vector(matrix, dim1*dim2)

  end subroutine sync_matrix
  !******
  !------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_vector_scalapac
  !  NAME
  !    sync_vector_scalapac
  !  SYNOPSIS
  subroutine sync_vector_scalapac(matrix, dim1, scalapack_comm)
    !  PURPOSE
    !    Synchronize a double vector for a ScaLAPACK communicator.
    ! USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim1
    real*8 :: matrix(dim1)
    integer :: scalapack_comm
    !  INPUTS
    !    o dim1 -- dimension of the array
    !    o matrix -- the array itself
    !    o scalapack_comm -- the ScaLAPACK communicator
    !  OUTPUT
    !    the array is synchronized over tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    call sync_vector(matrix, dim1, scalapack_comm)

  end subroutine sync_vector_scalapac
  !******
  !------------------------------------------------
  !****s* synchronize_mpi_basic/sync_matrix_complex
  !  NAME
  !    sync_matrix_complex
  !  SYNOPSIS
  subroutine sync_matrix_complex(matrix, dim1,dim2 )
    !  PURPOSE
    !    Synchronize a genral double complex matrix.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer:: dim1,dim2
    complex*16 :: matrix(dim1, dim2)
    !  INPUTS
    !    o dim1 -- leading dimension of the array
    !    o dim2 -- second dimension of the array
    !    o matrix -- the array itself
    !  OUTPUT
    !    the array is synchronized over tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    call sync_vector_complex(matrix,dim1*dim2)

  end subroutine sync_matrix_complex
  !******
  !---------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_one_integer_number
  !  NAME
  !    sync_one_integer_number
  !  SYNOPSIS
  subroutine sync_one_integer(m)
    !  PURPOSE
    !    Synchronize an integer
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer :: m
    !  INPUTS
    !    o m -- the integer to be synchronized
    !  OUTPUT
    !    the integer m is synchronized over all tasks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: mpierr, n

    if (.not.use_mpi) return
    n = 0
    call MPI_ALLREDUCE(m,n,1,MPI_INTEGER,MPI_SUM, &
    mpi_comm_global,mpierr )
    m = n
  end subroutine sync_one_integer
  !******
  !-------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_timing
  !  NAME
  !    sync_timing
  !  SYNOPSIS
  subroutine sync_timing( timing )
    !  PURPOSE
    !    Synchronize a given timing via max-operation.
    !  ARGUMENTS
    real*8 :: timing
    !  INPUTS
    !    o timing -- a timing data
    !  OUTPUT
    !    timing is substituted by a maximum over all tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !     locals
    real*8 :: timing_mpi

    call get_max_double(timing_mpi, timing)
    timing = timing_mpi

  end subroutine sync_timing
  !******
  !------------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_int_vector
  !  NAME
  !    sync_int_vector
  !  SYNOPSIS
  subroutine sync_int_vector(vector,n_dim)
    !  PURPOSE
    !    Synchronize and integer vector
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer n_dim
    integer ::  vector(n_dim)
    !  INPUTS
    !    o n_dim -- dimension of the vector
    !    o vector -- an integer vector
    !  OUTPUT
    !    the vector is synchronized over all tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !      local varibale
    integer :: vector_mpi(n_dim)
    integer :: mpierr

    if (.not.use_mpi) then
       return
    endif

    ! Debug
    ! write(use_unit,*) myid, 'N_DIM = ', n_dim
    ! write(use_unit,*) myid, 'VECTOR LENGTH = ', size(vector)
    ! write(use_unit,*) myid, 'VECTOR_MPI LENGTH = ', size(vector_mpi)

    call MPI_ALLREDUCE(vector, vector_mpi, &
    n_dim, MPI_INTEGER, &
    MPI_SUM, mpi_comm_global, mpierr)

    vector(:) = vector_mpi(:)
    return
  end subroutine sync_int_vector
  !******
  !------------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_real_number
  !  NAME
  !    sync_real_number
  !  SYNOPSIS
  subroutine sync_real_number(real_number)
    !  PURPOSE
    !    Synchronize a double precision real number
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    real*8  ::  real_number
    !  INPUTS
    !    o real_number -- the double to be synchronized
    !  OUTPUT
    !    the double real_number is synchronized over all tasks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !      local varibale
    real*8 ::  real_number_mpi
    integer :: mpierr

    if (.not.use_mpi) then
       return
    endif

    call MPI_ALLREDUCE(real_number, real_number_mpi, &
    1, MPI_DOUBLE_PRECISION, &
    MPI_SUM, mpi_comm_global, mpierr)

    real_number = real_number_mpi

  end subroutine sync_real_number

  !******
  !------------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_complex_number
  !  NAME
  !    sync_complex_number
  !  SYNOPSIS
  subroutine sync_complex_number(complex_number)
    !  PURPOSE
    !    Synchronize a double precision real number
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    complex(kind=8)  ::  complex_number
    !  INPUTS
    !    o complex_number -- the double to be synchronized
    !  OUTPUT
    !    the double complex_number is synchronized over all tasks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !      local varibale
    complex(kind=8) ::  complex_number_mpi
    integer :: mpierr

    if (.not.use_mpi) then
       return
    endif

    call MPI_ALLREDUCE(complex_number, complex_number_mpi, &
    1, MPI_DOUBLE_COMPLEX, &
    MPI_SUM, mpi_comm_global, mpierr)

    complex_number = complex_number_mpi

  end subroutine sync_complex_number
  !******
  !-------------------------------------------------------------
  !****s* synchronize_mpi_basic/get_max_double
  !  NAME
  !    get_max_double
  !  SYNOPSIS
  subroutine get_max_double( max_value, value_part )
    !  PURPOSE
    !    Find the maximal value of double numbers over tasks.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    real*8, intent(OUT) :: max_value
    real*8, intent(IN) :: value_part
    !  INPUTS
    !    o value_part -- the value in a task.
    !  OUTPUT
    !    o max_value -- the maximal value over all tasks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !     locals
    integer :: mpierr

    if (.not.use_mpi) then
       max_value = value_part
       return
    end if
    max_value = 0.0d0

    call MPI_ALLREDUCE(value_part, &
    max_value, 1, MPI_DOUBLE_PRECISION, &
    MPI_MAX, mpi_comm_global, mpierr)

  end subroutine get_max_double
  !******
  !-------------------------------------------------------------
  !****s* synchronize_mpi_basic/get_min_double
  !  NAME
  !    get_min_double
  !  SYNOPSIS
  subroutine get_min_double( min_value, value_part )
    !  PURPOSE
    !    Find the minimal value of double numbers over tasks.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    real*8, intent(OUT) :: min_value
    real*8, intent(IN) :: value_part
    !  INPUTS
    !    o value_part -- the value in a task.
    !  OUTPUT
    !    o min_value -- the minimal value over all tasks.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !     locals
    integer :: mpierr

    if (.not.use_mpi) then
       min_value = value_part
       return
    end if
    min_value = 0.0d0

    call MPI_ALLREDUCE(value_part, &
    min_value, 1, MPI_DOUBLE_PRECISION, &
    MPI_MIN, mpi_comm_global, mpierr)

  end subroutine get_min_double
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/sync_find_max
  !  NAME
  !    sync_find_max
  !  SYNOPSIS
  subroutine sync_find_max(number,max_number )
    !  PURPOSE
    !    Find the maximal integer over the tasks in the global communicator
    !    and broadcast this to all tasks in the communicator.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer, intent(IN) :: number
    integer, intent(OUT) :: max_number
    !  INPUTS
    !    o number -- value on the task
    !  OUTPUT
    !    o max_number -- maximum value of number over the tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer :: mpierr


    if (.not.use_mpi)then

       max_number =  number

    else

       max_number = 0

       call MPI_ALLREDUCE( number, max_number, &
       1 , MPI_INTEGER, &
       MPI_MAX, mpi_comm_global, mpierr)

    end if


  end subroutine sync_find_max
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/bcast_real
  !  NAME
  !    bcast_real
  !  SYNOPSIS
  subroutine bcast_real( value_real, from_here)
    !  PURPOSE
    !    Broadcast a single double number from a given task to all tasks in the
    !    global communicator.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    real*8::  value_real
    integer::  from_here
    !  INPUTS
    !    o value_real -- the value to be broadcast
    !    o from_here -- the index of the broadcasting task
    !  OUTPUT
    !    o value_real -- set on the receiving tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE

    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(value_real, 1, &
    MPI_DOUBLE_PRECISION, &
    from_here, mpi_comm_global, mpierr)


  end subroutine bcast_real
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/bcast_real_vector
  !  NAME
  !    bcast_real_vector
  !  SYNOPSIS
  subroutine bcast_real_vector( vector_real, dim, from_here)
    !  PURPOSE
    !    Broadcast a double precision array from a given task to all tasks in the
    !    global communicator.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    integer :: dim
    real*8 :: vector_real(dim)
    integer ::  from_here
    !  INPUTS
    !    o vector_real -- the array to be broadcast
    !    o dim -- number of elements to be broadcast
    !    o from_here -- the index of the broadcasting task
    !  OUTPUT
    !    o vector_real -- set on the receiving tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    ! SOURCE

    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(vector_real, dim, &
        MPI_DOUBLE_PRECISION, &
        from_here, mpi_comm_global, mpierr)


  end subroutine bcast_real_vector
   !******
   !----------------------------------------------------------------------
   !****s* synchronize_mpi_basic/bcast_integer
   !  NAME
   !    bcast_integer
   !  SYNOPSIS
   subroutine bcast_integer(value_integer, from_here)
   !  USES
      use mpi_tasks
      implicit none
   !  PURPOSE
   !    Broadcast a single integer number from a given task to all
   !    tasks in the global communicator.
   !  ARGUMENTS
      integer, intent(inout) :: value_integer
      integer, intent(in) ::  from_here
   !  INPUTS
   !    o value_integer -- the value to be broadcast
   !    o from_here -- the index of the broadcasting task
   !  OUTPUT
   !    o value_integer -- set on the receiving tasks
   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Development version, FHI-aims (2014).
   !  SOURCE

      integer :: mpierr

      if (.not.use_mpi) return

      call MPI_Bcast(value_integer, 1, &
            MPI_INTEGER, &
            from_here, mpi_comm_global, mpierr)

   end subroutine bcast_integer
   !******
   !----------------------------------------------------------------------
   !****s* synchronize_mpi_basic/bcast_integer_vector
   !  NAME
   !    bcast_integer_vector
   !  SYNOPSIS
   subroutine bcast_integer_vector( vector_integer, dim, from_here)
   !  USES
      use mpi_tasks
      implicit none
   !  PURPOSE
   !    Broadcast an integer type array from a given task to all tasks in the
   !    global communicator.
   !  ARGUMENTS
      integer, intent(in) :: dim
      integer, dimension(dim), intent(inout) :: vector_integer
      integer ::  from_here
   !  INPUTS
   !    o vector_integer -- the array to be broadcast
   !    o dim -- number of elements to be broadcast
   !    o from_here -- the index of the broadcasting task
   !  OUTPUT
   !    o vector_integer -- set on the receiving tasks
   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Development version, FHI-aims (2014).
   !  SOURCE

      integer :: mpierr

      if (.not.use_mpi) return

      call MPI_Bcast(vector_integer, dim, &
            MPI_INTEGER, &
            from_here, mpi_comm_global, mpierr)

   end subroutine bcast_integer_vector
  !******
  !----------------------------------------------------------------------
  !****s* synchronize_mpi_basic/bcast_logical
  !  NAME
  !    bcast_logical
  !  SYNOPSIS
  subroutine bcast_logical( value_logical, from_here)
    !  PURPOSE
    !    Broadcast a single logical from a given task to all tasks in the
    !    global communicator.
    !  USES
    use mpi_tasks
    implicit none
    !  ARGUMENTS
    logical::  value_logical
    integer::  from_here
    !  INPUTS
    !    o value_logical -- the value to be broadcast
    !    o from_here -- the index of the broadcasting task
    !  OUTPUT
    !    o value_logical -- set on the receiving tasks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(value_logical, 1, &
    MPI_LOGICAL, &
    from_here, mpi_comm_global, mpierr)


  end subroutine bcast_logical

!-------------------------------------------------------------
        subroutine send_complex_vector(data, dim, where )
          use mpi_tasks
          implicit none

          integer::dim, where, info
          complex*16:: data(dim)

          call mpi_send(data, dim, MPI_DOUBLE_COMPLEX, where, 0, mpi_comm_global, info)

        end subroutine send_complex_vector
!-------------------------------------------------------------
        subroutine send_real_vector(data, dim, where )
          use mpi_tasks
          implicit none

          integer::dim, where, info
          real*8:: data(dim)

          call mpi_send(data, dim, MPI_DOUBLE_PRECISION, where, 0, mpi_comm_global, info)

        end subroutine send_real_vector
!-------------------------------------------------------------
        subroutine send_integer_vector(data, dim, where )
          use mpi_tasks
          implicit none

          integer :: dim, where, info
          integer :: data(dim)

          call mpi_send(data, dim, MPI_INTEGER, where, 0, mpi_comm_global, info)

        end subroutine send_integer_vector
!-------------------------------------------------------------
        subroutine receive_complex_vector(data, dim, from )
          use mpi_tasks
          implicit none

          integer::dim, from, info
          complex*16:: data(dim)
          integer, dimension(MPI_STATUS_SIZE):: status

          call mpi_recv(data, dim, MPI_DOUBLE_COMPLEX, from, 0, mpi_comm_global, status, info)

        end subroutine receive_complex_vector
!-------------------------------------------------------------
        subroutine send_string(string,dim,where)
          use mpi_tasks
          implicit none

          integer :: dim, where, info
          character(len=dim) :: string

          call mpi_send(string, dim, MPI_CHARACTER, where, 0, mpi_comm_global, info)

        end subroutine send_string
!-------------------------------------------------------------
        subroutine receive_string(string,dim,from)
          use mpi_tasks
          implicit none

          integer :: dim, from, info
          character(len=dim) :: string
          integer, dimension(MPI_STATUS_SIZE):: status

          call mpi_recv(string, dim, MPI_CHARACTER, from, 0, mpi_comm_global, status, info)

        end subroutine receive_string
!-------------------------------------------------------------
        subroutine receive_real_vector(data, dim, from )
          use mpi_tasks
          implicit none

          integer::dim, from, info
          real*8:: data(dim)
          integer, dimension(MPI_STATUS_SIZE):: status

          call mpi_recv(data, dim, MPI_DOUBLE_PRECISION, from, 0, mpi_comm_global, status, info)

        end subroutine receive_real_vector
!-------------------------------------------------------------
        subroutine receive_integer_vector(data, dim, from )
          use mpi_tasks
          implicit none

          integer::dim, from, info
          integer:: data(dim)
          integer, dimension(MPI_STATUS_SIZE):: status

          call mpi_recv(data, dim, MPI_INTEGER, from, 0, mpi_comm_global, status, info)

        end subroutine receive_integer_vector
!------------------------------------------------------------
      subroutine bcast_complex_vector( data, dim, from_here)
      use mpi_tasks
      implicit none

      integer::  from_here, dim
      complex*16::  data(dim)
      integer :: mpierr

      if (.not.use_mpi) return

         call MPI_Bcast(data, dim, &
               MPI_DOUBLE_COMPLEX, &
              from_here, mpi_comm_global, mpierr)


       end subroutine bcast_complex_vector
!-------------------------------------------------------------
        subroutine isend_complex_vector(data, dim, where, request)
          use mpi_tasks
          implicit none

          integer::dim, where, info
          complex*16:: data(dim)
          integer, dimension(MPI_STATUS_SIZE) :: request

          call mpi_isend(data, dim, MPI_DOUBLE_COMPLEX, where, 0, mpi_comm_global, request, info)

        end subroutine isend_complex_vector
!-------------------------------------------------------------
        subroutine isend_wait(request)
          use mpi_tasks
          implicit none

          integer::info
          integer, dimension(MPI_STATUS_SIZE):: status, request

          call mpi_wait(request, status, info)

        end subroutine isend_wait

! MR: I am here adding a couple of routines that are probably duplicated, and
! are solely used for the pimd-wrapper broadcasting procedure.
  subroutine mp_bcast_int(val)
    use mpi_tasks
    implicit none

    integer:: val
    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(val, 1, MPI_INTEGER, 0, mpi_comm_global, mpierr)
  end subroutine

  subroutine mp_bcast_mat(arr)
    use mpi_tasks
    implicit none

    real*8 :: arr(:,:)
    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(arr, size(arr), MPI_DOUBLE_PRECISION, 0, mpi_comm_global, mpierr)
  end subroutine

  subroutine mp_bcast_arr(arr)
    use mpi_tasks
    implicit none

    real*8:: arr(:)
    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(arr, size(arr), MPI_DOUBLE_PRECISION, 0, mpi_comm_global, mpierr)
  end subroutine

  subroutine mp_bcast_char(arr)
    use mpi_tasks
    implicit none

    CHARACTER (len=*):: arr
    integer :: mpierr

    if (.not.use_mpi) return

    call MPI_Bcast(arr, len(arr), MPI_CHARACTER, 0, mpi_comm_global, mpierr)
  end subroutine
! MR: Here I finish adding these routines.
  !******
end module synchronize_mpi_basic
!******
