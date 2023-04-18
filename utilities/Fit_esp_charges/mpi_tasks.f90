! This module provides basic utilities for MPI calculations

module mpi_tasks

  implicit none
  
  integer :: n_tasks
  integer :: myid
  integer :: mpi_comm_global
  include 'mpif.h'
  
contains
  
  subroutine initialize_mpi()
    implicit none
    integer mpierr
    call MPI_INIT(mpierr)
    if (mpierr.ne.MPI_SUCCESS) then
       write(6,*) "Error in MPI initialization."
       write(6,*) "Exiting..."
       stop
    end if
    mpi_comm_global = MPI_COMM_WORLD
    call get_my_task()
  end subroutine initialize_mpi

  subroutine get_my_task()
    implicit none
    integer mpierr
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, mpierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierr)
    if (mpierr.ne.MPI_SUCCESS) then
       write(6,'(1X,A)') "* get_my_task() failed"
       write(6,'(1X,A)') "* Exiting..."
       stop
    end if
  end subroutine get_my_task
  
  subroutine finalize_mpi()
    implicit none
    integer mpierr
    call MPI_FINALIZE(mpierr)
  end subroutine finalize_mpi

  subroutine sync_vector_complex(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double complex vector over all tasks in
    !    a given MPI-communicator in manageable parts.
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



    if(present(mpi_comm)) then
       comm = mpi_comm
    else
       comm = MPI_COMM_WORLD
    endif

    allocate(temp_mpi(min(dim,max_len)),stat=i)

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

  subroutine sync_vector(vector, dim, mpi_comm)
    !  PURPOSE
    !    Synchronize an arbitrary double vector over all tasks in a given 
    !    MPI-communicator in manageable parts.
    implicit none
    !  ARGUMENTS
    integer:: dim
    real*8 :: vector(dim)
    integer :: mpi_comm
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

    integer:: i, len
    integer :: mpierr

    !     adjust if needed
    integer, parameter :: max_len = 1000000

    allocate(temp_mpi(min(dim,max_len)),stat=i)

    do i=1,dim,max_len

       if(dim-i+1 < max_len) then
          len = dim-i+1
       else
          len = max_len
       endif

       call MPI_ALLREDUCE(vector(i), temp_mpi, len, &
       MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, mpierr)
       vector(i:i-1+len) = temp_mpi(1:len)

    enddo

    deallocate(temp_mpi)

  end subroutine sync_vector

  subroutine sync_matrix(matrix, dim1,dim2, mpi_comm)
    !  PURPOSE
    !    Synchronize a general double matrix.
    !  ARGUMENTS
    integer:: dim1,dim2
    real*8 :: matrix(dim1, dim2)
    integer :: mpi_comm
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

    call sync_vector(matrix, dim1*dim2,mpi_comm)

  end subroutine sync_matrix

  subroutine sync_real_number(real_number,mpi_comm)
    !  PURPOSE
    !    Synchronize a double precision real number
    !  ARGUMENTS
    real*8  ::  real_number
    integer :: mpi_comm
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


    call MPI_ALLREDUCE(real_number, real_number_mpi, &
    1, MPI_DOUBLE_PRECISION, &
    MPI_SUM, mpi_comm, mpierr)

    real_number = real_number_mpi

  end subroutine sync_real_number
end module mpi_tasks
