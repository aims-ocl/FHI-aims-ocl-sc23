! This module provides basic utilities for MPI calculations

module mpi_tasks

  implicit none
  
  integer :: n_tasks
  integer :: myid
  
  include 'mpif.h'
  
contains
  
  subroutine initialize_mpi()
    implicit none
    integer mpierr
    call MPI_INIT(mpierr)
    if (mpierr.ne.MPI_SUCCESS) then
       write(use_unit,*) "Error in MPI initialization."
       write(use_unit,*) "Exiting..."
       stop
    end if
    call get_my_task()
  end subroutine initialize_mpi

  subroutine get_my_task()
    implicit none
    integer mpierr
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, mpierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpierr)
    if (mpierr.ne.MPI_SUCCESS) then
       write(use_unit,'(1X,A)') "* get_my_task() failed"
       write(use_unit,'(1X,A)') "* Exiting..."
       stop
    end if
  end subroutine get_my_task
  
  subroutine finalize_mpi()
    implicit none
    integer mpierr
    call MPI_FINALIZE(mpierr)
  end subroutine finalize_mpi
  
end module mpi_tasks
