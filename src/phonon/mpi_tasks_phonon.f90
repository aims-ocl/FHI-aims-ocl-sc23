! This module provides basic utilities for MPI calculations

module mpi_tasks_phonon

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
    call wait_for_all_tasks()
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
  
  subroutine wait_for_all_tasks()
    implicit none
    integer :: mpierr
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  end subroutine wait_for_all_tasks

  subroutine finalize_mpi()
    implicit none
    integer mpierr
    call wait_for_all_tasks()
    call MPI_FINALIZE(mpierr)
  end subroutine finalize_mpi
  
  subroutine sync_vector(vector, dim)
    implicit none
    integer:: dim, mpierr
    real*8, dimension(dim) :: vector, temp_mpi
    
    call MPI_ALLREDUCE(vector, temp_mpi, dim, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

    vector(:) = temp_mpi(:)

    if (mpierr.ne.0) then
       write(use_unit,*) "WARNING: MPI error on task ",myid
       stop
    end if

  end subroutine sync_vector

end module mpi_tasks_phonon
