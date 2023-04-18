module mpi_dummy

  ! Module that contains all 'mpi-routines' for use on a single-cpu computer. 
  ! should be linked to the NEB code

  integer :: n_tasks = 1
  integer :: myid = 0
  
contains

subroutine initialize_mpi()
implicit none
integer mpierr
mpierr = -1
call get_my_task()
end subroutine initialize_mpi

subroutine get_my_task()
implicit none
n_tasks = 1
myid    = 0
end subroutine get_my_task

subroutine wait_for_all_tasks()
  implicit none
end subroutine wait_for_all_tasks

subroutine finalize_mpi()
implicit none
end subroutine finalize_mpi

subroutine sync_vector(vector, dim)
  implicit none
  integer :: dim
  real*8, dimension(dim) :: vector
end subroutine sync_vector

end module mpi_dummy
