module mpi_dummy

  ! Module that contains all 'mpi-routines' for use on a single-cpu computer. 
  ! should be linked to the NEB code

  integer :: n_tasks
  integer :: myid
  
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


subroutine finalize_mpi()
implicit none
end subroutine finalize_mpi

end module mpi_dummy
