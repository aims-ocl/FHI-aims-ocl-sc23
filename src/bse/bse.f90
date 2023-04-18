subroutine bse( )
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  use mpi_tasks
  use timing
  use hartree_fock
  implicit none
  character(*), parameter :: func = 'bse'

  if(myid == 0) print*, '----------------------------------------------------------------------'
  if(myid == 0) print*, 'BSE calculation starts...'
  ! Time total
  if(n_tasks > 1) then
    call bse_parallel_wrapper()
  else
    call bse_serial_wrapper()
  end if
  if(myid == 0) print*, 'End subroutine bse'
end subroutine bse
