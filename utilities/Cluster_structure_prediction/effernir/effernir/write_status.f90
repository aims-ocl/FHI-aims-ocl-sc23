subroutine write_status(status)
  
  implicit none

  ! imported variables

  ! input
  integer :: status

  open(20, FILE="opt_data.out")

  select case(status)
     
  case(1)
     write (20,'(A)') "# still optimizing"
  case(2)
     write (20,'(A)') "# optimization too slow!!! aborted "
  case(3)
     write (20,'(A)') "# bh-run finished "
  case default

  end select

  close(20)

end subroutine write_status
