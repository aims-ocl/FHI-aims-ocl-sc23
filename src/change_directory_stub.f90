    subroutine change_directory(path, error)
      use localorb_io
      use mpi_tasks, only: aims_stop
      implicit none
      character(*), intent(in) :: path
      integer, intent(out) :: error
      character*200 :: message
      
      !TZ intent(out) type do not have an explicit value
      error = 1
      write (message,'(A)') "  | FORTRAN chdir() compiler dependent"
      call localorb_info(message)  
      write (message,'(A)') "  | C compiler needed for unified C chdir()"
      call localorb_info(message)  
      call aims_stop("Aims stops here", "change_directory")
    end subroutine change_directory
