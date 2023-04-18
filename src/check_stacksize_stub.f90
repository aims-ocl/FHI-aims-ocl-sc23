    subroutine check_stacksize()
      use localorb_io
      character*200 :: message

      write (message,'(A)') "  | Stacksize not measured: no C compiler"
      call localorb_info(message)  
    end  subroutine check_stacksize
