subroutine run_driver() 

  use localorb_io,only:use_unit
  
  write(use_unit,*) "* Attention. Your run has called a subroutine intended for "
  write(use_unit,*) "* the PIMD python wrapper."
  write(use_unit,*) "* This functionality must be compiled explicitly into the code. "
  write(use_unit,*) "* Either modify your control.in file, or build/use the code version "
  write(use_unit,*) "* including PIMD wrapper."
  write(use_unit,*) "* Stopping the code for now. "
  stop
      
end subroutine run_driver
!******

