!****h* FHI-aims/ifcore
!  NAME
!    ifcore -- Stub module for the Intel Fortran ifcore module
!  SYNOPSIS
module ifcore
!  PURPOSE
!    Serves as a stub module for ifcore, an Intel Fortran-specific module.
!  USES
    implicit none
!  AUTHOR
!    William Huhn
!  HISTORY
!    June 2018 - Added.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
contains
  subroutine tracebackqq(a_string, user_exit_code, status, eptr)

    use localorb_io, only: use_unit
    implicit none

    character(len=*), intent(in), optional :: a_string
    integer, intent(in), optional :: user_exit_code
    integer, intent(in), optional :: status
    integer, pointer, intent(in), optional :: eptr

    write(use_unit,*) "This is a stub TRACEBACKQQ call."
    write(use_unit,*) "You should not be here. Check your"
    write(use_unit,*) "setting on the print_stacktrace"
    write(use_unit,*) "variable, or include the ifcore"
    write(use_unit,*) "module during compilation"
    stop

  end subroutine tracebackqq
end module ifcore
