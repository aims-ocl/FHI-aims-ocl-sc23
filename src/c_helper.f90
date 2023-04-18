!****h* FHI-aims/c_helper
!*  NAME
!*    c_helper
!*  SYNOPSIS
module c_helper
!*  PURPOSE
!*    This module is a collection of data conversion functions for interfacing
!*    between Fortran and C.
!*  USES
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University)
!*  HISTORY
!*    February 2018 - Created.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is
!*    subject to the terms and conditions of the respective license
!*    agreement.
!*  SOURCE

  public :: fort_log_to_c_int
  public :: c_char_array_to_fort_str

contains

  ! Converts a Fortran logical variable into a C integer following the usual
  ! false-if-zero convention.
  pure function fort_log_to_c_int(fort_log) result(c_int_val)
    use, intrinsic :: iso_c_binding, only: c_int32_t
    implicit none

    logical, intent(in) :: fort_log
    integer(kind=c_int32_t) :: c_int_val

    if (.not. fort_log) then
      c_int_val = 0
    else
      c_int_val = 1
    end if
  end function fort_log_to_c_int

  ! Converts a C character array of known length to a Fortran (dynamic) string
  pure function c_char_array_to_fort_str(len_array, c_char_array) &
       result(fort_str)
    use, intrinsic :: iso_c_binding, only: c_char, C_NULL_CHAR
    implicit none

    integer, intent(in) :: len_array
    character(kind=c_char, len=1), dimension(len_array), intent(in) :: &
         c_char_array
    character(len=:), allocatable :: fort_str

    character(len=len_array) :: fort_str_temp
    integer :: i_char

    ! The following was shamelessly stolen from Stack Overflow (M. S. B.'s
    ! answer to https://stackoverflow.com/questions/8207997/calling-a-fortran-subroutine-from-c/8208960)
    fort_str_temp = " "
    do i_char = 1, len_array
      if (c_char_array(i_char) == C_NULL_CHAR) then
        exit
      else
        fort_str_temp(i_char:i_char) = c_char_array(i_char)
      end if
    end do
    fort_str = trim(adjustl(fort_str_temp))
  end function c_char_array_to_fort_str

end module
!******
