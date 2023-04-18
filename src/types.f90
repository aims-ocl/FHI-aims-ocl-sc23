!****h* FHI-aims/types
!  NAME
!    types
!  SYNOPSIS

module types

!  PURPOSE
!    provides kind constants for _basic_ fortran types, nothing else.
!
!  USES
   implicit none

   real*8, parameter, private :: the_aims_double_precision_type = 0.d0

   integer, parameter, public :: dp=kind(the_aims_double_precision_type)

end module
