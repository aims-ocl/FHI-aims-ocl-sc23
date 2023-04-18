!****f* FHI-aims/factorial
!  NAME
!    factorial
!  SYNOPSIS

real*8 elemental function factorial(n)

!  PURPOSE
!   calculate the factorial of n
!  INPUTS
!   o n -- number to be factorialed
!  OUTPUTS
!   o factorial  -- n!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  
  ! input
  integer, intent(in) :: n

  ! counter
  integer :: i_counter

  factorial = 1.d0
  do i_counter = 1, n, 1
     factorial = factorial * dble(i_counter)
  end do

end function factorial
!******
