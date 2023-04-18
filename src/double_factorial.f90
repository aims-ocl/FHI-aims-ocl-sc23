!****f* FHI-aims/double_factorial
!  NAME
!   double_factorial
!  SYNOPSIS

real*8 function double_factorial(n)

!  PURPOSE
!  The function calculates the double factorial n!! of n
!  using double real numbers
!
!  USES

  implicit none

!  ARGUMENTS

  integer :: n

!  INPUTS
!   o n -- number to be factorialed
!
!  OUTPUTS
!   o double_factorial  -- n!! = 1 * 3 * 5 * ... * n (for n odd)
  
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



  ! counter
  integer :: i_counter
  
  double_factorial = 1.d0
  do i_counter = modulo(n,2) , n, 2
     if (i_counter.eq.0) then 
        double_factorial = 1.d0
     else
        double_factorial = double_factorial * dble(i_counter)
     endif
  end do
  
end function double_factorial
!******	
