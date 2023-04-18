
!****f* FHI-aims/binomial_coefficient
!  NAME
!    binomial_coefficient
!  SYNOPSIS

real*8 elemental function binomial_coefficient(n, m)

!  PURPOSE
!  Small function to evaluate a binomial coefficient
!
!     /   \
!     | n |
!     |   |
!     | m |
!     \   / 
!
!
!  USES
!  ARGUMENTS

  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: m

!  INPUTS
!    o n, m -- coefficient parameters
!  OUTPUTS
!    o binomial_coefficient -- resultant value
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



  ! local variables
  integer :: nominator

  ! counter
  integer :: i_counter

  nominator = n
  binomial_coefficient = 1.d0

  if ((m .ge. 0) .and. (m .le. n)) then
     do i_counter = 1, m, 1
        binomial_coefficient = binomial_coefficient * dble(nominator) / dble(i_counter)
        nominator = nominator - 1
     end do
  else
     binomial_coefficient = 0.d0 
  end if
  

end function binomial_coefficient
!******
