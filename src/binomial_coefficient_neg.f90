!****f* FHI-aims/binomial_coefficient_neg
!  NAME
!    binomial_coefficient_neg
!  SYNOPSIS

real*8 function binomial_coefficient_neg(n, m)

!  PURPOSE
!   small function to evaluate a binomial coefficient
!
!    /   \
!    | n |
!    |   |
!    | m |
!    \   / 
!
!
!  USES
!  ARGUMENTS

  implicit none  
  integer, intent(in) :: n
  integer, intent(in) :: m

!  INPUTS
!    o n -- ??????????
!    o m -- ??????????
!  OUTPUT
!    o binomial_coefficient_neg -- ?????????
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
  
  !     counter
  integer :: i_counter
  
  binomial_coefficient_neg = 1.d0
  
  
  if(n.lt.0) then 
     
     nominator = abs(n) + m - 1
     
     if ((m .ge. 0) .and. (m .le. abs(n))) then
        do i_counter = 1, m, 1
           binomial_coefficient_neg = binomial_coefficient_neg * dble(nominator) / dble(i_counter)
           nominator = nominator - 1
        end do
        binomial_coefficient_neg = binomial_coefficient_neg * (-1)**m
     else
        binomial_coefficient_neg = 0.d0 
      endif
      
   else
      
      nominator = n
      
      if ((m .ge. 0) .and. (m .le. n)) then
         do i_counter = 1, m, 1
            binomial_coefficient_neg = binomial_coefficient_neg * dble(nominator) / dble(i_counter)
            nominator = nominator - 1
         end do
      else
         binomial_coefficient_neg = 0.d0 
      end if
   endif
   
 end function binomial_coefficient_neg
!******
