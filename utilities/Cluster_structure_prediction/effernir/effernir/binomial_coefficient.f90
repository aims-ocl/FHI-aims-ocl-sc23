! small function to evaluate a binomial coefficient
!
! /   \
! | n |
! |   |
! | m |
! \   / 
!
! R. Gehrke (2006)
!

real*8 function binomial_coefficient(n, m)

  implicit none
  
  ! imported variables
  integer, intent(in) :: n
  integer, intent(in) :: m

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
