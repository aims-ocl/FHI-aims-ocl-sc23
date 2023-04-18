! well, factorial of an integer...
!
! R.Gehrke (2006)

real*8 function factorial(n)

  implicit none
  
  ! input
  integer :: n

  ! counter
  integer :: i_counter

  factorial = 1.d0
  do i_counter = 1, n, 1
     factorial = factorial * dble(i_counter)
  end do

end function factorial
