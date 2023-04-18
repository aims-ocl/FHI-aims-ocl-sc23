!  Module vector encapsulates a three-dimensional vector class with
!  some basic functions
!    
!  Subroutines and functions:
!  * sum
!  * difference
!  * display
!  * copy
!
!  R.Gehrke (2005)

module vector_class
      
  type vector
     real*8 :: x
     real*8 :: y
     real*8 :: z
  end type vector
  
  interface operator (+)
     module procedure sum
  end interface
  
  interface operator (-)
     module procedure difference
  end interface

  interface assignment (=)
     module procedure copy
  end interface

  interface operator (*)
     module procedure rescale
  end interface

contains
  type (vector) function sum(a,b) 
  
    implicit none
  
!  imported variables

    type (vector), intent(in) :: a
    type (vector), intent(in) :: b
    
    sum%x = a%x + b%x
    sum%y = a%y + b%y
    sum%z = a%z + b%z
    
  end function sum

  type (vector) function difference(a, b) 
    
    implicit none
    
    !  imported variables
    
    type (vector), intent(in) :: a
    type (vector), intent(in) :: b
    
    difference%x = a%x - b%x
    difference%y = a%y - b%y
    difference%z = a%z - b%z
    
  end function difference

  subroutine copy(a,b)

    implicit none

    !  imported variables

    type (vector), intent(in)  :: b
    type (vector), intent(out) :: a
    
    a%x = b%x
    a%y = b%y
    a%z = b%z
    
  end subroutine copy
  
  subroutine display(a, unit)

    implicit none
    
!  imported variables

    type (vector), intent(in) :: a
    integer, intent(in), optional :: unit

!  local variables
    integer :: local_unit
    
    if (present(unit)) then
       local_unit = unit
    else
       local_unit = 6
    end if

    write(local_unit ,'(1X, E20.10, 1X, E20.10, 1X, E20.10)') &
    a%x, a%y, a%z

  end subroutine display

  type (vector) function rescale(this, c)

!  imported variables

    type (vector), intent(in) :: this
    real*8, intent(in) :: c
    
    rescale%x = this%x * c
    rescale%y = this%y * c
    rescale%z = this%z * c

  end function rescale

  real*8 function norm(this)
    
!  imported variables
    
  type (vector), intent(in) :: this

  norm = sqrt(this%x * this%x + this%y * this%y + this%z * this%z)

  end function norm
  
  real*8 function square_norm(this)
    
!  imported variables
    
    type (vector), intent(in) :: this
    
    square_norm = this%x * this%x + this%y * this%y + this%z * this%z

  end function square_norm

end module vector_class
