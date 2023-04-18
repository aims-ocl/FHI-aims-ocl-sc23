!  Module vector_array encapsulates an array of vectors with
!  some basic functions
!    
!  Subroutines and functions:
!  * allocate_array
!  * deallocate_array
!  * copy_array
!  * display_array
!
!  R.Gehrke (2005)

module vector_array_class

  use vector_class

  type vector_array
     
     type (vector), dimension(:), pointer :: coords
     
  end type vector_array
  
  interface assignment (=)
     module procedure copy_array
  end interface

  contains

    subroutine allocate_array(this, n_elements)

      implicit none

!  imported variables

      integer :: n_elements
      type (vector_array) :: this

      allocate(this%coords(n_elements))

    end subroutine allocate_array

    subroutine deallocate_array(this)
      
      implicit none

!  imported variables

      type (vector_array) :: this
      
      if (associated(this%coords)) then
         deallocate(this%coords)
      end if
!FIXME: size is NOT set to zero when memory is deallocated!!

    end subroutine deallocate_array

    subroutine copy_array(a,b)

      implicit none

!  imported variables

      type (vector_array), intent(out) :: a
      type (vector_array), intent(in)  :: b

!  local variables

      integer :: i_counter
      
      if (.not.associated(a%coords)) then
         allocate(a%coords(size(b%coords)))
      else if (size(a%coords) .ne. size(b%coords)) then
            deallocate(a%coords) 
            allocate(a%coords(size(b%coords)))
         end if
         
      do i_counter = 1, size(a%coords), 1
         a%coords(i_counter) = b%coords(i_counter)
      end do

    end subroutine copy_array

    subroutine display_array(this, unit)
      
      implicit none

!  imported variables

      type (vector_array), intent(in) :: this
      integer, intent(in), optional :: unit

!  local variables

      integer :: i_counter
      integer :: local_unit
      if (present(unit)) then
         local_unit = unit
      else
         local_unit = 6
      end if
      
      if (associated(this%coords)) then
         do i_counter = 1, size(this%coords), 1
            call display(this%coords(i_counter), local_unit)
         end do
      end if

    end subroutine display_array

    real*8 function array_norm(this)

      ! imported variables
      type (vector_array), intent(in) :: this

      ! local variables
      real*8 :: temp

      ! counter
      integer :: i_counter
      
      temp = 0.d0
      if (associated(this%coords)) then
         do i_counter = 1, size(this%coords), 1
            temp = temp + square_norm(this%coords(i_counter))
         end do
      end if

      array_norm = sqrt(temp)

    end function array_norm

end module vector_array_class
