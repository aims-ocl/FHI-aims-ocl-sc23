!****h* FHI-aims/quicksort_index
!  NAME
!    quicksort_index
!  SYNOPSIS

module quicksort_index

!  PURPOSE
!     This module contains subroutines to sort a list of keys according to
!     the so-called quicksort algorithm
!  AUTHOR
!     Markus Sinstein
!  HISTORY
!     Development version, FHI-aims (2014).
!  COPYRIGHT
!     Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!     e.V. Please note that any use of the "FHI-aims-Software" is subject to
!     the terms and conditions of the respective license agreement."
!  HISTORY
!     Development version, FHI-aims (2014).
!  SOURCE
   
   implicit none

   public :: dquicksort_indexlist
   
   private
   
contains


!****** 
!---------------------------------------------------------------------------
!****s* quicksort_index/dquicksort_indexlist
!  NAME
!    dquicksort_indexlist
!  SYNOPSIS

subroutine dquicksort_indexlist(dim, key, ind, first, last)

   implicit none
   
!  PURPOSE
!
!     The subroutine sorts a list of keys via the quicksort algorithm.
!     Sorting does not happen in-place but rather an index list is returned
!     which contains the indices of all keys in ascending order of their values.
!
!     The algorithm has been taken from
!     Numerical Recipes, 3rd Ed., pp. 430-431
!     and translated to Fortran.
!
!  USES
!     nothing
!  ARGUMENTS

   integer, intent(in) :: dim
   real*8, dimension(dim), intent(in) :: key
   integer, dimension(dim), intent(inout) :: ind
   integer, intent(in), optional :: first
   integer, intent(in), optional :: last

!  INPUTS
!  o dim -- number of entries in key and ind list
!  o key -- list of keys
!  o first -- and
!  o last -- mark the range in the index list ind where sorting
!            takes place. By default, the whole range is sorted.
!
!  OUTPUT
!  o ind -- list of key indices, sorted in ascending key values
!
!  AUTHOR
!     Markus Sinstein
!  COPYRIGHT
!     Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!     e.V. Please note that any use of the "FHI-aims-Software" is subject to
!     the terms and conditions of the respective license agreement."
!  HISTORY
!     Development version, FHI-aims (2014).
!  SOURCE

   ! set stack size
   integer, parameter :: NSTACK = 64
   integer :: istack(NSTACK)
   
   ! determine max size for insertion sort fallback
   integer, parameter :: ISORT = 7
   
   ! pivoting variables
   double precision :: keypiv
   integer :: indpiv
   
   integer :: l, r, i, j, k, jstack
   
   if ( present(first) ) then
      if (first.lt.1) then
         l = 1
      else
         l = first
      endif
   else
      l = 1    ! left starting point is first entry by default
   endif

   if ( present(last) ) then
      if (last.gt.dim) then
         r = dim
      else
         r = last
      endif
   else
      r = dim  ! right starting point is last entry by default
   endif
   
   jstack = 0
   
   outerloop: do
      if (r-l .lt. ISORT) then ! insertion sort
         do j = l+1, r
            indpiv = ind(j)
            keypiv = key(indpiv)
            do i = j-1, l, -1
               if ( key(ind(i)) .le. keypiv ) exit
               ind(i+1) = ind(i)
            end do
            ind(i+1) = indpiv
         end do
        
         if ( jstack .eq. 0 ) exit
         r = istack(jstack)
         l = istack(jstack-1)
         jstack = jstack-2
      else                    ! quicksort
         k = (l+r)/2
         call iswap(ind(k), ind(l+1))
         if ( key(ind(l)) .gt. key(ind(r)) ) call iswap( ind(l), ind(r) )
         if ( key(ind(l+1)) .gt. key(ind(r)) ) call iswap( ind(l+1), ind(r) )
         if ( key(ind(l)) .gt. key(ind(l+1)) ) call iswap( ind(l), ind(l+1) )
        
         i = l+1
         j = r
         indpiv = ind(i)
         keypiv = key(indpiv)
         do
            do
               i = i+1
               if ( key(ind(i)) .ge. keypiv ) exit
            end do
            do
               j = j-1
               if ( key(ind(j)) .le. keypiv ) exit
            end do
            if ( j .lt. i ) exit
            call iswap(ind(i), ind(j))
         end do
         ind(l+1) = ind(j)
         ind(j) = indpiv
         jstack = jstack + 2
        
         if ( jstack .gt. NSTACK ) then
            print *, "qs_indexing: stack size too small"
            return
         end if
        
         if ( r-i+1 .ge. j-1 ) then
            istack(jstack) = r
            istack(jstack-1) = i
            r = j-1
         else
            istack(jstack) = j-1
            istack(jstack-1) = l
            l = i
         end if
      end if
   
   end do outerloop
   
end subroutine dquicksort_indexlist


subroutine iswap(a, b)
   integer, intent(inout) :: a, b
   integer :: temp
   temp = a
   a = b
   b = temp
end subroutine iswap


end module quicksort_index

