!****s* FHI-aims/heapsort_general
!  NAME
!    heapsort_general
!  SYNOPSIS 
subroutine heapsort_general( array, dim1, dim2, icoord )
!  PURPOSE
!    Sorts a given array of double precision reals using a heapsort
!    algorithm. The sort keys are in the second dimension of the array.
!    at the row icoord.
  implicit none
!  ARGUMENTS
  integer :: dim1, dim2
  real*8, dimension(dim1,dim2) :: array
  integer :: icoord
!  INPUTS
!    o dim1 -- first dimension of the array
!    o dim2 -- second dimension of the array
!    o array -- array of elements to sort
!    o icoord -- the row to be used as a sort key
!  OUTPUT
!    the array is sorted on exit
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


  integer :: index, tmp

  call heapify()
  index = dim2
  do while (index > 1)
     call swap(1, index)
     index = index - 1
     call downheap(1, index)
  end do
!******
contains
!****s* heapsort_general/heapify
!  NAME
!    heapify
!  SYNOPSIS 
  subroutine heapify()
!  PURPOSE
!    Creates the original max-heap structure
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    the array is in max-heap order on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: v
    integer :: top

    top = INT(dble(dim2)/2.0)

    do v = top, 1, -1
       call downheap(v, dim2)
    end do

  end subroutine heapify
!******
!****s*  heapsort_general/downheap
!  NAME
!    downheap
!  SYNOPSIS 
  subroutine downheap( index, end )
!  PURPOSE
!    Restores the max-heap structure after extraction.
!  ARGUMENTS
    integer :: index, end
!  INPUTS
!    o index -- top of the head
!    o end -- end of the heap
!  OUTPUT
!    the max-heap order is restored on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: v, w

    v = index
    w = 2*v

    do while (w <= end)
       if (w+1 <= end) then
!!$          tmp = compare_function(array(:,w), array(:,w+1), icoord)
!!$          if (tmp == -1) then
          if (array(icoord,w) < array(icoord,w+1)) then
             w = w + 1
          end if
       end if
!!$       tmp = compare_function(array(:,w), array(:,v), icoord)
!!$       if (tmp == -1 .or. tmp == 0) then
       if (array(icoord,w) <= array(icoord,v)) then
          return
       else
          call swap(v ,w)
          v = w
          w = 2*v
       end if
    end do
    
  end subroutine downheap
!******
!****s*  heapsort_general/swap
!  NAME
!    swap
!  SYNOPSIS 
  subroutine swap(i,j)
!  PURPOSE
!    Swaps two entries in the heap array.
!  ARGUMENTS
    integer :: i,j
!  INPUTS
!    o i -- one of the indeces to swap
!    o j -- the other of the indeces to swap
!  OUTPUT
!    entries i and j are swapped on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    real*8 :: temp(dim1)

    temp = array(:,i)
    array(:,i) = array(:,j)
    array(:,j) = temp

  end subroutine swap
!******


!****s*  heapsort_general/compare_function
!  NAME
!   compare_function
!  SYNOPSIS

  integer function compare_function(a, b, icoord)

!  PURPOSE
!  Compares functions a and b
!
    implicit none
!  ARGUMENTS

   real*8 ::  a(*), b(*)
   integer :: icoord

!  INPUTS
!   o a, b -- compared functions
!   o icoord -- coordinate
!   
!  OUTPUT
!   o compare_function -- results
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




 

    if(a(icoord)<b(icoord)) then
       compare_function = -1
    else if (a(icoord)>b(icoord)) then
       compare_function = 1
    else
       compare_function = 0
    endif
  end function compare_function
!******
  
end subroutine heapsort_general


subroutine heapsort_general_integer( array, dim1, dim2, icoord )
!  PURPOSE
!    Sorts a given array of double precision reals using a heapsort
!    algorithm. The sort keys are in the second dimension of the array.
!    at the row icoord.
  implicit none
!  ARGUMENTS
  integer :: dim1, dim2
  integer, dimension(dim1,dim2) :: array
  integer :: icoord
!  INPUTS
!    o dim1 -- first dimension of the array
!    o dim2 -- second dimension of the array
!    o array -- array of elements to sort
!    o icoord -- the row to be used as a sort key
!  OUTPUT
!    the array is sorted on exit
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


  integer :: index, tmp

  call heapify()
  index = dim2
  do while (index > 1)
     call swap(1, index)
     index = index - 1
     call downheap(1, index)
  end do
!******
contains
!****s* heapsort_general/heapify
!  NAME
!    heapify
!  SYNOPSIS 
  subroutine heapify()
!  PURPOSE
!    Creates the original max-heap structure
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    the array is in max-heap order on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: v
    integer :: top

    top = INT(dble(dim2)/2.0)

    do v = top, 1, -1
       call downheap(v, dim2)
    end do

  end subroutine heapify
!******
!****s*  heapsort_general/downheap
!  NAME
!    downheap
!  SYNOPSIS 
  subroutine downheap( index, end )
!  PURPOSE
!    Restores the max-heap structure after extraction.
!  ARGUMENTS
    integer :: index, end
!  INPUTS
!    o index -- top of the head
!    o end -- end of the heap
!  OUTPUT
!    the max-heap order is restored on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: v, w

    v = index
    w = 2*v

    do while (w <= end)
       if (w+1 <= end) then
!!$          tmp = compare_function(array(:,w), array(:,w+1), icoord)
!!$          if (tmp == -1) then
          if (array(icoord,w) < array(icoord,w+1)) then
             w = w + 1
          end if
       end if
!!$       tmp = compare_function(array(:,w), array(:,v), icoord)
!!$       if (tmp == -1 .or. tmp == 0) then
       if (array(icoord,w) <= array(icoord,v)) then
          return
       else
          call swap(v ,w)
          v = w
          w = 2*v
       end if
    end do
    
  end subroutine downheap
!******
!****s*  heapsort_general/swap
!  NAME
!    swap
!  SYNOPSIS 
  subroutine swap(i,j)
!  PURPOSE
!    Swaps two entries in the heap array.
!  ARGUMENTS
    integer :: i,j
!  INPUTS
!    o i -- one of the indeces to swap
!    o j -- the other of the indeces to swap
!  OUTPUT
!    entries i and j are swapped on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: temp(dim1)

    temp = array(:,i)
    array(:,i) = array(:,j)
    array(:,j) = temp

  end subroutine swap
!******


!****s*  heapsort_general/compare_function
!  NAME
!   compare_function
!  SYNOPSIS

  integer function compare_function(a, b, icoord)

!  PURPOSE
!  Compares functions a and b
!
    implicit none
!  ARGUMENTS

   integer ::  a(*), b(*)
   integer :: icoord

!  INPUTS
!   o a, b -- compared functions
!   o icoord -- coordinate
!   
!  OUTPUT
!   o compare_function -- results
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




 

    if(a(icoord)<b(icoord)) then
       compare_function = -1
    else if (a(icoord)>b(icoord)) then
       compare_function = 1
    else
       compare_function = 0
    endif
  end function compare_function
!******
  
end subroutine heapsort_general_integer
