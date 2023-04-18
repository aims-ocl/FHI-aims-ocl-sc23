!****s* FHI-aims/heapsort
!  NAME
!    heapsort
!  SYNOPSIS 
subroutine heapsort( array, length )
!  PURPOSE
!    Sorts a given array of double precision reals using the heapsort
!    algorithm.
  implicit none
!  ARGUMENTS
  integer :: length
  real*8, dimension(length) :: array
!  INPUTS
!    o length -- number of elements to sort
!    o array -- array of elements to sort
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


  integer :: index

  call heapify()
  index = length
  do while (index > 1)
     call swap(1, index)
     index = index - 1
     call downheap(1, index)
  end do

!******
contains
!****s* heapsort/heapify
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

    top = INT(dble(length)/2.0)

    do v = top, 1, -1
       call downheap(v, length)
    end do

  end subroutine heapify
!******
!****s* heapsort/downheap
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
!    the array is restored in max-heap order on exit
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
          if (array(w+1) > array(w)) then
             w = w + 1
          end if
       end if
       if (array(v) >= array(w)) then
          return
       else
          call swap(v ,w)
          v = w
          w = 2*v
       end if
    end do
    
  end subroutine downheap
!******
!****s* heapsort/swap
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
!    entries i and j in the array are swapped on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    real*8 :: temp

    temp = array(i)
    array(i) = array(j)
    array(j) = temp

  end subroutine swap
!******
end subroutine heapsort
