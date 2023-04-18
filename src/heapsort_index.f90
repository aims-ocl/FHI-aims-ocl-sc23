!****s* FHI-aims/heapsort_index
!  NAME
!    heapsort_index
!  SYNOPSIS 
subroutine heapsort_index(iarray, rank, num, idim, new2old, reset)
  !  PURPOSE
  !    Sorts a given array of integer indices using a heapsort
  !    algorithm. The sort keys are in the second dimension of the iarray
  !    at the row idim.
  implicit none
  !  ARGUMENTS
  !f2py integer, intent(OUT) :: new2old(num)
  !f2py logical, intent(IN), optional :: reset = 1
  integer, intent(IN) :: rank, num
  integer, intent(IN) :: iarray(rank, num)
  integer, intent(IN) :: idim
  integer :: new2old(num)
  logical, intent(IN) :: reset
  !  INPUTS
  !    o rank -- first dimension of the iarray
  !    o num -- second dimension of the iarray (number of sort items)
  !    o iarray -- array of elements to sort
  !    o idim -- the row (first index) to be used as a sort key;
  !              at tie, the next cyclic index is used
  !    o new2old -- if (.not. reset), the initial permutation of iarray
  !    o reset -- reset new2old(i) to i?
  !  OUTPUT
  !    o new2old -- permutation iarray: i_sorted -> i_old
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

  integer :: i, n_heap, tmp

  if (reset) then
     do i = 1, num
        new2old(i) = i
     end do
  end if
  
  call heapify()
  n_heap = num
  do while (n_heap > 1)
     call swap(1, n_heap)
     n_heap = n_heap - 1
     call downheap(1, n_heap)
  end do

contains

  ! ---

  subroutine heapify()
    ! Creates original max-heap structure

    integer :: v
    integer :: top

    top = num / 2
    do v = top, 1, -1
       call downheap(v, num)
    end do

  end subroutine heapify

  ! ---

  subroutine downheap(index, n_heap)
    integer, intent(IN) :: index, n_heap
    ! Restores the max-heap structure after extraction

    integer :: v, w

    v = index
    w = 2*v
    do while (w <= n_heap)
       if (w+1 <= n_heap) then
          if (greater(w+1, w, .false.)) then
             w = w + 1
          end if
       end if
       if (greater(v, w, .true.)) then
          return    ! index element is already larger than child elements
       else
          ! index element is too small; swap with child ...
          call swap(v ,w)
          ! ... and recurse to grandchilds if applicable
          v = w
          w = 2*v
       end if
    end do

  end subroutine downheap

  ! ---

  logical function greater(k1, k2, greater_equal)
    implicit none
    integer, intent(IN) :: k1, k2
    logical, intent(IN) :: greater_equal
    integer :: iv1(rank), iv2(rank)
    integer :: i, ii

    iv1 = iarray(:, new2old(k1))
    iv2 = iarray(:, new2old(k2))
    do i = 0, rank-1
       ii = modulo(idim+i-1, rank) + 1
       if (iv1(ii) > iv2(ii)) then
          greater = .true.; return
       else if (iv1(ii) < iv2(ii)) then
          greater = .false.; return
       end if
    end do
    greater = greater_equal
  end function greater

  ! ---

  subroutine swap(i,j)
    integer, intent(IN) :: i,j
    real*8 :: temp1
    temp1 = new2old(i)
    new2old(i) = new2old(j)
    new2old(j) = temp1
  end subroutine swap

end subroutine heapsort_index
