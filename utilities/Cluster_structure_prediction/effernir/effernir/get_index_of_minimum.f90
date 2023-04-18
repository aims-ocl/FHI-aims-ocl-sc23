!  subroutine to search the energetically ordered list of minima for a 
!  candidate to add; if candidate does not exist yet, it is add 
!  at the corresponding position
!
!  R. Gehrke (2005)

subroutine get_index_of_minimum(list_of_geometries, new_minimum, index)

  use geometry_class
  use geometry_list
  use cluster

  implicit none

  include "constants.f90"

  ! imported variables

  ! input
  type (list_of_geometry) :: list_of_geometries
  type (node_of_geometry), target :: new_minimum

  ! output
  integer, intent(out) :: index

  ! local
  type (node_of_geometry), pointer :: list 
  logical :: found

  ! counter
  integer :: i_distance
  integer :: i_counter

  index = -1
  i_counter = 1

  list => list_of_geometries%head
  do while (associated(list))
     call compare_distances(list%inst_geometry, new_minimum%inst_geometry, diff_tolerance, found)
     if (found) then
        index = i_counter
        return
     end if
     list => list%next
     i_counter = i_counter + 1
  end do

end subroutine get_index_of_minimum
