!  subroutine to search the energetically ordered list of minima for a 
!  candidate to add; if candidate does not exist yet, it is add 
!  at the corresponding position
!
!  R. Gehrke (2005)

subroutine check_and_add_list_of_minima(list_of_geometries, new_minimum, found)

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
  logical, intent(out) :: found

  ! local
  type (node_of_geometry), pointer :: list 
  type (node_of_geometry), pointer :: previous_node 

  ! counter
  integer :: i_distance

  found = .false.
  list => list_of_geometries%head

  nullify(previous_node)
  do while (associated(list))
     call compare_distances(list%inst_geometry, new_minimum%inst_geometry, diff_tolerance, found)
!     write (*,*) "found?", found

!  structure equivalent? if yes, exit...
     if (found) then
        list%inst_geometry%n_occurence = list%inst_geometry%n_occurence + 1
        call update_energy_interval(list%inst_geometry, new_minimum%inst_geometry%energy)
        return
     end if
!  ...if no, check, whether energy of minimum to be checked is lower
!  if yes, add new minimum to linked list...

!     write (*,*) "check energies:", "new energy: ", new_minimum%inst_geometry%energy , &
!          "energy in current_node", list%inst_geometry%energy

     if (new_minimum%inst_geometry%energy .lt. list%inst_geometry%energy) then
        new_minimum%next => list
!  is new minimum energetically the lowest one? ...
        if (associated(previous_node)) then
           previous_node%next => new_minimum
        else
!  yes => add minimum to the beginning of list
           list_of_geometries%head => new_minimum
        end if
        call update_energy_interval(new_minimum%inst_geometry, new_minimum%inst_geometry%energy)
        return
     end if
!  ...if no, go ahead
     previous_node => list
     list          => list%next

  end do
!  ...energy of new minimum is the highest at the moment
!  so add to the end of list
  previous_node%next => new_minimum
  call update_energy_interval(new_minimum%inst_geometry, new_minimum%inst_geometry%energy)

end subroutine check_and_add_list_of_minima
