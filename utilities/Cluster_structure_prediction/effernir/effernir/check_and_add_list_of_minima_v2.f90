!  subroutine to search the energetically ordered list of minima for a 
!  candidate to add; if candidate does not exist yet, it is add 
!  at the corresponding position
!
!  R. Gehrke (2005)

subroutine check_and_add_list_of_minima_v2(list_of_geometries, new_minimum, found)

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

  ! first run to look if the minimum has already been found
  found = .false.
  list => list_of_geometries%head

  do while (associated(list))
     call compare_distances(list%inst_geometry, new_minimum%inst_geometry, diff_tolerance, found)
     
     if (found) then
        ! are energies also sufficiently similar?
        found = .false.
        if (abs(list%inst_geometry%energy - new_minimum%inst_geometry%energy) .lt. energy_tolerance) then
           found = .true.
        end if
     end if

     !  structure equivalent? if yes, exit...
     if (found) then
        list%inst_geometry%n_occurence = list%inst_geometry%n_occurence + 1
        call update_energy_interval(list%inst_geometry, new_minimum%inst_geometry%energy)
        return
     end if
     list => list%next
  end do

  if (.not.found) then

     ! minimum hasn't been found yet so add it
     ! according to energetic position

     list => list_of_geometries%head

     nullify(previous_node)
     do while (associated(list))

        ! check, whether energy of minimum to be checked is lower
        ! if yes, add new minimum to linked list...
        if (new_minimum%inst_geometry%energy .lt. list%inst_geometry%energy) then
           new_minimum%next => list
           ! is new minimum energetically the lowest one? ...
           if (associated(previous_node)) then
              previous_node%next => new_minimum
           else
              ! yes => add minimum to the beginning of list
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
     
  else
     write (6,*) "Internal inconsistency in check_and_add_list_of_minima."
     write (6,*) "* Aborting."
     stop
  end if

end subroutine check_and_add_list_of_minima_v2
