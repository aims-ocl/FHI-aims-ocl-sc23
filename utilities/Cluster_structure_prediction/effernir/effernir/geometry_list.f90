!  module geometry_list encapsulates a 
!  linked list of geometry objects
!
!  R. Gehrke (2005)
!

module geometry_list

  use geometry_class, g_allocate => allocate, g_deallocate => deallocate
  use control
  use cluster

  implicit none

  private

  public :: &
       list_of_geometry, &
       node_of_geometry, &
       write_energy_of_structures, &
       write_geometry_of_structures, &
       write_aims_files_of_structures, &
       get_number_of_minima, &
       analyze_efficiency_of_sampling, &
       compare_lowest_energy, &
       write_angular_files_of_structures, &
       get_coords, &
       increase_histogram_out
 
  include "constants.f90"

  type list_of_geometry
     type (node_of_geometry), pointer :: head
  end type list_of_geometry

  type node_of_geometry
     type (geometry) :: inst_geometry
     type (node_of_geometry), pointer :: next
  end type node_of_geometry
  
  contains

    subroutine write_energy_of_structures(this)
  
      ! imported variables

      ! input
      type (list_of_geometry) :: this

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter
      integer :: checksum
      real*8 :: global_minimum
      real*8 :: ref_energy
      real*8 :: average_compact_ratio

      i_counter = 1
      checksum = 0
      average_compact_ratio = 0.d0
      list => this%head
      global_minimum = list%inst_geometry%energy
      ref_energy = reference_energy()
      if (spin_polarized) then
         write (150, '(A,A)') "energy of minimum [eV] (total energy, diff to E_0, dissociation energy), total moment, hist_in", &
              ", hist_out, checksum, deviation, min, max, compact, found in loop"
      else
         write (150, '(A,A)') "energy of minimum [eV] (total energy, diff to E_0, dissociation energy), hist_in", &
              ", hist_out, checksum, deviation, min, max, compact, found in loop"
      end if
      do while (associated(list))
         checksum = checksum + list%inst_geometry%n_occurence
         average_compact_ratio = average_compact_ratio + list%inst_geometry%n_occurence * list%inst_geometry%compact_ratio
         if (spin_polarized) then
            write (150, '(A, I4, 1X, F14.6, 1X, F10.6, 1X, F10.6, 1X, F8.3, 1X, I6, 1X, I6, 1X, I6, 1X, F10.6, 1X, F14.6, 1X, F14.6, 1X, F8.4, 1X, I4)') &
                 "# ", i_counter, list%inst_geometry%energy, &
                 list%inst_geometry%energy - global_minimum, list%inst_geometry%energy - ref_energy, &
                 list%inst_geometry%total_moment, list%inst_geometry%n_occurence, list%inst_geometry%histogram_out, checksum, &
                 (list%inst_geometry%max_energy - list%inst_geometry%min_energy) , &
                 list%inst_geometry%min_energy, list%inst_geometry%max_energy, list%inst_geometry%compact_ratio, &
                 list%inst_geometry%found_in_loop
         else
            write (150, '(A, I4, 1X, F14.6, 1X, F10.6, 1X, F10.6, 1X, I6, 1X, I6, 1X, I6, 1X, F10.6, 1X, F14.6, 1X, F14.6, 1X, F8.4, 1X, I4)') &
                 "# ", i_counter, list%inst_geometry%energy, &
                 list%inst_geometry%energy - global_minimum, list%inst_geometry%energy - ref_energy, &
                 list%inst_geometry%n_occurence, list%inst_geometry%histogram_out, checksum, &
                 (list%inst_geometry%max_energy - list%inst_geometry%min_energy) , &
                 list%inst_geometry%min_energy, list%inst_geometry%max_energy, list%inst_geometry%compact_ratio, &
                 list%inst_geometry%found_in_loop
         end if
         i_counter = i_counter + 1
         list => list%next
      end do
      write (150,'(1X, A, I6)') "number of all occured minima: ", checksum
      write (150,'(1X, A, F8.4)') "average compact ratio: ", average_compact_ratio / checksum

    end subroutine write_energy_of_structures

    subroutine write_geometry_of_structures(this)

      ! imported variables

      ! input
      type (list_of_geometry) :: this

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter
      
      i_counter = 1
      list => this%head
      do while (associated(list))
         call write_molden_file(list%inst_geometry, i_counter)
         i_counter = i_counter + 1
         list => list%next
      end do
      
    end subroutine write_geometry_of_structures

    subroutine write_aims_files_of_structures(this)

      ! imported variables

      ! input
      type (list_of_geometry) :: this

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter
      real*8 :: global_minimum
      real*8 :: ref_energy

      i_counter = 1
      list => this%head
      global_minimum = list%inst_geometry%energy
      ref_energy = reference_energy()
      do while (associated(list))
         call write_aims_file(list%inst_geometry, i_counter, global_minimum, ref_energy)
         i_counter = i_counter + 1
         list => list%next
      end do
      
    end subroutine write_aims_files_of_structures

    subroutine write_angular_files_of_structures(this)

      ! imported variables
      
      ! input
      type (list_of_geometry) :: this
      
      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter

      i_counter = 1
      list => this%head
      do while (associated(list))
         call write_angular_distribution_function(list%inst_geometry, i_counter)
         i_counter = i_counter + 1
         list => list%next
      end do

    end subroutine write_angular_files_of_structures

    subroutine get_number_of_minima(this, delta_energy, number_of_minima)
    
      ! imported variables
      
      ! input
      type (list_of_geometry), intent(in) :: this
      real*8, intent(in) :: delta_energy
      
      ! output
      integer, intent(out) :: number_of_minima

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter
      real*8 :: min_energy
      
      i_counter = 0
      list => this%head
      min_energy = list%inst_geometry%energy
      do while (associated(list) .and. (list%inst_geometry%energy .lt. (min_energy + delta_energy/n_atoms)))
         i_counter = i_counter + 1
         list => list%next
      end do
      number_of_minima = i_counter

    end subroutine get_number_of_minima

    subroutine analyze_efficiency_of_sampling(this, number_of_minima, mean_value, standard_deviation)

      ! imported variables
      
      ! input
      type (list_of_geometry), intent(in) :: this
      
      ! output
      integer, intent(out) :: number_of_minima
      real*8, intent(out) :: mean_value
      real*8, intent(out) :: standard_deviation

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter
      real*8 :: energy_min

      i_counter = 0
      mean_value = 0.d0
      list => this%head
      energy_min = this%head%inst_geometry%energy
      do while (associated(list))
         if ((list%inst_geometry%energy - energy_min)*n_atoms .lt. energy_interval) then
            i_counter = i_counter + 1
            mean_value = mean_value + list%inst_geometry%n_occurence
         end if
         list => list%next
      end do
      number_of_minima = i_counter
      mean_value = mean_value / number_of_minima
      
      ! run a second time through linked list to calculate the standard_deviation of
      ! histogramm entries
      standard_deviation = 0.d0
      list => this%head
      do while (associated(list))
         if ((list%inst_geometry%energy - energy_min)*n_atoms .lt. energy_interval) then
            i_counter = i_counter + 1
            standard_deviation = standard_deviation + (list%inst_geometry%n_occurence - mean_value) * (list%inst_geometry%n_occurence - mean_value)
         end if
         list => list%next
      end do
      if (number_of_minima .gt. 1) then
         standard_deviation = sqrt(standard_deviation / (number_of_minima - 1))
      else
         standard_deviation = -1
      end if

    end subroutine analyze_efficiency_of_sampling

    subroutine compare_lowest_energy(this, target_energy, optimized)

      !  imported variables
      
      !  input
      type (list_of_geometry), intent(in) :: this
      real*8, intent(in) :: target_energy
      
      !  output
      logical, intent(out) :: optimized

      ! local variables
      real*8 :: energy
      real*8 :: ref_energy

      ref_energy = reference_energy()
      energy = this%head%inst_geometry%energy - ref_energy

      if (abs(energy - target_energy) .lt. 1e-2) then
         optimized = .true.
      else
         optimized = .false.
      endif
            
    end subroutine compare_lowest_energy

    subroutine get_coords(this, index, coords, start_sample)

      !  imported variables
      
      !  input
      type (list_of_geometry), intent(in) :: this
      integer, intent(in) :: index
      
      ! output
      type (vector_array) :: coords
      integer, intent(out) :: start_sample

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter
      
      i_counter = 0
      list => this%head
      do while (associated(list))
         i_counter = i_counter + 1
!         write (150,*) i_counter, index
!         call display_array(list%inst_geometry%coords, 150)
         if (i_counter .eq. index) then
            coords = list%inst_geometry%coords
            start_sample = list%inst_geometry%n_occurence
         end if
         list => list%next
      end do

    end subroutine get_coords

    subroutine increase_histogram_out(this, index)

      !  imported variables
      
      !  input
      type (list_of_geometry), intent(in) :: this
      integer, intent(in) :: index

      ! local variables
      type (node_of_geometry), pointer :: list
      integer :: i_counter

      list => this%head
      do i_counter = 2, index, 1
         list => list%next
      end do
      list%inst_geometry%histogram_out = list%inst_geometry%histogram_out + 1

    end subroutine increase_histogram_out

 end module geometry_list
