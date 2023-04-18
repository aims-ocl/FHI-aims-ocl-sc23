!  module geometry encapsulates all data and methods
!  to define the structure of an atomic cluster
!  
!  R. Gehrke (2005)
!

module geometry_class
  
  use cartesian_ylm
  use vector_array_class
  use cluster
  
  implicit none
  
  private
  
  include "constants.f90"
  
  public :: &
       allocate, &
       deallocate, &
       geometry, &
       get_distances, &
       get_angles, &
       write_molden_file, &
       write_aims_file, &
       write_angular_distribution_function, &
       update_energy_interval, &
       get_order_parameter_o4, &
       copy_geometry

  type geometry
     real*8 :: energy
     real*8 :: min_energy
     real*8 :: max_energy
     real*8 :: compact_ratio
     type (vector_array) :: coords
     real*8, dimension(:), pointer :: distance
     real*8, dimension(:,:), pointer :: distance_matrix
     real*8, dimension(:), pointer :: angle
     integer, dimension(:,:), pointer :: next_neighbours
     real*8, dimension(:), pointer :: spin_per_atom
     real*8, dimension(:), pointer :: charge_per_atom
     real*8 :: total_moment
     integer :: n_distances
     integer :: n_angle
     integer :: n_occurence
     integer :: histogram_out
     integer :: found_in_loop
     real*8 :: order_parameter_o4
  end type geometry

  interface assignment (=)
     module procedure copy_geometry
  end interface

  real*8, public :: diff_tolerance
  integer, public :: n_next_neighbours
  real*8, public :: r_cut_off
  real*8, public :: energy_tolerance
  
contains
  subroutine allocate(this)
    
    implicit none

    ! imported variables

    ! input
    type (geometry) :: this
    this%n_distances = n_atoms * (n_atoms - 1) / 2
    this%n_angle   = n_atoms * n_next_neighbours * (n_next_neighbours - 1) / 2
    allocate(this%distance(this%n_distances))
    allocate(this%angle(this%n_angle))
    allocate(this%next_neighbours(n_atoms, n_atoms - 1))
    allocate(this%distance_matrix(n_atoms, n_atoms))
    if (spin_polarized) then
       allocate(this%spin_per_atom(n_atoms))
    end if
    if (charged) then
       allocate(this%charge_per_atom(n_atoms))
    end if
    this%n_occurence    = 1
    this%histogram_out = 0
    this%min_energy = 0.d0
    this%max_energy = -1d9
    call allocate_array(this%coords, n_atoms)

  end subroutine allocate
  
  subroutine deallocate(this)
    
    ! imported variables
    
    ! input
    type (geometry) :: this
    
    deallocate(this%distance)
    deallocate(this%angle)
    deallocate(this%next_neighbours)
    deallocate(this%distance_matrix)
    if (spin_polarized) then
       deallocate(this%spin_per_atom)
    end if
    if (charged) then
       deallocate(this%charge_per_atom)
    end if
    call deallocate_array(this%coords)
    
  end subroutine deallocate
  
  subroutine get_distances(this)
    
    ! imported variables
    
    ! input
    type (geometry), intent(inout) :: this
    
    ! local variables
    real*8 :: max_distance
    real*8 :: min_distance
    real*8 :: temp_distance
    integer :: temp_neighbour
    integer :: max_index
    integer :: min_index
    integer :: min_counter

    ! counter
    integer :: i_atom
    integer :: i_atom_2
    integer :: i_atom_3
    integer :: i_distance
    integer :: i_distance_2
    integer :: n_distances
    integer :: i_counter
    integer :: i_counter_2
    
    ! first calculate all distances
    i_distance = 1
    do i_atom = 1, n_atoms, 1
       do i_atom_2 = i_atom + 1, n_atoms, 1
          this%distance(i_distance) = norm(this%coords%coords(i_atom) - this%coords%coords(i_atom_2))
          this%distance_matrix(i_atom  , i_atom_2) = this%distance(i_distance)
          this%distance_matrix(i_atom_2, i_atom)   = this%distance_matrix(i_atom, i_atom_2)
          i_distance = i_distance + 1
       end do
    end do
    ! then sort them
    do i_distance = 1, this%n_distances, 1
       max_distance = 0.d0
       max_index    = -1
       do i_distance_2 = i_distance, this%n_distances, 1
          if (this%distance(i_distance_2) .gt. max_distance) then
             max_distance = this%distance(i_distance_2)
             max_index    = i_distance_2
          end if 
       end do
       if (max_index .eq. -1) then
          write (6,*) "Zero coordinates!!."
          write (6,*) "* Aborting."
          stop
       end if
       temp_distance        = this%distance(max_index)
       this%distance(max_index)  = this%distance(i_distance)
       this%distance(i_distance) = temp_distance 
    end do

    this%compact_ratio = this%distance(1) / this%distance(this%n_distances)

!    write (6,'(1X, A)') "in get_distances()..."
!    do i_distance = 1, n_distances, 1
!       write (6,'(1X, A, I4, F8.3)') "distance # ", i_distance, this%distance(i_distance) * bohr
!    end do

    ! now create next_neighbours-list
    do i_atom = 1, n_atoms, 1
       ! create list for i_atom, therefore distances of distance_matrix need to be sorted
       ! first initialize next_neighbours list
       i_counter = 1
       do i_atom_2 = 1, n_atoms, 1
          if (i_atom_2 .ne. i_atom) then
             this%next_neighbours(i_atom, i_counter) = i_atom_2
             i_counter = i_counter + 1
          end if
       end do
 !      write (6,*) "initialization"
 !      write (6,*) "atom # ", i_atom
 !      do i_counter = 1, n_atoms - 1, 1
 !         write (6,*) this%next_neighbours(i_atom, i_counter)
 !      end do
       i_counter = 1
       do i_atom_2 = 1, n_atoms, 1
 !         write (6,*) "i_counter = ", i_counter, " i_atom_2 = ", i_atom_2
          if (i_atom_2 .ne. i_atom) then

             min_distance = 1e10

             i_counter_2 = i_counter
             do i_atom_3 = i_atom_2, n_atoms, 1
 !               write (6,*) "i_counter_2 = ", i_counter_2, " i_atom_3 = ", i_atom_3
                if (i_atom_3 .ne. i_atom) then
                   if (this%distance_matrix(i_atom, i_atom_3) .lt. min_distance) then
                      min_distance = this%distance_matrix(i_atom, i_atom_3)
                      min_index    = i_atom_3
                      min_counter  = i_counter_2
                   end if

                   i_counter_2 = i_counter_2 + 1
                end if
             end do
 !            write (6,*) "min_index = ", min_index, " min_counter = ", min_counter
             temp_distance = this%distance_matrix(i_atom, min_index)
             this%distance_matrix(i_atom, min_index) = this%distance_matrix(i_atom, i_atom_2)
             this%distance_matrix(i_atom, i_atom_2) = temp_distance

             temp_neighbour = this%next_neighbours(i_atom, min_counter)
             this%next_neighbours(i_atom, min_counter) = this%next_neighbours(i_atom, i_counter)
             this%next_neighbours(i_atom, i_counter) = temp_neighbour
             i_counter = i_counter + 1
          end if
       end do

!       write (6,*) "atom # ", i_atom
!       i_counter = 1
!       do i_atom_2 = 1, n_atoms, 1
       !          if (i_atom_2 .ne. i_atom) then
       
       !             write (6,*) this%next_neighbours(i_atom, i_counter), this%distance_matrix(i_atom, i_atom_2)

!             i_counter = i_counter + 1
!          end if
!       end do

    end do

  end subroutine get_distances
    
  subroutine get_angles(this)

    ! imported variables
    
    ! input
    type (geometry), intent(inout) :: this
      
    ! local variables
    real*8 :: max_angle
    real*8 :: temp_angle
    integer :: max_index
    type (vector) :: diff_a
    type (vector) :: diff_b
  
    ! counter
    integer :: i_atom
    integer :: i_atom_2
    integer :: i_atom_3
    integer :: i_angle
    integer :: i_angle_2
    integer :: n_angle
    
    ! first calculate all cosines of angles
    i_angle = 0
    do i_atom = 1, n_atoms, 1
       do i_atom_2 = 1, n_next_neighbours, 1
          do i_atom_3 = i_atom_2 + 1, n_next_neighbours, 1
             i_angle = i_angle + 1
             diff_a = this%coords%coords(this%next_neighbours(i_atom, i_atom_2)) - this%coords%coords(i_atom)
             diff_b = this%coords%coords(this%next_neighbours(i_atom, i_atom_3)) - this%coords%coords(i_atom)
             this%angle(i_angle) = (diff_a%x * diff_b%x + diff_a%y * diff_b%y + diff_a%z * diff_b%z) / &
                  (norm(diff_a) * norm(diff_b))
          end do
       end do
    end do
        
    ! then sort them
    do i_angle = 1, this%n_angle, 1
       max_angle = -1.01
       do i_angle_2 = i_angle, this%n_angle, 1
          if (this%angle(i_angle_2) .gt. max_angle) then
             max_angle = this%angle(i_angle_2)
             max_index   = i_angle_2
          end if
       end do
       temp_angle = this%angle(max_index)
       this%angle(max_index) = this%angle(i_angle)
       this%angle(i_angle) = temp_angle
    end do

  end subroutine get_angles

  subroutine get_order_parameter_o4(this, order_parameter)

    ! evaluates order parameter o4 according to Doye et al. JCP, 110(14), 6896(1999)

    ! imported variables
    
    ! input
    type (geometry), intent(inout) :: this

    ! output
    real*8, intent(out) :: order_parameter

    ! local variables
    real*8, dimension(3, n_atoms) :: dir_tab
    real*8, dimension(n_max_cartesian, 0:l_wave_max, n_atoms) :: cartesians
    real*8, dimension((l_wave_max)*(l_wave_max+1), n_atoms) :: ylm_tab
    real*8, dimension(-l_wave_max:l_wave_max) :: q_order
    integer :: lm_index
    integer :: n_bonds

    ! counter
    integer :: i_atom
    integer :: i_atom_2
    integer :: i_m

    n_bonds = 0
    q_order(:) = 0.d0
    order_parameter = 0.d0
    do i_atom = 1, n_atoms, 1
       write (6,*) "!!!!", i_atom, n_bonds
       do i_atom_2 = 1, n_atoms, 1
          dir_tab(1, i_atom_2) = this%coords%coords(i_atom)%x - this%coords%coords(i_atom_2)%x 
          dir_tab(2, i_atom_2) = this%coords%coords(i_atom)%y - this%coords%coords(i_atom_2)%y
          dir_tab(3, i_atom_2) = this%coords%coords(i_atom)%z - this%coords%coords(i_atom_2)%z
          ! exclude case i_atom eq i_atom_2 !!
          dir_tab(:, i_atom_2) = dir_tab(:, i_atom_2) / norm(this%coords%coords(i_atom) - this%coords%coords(i_atom_2))
       end do
       
       ! tabulate all ylms
       call evaluate_cartesians(dir_tab, cartesians)
       call tab_ylm_cartesian(cartesians, 4, ylm_tab)
       do i_atom_2 = i_atom + 1, n_atoms, 1
          ! check if distance lt cut off radius
 !         write (6,*) "xxx", norm(this%coords%coords(i_atom) - this%coords%coords(i_atom_2)), r_cut_off
          if (norm(this%coords%coords(i_atom) - this%coords%coords(i_atom_2)) .lt. r_cut_off) then
             n_bonds = n_bonds + 1
             do i_m = - 4, 4, 1
                lm_index = index_lm(i_m, 4)
                q_order(i_m) = q_order(i_m) + ylm_tab(lm_index, i_atom_2) 
                write (6,*) i_m, index_lm(i_m, 4), ylm_tab(lm_index, i_atom_2), q_order(i_m)
             end do
          end if
       end do

    end do
    do i_m = - 4, 4, 1
       write (6,*) i_m, q_order(i_m)
       order_parameter = order_parameter + q_order(i_m) * q_order(i_m)
    end do
    order_parameter = sqrt(order_parameter * pi4 / (2 * 4.d0 + 1)) / n_bonds
    this%order_parameter_o4 = order_parameter

  end subroutine get_order_parameter_o4

  subroutine write_molden_file(this, index)
      
    include "constants.f90"
      
    ! imported variables

    ! input
    type (geometry), intent(in) :: this
    integer :: index

    ! counters
    integer :: i_atom

    ! local variables
    character(len = 17) :: filename
    character(len =  4) :: number
    
    write(number, '(I4.4)') index
      
    number = adjustl(number)
  
    filename = "structure" // trim(number) // ".xyz"
  
    open (50, file = filename)

    write (50, '(2X, I4)') n_atoms
    write (50, *)
  
    do i_atom = 1, n_atoms, 1

       write (50,'(A, F10.6, 1X, F10.6, 1X, F10.6)') &
            species_name(species(i_atom)), &
            this%coords%coords(i_atom)%x, this%coords%coords(i_atom)%y, &
            this%coords%coords(i_atom)%z
       
    end do
      
    close(50)
      
  end subroutine write_molden_file

  subroutine write_aims_file(this, index, global_minimum, ref_energy)
      
    include "constants.f90"
      
    ! imported variables

    ! input
    type (geometry), intent(in) :: this
    integer, intent(in) :: index
    real*8, intent(in) :: global_minimum
    real*8, intent(in) :: ref_energy

    ! counters
    integer :: i_atom

    ! local variables
    character(len = 16) :: filename
    character(len =  4) :: number
    real*8 :: total_moment
    
    write(number, '(I4.4)') index
      
    number = adjustl(number)
  
    filename = "geometry.in." // trim(number)
  
    open (50, file = filename)

    write (50, '(A)') "# geometry-file created by effernir"
    write (50, '(A, F14.6)') "# total energy / N [eV] ", this%energy
    write (50, '(A, F10.6)') "# difference to assumed global minimum [eV] ", this%energy - global_minimum
    write (50, '(A, F10.6)') "# dissociation energy [eV] ", this%energy  - ref_energy
    write (50, '(A, F10.6)') "# order parameter o4 ", this%order_parameter_o4
    write (50, '(A, I4)') "# found in loop # ", this%found_in_loop

    if (spin_polarized) then
       total_moment = 0.d0
       do i_atom = 1, n_atoms, 1
          total_moment = total_moment + this%spin_per_atom(i_atom)
          write (50, '(A, I4, 1X, F10.6)') "# spin on at ", i_atom, this%spin_per_atom(i_atom)
       end do
       write (50, '(A, F10.6)') "# total moment N_up - N_down ", total_moment 
    end if

    do i_atom = 1, n_atoms, 1

       write (50,'(A, F10.6, 1X, F10.6, 1X, F10.6, 1X, A)') &
            "atom  ", this%coords%coords(i_atom)%x, this%coords%coords(i_atom)%y, &
            this%coords%coords(i_atom)%z, species_name(species(i_atom))
       
    end do
      
    close(50)
      
  end subroutine write_aims_file

  subroutine write_angular_distribution_function(this, index)
    
    include "constants.f90"
    
    ! imported variables
    
    ! input
    type (geometry), intent(in) :: this
    integer :: index

    ! counters
    integer :: i_angle

    ! local variables
    character(len = 13) :: filename
    character(len =  4) :: number
    real*8 :: delta_cosine = - 1e-2
    real*8 :: cosine
    real*8 :: value
    
    write(number, '(I4)') index
      
    number = adjustl(number)
  
    filename = "angle" // trim(number) // ".dat"
  
    open (50, file = filename)
    
    i_angle = 1
    do cosine= 1.0, - 1.0, delta_cosine
       value = 0.d0
       do while ((this%angle(i_angle) .gt. cosine + delta_cosine) .and. (i_angle .le. this%n_angle)) 
          value = value + 1.d0
          i_angle = i_angle + 1
       end do
       value = value / this%n_angle
       write (50, '(F10.4, 1X, F10.4)') cosine, value
    end do
      
    close(50)
 
  end subroutine write_angular_distribution_function

  subroutine update_energy_interval(this, energy)

! imported variables

! input
    type (geometry), intent(inout) :: this
    real*8, intent(in) :: energy

    if (energy .lt. this%min_energy) then
       this%min_energy = energy
    end if
    if (energy .gt. this%max_energy) then
       this%max_energy = energy
    end if

  end subroutine update_energy_interval

  subroutine copy_geometry(a,b)

    implicit none
    
    ! imported variables

    ! input
    type (geometry), intent(out) :: a
    type (geometry), intent(in) :: b

    ! counter
    integer :: i_distance
    integer :: i_angle
    integer :: i_next_neighbour
    integer :: i_atom
    integer :: i_atom_2

    a%energy = b%energy
    a%min_energy = b%min_energy
    a%max_energy = b%max_energy
    a%coords = b%coords
    a%n_distances = b%n_distances
    a%n_angle = b%n_angle
    a%n_occurence = b%n_occurence
    a%order_parameter_o4 = b%order_parameter_o4
    
    if (a%n_distances .ne. b%n_distances) then
       write (6,*) "Geometries to be copied have different number of distances!!"
       stop
    end if

    if (a%n_angle .ne. b%n_angle) then
       write (6,*) "Geometries to be copied have different number of angles!!"
       stop
    end if

    do i_distance = 1, a%n_distances, 1
       a%distance(i_distance) = b%distance(i_distance) 
    end do

    a%found_in_loop = b%found_in_loop

    do i_atom = 1, n_atoms, 1

       if (spin_polarized) then
          a%spin_per_atom(i_atom) = b%spin_per_atom(i_atom)
       end if

       do i_atom_2 = 1, n_atoms, 1
          a%distance_matrix(i_atom, i_atom_2) = b%distance_matrix(i_atom, i_atom_2)
       end do
    end do
    if (spin_polarized) then
       a%total_moment = b%total_moment
    end if
    do i_atom = 1, n_atoms, 1
       do i_atom_2 = 1, n_atoms - 1, 1
          a%next_neighbours(i_atom, i_atom_2) = b%next_neighbours(i_atom, i_atom_2)
       end do
    end do

    do i_angle = 1, a%n_angle, 1
       a%angle(i_angle) = b%angle(i_angle) 
    end do

  end subroutine copy_geometry

end module geometry_class
  
