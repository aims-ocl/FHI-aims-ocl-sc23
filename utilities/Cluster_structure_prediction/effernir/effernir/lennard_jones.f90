!  calculation of total energy of a cluster 
!  using a simple lennard-jones potential
!
!  R.Gehrke(2005)


module lennard_jones

use cluster
use vector_array_class

real*8, dimension(:,:), allocatable :: sigma
real*8, dimension(:,:), allocatable :: epsilon

contains
  
  subroutine initialize_lennard_jones()
    
    write (*,*) "Allocate memory for lennard jones data"
    if ((n_atoms .gt. 0) .and. (n_species .gt. 0)) then
       if (.not.allocated(sigma)) then
          allocate(sigma(n_species, n_species))
       end if
       if (.not.allocated(epsilon)) then
          allocate(epsilon(n_species, n_species))
       end if
    else
       write (*,*) "No atoms and species found yet."
       write(*,*) "* Aborting."
       stop
    end if

    !  initialize lennard-jones parameters due to Lorentz-Berthelot's mixing rules
    do i_species = 1, n_species, 1
       do i_species_2 = i_species, n_species, 1
          
          sigma(i_species, i_species_2) = &
               0.5 * (LJ_sigma(i_species) + LJ_sigma(i_species_2))
          sigma(i_species_2, i_species) = sigma(i_species, i_species_2)
          epsilon(i_species, i_species_2) = &
               sqrt(LJ_epsilon(i_species) * LJ_epsilon(i_species_2))
          epsilon(i_species_2, i_species) = epsilon(i_species, i_species_2)
          write (*,*) "sigma/epsilon = ", sigma(i_species, i_species_2), " ", epsilon(i_species, i_species_2)
       end do
    end do

  end subroutine initialize_lennard_jones

  subroutine deallocate_lennard_jones()

    implicit none

    if (allocated(sigma)) then
       deallocate(sigma)
    end if
    if (allocated(epsilon)) then
       deallocate(epsilon)
    end if
  
  end subroutine deallocate_lennard_jones

  real*8 function get_lennard_jones_energy(positions)

    implicit none

    include "constants.f90"
    
    !  imported variables

    !  input
    type (vector_array) :: positions

    !  local variables
    
    real*8 :: distance
    real*8 :: energy
    
    !  counters
    
    integer :: i_atom
    integer :: i_atom_2
    
    energy = 0.d0
    do i_atom = 1, n_atoms, 1
       do i_atom_2 = i_atom + 1, n_atoms, 1
           distance = norm(positions%coords(i_atom) - positions%coords(i_atom_2))
           energy = energy + 4 * epsilon(species(i_atom), species(i_atom_2)) * &
                ((sigma(species(i_atom), species(i_atom_2)) / distance) ** 12 - &
                (sigma(species(i_atom), species(i_atom_2)) / distance) ** 6)
        end do
     end do
    
    get_lennard_jones_energy = energy
    
  end function get_lennard_jones_energy
  
  real*8 function get_single_lj_energy(distance, i_atom, i_atom_2)

    implicit none

    !  imported variables

    !  input
    real*8 :: distance
    integer :: i_atom
    integer :: i_atom_2

    get_single_lj_energy = 4 * epsilon(species(i_atom), species(i_atom_2)) * &
         ((sigma(species(i_atom), species(i_atom_2)) / distance) ** 12 - &
         (sigma(species(i_atom), species(i_atom_2)) / distance) ** 6)
    
    
  end function get_single_lj_energy

  subroutine get_lj_forces_and_energy(positions, forces, energy)

    implicit none
    
    !  imported variables

    !  input
    type (vector_array), intent(in) :: positions

    !  output
    type (vector_array), intent(out) :: forces
    real*8, intent(out) :: energy

    !  local variables
    
    real*8 :: distance
    real*8 :: one_over_d_sqr
    real*8 :: repulsive
    real*8 :: attractive
    
    !  counters
    
    integer :: i_atom
    integer :: i_atom_2
    
    energy = 0.d0

    do i_atom = 1, n_atoms, 1

       forces%coords(i_atom)%x = 0.d0
       forces%coords(i_atom)%y = 0.d0
       forces%coords(i_atom)%z = 0.d0
       
       do i_atom_2 = 1, n_atoms, 1
          
          if (i_atom_2 .ne. i_atom) then
             distance = norm(positions%coords(i_atom) - positions%coords(i_atom_2))
             one_over_d_sqr = 1. / (distance * distance)
             repulsive  = (sigma(species(i_atom), species(i_atom_2)) / distance) ** 12
             attractive = (sigma(species(i_atom), species(i_atom_2)) / distance) ** 6

! divide by factor 2 to avoid double counting => 2 * epsilon instead of 4 * epsilon
             energy = energy + 2 * epsilon(species(i_atom), species(i_atom_2)) * &
                  (repulsive - attractive)

             forces%coords(i_atom)%x = forces%coords(i_atom)%x + epsilon(species(i_atom), species(i_atom_2)) * &
                  one_over_d_sqr * (48 * repulsive - 24 * attractive) * &
                  (positions%coords(i_atom)%x - positions%coords(i_atom_2)%x)
             forces%coords(i_atom)%y = forces%coords(i_atom)%y + epsilon(species(i_atom), species(i_atom_2)) * &
                  one_over_d_sqr * (48 * repulsive - 24 * attractive) * &
                  (positions%coords(i_atom)%y - positions%coords(i_atom_2)%y)
             forces%coords(i_atom)%z = forces%coords(i_atom)%z + epsilon(species(i_atom), species(i_atom_2)) * &
                  one_over_d_sqr * (48 * repulsive - 24 * attractive) * &
                  (positions%coords(i_atom)%z - positions%coords(i_atom_2)%z)
          end if

       end do
    end do

  end subroutine get_lj_forces_and_energy

end module lennard_jones
