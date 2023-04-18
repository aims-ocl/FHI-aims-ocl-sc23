!  Module basin-hopping encapsulates all necessary procedures
!  for a bh algorithm on clusters
!    
!  Subroutines and functions:
!
!  R.Gehrke (2005)
!

module basin_hopping
  
  use vector_array_class
  use cluster
  use control
  use geometry_class, g_allocate => allocate, g_deallocate => deallocate
  use geometry_list
  use dft_file
  use arch_specific
  use conjugate_gradient, cg_optimize => optimize, cg_allocate => allocate, cg_deallocate => deallocate
  use lennard_jones
  use pair_potential, pp_allocate => allocate, pp_deallocate => deallocate, pp_initialize => initialize

  implicit none

  private
  
  include "constants.f90"

  type (vector_array), public :: current_coords
  type (vector_array) :: try_coords
  type (vector_array) :: optimized_coords
  real*8, allocatable, dimension(:) :: try_spin_per_atom 
  real*8, allocatable, dimension(:) :: current_spin_per_atom 
  real*8, allocatable, dimension(:) :: optimized_spin_per_atom
  integer, public :: max_loop
  real*8, public  :: energy_threshold_theta
  real*8, public  :: energy_threshold_boltzmann
  logical, public :: reset
  real*8, public :: cavity_radius
  logical, public :: cavity
  character*30, public :: minimizer
  type (list_of_geometry) :: list_of_geometries
  real*8, public :: max_move
  character*30, public :: kind_of_move
  real*8, public :: max_move_increase
  real*8, public :: max_move_decrease
  logical, public :: move_adaptive
  character*30, public :: move_distribution
  real*8, public :: max_distance
  real*8, public :: target_energy
  logical, public :: optimize_energy
  integer, public :: target_number_of_minima
  integer, allocatable, dimension(:), public :: adsorption_atom
  logical, public :: adsorption_move
  logical, public :: pure_adsorption_move
  real*8, allocatable, dimension(:), public :: ad_radius_reduction
  logical, public :: use_rotation_move
  logical, public :: use_angular_move
  integer, public :: fraction_of_rotation_move
  real*8, public :: energy_ratio_ang_move
  real*8, public :: max_successful_rate
  real*8, public :: min_successful_rate
  character*30, public :: adapting_scheme
  character*30, public :: acceptance_reference
  integer, public :: n_adsorption_atom
  integer, public :: scf_loop_to_switch
  real*8, public :: switched_E_theta
  real*8, public :: switched_E_boltzmann
  real*8, public :: spin_step_width
  character*30, public :: bh_flavour
  integer, public :: n_starting_isomers
  integer, allocatable, dimension(:), public :: starting_index
  integer, public :: n_samples
 
  ! parameter for maximum rock adaptive scheme
  integer, public :: number_of_averaging
  integer, public :: number_of_coarse_averaging
  real*8, public  :: max_delta_max_move
  real*8, public  :: min_delta_max_move
  real*8, public  :: delta_max_move
  real*8, public  :: max_max_move
  real*8, public  :: kick_in_the_ass
  real*8, public  :: max_move_modifier_one
  real*8, public  :: max_move_modifier_two
  real*8, public  :: alpha_unsucc_max
  real*8, public :: alpha_highE_max
  real*8, public :: alpha_tolerance

  public :: allocate, deallocate, optimize, node_of_geometry, rebuild_linked_list, analyze_moves, sample_hop_matrix

  contains

    subroutine allocate()
      
      implicit none
      
      call allocate_array(optimized_coords, n_atoms)
      call allocate_array(current_coords, n_atoms)
      call allocate_array(try_coords, n_atoms)
      
      if (potential_flag .eq. 'LJ') then
         call initialize_lennard_jones()
      end if
      if (minimizer .eq. 'cg') then
         call cg_allocate()
      end if
      if (use_pair_pot) then
         call pp_initialize()
      end if
      if (spin_polarized) then
         allocate(try_spin_per_atom(n_atoms))
         allocate(current_spin_per_atom(n_atoms))
         allocate(optimized_spin_per_atom(n_atoms))
      end if

    end subroutine allocate

    subroutine deallocate()
    
      call deallocate_array(optimized_coords)
      call deallocate_array(current_coords)
      call deallocate_array(try_coords)
      if (potential_flag .eq. 'LJ') then
         call deallocate_lennard_jones()
      end if
      if (minimizer .eq. 'cg') then
         call cg_deallocate()
      end if
      if (use_pair_pot) then
         call pp_deallocate()
      end if
      if (allocated(adsorption_atom)) then
         deallocate(adsorption_atom)
      end if
      if (allocated(ad_radius_reduction)) then
         deallocate(ad_radius_reduction)
      end if
      if (allocated(try_spin_per_atom)) then
         deallocate(try_spin_per_atom)
      end if
      if (allocated(current_spin_per_atom)) then
         deallocate(current_spin_per_atom)
      end if
      if (allocated(optimized_spin_per_atom)) then
         deallocate(optimized_spin_per_atom)
      end if

    end subroutine deallocate

    subroutine initialize_positions()

      implicit none

      real*8 :: dummy
      character*10 :: arg1, arg2, arg3
      integer, dimension(8) :: values

      integer :: i_counter
      integer :: random_number
      integer :: i_atom

      current_coords = initial_coords
      try_coords     = current_coords

      if (spin_polarized) then
         do i_atom = 1, n_atoms, 1
            current_spin_per_atom(i_atom) = ini_spin_per_atom(i_atom)
            try_spin_per_atom (i_atom) = current_spin_per_atom(i_atom)
         end do
      end if

      ! initialize random number generator with milliseconds of real-time clock 

      call date_and_time(arg1, arg2, arg3, values)
    
      if (potential_flag .ne. 'external') then
         ! initialize random number generator with milliseconds of real-time clock 
         call date_and_time(arg1, arg2, arg3, values)
         seed_uni = values(8)
         call date_and_time(arg1, arg2, arg3, values)
         seed_poiss = values(7)
      end if

      ! initialize uniform RNG
      call srand(seed_uni)

      ! initialize poisson-distributed RNG
      ! FIXME: This one is very crappy!!!
      ! FIXME: Search for a better one
      call zufalli(seed_poiss)

      write (6,*) "uniform RNG test:"
      write (6,*) "seed=", seed_uni
      do i_counter = 1, 10, 1
         write (6,*) arch_rand()
      end do
      write (6,*) "poisson RNG test:"
      write (6,*) "seed=", seed_poiss
      do i_counter = 1, 10, 1
         call fische(1, 100d0, random_number)
         write (6,*) random_number
      end do
      
    end subroutine initialize_positions

    subroutine trial_move(square_of_move)
      
      implicit none

      ! imported variables
      
      ! output
      real*8, intent(out) :: square_of_move

      ! local variables
      type (vector) :: move
      integer :: i_atom
      logical :: reasonable
      logical :: outside
      logical :: dissociated
      real*8 :: theta
      real*8 :: phi
      real*8 :: radius
      integer :: random_number
      real*8 :: pair_energy_max

      integer, save :: average_random_atom = 0
      real*8, save :: average_move_x = 0
      real*8, save :: average_move_y = 0
      real*8, save :: average_move_z = 0
      integer, save :: move_counter = 0

      logical :: found

      ! counter
      integer :: i_adsorp

      reasonable = .false.
      write (150,*) "generating single-particle move..."
      do while (.not. reasonable)
         
         try_coords = current_coords

         ! choose the atom
         select case (kind_of_move)
            
         case ('worst-single')
            call get_highest_pair_energy(current_coords, pair_energy_max, i_atom)
            
         case default

            i_atom = int(n_atoms * arch_rand()) + 1

            found = .false.
            do i_adsorp = 1, n_adsorption_atom, 1
               if (i_atom .eq. adsorption_atom(i_adsorp)) then
                  found = .true.
               end if
            end do

            do while (adsorption_move .and. found)
               i_atom = int(n_atoms * arch_rand()) + 1

               found = .false.
               do i_adsorp = 1, n_adsorption_atom, 1
                  if (i_atom .eq. adsorption_atom(i_adsorp)) then
                     found = .true.
                  end if
               end do
               
            end do

         end select

         average_random_atom = average_random_atom + i_atom
         move_counter = move_counter + 1

         ! get spherical coordinates
         theta = pi * arch_rand()
         phi   = 2 * pi * arch_rand()

         select case (move_distribution)

         case ('uniform')
            radius = max_move * arch_rand()

         case ('poisson')
            call fische(1, 100d0 * max_move, random_number)
            radius = dble(random_number) / 100d0
    
         case default
            write (6,*) "Internal incosistency."
            write (6,*) "(trial_move)"
            write (6,*) "* Aborting."
            stop

         end select

         if (verbose) then
            write (150,*) "radius=", radius
            write (150,*) "theta=", theta
            write (150,*) "phi=", phi
         end if

         move%x = radius * sin(theta) * cos(phi)
         move%y = radius * sin(theta) * sin(phi)
         move%z = radius * cos(theta)
         
         average_move_x = average_move_x + move%x
         average_move_y = average_move_y + move%y
         average_move_z = average_move_z + move%z

         square_of_move = square_norm(move)
      
         try_coords%coords(i_atom) = &
              try_coords%coords(i_atom) + move
         
         call check_hard_sphere(try_coords, n_atoms, reasonable)
         if (.not.reasonable .and. verbose) then
            write (150,*) "Hard sphere criterium violated!!"
         end if

         if (cavity) then
            call check_cavity(try_coords, cavity_radius, outside)
            if (outside .and. verbose) then
               write (150,*) "Trial move put some atoms outside!!"
            end if
            reasonable = (reasonable .and. .not.(outside))
         end if
         call check_compound(try_coords, max_distance, dissociated)
         if (dissociated .and. verbose) then
            write (150,*) "Cluster structure dissociated from trial move!!"
         end if
         reasonable = (reasonable .and. .not.(dissociated))

         if (verbose) then
            write (150,*) "atom to move", i_atom
            write (150,*) "trial move:"
            call display(move, 150)
         end if
      end do
      
!      write (150,*) "average atom" , dble(average_random_atom) / dble(move_counter)
!      write (150,*) "average move x" , average_move_x / dble(move_counter)
!      write (150,*) "average move y" , average_move_y / dble(move_counter)
!      write (150,*) "average move z" , average_move_z / dble(move_counter)

    end subroutine trial_move

    subroutine trial_move_all()
      
      implicit none

      ! imported variables
      
      ! local variables
      type (vector) :: move
      logical :: reasonable
      logical :: outside
      logical :: dissociated
      real*8 :: theta
      real*8 :: phi
      real*8 :: radius
      integer :: random_number
      logical :: found

      ! functions
      real*4 :: rand

      ! counter
      integer :: i_atom
      integer :: i_adsorp

      reasonable = .false.
      write (150,*) "generating collective move..."
      do while (.not. reasonable)
         
         try_coords = current_coords
         if (verbose) then
            write (150,*) "trial move"
         end if
         do i_atom = 1, n_atoms, 1

            found = .false.
            do i_adsorp = 1, n_adsorption_atom, 1
               if (i_atom .eq. adsorption_atom(i_adsorp)) then
                  found = .true.
               end if
            end do

            if (.not.(adsorption_move .and. found)) then
!            move%x = 2 * max_move * (arch_rand() - 0.5d0)
!            move%y = 2 * max_move * (arch_rand() - 0.5d0)
!            move%z = 2 * max_move * (arch_rand() - 0.5d0)
            ! get spherical coordinates
               theta = pi * arch_rand()
               phi   = 2 * pi * arch_rand()
            
               select case (move_distribution)
               
               case ('uniform')
                  radius = max_move * arch_rand()
               
               case ('poisson')
                  call fische(1, 100d0 * max_move, random_number)
                  radius = dble(random_number) / 100d0
                  
               case default
                  write (6,*) "Internal incosistency."
                  write (6,*) "(trial_move)"
                  write (6,*) "* Aborting."
                  stop
               
               end select
            
               if (verbose) then
                  write (150,*) "radius=", radius
                  write (150,*) "theta=", theta
                  write (150,*) "phi=", phi
               end if
               move%x = radius * sin(theta) * cos(phi)
               move%y = radius * sin(theta) * sin(phi)
               move%z = radius * cos(theta)
               
               try_coords%coords(i_atom) = &
                    try_coords%coords(i_atom) + move

               if (verbose) then
                  call display(move, 150)
               end if

            end if

         end do

         call check_hard_sphere(try_coords, n_atoms, reasonable)
         if (.not.reasonable .and. verbose) then
            write (150,*) "Hard sphere criterium violated!!"
         end if
         
         if (cavity) then
            call check_cavity(try_coords, cavity_radius, outside)
            if (outside .and. verbose) then
               write (150,*) "Trial move put some atoms outside!!"
            end if
            reasonable = (reasonable .and. .not.(outside)) 
         end if
         call check_compound(try_coords, max_distance, dissociated)
         if (dissociated .and. verbose) then
            write (150,*) "Cluster structure dissociated from trial move!!"
         end if
         reasonable = (reasonable .and. .not.(dissociated))

      end do
      
    end subroutine trial_move_all

    subroutine trial_move_collective_vista()
      
      implicit none

      ! imported variables
      
      ! local variables
      type (vector) :: move
      logical :: reasonable
      logical :: outside
      logical :: dissociated
      real*8 :: theta
      real*8 :: phi
      real*8 :: radius
      integer :: random_number
      logical :: found

      ! functions
      real*4 :: rand

      ! counter
      integer :: i_atom
      integer :: i_adsorp
      integer :: i_move

      reasonable = .false.
      i_move = 0
      write (150,*) "generating collective move with hard sphere prerelaxation..."
         
      try_coords = current_coords
      if (verbose) then
         write (150,*) "trial move"
      end if

      reasonable = .false.
      do while (.not.reasonable)
         
         i_move = i_move + 1
         if (i_move .gt. 1000) then
            write (150,*) "unable to generate trial move."
            write (150,*) "* Abort."
            stop
         end if

         do i_atom = 1, n_atoms, 1
         
            found = .false.
            do i_adsorp = 1, n_adsorption_atom, 1
               if (i_atom .eq. adsorption_atom(i_adsorp)) then
                  found = .true.
               end if
            end do
            
            if (.not.(adsorption_move .and. found)) then
               
               ! move each atom as many times until it is not dissociated
               ! that way surface atoms are automatically moved towards the
               ! clusters
               dissociated = .true.
               do while (dissociated)

                  ! get spherical coordinates
                  theta = pi * arch_rand()
                  phi   = 2 * pi * arch_rand()
               
                  select case (move_distribution)
                  
                  case ('uniform')
                     radius = max_move * arch_rand()
                  
                  case ('poisson')
                     call fische(1, 100d0 * max_move, random_number)
                     radius = dble(random_number) / 100d0
                  
                  case default
                     write (6,*) "Internal incosistency."
                     write (6,*) "(trial_move)"
                     write (6,*) "* Aborting."
                     stop
                     
                  end select
                  
                  if (verbose) then
                     write (150,*) "radius=", radius
                     write (150,*) "theta=", theta
                     write (150,*) "phi=", phi
                  end if
                  
                  ! move each atom as many times until it is not dissociated
                  ! that way surface atoms are automatically moved towards the
                  ! clusters
                  
                  move%x = radius * sin(theta) * cos(phi)
                  move%y = radius * sin(theta) * sin(phi)
                  move%z = radius * cos(theta)
                  
                  try_coords%coords(i_atom) = &
                       try_coords%coords(i_atom) + move
                  
                  call check_single_dissociated(try_coords, max_distance, n_atoms, i_atom, dissociated)
                  
                  if (dissociated) then
                     ! reset coordinates
                     write (150,*) i_atom, "dissociated!"
                     try_coords%coords(i_atom) = &
                          try_coords%coords(i_atom) - move  
                  end if

               end do

               if (verbose) then
                  call display(move, 150)
               end if
               
            end if
            
         end do
         
         call pre_relax_hard_sphere(try_coords, hard_sphere_radius, n_atoms)
         
         call check_hard_sphere(try_coords, n_atoms, reasonable)
         if (.not.reasonable .and. verbose) then
            write (150,*) "Hard sphere criterium violated!!"
         end if
         
         if (cavity) then
            call check_cavity(try_coords, cavity_radius, outside)
            if (outside .and. verbose) then
               write (150,*) "Trial move put some atoms outside!!"
            end if
            reasonable = (reasonable .and. .not.(outside)) 
         end if
         call check_compound(try_coords, max_distance, dissociated)
         if (dissociated .and. verbose) then
            write (150,*) "Cluster structure dissociated from trial move!!"
         end if
         reasonable = (reasonable .and. .not.(dissociated))
         
      end do
      
    end subroutine trial_move_collective_vista
    
    subroutine trial_adsorption_move()
    
      ! put the selected adsorption atom randomly on a sphere centered
      ! in the center of mass of the cluster (without this atom)
      ! with the current distance of the adsorption atom to the center
      ! as radius

      !  local variables
      type (vector) :: center
      real*8 :: inv_n_atoms
      real*8 :: radius
      real*8 :: theta
      real*8 :: phi
      logical :: found
      real*8 :: min_radius

      !  counter
      integer :: i_atom
      integer :: i_adsorp

      write (150,*) "displacing adsorption atom..."

      inv_n_atoms = 1. / (n_atoms - n_adsorption_atom)

      center%x = 0; center%y = 0; center%z = 0;
      do i_atom = 1, n_atoms, 1

         found = .false.
         do i_adsorp = 1, n_adsorption_atom, 1
            if (i_atom .eq. adsorption_atom(i_adsorp)) then
               found = .true.
            end if
         end do
         
         if (.not. found) then
            center%x = center%x + try_coords%coords(i_atom)%x
            center%y = center%y + try_coords%coords(i_atom)%y
            center%z = center%z + try_coords%coords(i_atom)%z
         end if
      end do
      center = center * inv_n_atoms

      if (verbose) then
         write (150,*) "center of mass:"
         call display(center, 150)
      end if

      ! get closest distance of adsorption atom to center of mass of cluster
      min_radius = norm(try_coords%coords(adsorption_atom(1)) - center)
      do i_adsorp = 2, n_adsorption_atom, 1
         if (norm(try_coords%coords(adsorption_atom(i_adsorp)) - center) .lt. min_radius) then
            min_radius = norm(try_coords%coords(adsorption_atom(i_adsorp)) - center)
         end if
      end do

      if (verbose) then
         write (150,*) "minimum radius of adsorption atom:", min_radius
      end if

      do i_adsorp = 1, n_adsorption_atom, 1

         if (verbose) then
            write (150,*) "coordinates of adsorption atom:"
            call display(try_coords%coords(adsorption_atom(i_adsorp)), 150)
            
            write (150,*) "adsorption atom # ", adsorption_atom(i_adsorp)
         end if
         
         radius = ad_radius_reduction(i_adsorp) * min_radius
         
         theta = pi * arch_rand()
         phi   = 2 * pi * arch_rand()

         if (verbose) then
            write (150,*) "radius=", radius
            write (150,*) "theta=", theta
            write (150,*) "phi=", phi
         end if

         try_coords%coords(adsorption_atom(i_adsorp))%x = radius * sin(theta) * cos(phi) + center%x
         try_coords%coords(adsorption_atom(i_adsorp))%y = radius * sin(theta) * sin(phi) + center%y
         try_coords%coords(adsorption_atom(i_adsorp))%z = radius * cos(theta) + center%z
         
      end do

    end subroutine trial_adsorption_move

    subroutine trial_ini_spin()

      ! counter
      integer :: i_atom

      write (150,*) "random initial spin moments per atom..."

      do i_atom = 1 , n_atoms, 1
         try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom) +  2.d0 * spin_step_width * (arch_rand() - 0.5d0)
         if (try_spin_per_atom(i_atom) .ge. ((species_z(species(i_atom)) - ini_charge_per_atom(i_atom)) / 2.d0)) then
            try_spin_per_atom(i_atom) = (species_z(species(i_atom)) - ini_charge_per_atom(i_atom)) / 2.d0
         end if
         if (try_spin_per_atom(i_atom) .le. -((species_z(species(i_atom)) - ini_charge_per_atom(i_atom)) / 2.d0)) then
            try_spin_per_atom(i_atom) = -(species_z(species(i_atom)) - ini_charge_per_atom(i_atom)) / 2.d0
         end if
         write (150,'(A,I4,F10.4)') "atom # ", i_atom, try_spin_per_atom(i_atom)
      end do

    end subroutine trial_ini_spin
    
    subroutine rotation_move()

      implicit none

      ! imported variables

      ! local variables
      real*8 :: theta
      real*8 :: phi
      type(vector) :: plane
      real*8 :: angle
      real*8 :: distance_to_origin
      logical :: reasonable
      logical :: dissociated
      logical :: outside

      ! functions
      real*4 :: rand

      ! counter
      integer :: i_atom
      integer :: i_rotated
      integer :: i_move

      write (150,*) "rotation_move ..."

      try_coords = current_coords
      ! choose plane through origin randomly
      ! center of mass of the structure is in the origin
      theta = pi * arch_rand()
      phi   = 2 * pi * arch_rand()
      
      ! choose angle to move randomly
      angle = 2.d0 * pi * arch_rand() 
      
      plane%x = sin(theta) * cos(phi)
      plane%y = sin(theta) * sin(phi)
      plane%z = cos(theta)
      
      if (verbose) then
         write (150,*) "normal vector of plane: "
         call display(plane, 150)
         
         write (150,*) "angle: ", angle
      end if
      
      i_rotated = 0
      do i_atom = 1, n_atoms, 1
         distance_to_origin = try_coords%coords(i_atom)%x * plane%x + &
              try_coords%coords(i_atom)%y * plane%y + try_coords%coords(i_atom)%z * plane%z
         if (distance_to_origin .gt. 0) then
            ! we are on one side of the plane, so rotate this point around the axis
            call rotate(plane, angle, try_coords%coords(i_atom))
            i_rotated = i_rotated + 1
         end if
      end do
      write (150,*) i_rotated, " atoms rotated."
      
      call pre_relax_hard_sphere(try_coords, hard_sphere_radius, n_atoms)
      
      call check_hard_sphere(try_coords, n_atoms, reasonable)
      if (.not.reasonable) then
         write (150,*) "Hard sphere criterium violated!!"
      end if
      
    end subroutine rotation_move
    
    subroutine evaluate_acceptance(energy_old, current_groundstate, energy_new, accepted)
      
      implicit none

      ! imported variables
      
      ! input
      real*8, intent(inout) :: energy_old
      real*8, intent(in)    :: energy_new
      real*8, intent(in)    :: current_groundstate

      ! output
      logical, intent(out) :: accepted

      ! local variables
      real*8 :: random_number
      real*8 :: energy_reference

      ! functions
      real*4 :: rand

      ! counter
      integer :: i_atom

      select case (acceptance_reference)

      case('last_energy')
         energy_reference = energy_old
   
      case ('current_groundstate')
         energy_reference = current_groundstate
         
      case default
         write (6,*) "Internal inconsistency (evaluate_acceptance)."
         write (6,*) "* Aborting."
         stop

      end select

      write (150,*) "new structure accepted?"
      write (150,*) "diff energy [eV]:", (energy_new - energy_reference) / n_atoms
      if ((energy_new / n_atoms) .le. ((energy_reference / n_atoms) + energy_threshold_theta)) then
         accepted = .true.
      else 
         random_number = arch_rand()
         write (150,*) "random number", random_number
         write (150,*) "exp()", exp( - (((energy_new - energy_reference) / n_atoms) - energy_threshold_theta) / energy_threshold_boltzmann)
         if (random_number .le. exp( - (((energy_new - energy_reference) / n_atoms) - energy_threshold_theta) / energy_threshold_boltzmann)) then
            accepted = .true.
         else
            accepted = .false.
         end if
      end if
      write (150,*) "accepted? ", accepted
      if (accepted) then
         energy_old = energy_new
         
         if (spin_polarized) then
            do i_atom = 1, n_atoms, 1
               current_spin_per_atom(i_atom) = optimized_spin_per_atom(i_atom)
               try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom)
            end do
         end if
         
         if (reset) then
            ! don't make next trial move from optimized structure but
            ! from last generated trial structure which results in 
            ! a metropolis run on the deformed PES
            current_coords = try_coords
         else
            ! always start from optimized structures which results
            ! in a hopping between minima but we are not sampling
            ! equilibrium quantities that way
            current_coords = optimized_coords
            try_coords     = current_coords
         end if
      else
         ! well, just set the structure back to last initial configuration
         try_coords = current_coords
         
         if (spin_polarized) then
            do i_atom = 1, n_atoms, 1
               try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom)
            end do
         end if

      end if
    
    end subroutine evaluate_acceptance

    subroutine select_next_starting_structure(current_groundstate, found, new_minimum)

      implicit none

      ! imported variables

      ! input variables
      logical :: found
      type (node_of_geometry), target :: new_minimum
      real*8 :: current_groundstate

      ! local
      type (node_of_geometry), pointer :: list 
      real*8, allocatable, dimension(:) :: histogram 
      integer :: i_counter
      integer :: n_structures
      integer :: histogram_max
      real*8 :: sum
      real*8 :: inv_sum
      real*8 :: propability
      integer :: n_selected
      real*8 :: random_number
      integer :: number_of_minima
      
      write (150,*) "select starting structure according to histogram-weighted statistics..."

      call get_number_of_minima(list_of_geometries, 1d6, number_of_minima)

!      write (150,*) number_of_minima

      allocate(histogram(number_of_minima))

      ! if latest minimum is completely new then automatically choose this as starting point
      ! (and it is within the current energy interval of interest)
      if ((.not.found) .and. ((new_minimum%inst_geometry%energy - current_groundstate / n_atoms) .le. energy_threshold_theta)) then
         current_coords = new_minimum%inst_geometry%coords
         new_minimum%inst_geometry%histogram_out = new_minimum%inst_geometry%histogram_out + 1 
      else
         ! choose basin according to histogram-weighted statistics
!         write (150,*) "collect histogram entries..."
         i_counter = 0
         list => list_of_geometries%head
         
         histogram_max = 0
         do while ((associated(list)) .and. (i_counter .lt. number_of_minima))

            ! if out of energy window of interest, stop counting
            if ((list%inst_geometry%energy - current_groundstate / n_atoms) .gt. energy_threshold_theta) then
               exit
            end if
            i_counter = i_counter + 1

            histogram(i_counter) = list%inst_geometry%histogram_out
!            write (150,*) i_counter, histogram(i_counter)
            if (histogram(i_counter) .gt. histogram_max) then
               histogram_max = histogram(i_counter) 
            end if

            list => list%next
         end do
         n_structures = i_counter
!         write (150,*) "n_structures...", n_structures
         ! calculate ratio of propabilities to choose the basin according to the number of histogram entries
!         write (150,*) "calculate ratio of propabilities..."
         sum = 0.d0
         do i_counter = 1, n_structures, 1
            histogram(i_counter) = dble(histogram_max) / histogram(i_counter)
            sum = sum + histogram(i_counter)
!            write (150,*) i_counter, histogram(i_counter)
         end do
         inv_sum = 1.d0 / sum

         ! normalize...
!         write (150,*) "normalize and select structure..."
         random_number = arch_rand()
         write (150,*) random_number
         i_counter = 0
         propability = 0.d0
         do while ((propability .lt. random_number) .and. (i_counter .lt. n_structures))
            i_counter = i_counter + 1
            propability = propability + inv_sum * histogram(i_counter)  
            write (150,*) i_counter, propability, inv_sum * histogram(i_counter)
         end do
         n_selected = i_counter
         write (150,*) "selected: ", n_selected
         ! go once again through linked list to get the corresponding coordinates
!         write (150,*) "get coordinates..."
         list => list_of_geometries%head
         
         do i_counter = 2, n_selected, 1
            list => list%next
         end do
         current_coords = list%inst_geometry%coords
         if (verbose) then
            write (150,*) "new starting basin(current_coords):"
            call display_array(current_coords, 150)
         end if
         list%inst_geometry%histogram_out = list%inst_geometry%histogram_out + 1

      end if

      try_coords = current_coords
      deallocate(histogram)

    end subroutine select_next_starting_structure

    subroutine optimize(final_coords, final_spin_per_atom, status)
      
      implicit none

      ! imported variables
      
      ! output
      type (vector_array), intent(out) :: final_coords
      real*8 , dimension(n_atoms) :: final_spin_per_atom
      integer, intent(inout) :: status

      ! local variables
      real*8 :: energy_old
      real*8 :: energy_new
      real*8 :: square_of_move
      logical :: accepted
      logical :: found
      logical :: outside
      logical :: successful
      type (node_of_geometry), pointer :: new_minimum
      type (geometry) :: old_minimum
      real*8 :: conversion_factor
      type (vector_array) :: forces

      real*8 :: max_force
      integer :: i_line
      integer :: i_evaluations
      integer :: number_of_minima
      integer :: number_of_evaluations
      real*8 :: mean_number_of_evaluations
      logical :: optimized
      logical :: dissociated
      real*8 :: order_parameter_o4_new
      real*8 :: order_parameter_o4_old
      integer :: new_basin
      integer :: old_basin
      real*8 :: configurational_jump
      real*8 :: average_configurational_jump
      real*8 :: successful_rate
      real*8 :: acceptance_rate
      real*8 :: add_successful_rate
      real*8 :: add_average_configurational_jump
      real*8 :: average_configurational_jump_rate
      real*8 :: add_average_configurational_jump_rate
      real*8 :: current_groundstate
      real*8 :: mean_delta_energy
      real*8 :: add_mean_delta_energy
      real*8 :: add_mean_delta_energy_rate
      real*8, dimension(:), allocatable :: nn_atomic_energy

      ! parameter for maximum rock adaptive scheme
      integer :: alpha_succ
      integer :: alpha_unsucc
      integer :: alpha_highE
      real*8 :: alpha_succ_rate
      real*8 :: alpha_unsucc_rate
      real*8 :: alpha_highE_rate
      real*8 :: ini_delta_max_move
      real*8 :: alpha_succ_old
      integer :: atom_index

      ! counters
      integer :: i_accepted
      integer :: i_loop
      integer :: i_minima
      integer :: i_atom
      integer :: i_successful_add
      integer :: i_successful
      integer :: i_distance
      integer :: i_average
      integer :: i_average_loop

      call allocate_array(forces, n_atoms)
      allocate(nn_atomic_energy(n_atoms))

      i_accepted = 0
      i_successful = 0
      i_successful_add = 0
      i_minima   = 0
      i_line = 0
      i_evaluations = 0
      i_average = 0
      number_of_evaluations = 0
      add_successful_rate = 0
      add_average_configurational_jump = 0
      mean_delta_energy = 0.d0
      add_mean_delta_energy = 0.d0
      i_average_loop = 0

      ! maximum rock parameter
      alpha_succ = 0
      alpha_unsucc = 0
      alpha_highE = 0
      
      open (150, FILE="bh_log.out")

      allocate(list_of_geometries%head) 
      nullify(list_of_geometries%head%next)
      call g_allocate(list_of_geometries%head%inst_geometry)

      allocate(new_minimum)
      nullify(new_minimum%next)
      call g_allocate(new_minimum%inst_geometry)

      call g_allocate(old_minimum)

      call initialize_positions()
      call center_structure(current_coords)

      try_coords = current_coords
      ini_delta_max_move = delta_max_move
      alpha_succ_old = -1.d0

      write (150,*) "optimization of initial structure..."
      write (150,*) "doing bh-loop # 1"
      if (verbose) then
         write (150,*) "try_coords"
      end if

      status = 2 
      do while (status .eq. 2)

         if (verbose) then
            call display_array(try_coords, 150)
         end if
         select case (minimizer)
            
         case ('external')
            call get_dft_data(optimized_coords, energy_old, forces, optimized_spin_per_atom, status)
            
            if (status .eq. 1) then
               final_coords = try_coords
               
               if (spin_polarized) then
                  do i_atom = 1, n_atoms, 1
                     final_spin_per_atom(i_atom) = try_spin_per_atom(i_atom)
                  end do
               end if
               
               return
            end if
            
            ! if initial structure didn't converge, try a move
            if (status .eq. 2) then
               
               select case (kind_of_move)
               
               case ('single')
                  call trial_move(square_of_move)
                  
               case ('worst-single')
                  call trial_move(square_of_move)
                  
               case ('collective')
                  call trial_move_all()

               case ('collective-vista')
                  call trial_move_collective_vista()
                  
               case default
                  write (150,*) "Internal inconsistency in basin-hopping!"
                  write (150,*) "(trial moves)"
                  write (150,*) "* Aborting."
                  stop   
                  
               end select

               if (spin_polarized) then
                  call trial_ini_spin()
               end if

            end if

         case ('cg')
            call cg_optimize(try_coords, optimized_coords, energy_old, forces, max_force, i_line, i_evaluations, status)
            number_of_evaluations = number_of_evaluations + i_evaluations
            
         case ('magic')
            call nn_optimize(try_coords, optimized_coords, energy_old, nn_atomic_energy, status)

         case default
            write (150,*) "Internal inconsistency in basin-hopping!"
            write (150,*) "* Aborting."
            stop
            
         end select
         
      end do

      current_groundstate = energy_old

      if (verbose) then
         write (150,*) "optimized_coords"
         call display_array(optimized_coords, 150)
      end if

      list_of_geometries%head%inst_geometry%energy = energy_old / n_atoms
      list_of_geometries%head%inst_geometry%coords = optimized_coords
      list_of_geometries%head%inst_geometry%found_in_loop = 1
      list_of_geometries%head%inst_geometry%histogram_out = 1

      if (spin_polarized) then
         if (verbose) then
            write (150,*) "spin per atom"
         end if
         list_of_geometries%head%inst_geometry%total_moment = 0.d0
         do i_atom = 1, n_atoms, 1
            if (verbose) then
               write (150,*) optimized_spin_per_atom(i_atom)
            end if
            list_of_geometries%head%inst_geometry%spin_per_atom(i_atom) = optimized_spin_per_atom(i_atom)
            list_of_geometries%head%inst_geometry%total_moment = list_of_geometries%head%inst_geometry%total_moment + optimized_spin_per_atom(i_atom)
         end do
      end if

      call get_distances(list_of_geometries%head%inst_geometry)
      call get_angles(list_of_geometries%head%inst_geometry)
      call update_energy_interval(list_of_geometries%head%inst_geometry, energy_old / n_atoms)

      if (reset) then
         current_coords = try_coords
      else
         current_coords = optimized_coords
      end if
      try_coords = current_coords
      
      i_minima = i_minima + 1
      old_basin = 1
      average_configurational_jump = 0
      call copy_geometry(old_minimum, list_of_geometries%head%inst_geometry)
      write (150,*) "energy of head:", list_of_geometries%head%inst_geometry%energy
      write (150,*) "energy of old minimum:", old_minimum%energy
      if (verbose) then
         write (150,*) "current_coords"
         call display_array(current_coords, 150)
      end if

      ! count initial optimization as move # 1
      do i_loop = 2, max_loop, 1
         i_average_loop = i_average_loop + 1
         write (150,*) "doing bh-loop #", i_loop
         write (150,*) "------------------------"

         ! switch acceptance reference ?
         if (i_loop .eq. scf_loop_to_switch) then
            acceptance_reference = 'current_groundstate'
            energy_threshold_theta = switched_E_theta
            energy_threshold_boltzmann = switched_E_boltzmann
            write (150,'(1X,A)') "switching acceptance criteria to current_groundstate with..."
            write (150,'(1X,A,F8.4,A)') "energy_threshold_theta = ", energy_threshold_theta, " eV." 
            write (150,'(1X,A,F8.4,A)') "energy_threshold_boltzmann = ", energy_threshold_boltzmann, " eV."
         end if

         if ((use_rotation_move) .and. (mod(i_loop, fraction_of_rotation_move) .eq. 0)) then

            call rotation_move()

         else if (pure_adsorption_move) then
            
            call trial_adsorption_move()

         else if (use_angular_move) then

            call check_nn_atomic_energies(nn_atomic_energy, energy_ratio_ang_move, atom_index)
            if (atom_index .gt. 0) then
               call angular_move(atom_index)
            end if

         else
            select case (kind_of_move)
               
            case ('single')
               call trial_move(square_of_move)
               
            case ('worst-single')
               call trial_move(square_of_move)
               
            case ('collective')
               call trial_move_all()
              
            case ('collective-vista')
               call trial_move_collective_vista()
 
            case default
               write (150,*) "Internal inconsistency in basin-hopping!"
               write (150,*) "(trial moves)"
               write (150,*) "* Aborting."
               stop   
               
            end select
            
            if (adsorption_move) then
               call trial_adsorption_move()
            end if
            
         end if

         if (verbose) then
            write (150,*) "try_coords"
            call display_array(try_coords, 150)
         end if

         if (spin_polarized) then
            call trial_ini_spin()
         end if

         select case (minimizer)

         case ('external')
            call get_dft_data(optimized_coords, energy_new, forces, optimized_spin_per_atom, status)

         case ('cg')
            call cg_optimize(try_coords, optimized_coords, energy_new, forces, max_force, i_line, i_evaluations, status)
            number_of_evaluations =  number_of_evaluations + i_evaluations

         case ('magic')
            call nn_optimize(try_coords, optimized_coords, energy_new, nn_atomic_energy, status)

         case default
            write (150,*) "Internal inconsistency in basin-hopping!"
            write (150,*) "(optimizer)"
            write (150,*) "* Aborting."
            stop

         end select
 
         if (status .eq. 1) then
            ! need new energy and forces by aims
            call write_geometry_of_structures(list_of_geometries)
            call write_aims_files_of_structures(list_of_geometries)
            call write_angular_files_of_structures(list_of_geometries)

            final_coords = try_coords

            if (spin_polarized) then
               do i_atom = 1, n_atoms, 1
                  final_spin_per_atom(i_atom) = try_spin_per_atom(i_atom)
               end do
            end if

            return
         end if

         if (status .eq. 2) then
            write (150,*) "Something went wrong with the relaxation scheme."
            write (150,*) "Discarding trial move."
            try_coords = current_coords
            if (spin_polarized) then
               do i_atom = 1, n_atoms, 1
                  try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom)
               end do
            end if
            cycle
         end if
         
         ! check if maybe the cluster has dissociated during relaxation
         if (status .eq. 0) then
            call check_compound(optimized_coords, max_distance, dissociated)
            if (dissociated) then
               write (150,*) "Cluster has dissociated during local relaxation!"
               write (150,*) "Discarding trial move."
               try_coords = current_coords
               if (spin_polarized) then
                  do i_atom = 1, n_atoms, 1
                     try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom)
                  end do
               end if
               cycle
            end if
         end if

         if (verbose) then
            write (150,*) "optimized_coords"
            call display_array(optimized_coords, 150)
         end if

         if (status .eq. 0) then
            write (150,*) "energy_of_optimized_coords ", energy_new / n_atoms
         end if

         mean_delta_energy = mean_delta_energy + (energy_new - energy_old) * (energy_new - energy_old)

         ! keep track of current groundstate
         if (energy_new .lt. current_groundstate) then
            current_groundstate = energy_new
         end if

         new_minimum%inst_geometry%coords = optimized_coords
         new_minimum%inst_geometry%energy = energy_new / n_atoms
         new_minimum%inst_geometry%found_in_loop = i_loop
         new_minimum%inst_geometry%histogram_out = 0

         if (spin_polarized) then
            if (verbose) then
               write (150,*) "spin per atom"
            end if
            new_minimum%inst_geometry%total_moment = 0.d0
            do i_atom = 1, n_atoms, 1
               if (verbose) then
                  write (150,*) optimized_spin_per_atom(i_atom)
               end if
               new_minimum%inst_geometry%spin_per_atom(i_atom) = optimized_spin_per_atom(i_atom)
               new_minimum%inst_geometry%total_moment = new_minimum%inst_geometry%total_moment + optimized_spin_per_atom(i_atom)
            end do
         end if

         call get_distances(new_minimum%inst_geometry)
         call get_angles   (new_minimum%inst_geometry)
         call get_order_parameter_o4(new_minimum%inst_geometry, order_parameter_o4_new)

         call check_and_add_list_of_minima_v2(list_of_geometries, new_minimum, found)
         call get_index_of_minimum(list_of_geometries, new_minimum, new_basin)
         call get_number_of_minima(list_of_geometries, energy_interval, number_of_minima)
         if (verbose) then
            write (150,'(1X, A, I10, 1X, I10, 1X, L4)') "old basin / new basin / found:", old_basin, new_basin, found
         end if

         successful = .true.
         if (.not. found) then
            i_successful = i_successful + 1
            i_successful_add = i_successful_add + 1
         else
            if (new_basin .ne. old_basin) then
               i_successful = i_successful + 1
               i_successful_add = i_successful_add + 1
            else
               alpha_unsucc = alpha_unsucc + 1
               successful = .false.
            end if
         end if

         successful_rate = dble(i_successful) / dble(i_loop - 1)

         if (number_of_averaging .gt. 0) then
            if (mod((i_loop - 1), number_of_averaging) .eq. 0) then
               ! new averaging period starts, so start counting successful jumps from zero
               i_successful_add = 0
            else
               add_successful_rate = dble(i_successful_add) / dble(mod((i_loop - 1), number_of_averaging))
            end if
         end if

         call diff_distances(new_minimum%inst_geometry, old_minimum, configurational_jump)

         average_configurational_jump = average_configurational_jump + configurational_jump
         average_configurational_jump_rate = average_configurational_jump / dble(i_loop - 1)

         add_average_configurational_jump = add_average_configurational_jump + configurational_jump
         add_mean_delta_energy = add_mean_delta_energy + (energy_new - energy_old) * (energy_new - energy_old)
         write (150,*) "mean_delta_energy ", add_mean_delta_energy
         if (number_of_averaging .gt. 0) then
            if (mod((i_loop - 1), number_of_averaging) .eq. 0) then
               ! new averaging period starts, so start counting successful jumps from zero
               add_average_configurational_jump = 0.d0
               add_mean_delta_energy = 0.d0
            else
               add_average_configurational_jump_rate = add_average_configurational_jump / dble(mod((i_loop - 1), number_of_averaging))
               add_mean_delta_energy_rate = sqrt(add_mean_delta_energy / dble(mod((i_loop - 1), number_of_averaging))) / dble(n_atoms)
            end if
         end if

         select case (bh_flavour)

         case ('standard')
            call evaluate_acceptance(energy_old, current_groundstate, energy_new, accepted)

         case ('histogram-weighted')
            ! if new structure has been found then automatically select this one as next starting point
            ! since number_of_minima is set from last cycle, it could be increased by one
            ! it's only used for the dimension of the histogram-array within select_next_starting_structure
            ! so it doesn't matter here
            call select_next_starting_structure(current_groundstate, found, new_minimum)

         case default
            write (150,*) "Internal inconsistency in basin-hopping!"
            write (150,*) "(bh-flavour)"
            write (150,*) "* Aborting."
            stop   

         end select

         if (accepted) then
            ! check, whether successful or not
            if (successful) then
               alpha_succ = alpha_succ + 1
            end if
            i_accepted = i_accepted + 1
            order_parameter_o4_old = order_parameter_o4_new
            old_basin = new_basin
            call copy_geometry(old_minimum, new_minimum%inst_geometry)
            write (150,'(1X, A, I10, 1X, F20.10)') "current energy #", i_loop, energy_new / n_atoms
         else
            alpha_highE = alpha_highE + 1
            order_parameter_o4_new = order_parameter_o4_old
            write (150,'(1X, A, I10, 1X, F20.10)') "current energy #", i_loop, energy_old / n_atoms
         end if


         !! maximum rock adaptive scheme !!

         alpha_unsucc_rate = dble(alpha_unsucc) / dble(i_average_loop)
         alpha_succ_rate = dble(alpha_succ) / dble(i_average_loop)
         alpha_highE_rate = dble(alpha_highE) / dble(i_average_loop)
         
         if ((move_adaptive) .and.(adapting_scheme .eq. 'maximum_rock')) then

            ! coarse reset of maximum move if completely nonsense

            if (mod(i_average_loop, number_of_coarse_averaging).eq.0) then

               if (alpha_unsucc_rate .gt. alpha_unsucc_max) then
                  ! give it a kick in the ass and reset averaging
                  max_move = max_move + kick_in_the_ass
                  alpha_succ_old = -1.d0
                  write (150,'(1X,A,F8.4)') "maximum rock: reset max move due to kick-in-the-ass to ", max_move
                  delta_max_move = abs(ini_delta_max_move)
                  alpha_unsucc = 0
                  alpha_succ = 0
                  alpha_highE = 0
                  i_average_loop = 0
               end if
               if (alpha_highE_rate .gt. alpha_highE_max) then
                  ! give it a kick in the ass and reset averaging
                  max_move = max_move - kick_in_the_ass
                  alpha_succ_old = -1.d0
                  write (150,'(1X,A,F8.4)') "maximum rock: reset max move due to kick-in-the-ass to ", max_move
                  delta_max_move = -abs(ini_delta_max_move)
                  alpha_unsucc = 0
                  alpha_succ = 0
                  alpha_highE = 0
                  i_average_loop = 0
               end if
            end if
            
            if (i_average_loop .eq. number_of_averaging) then

               ! new averaging period starts, so start counting successful jumps from zero
               alpha_unsucc = 0
               alpha_succ = 0
               alpha_highE = 0
               
               if (alpha_succ_old .ge. 0) then
                  if (abs(alpha_succ_rate - alpha_succ_old) .lt. alpha_tolerance) then
                     ! move distance seems to make sense so vary in smaller steps
                     delta_max_move = delta_max_move * max_move_modifier_one
                     if (abs(delta_max_move) .lt. min_delta_max_move) then
                        delta_max_move = min_delta_max_move
                     end if
                  end if
                  if (alpha_succ_rate .gt. (alpha_succ_old + alpha_tolerance)) then
                     ! move distance is obviously corrected in the right direction
                     delta_max_move = delta_max_move * max_move_modifier_two
                     if (abs(delta_max_move) .gt. max_delta_max_move) then
                        delta_max_move = max_delta_max_move
                     end if
                  end if
                  if (alpha_succ_rate .lt. (alpha_succ_old - alpha_tolerance)) then
                     ! move distance significantly worse, so change direction
                     delta_max_move = - delta_max_move * max_move_modifier_one
                     if (abs(delta_max_move) .lt. min_delta_max_move) then
                        delta_max_move = min_delta_max_move
                     end if
                  end if
               end if
               
               max_move = max_move + delta_max_move
               write (150,'(A)') "-------------------------------------------"
               write (150,'(A,F8.4)') "reset move distance according to maximum rock scheme to ", max_move
               write (150,'(A,F8.4,1X,F8.4,1X,F8.4)') "delta_max_move, alpha_succ_old, alpha_succ_rate ", delta_max_move, alpha_succ_old, alpha_succ_rate
               write (150,'(A)') "-------------------------------------------"
               ! store old alpha_succ
               alpha_succ_old = alpha_succ_rate
               i_average_loop = 0
            end if
         end if
         
         acceptance_rate = dble(i_accepted) / dble(i_loop - 1)

         if (.not. found) then
            i_minima = i_minima + 1
            allocate(new_minimum)
            nullify(new_minimum%next)
            call g_allocate(new_minimum%inst_geometry)
         end if

         ! adapt moves if desired
         if (move_adaptive) then
            if (mod((i_loop - 1), number_of_averaging) .eq. 0) then
               select case(adapting_scheme)

               case ('successful_additive')
               
                  if (add_successful_rate .gt. max_successful_rate) then
                     max_move = max_move - max_move_decrease
                     write (150,'(1X, A, F10.4, A)') "Maximum move decreased to ", max_move, " Ang."
                  else if (add_successful_rate .lt. min_successful_rate) then
                     max_move = max_move + max_move_increase
                     write (150,'(1X, A, F10.4, A)') "Maximum move increased to ", max_move, " Ang."
                  end if
                  
               case ('maximum_rock')
                  continue
                  
               case default
                  
                  write (6,*) "Internal inconsistency (move_adaptive)."
                  write(6,*) "* Aborting."
                  stop
                  
               end select
            end if
         end if
         write (150,'(1X,A,F10.4)') "current maximum move [Ang]:", max_move
         write (150,'(1X,A)') "---------------------------------------------"
         write (150,'(1X,A,F8.4)') "alpha_unsucc  :", alpha_unsucc_rate
         write (150,'(1X,A,F8.4)') "alpha_highE   :", alpha_highE_rate
         write (150,'(1X,A,F8.4)') "alpha_succ    :", alpha_succ_rate
         write (150,'(1X,A)') "---------------------------------------------"
         write (150,'(1X,A,F8.4)') "consistency check :", alpha_unsucc_rate + alpha_highE_rate + alpha_succ_rate 
         write (150,'(1X,A)') "---------------------------------------------"
         write (150,'(1X,A,F10.4)') "current order parameter o4                            : ", order_parameter_o4_new
         write (150,'(1X,A,F10.4)') "current acceptance rate                               : ", acceptance_rate
         write (150,'(1X,A,F10.4)') "overall successful rate                               : ", successful_rate
         write (150,'(1X,A,E10.4)') "last configurational jump                             : ", configurational_jump
         write (150,'(1X,A,E10.4)') "overall average configurational jump                  : ", average_configurational_jump_rate
         write (150,'(1X,A,E10.4)') "overall average sqrt(diff energy square)              : ", sqrt(mean_delta_energy / dble(i_loop - 1)) / dble(n_atoms)
         if (number_of_averaging .gt. 0) then
            write (150,'(1X,A,F10.4)') "current successful rate for last averaging period         : ", add_successful_rate
            write (150,'(1X,A,E10.4)') "average configurational jump for last averaging period    : ", add_average_configurational_jump_rate
            write (150,'(1X,A,E10.4)') "average sqrt(diff energy square) for last averaging period: ", add_mean_delta_energy_rate
         end if

         write (150,'(1X,A,I4, 1X, I4)') "number of different minima found # ", i_loop, number_of_minima
         if (potential_flag .ne. 'external') then
            write (150,'(1X,A,I4, 1X, F6.2)') "mean number of energy evaluations per relaxation # ", i_loop, dble(number_of_evaluations) / dble(i_loop)
         end if

         call write_energy_of_structures(list_of_geometries)
         call center_structure(current_coords)
      
         if (verbose) then
            write (150,*) "current_coords (beware: centered!!)"
            call display_array(current_coords, 150)
         end if
         
         try_coords = current_coords
         call compare_lowest_energy(list_of_geometries, target_energy, optimized)
         if ((number_of_minima .eq. target_number_of_minima) .and. (optimized)) then
            ! if a certain number of minima that are desired have been identified then stop
            ! the run!!!
            exit
         end if
         if (optimize_energy) then
            if (optimized) then
               write (150,*)"target energy reached # ", i_loop
               optimize_energy = .false.
            end if
         end if

         if (minimizer .eq. 'magic') then
            call write_geometry_of_structures(list_of_geometries)
         end if

      end do
      call write_geometry_of_structures(list_of_geometries)
      call write_aims_files_of_structures(list_of_geometries)
      call write_angular_files_of_structures(list_of_geometries)
      final_coords = current_coords
      status = 3
      close (150)
      call deallocate_array(forces)
      if (allocated(nn_atomic_energy)) then
         deallocate(nn_atomic_energy)
      end if

    end subroutine optimize

    subroutine analyze_moves()

      ! local variables
      type (node_of_geometry), pointer :: new_minimum
      real*8 :: energy
      character*30 :: desc_str
      integer :: i_code
      logical :: eof
      integer :: index
      integer :: basin_old
      integer :: basin_new
      integer :: successful_moves
      integer :: number_of_minima
      real*8 :: mean_value
      real*8 :: standard_deviation
      logical :: discard
      logical :: rotation

      ! counter
      integer :: i_atom
      integer :: i_loop
      integer :: i_distance

      i_loop = -1

      allocate(new_minimum)
      nullify(new_minimum%next)
      call g_allocate(new_minimum%inst_geometry)

      open (10, FILE= "move.out")
      open (7, FILE= "bh_log.out")
      open (15, FILE= "lookup.out")
      open (150, FILE="bh_analysis.out")

      ! read first line
      read (7, *, iostat = i_code) desc_str
      
      if (i_code.ne.0) then
         eof = .true. 
      end if

      successful_moves = 0
      rotation = .false.
   
      do while (.not.eof)
         
         ! decode contents of current line
         select case(desc_str)
        
         case ('current_coords')

            do i_atom = 1, n_atoms, 1
               read (7,*) current_coords%coords(i_atom)%x, current_coords%coords(i_atom)%y, current_coords%coords(i_atom)%z  
            end do

            new_minimum%inst_geometry%coords = current_coords
            call get_distances(new_minimum%inst_geometry)
            call get_index_of_minimum(list_of_geometries, new_minimum, index)
            call increase_histogram_out(list_of_geometries, index)
!            write (6,*) "after current_coords"
!            do i_distance = 1, new_minimum%inst_geometry%n_distances, 1
!               write (6,*) i_distance, new_minimum%inst_geometry%distance(i_distance)
!            end do

            if (index .eq. -1) then
               write (6,*) "Internal inconsistency in basin-hopping."
               write (6,*) "after current_coords..."
               write (6,*) "Analyze_moves: structure could not be found in linked list."
               stop
            end if
            
            basin_old = index

            write (15,'(I4)') basin_old

         case ('try_coords')

            discard = .false.
            do i_atom = 1, n_atoms, 1
               read (7,*) try_coords%coords(i_atom)%x, try_coords%coords(i_atom)%y, try_coords%coords(i_atom)%z  
            end do

         case ('rotation_move')

            rotation = .true.

         case ('Discarding')

            write (15,'(I4)') -1

         case ('optimized_coords')

            i_loop = i_loop + 1
            do i_atom = 1, n_atoms, 1
               read (7,*) optimized_coords%coords(i_atom)%x, optimized_coords%coords(i_atom)%y, optimized_coords%coords(i_atom)%z  
            end do
            ! call display_array(optimized_coords)
            new_minimum%inst_geometry%coords = optimized_coords
            call get_distances(new_minimum%inst_geometry)
            call get_index_of_minimum(list_of_geometries, new_minimum, index)
!            write (6,*) "after optimized_coords"
!            do i_distance = 1, new_minimum%inst_geometry%n_distances, 1
!               write (6,*) i_distance, new_minimum%inst_geometry%distance(i_distance)
!            end do
            if (index .eq. -1) then
               write (6,*) "Internal inconsistency in basin-hopping."
               write (6,*) "after optimized_coords..."
               write (6,*) "Analyze_moves: structure could not be found in linked list."
               stop
            end if
            
            basin_new = index
            
            if (i_loop .gt. 0) then
               if (basin_new .ne. basin_old) then
                  successful_moves = successful_moves + 1
               end if
               ! i_loop eq. 0 is the initial relaxation, so this output doesn't make sense
               if (rotation) then
                  write (10, '(A,I4,A,I4,A,I4,A)') "bh-loop # ", i_loop, " minimum # ", basin_old, " -> minimum # ", basin_new, "   rotation"
               else
                  write (10, '(A,I4,A,I4,A,I4,A)') "bh-loop # ", i_loop, " minimum # ", basin_old, " -> minimum # ", basin_new, "   displacement"
               end if
               write (10, '(1X,I4)') n_atoms
               write (10,*)
               do i_atom = 1, n_atoms, 1
                  write (10,'(A, 3(F10.6,2X))') species_name(species(i_atom)), current_coords%coords(i_atom)%x, current_coords%coords(i_atom)%y, &
                       current_coords%coords(i_atom)%z
               end do
               
               write (10, '(1X,I4)') n_atoms
               write (10,*)
               do i_atom = 1, n_atoms, 1
                  write (10,'(A, 3(F10.6,2X))') species_name(species(i_atom)), try_coords%coords(i_atom)%x, try_coords%coords(i_atom)%y, try_coords%coords(i_atom)%z
               end do
               
               write (10, '(1X,I4)') n_atoms
               write (10,*)
               do i_atom = 1, n_atoms, 1
                  write (10,'(A, 3(F10.6,2X))') species_name(species(i_atom)), optimized_coords%coords(i_atom)%x, optimized_coords%coords(i_atom)%y, &
                       optimized_coords%coords(i_atom)%z
               end do
               write (10,*) i_loop, "# of successful moves: ", successful_moves, " (", successful_moves / dble(i_loop), ")" ;
            end if
            rotation = .false.

         case ('energy_of_optimized_coords')
   
            backspace(7)
     
            read (7,*) desc_str, energy

         end select
         
         ! read next line
         
         read (7,*,iostat = i_code) desc_str
         
         if (i_code.ne.0) then
            eof = .true.
         end if
         
      end do

      call analyze_efficiency_of_sampling(list_of_geometries, number_of_minima, mean_value, standard_deviation)

      call write_energy_of_structures(list_of_geometries)

      write (10,*) "# of different minima found: ", number_of_minima
      write (10,*) "mean_value of histogram: ", mean_value
      write (10,'(A,F8.4,A,F8.4)') "standard deviation: ", standard_deviation, " stan.dev./mean_value:", standard_deviation / mean_value
      write (10,*) "# of successful moves: ", successful_moves, " (", successful_moves / dble(i_loop), ")" 

      close (7)
      close (10)
      close (15)
      close (150)

    end subroutine analyze_moves
 
    subroutine rebuild_linked_list()

      ! local variables
      type (vector_array) :: forces
      real*8 :: energy
      type (node_of_geometry), pointer :: new_minimum
      integer :: status
      logical :: found

      ! counter
      integer :: i_loop

      open (150, FILE="bh_analysis.out")
      call allocate_array(forces, n_atoms)

      allocate(list_of_geometries%head) 
      nullify(list_of_geometries%head%next)
      call g_allocate(list_of_geometries%head%inst_geometry)

      allocate(new_minimum)
      nullify(new_minimum%next)
      call g_allocate(new_minimum%inst_geometry)

      call get_dft_data(optimized_coords, energy, forces, optimized_spin_per_atom, status)
      if (status .gt. 0) then
         return
      end if

      list_of_geometries%head%inst_geometry%energy = energy / n_atoms
      list_of_geometries%head%inst_geometry%coords = optimized_coords
      call get_distances(list_of_geometries%head%inst_geometry)
      call get_angles(list_of_geometries%head%inst_geometry)
      call update_energy_interval(list_of_geometries%head%inst_geometry, energy / n_atoms)

      do i_loop = 1, max_loop, 1
         call get_dft_data(optimized_coords, energy, forces, optimized_spin_per_atom, status) 
         
         if (status .eq. 1) then
            ! dft-file has ended here
            exit
         end if
         if (status .eq. 2) then
            ! Something went wrong with the relaxation scheme.
            ! Discarding this configuration
            cycle
         end if
         
         new_minimum%inst_geometry%coords = optimized_coords
         new_minimum%inst_geometry%energy = energy / n_atoms
         
         call get_distances(new_minimum%inst_geometry)
         call get_angles   (new_minimum%inst_geometry)
         call check_and_add_list_of_minima_v2(list_of_geometries, new_minimum, found)
         
         if (.not. found) then
            allocate(new_minimum)
            nullify(new_minimum%next)
            call g_allocate(new_minimum%inst_geometry)
         end if

      end do

      call write_energy_of_structures(list_of_geometries)
      call write_geometry_of_structures(list_of_geometries)
      call write_aims_files_of_structures(list_of_geometries)
      call write_angular_files_of_structures(list_of_geometries)
      call deallocate_array(forces) 

      close (150)

    end subroutine rebuild_linked_list
    
    subroutine sample_hop_matrix(final_coords, final_spin_per_atom, status)

      ! output
      type (vector_array), intent(out) :: final_coords
      real*8 , dimension(n_atoms) :: final_spin_per_atom
      integer, intent(inout) :: status

      type (vector_array), dimension(n_starting_isomers) :: coords

      ! counter
      integer :: i_starting_isomer
      integer :: i_sample
      integer :: i_loop
      integer :: i_atom

      ! local variables
      real*8 :: square_of_move
      logical :: found
      logical :: outside
      type (node_of_geometry), pointer :: new_minimum
      type (geometry) :: old_minimum
      real*8 :: conversion_factor
      type (vector_array) :: forces
      real*8 :: energy_new
      real*8 :: max_force
      integer :: i_line
      integer :: i_evaluations
      integer :: number_of_minima
      integer :: number_of_evaluations
      real*8 :: mean_number_of_evaluations
      logical :: optimized
      logical :: dissociated
      real*8 :: order_parameter_o4_new

      integer, dimension(n_starting_isomers) :: start_sample
      
      allocate(new_minimum)
      nullify(new_minimum%next)
      call g_allocate(new_minimum%inst_geometry)

      call allocate_array(forces, n_atoms)

      call initialize_positions()

      open (150, FILE="bh_log.out")

      ! first get coordinates of basins from which to jump
      do i_starting_isomer = 1, n_starting_isomers, 1
         call allocate_array(coords(i_starting_isomer), n_atoms)
         call get_coords(list_of_geometries, starting_index(i_starting_isomer), coords(i_starting_isomer), start_sample(i_starting_isomer))
         write (150,*) "coordinates of starting_isomer # ", starting_index(i_starting_isomer)
         call display_array(coords(i_starting_isomer), 150)
         write (150,*) start_sample(i_starting_isomer), " moves already done."
      end do

      do i_starting_isomer = 1, n_starting_isomers, 1

         write (150,*) "start from isomer #  ", starting_index(i_starting_isomer)

         do i_sample = start_sample(i_starting_isomer)+1, n_samples, 1

            write (150,*) "do move # ", i_sample
            ! set current_coord to starting coordinates
            current_coords = coords(i_starting_isomer)
            try_coords = current_coords
            
            write (150,*) "current_coords (beware: centered!!)"
            call display_array(current_coords, 150)
            
            if ((use_rotation_move) .and. (mod(i_loop, fraction_of_rotation_move) .eq. 0)) then
               
               call rotation_move()
               
            else if (pure_adsorption_move) then
               
               call trial_adsorption_move() 
               
            else
               select case (kind_of_move)
                  
               case ('single')
                  call trial_move(square_of_move)
                  
               case ('worst-single')
                  call trial_move(square_of_move)
                  
               case ('collective')
                  call trial_move_all()
                  
               case default
                  write (150,*) "Internal inconsistency in basin-hopping!"
                  write (150,*) "(trial moves)"
                  write (150,*) "* Aborting."
                  stop   
                  
               end select
               
               if (adsorption_move) then
                  call trial_adsorption_move()
               end if
               
            end if
         
            write (150,*) "try_coords"
            call display_array(try_coords, 150)

            select case (minimizer)
               
            case ('external')
               call get_dft_data(optimized_coords, energy_new, forces, optimized_spin_per_atom, status)

            case ('cg')
               call cg_optimize(try_coords, optimized_coords, energy_new, forces, max_force, i_line, i_evaluations, status)
               number_of_evaluations =  number_of_evaluations + i_evaluations
               
            case default
               write (150,*) "Internal inconsistency in basin-hopping!"
               write (150,*) "(optimizer)"
               write (150,*) "* Aborting."
               write (150,*) minimizer
               stop
               
            end select

            if (status .eq. 1) then
               ! need new energy and forces by aims
               call write_geometry_of_structures(list_of_geometries)
               call write_aims_files_of_structures(list_of_geometries)
               call write_angular_files_of_structures(list_of_geometries)
               
               final_coords = try_coords
               
               if (spin_polarized) then
                  do i_atom = 1, n_atoms, 1
                     final_spin_per_atom(i_atom) = try_spin_per_atom(i_atom)
                  end do
               end if
               
               return
            end if
            
            if (status .eq. 2) then
               write (150,*) "Something went wrong with the relaxation scheme."
               write (150,*) "Discarding trial move."
               try_coords = current_coords
               if (spin_polarized) then
                  do i_atom = 1, n_atoms, 1
                     try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom)
                  end do
               end if
               cycle
            end if

            ! check if maybe the cluster has dissociated during relaxation
            if (status .eq. 0) then
               call check_compound(optimized_coords, max_distance, dissociated)
               if (dissociated) then
                  write (150,*) "Cluster has dissociated during local relaxation!"
                  write (150,*) "Discarding trial move."
                  try_coords = current_coords
                  if (spin_polarized) then
                     do i_atom = 1, n_atoms, 1
                        try_spin_per_atom(i_atom) = current_spin_per_atom(i_atom)
                     end do
                  end if
                  cycle
               end if
            end if

            write (150,*) "optimized_coords"
            call display_array(optimized_coords, 150)
            
            if (status .eq. 0) then
               write (150,*) "energy_of_optimized_coords ", energy_new / n_atoms
            end if
            
            new_minimum%inst_geometry%coords = optimized_coords
            new_minimum%inst_geometry%energy = energy_new / n_atoms
            new_minimum%inst_geometry%found_in_loop = i_loop
            new_minimum%inst_geometry%histogram_out = 0
            
            if (spin_polarized) then
               write (150,*) "spin per atom"
               new_minimum%inst_geometry%total_moment = 0.d0
               do i_atom = 1, n_atoms, 1
                  write (150,*) optimized_spin_per_atom(i_atom)
                  new_minimum%inst_geometry%spin_per_atom(i_atom) = optimized_spin_per_atom(i_atom)
                  new_minimum%inst_geometry%total_moment = new_minimum%inst_geometry%total_moment + optimized_spin_per_atom(i_atom)
               end do
            end if
            
            call get_distances(new_minimum%inst_geometry)
            call get_angles   (new_minimum%inst_geometry)
            call get_order_parameter_o4(new_minimum%inst_geometry, order_parameter_o4_new)
            
            call check_and_add_list_of_minima_v2(list_of_geometries, new_minimum, found)

            if (.not. found) then
               allocate(new_minimum)
               nullify(new_minimum%next)
               call g_allocate(new_minimum%inst_geometry)
            end if
            
         end do
      end do
      final_coords = current_coords
      status = 3
      close (150)

    end subroutine sample_hop_matrix

    subroutine check_nn_atomic_energies(nn_atomic_energy, energy_ratio_ang_move, atom_index)

      implicit none

      ! imported variables
      
      ! input 
      real*8, intent(in) :: nn_atomic_energy(n_atoms)
      real*8, intent(in) :: energy_ratio_ang_move
      
      ! output
      integer, intent(out) :: atom_index
      
      ! local variables
      real*8 :: lowest_energy
      real*8 :: highest_energy

      ! counter
      integer :: i_atom

      atom_index = -1

      ! get lowest_energy
      lowest_energy = nn_atomic_energy(1)
      do i_atom = 2, n_atoms, 1
         if (nn_atomic_energy(i_atom) .lt. lowest_energy) then
            lowest_energy = nn_atomic_energy(i_atom)
         end if
      end do

      ! get highest_energy 
      highest_energy = nn_atomic_energy(1)
      do i_atom = 2, n_atoms, 1
         if (nn_atomic_energy(i_atom) .gt. highest_energy) then
            highest_energy = nn_atomic_energy(i_atom)
            atom_index = i_atom
         end if
      end do

      ! check ratio
      if ((atom_index .gt. 0) .and. (highest_energy / lowest_energy .gt. energy_ratio_ang_move)) then
         atom_index = -1
      end if
      if (verbose) then
         write (150,*) highest_energy, lowest_energy, highest_energy / lowest_energy, atom_index
      end if

    end subroutine check_nn_atomic_energies

    subroutine angular_move(atom_index)

      implicit none

      ! imported variables
      
      ! input 
      integer, intent(in) :: atom_index

      ! local variables
      type (vector) :: center
      real*8 :: inv_n_atoms
      real*8 :: distance
      real*8 :: max_distance
      real*8 :: theta
      real*8 :: phi

      ! counter
      integer :: i_atom

      write (150,*) "do angular move for atom # ", atom_index
      
      ! get center of mass
      inv_n_atoms = 1. / n_atoms
      center%x = 0; center%y = 0; center%z = 0;
      
      do i_atom = 1, n_atoms, 1
         center%x = center%x + try_coords%coords(i_atom)%x
         center%y = center%y + try_coords%coords(i_atom)%y
         center%z = center%z + try_coords%coords(i_atom)%z
      end do
      center%x = inv_n_atoms * center%x
      center%y = inv_n_atoms * center%y
      center%z = inv_n_atoms * center%z

      ! get largest distance in cluster
      max_distance = 0.d0
      do i_atom = 1, n_atoms, 1
         if (i_atom .ne. atom_index) then
            distance = norm(try_coords%coords(i_atom) - center)
            if (distance .gt. max_distance) then
               max_distance = distance
            end if
         end if
      end do

      ! determine angles randomly
      theta = pi * arch_rand()
      phi   = 2 * pi * arch_rand()

      try_coords%coords(atom_index)%x = max_distance * sin(theta) * cos(phi)
      try_coords%coords(atom_index)%y = max_distance * sin(theta) * sin(phi)
      try_coords%coords(atom_index)%z = max_distance * cos(theta)

    end subroutine angular_move

  end module basin_hopping
  
