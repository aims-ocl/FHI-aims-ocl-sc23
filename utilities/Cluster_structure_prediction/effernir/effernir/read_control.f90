!  Subroutine read_control provides all general parameters for run

subroutine read_control()
  
  use basin_hopping, bh_energy_threshold_theta => energy_threshold_theta, bh_energy_threshold_boltzmann => energy_threshold_boltzmann, &
       bh_cavity_radius => cavity_radius, bh_cavity => cavity, bh_reset => reset, bh_max_loop => max_loop, bh_minimizer => minimizer, &
       bh_max_move => max_move, bh_kind_of_move => kind_of_move, bh_max_move_increase => max_move_increase, &
       bh_max_move_decrease => max_move_decrease, bh_move_adaptive => move_adaptive, &
       bh_move_distribution => move_distribution, bh_max_distance => max_distance
  
  use cluster
  use control
  use geometry_class
  use conjugate_gradient
  use pair_potential
  
  implicit none
  
  include "constants.f90"
  
  ! local variables
  
  ! eof    : set when end of input is reached
  ! iostat : read statement error condition
  ! desc_str: line descriptor
  
  
  logical :: eof
  integer :: i_code
  
  character*30 :: desc_str
  
  logical :: flag_simulation_flag
  logical :: flag_bh_max_loop
  logical :: flag_bh_energy_threshold_theta
  logical :: flag_bh_energy_threshold_boltzmann
  logical :: flag_bh_reset
  logical :: flag_bh_cavity_radius
  logical :: flag_diff_tolerance
  logical :: flag_hard_sphere_radius
  logical :: flag_bh_max_move
  logical :: flag_bh_max_move_decrease
  logical :: flag_bh_max_move_increase
  logical :: flag_bh_minimizer
  logical :: flag_diff_norm
  logical :: flag_bh_kind_of_move
  logical :: flag_bh_move_adaptive
  logical :: flag_bh_move_distribution
  logical :: flag_bh_max_distance
  logical :: flag_n_next_neighbours
  logical :: flag_r_cut_off
  logical :: flag_energy_interval
  logical :: flag_max_successful_rate
  logical :: flag_min_successful_rate
  logical :: flag_number_of_averaging
  logical :: flag_acceptance_reference
  logical :: flag_switched_E_theta 
  logical :: flag_switched_E_boltzmann
  logical :: flag_spin_step_width
  logical :: flag_energy_tolerance
  logical :: flag_verbose

  logical :: found
  character*30 :: species_temp
  character*30 :: species_temp_2
  real*8 :: species_z_temp
  integer :: n_data_temp
  real*8 :: d_min_temp
  real*8 :: delta_d_temp
  real*8 :: energy_temp

  ! counter
  integer :: i_species
  integer :: i_species_2
  integer :: i_data
  integer :: i_pair_pot
  integer :: i_counter
  integer :: i_atomic_energy
  integer :: i_adsorp

  ! begin work

  if (potential_flag .eq. 'external') then

     write(6,*)
     write(6,'(A)') &
          "------------------------------------------------------------"
     write(6,'(10X, A)') "Reading file control.in"
     write(6,'(A)') &
          "------------------------------------------------------------"
     open(7, file = "control.in", status='OLD', iostat = i_code)
     
     eof = .false.
     
     i_species = 0
     
     if (i_code.ne.0) then
        write(*,*) "* Input file control.in not found."
        stop
     end if
     
     do while (.not. eof)
        read(7,*,iostat = i_code) desc_str
        if (i_code.ne.0) then
           eof = .true.
        else if (desc_str.eq."species") then
           backspace(7)
           i_species = i_species + 1
           read(7,*) desc_str, species_name(i_species)
        end if
     end do
     
     close(7)
     
  end if
  
  write(6,*)
  write(6,'(A)') &
       "------------------------------------------------------------"
  write(6,'(10X, A)') "Reading file control.in.opt"
  write(6,'(A)') &
       "------------------------------------------------------------"
  
  ! initialize
  
  eof = .false.      

  ! provide reasonable defaults
  
  flag_simulation_flag      = .false.
  flag_bh_max_loop          = .false.
  flag_bh_reset             = .false.
  flag_bh_cavity_radius     = .false.
  flag_diff_tolerance       = .false.
  flag_hard_sphere_radius   = .false.
  flag_bh_max_move          = .false.
  flag_bh_max_move_increase = .false.
  flag_bh_max_move_decrease = .false.
  flag_bh_minimizer         = .false.
  flag_diff_norm            = .false.
  flag_bh_kind_of_move      = .false.
  flag_bh_move_adaptive     = .false.
  flag_bh_move_distribution = .false.
  flag_bh_max_distance      = .false.
  flag_bh_energy_threshold_theta = .false.
  flag_bh_energy_threshold_boltzmann = .false.
  flag_n_next_neighbours = .false.
  flag_r_cut_off = .false.
  flag_energy_interval = .false.
  flag_max_successful_rate = .false.
  flag_min_successful_rate = .false.
  flag_number_of_averaging = .false.
  flag_acceptance_reference = .false.
  flag_switched_E_theta  = .false.
  flag_switched_E_boltzmann = .false.
  flag_spin_step_width = .false.
  flag_energy_tolerance = .false.
  flag_verbose = .false.

  target_number_of_minima = -1
  optimize_energy = .false.
  adsorption_move = .false.
  pure_adsorption_move = .false.
  use_rotation_move = .false.
  use_angular_move = .false.
  fraction_of_rotation_move = 1
  scf_loop_to_switch = -1
  spin_polarized = .false.
  charged = .false.
  bh_flavour = 'standard'

  ! open file
  
  i_species = 0
  i_pair_pot = 0
  i_atomic_energy = 0
  i_adsorp = 0

  open (7, FILE="control.in.opt")
  
  ! read first line
  
  read (7, *, iostat = i_code) desc_str
  
  if (i_code.ne.0) then
     write (*,*) "Empty input file control.in.opt "
     eof = .true. 
  end if
  
  do while (.not.eof)
     
     ! decode contents of current line
     
     if (desc_str(1:1) .eq. '#') then
        read (7, *, iostat = i_code) desc_str
        
        if (i_code.ne.0) then
           write(6,*)
           write(6,'(2X,A)') "Input file control.in.opt ends. "
           write(6,'(A)') & 
                "------------------------------------------------------------"
           eof = .true.
        end if
        cycle
     end if
     
     select case(desc_str)
        
     case ('simulation')
        
        backspace(7)
        
        read (7,*) desc_str, simulation_flag
        write (6,*) simulation_flag, " chosen as simulation method."
        
        if (simulation_flag .eq. 'sample-hop-matrix') then

           backspace(7)

           read (7,*) desc_str, desc_str, n_samples, n_starting_isomers
           allocate(starting_index(n_starting_isomers))
           backspace(7)
           read (7,*) (desc_str, i_counter = 1, 4), (starting_index(i_counter), i_counter=1, n_starting_isomers)
           
           write (6,*) "converge hopping matrix up to ", n_samples, " samples"
           write (6,*) "starting indizes:"
           do i_counter = 1, n_starting_isomers, 1
              write (6,'(I4,1X)',advance='no') starting_index(i_counter)
           end do
           write (6,*) 

        end if

        flag_simulation_flag = .true.

     case ('spin_polarized')
        
        backspace(7)
        
        read (7,*) desc_str, spin_polarized
        write (6,*) "Spin polarized calculations? ", spin_polarized

     case ('charged')
        
        backspace(7)
        
        read (7,*) desc_str, charged
        write (6,*) "Charged system? ", charged
      
     case ('spin_step_width')

        backspace(7)
        
        read (7,*) desc_str, spin_step_width
        write (6,*) "Spin per atom is randomly varied by a maximum amount of ", spin_step_width, " electrons."

        flag_spin_step_width = .true.

     case ('energy_tolerance') 

        backspace(7)
        
        read (7,*) desc_str, energy_tolerance
        write (6,*) "Energy tolerance to distinguish structures ", energy_tolerance, " eV/atom. "

        flag_energy_tolerance = .true.

     case ('potential')
        
        continue

     case ('diff_norm')
        
        backspace(7)
        
        read (7,*) desc_str, diff_norm
        select case (diff_norm)

        case ('average')
           write (6,*) "Average norm for measuring difference in structures chosen."
        case ('maximum')
           write (6,*) "Maximum norm for measuring difference in structures chosen."
        case default
           write (6,*) "Chosen norm not implemented."
           write (6,*) "* Aborting."
           stop
        end select
        
        flag_diff_norm = .true.
        
     case ('species')
        
        i_species = i_species + 1

        backspace(7)
        
        select case (potential_flag)

        case ('LJ')
           read (7,*) desc_str, species_name(i_species), LJ_epsilon(i_species), LJ_sigma(i_species)

        case ('NN')
           read (7,*) desc_str, species_name(i_species)

        case ('external')
           continue

        case default
           write (6,*) "Chosen interaction not implemented."
           write (6,*) "* Aborting."
           stop

        end select

     case ('atomic_energy')
        
        backspace(7)
        
        read (7,*) desc_str, species_temp, energy_temp, species_z_temp

        ! check whether we know this species
        found = .false.
        do i_species = 1, n_species, 1
           if (species_temp .eq. species_name(i_species)) then
              i_atomic_energy = i_atomic_energy  + 1
              found = .true.
              exit
           end if
        end do
        if (.not.found) then
           write(6,'(1X,A,A,A,A)') &
                "* Species ", species_temp, ", listed in ", &
                "control.in.opt, is not described in input file control.in."
           write(6,*) "* List of ", n_species, " species: "
           do i_species = 1, n_species, 1
              write(*,*) "* ",i_species,":",species_name(i_species)
           enddo
           write(6,*) "* Aborting."
           stop
        end if
        atomic_energy(i_species) = energy_temp
        species_z(i_species) = species_z_temp

     case ('target_energy')

        backspace(7)
            
        read (7,*) desc_str, target_energy
        
        write (6,'(2X, A, F10.4, A)') "target_energy for optimization: ", target_energy, " eV."
        
        optimize_energy = .true.

     case ('adsorption_atom')
        
        i_adsorp = i_adsorp + 1

        backspace(7)
        
        read (7,*) desc_str, adsorption_atom(i_adsorp), ad_radius_reduction(i_adsorp), pure_adsorption_move
        
        write (6,'(2X, A, I4, 1X, I4, A)') "index for adsorption atom # ", i_adsorp, adsorption_atom(i_adsorp), " ."
        write (6,'(2X, A, F8.4, A)') "reduction factor for radius: ", ad_radius_reduction(i_adsorp), " ."
        write (6,'(2X, A, L4, A)') "additional displacement moves?", .not.(pure_adsorption_move), " ."

        adsorption_move = .true.

     case ('energy_interval')

        backspace(7)
            
        read (7,*) desc_str, energy_interval
        
        write (6,'(2X, A, F10.4, A)') "energy_interval for analysis of sampling-efficiency: ", energy_interval, " eV."
        
        flag_energy_interval = .true.

     case ('target_number_of_minima')
        
        backspace(7)
        
        read (7,*) desc_str, target_number_of_minima 
        
        write (6,'(2X, A, I4)') "aimed number of minima within energy_interval above groundstate: ", target_number_of_minima
        
     case ('pair_pot')

        backspace(7)
        
        select case (potential_flag)

        case('external')
           
           read (7,*) desc_str, species_temp, species_temp_2, n_data_temp, d_min_temp, delta_d_temp
           ! check whether we know this species
           found = .false.
           do i_species = 1, n_species, 1
              if (species_temp .eq. species_name(i_species)) then
                 found = .true.
                 exit
              end if
           end do
           if (found) then
              found = .false.
              do i_species_2 = 1, n_species, 1
                 if (species_temp_2 .eq. species_name(i_species_2)) then
                    i_pair_pot = i_pair_pot + 1
                    found = .true.
                    exit
                 end if
              end do
           end if
           if (.not.found) then
              write(6,'(1X,A,A,A,A)') &
                   "* Species ", species_temp, ", listed in ", &
                   "control.in.opt, is not described in input file control.in."
              write(6,*) "* List of ", n_species, " species: "
              do i_species = 1, n_species, 1
                 write(*,*) "* ",i_species,":",species_name(i_species)
              enddo
              write(6,*) "* Aborting."
              stop
           end if

           n_data(i_species, i_species_2) = n_data_temp
           d_min(i_species, i_species_2) = d_min_temp
           delta_d(i_species, i_species_2) = delta_d_temp

           backspace(7)

           read (7,*) (desc_str, i_counter = 1, 6), (pot_par(i_data, i_species, i_species_2), i_data= 1, n_data(i_species, i_species_2))
           write (6,'(1X,A,A5,A3,A5,A,F4.2,A,F4.2)') "pair-potential[eV]: ", species_name(i_species), " - ", species_name(i_species_2), &
                " d_min[Ang]=", d_min(i_species, i_species_2), " delta_d[Ang]=", delta_d(i_species, i_species_2)
           write (6,'(E10.4, 1X)') (pot_par(i_data, i_species, i_species_2), i_data= 1, n_data(i_species, i_species_2))
   
           n_data(i_species_2, i_species) =  n_data(i_species, i_species_2)
           d_min(i_species_2, i_species) =  d_min(i_species, i_species_2)
           delta_d(i_species_2, i_species) = delta_d(i_species, i_species_2)
           pot_par(:, i_species_2, i_species) = pot_par(:, i_species, i_species_2)

        case default
           write (6,*) "For non-external potentials, pair-potential parameters are not used."
           continue

        end select

     case ('bh_max_loop')
        
        backspace(7)
        
        read (7,*) desc_str, bh_max_loop
        write (6,'(2X,A,I10)') &
             "maximum number of bh_loops : ", bh_max_loop
        
        flag_bh_max_loop = .true.  
        
        if (bh_max_loop .le. 0) then
           write (6,*) "What do you mean with a negative number of max_loops?"
           write (6,*) "* Aborting"
           stop
        end if
        
     case ('bh_energy_threshold_theta')
        
        backspace(7)
        
        read (7,*) desc_str, bh_energy_threshold_theta
        write (6,'(2X,A,F10.4)') &
             "Energy threshold theta for basin hopping [eV]: ", &
             bh_energy_threshold_theta
        
        flag_bh_energy_threshold_theta = .true.

     case ('bh_energy_threshold_boltzmann')
        
        backspace(7)
        
        read (7,*) desc_str, bh_energy_threshold_boltzmann
        write (6,'(2X,A,F10.4)') &
             "Energy threshold boltzmann for basin hopping [eV]: ", &
             bh_energy_threshold_boltzmann
        
        flag_bh_energy_threshold_boltzmann = .true. 

     case ('bh_reset')
        
        backspace(7)
        
        read (7,*) desc_str, bh_reset
        
        if (bh_reset) then
           write (6,*) "After each local relaxation, coordinates are reset to beginning of relaxation."
        else
           write (6,*) "After each local relaxation, coordinates of local minimum are kept."
        end if
        
        flag_bh_reset = .true.
        
     case ('bh_cavity_radius')
        
        backspace(7)
        
        read (7,*) desc_str, bh_cavity_radius
        
        write (6,'(2X,A,F10.4)') &
             "cavity_radius [Ang]: ", &
             bh_cavity_radius
        
        flag_bh_cavity_radius = .true.
        
        bh_cavity = .true.
        
        if (bh_cavity_radius .le. 0.d0) then
           write (6,*) "Why so complicated? Let's just take a positive cavity_radius."
           write (6,*) "Set bh_cavity_radius to - bh_cavity_radius"
           bh_cavity_radius = - bh_cavity_radius
        end if

     case ('bh_kind_of_move')

        backspace(7)
        
        read (7,*) desc_str, bh_kind_of_move
        
        select case (bh_kind_of_move)

        case ('single')
           write (6,'(2X,A)') "Single atomic moves chosen for sampling."
        case ('worst-single')
           write (6,'(2X,A)') "Single atomic moves of worst atom chosen for sampling."  

           select case (potential_flag)

           case ('external')
              use_pair_pot = .true.

           case default
              continue

           end select

        case ('collective')
           write (6,'(2X,A)') "Collective atomic moves chosen for sampling."

        case ('collective-vista')
           write (6,'(2X,A)') "Collective-vista atomic moves chosen for sampling."

        case default
           write (6,*) "Chosen move is not implemented."
           write (6,*) "* Aborting."
           stop
        end select
        
        flag_bh_kind_of_move = .true.
        
     case ('bh_max_move')
        
        backspace(7)
        
        read (7,*) desc_str, bh_max_move
        
        write (6,'(2X,A,F10.4)') &
             "max_move [Ang]: ", &
             bh_max_move
        
        flag_bh_max_move = .true.
        
        if (bh_max_move .le. 0.d0) then
           write (6,*) "Why so complicated? Let's just take a positive max_move."
           write (6,*) "Set bh_max_move to - bh_max_move"
           bh_max_move = - bh_max_move
        end if

     case ('bh_max_distance')
        
        backspace(7)
        
        read (7,*) desc_str, bh_max_distance
        
        write (6,'(2X,A,F10.4)') &
             "max_distance [Ang]: ", &
             bh_max_distance
        
        flag_bh_max_distance = .true.
        
        if (bh_max_distance .le. 0.d0) then
           write (6,*) "Negative distance doesn't make sense."
           write (6,*) "* Aborting."
           stop
        end if

     case ('rotation_move')

        backspace(7)
        
        read (7,*) desc_str, fraction_of_rotation_move
        
        write (6,'(2X,A,I4,A)') &
             "Perform rotation moves every : ",  fraction_of_rotation_move, &
             " . trial move."
        
        use_rotation_move = .true.

     case ('angular_move')

        backspace(7)
        
        read (7,*) desc_str, energy_ratio_ang_move
        
        write (6,'(2X,A,F6.2)') &
             "Perform angular moves for atoms with a ratio high/low atomic energy of : ", energy_ratio_ang_move 
        
        use_angular_move = .true.

     case ('bh_move_distribution')
        
        backspace(7)
        
        read (7,*) desc_str, bh_move_distribution
        
        select case (bh_move_distribution)

        case ('uniform')
           write (6,'(2X,A)') "Uniform distribution of moves chosen."

        case ('poisson')
           write (6,'(2X,A)') "Poisson distribution of moves chosen."

        case default
           write (6,*) "Chosen distribution of moves not implemented."
           write (6,*) "* Aborting."
           stop

        end select

        flag_bh_move_distribution = .true.
  
     case ('adapting_scheme')
        
        backspace(7)
        
        read (7,*) desc_str, adapting_scheme
        
        select case (adapting_scheme)

        case ('successful_additive')
           write (6,'(2X,A)') "Additively adapt moves according to successful rate."

        case ('maximum_rock')
           write (6,'(2X,A)') "Adapt moves according to maximizing successful moves (so that it rocks)."

           backspace(7)
           read (7,*) desc_str, desc_str, number_of_coarse_averaging, number_of_averaging, delta_max_move, kick_in_the_ass, max_max_move, &
                min_delta_max_move, max_delta_max_move,  alpha_unsucc_max, alpha_highE_max, alpha_tolerance, max_move_modifier_one, max_move_modifier_two

           write (6,'(3X,A,I4)') "number_of_coarse_averaging ", number_of_coarse_averaging
           write (6,'(3X,A,I4)') "number_of_averaging ", number_of_averaging
           write (6,'(3X,A,F8.4)') "delta_max_move ", delta_max_move
           write (6,'(3X,A,F8.4)') "kick_in_the_ass ", kick_in_the_ass
           write (6,'(3X,A,F8.4)') "max_max_move ", max_max_move
           write (6,'(3X,A,F8.4)') "min_delta_max_move ", min_delta_max_move
           write (6,'(3X,A,F8.4)') "max_delta_max_move ", max_delta_max_move
           write (6,'(3X,A,F8.4)') "alpha_unsucc_max ", alpha_unsucc_max
           write (6,'(3X,A,F8.4)') "alpha_highE_max ", alpha_highE_max
           write (6,'(3X,A,F8.4)') "alpha_tolerance ", alpha_tolerance
           write (6,'(3X,A,F8.4)') "max_move_modifier_one ", max_move_modifier_one
           write (6,'(3X,A,F8.4)') "max_move_modifier_two ", max_move_modifier_two
           write (6,'(2X,A)') "-------------------------------"

           flag_number_of_averaging = .true.

        case default
           write (6,*) "Chosen adapting scheme not implemented."
           write (6,*) "* Aborting."
           stop

        end select

        flag_bh_move_adaptive = .true.
        bh_move_adaptive = .true.

     case ('bh-flavour')
        
        backspace(7)
        
        read (7,*) desc_str, bh_flavour
        
        select case (bh_flavour)

        case ('standard')
           write (6,'(2X,A)') "Standard serial BH-run chosen."

        case ('histogram-weighted')
           write (6,'(2X,A)') "Choose starting basin according to histogram weighted statistics."

        case default
           write (6,*) "Chosen BH-flavour not implemented."
           write (6,*) "* Aborting."
           stop

        end select

     case ('acceptance_reference')
        
        backspace(7)
        
        read (7,*) desc_str, acceptance_reference
        
        select case (acceptance_reference)

        case ('last_energy')

           write (6,'(2X,A)') "Evaluate acceptance of new isomer according to last energy."
           write (6,'(2X,A)') "So get a markovian chain."
                      
        case ('current_groundstate')

           backspace(7)

           read (7,*) desc_str, desc_str, scf_loop_to_switch

           if (scf_loop_to_switch .gt. 0) then
              ! first start with last-energy
              acceptance_reference = 'last_energy'
           end if

           write (6,'(2X,A)') "Evaluate acceptance of new isomer according to energy difference to current groundstate."
           write (6,'(2X,A,I4)') "Switching after scf-loop # ", scf_loop_to_switch

        case default

           write (6,*) "Acceptance reference not supported."
           write (6,*) "* Aborting."
           stop
           
      end select

        flag_acceptance_reference = .true.

     case ('switched_E_theta')

        backspace(7)
        
        read (7,*) desc_str, switched_E_theta
        write (6,'(2X,A,F8.4,A)') "Switching E_theta to ", switched_E_theta, " eV."

        flag_switched_E_theta = .true.

     case ('switched_E_boltzmann')

        backspace(7)
        
        read (7,*) desc_str, switched_E_boltzmann
        write (6,'(2X,A,F8.4,A)') "Switching E_boltzmann to ", switched_E_boltzmann, " eV."

        flag_switched_E_boltzmann = .true.

     case ('max_successful_rate')
        
        backspace(7)
        
        read (7,*) desc_str, max_successful_rate
        
        write (6,'(2X,A,F10.4)') &
             "maximum successful rate: ", &
             max_successful_rate 
        
        flag_max_successful_rate = .true.
        
        if ((max_successful_rate .gt. 1.d0) .or. (max_successful_rate .lt. 0.d0)) then
           write (6,*) "Successful rate must lie between zero and one."
           write (6,*) "* Aborting."
           stop
        end if
        
     case ('min_successful_rate')
        
        backspace(7)
        
        read (7,*) desc_str, min_successful_rate
        
        write (6,'(2X,A,F10.4)') &
             "minimum successful rate: ", &
             min_successful_rate 
        
        flag_min_successful_rate = .true.
        
        if ((min_successful_rate .gt. 1.d0) .or. (min_successful_rate .lt. 0.d0)) then
           write (6,*) "Successful rate must lie between zero and one."
           write (6,*) "* Aborting."
           stop
        end if

     case ('number_of_averaging')
        
        backspace(7)
        
        read (7,*) desc_str, number_of_averaging
        
        write (6,'(2X,A,I6)') "number of moves for averaging: ", number_of_averaging

        flag_number_of_averaging = .true.

     case ('bh_max_move_increase')
        
        backspace(7)
        
        read (7,*) desc_str, bh_max_move_increase
        
        write (6,'(2X,A,F10.4)') &
             "increase of max_move [Ang]: ", &
             bh_max_move_increase
        
        flag_bh_max_move_increase = .true.
        
        if (bh_max_move_increase .le. 0.d0) then
           write (6,*) "Negative value for an increase does not make sense."
           write (6,*) "* Aborting."
           stop
        end if

     case ('bh_max_move_decrease')
        
        backspace(7)
        
        read (7,*) desc_str, bh_max_move_decrease
        
        write (6,'(2X,A,F10.4)') &
             "decrease of max_move [Ang]: ", &
             bh_max_move_decrease
        
        flag_bh_max_move_decrease = .true.
        
        if (bh_max_move_decrease .le. 0.d0) then
           write (6,*) "Negative value for a decrease does not make sense."
           write (6,*) "* Aborting."
           stop
        end if

     case ('bh_minimizer')
        
        backspace(7)
        
        read (7,*) desc_str, bh_minimizer
        
        flag_bh_minimizer = .true.

        select case (bh_minimizer)

        case ('external')
           write (6,'(2X,A)') "MC-simulation coupled to external code."

        case ('cg')
           write (6,'(2X,A)') "Using internal conjugate-gradient scheme."
           write (6,'(2X,A)') "Currently only in connection with LJ-potential."
           potential_flag = 'LJ'

           backspace(7)
           
           read (7,*) desc_str, desc_str, max_force_allowed, max_pr_force, max_cosine, trial_step_initial, max_line_ts, max_increase_ts, min_progress, max_progress, &
                step_increase, step_decrease, sw_prerelax, max_pr_ts, damping_factor, pr_step_decrease, max_n_cg_loops, max_line_steps

           write (6,'(4X,A)') "Using the following parameters for cg:"
           write (6,'(4X,A,E10.4)') "max_force_allowed [eV/A]: ", max_force_allowed
           write (6,'(4X,A,E10.4)') "max_pr_force [eV/A]     : ", max_pr_force
           write (6,'(4X,A,F8.4)')  "max_cosine              : ", max_cosine
           write (6,'(4X,A,E10.4)') "trial_step              : ", trial_step_initial
           write (6,'(4X,A,E10.4)') "max_line_ts             : ", max_line_ts
           write (6,'(4X,A,E10.4)') "max_increase_ts         : ", max_increase_ts
           write (6,'(4X,A,E10.4)') "min_progress            : ", min_progress
           write (6,'(4X,A,E10.4)') "max_progress            : ", max_progress
           write (6,'(4X,A,E10.4)') "step_increase           : ", step_increase
           write (6,'(4X,A,E10.4)') "step_decrease           : ", step_decrease
           write (6,'(4X,A,E10.4)') "sw_prerelax             : ", sw_prerelax
           write (6,'(4X,A,E10.4)') "max_pr_ts               : ", max_pr_ts
           write (6,'(4X,A,E10.4)') "damping_factor          : ", damping_factor
           write (6,'(4X,A,E10.4)') "pr_step_decrease        : ", pr_step_decrease
           write (6,'(4X,A,I6)')    "max_n_cg_loops          : ", max_n_cg_loops
           write (6,'(4X,A,I6)')    "max_line_steps          : ", max_line_steps
           write (6,*) 

        case ('magic')
           potential_flag = 'NN'
           write (6,'(2X,A)') "Using neural network based optimization by J.Behler for obtaining optimized structures."           

        case default
           write (6,*) "Chosen minimizer not implemented."
           write (6,*) "* Aborting."
           stop

        end select

     case ('diff_tolerance')
        
        backspace(7)
        
        read (7,*) desc_str, diff_tolerance
        
        write (6,'(2X, A, E10.4, A)') &
             "distance tolerance for equivalence of structures ", diff_tolerance, " ."
        
        flag_diff_tolerance = .true.  
        
        if (diff_tolerance .le. 0) then
           write (6,*) "Negative number does not make sense."
           write (6,*) "* Aborting"
           stop
        end if

     case ('hard_sphere_radius')
        
        backspace(7)
        
        read (7,*) desc_str, hard_sphere_radius
        hard_sphere_radius = hard_sphere_radius
        
        write (6,'(2X, A, E10.4, A)') &
             "hard sphere radius for pre-rejecting a move = ", hard_sphere_radius, " Ang."
        
        flag_hard_sphere_radius = .true.  
        
        if (hard_sphere_radius .lt. 0.d0) then
           write (6,*) "Negative number does not make sense."
           write (6,*) "* Aborting"
           stop
        end if

     case ('n_next_neighbours')
        
        backspace(7)
        
        read (7,*) desc_str, n_next_neighbours
        
        write (6,'(2X, A, I4, A)') &
             " # neighbours to be taken into account for angular distribution function = ", &
             n_next_neighbours, " ."
        
        flag_n_next_neighbours = .true.  
        
        if (n_next_neighbours .lt. 2) then
           write (6,*) "How do you want to define an angle with less than 3 atoms ?"
           write (6,*) "* Aborting"
           stop
        end if
        if (n_next_neighbours .gt. (n_atoms - 1)) then
           write (6,*) "Sorry, you just have ", n_atoms - 1, " atoms as neighbours."
           write (6,*) "* Aborting"
           stop
        end if

     case ('r_cut_off')
        
        backspace(7)
        
        read (7,*) desc_str, r_cut_off
        
        write (6,'(2X, A, I4, A)') &
             " cut of radius for order parameter O_4 = ", &
             r_cut_off, " ."
        
        flag_r_cut_off = .true.  
        
     case ('verbose')

        backspace(7)
        
        read (7,*) desc_str, verbose

        write (6,'(2X, A, L, A)') &
             " verbosity = ", &
             verbose, " ."
        
        flag_verbose = .true.  

     case default
        
        write (*,*) "Unknown descriptor ", desc_str, &
             " in file control.in."
        stop
        
     end select
     
     ! read next line
     
     read (7,*,iostat = i_code) desc_str
     
     if (i_code.ne.0) then
        write(6,*)
        write (6,'(2X,A)') "Input file control.in.opt ends. "
        write(6,'(A)') & 
             "------------------------------------------------------------"
        eof = .true.
     end if
     
  end do
  
  ! close control.in
  
  close(7)
  
  ! check whether relevant quantities have been set - if not, use defaults:
  
  if (.not. flag_simulation_flag) then
     simulation_flag = 'basin-hopping'
     write (6,'(2X,A)') "No simulation method selected. "
     write (6,'(2X,A)') "Make a basin-hopping simulation."
  end if
  if (.not. flag_diff_tolerance) then
     diff_tolerance = 1d-4
     write (6,'(2X,A)') "diff_tolerance not set."
     write (6,'(2X,A,F8.4)') "Using default: diff_tolerance = ", diff_tolerance
  end if
  if (.not. flag_hard_sphere_radius) then
     hard_sphere_radius = 0.0
     write (6,'(2X,A)') "hard sphere radius not set."
     write (6,'(2X,A,F8.4,A)') "Using default: hard sphere radius = ", hard_sphere_radius, " Ang."
  end if
  if (.not. flag_diff_norm) then
     diff_norm = 'average'
     write (6,'(2X,A)') "Average norm for measuring difference in structures chosen."
  end if
  if (.not. flag_n_next_neighbours) then
     n_next_neighbours = 2
     write (6,'(2X,A)') "n_next_neighbours not set."
     write (6,'(2X,A)') "Using default: n_next_neighbours = 2"
  end if
  if (.not. flag_r_cut_off) then
     r_cut_off = 1.391 ! see Doye et al. , JCP 110(14), 6896 (1999)
     write (6,'(2X,A)') "cut off radius for order parameter O_4 not set."
     write (6,'(2X,A,F8.4,A)') "Using default: ", r_cut_off, " Ang."
  end if
  if (.not.flag_energy_interval) then
     energy_interval = 1.0d0
     write (6,'(2X,A)') "energy interval for sampling-efficiency analysis not set."
     write (6,'(2X,A,F8.4,A)') "Using default: ", energy_interval, " eV."
  end if
  if (.not.flag_verbose) then
     verbose = .false.
     write (6,'(2X,A)') "verbosity not specified."
     write (6,'(2X,A,L,A)') "defaulting to ", verbose, " ."
  end if
  if ((simulation_flag .eq. 'basin-hopping') .or. (simulation_flag .eq. 'sample-hop-matrix')) then
     if (.not. flag_bh_energy_threshold_theta) then
        bh_energy_threshold_theta = 0.25e0
        write(6,'(2X,A)') &
             "Energy threshold theta for basin hopping has not been set."
        write(6,'(2X,A,F10.4,A)') &
             "Energy threshold theta = ", bh_energy_threshold_theta, " eV."
     end if
     if (.not. flag_bh_energy_threshold_boltzmann) then
        bh_energy_threshold_boltzmann = 0.25e0
        write(6,'(2X,A)') &
             "Energy threshold boltzmann for basin hopping has not been set."
        write(6,'(2X,A,F10.4,A)') &
             "Energy threshold boltzmann = ", bh_energy_threshold_boltzmann, " eV."
     end if
     if (.not. flag_bh_max_loop) then
        bh_max_loop = 100
        write (6,'(2X,A)') &
             "Maximum number of bh-loops not chosen."
        write (6,'(2X,A,I4)') &
             "Defaulting to:", bh_max_loop
     end if
     if (.not. flag_bh_reset) then
        bh_reset = .false.
        write (6,'(2X,A)') &
             "After each local relaxation, coordinates of local minimum are kept."
     end if
     if (.not. flag_bh_move_adaptive) then
        bh_move_adaptive = .false.
        write (6,'(2X,A)') &
             "No history feedback method for moves."
     end if
     if (.not.flag_bh_move_distribution) then
        bh_move_distribution = 'poisson'
        write (6,'(2X,A)') &
             "No move distribution chosen. Poisson distribution is assumed."
     end if
     if (.not.flag_bh_max_distance) then
        bh_max_distance = 3.d0
        write (6,'(2X,A)') "maximum distance not set."
        write (6,'(2X,A,F8.4,A)') "Using default: maximum distance = ", bh_max_distance, " Ang."
     end if

     if (.not. flag_bh_minimizer) then

        select case (potential_flag)
           
        case ('LJ')
           write (6,'(2X,A)') "For LJ-potentials, internal cg-optimizer is chosen."
           bh_minimizer = 'cg'
           max_force_allowed = 1d-3
           max_pr_force = 1d-1
           max_cosine = 1d-2
           trial_step = 1d-3
           max_line_ts = 0.1
           max_increase_ts = 0.2
           min_progress = 1.5
           max_progress = 10.
           step_increase = 2.
           step_decrease = 0.5
           sw_prerelax = 1d-2
           max_pr_ts = 0.1
           damping_factor = 0.0
           pr_step_decrease = 0.3
           max_n_cg_loops = 1000
           max_line_steps = 100
           write (6,*) "Using the following parameters for cg:"
           write (6,*) "max_force_allowed [eV/A]:", max_force_allowed
           write (6,*) "max_pr_force [eV/A]:", max_pr_force
           write (6,*) "max_cosine:", max_cosine
           write (6,*) "trial_step:", trial_step
           write (6,*) "max_line_ts:", max_line_ts
           write (6,*) "max_increase_ts:", max_increase_ts
           write (6,*) "min_progress:", min_progress
           write (6,*) "max_progress:", max_progress
           write (6,*) "step_increase:", step_increase
           write (6,*) "step_decrease:", step_decrease
           write (6,*) "sw_prerelax:", sw_prerelax
           write (6,*) "max_pr_ts:", max_pr_ts
           write (6,*) "damping_factor:", damping_factor
           write (6,*) "pr_step_decrease:", pr_step_decrease
           write (6,*) "max_n_cg_loops:", max_n_cg_loops
           write (6,*) "max_line_steps:", max_line_steps
   
        case ('external')
           bh_minimizer = 'external'
           continue
           
        case ('NN')
           bh_minimizer = 'magic'
           continue

        case default
           write (6,'(2X,A)') "Internal inconsistency."
           write (6,'(2X,A)') "(read_control/flag_bh_minimizer)"
           write (6,'(2X,A)') "* Aborting."
           
        end select
     end if
     if (.not. flag_bh_max_move) then
        bh_max_move = 0.5d0
        write (6,'(2X,A,F10.4,A)') &
             "No maximum move chosen. Defaulting to ", bh_max_move, " Ang."
     end if
     if (.not. flag_bh_kind_of_move) then
        bh_kind_of_move = 'single'
        write (6,'(2X,A)') &
             "Single atomic moves are used."
     end if
     if (spin_polarized .and. .not.flag_spin_step_width) then
        spin_step_width = 2.d0
        write (6,'(2X,A,F10.4,A)') "No step width for spin variation chosen. Defaulting to ", spin_step_width, " electrons per atom."
     end if
     if (.not.flag_energy_tolerance) then
        energy_tolerance = 0.01d0
        write (6,'(2X,A,F10.4,A)') "No energy tolerance to distinguish structures chosen. Defaulting to ", energy_tolerance,  " meV/atom ."
     end if

     if ((bh_move_adaptive) .and. (adapting_scheme .eq. 'successful_additive')) then

        if (.not. flag_bh_max_move_decrease) then
           select case (bh_kind_of_move)
              
           case ('single')
              
              bh_max_move_decrease = 0.25d0
              
           case ('collective')
              
              bh_max_move_decrease = 0.1d0
              
           end select
           write (6,'(2X,A,F10.4,A)') &
                "No maximum move decrease for adapting scheme chosen. Defaulting to ", bh_max_move_decrease, " Ang."
        end if

        if (.not. flag_bh_max_move_increase) then
           select case (bh_kind_of_move)
              
           case ('single')
              
              bh_max_move_increase = 0.25d0
              
           case ('collective')
              
              bh_max_move_increase = 0.1d0
              
           end select
           write (6,'(2X,A,F10.4,A)') &
                "No maximum move increase for adapting scheme chosen. Defaulting to ", bh_max_move_increase, " Ang."
        end if

        if (.not. flag_max_successful_rate) then
           max_successful_rate = 0.95d0

           write (6,'(2X,A,F10.4,A)') &
                "No maximum successful rate for adapting scheme chosen. Defaulting to ", max_successful_rate, " ."
        end if

        if (.not. flag_min_successful_rate) then
           min_successful_rate = 0.6d0

           write (6,'(2X,A,F10.4,A)') &
                "No minimum successful rate for adapting scheme chosen. Defaulting to ", min_successful_rate, " ."
        end if
        if ((.not. flag_number_of_averaging) .and. (bh_move_adaptive)) then
        write (6,'(2X,A)') "Number of moves for averaging efficiency indicators required for adaptive scheme."
        write (6,'(2X,A)') "* Abort."
        stop
     end if

     end if
     
     if (.not. flag_number_of_averaging) then
        write (6,'(2X,A)') "Just overall averaging of efficiency indicators."
        number_of_averaging = -1 
     end if
     if (.not. flag_acceptance_reference) then
        write (6,'(2X,A)') "Acceptance reference not chosen. Defaulting to last energy."
        acceptance_reference = 'last_energy'
     end if
     if (scf_loop_to_switch .gt. 0) then
        if (.not.flag_switched_E_theta) then
           switched_E_theta = bh_energy_threshold_theta
        end if
        if (.not.flag_switched_E_boltzmann) then
           switched_E_boltzmann = bh_energy_threshold_boltzmann
        end if
     end if
  end if
  
  if ((potential_flag .eq. 'NN') .and. (bh_minimizer .ne. 'magic')) then
     write (6,'(2X,A)') "neural network PES only in connection with external optimizer."
     write (6,'(2X,A)') "Switching to it."
     bh_minimizer = 'magic'
  end if

  if ((potential_flag .eq. 'NN') .and. (simulation_flag .ne. 'basin-hopping')) then
     write (6,'(2X,A)') "neural network PES only in connection with basin-hopping optimization."
     write (6,'(2X,A)') "Switching to it."
     simulation_flag = 'basin-hopping'
  end if

  if ((use_angular_move) .and. ((potential_flag .ne. 'NN') .or. (simulation_flag .ne. 'basin-hopping'))) then
     write (6,'(2X,A)') "angular moves only in connection with neural network PES and basin-hopping optimization."
     write (6,'(2X,A)') "* Abort."
     stop
  end if

  if ((potential_flag .eq. 'LJ') .and. (bh_minimizer .eq. 'external')) then
     write (6,'(2X,A)') "LJ-potentials only in connection with internal cg-optimizer."
     write (6,'(2X,A)') "Switching to it."
     bh_minimizer = 'cg'
     max_force_allowed = 1d-3
     max_pr_force = 1d-1
     max_cosine = 1d-2
     trial_step = 1d-3
     max_line_ts = 0.1
     max_increase_ts = 0.2
     min_progress = 1.5
     max_progress = 10.
     step_increase = 2.
     step_decrease = 0.5
     sw_prerelax = 1d-2
     max_pr_ts = 0.1
     damping_factor = 0.0
     pr_step_decrease = 0.5
     max_n_cg_loops = 1000
     max_line_steps = 100
     write (6,*) "Using the following parameters for cg:"
     write (6,*) "max_force_allowed [eV/A]:", max_force_allowed
     write (6,*) "max_pr_force [eV/A]:", max_pr_force
     write (6,*) "max_cosine:", max_cosine
     write (6,*) "trial_step:", trial_step
     write (6,*) "max_line_ts:", max_line_ts
     write (6,*) "max_increase_ts:", max_increase_ts
     write (6,*) "min_progress:", min_progress
     write (6,*) "max_progress:", max_progress
     write (6,*) "step_increase:", step_increase
     write (6,*) "step_decrease:", step_decrease
     write (6,*) "sw_prerelax:", sw_prerelax
     write (6,*) "max_pr_ts:", max_pr_ts
     write (6,*) "damping_factor:", damping_factor
     write (6,*) "pr_step_decrease:", pr_step_decrease
  end if
  
  if (use_pair_pot .and. (i_pair_pot .ne. (n_species * (n_species +1) / 2))) then
     write (6,*) "You asked for pair-potentials but didn't provide all"
     write (6,*) "required data."
     write (6,*) "* Aborting."
     stop
  end if
  
  if (i_atomic_energy .ne. n_species) then
     write (6,*) "Number of atomic energies doesn't fit the number of species."
     write (6,*) "* Aborting."
     stop
  end if
  
  return
  
end subroutine read_control
   
