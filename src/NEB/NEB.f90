!    usage: 
!    NEB.x <control file> <jmol input file> <jmol output file>
program NEB
  use mpi_routines
  implicit none
  character*100 :: control_name, jmol_in_name, jmol_out_name, keyword, save_data_name, method, current_method
  character*100 :: inputbuf
  integer :: io_status, n_atoms, n_images, n_coords, i_coords, i_coords2, i_atom, i_buf, i_images
  integer :: imax_image, imax_atom, imax_coord, imax_image_before, imax_atom_before, imax_coord_before
  integer :: transition_image
  real*8  :: spring_constant, start_energy, end_energy, max_atomic_move, min_line_step, line_step
  real*8  :: line_step_reduce, max_coordinate_force, force_convergence_criterion, object_function
  real*8  :: stored_object_function, object_function_tolerance, spring_constant_min, spring_constant_max
  real*8  :: switch_to_CINEB_force, reference_energy, max_force_before, sum_of_energies, energy_barrier
  real*8  :: start_BFGS_force, force_times_search_dir, transition_max_change, transition_RMS_change, transition_max_force
  real*8  :: forcenew_magnitude, forceold_magnitude, force_angle, ddot, max_line_search_angle, max_line_step_increase
  logical :: save_data_exists, ini_relaxation, initial_hessian, eof, object_function_exists
  logical :: use_variable_spring_constants, use_CINEB_method, line_search_on
  ! coordinate arrays: allocatable
  real*8,      allocatable, dimension(:)     :: coords
  real*8,      allocatable, dimension(:)     :: stored_coords
  real*8,      allocatable, dimension(:)     :: total_forces
  real*8,      allocatable, dimension(:)     :: forces_keep
  real*8,      allocatable, dimension(:)     :: stored_forces
  real*8,      allocatable, dimension(:)     :: search_direction
  real*8,      allocatable, dimension(:,:)   :: hessian
  real*8,      allocatable, dimension(:)     :: energies
  real*8,      allocatable, dimension(:)     :: start_coords, end_coords
  real*8,      allocatable, dimension(:)     :: start_forces, end_forces
  real*8,      allocatable, dimension(:)     :: delta_R
  real*8,      allocatable, dimension(:)     :: spring_constants
  real*8,      allocatable, dimension(:,:)   :: centre_of_mass
  real*8,      allocatable, dimension(:,:)   :: net_force  
  real*8,      allocatable, dimension(:,:,:) :: relative_coords
  character*2, allocatable, dimension(:)     :: species_names

  !---------------------------------------------------------------------------------------------------
  ! MPI initialization ... 
  call initialize_mpi()

  if (myid.eq.0) then
     write(use_unit,*) '------------------------------------------------------------'
     write(use_unit,*) '     Entering NEB force update '
     write(use_unit,*) '------------------------------------------------------------'

     !-----------------------------------------------------------------------------------------------------
     ! read arguments to figure out proper control- and movie name files ... 
     call getarg(1,inputbuf)
     read(inputbuf,*) control_name
     call getarg(2,inputbuf)
     read(inputbuf,*) jmol_in_name
     call getarg(3,inputbuf)
     read(inputbuf,*) jmol_out_name
     
     !-----------------------------------------------------------------------------------------------------
     max_atomic_move               = 0.1d0 ! defaults, unless something better is given in the control file 
     min_line_step                 = 0.1d0
     line_step_reduce              = 0.5d0
     force_convergence_criterion   = 1d-2
     object_function_tolerance     = 1d0
     method                        = 'no_method_specified'
     object_function_exists        = .false.
     spring_constant_min           = 0d0
     spring_constant_max           = 0d0
     use_variable_spring_constants = .false.
     spring_constant               = 0d0
     switch_to_CINEB_force         = 0d0
     start_BFGS_force              = 0d0
     line_search_on                = .false.
     use_CINEB_method              = .false. 
     max_line_search_angle         = 20d0
     max_line_step_increase        = -1d0
     open(unit = 7, file = control_name, status = 'unknown')
     io_status = 0
     write(use_unit,*) 'Reading control file ', control_name
     do while (io_status.eq.0)
        read(unit=7,fmt=*,iostat=io_status) keyword
        if (io_status.eq.0) then
           if (keyword(1:1).eq.'#') then
           else if (keyword.eq.'save_data') then
              backspace(7)
              read(7,*) keyword, save_data_name
              write(use_unit,'(2X,2A)') 'Obtaining previously saved data from ', save_data_name
           else if (keyword.eq.'n_atoms') then
              backspace(7)
              read(7,*) keyword, n_atoms
              write(use_unit,'(2X,A,I8)') 'n_atoms            = ', n_atoms
           else if (keyword.eq.'n_images') then
              backspace(7)
              read(7,*) keyword, n_images
              write(use_unit,'(2X,A,I8)') 'n_images           = ', n_images
           else if (keyword.eq.'spring_constant') then
              backspace(7)
              read(7,*) keyword, spring_constant
              write(use_unit,'(2X,A,D14.6)') 'spring_constant = ', spring_constant
           else if (keyword.eq.'use_variable_spring_constants') then
              backspace(7)
              read(7,*) keyword, use_variable_spring_constants
              if (use_variable_spring_constants) then
                 write(use_unit,*) 'Using variable spring constants with a linear dependence on the energy.'
              else
                 write(use_unit,*) 'Explicit request not to use variable spring constant.'
              end if
           else if (keyword.eq.'spring_constant_max') then
              backspace(7) 
              read(7,*) keyword, spring_constant_max
              write(use_unit,*) 'Maximal value for range of spring constants: ', spring_constant_max
           else if (keyword.eq.'spring_constant_min') then
              backspace(7)
              read(7,*) keyword, spring_constant_min
              write(use_unit,*) 'Minimal value for range of spring constants:', spring_constant_min
           else if (keyword.eq.'max_atomic_move') then
              backspace(7)
              read(7,*) keyword, max_atomic_move
              write(use_unit,'(2X,A,D14.6)') 'max_atomic_move = ', max_atomic_move
           else if (keyword.eq.'min_line_step') then
              backspace(7)
              read(7,*) keyword, min_line_step
              write(use_unit,'(2X,A,D14.6)') 'min_line_step   = ', min_line_step
           else if (keyword.eq.'max_line_step_increase') then
              backspace(7)
              read(7,*) keyword, max_line_step_increase
              write(use_unit,'(2X,A,D14.6)') 'max_line_step_increase = ', max_line_step_increase
           else if (keyword.eq.'line_step_reduce') then
              backspace(7)
              read(7,*) keyword, line_step_reduce
              write(use_unit,'(2X,A,D14.6)') 'line_step_reduce   = ', line_step_reduce
           else if (keyword.eq.'force_convergence_criterion') then
              backspace(7)
              read(7,*) keyword, force_convergence_criterion
              ! write what was found!
           else if (keyword.eq.'object_function_tolerance') then
              backspace(7)
              read(7,*) keyword, object_function_tolerance
              ! describe what was found!
           else if (keyword.eq.'method') then
              backspace(7)
              read(7,*) keyword, method
              if (method.eq.'PEB') then
                 write(use_unit,*) 'Using the plain elastic band method. '
                 write(use_unit,*) 'Make sure you have a valid spring constant!'
                 object_function_exists = .true.
              else if (method.eq.'NEB') then
                 write(use_unit,*) 'Using the nudged elastic band.'
              else 
                 write(use_unit,*) '* WARNING: Method ',method,' not known.'
                 write(use_unit,*) '*          Aborting execution.'
                 stop
              end if
           else if (keyword.eq.'switch_to_CI-NEB') then
              backspace(7) 
              read(7,*) keyword, switch_to_CINEB_force
              write(use_unit,*) 'Switching to climbing-image NEB after reaching force convergence of ', switch_to_CINEB_force
           else if (keyword.eq.'start_BFGS') then
              backspace(7)
              read(7,*) keyword, start_BFGS_force
              write(use_unit,*) 'Starting BFGS algorithm after reaching a force convergence of ', start_BFGS_force
           else if (keyword.eq.'max_line_search_angle') then
              backspace(7)
              read(7,*) keyword, max_line_search_angle
              write(use_unit,'(2X,2A,F6.2,A)') 'Continuing line searches only if two successive forces don not ', &
                   'differ by more than ',max_line_search_angle, ' degrees.'
           else if ((keyword.eq.'n_iteration_start')             &
                .or.(keyword.eq.'n_max_iteration')               &
                .or.(keyword.eq.'input_template')                &
                .or.(keyword.eq.'aims_control_file_template')    &
                .or.(keyword.eq.'aims_geometry_file_template')   &
                .or.(keyword.eq.'aims_restart_string')           &
                .or.(keyword.eq.'dft_exe')                       &
                .or.(keyword.eq.'dft_code')                      &
                .or.(keyword.eq.'castep_control_file_template')  &
                .or.(keyword.eq.'castep_geometry_file_template') &
                .or.(keyword.eq.'castep_restart_string')         &
                )then
              ! do nothing as these are keywords specific to the PERL driver and should be ignored in this code
           else
              write(use_unit,*) '* WARNING: keyword ',keyword,' not known'
              write(use_unit,*) '*          Aborting execution.'
              stop
           end if
        end if
     end do
     close(unit = 7)
     if (method.eq.'no_method_specified') then
        write(use_unit,*) '* WARNING: no method flag specified explicitly,'
        write(use_unit,*) '*          defaulting to FEB.'
        method = 'NEB'
     end if
     if (use_variable_spring_constants.and.((spring_constant_max.eq.0d0).or.(spring_constant_min.eq.0d0))) then
        write(use_unit,*) '* WARNING: requested the use of variable spring constants but did not specify a proper range!'
        write(use_unit,*) '* Aborting.'
        stop
     end if
     if ((.not.use_variable_spring_constants).and.(spring_constant.eq.0d0)) then
        write(use_unit,*) '* WARNING: no proper spring constant specified.'
        write(use_unit,*) '*          Aborting.'
        stop
     end if

     !---------------------------------------------------------------------------------------------------
     ! set coordinate data and allocate working arrays
     n_coords = 3*n_images*n_atoms
     allocate(coords(n_coords))          
     allocate(stored_coords(n_coords))   
     allocate(total_forces(n_coords))    
     allocate(forces_keep(n_coords))
     allocate(stored_forces(n_coords))   
     allocate(search_direction(n_coords))
     allocate(hessian(n_coords,n_coords))
     allocate(energies(n_images)) 
     allocate(spring_constants(n_images))
     allocate(species_names(n_atoms))
     allocate(start_coords(3*n_atoms))
     allocate(start_forces(3*n_atoms))
     allocate(end_coords(3*n_atoms))
     allocate(end_forces(3*n_atoms))
     allocate(delta_R(0:n_images))
     allocate(centre_of_mass(3,n_images))
     allocate(net_force(3,n_images))
     allocate(relative_coords(3,n_atoms,n_images))
     
     !---------------------------------------------------------------------------------------------------
     ! read previous BFGS information from binary file, if it exists, and if it makes sense.
     inquire(FILE=save_data_name,EXIST=save_data_exists)
     if (save_data_exists) then
        write(use_unit,*) 'Reading saved Hessian data from file ', save_data_name
        open(unit = 7, file = save_data_name, status = 'unknown', form = 'unformatted')
        do i_coords = 1, n_coords
           read(7) stored_coords(i_coords)
           read(7) stored_forces(i_coords)
           read(7) search_direction(i_coords)
           do i_coords2 = 1, n_coords
              read(7) hessian(i_coords,i_coords2)
           end do
        end do
        read(7) line_step
        read(7) stored_object_function
        read(7) use_CINEB_method
        read(7) line_search_on
        read(7) force_times_search_dir
        close(unit = 7)
        ini_relaxation   = .false.
        initial_hessian  = .false.
     else
        ! initialize BFGS information if file does not exist or data does not make sense. 
        write(use_unit,*) 'Initializing BFGS information'
        stored_coords(:)       = 0d0
        stored_forces(:)       = 0d0
        search_direction(:)    = 0d0
        hessian(:,:)           = 0d0
        stored_object_function = 0d0
        do i_coords = 1, n_coords
           hessian(i_coords,i_coords) = 1d0
        end do
        line_step              = 1d0
        ini_relaxation         = .true.
        initial_hessian        = .true.
     end if
     
     !---------------------------------------------------------------------------------------------------
     ! Read geometry and energy information
     write(use_unit,*)
     write(use_unit,'(2X,2A)') 'Reading geometry information from file ', jmol_in_name
     ! read force, geometry, and energies from jmol-type movie 
     open(unit = 7, file = jmol_in_name, status = 'unknown')
     ! read starting configuration
     read(7,*) i_buf
     read(7,*) keyword, keyword, keyword, keyword, keyword, start_energy, keyword, current_method
     i_coords = 1
     do i_atom = 1, n_atoms 
        read(7,*) species_names(i_atom), (start_coords(i_coords+i_buf), i_buf = 0, 2), (start_forces(i_coords+i_buf), i_buf = 0, 2)
        i_coords = i_coords + 3
     end do
     if (current_method.eq.'CI-NEB') use_CINEB_method = .true.

     ! read n_images images
     i_coords = 1
     do i_images = 1, n_images
        read(7,*) i_buf
        read(7,*) keyword, keyword, keyword, keyword, keyword, energies(i_images)
        do i_atom = 1, n_atoms 
           read(7,*) species_names(i_atom), (coords(i_coords+i_buf), i_buf = 0, 2), (total_forces(i_coords+i_buf), i_buf = 0, 2)
           i_coords = i_coords + 3
        end do
     end do
     ! read final configuration
     read(7,*) i_buf
     read(7,*) keyword, keyword, keyword, keyword, keyword, end_energy
     i_coords = 1
     do i_atom = 1, n_atoms 
        read(7,*) species_names(i_atom), (end_coords(i_coords+i_buf), i_buf = 0, 2), (end_forces(i_coords+i_buf), i_buf = 0, 2)
        i_coords = i_coords + 3
     end do
     close(unit = 7)
  
     !---------------------------------------------------------------------------------------------------
     ! calculate and output difference vectors between the images: before BFGS
     call get_difference_magnitude(n_atoms,n_images,start_coords,coords,end_coords,delta_R)
     write(use_unit,'(2X,A)') 'Magnitudes of image-image vectors before geometry update:'
     write(use_unit,'(2X,A,F10.6)') ' | start    -> image  1: ',delta_R(0)
     do i_images = 1, n_images - 1
        write(use_unit,'(2X,A,I2,A,I2,A,F10.6)') ' | image ',i_images,' -> image ',i_images+1,': ',delta_R(i_images)
     end do
     write(use_unit,'(2X,A,I2,A,F10.6)') ' | image ',n_images,' -> end     : ',delta_R(i_images)
     

     write(use_unit,*) 
     write(use_unit,*) 'Cleaning translational force components from each input iteration.'
     call clean_translational_forces(n_atoms,n_images,total_forces)
     write(use_unit,*) 

!!$    !---------------------------------------------------------------------------------------------------
!!$    ! Clean translational forces from each image ... 
!!$     write(use_unit,*)
!!$     write(use_unit,'(2X,A)') 'Cleaning translational components out of input forces.'
!!$     call get_net_force(n_atoms,n_images,total_forces,net_force)
!!$     write(use_unit,*)
!!$     write(use_unit,'(2X,A)')        'Net force on image BEFORE updating forces, in eV/A '
!!$     write(use_unit,'(2X,A)')        '|     x       y        z'
!!$     do i_images = 1, n_images
!!$        write(use_unit,'(2X,A,3F10.6)') '|',(net_force(i_coords,i_images), i_coords = 1, 3)
!!$     end do
!!$     call clean_translational_forces(n_atoms,n_images,total_forces)
!!$     call get_net_force(n_atoms,n_images,total_forces,net_force)
!!$     write(use_unit,*)
!!$     write(use_unit,'(2X,A)')        'Net force on image BEFORE updating forces, after cleaning, in eV/A '
!!$     write(use_unit,'(2X,A)')        '|     x       y        z'
!!$     do i_images = 1, n_images
!!$        write(use_unit,'(2X,A,3F10.6)') '|',(net_force(i_coords,i_images), i_coords = 1, 3)
!!$     end do
!!$
!!$
!!$     write(use_unit,*)

     !---------------------------------------------------------------------------------------------------
     ! calculate and output NEB 'object function' - the area under the path via trapezoid rule
     ! FIXME: This should technically go into the NEB subroutine if it ends iup being 
     ! used, if only as an input option for particularly straight input paths, technically 
     ! it should not be necessary
     reference_energy = min(start_energy,end_energy)
     if (method.eq.'NEB') then
        ! start -> first image
        object_function = delta_R(0)*(energies(1)+start_energy-2d0*reference_energy)/2d0
        ! intermediate spacings
        do i_images = 2, n_images
           object_function = object_function + delta_R(i_images-1)*(energies(i_images)+energies(i_images-1)-2d0*reference_energy)/2d0
!           sum_of_energies = sum_of_energies + energies(n_images)-reference_energy
        end do
        ! last image -> end
        object_function = object_function + delta_R(n_images)*(end_energy+energies(n_images)-2d0*reference_energy)/2d0
        write(use_unit,*) 'The object function for this iteration is: ', object_function 
     end if
     
     !---------------------------------------------------------------------------------------------------
     ! calculate the sum of DFT energies for reference with the object function contributions
     ! also obtain the highest energy difference to the reference (i.e. the final energy!) and which image this 
     ! corresponds to!
     sum_of_energies  = energies(1)
     energy_barrier   = max(start_energy,end_energy)
     transition_image = 1
     do i_images = 2, n_images
        sum_of_energies = sum_of_energies + energies(i_images)
        if (energy_barrier.lt.energies(i_images)) then
           energy_barrier   = energies(i_images)
           transition_image = i_images
        end if
     end do
     sum_of_energies = sum_of_energies - dble(n_images)*reference_energy
     energy_barrier  = energy_barrier - reference_energy

     !---------------------------------------------------------------------------------------------------
     ! calculate spring constants
     if (use_variable_spring_constants) then
        call obtain_variable_spring_constants(n_images,            &  
                                              spring_constants,    &
                                              start_energy,        &
                                              end_energy,          &
                                              energies,            &
                                              spring_constant_max, &
                                              spring_constant_min)
     else
        spring_constants = spring_constant
     end if

     !---------------------------------------------------------------------------------------------------
     ! determine the maximum force on any image coordinate to compare with value after NEB
     max_force_before = 0d0
     do i_images = 1, n_images
        do i_atom = 1, n_atoms
           do i_coords = 1, 3
              if (max_force_before.lt.abs(total_forces((i_images-1)*3*n_atoms+(i_atom-1)*3+i_coords))) then
                 max_force_before  =  abs(total_forces((i_images-1)*3*n_atoms+(i_atom-1)*3+i_coords))
                 imax_image_before = i_images
                 imax_atom_before  = i_atom
                 imax_coord_before = i_coords
              end if
           end do
        end do
     end do
     
     !---------------------------------------------------------------------------------------------------
     ! calculate method-dependent forces  

     forces_keep(:) = total_forces(:)

     if (method.eq.'PEB') then
        if (.not.use_CINEB_method) then
           !   Pure elastic band
           write(use_unit,*) 'Adjusting forces with pure elastic band method.'
           call update_force_PEB(n_atoms,n_images,start_coords,coords,end_coords,start_forces,total_forces,end_forces,start_energy,energies,end_energy,spring_constants,object_function)
        else
           write(use_unit,*) 'Adjusting forces with the climbing image nudged elastic band method.'
           call update_force_CINEB(n_atoms,n_images,start_coords,coords,end_coords,start_forces,total_forces,end_forces,start_energy,energies,end_energy,spring_constants)
        end if

     else if (method.eq.'NEB') then
        if (.not.use_CINEB_method) then
           write(use_unit,*) 'Adjusting forces with nudged elastic band method.'
           call update_force_NEB(n_atoms,n_images,start_coords,coords,end_coords,start_forces,total_forces,end_forces,start_energy,energies,end_energy,spring_constants)
        else
           write(use_unit,*) 'Adjusting forces with the climbing image nudged elastic band method.'
           call update_force_CINEB(n_atoms,n_images,start_coords,coords,end_coords,start_forces,total_forces,end_forces,start_energy,energies,end_energy,spring_constants)
        end if
     else if (method.eq.'CI-NEB') then
        write(use_unit,*) 'Adjusting forces with the climbing image nudged elastic band method.'
        call update_force_CINEB(n_atoms,n_images,start_coords,coords,end_coords,start_forces,total_forces,end_forces,start_energy,energies,end_energy,spring_constants)
     end if
     
     !---------------------------------------------------------------------------------------------------
     ! determine the maximum force on the NEB (not absolute on image!) - for convergence criterion
     max_coordinate_force = 0d0
     do i_images = 1, n_images
        do i_atom = 1, n_atoms
           do i_coords = 1, 3
              if (max_coordinate_force.lt.abs(total_forces((i_images-1)*3*n_atoms+(i_atom-1)*3+i_coords))) then
                 max_coordinate_force = abs(total_forces((i_images-1)*3*n_atoms+(i_atom-1)*3+i_coords))
                 imax_image = i_images
                 imax_atom  = i_atom
                 imax_coord = i_coords
              end if
           end do
        end do
     end do
     
     ! check if one should switch to CI-NEB:
     if (      (switch_to_CINEB_force.gt.0d0)                  &  ! does the option exist at all?
          .and.(max_coordinate_force.lt.switch_to_CINEB_force) &  ! is force convergence reached?
          .and.(.not.use_CINEB_method)                         &  ! are we already using CI-NEB?
          .and.(.not.line_search_on)                           &  ! are we NOT in the middle of a line-search (that should be completed first ... )?
          ) then
        use_CINEB_method = .true.
        write(use_unit,*) 'Initial force convergence reached. Switching to climbing image NEB.'
        write(use_unit,*) 'Recalculating the force adjustment with CI-NEB method.'
        total_forces(:) = forces_keep(:)
        call update_force_CINEB(n_atoms,n_images,start_coords,coords,end_coords,start_forces,total_forces,end_forces,start_energy,energies,end_energy,spring_constants)
        
        ! redetermine maximum coordinate force for proper convergence checking
        max_coordinate_force = 0d0
        do i_images = 1, n_images
           do i_atom = 1, n_atoms
              do i_coords = 1, 3
                 if (max_coordinate_force.lt.abs(total_forces((i_images-1)*3*n_atoms+(i_atom-1)*3+i_coords))) then
                    max_coordinate_force = abs(total_forces((i_images-1)*3*n_atoms+(i_atom-1)*3+i_coords))
                    imax_image = i_images
                    imax_atom  = i_atom
                    imax_coord = i_coords
                 end if
              end do
           end do
        end do
     end if
     

     !---------------------------------------------------------------------------------------------------
     ! calculate the angle between two successive iteration's forces     
     forcenew_magnitude = ddot(3*n_atoms, total_forces,1, total_forces,1)
     forceold_magnitude = ddot(3*n_atoms,stored_forces,1,stored_forces,1)
     force_angle        = ddot(3*n_atoms, total_forces,1,stored_forces,1)
     if (forcenew_magnitude*forceold_magnitude.gt.0d0) then
        force_angle        = 90d0*acos(abs(force_angle/sqrt(forcenew_magnitude*forceold_magnitude)))/asin(1d0)
     end if
     write(use_unit,'(2X,A,F14.8)') 'Angle between forces of last two iterations: ', force_angle        
     if (line_search_on.and.(force_angle.gt.max_line_search_angle)) then
        write(use_unit,'(2X,A)') '| This angle is too large to continue a line search.'
        line_search_on = .false.
     end if
     write(use_unit,*)
     write(use_unit,'(2X,A,F14.8)') 'Maximum force on original  image coordinate: ', max_force_before
     write(use_unit,'(2X,A,3I6)')   '               on (image, atom, coordinate): ', imax_image_before, imax_atom_before, imax_coord_before
     write(use_unit,'(2X,A,F14.8)') 'Maximum force on corrected image coordinate: ', max_coordinate_force
     write(use_unit,'(2X,A,3I6)')   '               on (image, atom, coordinate): ', imax_image, imax_atom, imax_coord
     write(use_unit,'(2X,A,F14.8)') '                Force convergence criterion: ', force_convergence_criterion
     if (object_function_exists) then
        write(use_unit,'(2X,A,F14.8)') '                   Value of object function: ', object_function
     end if
     write(use_unit,*)
     write(use_unit,'(2X,A,F14.8)') '                        Sum of DFT energies: ', sum_of_energies
     write(use_unit,'(2X,A,F14.8)') '                  Maximal energy difference: ', energy_barrier
     write(use_unit,*) 
     


!!$     !------DEBUG----------------------DEBUG-------------------------------------------------------------
!!$     ! calculate the net force of the incoming images, should be set to zero if significant ... 
!!$     write(use_unit,'(2X,A)') 'Cleaning translational components out of corrected forces.'
!!$     call get_net_force(n_atoms,n_images,total_forces,net_force)
!!$     write(use_unit,*)
!!$     write(use_unit,'(2X,A)')        'Net force on image AFTER updating forces, after cleaning, in eV/A '
!!$     write(use_unit,'(2X,A)')        '|     x       y        z'
!!$     do i_images = 1, n_images
!!$        write(use_unit,'(2X,A,3F10.6)') '|',(net_force(i_coords,i_images), i_coords = 1, 3)
!!$     end do
!!$     call clean_translational_forces(n_atoms,n_images,total_forces)
!!$     call get_net_force(n_atoms,n_images,total_forces,net_force)
!!$     write(use_unit,*)
!!$     write(use_unit,'(2X,A)')        'Net force on image AFTER updating forces, after cleaning, in eV/A '
!!$     write(use_unit,'(2X,A)')        '|     x       y        z'
!!$     do i_images = 1, n_images
!!$        write(use_unit,'(2X,A,3F10.6)') '|',(net_force(i_coords,i_images), i_coords = 1, 3)
!!$     end do
!!$
!!$     !------DEBUG----------------------DEBUG-------------------------------------------------------------



     if (max_coordinate_force.lt.force_convergence_criterion) then
        write(use_unit,*) 'Algorithm is converged. '
        write(use_unit,'(2X,3A)') ' File ', jmol_in_name, ' contains the final path, forces, and energies.'
     else
        !---------------------------------------------------------------------------------------------------
        ! call BFGS routine to get updated iteration 
        call BFGS_n_coords(n_coords,                  &
                           coords,                    &
                           total_forces,              &
                           ini_relaxation,            &
                           initial_hessian,           & 
                           stored_forces,             &
                           stored_coords,             &
                           search_direction,          &
                           hessian,                   &   
                           max_atomic_move,           &
                           min_line_step,             &
                           line_step_reduce,          &
                           line_step,                 &
                           object_function_exists,    &
                           object_function,           &
                           stored_object_function,    &
                           object_function_tolerance, &
                           line_search_on,            & 
                           force_times_search_dir,    &
                           max_line_step_increase)
        
        !---------------------------------------------------------------------------------------------------
        ! calculate the change in the transition state image before and after the geometry update, 
        ! also calculate and output the maximal force component on transition image - all else does not matter
        call get_difference_image(n_atoms,n_images,transition_image,coords,stored_coords,transition_max_change,transition_RMS_change)
        call get_max_force_image(n_atoms,n_images,transition_image,total_forces,transition_max_force)
        write(use_unit,*)
        write(use_unit,'(2X,A)') 'Information on transition state image:'
        write(use_unit,'(2X,A,F10.6)') ' | Maximal coordinate change in this iteration : ', transition_max_change
        write(use_unit,'(2X,A,F10.6)') ' |     RMS coordinate change in this iteration : ', transition_RMS_change
        write(use_unit,'(2X,A,F10.6)') ' |                     Maximal force component : ', transition_max_force
        write(use_unit,*)

        !---------------------------------------------------------------------------------------------------
        ! calculate and output difference vectors between the images: after BFGS
        write(use_unit,*)
        call get_difference_magnitude(n_atoms,n_images,start_coords,coords,end_coords,delta_R)
        write(use_unit,'(2X,A)') 'Magnitudes of image-image vectors after geometry update:'
        write(use_unit,'(2X,A,F10.6)') ' | start    -> image  1: ',delta_R(0)
        do i_images = 1, n_images - 1
           write(use_unit,'(2X,A,I2,A,I2,A,F10.6)') ' | image ',i_images,' -> image ',i_images+1,': ',delta_R(i_images)
        end do
        write(use_unit,'(2X,A,I2,A,F10.6)') ' | image ',n_images,' -> end     : ',delta_R(i_images)
     end if
     
     !---------------------------------------------------------------------------------------------------
     ! Check if the force criterion for switching to the BFGS has been reached
     ! reset Hessian to unity if not.
     if ((max_coordinate_force.gt.start_BFGS_force).and.(.not.line_search_on)) then 
        write(use_unit,*) 'Force convergence for starting BFGS not yet reached. Resetting Hessian. '
        hessian(:,:)           = 0d0
        stored_object_function = 0d0
        do i_coords = 1, n_coords
           hessian(i_coords,i_coords) = 1d0
        end do
        line_step              = 1d0        
     end if

     ! check if BFGS is being used and switch off the line search if that's the case
     if (max_coordinate_force.lt.start_BFGS_force) line_search_on = .false.

     !---------------------------------------------------------------------------------------------------
     ! output for next iteration

     if (use_CINEB_method) method = 'CI-NEB'

     write(use_unit,*)
     write(use_unit,'(2X,2A)') 'Writing geometries to file ', jmol_out_name
     ! write new iteration jmol-type movie
     open(unit = 7, file = jmol_out_name, status = 'unknown')
     ! write starting configuration
     write(7,'(I8)') n_atoms
     write(7,'(A,E18.8E4,2A)') 'START - Maximum force component in previous iteration = ', max_coordinate_force, ' method: ', method
     i_coords = 1
     do i_atom = 1, n_atoms 
        write(7,'(A,3E20.6)') species_names(i_atom), (start_coords(i_coords+i_buf), i_buf = 0, 2)
        i_coords = i_coords + 3
     end do
     ! write n_images images
     i_coords = 1
     do i_images = 1, n_images
        write(7,'(I8)') n_atoms
        write(7,'(A,I8)') 'image ',i_images
        do i_atom = 1, n_atoms 
           write(7,'(A,3E20.6)') species_names(i_atom), (coords(i_coords+i_buf), i_buf = 0, 2)
           i_coords = i_coords + 3
        end do
     end do
     ! writefinal configuration
     write(7,'(I8)') n_atoms
     write(7,'(A)') 'END '
     i_coords = 1
     do i_atom = 1, n_atoms 
        write(7,'(A,3E20.6)') species_names(i_atom), (end_coords(i_coords+i_buf), i_buf = 0, 2)
        i_coords = i_coords + 3
     end do
     close(unit = 7)
  
     ! write relaxation data for next iteration to file ... in binary format
     open(unit = 7, file = save_data_name, status = 'unknown', form = 'unformatted')
     do i_coords = 1, n_coords
        write(7) stored_coords(i_coords)
        write(7) stored_forces(i_coords)
        write(7) search_direction(i_coords)
        do i_coords2 = 1, n_coords
           write(7) hessian(i_coords,i_coords2)
        end do
     end do
     write(7) line_step
     write(7) stored_object_function
     write(7) use_CINEB_method
     write(7) line_search_on
     write(7) force_times_search_dir
     close(unit = 7)

     !---------------------------------------------------------------------------------------------------
     ! deallocate all variables
     deallocate(coords)
     deallocate(stored_coords)
     deallocate(total_forces)
     deallocate(forces_keep)
     deallocate(stored_forces)
     deallocate(search_direction)
     deallocate(hessian)
     deallocate(energies)
     deallocate(species_names)
     deallocate(start_forces)
     deallocate(start_coords)
     deallocate(end_forces)
     deallocate(end_coords)
     deallocate(delta_R)
     deallocate(spring_constants)
     deallocate(centre_of_mass)
     deallocate(net_force)
     deallocate(relative_coords)
     write(use_unit,*) '------------------------------------------------------------'
     write(use_unit,*) '     Exiting NEB searcher '
     write(use_unit,*) '------------------------------------------------------------'
  end if
  call finalize_mpi()
end program NEB
