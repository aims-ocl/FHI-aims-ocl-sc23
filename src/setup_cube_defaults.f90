!****s* FHI-aims/setup_cube_defaults
!  NAME
!  setup_cube_defaults.f90
!  SYNOPSIS
     subroutine setup_cube_defaults()
!!  PURPOSE
! If cube files are requested, but
! not further specified, sensible default values
! will be employed. The goal is to thus 
! make generation of these files more
! convenient and also more fail-save. 
!

  !  INPUTS
  !   none
  !  OUTPUT
  !    cube_origin(n_cube), cube_edge_steps(n_cube), cube_edge_unit(n_cube) 
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


 
!  USES
  use geometry
  use constants
  use plot
  use runtime_choices
  use mpi_tasks, only: aims_stop, myid
  use physics
  use dimensions, only: flag_out_elec_real 


 implicit none

!local variables
 real*8, dimension(3) :: geom_min  !minimal values of x,y,z geometry
 real*8, dimension(3) :: geom_max  !minimal values of x,y,z geometry
 real*8, dimension(3,n_atoms) :: coord_temp !temporary coordinates
 real*8 :: vec_len !vector length
 real*8 :: slab_len !length of the slab
! real*8 :: step_length
 integer :: cube_edge_number
 character(LEN=5) :: file_suffix !suffix of the file name
 character(LEN=120) :: info_str
 logical :: cube_exists = .false.
 integer :: cube_increment
 character(LEN=120) :: original_name

 !counter
 integer :: i_cube !(i-th cube)
 integer :: i_coord !(1-3: x,y,z)
 integer :: i_counter
 integer :: i_edge

 logical :: slab_calculation 



!!Begin work

!Default for cube divisor
cube_divisor(0) = 10

!Initialize variables
    geom_min(1:3) = 1d10 !Initializing with 0 causes wrong results if only negative coordinates are supplied
    geom_max(1:3) = -1d10 !Like above for positive coordinates
    slab_calculation = use_dipole_correction !Is there a better telltale variable?

!Copy the coordinates so we can modify them
    coord_temp(1:3,1:n_atoms) = coords(1:3,1:n_atoms)
!For slab calculation, shift everything below the vacuum level
    if (slab_calculation) then
        do i_coord = 1,n_atoms,1
           if (coord_temp(3,i_coord).gt.vacuum_z_level) then
               coord_temp(1,i_coord) = coord_temp(1,i_coord)-lattice_vector(1,3)
               coord_temp(2,i_coord) = coord_temp(2,i_coord)-lattice_vector(2,3)
               coord_temp(3,i_coord) = coord_temp(3,i_coord)-lattice_vector(3,3)
           endif
        enddo
    endif

!Tell the user it's cube time now
      if (myid.eq.0) then
        call localorb_info("")
        call localorb_info("  ------------------------------------------------------------")
        call localorb_info("                      Starting Cube Output                    ")
        call localorb_info("  ------------------------------------------------------------")
        call localorb_info("")
        call localorb_info("  The following cube files will be created:")
        call localorb_info("")
        call localorb_info( "  --------------------------------------")
      endif

  coord_temp = coords
! For slab calculation, shift the coords around somewhat
  

!Start employing default values if necessary
!Do some analysis of the geometry in case defaults are required
    do i_coord = 1,n_atoms,1 !Determine min and max values of x,y,z
         geom_min(1) = min(geom_min(1),coord_temp(1,i_coord))
         geom_max(1) = max(geom_max(1),coord_temp(1,i_coord))
         geom_min(2) = min(geom_min(2),coord_temp(2,i_coord))
         geom_max(2) = max(geom_max(2),coord_temp(2,i_coord))
         geom_min(3) = min(geom_min(3),coord_temp(3,i_coord))
         geom_max(3) = max(geom_max(3),coord_temp(3,i_coord))
    enddo

!Calculate default for cube origin
          if(n_periodic.eq.0) then !set cube origin to center of system
             do i_coord = 1,3,1
                 cube_origin(i_coord,0)= &
                 0.5*(geom_max(i_coord)+geom_min(i_coord)) !Take average of min and max
             enddo
          elseif(n_periodic.eq.3) then
             !if (flag_out_elec_real) then 
             !   cube_origin(:,0) = 0.0d0 
             !else
               do i_coord = 1,3,1
                 !cube_origin(i_coord,0) = 0 !Set to origin.  
                 !OTH: As of Oct 2017, changing the default to center of unit cell
                 cube_origin(1,0) = (lattice_vector(1,1)+lattice_vector(1,2)+lattice_vector(1,3)) * 0.5
                 cube_origin(2,0) = (lattice_vector(2,1)+lattice_vector(2,2)+lattice_vector(2,3)) * 0.5
                 cube_origin(3,0) = (lattice_vector(3,1)+lattice_vector(3,2)+lattice_vector(3,3)) * 0.5
                 if (slab_calculation) cube_origin(3,0) =  0.5*(geom_max(3)+geom_min(3))
               enddo
            !endif
          else !Included for future support. Not sure what would be a sensible default, as that depends on how AIMS will handle n_periodic 1 and 2
                    info_str= "No default origin for cube available &
                           &   for periodic systems with less than 3 dimensions. &
                           &   Please set manually "
                    if (myid.eq.0) call localorb_info(info_str)
         endif
!Set default for for cube edges and steps
        if(n_periodic.eq.0) then
           do i_coord = 1,3,1
                cube_edge_unit(1,i_coord,0)=0.0  !Make sure everything is zero
                cube_edge_unit(2,i_coord,0)=0.0
                cube_edge_unit(3,i_coord,0)=0.0
!                cube_edge_unit(i_coord,i_coord,0)=0.10
                cube_edge_unit(i_coord,i_coord,0)=0.10/bohr !and then set diagonal terms to 0.1A 
                cube_edge_steps(i_coord,0)=(geom_max(i_coord)-geom_min(i_coord)+14)/cube_edge_unit(i_coord,i_coord,0) 
            enddo
         elseif(n_periodic.eq.3) then
            !OTH: As of Oct 2017, change of default to 100 steps along each
            !lattice vector. This avoids rounding errors and always leads to a
            !sensible size of the cube files
            if (flag_out_elec_real) then 
               cube_edge_steps(1,0) = step_x 
               cube_edge_steps(2,0) = step_y
               cube_edge_steps(3,0) = step_z
               do i_edge = 1,3,1
                  do i_coord = 1,3,1 
                     cube_edge_unit(i_coord,i_edge,0) = lattice_vector(i_coord,i_edge)/ cube_edge_steps(i_edge,0)
                  enddo
               enddo
               if (slab_calculation .and. (vacuum_length .ne. 0)) then    !Exception for slab calculations - only plot region near the slab 
                slab_len=geom_max(3)-geom_min(3)+vacuum_length*1.8897259886    
                cube_edge_unit(3,3,0) = slab_len/cube_edge_steps(3,0)
                !cube_edge_steps(3,0)=floor(slab_len/(0.1/bohr))
                !step_length=vec_len/cube_edge_steps(3,0)
                !cube_edge_unit(3,3,0)=0.1/bohr ! step_length
               endif
         
            else
               do i_edge =1,3,1
               cube_edge_steps(i_edge,0) = 100 
                  do i_coord = 1,3,1
                     cube_edge_unit(i_coord,i_edge,0) = lattice_vector(i_coord,i_edge) / cube_edge_steps(i_edge,0)
!                  cube_edge_unit(i_coord,i_edge,0) =  lattice_vector(i_coord,i_edge) / cube_edge_steps(i_edge,0)/bohr !I shouldn't have to use bohr, not sure why
                  enddo
               enddo
            
            !do i_edge = 1,3,1
            !     vec_len=sqrt(lattice_vector(1,i_edge)*lattice_vector(1,i_edge) &
            !       + lattice_vector(2,i_edge)*lattice_vector(2,i_edge) &
            !       +lattice_vector(3,i_edge)*lattice_vector(3,i_edge))
!
!                   cube_edge_steps(i_edge,0)=floor(vec_len/(0.1/bohr))
!                   step_length=vec_len/cube_edge_steps(i_edge,0)
!               do i_coord = 1,3,1 !!Norm to 0.1A length, keep direction of lattice-vector
!                       cube_edge_unit(i_coord,i_edge,0)=lattice_vector(i_coord,i_edge)*step_length/vec_len
!                enddo
!           enddo

           if (slab_calculation) then    !Exception for slab calculations - only plot region near the slab 
                slab_len=geom_max(3)-geom_min(3)+14
                cube_edge_steps(3,0)=floor(slab_len/(0.1/bohr))
                !step_length=vec_len/cube_edge_steps(3,0)
                cube_edge_unit(3,3,0)=0.1/bohr ! step_length
          endif
         endif

        else !Included for future support. Not sure what would be a sensible default, as that depends on how AIMS will handle n_periodic 1 and 2
          info_str =  "No default cube edges for cube available &
                  &   for periodic systems with less than 3 dimensions. &
                  &  Please set manually "
          if (myid.eq.0) call localorb_info(info_str) 
        endif

!Now provide all cubes with grids and origin. If no new grid is specified
!use the previous one
    do i_cube = 1,n_cube,1

       ! Check if HOMO or LUMO reference has been requested
       ! This way of indexing is bad; it will fail for metallic systems. It also
       ! requires a double -> integer conversion, which is unsettling.
       if (cube_fermi_ref(i_cube)) then
         if (out_cube_soc) then
           cube_index(i_cube) = cube_index(i_cube) + int(n_electrons)
         else 
           cube_index(i_cube) = cube_index(i_cube) + int(ceiling(n_electrons/2.0d0))
         end if
       end if

     !State the cube file
     write(info_str,*) '   Cube file: ', i_cube
     call localorb_info(info_str)

!If no origin has been set, use value of the previous cube
!Since the loop starts with "1", it will start by refering to the
!the default values defined above if nothing else is provided
      if (.not.flag_origin(i_cube)) then
                 do i_coord = 1,3,1
                   cube_origin(i_coord,i_cube)= &
                   cube_origin(i_coord,i_cube-1)
                 enddo
      endif


!Output the effectively used origin of the cube file
      if (myid.eq.0) then
              write(info_str,'(2X,A,3(1X,F11.4))') & 
!              & "| Origin : ", &
               & "  cube origin ", &
              & (cube_origin(i_coord,i_cube)*bohr, i_coord=1,3,1)
              call localorb_info(info_str)
      endif

!If no edge has been set, use the edge of the previous cube file
     if(.not.flag_edge(i_cube)) then
           do i_coord = 1,3,1
             do i_edge = 1,3,1
                cube_edge_unit(i_edge,i_coord,i_cube)=cube_edge_unit(i_edge,i_coord,i_cube-1) 
                cube_edge_steps(i_edge,i_cube)=cube_edge_steps(i_edge,i_cube-1)
             enddo
           enddo
      endif

!Output the effectively used Grid direction and number of steps
        if (myid.eq.0) then
          do i_counter=1,3,1
!              write(info_str,'(2X,A,I3,A,I6,A,3(1X,F11.4))') & 
               write(info_str,'(2X,A,I0,3(1X,F11.4))') &
!         &     "| Grid step direction", i_counter, ": ", &
          &     "  cube edge ", &
               cube_edge_steps(i_counter, i_cube), &
 !             cube_edge_steps(i_counter,i_cube), " steps of ", &
                
              (cube_edge_unit(i_coord,i_counter,i_cube)*bohr, i_coord=1,3,1)
             call localorb_info(info_str)
          enddo
       endif

!Check if cube divisor has been requested
        if (flag_divisor(i_cube).eqv..false.) then
           if(cube_type_needs_densmat(i_cube).and.packed_matrix_format==PM_index) then !safeguard
             cube_divisor(i_cube) = 1
           else
             cube_divisor(i_cube) = cube_divisor(i_cube-1)
           endif
        endif
        if (flag_divisor(i_cube).eqv..true.) then !Show that only if explicitely asked for
             write(info_str,'(2X,A,I0)') ' cube divisor ', cube_divisor(i_cube)
             if (myid.eq.0) call localorb_info(info_str)
        endif

!If no file name is set, give a default file name. 
!Check for format
      select case(cube_format(i_cube))
             case('xsf') 
                file_suffix='.xsf'
             case default 
              file_suffix='.cube'  !Assume that cube_format can only hold sensible choices - nonsense input is checked earlier
      end select

      if(myid.eq.0) then
        !cube_filename(i_cube)=cube_inputname(i_cube)
        if(cube_filename(i_cube)=="") then 
              if (cube_type(i_cube).eq.'total_density') then

                 if(cube_spin(1,i_cube).eq.1.and.cube_spin(2,i_cube).eq.1) then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_total_density", file_suffix
                 elseif (cube_spin(1,i_cube).eq.1.and.cube_spin(2,i_cube).eq.0) then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_density_spin_up", file_suffix
                 elseif (cube_spin(1,i_cube).eq.0.and.cube_spin(2,i_cube).eq.1) then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_density_spin_down", file_suffix
                 else
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_density_user_defined_spin_mask", file_suffix
                 endif
             else if (cube_type(i_cube).eq.'spin_density') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_spin_density", file_suffix
             else if (cube_type(i_cube).eq.'delta_density') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_delta_density", file_suffix
             else if (cube_type(i_cube).eq.'dielec_func') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_dielec_func", file_suffix
             else if (cube_type(i_cube).eq.'ion_dens') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_ion_dens", file_suffix             
             else if (cube_type(i_cube).eq.'delta_v') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                           "cube_", i_cube, "_delta_v", file_suffix                               
             else if (cube_type(i_cube).eq.'eigenstate'.or. &
                      cube_type(i_cube).eq.'eigenstate_non_soc' .or. &
                      cube_type(i_cube).eq.'eigenstate_imag' .or. &
                      cube_type(i_cube).eq.'eigenstate_density' .or. &
                      cube_type(i_cube).eq.'xc_potential') then
	     
                    if (cube_type(i_cube).eq.'eigenstate') then
                         if(n_periodic .gt. 0) then 
                            write (unit=cube_filename(i_cube),fmt='(A,I3.3,A,I5.5,A,I1.1,A,I4.4,A)') &
                           "cube_", i_cube, '_eigenstate_',cube_index(i_cube),"_spin_", &
                           cube_state(1,i_cube), "_k_point_", cube_state(2,i_cube), file_suffix
                         else
                            cube_state(2,i_cube)=1
                            write (unit=cube_filename(i_cube),fmt='(A,I3.3, A,I5.5,A,I1.1,A)') &
                            "cube_", i_cube, '_eigenstate_',cube_index(i_cube),"_spin_",&
                            cube_state(1,i_cube), file_suffix
                        endif
                    elseif (cube_type(i_cube).eq.'eigenstate_non_soc') then
                         if(n_periodic .gt. 0) then 
                            write (unit=cube_filename(i_cube),fmt='(A,I3.3,A,I5.5,A,I1.1,A,I4.4,A)') &
                           "cube_", i_cube, '_eigenstate_non_soc_',cube_index(i_cube),"_spin_", &
                           cube_state(1,i_cube), "_k_point_", cube_state(2,i_cube), file_suffix
                         else
                            cube_state(2,i_cube)=1
                            write (unit=cube_filename(i_cube),fmt='(A,I3.3,A,I5.5,A,I1.1,A)') &
                            "cube_", i_cube, '_eigenstate_non_soc_',cube_index(i_cube),"_spin_",&
                            cube_state(1,i_cube), file_suffix
                        endif
                    elseif(cube_type(i_cube).eq.'eigenstate_imag') then 
                        if(n_periodic .gt. 0) then !just for consistency, doesn't make sense, because eigenstates are real in the cluster case
                            write (unit=cube_filename(i_cube),fmt='(A,I3.3,A,I5.5,A,I1.1,A,I4.4,A)') &
                           "cube_", i_cube, '_eigenstate_imag_',cube_index(i_cube),"_spin_", &
                           cube_state(1,i_cube), "_k_point_", cube_state(2,i_cube), file_suffix
                         else
                            cube_state(2,i_cube)=1
                            write (unit=cube_filename(i_cube),fmt='(A,I3.3, A,I5.5,A,I1.1,A)') &
                            "cube_", i_cube, '_eigenstate_imag_',cube_index(i_cube),"_spin_",&
                            cube_state(1,i_cube), file_suffix
                        endif
                    elseif(cube_type(i_cube).eq.'xc_potential') then 
                           cube_state(2,i_cube)=1
                           write (unit=cube_filename(i_cube),fmt='(A,I3.3,A,I5.5,A,I1.1,A)') &
                           "cube_", i_cube, '_xc_potential_spin_', &
                           cube_state(1,i_cube), file_suffix
                    else
                           if(cube_stm(3,i_cube).lt.0) then
                               if(n_periodic .gt. 0) then
                                   write (unit=cube_filename(i_cube),fmt='(A, I3.3, A,I5.5,A,I1.1,A,I4.4,A)') &
                                   "cube_", i_cube, '_eigenstate_density_', &
                                  cube_index(i_cube),"_spin_", &
                                 cube_state(1,i_cube), "_k_point_", cube_state(2,i_cube),file_suffix
                               else
                                  cube_state(2,i_cube)=1
                                  write (unit=cube_filename(i_cube),fmt='(A, I3.3, A,I5.5,A,I1.1,A)') &
                                  "cube_", i_cube, '_eigenstate_density_', &
                                  cube_index(i_cube),"_spin_", cube_state(1,i_cube), file_suffix
                               endif
                           endif
                    endif
            elseif(cube_type(i_cube).eq.'stm') then 
                    write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,I2.2,A)') &
                    "cube_", i_cube , &
                   '_stm_',cube_index(i_cube),file_suffix
            elseif(cube_type(i_cube).eq.'potential') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                    "cube_", i_cube , &
                   '_potential',file_suffix
            elseif(cube_type(i_cube).eq.'hartree_potential') then
                  if (flag_out_elec_real) then 
                     write(unit=cube_filename(i_cube),fmt='(A)') "realspace_ESP.out"
                  else
                     write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                    "cube_", i_cube , &
                   '_hartree_potential',file_suffix
                  endif
            elseif(cube_type(i_cube).eq.'xc_potential') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                    "cube_", i_cube , &
                   '_xc_potential',file_suffix
            elseif(cube_type(i_cube).eq.'long_range_potential') then
                  write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                    "cube_", i_cube , &
                   '_long_range_potential',file_suffix

            elseif(cube_type(i_cube).eq.'elf') then
                   write(unit=cube_filename(i_cube),fmt='(A,I3.3,A,A)') &
                    "cube_", i_cube , &
                   '_elf',file_suffix
            else
                 call aims_stop('Cube type not recognized. Stopping')
            endif
        endif
      if (overwrite_existing_cube_files.eqv..false.) then
         INQUIRE(FILE=cube_filename(i_cube), EXIST=cube_exists)
         original_name=cube_filename(i_cube)
         cube_increment=1
         do while (cube_exists.eqv..true.) 
            write(cube_filename(i_cube),fmt='(A,I0)') trim(original_name), cube_increment
            INQUIRE(FILE=cube_filename(i_cube), EXIST=cube_exists)
            cube_increment=cube_increment+1
            if (cube_increment==100) then
               write(info_str,*) '100 attempts to generate new cube name unsucessful'
               call localorb_info(info_str)
               stop
            endif
         enddo
      endif

          write(info_str,*) '   cube filename: ', cube_filename(i_cube)
          call localorb_info(info_str)
          call localorb_info( "  --------------------------------------")
      endif

!End of setting defaults
!
!!Check for estimated size of cube-file
!!This is only a crude check - stop only if defaults are used
!!TODO provide keyword for maximum size of cube file
         cube_edge_number=cube_edge_steps(1,i_cube)*cube_edge_steps(2,i_cube)*cube_edge_steps(3,i_cube)
         if ((.not.any(flag_edge(1:i_cube))).and. &
             cube_edge_number.gt.cube_default_size_safeguard) then
             write(use_unit,*) 'Size of cube file estimated to be too large (~', cube_edge_number/8, ' Bytes)'
             write(use_unit,*) 'To make sure you are not unintentionally clogging up the hard dis, we stop here.'
             write(use_unit,*) 'In order to resolve this issue, please either '
             write(use_unit,*) '  - set up the cube edges manually (see manual), or '
             write(use_unit,*) '  - use the keyword cube_default_size_safeguard (size)'
             write(use_unit,*) '  - to increase the allowed size for cube files using the default setting '
             write(use_unit,*) 'Apologies for the inconvenience.'
             call aims_stop('Safeguard stop when setting up cube files. Please read the error message above.')
          endif


    enddo
          end subroutine setup_cube_defaults
! end work
