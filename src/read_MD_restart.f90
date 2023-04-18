!------------------------------------------------------------------------------
!****s* FHI-aims/read_MD_restart
!  NAME
!    read_MD_restart
!  SYNOPSIS
subroutine read_MD_restart
!  PURPOSE
!     read restart information from file specified in runtime_choices. This is only done if 
!     requested. Note that this routine overwrites positions and velocities (after checking 
!     for the usefulness of the restart file) - this makes restart all the more simple as 
!     one would not even have to dig out the last known positions and velocities; The latter
!     are completely unimportant anyway as the integration algorithms give out velocities
!     at 1/2 time steps while the output is really in full time steps. 
!  USES
  use mpi_tasks
  use synchronize_mpi
  use synchronize_mpi_basic
  use runtime_choices
  use dimensions
  use physics
  use geometry
  use localorb_io
  use timing
  use constants 
  use species_data
  use relaxation
  use molecular_dynamics
  implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  logical :: MD_restart_file_exists
  integer :: n_atoms_internal, i_atom, i_periodic, i_coord, i, i_s
  real*8  :: s_read, s_read_last, MD_gle_ns_restart
  real*8, dimension(n_atoms)   :: masses
  real*8, dimension(n_atoms)   :: charges 
  real*8 :: plumed_eunit
  character*400 :: info_str,desc_str
  integer  :: i_code,coords_counter,v_half_counter,v_last_counter,r_last_counter
  integer :: gs_counter, gt_counter, gp_atom_counter, gp_s_counter!, ngp_atom_counter, ngp_s_counter
  logical  :: eof
  integer:: cell_type

  if(myid == 0)then
     inquire(FILE=MD_restart_file,EXIST=MD_restart_file_exists)
  end if
  call bcast_logical(MD_restart_file_exists, 0)
  if(.not. MD_restart_file_exists) then 
     write(info_str,'(2X,2A)') 'Could not find MD restart file: ', trim(MD_restart_file)
     call localorb_info(info_str)
     write(info_str,'(2X,A)') 'Returning to default initialization.'
     call localorb_info(info_str)
     return
  else
     write(info_str,'(2X,2A)') 'Attempting to restart MD using MD restart file: ', trim(MD_restart_file)
     call localorb_info(info_str)
  end if

  ! must initialize correctly as we only read on process 0 and then synchronize via allreduce
  MD_stepcount         = 0
  MD_force_evaluations = 0
  MD_high_order_i = 0
  ! if GLE thermostat, must allocate matrices
  if(MD_ensemble.eq.'GLE_thermostat')then
     call allocate_gle()
  endif
  if (myid.eq.0) then
     open(file = MD_restart_file, unit = 88, status = 'old', action='read')
     !read first line
     read (88,*,iostat=i_code) desc_str
     if (i_code.ne.0) then
        if (myid.eq.0) then
           write(use_unit,*) "* Attention - empty restart file ",trim(MD_restart_file),". "
        end if
        eof = .true.
        stop
     end if

     eof = .false.
     coords_counter = 0
     v_half_counter = 0
     v_last_counter = 0
     r_last_counter = 0
     gp_s_counter = 0
     gp_atom_counter=1
     gt_counter=0
     gs_counter=0
     ! Parsing: Please note that certain variables have different names in read/write subroutine
     ! write_MD		     |||    read_MD
     ! n_atoms               -->    n_atoms_internal
     ! coords                -->    coords
     ! v_half_write          -->    v_half
     ! v_last_write          -->    v_last
     ! r_last                -->    r_last
     ! s_write               -->    s_read
     ! s_NP_half             -->    s_NP_half
     ! s_last_write          -->    s_read_last
     ! s_dot_NP_half         -->    s_dot_NP_half
     ! s_dot_NP_last         -->    s_dot_NP_last
     ! tsystem               -->    tsystem
     ! tsystem_last          -->    tsystem_last
     ! MD_H0                 -->    MD_H0
     ! MD_Epot_last          -->    MD_Epot_last
     ! MD_stepcount          -->    MD_stepcount
     ! MD_high_order_i       -->    MD_high_order_i
     ! MD_force_evaluations  -->    MD_force_evaluations
     do while (.not.eof)
       if (desc_str(1:1).eq.'#') then
           continue
       ! n_atoms
       elseif (desc_str.eq.'n_atoms') then
         backspace(88)
         read (88,*) desc_str, n_atoms_internal
         if (n_atoms_internal.ne.n_atoms) then
           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
           write(use_unit,*) '* has the wrong number of atoms and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
       ! coords
       elseif (desc_str.eq.'coords') then
         coords_counter = coords_counter + 1
         if (coords_counter .gt. n_atoms) then
           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
           write(use_unit,*) '* includes too many coords tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(88)
         read (88,*) desc_str, ( coords(i_coord,coords_counter), i_coord=1,3,1)
       ! v_half_write
       elseif (desc_str.eq.'v_half') then
         v_half_counter = v_half_counter + 1
         if (v_half_counter .gt. n_atoms) then
           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
           write(use_unit,*) '* includes too many v_half tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(88)
         read (88,*) desc_str, ( v_half(i_coord,v_half_counter), i_coord=1,3,1)
       ! v_last_write
       elseif (desc_str.eq.'v_last') then
         v_last_counter = v_last_counter + 1
         if (v_last_counter .gt. n_atoms) then
           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
           write(use_unit,*) '* includes too many v_last tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(88)
         read (88,*) desc_str, ( v_last(i_coord,v_last_counter), i_coord=1,3,1)
       ! r_last
       elseif (desc_str.eq.'r_last') then
         r_last_counter = r_last_counter + 1
         if (r_last_counter .gt. n_atoms) then
           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
           write(use_unit,*) '* includes too many r_last tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(88)
         read (88,*) desc_str, ( r_last(i_coord,r_last_counter), i_coord=1,3,1)
       ! s_write
       elseif (desc_str.eq.'s_write') then
         backspace(88)
         read (88,*) desc_str, s_read 
       ! s_NP_half
       elseif (desc_str.eq.'s_NP_half') then
         backspace(88)
         read (88,*) desc_str, s_NP_half
       ! s_last_write
       elseif (desc_str.eq.'s_last_write') then
         backspace(88)
         read (88,*) desc_str, s_read_last
       ! s_dot_NP_half
       elseif (desc_str.eq.'s_dot_NP_half') then
         backspace(88)
         read (88,*) desc_str, s_dot_NP_half
       ! s_dot_NP_last
       elseif (desc_str.eq.'s_dot_NP_last') then
         backspace(88)
         read (88,*) desc_str, s_dot_NP_last
       ! tsystem
       elseif (desc_str.eq.'tsystem_half') then
         backspace(88)
         read (88,*) desc_str, tsystem
       ! tsystem_last
       elseif (desc_str.eq.'tsystem_last') then
         backspace(88)
         read (88,*) desc_str, tsystem_last
       ! MD_H0
       elseif (desc_str.eq.'MD_H0') then
         backspace(88)
         read (88,*) desc_str, MD_H0
       ! MD_Epot_last
       elseif (desc_str.eq.'MD_Epot_last') then
         backspace(88)
         read (88,*) desc_str, MD_Epot_last
       ! BDP_conint (BDP pseudo-Hamilt.)
       elseif (desc_str.eq.'BDP_psuedo-Hamilt.') then
         backspace(88)
         read (88,*) desc_str, BDP_conint
       ! MD_stepcount
       elseif (desc_str.eq.'MD_stepcount') then
         backspace(88)
         read (88,*) desc_str, MD_stepcount
       ! MD_high_order_i
       elseif (desc_str.eq.'MD_high_order_i') then
         backspace(88)
         read (88,*) desc_str, MD_high_order_i
       ! MD_force_evaluations
       elseif (desc_str.eq.'MD_force_evaluations') then
         backspace(88)
         read (88,*) desc_str, MD_force_evaluations
       elseif (desc_str.eq.'GLE_pseudo') then
         backspace(88)
         read (88,*) desc_str, langham
       elseif (desc_str.eq.'MD_gle_ns') then
         backspace(88)
         read (88,*) desc_str, MD_gle_ns_restart
         ! sanity check
         if (MD_gle_ns_restart.ne.MD_gle_ns) then
           write(use_unit,*) '* Warning: In the restart file ',trim(MD_restart_file)
           write(use_unit,*) '* the thermostat parameter is not the same as in control.in'
           write(use_unit,*) '* This cannot belong to the current system. Please correct.'
           write(use_unit,*) 'Aborting execution'
           stop
         endif
       elseif (desc_str.eq.'gp') then
         gp_s_counter=gp_s_counter+1
         if (gp_s_counter .gt. MD_gle_ns+1) then
            gp_atom_counter=gp_atom_counter+1
            gp_s_counter=1
         endif
         if (gp_s_counter*gp_atom_counter .gt. n_atoms*(MD_gle_ns+1)) then
           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
           write(use_unit,*) '* includes too many gp tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(88)
         read (88,*) desc_str, (gp(i_coord, gp_atom_counter, gp_s_counter) , i_coord=1,3,1)
!   ngp does not need to be reinitialized IMHO... but I leave the code here for now
!       elseif (desc_str.eq.'ngp') then
!         ngp_s_counter=ngp_s_counter+1
!         if (ngp_s_counter .gt. MD_gle_ns+1) then
!            ngp_atom_counter=ngp_atom_counter+1
!            ngp_s_counter=1
!         endif
!         if (ngp_s_counter*ngp_atom_counter .gt. n_atoms*(MD_gle_ns+1)) then
!           write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
!           write(use_unit,*) '* includes too many ngp tags and can impossibly belong'
!           write(use_unit,*) '* to the current system.'
!           write(use_unit,*) 'Aborting execution'
!           stop
!         end if
!         backspace(88)
!         read (88,*) desc_str, (ngp(i_coord, ngp_atom_counter, ngp_s_counter) , i_coord=1,3,1)
       elseif (desc_str.eq.'gT') then
            gt_counter=gt_counter+1
            backspace(88)
            read (88,*) desc_str, (gT(gt_counter, i_s) , i_s=1,MD_gle_ns+1,1)
       elseif (desc_str.eq.'gS') then
            gs_counter=gs_counter+1
            backspace(88)
            read (88,*) desc_str, (gS(gs_counter, i_s) , i_s=1,MD_gle_ns+1,1)
       else
          if (myid.eq.0) then
             write(use_unit,*) "Unknown descriptor ", desc_str, " in file ",trim(MD_restart_file),"."
          end if
          stop
       end if
       ! read next line
       read (88,*,iostat = i_code) desc_str
       if (i_code.ne.0) then
          eof = .true.
       end if
     end do
     close(unit=88)

     MD_successful_restart_read = .true.
     write(info_str,'(2X,2A)') 'Successfully read all MD restart data from file ', trim(MD_restart_file)
     call localorb_info(info_str)

     if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
        do i_atom = 1, n_atoms
          v_half(:,i_atom) = species_m(species(i_atom))*v_half(:,i_atom)
          v_last(:,i_atom) = species_m(species(i_atom))*v_last(:,i_atom)
        end do
     end if

     s_NP      = s_read
     s_NP_last = s_read_last
     if (MD_ensemble.eq.'NVT_nose-hoover') then
        s_NP      = dlog(s_read)
        s_NP_last = dlog(s_read_last)
     end if

     ! print initial geometry from restart file
         write(use_unit,'(2X,A)') "Input geometry:"

         if (n_periodic.gt.0) then
            write(use_unit,'(2X,A)') "| Unit cell: "
            do i_periodic = 1, n_periodic, 1
               write(use_unit,'(2X,A,3(2X,F15.6))') "|", &
                    (lattice_vector(i_coord,i_periodic)*bohr, &
                    i_coord=1,3,1)
            enddo
         else
            write(use_unit,'(2X,A)') "| No unit cell requested."
         end if

         write(use_unit,'(2X,A)') "| Atomic structure: "
         write(use_unit,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Atom ","x [A]","y [A]","z [A]"
         do i_atom = 1, n_occ_atoms, 1
            write(use_unit,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
                 "|",i_atom, ": Species ", &
                 species_name(species(i_atom)), &
                 (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
         enddo

         if (n_occ_atoms.lt.n_atoms) then
           write(use_unit,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Ghost","x [A]","y [A]","z [A]"
           do i_atom = n_occ_atoms+1, n_atoms, 1
             write(use_unit,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
                 "|",i_atom, ": Species ", &
                 species_name(species(i_atom)), &
                 (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
           enddo
         endif
  end if ! end myid=0

  if (MD_time_restart) then
     tsystem      = 0d0
     tsystem_last = 0d0
     MD_stepcount = 0
     MD_force_evaluations = 0
  end if
  call bcast_logical(MD_successful_restart_read,0)
  if (MD_successful_restart_read) then
     call broadcast_MD_velocities(coords,0)
     call broadcast_MD_velocities(v_half,0)
     call broadcast_MD_velocities(v_last,0)
     call broadcast_MD_velocities(r_last,0)
     call bcast_real(s_NP         ,0)
     call bcast_real(s_NP_last    ,0)
     call bcast_real(s_dot_NP_half,0)
     call bcast_real(s_dot_NP_last,0)
     call bcast_real(tsystem      ,0)
     call bcast_real(tsystem_last ,0)
     call bcast_real(MD_H0        ,0)
     call bcast_real(MD_Epot_last ,0)
     call bcast_real(BDP_conint   ,0)
     call sync_integer(MD_stepcount)
     call sync_integer(MD_high_order_i)
     call sync_integer(MD_force_evaluations)
     if(MD_ensemble.eq.'GLE_thermostat')then
        do i_s=1, MD_gle_ns+1, 1
          call broadcast_MD_velocities(gp(:,:,i_s),0)
        enddo
        call mp_bcast_mat(gT)
        call mp_bcast_mat(gS) 
     endif
  else 
     ! If we were not successful, we should never have reached this point in the code.
     ! Therefore, stop execution.
     write(info_str,'(2X,4A)') 'MD restart from file ',trim(MD_restart_file),' not successful. ',&
                               'At this point, partial (inconsistent) data may have been read.'
     call aims_stop(info_str)
  end if

  if (MD_use_schedule) then
       ! FIXME: Intelligent output about the initialization!!
       MD_schedule_step   = 1 
       MD_ensemble        = MD_schedule_ensemble      (MD_schedule_step)
       MD_time            = MD_schedule_time          (MD_schedule_step) 
       MD_temperature     = MD_schedule_temperature   (MD_schedule_step)
       MD_tau_berendsen   = MD_schedule_tau_berendsen (MD_schedule_step)
       MD_Q_NP            = MD_schedule_Q             (MD_schedule_step)
       NVE_damping_factor = MD_schedule_damping_factor(MD_schedule_step)
       MD_nu_andersen     = MD_schedule_nu_andersen   (MD_schedule_step)
       MD_tau_BDP         = MD_schedule_tau_BDP       (MD_schedule_step)
       write(info_str,'(2X,2A)') 'Restart MD schedule in ensemble: ',trim(MD_ensemble)
       call localorb_info(info_str)
  end if

  if((MD_ensemble.eq.'NVT_andersen').or.(MD_ensemble.eq.'NVT_parrinello').or.(MD_ensemble.eq.'GLE_thermostat')) then 
!the RNG will be initialize after restart, the alternative is to save the RNG status in the restart
     call initialize_RNG ()
     MD_RNG_firstcall = .true.
  endif

    if(plumed_plugin) then
	if(myid.eq.0) then
		do i_atom = 1, n_atoms
			masses(i_atom) = species_m(species(i_atom))
			charges(i_atom) = 0.d0
		end do
		!call init_metadyn_(n_atoms, MD_tstep, masses, trim(adjustl(plumed_file))//char(0))
                ! verison 1.3
                if ( n_periodic .ge.1) then
                    cell_type=1 
                else
                    cell_type=0 
                endif 
                plumed_eunit=1.d0
                call init_metadyn_(n_atoms, MD_tstep, masses, charges, &
                        cell_type, plumed_eunit, &
                        trim(adjustl(plumed_file))//char(0))  
                !write(use_unit,*) "fhi atoms MD_tstep plumed_file",n_atoms, MD_tstep,plumed_file
		!write(use_unit,*) "fhi masses", masses
	endif
    endif

    ! find out number of mobile atoms, if some are restricted 
    if (use_relaxation_constraints) then 
       n_atoms_MD = 0
       constrain_MD(:) = constrain_relaxation(:)
       do i_atom = 1, n_atoms
          if (.not.constrain_MD(i_atom)) then
             n_atoms_MD = n_atoms_MD + 1 
          end if
       end do
       ! reset number of degrees of freedom
       MD_g_DOF = 3*n_atoms_MD
    else
       n_atoms_MD = n_atoms
    end if
!DEBUG
!do i_atom=1, n_atoms
!  do i_s=1, MD_gle_ns+1
!    write(use_unit,*) 'gp in restart', i_s, gp(:, i_atom, i_s)
!    write(use_unit,*) 'ngp in restart', i_s, ngp(:, i_atom, i_s)
!  enddo
!enddo
! END DEBUG

end subroutine read_MD_restart
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/read_MD_restart_binary
!  NAME
!    read_MD_restart_binary
!  SYNOPSIS
subroutine read_MD_restart_binary ![OBSOLETE]
!  PURPOSE
!     read restart information from file specified in runtime_choices. This is only done if 
!     requested. Note that this routine overwrites positions and velocities (after checking 
!     for the usefulness of the restart file) - this makes restart all the more simple as 
!     one would not even have to dig out the last known positions and velocities; The latter
!     are completely unimportant anyway as the integration algorithms give out velocities
!     at 1/2 time steps while the output is really in full time steps. 
!    [OBSOLETE] A new read_MD_restart function is provided that uses a human
!    readable, cross-architecture-compatible format
!  USES
  use mpi_tasks
  use synchronize_mpi
  use runtime_choices
  use dimensions
  use physics
  use geometry
  use localorb_io
  use timing
  use constants 
  use species_data
  use relaxation
  use molecular_dynamics
  implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
  logical :: MD_restart_file_exists
  integer :: n_atoms_internal, i_atom, i_periodic, i_coord, i
  real*8  :: s_read, s_read_last
  real*8, dimension(n_atoms)   :: masses
  real*8, dimension(n_atoms)   :: charges 
  real*8 :: plumed_eunit   
  character*100 :: info_str
  integer :: cell_type

  if(myid == 0)then
     inquire(FILE=MD_restart_file,EXIST=MD_restart_file_exists)
  end if
  call bcast_logical(MD_restart_file_exists, 0)
  if(.not. MD_restart_file_exists) then 
     write(info_str,'(2X,2A)') 'Could not find MD restart file: ', trim(MD_restart_file)
     call localorb_info(info_str)
     write(info_str,'(2X,A)') 'Returning to default initialization.'
     call localorb_info(info_str)
     return
  end if
  MD_stepcount         = 0
  MD_force_evaluations = 0
  if (myid.eq.0) then
     open(file = MD_restart_file, unit = 88, status = 'old', form = 'unformatted', action='read')
     read(88) n_atoms_internal
     if (n_atoms_internal.ne.n_atoms) then
        write(use_unit,*) '* Warning: The restart file ',trim(MD_restart_file)
        write(use_unit,*) '* has the wrong number of atoms and can impossibly belong'
        write(use_unit,*) '* to the current system.'
        write(use_unit,*) 'Aborting execution'
        stop
     end if
     read(88) coords
     read(88) v_half
     read(88) v_last
     read(88) r_last
     read(88) s_read
     read(88) s_NP_half
     read(88) s_read_last
     read(88) s_dot_NP_half
     read(88) s_dot_NP_last
     read(88) tsystem
     read(88) tsystem_last
     read(88) MD_H0
     read(88) MD_Epot_last
     read(88) BDP_conint
     read(88) MD_stepcount
     read(88) MD_high_order_i
     read(88) MD_force_evaluations
     close(unit=88)

     if ((MD_ensemble.eq.'NVT_nose-poincare').or.(MD_ensemble.eq.'NVT_nose-hoover')) then
     	do i_atom = 1, n_atoms
        	v_half(:,i_atom) = species_m(species(i_atom))*v_half(:,i_atom)
		v_last(:,i_atom) = species_m(species(i_atom))*v_last(:,i_atom)
     	end do
     end if	

     s_NP      = s_read
     s_NP_last = s_read_last
     if (MD_ensemble.eq.'NVT_nose-hoover') then
        s_NP      = dlog(s_read)
        s_NP_last = dlog(s_read_last)
     end if
     MD_successful_restart_read = .true.
     write(info_str,'(2X,2A)') 'Successfully read all MD restart data from file ', trim(MD_restart_file)
     call localorb_info(info_str)
! print initial geometry from restart file
         write(use_unit,'(2X,A)') "Input geometry:"

         if (n_periodic.gt.0) then
            write(use_unit,'(2X,A)') "| Unit cell: "
            do i_periodic = 1, n_periodic, 1
               write(use_unit,'(2X,A,3(2X,F15.6))') "|", &
                    (lattice_vector(i_coord,i_periodic)*bohr, &
                    i_coord=1,3,1)
            enddo
         else
            write(use_unit,'(2X,A)') "| No unit cell requested."
         end if

         write(use_unit,'(2X,A)') "| Atomic structure: "
         write(use_unit,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Atom ","x [A]","y [A]","z [A]"
         do i_atom = 1, n_occ_atoms, 1
            write(use_unit,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
                 "|",i_atom, ": Species ", &
                 species_name(species(i_atom)), &
                 (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
         enddo

         if (n_occ_atoms.lt.n_atoms) then
           write(use_unit,'(2X,A1,7X,A,15X,A,12X,A,12X,A)') &
              "|","Ghost","x [A]","y [A]","z [A]"
           do i_atom = n_occ_atoms+1, n_atoms, 1
             write(use_unit,'(2X,A1,I5,A,A2,3(2X,F15.6))') &
                 "|",i_atom, ": Species ", &
                 species_name(species(i_atom)), &
                 (coords(i_coord,i_atom)*bohr, i_coord=1,3,1)
           enddo
         endif
  else 
     MD_successful_restart_read = .false.     
     write(info_str,'(2X,3A)') 'MD restart from file ',trim(MD_restart_file),' not succesful.'
     call localorb_info(info_str)
     write(info_str,'(2X,A)') 'Using regular initialization for molecular dynamics run.'
     call localorb_info(info_str)
  end if
  if (MD_time_restart) then
     tsystem      = 0d0
     tsystem_last = 0d0
     MD_stepcount = 0
     MD_force_evaluations = 0
  end if
  call bcast_logical(MD_successful_restart_read,0)
  if (MD_successful_restart_read) then
     call broadcast_MD_velocities(coords,0)
     call broadcast_MD_velocities(v_half,0)
     call broadcast_MD_velocities(v_last,0)
     call broadcast_MD_velocities(r_last,0)
     call bcast_real(s_NP         ,0)
     call bcast_real(s_NP_last    ,0)
     call bcast_real(s_dot_NP_half,0)
     call bcast_real(s_dot_NP_last,0)
     call bcast_real(tsystem      ,0)
     call bcast_real(tsystem_last ,0)
     call bcast_real(MD_H0        ,0)
     call bcast_real(MD_Epot_last ,0)
     call bcast_real(BDP_conint   ,0)
     call sync_integer(MD_stepcount)
     call sync_integer(MD_high_order_i)
     call sync_integer(MD_force_evaluations)
  end if

  if (MD_use_schedule) then
       ! FIXME: Intelligent output about the initialization!!
       MD_schedule_step   = 1 
       MD_ensemble        = MD_schedule_ensemble      (MD_schedule_step)
       MD_time            = MD_schedule_time          (MD_schedule_step) 
       MD_temperature     = MD_schedule_temperature   (MD_schedule_step)
       MD_tau_berendsen   = MD_schedule_tau_berendsen (MD_schedule_step)
       MD_Q_NP            = MD_schedule_Q             (MD_schedule_step)
       NVE_damping_factor = MD_schedule_damping_factor(MD_schedule_step)
       MD_nu_andersen     = MD_schedule_nu_andersen   (MD_schedule_step)
       MD_tau_BDP         = MD_schedule_tau_BDP       (MD_schedule_step)
       write(info_str,'(2X,2A)') 'Restart MD schedule in ensemble: ',trim(MD_ensemble)
       call localorb_info(info_str)
  end if

  if((MD_ensemble.eq.'NVT_andersen').or.(MD_ensemble.eq.'NVT_parrinello')) then 
!the RNG will be initialize after restart, the alternative is to save the RNG status in the restart
     call initialize_RNG ()
     MD_RNG_firstcall = .true.
  endif

    if(plumed_plugin) then
	if(myid.eq.0) then
		do i_atom = 1, n_atoms
			masses(i_atom) = species_m(species(i_atom))
			charges(i_atom) = 0.d0 
		end do
		!call init_metadyn_(n_atoms, MD_tstep, masses, trim(adjustl(plumed_file))//char(0))
                if ( n_periodic .ge.1) then
                    cell_type=1 
                else
                    cell_type=0 
                endif 
                plumed_eunit=1.d0
                call init_metadyn_(n_atoms, MD_tstep, masses, charges, &
                        cell_type, plumed_eunit, &
                        trim(adjustl(plumed_file))//char(0))  
                !write(use_unit,*) "fhi atoms MD_tstep plumed_file",n_atoms, MD_tstep,plumed_file
		!write(use_unit,*) "fhi masses", masses
	endif
    endif

    ! find out number of mobile atoms, if some are restricted 
    if (use_relaxation_constraints) then 
       n_atoms_MD = 0
       constrain_MD(:) = constrain_relaxation(:)
       do i_atom = 1, n_atoms
          if (.not.constrain_MD(i_atom)) then
             n_atoms_MD = n_atoms_MD + 1 
          end if
       end do
       ! reset number of degrees of freedom
       MD_g_DOF = 3*n_atoms_MD
    else
       n_atoms_MD = n_atoms
    end if

end subroutine read_MD_restart_binary
!******
