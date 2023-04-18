!------------------------------------------------------------------------------
!****s* FHI-aims/read_PIMD_restart
!  NAME
!    read_PIMD_restart
!  SYNOPSIS
subroutine read_PIMD_restart
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
  use pi_molecular_dynamics
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
  logical :: PIMD_restart_file_exists
  integer :: n_atoms_internal, i_atom, i_bead, i_periodic, i_coord, i, n_beads_internal
  real*8  :: s_read, s_read_last
  real*8, dimension(n_atoms)   :: masses
  character*120 :: info_str,desc_str
  integer  :: i_code,coords_counter,v_half_counter,v_last_counter,r_last_counter,bead_counter
  logical  :: eof

  if(myid == 0)then
     inquire(FILE=PIMD_restart_file,EXIST=PIMD_restart_file_exists)
  end if
  call bcast_logical(PIMD_restart_file_exists, 0)
  if(.not. PIMD_restart_file_exists) then 
     write(info_str,'(2X,2A)') 'Could not find PIMD restart file: ', trim(PIMD_restart_file)
     call localorb_info(info_str)
     write(info_str,'(2X,A)') 'Returning to default initialization.'
     call localorb_info(info_str)
     return
  else
     write(info_str,'(2X,2A)') 'Attempting to restart MD using PIMD restart file: ', trim(PIMD_restart_file)
     call localorb_info(info_str)
  end if

  ! must initialize correctly as we only read on process 0 and then synchronize via allreduce
  PIMD_stepcount         = 0
  PIMD_force_evaluations = 0

  if (myid.eq.0) then
     open(file = PIMD_restart_file, unit = 89, status = 'old', action='read')
     !read first line
     read (89,*,iostat=i_code) desc_str
     if (i_code.ne.0) then
        if (myid.eq.0) then
           write(use_unit,*) "* Attention - empty restart file ",trim(PIMD_restart_file),". "
        end if
        eof = .true.
        stop
     end if

     eof = .false.
     coords_counter = 0
     v_half_counter = 0
     v_last_counter = 0
     r_last_counter = 0
     bead_counter   = 0
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
         backspace(89)
         read (89,*) desc_str, n_atoms_internal
         if (n_atoms_internal.ne.n_atoms) then
           write(use_unit,*) '* Warning: The restart file ',trim(PIMD_restart_file)
           write(use_unit,*) '* has the wrong number of atoms and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
       ! n_beads
       elseif (desc_str.eq.'n_beads') then
         backspace(89)
         read (89,*) desc_str, n_beads_internal
         if (n_beads_internal.ne.n_beads) then
           write(use_unit,*) '* Warning: The restart file ',trim(PIMD_restart_file)
           write(use_unit,*) '* has the wrong number of beads and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
       ! coords_beads
       elseif (desc_str.eq.'coords_beads') then
         coords_counter = coords_counter + 1
         if (coords_counter .gt. n_atoms*n_beads) then
           write(use_unit,*) '* Warning: The restart file ',trim(PIMD_restart_file)
           write(use_unit,*) '* includes too many coords tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(89)
         i_bead = mod(coords_counter,n_beads)
         if (i_bead.eq.0) i_bead = n_beads
         i_atom = int((coords_counter - i_bead)/n_beads) + 1
         read (89,*) desc_str, ( coords_beads(i_coord,i_bead,i_atom), i_coord=1,3,1)
       ! v_beads_half
       elseif (desc_str.eq.'v_beads_half') then
         v_half_counter = v_half_counter + 1
         if (v_half_counter .gt. n_atoms*n_beads) then
           write(use_unit,*) '* Warning: The restart file ',trim(PIMD_restart_file)
           write(use_unit,*) '* includes too many v_half tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(89)
         i_bead = mod(v_half_counter,n_beads)
         if (i_bead.eq.0) i_bead = n_beads
         i_atom = int((v_half_counter - i_bead)/n_beads) + 1
         read (89,*) desc_str, ( v_beads_half(i_coord,i_bead,i_atom), i_coord=1,3,1)
       ! v_beads_last
       elseif (desc_str.eq.'v_beads_last') then
         v_last_counter = v_last_counter + 1
         if (v_last_counter .gt. n_atoms*n_beads) then
           write(use_unit,*) '* Warning: The restart file ',trim(PIMD_restart_file)
           write(use_unit,*) '* includes too many v_last tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(89)
         i_bead = mod(v_last_counter,n_beads)
         if (i_bead.eq.0) i_bead = n_beads
         i_atom = int((v_last_counter - i_bead)/n_beads) + 1
         read (89,*) desc_str, ( v_beads_last(i_coord,i_bead,i_atom), i_coord=1,3,1)
       ! r_last
       elseif (desc_str.eq.'r_beads_last') then
         r_last_counter = r_last_counter + 1
         if (r_last_counter .gt. n_atoms*n_beads) then
           write(use_unit,*) '* Warning: The restart file ',trim(PIMD_restart_file)
           write(use_unit,*) '* includes too many r_last tags and can impossibly belong'
           write(use_unit,*) '* to the current system.'
           write(use_unit,*) 'Aborting execution'
           stop
         end if
         backspace(89)
         i_bead = mod(r_last_counter,n_beads)
         if (i_bead.eq.0) i_bead = n_beads
         i_atom = int((r_last_counter - i_bead)/n_beads) + 1
         read (89,*) desc_str, ( r_beads_last(i_coord,i_bead,i_atom), i_coord=1,3,1)
       ! s_NP_beads
       elseif (desc_str.eq.'s_NP_beads') then
         backspace(89)
         read (89,*) desc_str, s_NP_beads
       ! s_NP_beads_half
       elseif (desc_str.eq.'s_NP_beads_half') then
         backspace(89)
         read (89,*) desc_str, s_NP_beads_half
       ! s_NP_beads_last
       elseif (desc_str.eq.'s_NP_beads_last') then
         backspace(89)
         read (89,*) desc_str, s_NP_beads_last
       ! s_dot_NP_beads_half
       elseif (desc_str.eq.'s_dot_NP_beads_half') then
         backspace(89)
         read (89,*) desc_str, s_dot_NP_beads_half
       ! s_dot_NP_beads_last
       elseif (desc_str.eq.'s_dot_NP_beads_last') then
         backspace(89)
         read (89,*) desc_str, s_dot_NP_beads_last
       ! tsystem_beads
       elseif (desc_str.eq.'tsystem_beads') then
         backspace(89)
         read (89,*) desc_str, tsystem_beads
       ! tsystem_beads_last
       elseif (desc_str.eq.'tsystem_beads_last') then
         backspace(89)
         read (89,*) desc_str, tsystem_beads_last
       ! PIMD_H0
       elseif (desc_str.eq.'PIMD_H0') then
         backspace(89)
         read (89,*) desc_str, PIMD_H0
       ! PIMD_Epot_last
       elseif (desc_str.eq.'PIMD_Epot_last') then
         backspace(89)
         read (89,*) desc_str, PIMD_Epot_last
       ! PIMD_stepcount
       elseif (desc_str.eq.'PIMD_stepcount') then
         backspace(89)
         read (89,*) desc_str, PIMD_stepcount
       ! PIMD_force_evaluations
       elseif (desc_str.eq.'PIMD_force_evaluations') then
         backspace(89)
         read (89,*) desc_str, PIMD_force_evaluations
       else
          if (myid.eq.0) then
             write(use_unit,*) "Unknown descriptor ", desc_str, " in file ",trim(PIMD_restart_file),"."
          end if
          stop
       end if
       ! read next line
       read (89,*,iostat = i_code) desc_str
       if (i_code.ne.0) then
          eof = .true.
       end if
     end do
     close(unit=89)

     PIMD_successful_restart_read = .true.
     write(info_str,'(2X,2A)') 'Successfully read all PIMD restart data from file ', trim(PIMD_restart_file)
     call localorb_info(info_str)

     if ((PIMD_ensemble.eq.'NVT_nose-poincare').or.(PIMD_ensemble.eq.'NVT_nose-hoover')) then
        do i_atom = 1, n_atoms
          do i_bead = 1, n_beads
            v_beads_half(:,i_bead,i_atom) = species_m(species(i_atom))*v_beads_half(:,i_bead,i_atom)
            v_beads_last(:,i_bead,i_atom) = species_m(species(i_atom))*v_beads_last(:,i_bead,i_atom)
          end do
        end do
     end if

     if (MD_ensemble.eq.'NVT_nose-hoover') then
        s_NP_beads      = dlog(s_NP_beads)
        s_NP_beads_last = dlog(s_NP_beads_last)
     end if

  end if ! end myid=0

  if (PIMD_time_restart) then
     tsystem_beads      = 0d0
     tsystem_beads_last = 0d0
     PIMD_stepcount = 0
     PIMD_force_evaluations = 0
  end if

  do i_atom = 1, n_atoms
    coords(:,i_atom) = coords_beads(:,1,i_atom)
  end do
 
  call bcast_logical(PIMD_successful_restart_read,0)
  if (PIMD_successful_restart_read) then
     call broadcast_MD_velocities(coords, 0 )
     call broadcast_PIMD_velocities(coords_beads,0)
     call broadcast_PIMD_velocities(v_beads_half,0)
     call broadcast_PIMD_velocities(v_beads_last,0)
     call broadcast_PIMD_velocities(r_beads_last,0)
     call bcast_real(s_NP_beads          ,0)
     call bcast_real(s_NP_beads_last     ,0)
     call bcast_real(s_dot_NP_beads_half ,0)
     call bcast_real(s_dot_NP_beads_last ,0)
     call bcast_real(tsystem_beads       ,0)
     call bcast_real(tsystem_beads_last  ,0)
     call bcast_real(PIMD_H0             ,0)
     call bcast_real(PIMD_Epot_last      ,0)
     call sync_integer(PIMD_stepcount)
     call sync_integer(PIMD_force_evaluations)
  else 
     ! If we were not successful, we should never have reached this point in the code.
     ! Therefore, stop execution.
     write(info_str,'(2X,4A)') 'PIMD restart from file ',trim(PIMD_restart_file),' not successful. ',&
                               'At this point, partial (inconsistent) data may have been read.'
     call aims_stop(info_str)
  end if

end subroutine read_PIMD_restart
