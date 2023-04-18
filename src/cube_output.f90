!****s* FHI-aims/cube_output
!  NAME
!    cube_output
!  SYNOPSIS
subroutine cube_output ( )
!  PURPOSE
!  Wrapper around cube output subroutines to allow for tracking of overall timing 
!  and memory tracking (should that ever become an issue)
!
!  USES
  use runtime_choices, only : use_old_cube_routine
  use timing,          only : get_times, get_timestamps, &
                              tot_time_cube_output, tot_clock_time_cube_output
  use localorb_io,     only : localorb_info
  implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    William Huhn (Duke University)
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    September 2017 - Created 
!  SOURCE
  real*8        :: time_cube_output = 0.0d0, clock_time_cube_output = 0.0d0
  character*100 :: info_str

  character(*), parameter :: func = 'cube_output'

  call get_timestamps( time_cube_output, clock_time_cube_output )
  ! if(n_periodic>1)then
  ! If cube files are requested, setup default values
  ! Must be done after reading geometry
  ! Possibly, a better place would be near read_geometry, so that possible 
  ! input errors resulting in a stop of the code
  ! are caught there
  call setup_cube_defaults ( )
  if (use_old_cube_routine) then
    call output_cube_files_p1 ()
  else
    call output_cube_files_p2 ()
  endif
  call get_times( time_cube_output, clock_time_cube_output, &
                  tot_time_cube_output, tot_clock_time_cube_output, .true. )

  write(info_str, *)  " Cube output                                             :  max(cpu_time)    wall_clock(cpu1)"
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Total time for cube output                            :", &
       tot_time_cube_output, tot_clock_time_cube_output
  call localorb_info( info_str )
  write(info_str, '(A)') "------------------------------------------------------------"
  call localorb_info( info_str )

end subroutine cube_output
!******
