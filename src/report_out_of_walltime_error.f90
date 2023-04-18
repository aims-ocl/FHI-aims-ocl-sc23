!****s* FHI-aims/report_out_of_walltime_error
!  NAME
!    report_out_of_walltime_error
!  SYNOPSIS
        subroutine report_out_of_walltime_error()
!  PURPOSE
!    handles the warning output and gives the final atomic structure for the case where the system ran out of wall time
!  USES
        use geometry
        use species_data
        use dimensions
        use constants
        use runtime_choices
        use localorb_io
        implicit none
!  ARGUMENTS
!    none
!  INPUT
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
          integer :: i_atom, i_coord
          character*150 :: info_str

          write(info_str,'(A)') '*** WARNING: FHI-aims is terminating due to walltime restrictions'
          call localorb_info( info_str )
          ! did we request relaxation along the way???
          if (relax_mode.ne.RELAX_OFF) then
             write(info_str,'(A)') '*** Last used atomic structure (not necessarily the equilibrium result):'
             call localorb_info( info_str )
             write(info_str,'(A,15X,A,13X,A,13X,A,7X,A)') &
                  "***","x [A]","y [A]","z [A]", "Atom"
             call localorb_info( info_str )
             do i_atom = 1, n_atoms, 1
                write(info_str,'(A,2X,A,3(2X,F16.8),2X,A)') &
                     "***","atom ", &
                     (coords(i_coord,i_atom)*bohr, i_coord=1,3,1), &
                     species_name(species(i_atom))
                call localorb_info( info_str )
             enddo
          end if
        end subroutine report_out_of_walltime_error
!******
