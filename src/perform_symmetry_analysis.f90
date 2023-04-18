! BL 2015
!
! Subroutine checks that all CPUs work on the same lattice and
! stress tensors. It occurs that for some reasons there is a bit flip.

 subroutine perform_symmetry_analysis ()

   use runtime_choices
   use dimensions
   use localorb_io
   use mpi_utilities
   use spglib_symmetry

   implicit none

   character*200 :: info_str

    ! Spacegroup output
    if (n_periodic == 3 .and. use_spglib .and. use_symmetry_analysis) then
      call write_symm_info()
    end if
    if (.not.(use_spglib .and. n_periodic == 3) .and. &
          (use_symmetry_analysis .and. set_symmetry_analysis)) then
      write(info_str,'(2X,A)') &
            "*WARNING: You have asked for a symmetry analysis."
      call localorb_info( info_str )
      write(info_str,'(2X,A)') &
            "*WARNING: However, either you have not compiled spglib"
      call localorb_info( info_str )
      write(info_str,'(2X,A)') &
            "*WARNING: or specified a cluster model."
      call localorb_info( info_str )
      write(info_str,'(2X,A)') &
            "*WARNING: No symmetry output can be provided this time."
      call localorb_info( info_str )
    end if
   
 end subroutine perform_symmetry_analysis
