!****s* FHI-aims/read_plot_band
!  NAME
!   read_plot_band
!  SYNOPSIS

  subroutine read_plot_band(inputline)
!  USES
    use dimensions
    use localorb_io
    use mpi_tasks, only: myid
    use runtime_choices
!  PURPOSE
!   The subroutine reads the ouput information from control.in for the band points.
!   This subroutine is called from read_control
!  INPUTS
!    none
!  OUTPUT
!    none
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

    implicit none

    character(*),intent(in) :: inputline
    character*30 desc_str
    character*100 :: info_str

    if( .not.  plot_band_allocated)then

       allocate(band_begin(n_plot_band,3))
       allocate(band_end(n_plot_band,3))
       allocate(n_points_in_band(n_plot_band ))
       plot_band_allocated = .true.
       band_initialized = 0

    end if

    band_initialized = band_initialized + 1

    write(info_str,*)
    call localorb_info(info_str, use_unit, '(A)')
    write(info_str,'(A,I5)') '  Plot band',band_initialized
    call localorb_info(info_str, use_unit, '(A)')


    read(inputline,*) desc_str, desc_str, band_begin(band_initialized,1), band_begin(band_initialized,2), &
         band_begin(band_initialized,3),&
         band_end(band_initialized,1),   band_end(band_initialized,2), band_end(band_initialized,3), &
         n_points_in_band(band_initialized )





    write(info_str,'(A,3F10.6)')  '  | begin', band_begin(band_initialized,:)
    call localorb_info(info_str, use_unit,'(A)')
    write(info_str,'(A,3F10.6)')  '  | end  ', band_end(band_initialized,:)
    call localorb_info(info_str, use_unit,'(A)')
    write(info_str,'(A,I5)')   '  | number of points:', n_points_in_band(band_initialized)
    call localorb_info(info_str, use_unit,'(A)')

   if (n_points_in_band(band_initialized) .le. 1 ) then
       if (myid.eq.0) then
            write(info_str,'(1X,A)')         "* Warning:  "
            call localorb_info(info_str, use_unit)
            write(info_str,'(1X,A,I5,A)')    "* In band number", band_initialized, "   you have requested less than two k-points."
            call localorb_info(info_str, use_unit)
            write(info_str,'(1X,A)')         "* Since every band has at least a starting and an end point,"
            call localorb_info(info_str, use_unit)
            write(info_str,'(1X,A)')         "* Please request two or more points to be calculated."
            call localorb_info(info_str, use_unit)
       end if
   stop
   end if

  end subroutine read_plot_band
!******
