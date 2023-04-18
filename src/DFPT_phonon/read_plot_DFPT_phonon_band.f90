!****s* FHI-aims/DFPT_phonon/read_plot_DFPT_phonon_band
!  NAME
!   read_plot_DFPT_phonon_band
!  SYNOPSIS

  subroutine read_plot_DFPT_phonon_band(inputline)
!  USES
    use dimensions, only : n_plot_band_DFPT_phonon
    use localorb_io
    use mpi_tasks, only: aims_stop_coll, myid
    use runtime_choices, only : DFPT_phonon_band_begin, DFPT_phonon_band_end, n_points_in_band_DFPT_phonon
!  PURPOSE
!   The subroutine reads the ouput information from control.in for the DFTP_phonon_band points.
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
    logical, save :: plot_DFPT_phonon_band_allocated = .false.
    integer ,save :: i_band = 0

    if( .not.  plot_DFPT_phonon_band_allocated)then
      allocate(DFPT_phonon_band_begin(n_plot_band_DFPT_phonon,3))
      allocate(DFPT_phonon_band_end(n_plot_band_DFPT_phonon,3))
      allocate(n_points_in_band_DFPT_phonon(n_plot_band_DFPT_phonon ))

      plot_DFPT_phonon_band_allocated = .true.
      i_band = 0
    endif

    i_band = i_band + 1

    write(info_str,*)
    call localorb_info(info_str, use_unit, '(A)')
    write(info_str,'(A,I5)') '  Plot DFPT_phonon_band', i_band
    call localorb_info(info_str, use_unit, '(A)')

    read(inputline,*) desc_str, desc_str, DFPT_phonon_band_begin(i_band,1), &
                      DFPT_phonon_band_begin(i_band,2), &
                      DFPT_phonon_band_begin(i_band,3), &
                      DFPT_phonon_band_end(i_band,1),   &
                      DFPT_phonon_band_end(i_band,2),   &
                      DFPT_phonon_band_end(i_band,3),   &
                      n_points_in_band_DFPT_phonon(i_band )


    write(info_str,'(A,3F10.6)')  '  | begin', DFPT_phonon_band_begin(i_band,:)
    call localorb_info(info_str, use_unit,'(A)')
    write(info_str,'(A,3F10.6)')  '  | end  ', DFPT_phonon_band_end(i_band,:)
    call localorb_info(info_str, use_unit,'(A)')
    write(info_str,'(A,I5)')   '  | number of points:', n_points_in_band_DFPT_phonon(i_band)
    call localorb_info(info_str, use_unit,'(A)')

   if (n_points_in_band_DFPT_phonon(i_band) .le. 1 ) then
       if (myid.eq.0) then
            write(info_str,'(1X,A)')         "* Warning:  "
            call localorb_info(info_str, use_unit)
            write(info_str,'(1X,A,I5,A)')    "* In band number", i_band, "   you have requested less than two k-points."
            call localorb_info(info_str, use_unit)
            write(info_str,'(1X,A)')         "* Since every band has at least a starting and an end point,"
            call localorb_info(info_str, use_unit)
            write(info_str,'(1X,A)')         "* Please request two or more points to be calculated."
            call localorb_info(info_str, use_unit)
       end if
    call aims_stop_coll("More q point please ",'DFPT_phonon/read_plot_DFPT_phonon_band')
   end if

  end subroutine read_plot_DFPT_phonon_band
!******
