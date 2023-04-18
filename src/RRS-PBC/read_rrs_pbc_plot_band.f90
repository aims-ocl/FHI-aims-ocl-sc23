!****s* FHI-aims/read_rrs_pbc_plot_band
!  NAME
!   read_rrs_pbc_plot_band
!  SYNOPSIS

  subroutine read_rrs_pbc_plot_band(inputline)
!  USES
    use dimensions
    use runtime_choices
    use localorb_io
    use mpi_tasks
!  PURPOSE
!   The subroutine reads the ouput information from control.in for the band poitns.
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

    character(*),intent(in) :: inputline
    character*30 desc_str
    character*100 :: info_str
    logical, save :: rrs_pbc_plot_band_allocated = .false.
    integer, save :: rrs_pbc_band_initialized = 0

    !if( .not. rrs_pbc_plot_band_allocated )then
    if( .not. allocated(rrs_pbc_band_begin) )then
       allocate(rrs_pbc_band_begin(rrs_pbc_n_plot_band,3))
       allocate(rrs_pbc_band_end(rrs_pbc_n_plot_band,3))
       allocate(rrs_pbc_n_points_in_band(rrs_pbc_n_plot_band ))
       rrs_pbc_band_initialized = 0
    end if

    rrs_pbc_band_initialized = rrs_pbc_band_initialized + 1

    write(info_str,*) 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(A,I5)') '  Plot RRS-PBC band',rrs_pbc_band_initialized
    call localorb_info(info_str,use_unit,'(A)')    


    read(inputline,*) desc_str,desc_str,  rrs_pbc_band_begin(rrs_pbc_band_initialized,1), &
         rrs_pbc_band_begin(rrs_pbc_band_initialized,2), &
         rrs_pbc_band_begin(rrs_pbc_band_initialized,3),&
         rrs_pbc_band_end(rrs_pbc_band_initialized,1), &
         rrs_pbc_band_end(rrs_pbc_band_initialized,2), &
         rrs_pbc_band_end(rrs_pbc_band_initialized,3), &
         rrs_pbc_n_points_in_band(rrs_pbc_band_initialized )


    write(info_str,'(A,3F10.6)')  '  | begin', rrs_pbc_band_begin(rrs_pbc_band_initialized,:)
    call localorb_info(info_str,use_unit,'(A)')    
    write(info_str,'(A,3F10.6)')  '  | end  ', rrs_pbc_band_end(rrs_pbc_band_initialized,:)
    call localorb_info(info_str,use_unit,'(A)')    
    write(info_str,'(A,I5)')   '  | number of points:', rrs_pbc_n_points_in_band(rrs_pbc_band_initialized)
    call localorb_info(info_str,use_unit,'(A)')    

   if (rrs_pbc_n_points_in_band(rrs_pbc_band_initialized) .le. 1 ) then
       if (myid.eq.0) then
            write(use_unit,'(1X,A)')         "* Warning:  "
            write(use_unit,'(1X,A,I5,A)')    "* In band number", rrs_pbc_band_initialized, &
                                      " you have requested less than two k-points."
            write(use_unit,'(1X,A)')         "* Since every band has at least a starting and an end point,"
            write(use_unit,'(1X,A)')         "* Please request two or more points to be calculated."
       end if
   stop
   end if 

  end subroutine read_rrs_pbc_plot_band
!******

