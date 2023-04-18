!****s* FHI-aims/read_plot_band
!  NAME
!   read_plot_band
!  SYNOPSIS

  subroutine read_plot_band_during_scf(inputline)
!  USES
    use dimensions
    use runtime_choices
    use localorb_io
!  PURPOSE
!   The subroutine reads the ouput information from control.in for the band points.
!   This subroutine is called from read_control.f90 .
!
!   The infrastructure (arrays etc) used in this routine must be different from
!   the original 'output band' version, as both options should work side by side
!   without disturbing one another.
!
!   Notice that, in this case, we do not need to store the number of points requested along
!   each band. The reason is that we are aiming to find all those k-points of the normal
!   k_grid that are located along the requested line in reciprocal space. That number of
!   available k-points is fixed by the specified k_grid , and could even be zero. More
!   commonly, it could be 1 (in Gamma-only calculations ...)
!
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
    logical, save :: plot_band_scf_allocated = .false.
    integer, save :: band_scf_initialized = 0

    if( .not.  plot_band_scf_allocated)then

       allocate(band_scf_begin(n_plot_band_scf,3))
       allocate(band_scf_end(n_plot_band_scf,3))
       plot_band_scf_allocated = .true.
       band_scf_initialized = 0

    end if

    band_scf_initialized = band_scf_initialized + 1

    write(info_str,*) 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(A,I5)') '# Plot band during s.c.f.:',band_scf_initialized
    call localorb_info(info_str,use_unit,'(A)')    

    read(inputline,*) desc_str,desc_str,  band_scf_begin(band_scf_initialized,1), band_scf_begin(band_scf_initialized,2), &
         band_scf_begin(band_scf_initialized,3),&
         band_scf_end(band_scf_initialized,1),   band_scf_end(band_scf_initialized,2), band_scf_end(band_scf_initialized,3)

    write(info_str,'(A,3F10.6)')  '  | begin', band_scf_begin(band_scf_initialized,:)
    call localorb_info(info_str,use_unit,'(A)')    
    write(info_str,'(A,3F10.6)')  '  | end  ', band_scf_end(band_scf_initialized,:)
    call localorb_info(info_str,use_unit,'(A)')    

  end subroutine read_plot_band_during_scf
!******
