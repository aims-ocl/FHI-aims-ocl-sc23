!****s* FHI-aims/read_plot_dielectric
!  NAME
!   read_plot_dielectric 
!  SYNOPSIS

  subroutine read_plot_dielectric(inputline)
!  USES
    use dimensions
    use runtime_choices
    use localorb_io
!  PURPOSE
!   The subroutine reads the ouput information from control.in for the dielectric output .
!   This subroutine is called from read_control.in 
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
    character*30  :: desc_str
    character*30 :: direction1
    character*30 :: direction2
    character*100 :: info_str
    real*8 ::  broaden_para_temp = 0 
    if ( n_plot_dielectric ==0 ) then 
     ! have not set "output dielectric broaden_type broaden_para direc1 direc 2"
     ! automatically output all the dielectric tensor component and diagonal
     ! part's absorption coefficients. 
     ! The rest are the default values for output dielectric
   
         n_plot_dielectric = 6  ! output all 6 independent dielectric tensor  
      if( .not.  plot_dielectric_allocated)then

         allocate(broaden_type(n_plot_dielectric))
         allocate(broaden_para(n_plot_dielectric))
         allocate(direc_1(n_plot_dielectric))
         allocate(direc_2(n_plot_dielectric))
         plot_dielectric_allocated = .true. 

      end if
      
      !default broadening paramter  Lorentzian 0.1 eV 
      broaden_type(:) = die_broaden_type
      broaden_para(:) = die_broaden_para
      omega_min = die_broaden_para/10.0    ! the min value for the omega are
                                          ! setting to 1/10.0 of the broadening
                                          ! parameter 

      ! calculate all direction momentem matrix (x, y, z) 
      calc_px_dielectric = .true.
      calc_py_dielectric = .true.
      calc_pz_dielectric = .true.
      ! direction settings for all the dielectric tensor component 
      ! x  x direction   
      direc_1(1) = 1   
      direc_2(1) = 1    
      ! y y direction  
      direc_1(2) = 2 
      direc_2(2) = 2
      ! z z direction 
      direc_1(3) = 3
      direc_2(3) = 3
      ! x y direction 
      direc_1(4) = 1 
      direc_2(4) = 2
      ! x z direction 
      direc_1(5) = 1
      direc_2(5) = 3 
      ! y z direction 
      direc_1(6) = 2
      direc_2(6) = 3 
   elseif (flag_dielectric_test) then  
      ! have set "output dielectric broaden_type broaden_para direc1 direc
      ! accordint to the number of  "output dielectric" line, output different
      ! broadening settings for different directions        
      ! the omega value are set according to the smallest broadening parameter         

      if( .not.  plot_dielectric_allocated)then

         allocate(broaden_type(n_plot_dielectric))
         allocate(broaden_para(n_plot_dielectric))
         allocate(direc_1(n_plot_dielectric))
         allocate(direc_2(n_plot_dielectric))
         plot_dielectric_allocated = .true.
         dielectric_initialized = 0

      end if

      dielectric_initialized = dielectric_initialized + 1
   
      write(info_str,*) 
      call localorb_info(info_str, use_unit, '(A)')
      write(info_str,'(A,I5)') '  Plot dielectric direction',dielectric_initialized
      call localorb_info(info_str, use_unit, '(A)')    

      read(inputline,*) desc_str, desc_str, die_method, broaden_para(dielectric_initialized),&
            direction1,direction2
   
      if (die_method == "gaussian") then
         broaden_type(dielectric_initialized) = 1 
      elseif (die_method == "lorentzian") then 
         broaden_type(dielectric_initialized) = 2 
      endif
      if (direction1 == "x") then 
         direc_1(dielectric_initialized) = 1 
         calc_px_dielectric = .true.
      elseif (direction1 == "y") then 
         direc_1(dielectric_initialized) = 2
         calc_py_dielectric = .true.
      elseif (direction1 == "z") then 
         direc_1(dielectric_initialized) = 3 
         calc_pz_dielectric = .true. 
      endif
      if (direction2 == "x") then 
         direc_2(dielectric_initialized) = 1
         calc_px_dielectric = .true.
      elseif (direction2 == "y") then 
         direc_2(dielectric_initialized) = 2
         calc_py_dielectric = .true.
      elseif (direction2 == "z") then 
         direc_2(dielectric_initialized) = 3
         calc_pz_dielectric = .true. 
      endif
      if (dielectric_initialized == 1) then 
         broaden_para_temp = broaden_para(dielectric_initialized)
      elseif (broaden_para_temp > broaden_para(dielectric_initialized)) then 
         broaden_para_temp = broaden_para(dielectric_initialized)
      endif  
      omega_min = broaden_para_temp/10.0 



   endif
   
!    write(info_str,'(A,3F10.6)')  '  | begin', band_begin(band_initialized,:)
!    call localorb_info(info_str, use_unit,'(A)')    
!    write(info_str,'(A,3F10.6)')  '  | end  ', band_end(band_initialized,:)
!    call localorb_info(info_str, use_unit,'(A)')    
!    write(info_str,'(A,I5)')   '  | number of points:', n_points_in_band(band_initialized)
!    call localorb_info(info_str, use_unit,'(A)')    

!!   if (n_points_in_band(band_initialized) .le. 1 ) then
!       if (myid.eq.0) then
!            write(info_str,'(1X,A)')         "* Warning:  "
!            call localorb_info(info_str, use_unit)    
!            write(info_str,'(1X,A,I5,A)')    "* In band number", band_initialized, "   you have requested less than two k-points."
!            call localorb_info(info_str, use_unit)    
!            write(info_str,'(1X,A)')         "* Since every band has at least a starting and an end point,"
!            call localorb_info(info_str, use_unit)    
!            write(info_str,'(1X,A)')         "* Please request two or more points to be calculated."
!            call localorb_info(info_str, use_unit)    
!       end if
!   stop
!   end if 

  end subroutine read_plot_dielectric
!******
