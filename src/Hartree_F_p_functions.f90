!****h* FHI-aims/Hartree_F_p_functions
!  NAME
!    Hartree_F_p_functions
!  SYNOPSIS

module Hartree_F_p_functions

  !  PURPOSE
  !  The module includes fonctions for so called Fp functions. The functions are used
  !  in periodic systems for the far distance Hartree potential.
  !
  !  The method and the equations are taken from :
  !  Bernard Delley, J. Phys. Chem 1996, 100, 6107-6110.
  !
  !  USES

  implicit none
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

  real*8,dimension(0:21),private:: Kp
  real*8, dimension(0:17) :: F_1_div_r_coeff

  real*8, dimension(:,:,:), allocatable :: Fp_function_spline   ! (1:lmax,n_max_spline,Fp_max_grid)
  real*8, dimension(:,:,:), allocatable :: Fpc_function_spline   ! (1:lmax,n_max_spline,Fp_max_grid)

  ! not cleaned up
  real*8 ::Ewald_radius_to(11), inv_Ewald_radius_to(2)
  real*8:: P_erfc_4(6), P_erfc_5(7), P_erfc_6(8)


  real*8 :: Fp_grid_min, Fp_grid_inc, Fp_grid_max
  integer :: Fp_max_grid, lmax_Fp

  !******


contains
  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/initialize_F_p_functions
  !  NAME
  !   initialize_F_p_functions
  !  SYNOPSIS

  subroutine initialize_F_p_functions

    !  PURPOSE
    !  The subroutine initializes the Hartree_F_p_functions module variables.
    !
    !  USES

    use dimensions, only: l_pot_max, n_max_spline, n_periodic, n_species
    use runtime_choices, only: Ewald_radius, use_hartree_non_periodic_ewald
    use constants, only: sqrt_pi
    use species_data, only: l_hartree
    use grids, only: r_grid_min, r_grid_max, r_grid_inc, invert_log_grid
    use spline, only: cubic_spline
    use localorb_io, only : use_unit

    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    implicit none
    integer:: i, i_species, i_grid, i_l, n_grid_switch
    real*8 :: factorial, rlog
    real*8, dimension(:,:), allocatable :: Fp_values, Fpc_values, aux_spl

!!$    ! DEBUG:
!!$    real*8, dimension(20) :: orig, splined

    do i=1,11
       Ewald_radius_to(i) = Ewald_radius**i
    end do
    do i=1,2
       inv_Ewald_radius_to(i) = 1.d0/ (Ewald_radius**i)
    end do

    P_erfc_4(1)= 210.0d0 * Ewald_radius_to(6)
    P_erfc_4(2)= 140.0d0 * Ewald_radius_to(4) 
    P_erfc_4(3)= 56.0d0 * Ewald_radius_to(2) 
    P_erfc_4(4)= 16.0d0 
    P_erfc_4(5)= 105.0d0 * Ewald_radius_to(7) *sqrt_pi
    P_erfc_4(6)= Ewald_radius_to(7)  * sqrt_pi


    P_erfc_5(1) = 2* 945.0d0 * Ewald_radius_to(8)
    P_erfc_5(2) = 2* 630.0d0 * Ewald_radius_to(6) 
    P_erfc_5(3) = 2*252.0d0 * Ewald_radius_to(4)
    P_erfc_5(4) = 2*72.0d0 * Ewald_radius_to(2) 
    P_erfc_5(5) = 2*16.0d0 
    P_erfc_5(6) = 945.0d0 * Ewald_radius_to(9) * sqrt_pi 
    P_erfc_5(7) = Ewald_radius_to(9) * sqrt_pi

    P_erfc_6(1) = 20790.0d0 * Ewald_radius_to(10)
    P_erfc_6(2) = 13860.0d0 * Ewald_radius_to(8)
    P_erfc_6(3) = 5544.0d0 * Ewald_radius_to(6)
    P_erfc_6(4) =  1584.0d0 * Ewald_radius_to(4)
    P_erfc_6(5) =  352.0d0 *  Ewald_radius_to(2)
    P_erfc_6(6) =  64.0d0
    P_erfc_6(7) = 10395.0d0  * Ewald_radius_to(11) * sqrt_pi
    P_erfc_6(8) = Ewald_radius_to(11)  * sqrt_pi



    Kp(0) = 2/sqrt_pi/(Ewald_radius)
    factorial = 1
    do i=1,20

       factorial = factorial*i

       Kp(i) = 2/sqrt_pi *(-1)**i/ (factorial * (2*i+1)) /(Ewald_radius)**(2*i+1)

    end do

    if(l_pot_max>17)then
       write(use_unit,*)'Error: Fp functions are not implemented so high angular l values!'
       write(use_unit,*)'Maximum is now 17.'
       stop
    end if


    ! set F_1_div_r_coeff lookup table
    F_1_div_r_coeff(0)  =  1.0d0
    F_1_div_r_coeff(1)  = -1.0d0
    F_1_div_r_coeff(2)  =  3.0d0
    F_1_div_r_coeff(3)  = -15.0d0
    F_1_div_r_coeff(4)  =  105.0d0
    F_1_div_r_coeff(5)  = -945.0d0
    F_1_div_r_coeff(6)  =  10395.0d0
    F_1_div_r_coeff(7)  = -135135.0d0
    F_1_div_r_coeff(8)  =  2027025.0d0
    F_1_div_r_coeff(9)  = -34459425.0d0
    F_1_div_r_coeff(10) =  654729075.0d0
    F_1_div_r_coeff(11) = -13749310575.0d0
    F_1_div_r_coeff(12) =  316234143225.0d0
    F_1_div_r_coeff(13) = -7905853580625.0d0
    F_1_div_r_coeff(14) =  213458046676875.0d0
    F_1_div_r_coeff(15) = -6190283353629375.0d0
    F_1_div_r_coeff(16) =  191898783962510625.0d0
    F_1_div_r_coeff(17) = -6332659870762850625.0d0

    ! allocate and calculate spline array for (periodic) Fp_function_spline 
    if (  n_periodic > 0  .or.  use_hartree_non_periodic_ewald  ) then

       ! determine maximal l-value & the maximal logarithmic grid settings
       lmax_Fp = 0
       Fp_grid_min =  1d6
       Fp_grid_max = -1d6
       Fp_grid_inc =  1d6
       do i_species = 1, n_species
          lmax_Fp     = max(l_hartree(i_species)+1,lmax_Fp) ! can't be more than l_hartree + 1 
          Fp_grid_min = min(r_grid_min(i_species),Fp_grid_min)
          Fp_grid_max = max(r_grid_max(i_species),Fp_grid_max)
          Fp_grid_inc = min(r_grid_inc(i_species),Fp_grid_inc)
       end do
       
       ! calculate the number of grid poins in the logarithmic Fp-grid, borrowed from subroutine parse_species in dimension.f90
       Fp_max_grid = int ( log(Fp_grid_max/Fp_grid_min) / log(Fp_grid_inc) ) + 1
       
       if (.not.allocated(Fp_function_spline))  allocate(Fp_function_spline (0:lmax_Fp,n_max_spline,Fp_max_grid))
       if (.not.allocated(Fpc_function_spline)) allocate(Fpc_function_spline(0:lmax_Fp,n_max_spline,Fp_max_grid))
       if (.not.allocated(Fp_values))           allocate(Fp_values (0:lmax_Fp,Fp_max_grid))
       if (.not.allocated(Fpc_values))          allocate(Fpc_values(0:lmax_Fp,Fp_max_grid))
       if (.not.allocated(aux_spl))             allocate(aux_spl(n_max_spline,Fp_max_grid))
       
       ! calculate the splines for all species & all grid points 
       rlog = Fp_grid_min
       do i_grid = 1, Fp_max_grid
          ! calculate all Fp functions for this particular radius
          call F_erf_original (Fp_values (:,i_grid),rlog,lmax_Fp)
          call F_erfc_original(Fpc_values(:,i_grid),rlog,lmax_Fp)
          rlog = rlog * Fp_grid_inc ! next radius, please
       end do

       ! calculate splines for all of these functions
       ! FIXME: The current splines appear to be accurate to within 10^-10, if this is not
       !        enough, one can think of using either a denser log grid or hermite splines. 
       !        tried the latter, but I could not figure out the proper derivatives in the 
       !        time I had available. (FH, 04/2009)
       Fp_function_spline  = 0d0
       Fpc_function_spline = 0d0 
       n_grid_switch       = invert_log_grid(25d0*Ewald_radius,Fp_grid_min,Fp_grid_inc) + 1
       n_grid_switch       = min(n_grid_switch,Fp_max_grid)

       do i_l = 0, lmax_Fp
          call cubic_spline(Fp_values(i_l,:),Fp_max_grid,aux_spl)
          do i_grid = 1, Fp_max_grid
             Fp_function_spline(i_l,:,i_grid) = aux_spl(:,i_grid)
          end do
          call cubic_spline(Fpc_values(i_l,:),Fp_max_grid,aux_spl)
          do i_grid = 1, n_grid_switch
             Fpc_function_spline(i_l,:,i_grid) = aux_spl(:,i_grid)
          end do
       end do
       
       ! deallocate stuff that was only needed for this particular initialization
       if (allocated(Fp_values )) deallocate(Fp_values )
       if (allocated(Fpc_values)) deallocate(Fpc_values)
       if (allocated(aux_spl))    deallocate(aux_spl   )


!!$       ! DEBUG: output the comparison between initial and splined data
!!$       do rlog = 1d-2, 100d0, 1d-1
!!$          call F_erf_original(splined,rlog,14)
!!$!          call F_erf_spline(splined,rlog,14)
!!$          call F_erfc_original(orig,rlog,14)
!!$          if (myid.eq.0)write(use_unit, '(A,3E18.8E4)') "DEBUG: ", rlog,orig(8),splined(8)
!!$       end do

    end if   ! splining of Fp functions!

  end subroutine initialize_F_p_functions
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/clean_F_p_functions
  !  NAME
  !   clean_F_p_functions
  !  SYNOPSIS
  subroutine clean_F_p_functions
    !  PURPOSE
    !  Deallocate data
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    implicit none
    if (allocated(Fp_function_spline))  deallocate(Fp_function_spline)
    if (allocated(Fpc_function_spline)) deallocate(Fpc_function_spline)
  end subroutine clean_F_p_functions
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/F_erf
  !  NAME
  !   F_erf
  !  SYNOPSIS
  subroutine F_erf(F,r,p)
    ! PURPOSE
    !    provides an interface to two versions of Delley's Fp functions (1996)
    ! USES
    use runtime_choices, only: hartree_fp_function_splines
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    integer :: p
    real*8, dimension(0:p) :: F
    real*8 :: r
    if (hartree_fp_function_splines) then
       call F_erf_spline(F,r,p)
    else
       call F_erf_original(F,r,p)
    end if
  end subroutine F_erf
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/F_erfc
  !  NAME
  !   F_erfc
  !  SYNOPSIS
  subroutine F_erfc(F,r,p)
    ! PURPOSE
    !   provides splined/nonsplined interfaces to the Hartree Fp functions
    ! USES
    use runtime_choices, only: hartree_fp_function_splines
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    integer :: p
    real*8, dimension(0:p) :: F
    real*8 :: r
    if (hartree_fp_function_splines) then
       call F_erfc_spline(F,r,p)
    else
       call F_erfc_original(F,r,p)
    end if
  end subroutine F_erfc
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/F_erf_original 
  !  NAME
  !   F_erf_original 
  !  SYNOPSIS
  subroutine F_erf_original(F,r,p)
    !  PURPOSE
    !    The original (direct) calculation of the Fp functions
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    integer :: p
    real*8, dimension(0:p) :: F
    real*8 :: r!, F_erf_single
!    do l = 0, p
       call F_erf_table_original(F(0:p), r,p)
!    end do
  end subroutine F_erf_original
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/F_erfc_original
  !  NAME
  !    F_erfc_original
  !  SYNOPSIS
  subroutine F_erfc_original(F,r,p)
    !  PURPOSE
    !    original calculation of the Fpc functions
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    integer :: p
    real*8, dimension(0:p) :: F
    real*8 :: r!, F_erfc_single
!    do l = 0, p
    call F_erfc_table_original( F(0:p), r,p)
!    end do
  end subroutine F_erfc_original
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/F_erf_spline
  !  NAME
  !   F_erf_spline
  !  SYNOPSIS
  subroutine F_erf_spline(F,r,p)
    !  PURPOSE
    !    provides vectorizable spline call to Hartree Fp functions
    !  USES
    use grids, only: invert_log_grid
    use spline, only: spline_vector
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    integer :: p
    real*8, dimension(0:p) :: F
    real*8 :: r, rlog
    ! calculate the grid index on the logarithmic grid
    rlog = invert_log_grid(r,Fp_grid_min,Fp_grid_inc)
    ! call cubic spline routine
    call spline_vector(rlog,Fp_function_spline,Fp_max_grid,lmax_Fp+1,Fp_max_grid,p+1,F)
  end subroutine F_erf_spline
  !******

  !--------------------------------------------------------------
  !****s* Hartree_F_p_functions/F_erfc_spline
  !  NAME
  !   F_erfc_spline
  !  SYNOPSIS
  subroutine F_erfc_spline(F,r,p)
    ! PURPOSE
    !   spline call for the hartree fpc functions
    ! USES
    use grids, only: invert_log_grid
    use spline, only: spline_vector
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    integer :: p
    real*8, dimension(0:p) :: F
    real*8 :: r, rlog
    ! calculate the grid index on the logarithmic grid
    rlog = invert_log_grid(r,Fp_grid_min,Fp_grid_inc)
    ! call cubic spline routine
    call spline_vector(rlog,Fpc_function_spline,Fp_max_grid,lmax_Fp+1,Fp_max_grid,p+1,F)
  end subroutine F_erfc_spline
  !******

  !--------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erf_single_spline
  !  NAME
  !   F_erf_single_spline
  !  SYNOPSIS
  function F_erfc_single_spline(r,p)
    ! PURPOSE
    !  single splined F-p-function
    ! USES
    use grids, only: invert_log_grid
    use spline, only: val_spline
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    real*8 :: F_erfc_single_spline, r, rlog
    integer :: p
    rlog = invert_log_grid(r,Fp_grid_min,Fp_grid_inc)
    F_erfc_single_spline = val_spline(rlog,Fpc_function_spline(p,:,:),Fp_max_grid)
  end function F_erfc_single_spline
  !******

  !--------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erfc_single_spline
  !  NAME
  !   F_erfc_single_spline
  !  SYNOPSIS
  function F_erf_single_spline(r,p)
    ! PURPOSE
    !  single splined F-p-function
    ! USES
    use grids, only: invert_log_grid
    use spline, only: val_spline
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    real*8 :: F_erf_single_spline, r, rlog
    integer :: p
    rlog = invert_log_grid(r,Fp_grid_min,Fp_grid_inc)
    F_erf_single_spline = val_spline(rlog,Fp_function_spline(p,:,:),Fp_max_grid)
  end function F_erf_single_spline
  !******

  !--------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erf_single
  !  NAME
  !   F_erf_single
  !  SYNOPSIS
  function F_erf_single(r,p)
    ! PURPOSE
    !  interface between splined and non-splined functions for the calculation
    !  of a single F-p-function
    ! USES
    use runtime_choices, only: hartree_fp_function_splines
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    real*8  :: F_erf_single, r
    integer :: p
    if (hartree_fp_function_splines) then
       F_erf_single = F_erf_single_spline(r,p)
    else
       F_erf_single = F_erf_single_original(r,p)
    end if    
  end function F_erf_single
  !******

  !--------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erfc_single
  !  NAME
  !   F_erfc_single
  !  SYNOPSIS
  function F_erfc_single(r,p)
    ! PURPOSE
    !  interface between splined and non-splined functions for the calculation
    !  of a single F-p-function
    ! USES
    use runtime_choices, only: hartree_fp_function_splines
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    implicit none
    real*8  :: F_erfc_single, r
    integer :: p
    if (hartree_fp_function_splines) then
       F_erfc_single = F_erfc_single_spline(r,p)
    else
       F_erfc_single = F_erfc_single_original(r,p)
    end if        
  end function F_erfc_single
  !******

  !--------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_1_div_r
  !  NAME
  !   F_1_div_r
  !  SYNOPSIS

  function F_1_div_r( r, p )

    !  PURPOSE
    !  The function calculates Fp function for the cluster systems.
    !  This function is replaced by direct and farter loops in module hartree_potential_real_p0.
    !  However it is also usefull to have this for testing purposes.
    !
    implicit none
    !  ARGUMENTS

    real*8:: r,F_1_div_r
    integer:: p

    !  INPUTS
    !    o r -- distance
    !    o p -- order of function
    !  OUTPUT
    !    o F_1_div_r -- value of Fp function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    select case(p)

    case(0)

       F_1_div_r = 1/r

    case(1)

       F_1_div_r = 1/r**3

    case(2)

       F_1_div_r = 1/r**5

    case(3)

       F_1_div_r = 1/r**7

    case(4)

       F_1_div_r = 1/r**9

    case(5)

       F_1_div_r = 1/r**11

    case(6)

       F_1_div_r = 1/r**13

    case(7)

       F_1_div_r = 1/r**15

    case(8)

       F_1_div_r = 1/r**17

    case(9)

       F_1_div_r = 1/r**19

    case(10)

       F_1_div_r = 1/r**21

    case(11)

       F_1_div_r = 1/r**23

    case(12)

       F_1_div_r = 1/r**25

    case(13)

       F_1_div_r = 1/r**27

    case(14)

       F_1_div_r = 1/r**29

    case(15)

       F_1_div_r = 1/r**31

    case(16)

       F_1_div_r = 1/r**33

    case(17)

       F_1_div_r = 1/r**35

!!$    case(18)
!!$       
!!$       F_1_div_r = 221643095476699771875/r**37
!!$
!!$    case(19)
!!$       
!!$       F_1_div_r = -8200794532637891559375/r**39
!!$
!!$    case(20)
!!$       
!!$       F_1_div_r = 319830986772877770815625/r**41

    end select



  end function F_1_div_r
  !******
  !-------------------------------------------------------
  !****f* Hartree_F_p_functions/F_1_div_r_2
  !  NAME
  !   F_1_div_r_2
  !  SYNOPSIS

  function F_1_div_r_2( r, p )

    !  PURPOSE
    !  The function calculates Fp function for the cluster systems.
    !  This function is replaced by direct and farter loops in module hartree_potential_real_p0.
    !  However it is also usefull to have this for testing purposes.
    !  Note: this is updated version of F_1_div_r
    !
    implicit none
    !  ARGUMENTS

    real*8:: r,F_1_div_r_2
    integer:: p

    !  INPUTS
    !    o r -- distance
    !    o p -- order of function
    !  OUTPUT
    !    o F_1_div_r_2 -- value of Fp function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE







    select case(p)

    case(0)

       F_1_div_r_2 = 1.0d0/r

    case(1)

       F_1_div_r_2 = -1.0d0/r**3

    case(2)

       F_1_div_r_2 = 3.0d0/r**5

    case(3)

       F_1_div_r_2 = -15.0d0/r**7

    case(4)

       F_1_div_r_2 = 105.0d0/r**9

    case(5)

       F_1_div_r_2 = -945.0d0/r**11

    case(6)

       !       F_1_div_r_2 = 10395/r**13
       F_1_div_r_2 = ( 2.036978755158300/r) **13

    case(7)

       F_1_div_r_2 = -135135.0d0/r**15

    case(8)

       F_1_div_r_2 = 2027025.0d0/r**17

    case(9)

       F_1_div_r_2 = -34459425.0d0/r**19

    case(10)

       F_1_div_r_2 = 654729075.0d0/r**21

    case(11)

       F_1_div_r_2 = -13749310575.0d0/r**23

    case(12)

       F_1_div_r_2 = 316234143225.0d0/r**25

    case(13)

       F_1_div_r_2 = -7905853580625.0d0/r**27

    case(14)

       F_1_div_r_2 = 213458046676875.0d0/r**29

    case(15)

       F_1_div_r_2 = -6190283353629375.0d0/r**31

    case(16)

       F_1_div_r_2 = 191898783962510625.0d0/r**33

    case(17)

       F_1_div_r_2 = -6332659870762850625.0d0/r**35

!!$    case(18)
!!$       
!!$       F_1_div_r_2 = 221643095476699771875/r**37
!!$
!!$    case(19)
!!$       
!!$       F_1_div_r_2 = -8200794532637891559375/r**39
!!$
!!$    case(20)
!!$       
!!$       F_1_div_r_2 = 319830986772877770815625/r**41

    end select



  end function F_1_div_r_2

  !---------------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erfc_table_original
  !  NAME
  !   F_erfc_single_original
  !  SYNOPSIS

  subroutine F_erfc_table_original( F_table, r,p_max)

    !  PURPOSE
    !  The function calculates Fp erfc function for the periodic systems.
    !  This is Fp function for regions outside the multipole radius.
    !
    use arch_specific, only: arch_erfc
    use constants, only: sqrt_pi
    use runtime_choices, only : Ewald_radius
    implicit none
    !  ARGUMENTS
    
    integer :: p_max
    real*8:: F_table(0:p_max)
    real*8::  r


    !  INPUTS
    !    o r -- distance
    !    o p -- order of function
    !  OUTPUT
    !    o F_1_div_r -- value of Fp function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    !    real*8:: n=100
    integer:: p
    real*8:: rto2, rto(0:p_max)
    real*8:: erfc_r, exp_m, exp_p
   
    real*8, external :: ddot


    F_table  = 0.0d0
!    return
    if(r > 25*Ewald_radius)then
       return
    end if

    rto2 = r*r

    rto(0) = r
    do p=1,p_max,1
       rto(p) = rto(p-1)*rto2
    end do


    do p=0,p_max,1


       select case(p)

       case(0)

          erfc_r = arch_erfc(r*inv_Ewald_radius_to(1))

          F_table(p)  = erfc_r /r

       case(1)

!          F_table(p) =  - (((2.0d0 * exp(- (r**2 /Ewald_radius**2 )) * r ) & 
!               / (Ewald_radius * sqrt_pi ) + Arch_erfc(r /Ewald_radius) ) /r**3 )

          exp_m = exp(- (rto2 *inv_Ewald_radius_to(2) ))

          F_table(p) =  - (((2.0d0 * exp_m * r ) & 
               / (Ewald_radius * sqrt_pi ) + erfc_r ) /rto(1) )

       case(2)


!          F_table(p) = &
!               ( ( exp(-(r**2 /Ewald_radius**2 ) ) * ((6 * Ewald_radius**2 * r + 4 * r**3 + 3 *  exp(r**2 /Ewald_radius**2 ) * &
!               Ewald_radius**3 * sqrt_pi * arch_erfc(r/Ewald_radius)) ) ) / (Ewald_radius**3 * sqrt_pi* r**5 ) ) 

          exp_p = 1.d0/exp_m

          F_table(p) = &
               ( ( exp_m * ((6 * Ewald_radius_to(2) * r + 4 * rto(1) + 3 * exp_p * &
               Ewald_radius_to(3) * sqrt_pi * erfc_r  ) ) ) / (Ewald_radius_to(3) * sqrt_pi* rto(2) ) ) 


       case(3)

!          F_table(p) = &
!               -(( exp (- (r**2 /Ewald_radius**2 ) )  * ((30 * Ewald_radius**4 * r +  20 * Ewald_radius_to(2) * r**3 + 8 * r**5 +  &
!               15  * exp(r**2 /Ewald_radius**2 ) * Ewald_radius**5  * sqrt_pi * arch_erfc( r /Ewald_radius)) ) ) &
!               / (Ewald_radius**5 *  sqrt_pi * r**7 ))

          F_table(p) = &
               -(( exp_m  * ((30 * Ewald_radius_to(4) * r +  20 * Ewald_radius_to(2) * rto(1) + 8 * rto(2) +  &
               15  * exp_p * Ewald_radius_to(5)  * sqrt_pi * erfc_r    ) ) ) &
               / (Ewald_radius_to(5) *  sqrt_pi * rto(3) ))


       case(4)

!          F_table(p) = &
!               ( ( exp(- (r**2 /Ewald_radius**2 ) ) * ((210.0d0 * Ewald_radius**6 * r + 140.0d0 * Ewald_radius**4 * r**3 &  
!               + 56.0d0 * Ewald_radius**2 * r**5 + 16.0d0 * r**7 + 105.0d0 * exp(r**2 /Ewald_radius**2 ) * Ewald_radius**7 &
!               *sqrt_pi* arch_erfc(r/Ewald_radius)) ) ) / (Ewald_radius**7  * sqrt_pi * r**9 ) )

           F_table(p) = ( exp_m  * ( &
                ddot(4,P_erfc_4(1:4),1, rto(0:3),1) +  P_erfc_4(5) * exp_p * erfc_r )) & 
                / ( P_erfc_4(6) * rto(4) ) 


       case(5)

!          F_table(p) = &
!               ( (- ( ( exp(- (r**2 /Ewald_radius**2 )) *  ((2.0d0 *  ((945.0d0 * Ewald_radius**8 * r + &
!               630.0d0 * Ewald_radius**6 * r**3 + 252.0d0 * Ewald_radius**4 * r**5 + 72.0d0 * Ewald_radius**2 * r**7 + & 
!               16.0d0 * r**9) ) +  945.0d0 * exp(r**2 /Ewald_radius**2) * &
!               Ewald_radius**9 * sqrt_pi * arch_erfc( r /Ewald_radius)))) / (Ewald_radius**9 * sqrt_pi * r**11 ) ) ) )

          F_table(p) =  - ( exp_m  * &
               ( ddot(5,P_erfc_5(1:5),1, rto(0:4),1) +   P_erfc_5(6) * exp_p * erfc_r  )) &
               / (  P_erfc_5(7) * rto(5) ) 


       case(6)

!         F_table(p) = &
!               ( ( exp (-(r**2 /Ewald_radius**2 )) * ((20790.0d0 * Ewald_radius**10 * r + 13860.0d0 * Ewald_radius**8  & 
!               * r**3 + 5544.0d0 * Ewald_radius**6 * r**5 + 1584.0d0 * Ewald_radius**4 * r**7 + &
!               352.0d0 *  Ewald_radius**2 * r**9 + 64.0d0 * r**11 + &  
!               10395.0d0 * exp(r**2 /Ewald_radius**2 ) * Ewald_radius**11 * &
!               sqrt_pi *  arch_erfc(r/Ewald_radius)) ) ) / (Ewald_radius**11  * sqrt_pi *  r**13 ) )

          F_table(p) = ( exp_m * &
               ( ddot(6, P_erfc_6(1:6),1, rto(0:5),1) + P_erfc_6(7) * exp_p * erfc_r )) &
               / (P_erfc_6(8) *  rto(6) ) 


       case(7)


          F_table(p) = &
               ( (- ( ( exp (- (r**2 /Ewald_radius**2 ) )  * ((270270.0d0 * Ewald_radius**12 *  r + &
               180180.0d0 * Ewald_radius**10 * r**3 + 72072.0d0 * Ewald_radius**8 * r**5 + 20592.0d0 * Ewald_radius**6 * r**7 + &
               4576.0d0 * Ewald_radius**4 * r**9 + 832.0d0 * Ewald_radius**2 * r**11 + 128.0d0 * r**13 + & 
               135135.0d0 *  exp (r**2 /Ewald_radius**2 ) * Ewald_radius**13 * &
               sqrt_pi * Arch_erfc(r /Ewald_radius)) ) ) / (Ewald_radius**13 * sqrt_pi * r**15 ) ) ) )

       case(8)

          F_table(p) = &
               ( ( exp (- (r**2 /Ewald_radius**2 ) ) * ((2 * ((2027025.0d0 * Ewald_radius**14 * r + &  
               1351350.0d0 * Ewald_radius**12 * r**3 + 540540.0d0 * Ewald_radius**10 * r**5 + &
               154440.0d0 * Ewald_radius**8 * r**7 + 34320.0d0 * Ewald_radius**6 *& 
               r**9 + 6240.0d0 * Ewald_radius**4 * r**11 + 960.0d0 * Ewald_radius**2 * r**13 + 128.0d0 * r**15) ) + 2027025.0d0 * &  
               exp (r**2 /Ewald_radius**2 ) * Ewald_radius**15 *  sqrt_pi * arch_erfc(r/Ewald_radius)) ) ) &
               / (Ewald_radius**15  * sqrt_pi * r**17 ) )

       case(9)

          F_table(p) = &
               ( (- ( ( exp(- (r**2 /Ewald_radius**2 ) ) * ((68918850.0d0 * Ewald_radius**16 * r + & 
               45945900.0d0 *Ewald_radius**14 * r**3 + 18378360.0d0 * Ewald_radius**12 * r**5 + &
               5250960.0d0 * Ewald_radius**10 * r**7 + 1166880.0d0 * Ewald_radius**8 * r**9 +  &
               212160.0d0 * Ewald_radius**6  *r**11 + 32640.0d0 * Ewald_radius**4 * r**13 + 4352.0d0 *  Ewald_radius**2*  r**15 + & 
               512.0d0 *  r**17 + &
               34459425.0d0 *exp (r**2 /Ewald_radius**2 ) * Ewald_radius**17  * &
               sqrt_pi * Arch_erfc(r /Ewald_radius)) ) ) / (Ewald_radius**17  * sqrt_pi*  r**19 ) ) ) )

       case(10)

          F_table(p) = &
               ( ( (1 / (Ewald_radius**19 * sqrt_pi * r**21 ) )* (( exp(- (r**2 /Ewald_radius**2 ) ) *&
               ((1309458150.0d0 * Ewald_radius**18 * r + 872972100.0d0 * Ewald_radius**16 * r**3 + & 
               349188840.0d0 * Ewald_radius**14 * r**5 + &
               99768240.0d0 *Ewald_radius**12 * r**7 + 22170720.0d0 * Ewald_radius**10 * r**9 + &
               4031040.0d0 * Ewald_radius**8 * r**11 + 620160.0d0 * Ewald_radius**6 * r**13 + & 
               82688.0d0 * Ewald_radius**4 * r**15 + &
               9728.0d0 * Ewald_radius**2 * r**17 + 1024.0d0 * r**19 + &
               654729075.0d0 *  exp (r**2 /Ewald_radius**2 ) * Ewald_radius**19 *  sqrt_pi * Arch_erfc(r /Ewald_radius)) )) ) ) )

       case(11)


          F_table(p) = &
               ( ( (1 / (Ewald_radius**21  * sqrt_pi * r**23 ) )* (( exp (- (r**2 /Ewald_radius**2 ) ) * &
               (( (-2.0d0)  * ((13749310575.0d0 * Ewald_radius**20 * r + 9166207050.0d0 * Ewald_radius**18*  r**3 + &
               3666482820.0d0 * Ewald_radius**16*  r**5 + 1047566520.0d0 * Ewald_radius**14*  r**7 + &
               232792560.0d0 * Ewald_radius**12 * r**9 + 42325920.0d0 * Ewald_radius**10*  r**11 + & 
               6511680.0d0 *  Ewald_radius**8  *r**13 + 868224.0d0 * Ewald_radius**6 * r**15 + &
               102144.0d0 * Ewald_radius**4 * r**17 + 10752.0d0 *  Ewald_radius**2 * r**19 + 1024.0d0 * r**21) ) - &
               13749310575.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**21  * sqrt_pi*  Arch_erfc(r /Ewald_radius)) )) ) ) )

       case(12)

          F_table(p) = &
               ( ( (1 / (Ewald_radius**23  * sqrt_pi * r**25 ) )* (( exp (- (r**2 /Ewald_radius**2 ) ) *&
               ((2 *  ((316234143225.0d0 * Ewald_radius**22 * r + 210822762150.0d0 * Ewald_radius**20 * r**3 + &
               84329104860.0d0 * Ewald_radius**18 * r**5 + 24094029960.0d0 * Ewald_radius**16 * r**7 + &
               5354228880.0d0 * Ewald_radius**14 * r**9 + 973496160.0d0 * Ewald_radius**12 * r**11 + &
               149768640.0d0 * Ewald_radius**10 * r**13 + 19969152.0d0 * Ewald_radius**8 * r**15 + &
               2349312.0d0 * Ewald_radius**6 * r**17 + 247296.0d0 * Ewald_radius**4 * r**19 + &
               23552.0d0 * Ewald_radius**2 * r**21 + 2048 * r**23) ) + &
               316234143225.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**23  * & 
               sqrt_pi * Arch_erfc(r /Ewald_radius)) )) ) ) )

       case(13)

          F_table(p) = &
               ( ( (1 / (Ewald_radius**25 *  sqrt_pi * r**27 ) ) *(( exp (- (r**2 /Ewald_radius**2 ) ) * &
               (( (-2 ) *  ((7905853580625.0d0 * Ewald_radius**24 * r + 5270569053750.0d0 * Ewald_radius**22*  r**3 + &
               2108227621500.0d0 * Ewald_radius**20 * r**5 + 602350749000.0d0 * Ewald_radius**18 * r**7 + &
               133855722000.0d0 * Ewald_radius**16 * r**9 + 24337404000.0d0 * Ewald_radius**14 * r**11 + &
               3744216000.0d0 * Ewald_radius**12 * r**13 + 499228800.0d0 * Ewald_radius**10 * r**15 + &
               58732800.0d0 * Ewald_radius**8 * r**17 + 6182400.0d0 * Ewald_radius**6*  r**19 + &
               588800.0d0 * Ewald_radius**4 * r**21 + 51200.0d0 * Ewald_radius**2 * r**23 + 4096.0d0 * r**25) ) - &
               7905853580625.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**25 * & 
               sqrt_pi *  arch_erfc(r/Ewald_radius)) )) ) ) )

       case(14)

          F_table(p) = &
               ( ( (1 / (Ewald_radius**27 *  sqrt_pi * r**29 ) )* (( exp (- (r**2 /Ewald_radius**2 ) ) * &
               ((2 *  ((213458046676875.0d0 * Ewald_radius**26 * r + 142305364451250.0d0 * Ewald_radius**24 * r**3 + &
               56922145780500.0d0 * Ewald_radius**22 * r**5 + 16263470223000.0d0 * Ewald_radius**20 * r**7 + &
               3614104494000.0d0 * Ewald_radius**18 * r**9 + 657109908000.0d0 * Ewald_radius**16 * r**11 + &
               101093832000.0d0 * Ewald_radius**14 * r**13 + 13479177600.0d0 * Ewald_radius**12 * r**15 + &
               1585785600.0d0 * Ewald_radius**10 * r**17 + 166924800.0d0 * Ewald_radius**8 * r**19 + &
               15897600.0d0 * Ewald_radius**6 * r**21 + 1382400.0d0 * Ewald_radius**4 * r**23 + &
               110592.0d0 * Ewald_radius**2 * r**25 + 8192 * r**27) ) + &
               213458046676875.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**27  * & 
               sqrt_pi* arch_erfc(r/Ewald_radius)) )) ) ) )


       end select


    end do

end subroutine F_erfc_table_original
  !******
  !------------------------------------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erf_single_original
  !  NAME
  !    F_erf_single_original
  !  SYNOPSIS

  subroutine F_erf_table_original(F_erf_table, r, p_max)

    !  PURPOSE
    !  The function calculates Fp erf function for the periodic systems.
    !  This is Fp function for regions inside the multipole radius.
    !  Close origin, we have to use polynomial expansion, because of the
    !  numberical problems.
    !

    use arch_specific, only: arch_erf
    use constants, only: sqrt_pi
    use runtime_choices, only : Ewald_radius
    implicit none
    !  ARGUMENTS

    integer:: p_max
    real*8:: F_erf_table(0:p_max), r


    !  INPUTS
    !    o r -- distance
    !    o p -- order of function
    !  OUTPUT
    !    o F_1_div_r -- value of Fp function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8,parameter:: n=100.0d0, epsi = 1e-8
    integer:: i,p


    F_erf_table = 0.0d0
!    return
    
    do p=0,p_max,1

       if( r   > (Ewald_radius) .or. p==0)then

          if ((r > 25*Ewald_radius).and.(p>0)) then
             F_erf_table(p) =  F_1_div_r_2( r, p )
          else


             select case(p)

       case(0)


          if(r < 1e-5)then
             do i=0,8
                F_erf_table(p) =  F_erf_table(p) + Kp(i) * r**(2*i)
             end do
          else
             F_erf_table(p) = (arch_erf(r/Ewald_radius)/r)
          end if


       case(1)

          F_erf_table(p) = ((2.0d0 * exp(-(r**2/Ewald_radius**2))* r) /  (Ewald_radius *sqrt_pi) - arch_erf(r/Ewald_radius))/r**3

       case(2)

          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) * (((-6.0d0) * Ewald_radius**2* r - 4.0d0 * r**3  &
               + 3.0d0 * exp(r**2/Ewald_radius**2) * Ewald_radius**3 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               / (Ewald_radius**3 * sqrt_pi * r**5))

       case(3)

          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) *((30.0d0 * Ewald_radius**4 *r & 
               + 20.0d0 * Ewald_radius**2 *r**3 + 8.0d0 *r**5 &
               - 15.0d0 * exp(r**2/Ewald_radius**2) *Ewald_radius**5 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               / (Ewald_radius**5 *sqrt_pi * r**7))

       case(4)

          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) *(((-2.0d0)* ((105.0d0 * Ewald_radius**6 *r + &
               70.0d0 * Ewald_radius**4* r**3 + 28.0d0 *Ewald_radius**2 *r**5 + 8.0d0 *r**7)) &
               + 105.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**7 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               /(Ewald_radius**7 *sqrt_pi *r**9+epsi))

       case(5)

          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2))* ((2.0d0 *((945.0d0 *Ewald_radius**8 *r & 
               + 630.0d0 *Ewald_radius**6 *r**3 &
               + 252.0d0 * Ewald_radius**4 *r**5 + 72.0d0 *Ewald_radius**2 *r**7 + 16.0d0 *r**9)) &
               - 945.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**9 * & 
               sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**9 *sqrt_pi *r**11))

       case(6)

          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) *(((-2.0d0) *((10395.0d0 *Ewald_radius**10 *r & 
               + 6930.0d0 *Ewald_radius**8 *r**3 &
               + 2772.0d0 *Ewald_radius**6 *r**5 & 
               + 792.0d0 *Ewald_radius**4 *r**7 + 176.0d0 *Ewald_radius**2 *r**9 + 32.0d0 *r**11)) &
               + 10395.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**11 * &
               sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**11 *sqrt_pi *r**13))

       case(7)

          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) *((270270.0d0 *Ewald_radius**12 *r & 
               + 180180.0d0 *Ewald_radius**10 *r**3 &
               + 72072.0d0 *Ewald_radius**8 *r**5 + 20592.0d0 *Ewald_radius**6 *r**7 &
               + 4576.0d0 *Ewald_radius**4 *r**9 + 832.0d0 *Ewald_radius**2 *r**11 + & 
               128.0d0 *r**13 - 135135.0d0 *exp(r**2/Ewald_radius**2) &
               *Ewald_radius**13 *sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**13 *sqrt_pi *r**15))

       case(8)
          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) &
               * (((-2.0d0) *((2027025.0d0 *Ewald_radius**14 *r + 1351350.0d0 *Ewald_radius**12 *r**3 &
               + 540540.0d0 *Ewald_radius**10 *r**5 + 154440.0d0 *Ewald_radius**8 *r**7 + 34320.0d0 *Ewald_radius**6 &
               *r**9 + 6240.0d0 *Ewald_radius**4 *r**11 + 960.0d0 *Ewald_radius**2 *r**13 + 128.0d0 *r**15)) + &
               2027025.0d0 *exp(r**2/Ewald_radius**2) &
               *Ewald_radius**15 *sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**15 *sqrt_pi *r**17))

       case(9)
          F_erf_table(p) = ((exp(-(r**2/Ewald_radius**2)) *((68918850.0d0 *Ewald_radius**16 *r & 
               + 45945900.0d0 *Ewald_radius**14 *r**3 &
               + 18378360.0d0 *Ewald_radius**12 *r**5 + 5250960.0d0 *Ewald_radius**10 *r**7 &
               + 1166880.0d0 *Ewald_radius**8 *r**9 + 212160.0d0 *Ewald_radius**6 *r**11 + &
               32640.0d0 *Ewald_radius**4 *r**13 + 4352.0d0 *Ewald_radius**2 *r**15 + 512.0d0 *r**17 &
               - 34459425.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**17 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               /(Ewald_radius**17 *sqrt_pi *r**19))

       case(10)
          F_erf_table(p) = (((1.0d0/(Ewald_radius**19 *sqrt_pi *r**21))*((exp(-(r**2/Ewald_radius**2)) * &
               (((-2.0d0) *((654729075.0d0 *Ewald_radius**18 *r + 436486050.0d0 *Ewald_radius**16 *r**3 + &
               174594420.0d0 *Ewald_radius**14 *r**5 + 49884120.0d0 *Ewald_radius**12 *r**7 & 
               + 11085360.0d0 *Ewald_radius**10 *r**9 &
               + 2015520.0d0 *Ewald_radius**8 *r**11 + 310080.0d0 *Ewald_radius**6 *r**13 &
               + 41344.0d0 *Ewald_radius**4 *r**15 + 4864.0d0 *Ewald_radius**2 *r**17 + 512.0d0 *r**19)) + &
               654729075.0d0 *exp(r**2/Ewald_radius**2) &
               *Ewald_radius**19 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(11)
          F_erf_table(p) = (((1.0d0/(Ewald_radius**21 *sqrt_pi *r**23))*((exp(-(r**2/Ewald_radius**2)) * &
               ((2.0d0 *((13749310575.0d0 *Ewald_radius**20 *r + 9166207050.0d0 *Ewald_radius**18 *r**3 &
               + 3666482820.0d0 *Ewald_radius**16 *r**5 + 1047566520.0d0 *Ewald_radius**14 *r**7 + &
               232792560.0d0 *Ewald_radius**12 *r**9 + 42325920.0d0 *Ewald_radius**10 *r**11 &
               +  6511680.0d0 *Ewald_radius**8 *r**13 + 868224.0d0 *Ewald_radius**6 *r**15 + &
               102144.0d0 *Ewald_radius**4 *r**17 + 10752.0d0 *Ewald_radius**2 *r**19 + 1024.0d0 *r**21)) &
               - 13749310575.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**21 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(12)
          F_erf_table(p) = (((1/(Ewald_radius**23 *sqrt_pi *r**25))*((exp(-(r**2/Ewald_radius**2)) *(((-2.0d0) &
               *((316234143225.0d0 *Ewald_radius**22 *r + 210822762150.0d0 *Ewald_radius**20 *r**3 + &
               84329104860.0d0 *Ewald_radius**18 *r**5 + 24094029960.0d0 *Ewald_radius**16 *r**7 + &
               5354228880.0d0 *Ewald_radius**14 *r**9 + 973496160.0d0 *Ewald_radius**12 *r**11 + &
               149768640.0d0 *Ewald_radius**10 *r**13 + 19969152.0d0 *Ewald_radius**8 *r**15 + & 
               2349312.0d0 *Ewald_radius**6 *r**17 + &
               247296.0d0 *Ewald_radius**4 *r**19 + &
               23552.0d0 * Ewald_radius**2 *r**21 + 2048.0d0 *r**23)) + &
               316234143225.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**23 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(13)
          F_erf_table(p) =(((1/(Ewald_radius**25 *sqrt_pi *r**27))*((exp(-(r**2/Ewald_radius**2)) &
               *((2.0d0 *((7905853580625.0d0 *Ewald_radius**24 *r + 5270569053750.0d0 *Ewald_radius**22 *r**3 &
               + 2108227621500.0d0 *Ewald_radius**20 *r**5 + 602350749000.0d0 *Ewald_radius**18 *r**7 + &
               133855722000.0d0 *Ewald_radius**16 *r**9 + 24337404000.0d0 *Ewald_radius**14 *r**11 &
               + 3744216000.0d0 *Ewald_radius**12 *r**13 + 499228800.0d0 *Ewald_radius**10 *r**15 & 
               + 58732800.0d0 *Ewald_radius**8 *r**17 &
               + 6182400.0d0 *Ewald_radius**6 *r**19 + 588800.0d0 *Ewald_radius**4 *r**21 &
               + 51200.0d0 *Ewald_radius**2 *r**23 + 4096.0d0 *r**25)) - &
               7905853580625.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**25 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(14)
          F_erf_table(p) = (((1.0d0/(Ewald_radius**27 *sqrt_pi* r**29))*((exp(-(r**2/Ewald_radius**2))* & 
               (((-2.0d0) *((213458046676875.0d0 *Ewald_radius**26 *r &
               + 142305364451250.0d0 *Ewald_radius**24 *r**3 + &
               56922145780500.0d0 *Ewald_radius**22 *r**5 + 16263470223000.0d0 *Ewald_radius**20 *r**7 + &
               3614104494000.0d0 *Ewald_radius**18 *r**9 + 657109908000.0d0 *Ewald_radius**16 *r**11 + &
               101093832000.0d0 *Ewald_radius**14 *r**13 + 13479177600.0d0 *Ewald_radius**12 *r**15 + &
               1585785600.0d0 *Ewald_radius**10 *r**17 + 166924800.0d0 *Ewald_radius**8 *r**19 + &
               15897600.0d0 *Ewald_radius**6 *r**21 + 1382400.0d0 *Ewald_radius**4 *r**23 + &
               110592.0d0 *Ewald_radius**2 *r**25 + 8192.0d0 *r**27)) + &
               213458046676875.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**27 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(15)
          F_erf_table(p) = (((1.0d0/(Ewald_radius**29 *sqrt_pi *r**31))*((exp(-(r**2/Ewald_radius**2))* &
               ((12380566707258750.0d0 *Ewald_radius**28 *r &
               + 8253711138172500.0d0 *Ewald_radius**26 *r**3 + &
               3301484455269000.0d0 *Ewald_radius**24 *r**5 + 943281272934000.0d0 *Ewald_radius**22 *r**7 + 209618060652000.0d0 &
               *Ewald_radius**20 *r**9 + 38112374664000.0d0 *Ewald_radius**18 *r**11 + &
               5863442256000.0d0 *Ewald_radius**16 *r**13 + 781792300800.0d0 *Ewald_radius**14 *r**15 + &
               91975564800.0d0 *Ewald_radius**12 *r**17 + 9681638400.0d0 *Ewald_radius**10 *r**19 + &
               922060800.0d0 *Ewald_radius**8 *r**21 + 80179200.0d0 *Ewald_radius**6 *r**23 + &
               6414336.0d0 *Ewald_radius**4 *r**25 + 475136.0d0 *Ewald_radius**2 *r**27 + 32768.0d0 *r**29 - &
               6190283353629375.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**29 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(16)
          F_erf_table(p) = (((1.0d0/(Ewald_radius**31 *sqrt_pi *r**33))*((exp(-(r**2/Ewald_radius**2))* &
               (((-2.0d0) *((191898783962510625.0d0 *Ewald_radius**30 *r + &
               127932522641673750.0d0 *Ewald_radius**28 *r**3 + 51173009056669500.0d0 *Ewald_radius**26 *r**5 + &
               14620859730477000.0d0 *Ewald_radius**24 *r**7 + &
               3249079940106000.0d0 *Ewald_radius**22 *r**9 + 590741807292000.0d0 *Ewald_radius**20 *r**11 + &
               90883354968000.0d0 *Ewald_radius**18 *r**13 &
               + 12117780662400.0d0 *Ewald_radius**16 *r**15 + &
               1425621254400.0d0 *Ewald_radius**14 *r**17 + 150065395200.0d0 *Ewald_radius**12 *r**19 + &
               14291942400.0d0 *Ewald_radius**10 *r**21 &
               + 1242777600.0d0 *Ewald_radius**8 *r**23 + &
               99422208.0d0 *Ewald_radius**6 *r**25 + 7364608.0d0 *Ewald_radius**4 *r**27 + &
               507904.0d0 *Ewald_radius**2 *r**29 + 32768.0d0 *r**31)) + &
               191898783962510625.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**31 *sqrt_pi  *arch_erf(r/Ewald_radius)))))))

       end select
    end if
    else


       F_erf_table(p) = 0.0d0
!       return
       select case(p)

       case(0)

          do i=0,8             
             F_erf_table(p) =  F_erf_table(p) + Kp(i) * r**(2*i)
          end do

       case(1)

          do i=0,8             
             F_erf_table(p) =   F_erf_table(p) + 2*(i+1) * Kp(i+1) * r**(2*i)
          end do

       case(2)

          do i=0,8             
             F_erf_table(p) =  F_erf_table(p) +  (2*i+2)* (2*i+4)* Kp(i+2) * r**(2*i)
          end do

       case(3)

          do i=0,8             
             F_erf_table(p)=  F_erf_table(p) +  (2*i+2)* (2*i+4)*(2*i+6)* Kp(i+3) * r**(2*i)
          end do

       case(4)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)* Kp(i+4) * r**(2*i)
          end do

       case(5)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*Kp(i+5) * r**(2*i)
          end do

       case(6)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*Kp(i+6) * r**(2*i)
          end do

       case(7)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*Kp(i+7) * r**(2*i)
          end do

       case(8)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) +  & 
             (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*(2*i+16)* Kp(i+8) * r**(2*i)
          end do

       case(9)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) +  & 
             (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*(2*i+16)*(2*i+18)* Kp(i+9) * r**(2*i)
          end do

       case(10)

          do i=0,8
             F_erf_table(p) =  F_erf_table(p) + & 
             (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*(2*i+16)*(2*i+18)*(2*i+18) &
                  *Kp(i+10) * r**(2*i)
          end do

       end select

    end if
 end do  ! p

end subroutine F_erf_table_original
  !******
  !-----------------------------------------------------------------------------
 !---------------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erfc_single_original
  !  NAME
  !   F_erfc_single_original
  !  SYNOPSIS

  function F_erfc_single_original( r,p)

    !  PURPOSE
    !  The function calculates Fp erfc function for the periodic systems.
    !  This is Fp function for regions outside the multipole radius.
    !
    use arch_specific, only: arch_erfc
    use constants, only: sqrt_pi
    use runtime_choices, only : Ewald_radius
    implicit none
    !  ARGUMENTS
    
    real*8:: F_erfc_single_original, r
    integer :: p

    !  INPUTS
    !    o r -- distance
    !    o p -- order of function
    !  OUTPUT
    !    o F_1_div_r -- value of Fp function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    !    real*8:: n=100

    if(r > 25*Ewald_radius)then

       F_erfc_single_original  = 0.0

    else

       select case(p)

       case(0)

          F_erfc_single_original  = arch_erfc(r/Ewald_radius)/r

       case(1)

          F_erfc_single_original =  - (((2.0d0 * exp(- (r**2 /Ewald_radius**2 )) * r ) & 
               / (Ewald_radius * sqrt_pi ) + Arch_erfc(r /Ewald_radius) ) /r**3 )

       case(2)



          F_erfc_single_original = &
               ( ( exp(-(r**2 /Ewald_radius**2 ) ) * ((6 * Ewald_radius**2 * r + 4 * r**3 + 3 *  exp(r**2 /Ewald_radius**2 ) * &
               Ewald_radius**3 * sqrt_pi * arch_erfc(r/Ewald_radius)) ) ) / (Ewald_radius**3 * sqrt_pi* r**5 ) ) 


       case(3)

          F_erfc_single_original = &
               -(( exp (- (r**2 /Ewald_radius**2 ) )  * ((30 * Ewald_radius**4 * r +  20 * Ewald_radius**2 * r**3 + 8 * r**5 +  &
               15  * exp(r**2 /Ewald_radius**2 ) * Ewald_radius**5  * sqrt_pi * arch_erfc( r /Ewald_radius)) ) ) &
               / (Ewald_radius**5 *  sqrt_pi * r**7 ))


       case(4)


          F_erfc_single_original = &
               ( ( exp(- (r**2 /Ewald_radius**2 ) ) * ((210.0d0 * Ewald_radius**6 * r + 140.0d0 * Ewald_radius**4 * r**3 &  
               + 56.0d0 * Ewald_radius**2 * r**5 + 16.0d0 * r**7 + 105.0d0 * exp(r**2 /Ewald_radius**2 ) * Ewald_radius**7 &
               *sqrt_pi* arch_erfc(r/Ewald_radius)) ) ) / (Ewald_radius**7  * sqrt_pi * r**9 ) )


       case(5)

          F_erfc_single_original = &
               ( (- ( ( exp(- (r**2 /Ewald_radius**2 )) *  ((2.0d0 *  ((945.0d0 * Ewald_radius**8 * r + &
               630.0d0 * Ewald_radius**6 * r**3 + 252.0d0 * Ewald_radius**4 * r**5 + 72.0d0 * Ewald_radius**2 * r**7 + & 
               16.0d0 * r**9) ) +  945.0d0 * exp(r**2 /Ewald_radius**2) * &
               Ewald_radius**9 * sqrt_pi * arch_erfc( r /Ewald_radius)))) / (Ewald_radius**9 * sqrt_pi * r**11 ) ) ) )

       case(6)

          F_erfc_single_original = &
               ( ( exp (-(r**2 /Ewald_radius**2 )) * ((20790.0d0 * Ewald_radius**10 * r + 13860.0d0 * Ewald_radius**8  & 
               * r**3 + 5544.0d0 * Ewald_radius**6 * r**5 + 1584.0d0 * Ewald_radius**4 * r**7 + &
               352.0d0 *  Ewald_radius**2 * r**9 + 64.0d0 * r**11 + &  
               10395.0d0 * exp(r**2 /Ewald_radius**2 ) * Ewald_radius**11 * &
               sqrt_pi *  arch_erfc(r/Ewald_radius)) ) ) / (Ewald_radius**11  * sqrt_pi *  r**13 ) )

       case(7)


          F_erfc_single_original = &
               ( (- ( ( exp (- (r**2 /Ewald_radius**2 ) )  * ((270270.0d0 * Ewald_radius**12 *  r + &
               180180.0d0 * Ewald_radius**10 * r**3 + 72072.0d0 * Ewald_radius**8 * r**5 + 20592.0d0 * Ewald_radius**6 * r**7 + &
               4576.0d0 * Ewald_radius**4 * r**9 + 832.0d0 * Ewald_radius**2 * r**11 + 128.0d0 * r**13 + & 
               135135.0d0 *  exp (r**2 /Ewald_radius**2 ) * Ewald_radius**13 * &
               sqrt_pi * Arch_erfc(r /Ewald_radius)) ) ) / (Ewald_radius**13 * sqrt_pi * r**15 ) ) ) )

       case(8)

          F_erfc_single_original = &
               ( ( exp (- (r**2 /Ewald_radius**2 ) ) * ((2 * ((2027025.0d0 * Ewald_radius**14 * r + &  
               1351350.0d0 * Ewald_radius**12 * r**3 + 540540.0d0 * Ewald_radius**10 * r**5 + &
               154440.0d0 * Ewald_radius**8 * r**7 + 34320.0d0 * Ewald_radius**6 *& 
               r**9 + 6240.0d0 * Ewald_radius**4 * r**11 + 960.0d0 * Ewald_radius**2 * r**13 + 128.0d0 * r**15) ) & 
               + 2027025.0d0 * &  
               exp (r**2 /Ewald_radius**2 ) * Ewald_radius**15 *  sqrt_pi * arch_erfc(r/Ewald_radius)) ) ) &
               / (Ewald_radius**15  * sqrt_pi * r**17 ) )

       case(9)

          F_erfc_single_original = &
               ( (- ( ( exp(- (r**2 /Ewald_radius**2 ) ) * ((68918850.0d0 * Ewald_radius**16 * r + & 
               45945900.0d0 *Ewald_radius**14 * r**3 + 18378360.0d0 * Ewald_radius**12 * r**5 + &
               5250960.0d0 * Ewald_radius**10 * r**7 + 1166880.0d0 * Ewald_radius**8 * r**9 +  &
               212160.0d0 * Ewald_radius**6  *r**11 + 32640.0d0 * Ewald_radius**4 * r**13 + & 
               4352.0d0 *  Ewald_radius**2*  r**15 + 512.0d0 *  r**17 + &
               34459425.0d0 *exp (r**2 /Ewald_radius**2 ) * Ewald_radius**17  * &
               sqrt_pi * Arch_erfc(r /Ewald_radius)) ) ) / (Ewald_radius**17  * sqrt_pi*  r**19 ) ) ) )

       case(10)

          F_erfc_single_original = &
               ( ( (1 / (Ewald_radius**19 * sqrt_pi * r**21 ) )* (( exp(- (r**2 /Ewald_radius**2 ) ) *&
               ((1309458150.0d0 * Ewald_radius**18 * r + 872972100.0d0 * Ewald_radius**16 * r**3 + & 
               349188840.0d0 * Ewald_radius**14 * r**5 + &
               99768240.0d0 *Ewald_radius**12 * r**7 + 22170720.0d0 * Ewald_radius**10 * r**9 + &
               4031040.0d0 * Ewald_radius**8 * r**11 + 620160.0d0 * Ewald_radius**6 * r**13 + & 
               82688.0d0 * Ewald_radius**4 * r**15 + 9728.0d0 * Ewald_radius**2 * r**17 + 1024.0d0 * r**19 + &
               654729075.0d0 *  exp (r**2 /Ewald_radius**2 ) * Ewald_radius**19 *  sqrt_pi * Arch_erfc(r /Ewald_radius)) )) ) ) )

       case(11)


          F_erfc_single_original = &
               ( ( (1 / (Ewald_radius**21  * sqrt_pi * r**23 ) )* (( exp (- (r**2 /Ewald_radius**2 ) ) * &
               (( (-2.0d0)  * ((13749310575.0d0 * Ewald_radius**20 * r + 9166207050.0d0 * Ewald_radius**18*  r**3 + &
               3666482820.0d0 * Ewald_radius**16*  r**5 + 1047566520.0d0 * Ewald_radius**14*  r**7 + &
               232792560.0d0 * Ewald_radius**12 * r**9 + 42325920.0d0 * Ewald_radius**10*  r**11 + & 
               6511680.0d0 *  Ewald_radius**8  *r**13 + 868224.0d0 * Ewald_radius**6 * r**15 + &
               102144.0d0 * Ewald_radius**4 * r**17 + 10752.0d0 *  Ewald_radius**2 * r**19 + 1024.0d0 * r**21) ) - &
               13749310575.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**21  * sqrt_pi*  Arch_erfc(r /Ewald_radius)) )) ) ) )

       case(12)

          F_erfc_single_original = &
               ( ( (1 / (Ewald_radius**23  * sqrt_pi * r**25 ) )* (( exp (- (r**2 /Ewald_radius**2 ) ) *&
               ((2 *  ((316234143225.0d0 * Ewald_radius**22 * r + 210822762150.0d0 * Ewald_radius**20 * r**3 + &
               84329104860.0d0 * Ewald_radius**18 * r**5 + 24094029960.0d0 * Ewald_radius**16 * r**7 + &
               5354228880.0d0 * Ewald_radius**14 * r**9 + 973496160.0d0 * Ewald_radius**12 * r**11 + &
               149768640.0d0 * Ewald_radius**10 * r**13 + 19969152.0d0 * Ewald_radius**8 * r**15 + &
               2349312.0d0 * Ewald_radius**6 * r**17 + 247296.0d0 * Ewald_radius**4 * r**19 + &
               23552.0d0 * Ewald_radius**2 * r**21 + 2048 * r**23) ) + &
               316234143225.0d0 * exp (r**2 /Ewald_radius**2 ) * & 
               Ewald_radius**23  * sqrt_pi * Arch_erfc(r /Ewald_radius)) )) ) ) )

       case(13)

          F_erfc_single_original = &
               ( ( (1 / (Ewald_radius**25 *  sqrt_pi * r**27 ) ) *(( exp (- (r**2 /Ewald_radius**2 ) ) * &
               (( (-2 ) *  ((7905853580625.0d0 * Ewald_radius**24 * r + 5270569053750.0d0 * Ewald_radius**22*  r**3 + &
               2108227621500.0d0 * Ewald_radius**20 * r**5 + 602350749000.0d0 * Ewald_radius**18 * r**7 + &
               133855722000.0d0 * Ewald_radius**16 * r**9 + 24337404000.0d0 * Ewald_radius**14 * r**11 + &
               3744216000.0d0 * Ewald_radius**12 * r**13 + 499228800.0d0 * Ewald_radius**10 * r**15 + &
               58732800.0d0 * Ewald_radius**8 * r**17 + 6182400.0d0 * Ewald_radius**6*  r**19 + &
               588800.0d0 * Ewald_radius**4 * r**21 + 51200.0d0 * Ewald_radius**2 * r**23 + 4096.0d0 * r**25) ) - &
               7905853580625.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**25 *  & 
               sqrt_pi *  arch_erfc(r/Ewald_radius)) )) ) ) )

       case(14)

          F_erfc_single_original = &
               ( ( (1 / (Ewald_radius**27 *  sqrt_pi * r**29 ) )* (( exp (- (r**2 /Ewald_radius**2 ) ) * &
               ((2 *  ((213458046676875.0d0 * Ewald_radius**26 * r + 142305364451250.0d0 * Ewald_radius**24 * r**3 + &
               56922145780500.0d0 * Ewald_radius**22 * r**5 + 16263470223000.0d0 * Ewald_radius**20 * r**7 + &
               3614104494000.0d0 * Ewald_radius**18 * r**9 + 657109908000.0d0 * Ewald_radius**16 * r**11 + &
               101093832000.0d0 * Ewald_radius**14 * r**13 + 13479177600.0d0 * Ewald_radius**12 * r**15 + &
               1585785600.0d0 * Ewald_radius**10 * r**17 + 166924800.0d0 * Ewald_radius**8 * r**19 + &
               15897600.0d0 * Ewald_radius**6 * r**21 + 1382400.0d0 * Ewald_radius**4 * r**23 + &
               110592.0d0 * Ewald_radius**2 * r**25 + 8192 * r**27) ) + &
               213458046676875.0d0 * exp (r**2 /Ewald_radius**2 ) * Ewald_radius**27  * & 
               sqrt_pi* arch_erfc(r/Ewald_radius)) )) ) ) )


       end select

    end if

  end function F_erfc_single_original
  !******
  !------------------------------------------------------------------------------------------
  !****f* Hartree_F_p_functions/F_erf_single_original
  !  NAME
  !    F_erf_single_original
  !  SYNOPSIS

  function F_erf_single_original(r, p)

    !  PURPOSE
    !  The function calculates Fp erf function for the periodic systems.
    !  This is Fp function for regions inside the multipole radius.
    !  Close origin, we have to use polynomial expansion, because of the
    !  numberical problems.
    !

    use arch_specific, only: arch_erf
    use constants, only: sqrt_pi
    use runtime_choices, only : Ewald_radius
    implicit none
    !  ARGUMENTS

    real*8:: F_erf_single_original, r
    integer:: p

    !  INPUTS
    !    o r -- distance
    !    o p -- order of function
    !  OUTPUT
    !    o F_1_div_r -- value of Fp function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8,parameter:: n=100.0d0, epsi = 1e-8
    integer:: i

    F_erf_single_original = 0.0d0
    
    if ((r > 25*Ewald_radius).and.(p>0)) then
       F_erf_single_original =  F_1_div_r_2( r, p )

    else if( r   > (Ewald_radius) .or. p==0)then

       select case(p)

       case(0)


          if(r < 1e-5)then
             do i=0,8
                F_erf_single_original =  F_erf_single_original + Kp(i) * r**(2*i)
             end do
          else
             F_erf_single_original = (arch_erf(r/Ewald_radius)/r)
          end if


       case(1)

          F_erf_single_original = ((2.0d0 * exp(-(r**2/Ewald_radius**2))* r) /  (Ewald_radius *sqrt_pi) & 
                                   - arch_erf(r/Ewald_radius))/r**3

       case(2)

          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) * (((-6.0d0) * Ewald_radius**2* r - 4.0d0 * r**3  &
               + 3.0d0 * exp(r**2/Ewald_radius**2) * Ewald_radius**3 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               / (Ewald_radius**3 * sqrt_pi * r**5))

       case(3)

          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) *((30.0d0 * Ewald_radius**4 *r + 20.0d0 * Ewald_radius**2 *r**3 & 
               + 8.0d0 *r**5 - 15.0d0 * exp(r**2/Ewald_radius**2) *Ewald_radius**5 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               / (Ewald_radius**5 *sqrt_pi * r**7))

       case(4)

          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) *(((-2.0d0)* ((105.0d0 * Ewald_radius**6 *r + &
               70.0d0 * Ewald_radius**4* r**3 + 28.0d0 *Ewald_radius**2 *r**5 + 8.0d0 *r**7)) &
               + 105.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**7 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               /(Ewald_radius**7 *sqrt_pi *r**9+epsi))

       case(5)

          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2))* ((2.0d0 *((945.0d0 *Ewald_radius**8 *r & 
               + 630.0d0 *Ewald_radius**6 *r**3 &
               + 252.0d0 * Ewald_radius**4 *r**5 + 72.0d0 *Ewald_radius**2 *r**7 + 16.0d0 *r**9)) &
               - 945.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**9 * & 
               sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**9 *sqrt_pi *r**11))

       case(6)

          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) *(((-2.0d0) *((10395.0d0 *Ewald_radius**10 *r & 
               + 6930.0d0 *Ewald_radius**8 *r**3 + 2772.0d0 *Ewald_radius**6 *r**5 & 
               + 792.0d0 *Ewald_radius**4 *r**7 + 176.0d0 *Ewald_radius**2 *r**9 + 32.0d0 *r**11)) &
               + 10395.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**11 * &
               sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**11 *sqrt_pi *r**13))

       case(7)

          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) *((270270.0d0 *Ewald_radius**12 *r & 
               + 180180.0d0 *Ewald_radius**10 *r**3 &
               + 72072.0d0 *Ewald_radius**8 *r**5 + 20592.0d0 *Ewald_radius**6 *r**7 &
               + 4576.0d0 *Ewald_radius**4 *r**9 + 832.0d0 *Ewald_radius**2 *r**11 + & 
               128.0d0 *r**13 - 135135.0d0 *exp(r**2/Ewald_radius**2) &
               *Ewald_radius**13 *sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**13 *sqrt_pi *r**15))

       case(8)
          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) &
               * (((-2.0d0) *((2027025.0d0 *Ewald_radius**14 *r + 1351350.0d0 *Ewald_radius**12 *r**3 &
               + 540540.0d0 *Ewald_radius**10 *r**5 + 154440.0d0 *Ewald_radius**8 *r**7 + 34320.0d0 *Ewald_radius**6 &
               *r**9 + 6240.0d0 *Ewald_radius**4 *r**11 + 960.0d0 *Ewald_radius**2 *r**13 + 128.0d0 *r**15)) + &
               2027025.0d0 *exp(r**2/Ewald_radius**2) &
               *Ewald_radius**15 *sqrt_pi *arch_erf(r/Ewald_radius))))/(Ewald_radius**15 *sqrt_pi *r**17))

       case(9)
          F_erf_single_original = ((exp(-(r**2/Ewald_radius**2)) *((68918850.0d0 *Ewald_radius**16 *r & 
               + 45945900.0d0 *Ewald_radius**14 *r**3 &
               + 18378360.0d0 *Ewald_radius**12 *r**5 + 5250960.0d0 *Ewald_radius**10 *r**7 &
               + 1166880.0d0 *Ewald_radius**8 *r**9 + 212160.0d0 *Ewald_radius**6 *r**11 + &
               32640.0d0 *Ewald_radius**4 *r**13 + 4352.0d0 *Ewald_radius**2 *r**15 + 512.0d0 *r**17 &
               - 34459425.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**17 *sqrt_pi *arch_erf(r/Ewald_radius)))) &
               /(Ewald_radius**17 *sqrt_pi *r**19))

       case(10)
          F_erf_single_original = (((1.0d0/(Ewald_radius**19 *sqrt_pi *r**21))*((exp(-(r**2/Ewald_radius**2)) * &
               (((-2.0d0) *((654729075.0d0 *Ewald_radius**18 *r + 436486050.0d0 *Ewald_radius**16 *r**3 + &
               174594420.0d0 *Ewald_radius**14 *r**5 + 49884120.0d0 *Ewald_radius**12 *r**7 + 11085360.0d0 *Ewald_radius**10 *r**9 &
               + 2015520.0d0 *Ewald_radius**8 *r**11 + 310080.0d0 *Ewald_radius**6 *r**13 &
               + 41344.0d0 *Ewald_radius**4 *r**15 + 4864.0d0 *Ewald_radius**2 *r**17 + 512.0d0 *r**19)) + &
               654729075.0d0 *exp(r**2/Ewald_radius**2) &
               *Ewald_radius**19 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(11)
          F_erf_single_original = (((1.0d0/(Ewald_radius**21 *sqrt_pi *r**23))*((exp(-(r**2/Ewald_radius**2)) * &
               ((2.0d0 *((13749310575.0d0 *Ewald_radius**20 *r + 9166207050.0d0 *Ewald_radius**18 *r**3 &
               + 3666482820.0d0 *Ewald_radius**16 *r**5 + 1047566520.0d0 *Ewald_radius**14 *r**7 + &
               232792560.0d0 *Ewald_radius**12 *r**9 + 42325920.0d0 *Ewald_radius**10 *r**11 &
               +  6511680.0d0 *Ewald_radius**8 *r**13 + 868224.0d0 *Ewald_radius**6 *r**15 + &
               102144.0d0 *Ewald_radius**4 *r**17 + 10752.0d0 *Ewald_radius**2 *r**19 + 1024.0d0 *r**21)) &
               - 13749310575.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**21 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(12)
          F_erf_single_original = (((1/(Ewald_radius**23 *sqrt_pi *r**25))*((exp(-(r**2/Ewald_radius**2)) *(((-2.0d0) &
               *((316234143225.0d0 *Ewald_radius**22 *r + 210822762150.0d0 *Ewald_radius**20 *r**3 + &
               84329104860.0d0 *Ewald_radius**18 *r**5 + 24094029960.0d0 *Ewald_radius**16 *r**7 + &
               5354228880.0d0 *Ewald_radius**14 *r**9 + 973496160.0d0 *Ewald_radius**12 *r**11 + &
               149768640.0d0 *Ewald_radius**10 *r**13 + 19969152.0d0 *Ewald_radius**8 *r**15 + & 
               2349312.0d0 *Ewald_radius**6 *r**17 + &
               247296.0d0 *Ewald_radius**4 *r**19 + &
               23552.0d0 * Ewald_radius**2 *r**21 + 2048.0d0 *r**23)) + &
               316234143225.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**23 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(13)
          F_erf_single_original =(((1/(Ewald_radius**25 *sqrt_pi *r**27))*((exp(-(r**2/Ewald_radius**2)) &
               *((2.0d0 *((7905853580625.0d0 *Ewald_radius**24 *r + 5270569053750.0d0 *Ewald_radius**22 *r**3 &
               + 2108227621500.0d0 *Ewald_radius**20 *r**5 + 602350749000.0d0 *Ewald_radius**18 *r**7 + &
               133855722000.0d0 *Ewald_radius**16 *r**9 + 24337404000.0d0 *Ewald_radius**14 *r**11 &
               + 3744216000.0d0 *Ewald_radius**12 *r**13 + 499228800.0d0 *Ewald_radius**10 *r**15 & 
               + 58732800.0d0 *Ewald_radius**8 *r**17 &
               + 6182400.0d0 *Ewald_radius**6 *r**19 + 588800.0d0 *Ewald_radius**4 *r**21 &
               + 51200.0d0 *Ewald_radius**2 *r**23 + 4096.0d0 *r**25)) - &
               7905853580625.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**25 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(14)
          F_erf_single_original = (((1.0d0/(Ewald_radius**27 *sqrt_pi* r**29))*((exp(-(r**2/Ewald_radius**2))* & 
               (((-2.0d0) *((213458046676875.0d0 *Ewald_radius**26 *r &
               + 142305364451250.0d0 *Ewald_radius**24 *r**3 + &
               56922145780500.0d0 *Ewald_radius**22 *r**5 + 16263470223000.0d0 *Ewald_radius**20 *r**7 + &
               3614104494000.0d0 *Ewald_radius**18 *r**9 + 657109908000.0d0 *Ewald_radius**16 *r**11 + &
               101093832000.0d0 *Ewald_radius**14 *r**13 + 13479177600.0d0 *Ewald_radius**12 *r**15 + &
               1585785600.0d0 *Ewald_radius**10 *r**17 + 166924800.0d0 *Ewald_radius**8 *r**19 + &
               15897600.0d0 *Ewald_radius**6 *r**21 + 1382400.0d0 *Ewald_radius**4 *r**23 + &
               110592.0d0 *Ewald_radius**2 *r**25 + 8192.0d0 *r**27)) + &
               213458046676875.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**27 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(15)
          F_erf_single_original = (((1.0d0/(Ewald_radius**29 *sqrt_pi *r**31))*((exp(-(r**2/Ewald_radius**2))* &
               ((12380566707258750.0d0 *Ewald_radius**28 *r &
               + 8253711138172500.0d0 *Ewald_radius**26 *r**3 + &
               3301484455269000.0d0 *Ewald_radius**24 *r**5 + 943281272934000.0d0 *Ewald_radius**22 *r**7 + 209618060652000.0d0 &
               *Ewald_radius**20 *r**9 + 38112374664000.0d0 *Ewald_radius**18 *r**11 + &
               5863442256000.0d0 *Ewald_radius**16 *r**13 + 781792300800.0d0 *Ewald_radius**14 *r**15 + &
               91975564800.0d0 *Ewald_radius**12 *r**17 + 9681638400.0d0 *Ewald_radius**10 *r**19 + &
               922060800.0d0 *Ewald_radius**8 *r**21 + 80179200.0d0 *Ewald_radius**6 *r**23 + &
               6414336.0d0 *Ewald_radius**4 *r**25 + 475136.0d0 *Ewald_radius**2 *r**27 + 32768.0d0 *r**29 - &
               6190283353629375.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**29 *sqrt_pi *arch_erf(r/Ewald_radius)))))))

       case(16)
          F_erf_single_original = (((1.0d0/(Ewald_radius**31 *sqrt_pi *r**33))*((exp(-(r**2/Ewald_radius**2))* &
               (((-2.0d0) *((191898783962510625.0d0 *Ewald_radius**30 *r + &
               127932522641673750.0d0 *Ewald_radius**28 *r**3 + 51173009056669500.0d0 *Ewald_radius**26 *r**5 + &
               14620859730477000.0d0 *Ewald_radius**24 *r**7 + &
               3249079940106000.0d0 *Ewald_radius**22 *r**9 + 590741807292000.0d0 *Ewald_radius**20 *r**11 + &
               90883354968000.0d0 *Ewald_radius**18 *r**13 &
               + 12117780662400.0d0 *Ewald_radius**16 *r**15 + &
               1425621254400.0d0 *Ewald_radius**14 *r**17 + 150065395200.0d0 *Ewald_radius**12 *r**19 + &
               14291942400.0d0 *Ewald_radius**10 *r**21 &
               + 1242777600.0d0 *Ewald_radius**8 *r**23 + &
               99422208.0d0 *Ewald_radius**6 *r**25 + 7364608.0d0 *Ewald_radius**4 *r**27 + &
               507904.0d0 *Ewald_radius**2 *r**29 + 32768.0d0 *r**31)) + &
               191898783962510625.0d0 *exp(r**2/Ewald_radius**2) *Ewald_radius**31 *sqrt_pi  *arch_erf(r/Ewald_radius)))))))

       end select

    else


       F_erf_single_original = 0.0d0

       select case(p)

       case(0)

          do i=0,8             
             F_erf_single_original =  F_erf_single_original + Kp(i) * r**(2*i)
          end do

       case(1)

          do i=0,8             
             F_erf_single_original =   F_erf_single_original + 2*(i+1) * Kp(i+1) * r**(2*i)
          end do

       case(2)

          do i=0,8             
             F_erf_single_original =  F_erf_single_original +  (2*i+2)* (2*i+4)* Kp(i+2) * r**(2*i)
          end do

       case(3)

          do i=0,8             
             F_erf_single_original=  F_erf_single_original +  (2*i+2)* (2*i+4)*(2*i+6)* Kp(i+3) * r**(2*i)
          end do

       case(4)

          do i=0,8
             F_erf_single_original =  F_erf_single_original +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)* Kp(i+4) * r**(2*i)
          end do

       case(5)

          do i=0,8
             F_erf_single_original =  F_erf_single_original +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*Kp(i+5) * r**(2*i)
          end do

       case(6)

          do i=0,8
             F_erf_single_original =  F_erf_single_original +  (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*Kp(i+6) * r**(2*i)
          end do

       case(7)

          do i=0,8
             F_erf_single_original =  F_erf_single_original +  & 
                (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*Kp(i+7) * r**(2*i)
          end do

       case(8)

          do i=0,8
             F_erf_single_original =  F_erf_single_original + & 
                (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*(2*i+16)* Kp(i+8) * r**(2*i)
          end do

       case(9)

          do i=0,8
             F_erf_single_original =  F_erf_single_original + & 
                (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*(2*i+16)*(2*i+18)* Kp(i+9) * r**(2*i)
          end do

       case(10)

          do i=0,8
             F_erf_single_original =  F_erf_single_original + & 
                 (2*i+2)* (2*i+4)*(2*i+6)*(2*i+8)*(2*i+10)*(2*i+12)*(2*i+14)*(2*i+16)*(2*i+18)*(2*i+18) &
                  *Kp(i+10) * r**(2*i)
          end do

       end select

    end if

  end function F_erf_single_original
  !******




end module Hartree_F_p_functions
