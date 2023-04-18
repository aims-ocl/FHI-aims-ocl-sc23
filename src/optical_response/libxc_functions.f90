!****** FHI-aims/initialise_libxc
!  NAME
!    initialise_libxc
!  SYNOPSIS
!    Initialises and destroys LibXC, to make sure it works OK.

     subroutine initialise_libxc() 

!  PURPOSE
!    Seems to be an outdated process, moved from read_control. 
!    Hopefully this can just be removed completely? Needs testing
!
!  USES

     use xc_f03_lib_m

     implicit none

! ARGUMENTS

     type(xc_f03_func_t)        :: xc_func
     type(xc_f03_func_info_t)   :: xc_info

!  INPUTS
!
!  OUTPUT
!
!  AUTHOR
!    AJL, UCL. 2016. Contact: a.logsdail@ucl.ac.uk
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
!

  call xc_f03_func_init(xc_func, 1, XC_UNPOLARIZED)
  call xc_f03_func_end(xc_func)

end subroutine initialise_libxc

!****** FHI-aims/check_libxc_func
!  NAME
!    check_libxc_func
!  SYNOPSIS
!    Goes through the process of making sure the choices in LibXC are reasonable

  subroutine check_libxc_func(selected_functional, tddft, functional_type, hybrid_percentage)

!  PURPOSE
!    Make sure the choices are reasonable and print out associated literature
!
!  USES

  use localorb_io
  use mpi_tasks, only: aims_stop_coll
  use runtime_choices
  use xc_f03_lib_m
  implicit none

! ARGUMENTS

  integer, intent(in)           :: selected_functional
  logical, intent(in)           :: tddft
  integer, intent(out)          :: functional_type
  real*8,  intent(out)          :: hybrid_percentage

!  INPUTS
!  o selected_functional: the ID of the functional of interest
!  o tddft: logical to check if we are doing TDDFT and then filter
!           out unsuitable functionals.
!
!  OUTPUT
!
!  AUTHOR
!    Originally Jan Kloppenberg, ~2013.
!    Updated by AJL, UCL. 2016. Contact: a.logsdail@ucl.ac.uk
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
!  LOCAL VARIABLES

  integer                       :: i, j
  type(xc_f03_func_t)           :: xc_func
  type(xc_f03_func_info_t)      :: xc_info 
  character*128                 :: s1_info, s2_info, s3_info, info_str

! SOURCE

  ! Initialise
  call xc_f03_func_init(xc_func, selected_functional, XC_UNPOLARIZED)
  xc_info = xc_f03_func_get_info(xc_func)

  ! Store this so we can compare outside this subroutine
  functional_type = xc_f03_func_info_get_kind(xc_info)

  ! Print out info on type
  select case(functional_type)
  case(XC_EXCHANGE)
    write(s3_info, '(a)') 'Exchange'
    if (tddft) libxc_tddft_chosen_x = .true.
  case(XC_CORRELATION)
    write(s3_info, '(a)') 'Correlation'
    if (tddft) libxc_tddft_chosen_c = .true.
  case(XC_EXCHANGE_CORRELATION)
    write(s3_info, '(a)') 'Exchange-correlation'
    if (tddft) then
      libxc_tddft_chosen_x = .true.
      libxc_tddft_chosen_c = .true.
    endif
  end select

  ! Print out family
  s1_info = xc_f03_func_info_get_name(xc_info)
  select case(xc_f03_func_info_get_family(xc_info))
    case (XC_FAMILY_LDA);      write(s2_info,'(a)') "LDA"
    case (XC_FAMILY_GGA);      write(s2_info,'(a)') "GGA"
    case (XC_FAMILY_HYB_GGA);  write(s2_info,'(a)') "Hybrid GGA"
    case (XC_FAMILY_MGGA);     write(s2_info,'(a)') "MGGA"
    case (XC_FAMILY_HYB_MGGA); write(s2_info,'(a)') "Hybrid MGGA"
    case (XC_FAMILY_LCA);      write(s2_info,'(a)') "LCA"
  end select

  ! Write out references
  write(info_str, '(2X,2a)') "Description: ", trim(s1_info)!, ' (', trim(s3_info), '-', trim(s2_info), ')'
  call localorb_info(info_str)
  write(info_str, '(2X,5a)') "Type: ", trim(s3_info), " (", trim(s2_info), ")"
  call localorb_info(info_str)
  write(info_str, '(2X,1a)') "Reference to Literature:"
  call localorb_info(info_str)

  ! Need two indices - if someone can tidy this to one then I'd be grateful!
  i = 0
  j = 0
  do while(i.gt.-1)
    j = j + 1
    ! AJL, Nov 2017. It'd be nice to make this flexible.... though we'd need some stubs
    !s1_info = xc_f03_func_info_get_refs(xc_info, i) ! LibXC v3.0.0
    s1_info = xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(xc_info, i)) ! LibXC 4.0.2
    write(info_str, '(2X,a,i1,2a)') '[', j, '] ', trim(s1_info)
    call localorb_info(info_str)
  enddo

  ! Check for, and get if neccesary, the Hybrid Exact Exchange
  if ( .not.tddft.and. &
  ( xc_f03_func_info_get_family(xc_info).eq.XC_FAMILY_HYB_GGA ) .or. &
  ( xc_f03_func_info_get_family(xc_info).eq.XC_FAMILY_HYB_MGGA )) then
    ! call xc_f03_hyb_exx_coef(xc_info, hybrid_percentage)
    hybrid_percentage = xc_f03_hyb_exx_coef(xc_func);
    ! write(use_unit,*) 'HYBRID COEF: ', hybrid_percentage
  endif

  ! Catch error as TDDFT only works with LDA
  if ( tddft.and.(xc_f03_func_info_get_family(xc_info).ne.XC_FAMILY_LDA) ) then
    write(use_unit,'(4x,a)') '*** At the moment, TDDFT calculations only support LDA functionals.'
    call aims_stop_coll('Please select only LDA functionals from libxc.', 'check_libxc_func')
  endif

  ! Terminate
  call xc_f03_func_end(xc_func)
end subroutine check_libxc_func

subroutine check_tddft_fxc_func(selected_functional)
  use localorb_io
  use mpi_tasks, only: aims_stop
!  use xc_f03_lib_m
  implicit none

  integer, intent(in)           :: selected_functional

  integer                       :: k_cnt
  character*128			:: info_str
!  type(xc_f03_func_t)           :: xc_func

  integer, parameter            :: non_working_fxc(44) = (/ 16, 18, 22, 23, 50, 51, 131, 132, 134, &
                                       145, 146, 160, 165, 166, 174, 175, 182, 401, 402, 403, 404, &
                                       405, 406, 409, 411, 412, 415, 416, 417, 418, 419, 201, 202, &
                                       203, 204, 205, 206, 207, 208, 209, 210, 211, 301, 302 /)

!  Not actually needed here. AJL
!  call xc_f03_func_init(xc_func, selected_functional, XC_UNPOLARIZED)

  do k_cnt=1, size(non_working_fxc)
   if(selected_functional==non_working_fxc(k_cnt)) then
    write(info_str,'(2X,a)') "You have selected a functional for TDDFT that does not provide an f_xc."
    call localorb_info(info_str)
    write(info_str,'(2X,a)') "Please choose one that does and restart."
    call localorb_info(info_str)
    call aims_stop('Invalid f_xc functional chosen!')
   endif
  enddo

!  call xc_f03_func_end(xc_func)
end subroutine check_tddft_fxc_func

subroutine check_funcs_consistency()
  use localorb_io
  use mpi_tasks, only: aims_stop
  use runtime_choices
  implicit none

  character*128	:: info_str

  if(.not.libxc_tddft_chosen_x) then
    write(info_str,'(2x,a)') 'You are missing the exchange part from libxc for TDDFT.'
    call localorb_info(info_str)
    call aims_stop('Cannot continue without exchange for TDDFT')
  endif
  if(.not.libxc_tddft_chosen_c) then
    write(info_str,'(2x,a)') 'You are missing the correlation part from libxc for TDDFT.'
    call localorb_info(info_str)
    call aims_stop('Cannot continue without correlation for TDDFT')
  endif
end subroutine check_funcs_consistency

! Should now be redundant
!subroutine determine_func_number(read_str, number_id)
!  use xc_f03_lib_m
!  implicit none
!  
!  character*40, intent(in)	:: read_str
!  integer, intent(out)		:: number_id
!
!  select case(read_str)
!  
!  case('XC_LDA_X')
!    number_id = XC_LDA_X !  Slater Exchange
!  case('XC_LDA_C_WIGNER')
!    number_id =  XC_LDA_C_WIGNER !  Wigner parametrization
!  case('XC_LDA_C_RPA')
!    number_id = XC_LDA_C_RPA !  Random Phase Approximation
!  case('XC_LDA_C_HL')
!    number_id = XC_LDA_C_HL !  Hedin & Lundqvist
!  case('XC_LDA_C_GL')
!    number_id = XC_LDA_C_GL !  Gunnarson & Lundqvist
!  case('XC_LDA_C_XALPHA')
!    number_id = XC_LDA_C_XALPHA !  Slater Xalpha
!  case('XC_LDA_C_VWN')
!    number_id = XC_LDA_C_VWN !  Vosko, Wilk, & Nussair
!  case('XC_LDA_C_VWN_RPA')
!    number_id = XC_LDA_C_VWN_RPA !  Vosko, Wilk, & Nussair (RPA)
!  case('XC_LDA_C_PZ')
!    number_id = XC_LDA_C_PZ !  Perdew & Zunger
!  case('XC_LDA_C_PZ_MOD')
!    number_id = XC_LDA_C_PZ_MOD !  Perdew & Zunger (Modified)
!  case('XC_LDA_C_OB_PZ')
!    number_id = XC_LDA_C_OB_PZ !  Ortiz & Ballone (PZ)
!  case('XC_LDA_C_PW')
!    number_id = XC_LDA_C_PW !  Perdew & Wang
!  case('XC_LDA_C_PW_MOD')
!    number_id = XC_LDA_C_PW_MOD !  Perdew & Wang (Modified)
!  case('XC_LDA_C_OB_PW')
!    number_id = XC_LDA_C_OB_PW !  Ortiz & Ballone (PW)
!  case('XC_LDA_C_2D_AMGB')
!    number_id = XC_LDA_C_2D_AMGB !  Attacalite et al
!  case('XC_LDA_C_2D_PRM')
!    number_id = XC_LDA_C_2D_PRM !  Pittalis, Rasanen & Marques correlation in 2D
!  case('XC_LDA_C_vBH')
!    number_id = XC_LDA_C_vBH !  von Barth & Hedin
!  case('XC_LDA_C_1D_CSC')
!    number_id = XC_LDA_C_1D_CSC !  Casula, Sorella, and Senatore 1D correlation
!  case('XC_LDA_X_2D')
!    number_id = XC_LDA_X_2D !  Exchange in 2D
! case('XC_LDA_XC_TETER93')
!    number_id = XC_LDA_XC_TETER93 !  Teter 93 parametrization
!  case('XC_LDA_X_1D')
!    number_id = XC_LDA_X_1D !  Exchange in 1D
!  case('XC_LDA_C_ML1')
!    number_id = XC_LDA_C_ML1 !  Modified LSD (version 1) of Proynov and Salahub
!  case('XC_LDA_C_ML2')
!    number_id = XC_LDA_C_ML2 !  Modified LSD (version 2) of Proynov and Salahub
!  case('XC_LDA_C_GOMBAS')
!    number_id = XC_LDA_C_GOMBAS !  Gombas parametrization
!  case('XC_LDA_C_PW_RPA')
!    number_id = XC_LDA_C_PW_RPA !  Perdew & Wang fit of the RPA
!  case('XC_LDA_K_TF')
!    number_id = XC_LDA_K_TF !  Thomas-Fermi kinetic energy functional
!  case('XC_LDA_K_LP')
!    number_id = XC_LDA_K_LP !  Lee and Parr Gaussian ansatz
!  case('XC_GGA_X_PBE')
!    number_id = XC_GGA_X_PBE !  Perdew, Burke & Ernzerhof exchange
!  case('XC_GGA_X_PBE_R')
!    number_id = XC_GGA_X_PBE_R !  Perdew, Burke & Ernzerhof exchange (revised)
!  case('XC_GGA_X_B86')
!    number_id = XC_GGA_X_B86 !  Becke 86 Xalfa,beta,gamma
!  case('XC_GGA_X_HERMAN')
!    number_id = XC_GGA_X_HERMAN !  Herman et al original GGA
!  case('XC_GGA_X_B86_MGC')
!    number_id = XC_GGA_X_B86_MGC !  Becke 86 Xalfa,beta,gamma (with mod. grad. correction)
!  case('XC_GGA_X_B88')
!    number_id = XC_GGA_X_B88 !  Becke 88
!  case('XC_GGA_X_G96')
!    number_id = XC_GGA_X_G96 !  Gill 96
!  case('XC_GGA_X_PW86')
!    number_id = XC_GGA_X_PW86 !  Perdew & Wang 86
!  case('XC_GGA_X_PW91')
!    number_id = XC_GGA_X_PW91 !  Perdew & Wang 91
!  case('XC_GGA_X_OPTX')
!    number_id = XC_GGA_X_OPTX !  Handy & Cohen OPTX 01
!  case('XC_GGA_X_DK87_R1')
!    number_id = XC_GGA_X_DK87_R1 !  dePristo & Kress 87 (version R1)
!  case('XC_GGA_X_DK87_R2')
!    number_id = XC_GGA_X_DK87_R2 !  dePristo & Kress 87 (version R2)
!  case('XC_GGA_X_LG93')
!    number_id = XC_GGA_X_LG93 !  Lacks & Gordon 93
!  case('XC_GGA_X_FT97_A')
!    number_id = XC_GGA_X_FT97_A !  Filatov & Thiel 97 (version A)
!  case('XC_GGA_X_FT97_B')
!    number_id = XC_GGA_X_FT97_B !  Filatov & Thiel 97 (version B)
!  case('XC_GGA_X_PBE_SOL')
!    number_id = XC_GGA_X_PBE_SOL !  Perdew, Burke & Ernzerhof exchange (solids)
!  case('XC_GGA_X_RPBE')
!    number_id = XC_GGA_X_RPBE !  Hammer, Hansen & Norskov (PBE-like)
!  case('XC_GGA_X_WC')
!    number_id = XC_GGA_X_WC !  Wu & Cohen
!  case('XC_GGA_X_mPW91')
!    number_id = XC_GGA_X_mPW91 !  Modified form of PW91 by Adamo & Barone
!  case('XC_GGA_X_AM05')
!    number_id = XC_GGA_X_AM05 !  Armiento & Mattsson 05 exchange
!  case('XC_GGA_X_PBEA')
!    number_id = XC_GGA_X_PBEA !  Madsen (PBE-like)
!  case('XC_GGA_X_MPBE')
!    number_id = XC_GGA_X_MPBE !  Adamo & Barone modification to PBE
!  case('XC_GGA_X_XPBE')
!    number_id = XC_GGA_X_XPBE !  xPBE reparametrization by Xu & Goddard
!  case('XC_GGA_X_2D_B86_MGC')
!    number_id = XC_GGA_X_2D_B86_MGC !  Becke 86 MGC for 2D systems
!  case('XC_GGA_X_BAYESIAN')
!    number_id = XC_GGA_X_BAYESIAN !  Bayesian best fit for the enhancement factor
!  case('XC_GGA_X_PBE_JSJR')
!    number_id = XC_GGA_X_PBE_JSJR !  JSJR reparametrization by Pedroza, Silva & Capelle
!  case('XC_GGA_X_2D_B88')
!    number_id = XC_GGA_X_2D_B88 !  Becke 88 in 2D
!  case('XC_GGA_X_2D_B86')
!    number_id = XC_GGA_X_2D_B86 !  Becke 86 Xalfa,beta,gamma
!  case('XC_GGA_X_2D_PBE')
!    number_id = XC_GGA_X_2D_PBE !  Perdew, Burke & Ernzerhof exchange in 2D
!  case('XC_GGA_C_PBE')
!    number_id = XC_GGA_C_PBE !  Perdew, Burke & Ernzerhof correlation
!  case('XC_GGA_C_LYP')
!    number_id = XC_GGA_C_LYP !  Lee, Yang & Parr
!  case('XC_GGA_C_P86')
!    number_id = XC_GGA_C_P86 !  Perdew 86
!  case('XC_GGA_C_PBE_SOL')
!    number_id = XC_GGA_C_PBE_SOL !  Perdew, Burke & Ernzerhof correlation SOL
!  case('XC_GGA_C_PW91')
!    number_id = XC_GGA_C_PW91 !  Perdew & Wang 91
!  case('XC_GGA_C_AM05')
!    number_id = XC_GGA_C_AM05 !  Armiento & Mattsson 05 correlation
!  case('XC_GGA_C_XPBE')
!    number_id = XC_GGA_C_XPBE !  xPBE reparametrization by Xu & Goddard
!  case('XC_GGA_C_LM')
!    number_id = XC_GGA_C_LM !  Langreth and Mehl correlation
!  case('XC_GGA_C_PBE_JRGX')
!    number_id = XC_GGA_C_PBE_JRGX !  JRGX reparametrization by Pedroza, Silva & Capelle
!  case('XC_GGA_X_OPTB88_VDW')
!    number_id = XC_GGA_X_OPTB88_VDW !  Becke 88 reoptimized to be used with vdW functional of Dion et al
!  case('XC_GGA_X_PBEK1_VDW')
!    number_id = XC_GGA_X_PBEK1_VDW !  PBE reparametrization for vdW
!  case('XC_GGA_X_OPTPBE_VDW')
!    number_id = XC_GGA_X_OPTPBE_VDW !  PBE reparametrization for vdW
!  case('XC_GGA_X_RGE2')
!    number_id = XC_GGA_X_RGE2 !  Regularized PBE
!  case('XC_GGA_C_RGE2')
!    number_id = XC_GGA_C_RGE2 !  Regularized PBE
!  case('XC_GGA_X_RPW86')
!    number_id = XC_GGA_X_RPW86 !  refitted Perdew & Wang 86
!  case('XC_GGA_X_KT1')
!    number_id = XC_GGA_X_KT1 !  Keal and Tozer version 1
!  case('XC_GGA_XC_KT2')
!    number_id = XC_GGA_XC_KT2 !  Keal and Tozer version 2
!  case('XC_GGA_C_WL')
!    number_id = XC_GGA_C_WL !  Wilson & Levy
!  case('XC_GGA_C_WI')
!    number_id = XC_GGA_C_WI !  Wilson & Ivanov
!  case('XC_GGA_X_MB88')
!    number_id = XC_GGA_X_MB88 !  Modified Becke 88 for proton transfer
!  case('XC_GGA_X_SOGGA')
!    number_id = XC_GGA_X_SOGGA !  Second-order generalized gradient approximation
!  case('XC_GGA_X_SOGGA11')
!    number_id = XC_GGA_X_SOGGA11 !  Second-order generalized gradient approximation 2011
!  case('XC_GGA_C_SOGGA11')
!    number_id = XC_GGA_C_SOGGA11 !  Second-order generalized gradient approximation 2011
!  case('XC_GGA_C_WI0')
!    number_id = XC_GGA_C_WI0 !  Wilson & Ivanov initial version
!  case('XC_GGA_XC_TH1')
!    number_id = XC_GGA_XC_TH1 !  Tozer and Handy v. 1
!  case('XC_GGA_XC_TH2')
!    number_id = XC_GGA_XC_TH2 !  Tozer and Handy v. 2
!  case('XC_GGA_XC_TH3')
!    number_id = XC_GGA_XC_TH3 !  Tozer and Handy v. 3
!  case('XC_GGA_XC_TH4')
!    number_id = XC_GGA_XC_TH4 !  Tozer and Handy v. 4
!  case('XC_GGA_X_C09X')
!    number_id = XC_GGA_X_C09X !  C09x to be used with the VdW of Rutgers-Chalmers
!  case('XC_GGA_C_SOGGA11_X')
!    number_id = XC_GGA_C_SOGGA11_X !  To be used with hyb_gga_x_SOGGA11-X
!  case('XC_GGA_X_LB')
!    number_id = XC_GGA_X_LB !  van Leeuwen & Baerends
!  case('XC_GGA_XC_HCTH_93')
!    number_id = XC_GGA_XC_HCTH_93 !  HCTH functional fitted to  93 molecules
!  case('XC_GGA_XC_HCTH_120')
!    number_id = XC_GGA_XC_HCTH_120 !  HCTH functional fitted to 120 molecules
!  case('XC_GGA_XC_HCTH_147')
!    number_id = XC_GGA_XC_HCTH_147 !  HCTH functional fitted to 147 molecules
!  case('XC_GGA_XC_HCTH_407')
!    number_id = XC_GGA_XC_HCTH_407 !  HCTH functional fitted to 147 molecules
!  case('XC_GGA_XC_EDF1')
!    number_id = XC_GGA_XC_EDF1 !  Empirical functionals from Adamson, Gill, and Pople
!  case('XC_GGA_XC_XLYP')
!    number_id = XC_GGA_XC_XLYP !  XLYP functional
!! These no longer exist in LibXC v3.0.0
!!  case('XC_GGA_XC_B97')
!!    number_id = XC_GGA_XC_B97  !  Becke 97
!!  case('XC_GGA_XC_B97_1')
!!    number_id = XC_GGA_XC_B97_1  !  Becke 97-1
!!  case('XC_GGA_XC_B97_2')
!!    number_id = XC_GGA_XC_B97_2 !  Becke 97-2
!  case('XC_GGA_XC_B97_D')
!    number_id = XC_GGA_XC_B97_D !  Grimme functional to be used with C6 vdW term
!!  case('XC_GGA_XC_B97_K')
!!    number_id = XC_GGA_XC_B97_K !  Boese-Martin for Kinetics
!!  case('XC_GGA_XC_B97_3')
!!    number_id = XC_GGA_XC_B97_3 !  Becke 97-3
!  case('XC_GGA_XC_PBE1W')
!    number_id = XC_GGA_XC_PBE1W !  Functionals fitted for water
!  case('XC_GGA_XC_MPWLYP1W')
!    number_id = XC_GGA_XC_MPWLYP1W !  Functionals fitted for water
!  case('XC_GGA_XC_PBELYP1W')
!    number_id = XC_GGA_XC_PBELYP1W !  Functionals fitted for water
!!  case('XC_GGA_XC_SB98_1a')
!!    number_id = XC_GGA_XC_SB98_1a !  Schmider-Becke 98 parameterization 1a
!!  case('XC_GGA_XC_SB98_1b')
!!    number_id = XC_GGA_XC_SB98_1b !  Schmider-Becke 98 parameterization 1b
!!  case('XC_GGA_XC_SB98_1c')
!!    number_id = XC_GGA_XC_SB98_1c !  Schmider-Becke 98 parameterization 1c
!!  case('XC_GGA_XC_SB98_2a')
!!    number_id = XC_GGA_XC_SB98_2a !  Schmider-Becke 98 parameterization 2a
!!  case('XC_GGA_XC_SB98_2b')
!!    number_id = XC_GGA_XC_SB98_2b !  Schmider-Becke 98 parameterization 2b
!!  case('XC_GGA_XC_SB98_2c')
!!    number_id = XC_GGA_XC_SB98_2c  !  Schmider-Becke 98 parameterization 2c
!  case('XC_GGA_X_LBM')
!    number_id = XC_GGA_X_LBM !  van Leeuwen & Baerends modified
!  case('XC_GGA_X_OL2')
!    number_id = XC_GGA_X_OL2 !  Exchange form based on Ou-Yang and Levy v.2
!  case('XC_GGA_X_APBE')
!    number_id = XC_GGA_X_APBE !  mu fixed from the semiclassical neutral atom
!  case('XC_GGA_K_APBE')
!    number_id = XC_GGA_K_APBE !  mu fixed from the semiclassical neutral atom
!  case('XC_GGA_C_APBE')
!    number_id = XC_GGA_C_APBE !  mu fixed from the semiclassical neutral atom
!  case('XC_GGA_K_TW1')
!    number_id = XC_GGA_K_TW1 !  Tran and Wesolowski set 1 (Table II)
!  case('XC_GGA_K_TW2')
!    number_id = XC_GGA_K_TW2 !  Tran and Wesolowski set 2 (Table II)
!  case('XC_GGA_K_TW3')
!    number_id = XC_GGA_K_TW3 !  Tran and Wesolowski set 3 (Table II)
!  case('XC_GGA_K_TW4')
!    number_id = XC_GGA_K_TW4 !  Tran and Wesolowski set 4 (Table II)
!  case('XC_GGA_X_HTBS')
!    number_id = XC_GGA_X_HTBS !  Haas, Tran, Blaha, and Schwarz
!  case('XC_GGA_X_AIRY')
!    number_id = XC_GGA_X_AIRY !  Constantin et al based on the Airy gas
!  case('XC_GGA_X_LAG')
!    number_id = XC_GGA_X_LAG !  Local Airy Gas
!  case('XC_GGA_XC_MOHLYP')
!    number_id = XC_GGA_XC_MOHLYP !  Functional for organometallic chemistry
!  case('XC_GGA_XC_MOHLYP2')
!    number_id = XC_GGA_XC_MOHLYP2 !  Functional for barrier heights
!  case('XC_GGA_XC_TH_FL')
!    number_id = XC_GGA_XC_TH_FL !  Tozer and Handy v. FL
!  case('XC_GGA_XC_TH_FC')
!    number_id = XC_GGA_XC_TH_FC !  Tozer and Handy v. FC
!  case('XC_GGA_XC_TH_FCFO')
!    number_id = XC_GGA_XC_TH_FCFO !  Tozer and Handy v. FCFO
!  case('XC_GGA_XC_TH_FCO')
!    number_id = XC_GGA_XC_TH_FCO  !  Tozer and Handy v. FCO
!  case('XC_GGA_K_VW')
!    number_id = XC_GGA_K_VW !  von Weiszaecker functional
!  case('XC_GGA_K_GE2')
!    number_id = XC_GGA_K_GE2 !  Second-order gradient expansion (l = 1/9)
!  case('XC_GGA_K_GOLDEN')
!    number_id = XC_GGA_K_GOLDEN !  TF-lambda-vW form by Golden (l = 13/45)
!  case('XC_GGA_K_YT65')
!    number_id = XC_GGA_K_YT65 !  TF-lambda-vW form by Yonei and Tomishima (l = 1/5)
!  case('XC_GGA_K_BALTIN')
!    number_id = XC_GGA_K_BALTIN !  TF-lambda-vW form by Baltin (l = 5/9)
!  case('XC_GGA_K_LIEB')
!    number_id = XC_GGA_K_LIEB !  TF-lambda-vW form by Lieb (l = 0.185909191)
!  case('XC_GGA_K_ABSR1')
!    number_id = XC_GGA_K_ABSR1 !  gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]
!  case('XC_GGA_K_ABSR2')
!    number_id = XC_GGA_K_ABSR2 !  gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]
!  case('XC_GGA_K_GR')
!    number_id = XC_GGA_K_GR !  gamma-TFvW form by Gázquez and Robles
!  case('XC_GGA_K_LUDENA')
!    number_id = XC_GGA_K_LUDENA !  gamma-TFvW form by Ludeña
!  case('XC_GGA_K_GP85')
!    number_id = XC_GGA_K_GP85 !  gamma-TFvW form by Ghosh and Parr
!  case('XC_GGA_K_PEARSON')
!    number_id = XC_GGA_K_PEARSON !  Pearson
!  case('XC_GGA_K_OL1')
!    number_id = XC_GGA_K_OL1 !  Ou-Yang and Levy v.1
!  case('XC_GGA_K_OL2')
!    number_id = XC_GGA_K_OL2 !  Ou-Yang and Levy v.2
!  case('XC_GGA_K_FR_B88')
!    number_id = XC_GGA_K_FR_B88 !  Fuentealba & Reyes (B88 version)
!  case('XC_GGA_K_FR_PW86')
!    number_id = XC_GGA_K_FR_PW86 !  Fuentealba & Reyes (PW86 version)
!  case('XC_GGA_K_DK')
!    number_id = XC_GGA_K_DK !  DePristo and Kress
!  case('XC_GGA_K_PERDEW')
!    number_id = XC_GGA_K_PERDEW !  Perdew
!  case('XC_GGA_K_VSK')
!    number_id = XC_GGA_K_PERDEW !  Vitos, Skriver, and Kollar
!  case('XC_GGA_K_VJKS')
!    number_id = XC_GGA_K_VJKS !  Vitos, Johansson, Kollar, and Skriver
!  case('XC_GGA_K_ERNZERHOF')
!    number_id = XC_GGA_K_ERNZERHOF !  Ernzerhof
!  case('XC_GGA_K_LC94')
!    number_id = XC_GGA_K_LC94 !  Lembarki & Chermette
!  case('XC_GGA_K_LLP')
!    number_id = XC_GGA_K_LLP !  Lee, Lee & Parr
!  case('XC_GGA_K_THAKKAR')
!    number_id = XC_GGA_K_THAKKAR !  Thakkar 1992
!  case('XC_HYB_GGA_XC_B3PW91')
!    number_id = XC_HYB_GGA_XC_B3PW91 !  The original hybrid proposed by Becke
!  case('XC_HYB_GGA_XC_B3LYP')
!    number_id = XC_HYB_GGA_XC_B3LYP !  The (in)famous B3LYP
!  case('XC_HYB_GGA_XC_B3P86')
!    number_id = XC_HYB_GGA_XC_B3P86 !  Perdew 86 hybrid similar to B3PW91
!  case('XC_HYB_GGA_XC_O3LYP')
!    number_id = XC_HYB_GGA_XC_O3LYP !  hybrid using the optx functional
!  case('XC_HYB_GGA_XC_mPW1K')
!    number_id = XC_HYB_GGA_XC_mPW1K !  mixture of mPW91 and PW91 optimized for kinetics
!  case('XC_HYB_GGA_XC_PBEH')
!    number_id = XC_HYB_GGA_XC_PBEH !  aka PBE0 or PBE1PBE
!  case('XC_HYB_GGA_XC_B97')
!    number_id = XC_HYB_GGA_XC_B97 !  Becke 97
!  case('XC_HYB_GGA_XC_B97_1')
!    number_id = XC_HYB_GGA_XC_B97_1 !  Becke 97-1
!  case('XC_HYB_GGA_XC_B97_2')
!    number_id = XC_HYB_GGA_XC_B97_2 !  Becke 97-2
!  case('XC_HYB_GGA_XC_X3LYP')
!    number_id = XC_HYB_GGA_XC_X3LYP !  maybe the best hybrid
!  case('XC_HYB_GGA_XC_B1WC')
!    number_id = XC_HYB_GGA_XC_B1WC !  Becke 1-parameter mixture of WC and PBE
!  case('XC_HYB_GGA_XC_B97_K')
!    number_id = XC_HYB_GGA_XC_B97_K !  Boese-Martin for Kinetics
!  case('XC_HYB_GGA_XC_B97_3')
!    number_id = XC_HYB_GGA_XC_B97_3 !  Becke 97-3
!  case('XC_HYB_GGA_XC_mPW3PW')
!    number_id = XC_HYB_GGA_XC_mPW3PW !  mixture with the mPW functional
!  case('XC_HYB_GGA_XC_B1LYP')
!    number_id = XC_HYB_GGA_XC_B1LYP !  Becke 1-parameter mixture of B88 and LYP
!  case('XC_HYB_GGA_XC_B1PW91')
!    number_id = XC_HYB_GGA_XC_B1PW91 !  Becke 1-parameter mixture of B88 and PW91
!  case('XC_HYB_GGA_XC_mPW1PW')
!    number_id = XC_HYB_GGA_XC_mPW1PW !  Becke 1-parameter mixture of mPW91 and PW91
!  case('XC_HYB_GGA_XC_mPW3LYP')
!    number_id = XC_HYB_GGA_XC_mPW3LYP !  mixture of mPW and LYP
!  case('XC_HYB_GGA_XC_SB98_1a')
!    number_id = XC_HYB_GGA_XC_SB98_1a !  Schmider-Becke 98 parameterization 1a
!  case('XC_HYB_GGA_XC_SB98_1b')
!    number_id = XC_HYB_GGA_XC_SB98_1b !  Schmider-Becke 98 parameterization 1b
!  case('XC_HYB_GGA_XC_SB98_1c')
!    number_id = XC_HYB_GGA_XC_SB98_1c !  Schmider-Becke 98 parameterization 1c
!  case('XC_HYB_GGA_XC_SB98_2a')
!    number_id = XC_HYB_GGA_XC_SB98_2a !  Schmider-Becke 98 parameterization 2a
!  case('XC_HYB_GGA_XC_SB98_2b')
!    number_id = XC_HYB_GGA_XC_SB98_2b !  Schmider-Becke 98 parameterization 2b
!  case('XC_HYB_GGA_XC_SB98_2c')
!    number_id = XC_HYB_GGA_XC_SB98_2c !  Schmider-Becke 98 parameterization 2c
!  case('XC_HYB_GGA_X_SOGGA11_X')
!    number_id = XC_HYB_GGA_X_SOGGA11_X !  Hybrid based on SOGGA11 form
!  case('XC_MGGA_X_LTA')
!    number_id = XC_MGGA_X_LTA  !  Local tau approximation of Ernzerhof & Scuseria
!  case('XC_MGGA_X_TPSS')
!    number_id = XC_MGGA_X_TPSS !  Perdew, Tao, Staroverov & Scuseria exchange
!  case('XC_MGGA_X_M06L')
!!    number_id = XC_MGGA_X_M06L !  Zhao, Truhlar exchange
!    number_id = XC_MGGA_X_M06_L ! Updated in version 3.0.0
!  case('XC_MGGA_X_GVT4')
!    number_id = XC_MGGA_X_GVT4 !  GVT4 from Van Voorhis and Scuseria (exchange part)
!  case('XC_MGGA_X_TAU_HCTH')
!    number_id = XC_MGGA_X_TAU_HCTH !  tau-HCTH from Boese and Handy
!  case('XC_MGGA_X_BR89')
!    number_id = XC_MGGA_X_BR89 !  Becke-Roussel 89
!  case('XC_MGGA_X_BJ06')
!    number_id = XC_MGGA_X_BJ06 !  Becke & Johnson correction to Becke-Roussel 89
!  case('XC_MGGA_X_TB09')
!    number_id = XC_MGGA_X_TB09 !  Tran & Blaha correction to Becke & Johnson
!  case('XC_MGGA_X_RPP09')
!    number_id = XC_MGGA_X_RPP09 !  Rasanen, Pittalis, and Proetto correction to Becke & Johnson
!  case('XC_MGGA_X_2D_PRHG07')
!    number_id = XC_MGGA_X_2D_PRHG07 !  Pittalis, Rasanen, Helbig, Gross Exchange Functional
!  case('XC_MGGA_X_2D_PRHG07_PRP10')
!    number_id = XC_MGGA_X_2D_PRHG07_PRP10 !  PRGH07 with PRP10 correction
!  case('XC_MGGA_C_TPSS')
!    number_id = XC_MGGA_C_TPSS !  Perdew, Tao, Staroverov & Scuseria correlation
!  case('XC_MGGA_C_VSXC')
!    number_id = XC_MGGA_C_VSXC !  VSxc from Van Voorhis and Scuseria (correlation part)
!  end select
!end subroutine determine_func_number
