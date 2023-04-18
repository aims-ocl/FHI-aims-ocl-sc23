! LibXC Stubs. Updated by AJL, August 2016

module xc_f03_lib_m
  
  use localorb_io
  use mpi_tasks
  
  implicit none

  
  integer, parameter, public :: &
    XC_UNPOLARIZED = 1, & ! Spin unpolarized
    XC_POLARIZED = 2 ! Spin polarized

    integer, parameter, public :: &
    XC_NON_RELATIVISTIC = 0, & ! Functional includes or not relativistic
    XC_RELATIVISTIC = 1 ! corrections. Only available in some functionals.

  ! Kinds
  integer, parameter, public :: &
    XC_EXCHANGE = 0, &
    XC_CORRELATION = 1, &
    XC_EXCHANGE_CORRELATION = 2, &
    XC_KINETIC = 3

  ! Families of xc functionals
  integer, parameter, public :: &
    XC_FAMILY_UNKNOWN = -1, &
    XC_FAMILY_NONE = 0, &
    XC_FAMILY_LDA = 1, &
    XC_FAMILY_GGA = 2, &
    XC_FAMILY_MGGA = 4, &
    XC_FAMILY_LCA = 8, &
    XC_FAMILY_OEP = 16, &
    XC_FAMILY_HYB_GGA = 32, &
    XC_FAMILY_HYB_MGGA = 64

  ! Stubs for commonly used functionals
  integer :: XC_LDA_X, XC_LDA_C_PW, XC_GGA_X_PBE, XC_GGA_C_PBE, XC_LDA_C_VWN, &
       & XC_GGA_X_PW91, XC_GGA_C_PW91, XC_GGA_X_RPBE, XC_LDA_C_PZ, &
       & XC_GGA_X_B88, XC_GGA_C_LYP, XC_HYB_GGA_XC_B3LYP

  type :: xc_f03_func_info_t
    private
    integer, pointer :: ptr
  end type xc_f03_func_info_t

  type :: xc_f03_func_t
    private
    integer, pointer :: ptr 
  end type xc_f03_func_t

  type :: xc_f03_func_reference_t
    private
    integer, pointer :: ptr
  end type xc_f03_func_reference_t

    contains
 
  ! Subroutines
  subroutine xc_f03_func_init(p, functional, nspin)
    use localorb_io
    use mpi_tasks
    implicit none

    type(xc_f03_func_t), intent(inout) :: p
    integer, intent(in) :: functional
    integer, intent(in) :: nspin

    ! Output an error
    character*128              :: info_str
    write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
    call localorb_info(info_str)
    write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
    call localorb_info(info_str)
    call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
  end subroutine xc_f03_func_init

  subroutine xc_f03_lda_vxc(p, np, rho, vrho)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*)
    real*8 :: vrho(*)
    character*128 :: info_str
    call localorb_info('***** This is a call to libxc_stubs and should not &
         &have happened!')
    call aims_stop(' General STOP. Check compile time variables. &
         &This is a dead end here.')
  end subroutine xc_f03_lda_vxc

  subroutine xc_f03_lda_exc_vxc(p, np, rho, zk, vrho)
    use localorb_io
    use mpi_tasks
    implicit none

    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*)
    real*8 :: zk(*), vrho(*)
  
    ! Output an error
    character*128              :: info_str
    write(info_str,'(a)') '***** This is a call to libxc_stubs and should not have happened!'
    call localorb_info(info_str)
    call aims_stop(' General STOP. Check compile time variables. This is a dead end here.')
  end subroutine xc_f03_lda_exc_vxc

  subroutine xc_f03_lda_fxc(p, np, rho, v2rho2)
    use localorb_io
    use mpi_tasks
    implicit none

    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*)
    real*8 :: v2rho2(*)

    ! Output an error
    character*128              :: info_str
    write(info_str,'(a)') '***** This is a call to libxc_stubs and should not have happened!'
    call localorb_info(info_str)
    call aims_stop(' General STOP. Check compile time variables. This is a dead end here.')
  end subroutine xc_f03_lda_fxc

  subroutine xc_f03_gga_vxc(p, np, rho, sigma, vrho, vsigma)
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*), sigma(*)
    real*8 :: vrho(*), vsigma(*)
    call localorb_info('***** This is a call to libxc_stubs and should not &
         &have happened!')
    call aims_stop(' General STOP. Check compile time variables. This is a &
         &dead end here.')
  end subroutine xc_f03_gga_vxc

  subroutine xc_f03_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
    use localorb_io
    use mpi_tasks
    implicit none

    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*), sigma(*)
    real*8 :: zk(*), vrho(*), vsigma(*)

    ! Output an error
    character*128              :: info_str
    write(info_str,'(a)') '***** This is a call to libxc_stubs and should not have happened!'
    call localorb_info(info_str)
    call aims_stop(' General STOP. Check compile time variables. This is a dead end here.')
  end subroutine xc_f03_gga_exc_vxc

  subroutine xc_f03_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    implicit none
    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*), sigma(*)
    real*8 :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    call localorb_info('***** This is a call to libxc_stubs and should not &
         &have happened!')
    call aims_stop(' General STOP. Check compile time variables. This is a &
         &dead end here.')
  end subroutine xc_f03_gga_fxc

  subroutine xc_f03_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    use localorb_io
    use mpi_tasks
    implicit none

    type(xc_f03_func_t), intent(in) :: p
    integer, intent(in) :: np
    real*8, intent(in) :: rho(*), sigma(*), lapl(*), tau(*)
    real*8 :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)

    ! Output an error
    character*128              :: info_str
    write(info_str,'(a)') '***** This is a call to libxc_stubs and should not have happened!'
    call localorb_info(info_str)
    call aims_stop(' General STOP. Check compile time variables. This is a dead end here.')
  end subroutine xc_f03_mgga_exc_vxc

  subroutine xc_f03_func_end(p)
    use localorb_io
    use mpi_tasks
    implicit none

    type(xc_f03_func_t), intent(inout) :: p

    ! Output an error
    character*128              :: info_str
    write(info_str,'(a)') '***** This is a call to libxc_stubs and should not have happened!'
    call localorb_info(info_str)
    call aims_stop(' General STOP. Check compile time variables. This is a dead end here.')
  end subroutine xc_f03_func_end

  ! Functions
  integer function xc_f03_family_from_id(id, family, number)
    integer, intent(in) :: id
    integer, optional, target :: family, number
    xc_f03_family_from_id = 0
  end function xc_f03_family_from_id

  integer function xc_f03_functional_get_number(func_string) result(number)
    character(len=*), intent(in) :: func_string
    number = 0
  end function xc_f03_functional_get_number

  type(xc_f03_func_info_t) function xc_f03_func_get_info(p) result(info)
    type(xc_f03_func_t), intent(in) :: p
    type(xc_f03_func_info_t) :: dummy
    info = dummy
  end function xc_f03_func_get_info

  integer function xc_f03_func_info_get_family(info) result(family)
    type(xc_f03_func_info_t), intent(in) :: info
    family = 0
  end function xc_f03_func_info_get_family

  integer function xc_f03_func_info_get_kind(info) result(kind)
    type(xc_f03_func_info_t), intent(in) :: info
    kind = 0
  end function xc_f03_func_info_get_kind

  character(len=128) function xc_f03_func_info_get_name(info) result(name)
    type(xc_f03_func_info_t), intent(in) :: info
    name = ""
  end function xc_f03_func_info_get_name

! AJL, November 2017
! Retired as of v3.0.0
!  character(len=120) function xc_f03_func_info_get_refs(info, number) result(ref)
!    type(xc_f03_func_info_t), intent(in) :: info
!    integer :: number ! number of the reference. Must be 0 in the first call
!    ref = ""
!  end function xc_f03_func_info_get_refs

! New as of V4.0.2
  type(xc_f03_func_reference_t) function xc_f03_func_info_get_references(info,number) result(reference)
    type(xc_f03_func_info_t), intent(in) :: info
    integer, intent(inout) :: number ! number of the reference. Must be 0 in the first call
    type(xc_f03_func_reference_t) :: dummy
    reference = dummy
  end function xc_f03_func_info_get_references

  character(len=120) function xc_f03_func_reference_get_ref(reference) result(ref)
    type(xc_f03_func_reference_t), intent(in) :: reference
    ref = ""
  end function xc_f03_func_reference_get_ref
! End, AJL

  real*8 function xc_f03_hyb_exx_coef(p) result(coef)
    type(xc_f03_func_t), intent(in) :: p
    coef = 0.d0
  end function xc_f03_hyb_exx_coef

end module xc_f03_lib_m

module xc_f90_types_m



  integer, public, parameter :: xc_f90_kind = selected_real_kind(14)


  type xc_f90_pointer_t
    private
    integer, pointer :: buffer
  end type xc_f90_pointer_t

end module xc_f90_types_m

module libxc_funcs_m
  implicit none

  public

  integer, parameter :: XC_LDA_X                       =   1  ! Exchange
  integer, parameter :: XC_LDA_C_WIGNER                =   2  ! Wigner parametrization
  integer, parameter :: XC_LDA_C_RPA                   =   3  ! Random Phase Approximation
  integer, parameter :: XC_LDA_C_HL                    =   4  ! Hedin & Lundqvist
  integer, parameter :: XC_LDA_C_GL                    =   5  ! Gunnarson & Lundqvist
  integer, parameter :: XC_LDA_C_XALPHA                =   6  ! Slater Xalpha
  integer, parameter :: XC_LDA_C_VWN                   =   7  ! Vosko, Wilk, & Nusair (5)
  integer, parameter :: XC_LDA_C_VWN_RPA               =   8  ! Vosko, Wilk, & Nusair (RPA)
  integer, parameter :: XC_LDA_C_PZ                    =   9  ! Perdew & Zunger
  integer, parameter :: XC_LDA_C_PZ_MOD                =  10  ! Perdew & Zunger (Modified)
  integer, parameter :: XC_LDA_C_OB_PZ                 =  11  ! Ortiz & Ballone (PZ)
  integer, parameter :: XC_LDA_C_PW                    =  12  ! Perdew & Wang
  integer, parameter :: XC_LDA_C_PW_MOD                =  13  ! Perdew & Wang (Modified)
  integer, parameter :: XC_LDA_C_OB_PW                 =  14  ! Ortiz & Ballone (PW)
  integer, parameter :: XC_LDA_C_2D_AMGB               =  15  ! Attaccalite et al
  integer, parameter :: XC_LDA_C_2D_PRM                =  16  ! Pittalis, Rasanen & Marques correlation in 2D
  integer, parameter :: XC_LDA_C_vBH                   =  17  ! von Barth & Hedin
  integer, parameter :: XC_LDA_C_1D_CSC                =  18  ! Casula, Sorella, and Senatore 1D correlation
  integer, parameter :: XC_LDA_X_2D                    =  19  ! Exchange in 2D
  integer, parameter :: XC_LDA_XC_TETER93              =  20  ! Teter 93 parametrization
  integer, parameter :: XC_LDA_X_1D                    =  21  ! Exchange in 1D
  integer, parameter :: XC_LDA_C_ML1                   =  22  ! Modified LSD (version 1) of Proynov and Salahub
  integer, parameter :: XC_LDA_C_ML2                   =  23  ! Modified LSD (version 2) of Proynov and Salahub
  integer, parameter :: XC_LDA_C_GOMBAS                =  24  ! Gombas parametrization
  integer, parameter :: XC_LDA_C_PW_RPA                =  25  ! Perdew & Wang fit of the RPA
  integer, parameter :: XC_LDA_C_1D_LOOS               =  26  ! P-F Loos correlation LDA
  integer, parameter :: XC_LDA_C_RC04                  =  27  ! Ragot-Cortona
  integer, parameter :: XC_LDA_C_VWN_1                 =  28  ! Vosko, Wilk, & Nusair (1)
  integer, parameter :: XC_LDA_C_VWN_2                 =  29  ! Vosko, Wilk, & Nusair (2)
  integer, parameter :: XC_LDA_C_VWN_3                 =  30  ! Vosko, Wilk, & Nusair (3)
  integer, parameter :: XC_LDA_C_VWN_4                 =  31  ! Vosko, Wilk, & Nusair (4)
  integer, parameter :: XC_LDA_XC_ZLP                  =  43  ! Zhao, Levy & Parr, Eq. (20)
  integer, parameter :: XC_LDA_K_TF                    =  50  ! Thomas-Fermi kinetic energy functional
  integer, parameter :: XC_LDA_K_LP                    =  51  ! Lee and Parr Gaussian ansatz
  integer, parameter :: XC_LDA_XC_KSDT                 = 259  ! Karasiev et al. parametrization
  integer, parameter :: XC_GGA_X_GAM                   =  32  ! GAM functional from Minnesota
  integer, parameter :: XC_GGA_C_GAM                   =  33  ! GAM functional from Minnesota
  integer, parameter :: XC_GGA_X_HCTH_A                =  34  ! HCTH-A
  integer, parameter :: XC_GGA_X_EV93                  =  35  ! Engel and Vosko
  integer, parameter :: XC_GGA_X_BGCP                  =  38  ! Burke, Cancio, Gould, and Pittalis
  integer, parameter :: XC_GGA_C_BGCP                  =  39  ! Burke, Cancio, Gould, and Pittalis
  integer, parameter :: XC_GGA_X_LAMBDA_OC2_N          =  40  ! lambda_OC2(N) version of PBE
  integer, parameter :: XC_GGA_X_B86_R                 =  41  ! Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
  integer, parameter :: XC_GGA_X_LAMBDA_CH_N           =  44  ! lambda_CH(N) version of PBE
  integer, parameter :: XC_GGA_X_LAMBDA_LO_N           =  45  ! lambda_LO(N) version of PBE
  integer, parameter :: XC_GGA_X_HJS_B88_V2            =  46  ! HJS screened exchange corrected B88 version
  integer, parameter :: XC_GGA_C_Q2D                   =  47  ! Chiodo et al
  integer, parameter :: XC_GGA_X_Q2D                   =  48  ! Chiodo et al
  integer, parameter :: XC_GGA_X_PBE_MOL               =  49  ! Del Campo, Gazquez, Trickey and Vela (PBE-like)
  integer, parameter :: XC_GGA_K_TFVW                  =  52  ! Thomas-Fermi plus von Weiszaecker correction
  integer, parameter :: XC_GGA_K_REVAPBEINT            =  53  ! interpolated version of REVAPBE
  integer, parameter :: XC_GGA_K_APBEINT               =  54  ! interpolated version of APBE
  integer, parameter :: XC_GGA_K_REVAPBE               =  55  ! revised APBE
  integer, parameter :: XC_GGA_X_AK13                  =  56  ! Armiento & Kuemmel 2013
  integer, parameter :: XC_GGA_K_MEYER                 =  57  ! Meyer,  Wang, and Young
  integer, parameter :: XC_GGA_X_LV_RPW86              =  58  ! Berland and Hyldgaard
  integer, parameter :: XC_GGA_X_PBE_TCA               =  59  ! PBE revised by Tognetti et al
  integer, parameter :: XC_GGA_X_PBEINT                =  60  ! PBE for hybrid interfaces
  integer, parameter :: XC_GGA_C_ZPBEINT               =  61  ! spin-dependent gradient correction to PBEint
  integer, parameter :: XC_GGA_C_PBEINT                =  62  ! PBE for hybrid interfaces
  integer, parameter :: XC_GGA_C_ZPBESOL               =  63  ! spin-dependent gradient correction to PBEsol
  integer, parameter :: XC_GGA_XC_OPBE_D               =  65  ! oPBE_D functional of Goerigk and Grimme
  integer, parameter :: XC_GGA_XC_OPWLYP_D             =  66  ! oPWLYP-D functional of Goerigk and Grimme
  integer, parameter :: XC_GGA_XC_OBLYP_D              =  67  ! oBLYP-D functional of Goerigk and Grimme
  integer, parameter :: XC_GGA_X_VMT84_GE              =  68  ! VMT{8,4} with constraint satisfaction with mu = mu_GE
  integer, parameter :: XC_GGA_X_VMT84_PBE             =  69  ! VMT{8,4} with constraint satisfaction with mu = mu_PBE
  integer, parameter :: XC_GGA_X_VMT_GE                =  70  ! Vela, Medel, and Trickey with mu = mu_GE
  integer, parameter :: XC_GGA_X_VMT_PBE               =  71  ! Vela, Medel, and Trickey with mu = mu_PBE
  integer, parameter :: XC_GGA_C_N12_SX                =  79  ! N12-SX functional from Minnesota
  integer, parameter :: XC_GGA_C_N12                   =  80  ! N12 functional from Minnesota
  integer, parameter :: XC_GGA_X_N12                   =  82  ! N12 functional from Minnesota
  integer, parameter :: XC_GGA_C_REGTPSS               =  83  ! Regularized TPSS correlation (ex-VPBE)
  integer, parameter :: XC_GGA_C_OP_XALPHA             =  84  ! one-parameter progressive functional (XALPHA version)
  integer, parameter :: XC_GGA_C_OP_G96                =  85  ! one-parameter progressive functional (G96 version)
  integer, parameter :: XC_GGA_C_OP_PBE                =  86  ! one-parameter progressive functional (PBE version)
  integer, parameter :: XC_GGA_C_OP_B88                =  87  ! one-parameter progressive functional (B88 version)
  integer, parameter :: XC_GGA_C_FT97                  =  88  ! Filatov & Thiel correlation
  integer, parameter :: XC_GGA_C_SPBE                  =  89  ! PBE correlation to be used with the SSB exchange
  integer, parameter :: XC_GGA_X_SSB_SW                =  90  ! Swarta, Sola and Bickelhaupt correction to PBE
  integer, parameter :: XC_GGA_X_SSB                   =  91  ! Swarta, Sola and Bickelhaupt
  integer, parameter :: XC_GGA_X_SSB_D                 =  92  ! Swarta, Sola and Bickelhaupt dispersion
  integer, parameter :: XC_GGA_XC_HCTH_407P            =  93  ! HCTH/407+
  integer, parameter :: XC_GGA_XC_HCTH_P76             =  94  ! HCTH p=7/6
  integer, parameter :: XC_GGA_XC_HCTH_P14             =  95  ! HCTH p=1/4
  integer, parameter :: XC_GGA_XC_B97_GGA1             =  96  ! Becke 97 GGA-1
  integer, parameter :: XC_GGA_C_HCTH_A                =  97  ! HCTH-A
  integer, parameter :: XC_GGA_X_BPCCAC                =  98  ! BPCCAC (GRAC for the energy)
  integer, parameter :: XC_GGA_C_REVTCA                =  99  ! Tognetti, Cortona, Adamo (revised)
  integer, parameter :: XC_GGA_C_TCA                   = 100  ! Tognetti, Cortona, Adamo
  integer, parameter :: XC_GGA_X_PBE                   = 101  ! Perdew, Burke & Ernzerhof exchange
  integer, parameter :: XC_GGA_X_PBE_R                 = 102  ! Perdew, Burke & Ernzerhof exchange (revised)
  integer, parameter :: XC_GGA_X_B86                   = 103  ! Becke 86 Xalpha,beta,gamma
  integer, parameter :: XC_GGA_X_HERMAN                = 104  ! Herman et al original GGA
  integer, parameter :: XC_GGA_X_B86_MGC               = 105  ! Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
  integer, parameter :: XC_GGA_X_B88                   = 106  ! Becke 88
  integer, parameter :: XC_GGA_X_G96                   = 107  ! Gill 96
  integer, parameter :: XC_GGA_X_PW86                  = 108  ! Perdew & Wang 86
  integer, parameter :: XC_GGA_X_PW91                  = 109  ! Perdew & Wang 91
  integer, parameter :: XC_GGA_X_OPTX                  = 110  ! Handy & Cohen OPTX 01
  integer, parameter :: XC_GGA_X_DK87_R1               = 111  ! dePristo & Kress 87 (version R1)
  integer, parameter :: XC_GGA_X_DK87_R2               = 112  ! dePristo & Kress 87 (version R2)
  integer, parameter :: XC_GGA_X_LG93                  = 113  ! Lacks & Gordon 93
  integer, parameter :: XC_GGA_X_FT97_A                = 114  ! Filatov & Thiel 97 (version A)
  integer, parameter :: XC_GGA_X_FT97_B                = 115  ! Filatov & Thiel 97 (version B)
  integer, parameter :: XC_GGA_X_PBE_SOL               = 116  ! Perdew, Burke & Ernzerhof exchange (solids)
  integer, parameter :: XC_GGA_X_RPBE                  = 117  ! Hammer, Hansen & Norskov (PBE-like)
  integer, parameter :: XC_GGA_X_WC                    = 118  ! Wu & Cohen
  integer, parameter :: XC_GGA_X_MPW91                 = 119  ! Modified form of PW91 by Adamo & Barone
  integer, parameter :: XC_GGA_X_AM05                  = 120  ! Armiento & Mattsson 05 exchange
  integer, parameter :: XC_GGA_X_PBEA                  = 121  ! Madsen (PBE-like)
  integer, parameter :: XC_GGA_X_MPBE                  = 122  ! Adamo & Barone modification to PBE
  integer, parameter :: XC_GGA_X_XPBE                  = 123  ! xPBE reparametrization by Xu & Goddard
  integer, parameter :: XC_GGA_X_2D_B86_MGC            = 124  ! Becke 86 MGC for 2D systems
  integer, parameter :: XC_GGA_X_BAYESIAN              = 125  ! Bayesian best fit for the enhancement factor
  integer, parameter :: XC_GGA_X_PBE_JSJR              = 126  ! JSJR reparametrization by Pedroza, Silva & Capelle
  integer, parameter :: XC_GGA_X_2D_B88                = 127  ! Becke 88 in 2D
  integer, parameter :: XC_GGA_X_2D_B86                = 128  ! Becke 86 Xalpha,beta,gamma
  integer, parameter :: XC_GGA_X_2D_PBE                = 129  ! Perdew, Burke & Ernzerhof exchange in 2D
  integer, parameter :: XC_GGA_C_PBE                   = 130  ! Perdew, Burke & Ernzerhof correlation
  integer, parameter :: XC_GGA_C_LYP                   = 131  ! Lee, Yang & Parr
  integer, parameter :: XC_GGA_C_P86                   = 132  ! Perdew 86
  integer, parameter :: XC_GGA_C_PBE_SOL               = 133  ! Perdew, Burke & Ernzerhof correlation SOL
  integer, parameter :: XC_GGA_C_PW91                  = 134  ! Perdew & Wang 91
  integer, parameter :: XC_GGA_C_AM05                  = 135  ! Armiento & Mattsson 05 correlation
  integer, parameter :: XC_GGA_C_XPBE                  = 136  ! xPBE reparametrization by Xu & Goddard
  integer, parameter :: XC_GGA_C_LM                    = 137  ! Langreth and Mehl correlation
  integer, parameter :: XC_GGA_C_PBE_JRGX              = 138  ! JRGX reparametrization by Pedroza, Silva & Capelle
  integer, parameter :: XC_GGA_X_OPTB88_VDW            = 139  ! Becke 88 reoptimized to be used with vdW functional of Dion et al
  integer, parameter :: XC_GGA_X_PBEK1_VDW             = 140  ! PBE reparametrization for vdW
  integer, parameter :: XC_GGA_X_OPTPBE_VDW            = 141  ! PBE reparametrization for vdW
  integer, parameter :: XC_GGA_X_RGE2                  = 142  ! Regularized PBE
  integer, parameter :: XC_GGA_C_RGE2                  = 143  ! Regularized PBE
  integer, parameter :: XC_GGA_X_RPW86                 = 144  ! refitted Perdew & Wang 86
  integer, parameter :: XC_GGA_X_KT1                   = 145  ! Keal and Tozer version 1
  integer, parameter :: XC_GGA_XC_KT2                  = 146  ! Keal and Tozer version 2
  integer, parameter :: XC_GGA_C_WL                    = 147  ! Wilson & Levy
  integer, parameter :: XC_GGA_C_WI                    = 148  ! Wilson & Ivanov
  integer, parameter :: XC_GGA_X_MB88                  = 149  ! Modified Becke 88 for proton transfer
  integer, parameter :: XC_GGA_X_SOGGA                 = 150  ! Second-order generalized gradient approximation
  integer, parameter :: XC_GGA_X_SOGGA11               = 151  ! Second-order generalized gradient approximation 2011
  integer, parameter :: XC_GGA_C_SOGGA11               = 152  ! Second-order generalized gradient approximation 2011
  integer, parameter :: XC_GGA_C_WI0                   = 153  ! Wilson & Ivanov initial version
  integer, parameter :: XC_GGA_XC_TH1                  = 154  ! Tozer and Handy v. 1
  integer, parameter :: XC_GGA_XC_TH2                  = 155  ! Tozer and Handy v. 2
  integer, parameter :: XC_GGA_XC_TH3                  = 156  ! Tozer and Handy v. 3
  integer, parameter :: XC_GGA_XC_TH4                  = 157  ! Tozer and Handy v. 4
  integer, parameter :: XC_GGA_X_C09X                  = 158  ! C09x to be used with the VdW of Rutgers-Chalmers
  integer, parameter :: XC_GGA_C_SOGGA11_X             = 159  ! To be used with HYB_GGA_X_SOGGA11_X
  integer, parameter :: XC_GGA_X_LB                    = 160  ! van Leeuwen & Baerends
  integer, parameter :: XC_GGA_XC_HCTH_93              = 161  ! HCTH functional fitted to  93 molecules
  integer, parameter :: XC_GGA_XC_HCTH_120             = 162  ! HCTH functional fitted to 120 molecules
  integer, parameter :: XC_GGA_XC_HCTH_147             = 163  ! HCTH functional fitted to 147 molecules
  integer, parameter :: XC_GGA_XC_HCTH_407             = 164  ! HCTH functional fitted to 407 molecules
  integer, parameter :: XC_GGA_XC_EDF1                 = 165  ! Empirical functionals from Adamson, Gill, and Pople
  integer, parameter :: XC_GGA_XC_XLYP                 = 166  ! XLYP functional
  integer, parameter :: XC_GGA_XC_B97_D                = 170  ! Grimme functional to be used with C6 vdW term
  integer, parameter :: XC_GGA_XC_PBE1W                = 173  ! Functionals fitted for water
  integer, parameter :: XC_GGA_XC_MPWLYP1W             = 174  ! Functionals fitted for water
  integer, parameter :: XC_GGA_XC_PBELYP1W             = 175  ! Functionals fitted for water
  integer, parameter :: XC_GGA_X_LBM                   = 182  ! van Leeuwen & Baerends modified
  integer, parameter :: XC_GGA_X_OL2                   = 183  ! Exchange form based on Ou-Yang and Levy v.2
  integer, parameter :: XC_GGA_X_APBE                  = 184  ! mu fixed from the semiclassical neutral atom
  integer, parameter :: XC_GGA_K_APBE                  = 185  ! mu fixed from the semiclassical neutral atom
  integer, parameter :: XC_GGA_C_APBE                  = 186  ! mu fixed from the semiclassical neutral atom
  integer, parameter :: XC_GGA_K_TW1                   = 187  ! Tran and Wesolowski set 1 (Table II)
  integer, parameter :: XC_GGA_K_TW2                   = 188  ! Tran and Wesolowski set 2 (Table II)
  integer, parameter :: XC_GGA_K_TW3                   = 189  ! Tran and Wesolowski set 3 (Table II)
  integer, parameter :: XC_GGA_K_TW4                   = 190  ! Tran and Wesolowski set 4 (Table II)
  integer, parameter :: XC_GGA_X_HTBS                  = 191  ! Haas, Tran, Blaha, and Schwarz
  integer, parameter :: XC_GGA_X_AIRY                  = 192  ! Constantin et al based on the Airy gas
  integer, parameter :: XC_GGA_X_LAG                   = 193  ! Local Airy Gas
  integer, parameter :: XC_GGA_XC_MOHLYP               = 194  ! Functional for organometallic chemistry
  integer, parameter :: XC_GGA_XC_MOHLYP2              = 195  ! Functional for barrier heights
  integer, parameter :: XC_GGA_XC_TH_FL                = 196  ! Tozer and Handy v. FL
  integer, parameter :: XC_GGA_XC_TH_FC                = 197  ! Tozer and Handy v. FC
  integer, parameter :: XC_GGA_XC_TH_FCFO              = 198  ! Tozer and Handy v. FCFO
  integer, parameter :: XC_GGA_XC_TH_FCO               = 199  ! Tozer and Handy v. FCO
  integer, parameter :: XC_GGA_C_OPTC                  = 200  ! Optimized correlation functional of Cohen and Handy
  integer, parameter :: XC_GGA_C_PBELOC                = 246  ! Semilocal dynamical correlation
  integer, parameter :: XC_GGA_XC_VV10                 = 255  ! Vydrov and Van Voorhis
  integer, parameter :: XC_GGA_C_PBEFE                 = 258  ! PBE for formation energies
  integer, parameter :: XC_GGA_C_OP_PW91               = 262  ! one-parameter progressive functional (PW91 version)
  integer, parameter :: XC_GGA_X_PBEFE                 = 265  ! PBE for formation energies
  integer, parameter :: XC_GGA_X_CAP                   = 270  ! Correct Asymptotic Potential
  integer, parameter :: XC_GGA_K_VW                    = 500  ! von Weiszaecker functional
  integer, parameter :: XC_GGA_K_GE2                   = 501  ! Second-order gradient expansion (l = 1/9)
  integer, parameter :: XC_GGA_K_GOLDEN                = 502  ! TF-lambda-vW form by Golden (l = 13/45)
  integer, parameter :: XC_GGA_K_YT65                  = 503  ! TF-lambda-vW form by Yonei and Tomishima (l = 1/5)
  integer, parameter :: XC_GGA_K_BALTIN                = 504  ! TF-lambda-vW form by Baltin (l = 5/9)
  integer, parameter :: XC_GGA_K_LIEB                  = 505  ! TF-lambda-vW form by Lieb (l = 0.185909191)
  integer, parameter :: XC_GGA_K_ABSP1                 = 506  ! gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]
  integer, parameter :: XC_GGA_K_ABSP2                 = 507  ! gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]
  integer, parameter :: XC_GGA_K_GR                    = 508  ! gamma-TFvW form by Gazquez and Robles
  integer, parameter :: XC_GGA_K_LUDENA                = 509  ! gamma-TFvW form by Ludena
  integer, parameter :: XC_GGA_K_GP85                  = 510  ! gamma-TFvW form by Ghosh and Parr
  integer, parameter :: XC_GGA_K_PEARSON               = 511  ! Pearson
  integer, parameter :: XC_GGA_K_OL1                   = 512  ! Ou-Yang and Levy v.1
  integer, parameter :: XC_GGA_K_OL2                   = 513  ! Ou-Yang and Levy v.2
  integer, parameter :: XC_GGA_K_FR_B88                = 514  ! Fuentealba & Reyes (B88 version)
  integer, parameter :: XC_GGA_K_FR_PW86               = 515  ! Fuentealba & Reyes (PW86 version)
  integer, parameter :: XC_GGA_K_DK                    = 516  ! DePristo and Kress
  integer, parameter :: XC_GGA_K_PERDEW                = 517  ! Perdew
  integer, parameter :: XC_GGA_K_VSK                   = 518  ! Vitos, Skriver, and Kollar
  integer, parameter :: XC_GGA_K_VJKS                  = 519  ! Vitos, Johansson, Kollar, and Skriver
  integer, parameter :: XC_GGA_K_ERNZERHOF             = 520  ! Ernzerhof
  integer, parameter :: XC_GGA_K_LC94                  = 521  ! Lembarki & Chermette
  integer, parameter :: XC_GGA_K_LLP                   = 522  ! Lee, Lee & Parr
  integer, parameter :: XC_GGA_K_THAKKAR               = 523  ! Thakkar 1992
  integer, parameter :: XC_GGA_X_WPBEH                 = 524  ! short-range version of the PBE
  integer, parameter :: XC_GGA_X_HJS_PBE               = 525  ! HJS screened exchange PBE version
  integer, parameter :: XC_GGA_X_HJS_PBE_SOL           = 526  ! HJS screened exchange PBE_SOL version
  integer, parameter :: XC_GGA_X_HJS_B88               = 527  ! HJS screened exchange B88 version
  integer, parameter :: XC_GGA_X_HJS_B97X              = 528  ! HJS screened exchange B97x version
  integer, parameter :: XC_GGA_X_ITYH                  = 529  ! short-range recipe for exchange GGA functionals
  integer, parameter :: XC_GGA_X_SFAT                  = 530  ! short-range recipe for exchange GGA functionals
  integer, parameter :: XC_HYB_GGA_X_N12_SX            =  81  ! N12-SX functional from Minnesota
  integer, parameter :: XC_HYB_GGA_XC_B97_1p           = 266  ! version of B97 by Cohen and Handy
  integer, parameter :: XC_HYB_GGA_XC_B3PW91           = 401  ! The original (ACM) hybrid of Becke
  integer, parameter :: XC_HYB_GGA_XC_B3LYP            = 402  ! The (in)famous B3LYP
  integer, parameter :: XC_HYB_GGA_XC_B3P86            = 403  ! Perdew 86 hybrid similar to B3PW91
  integer, parameter :: XC_HYB_GGA_XC_O3LYP            = 404  ! hybrid using the optx functional
  integer, parameter :: XC_HYB_GGA_XC_mPW1K            = 405  ! mixture of mPW91 and PW91 optimized for kinetics
  integer, parameter :: XC_HYB_GGA_XC_PBEH             = 406  ! aka PBE0 or PBE1PBE
  integer, parameter :: XC_HYB_GGA_XC_B97              = 407  ! Becke 97
  integer, parameter :: XC_HYB_GGA_XC_B97_1            = 408  ! Becke 97-1
  integer, parameter :: XC_HYB_GGA_XC_B97_2            = 410  ! Becke 97-2
  integer, parameter :: XC_HYB_GGA_XC_X3LYP            = 411  ! hybrid by Xu and Goddard
  integer, parameter :: XC_HYB_GGA_XC_B1WC             = 412  ! Becke 1-parameter mixture of WC and PBE
  integer, parameter :: XC_HYB_GGA_XC_B97_K            = 413  ! Boese-Martin for Kinetics
  integer, parameter :: XC_HYB_GGA_XC_B97_3            = 414  ! Becke 97-3
  integer, parameter :: XC_HYB_GGA_XC_MPW3PW           = 415  ! mixture with the mPW functional
  integer, parameter :: XC_HYB_GGA_XC_B1LYP            = 416  ! Becke 1-parameter mixture of B88 and LYP
  integer, parameter :: XC_HYB_GGA_XC_B1PW91           = 417  ! Becke 1-parameter mixture of B88 and PW91
  integer, parameter :: XC_HYB_GGA_XC_mPW1PW           = 418  ! Becke 1-parameter mixture of mPW91 and PW91
  integer, parameter :: XC_HYB_GGA_XC_MPW3LYP          = 419  ! mixture of mPW and LYP
  integer, parameter :: XC_HYB_GGA_XC_SB98_1a          = 420  ! Schmider-Becke 98 parameterization 1a
  integer, parameter :: XC_HYB_GGA_XC_SB98_1b          = 421  ! Schmider-Becke 98 parameterization 1b
  integer, parameter :: XC_HYB_GGA_XC_SB98_1c          = 422  ! Schmider-Becke 98 parameterization 1c
  integer, parameter :: XC_HYB_GGA_XC_SB98_2a          = 423  ! Schmider-Becke 98 parameterization 2a
  integer, parameter :: XC_HYB_GGA_XC_SB98_2b          = 424  ! Schmider-Becke 98 parameterization 2b
  integer, parameter :: XC_HYB_GGA_XC_SB98_2c          = 425  ! Schmider-Becke 98 parameterization 2c
  integer, parameter :: XC_HYB_GGA_X_SOGGA11_X         = 426  ! Hybrid based on SOGGA11 form
  integer, parameter :: XC_HYB_GGA_XC_HSE03            = 427  ! the 2003 version of the screened hybrid HSE
  integer, parameter :: XC_HYB_GGA_XC_HSE06            = 428  ! the 2006 version of the screened hybrid HSE
  integer, parameter :: XC_HYB_GGA_XC_HJS_PBE          = 429  ! HJS hybrid screened exchange PBE version
  integer, parameter :: XC_HYB_GGA_XC_HJS_PBE_SOL      = 430  ! HJS hybrid screened exchange PBE_SOL version
  integer, parameter :: XC_HYB_GGA_XC_HJS_B88          = 431  ! HJS hybrid screened exchange B88 version
  integer, parameter :: XC_HYB_GGA_XC_HJS_B97X         = 432  ! HJS hybrid screened exchange B97x version
  integer, parameter :: XC_HYB_GGA_XC_CAM_B3LYP        = 433  ! CAM version of B3LYP
  integer, parameter :: XC_HYB_GGA_XC_TUNED_CAM_B3LYP  = 434  ! CAM version of B3LYP tuned for excitations
  integer, parameter :: XC_HYB_GGA_XC_BHANDH           = 435  ! Becke half-and-half
  integer, parameter :: XC_HYB_GGA_XC_BHANDHLYP        = 436  ! Becke half-and-half with B88 exchange
  integer, parameter :: XC_HYB_GGA_XC_MB3LYP_RC04      = 437  ! B3LYP with RC04 LDA
  integer, parameter :: XC_HYB_GGA_XC_MPWLYP1M         = 453  ! MPW with 1 par. for metals/LYP
  integer, parameter :: XC_HYB_GGA_XC_REVB3LYP         = 454  ! Revised B3LYP
  integer, parameter :: XC_HYB_GGA_XC_CAMY_BLYP        = 455  ! BLYP with yukawa screening
  integer, parameter :: XC_HYB_GGA_XC_PBE0_13          = 456  ! PBE0-1/3
  integer, parameter :: XC_HYB_GGA_XC_B3LYPs           = 459  ! B3LYP* functional
  integer, parameter :: XC_HYB_GGA_XC_WB97             = 463  ! Chai and Head-Gordon
  integer, parameter :: XC_HYB_GGA_XC_WB97X            = 464  ! Chai and Head-Gordon
  integer, parameter :: XC_HYB_GGA_XC_LRC_WPBEH        = 465  ! Long-range corrected functional by Rorhdanz et al
  integer, parameter :: XC_HYB_GGA_XC_WB97X_V          = 466  ! Mardirossian and Head-Gordon
  integer, parameter :: XC_HYB_GGA_XC_LCY_PBE          = 467  ! PBE with yukawa screening
  integer, parameter :: XC_HYB_GGA_XC_LCY_BLYP         = 468  ! BLYP with yukawa screening
  integer, parameter :: XC_HYB_GGA_XC_LC_VV10          = 469  ! Vydrov and Van Voorhis
  integer, parameter :: XC_HYB_GGA_XC_CAMY_B3LYP       = 470  ! B3LYP with Yukawa screening
  integer, parameter :: XC_HYB_GGA_XC_WB97X_D          = 471  ! Chai and Head-Gordon
  integer, parameter :: XC_HYB_GGA_XC_HPBEINT          = 472  ! hPBEint
  integer, parameter :: XC_HYB_GGA_XC_LRC_WPBE         = 473  ! Long-range corrected functional by Rorhdanz et al
  integer, parameter :: XC_HYB_GGA_XC_B3LYP5           = 475  ! B3LYP with VWN functional 5 instead of RPA
  integer, parameter :: XC_HYB_GGA_XC_EDF2             = 476  ! Empirical functional from Lin, George and Gill
  integer, parameter :: XC_HYB_GGA_XC_CAP0             = 477  ! Correct Asymptotic Potential hybrid
  integer, parameter :: XC_MGGA_C_DLDF                 =  37  ! Dispersionless Density Functional
  integer, parameter :: XC_MGGA_XC_ZLP                 =  42  ! Zhao, Levy & Parr, Eq. (21)
  integer, parameter :: XC_MGGA_XC_OTPSS_D             =  64  ! oTPSS_D functional of Goerigk and Grimme
  integer, parameter :: XC_MGGA_C_CS                   =  72  ! Colle and Salvetti
  integer, parameter :: XC_MGGA_C_MN12_SX              =  73  ! Worker for MN12-SX functional
  integer, parameter :: XC_MGGA_C_MN12_L               =  74  ! MN12-L functional from Minnesota
  integer, parameter :: XC_MGGA_C_M11_L                =  75  ! M11-L functional from Minnesota
  integer, parameter :: XC_MGGA_C_M11                  =  76  ! Worker for M11 functional
  integer, parameter :: XC_MGGA_C_M08_SO               =  77  ! Worker for M08-SO functional
  integer, parameter :: XC_MGGA_C_M08_HX               =  78  ! Worker for M08-HX functional
  integer, parameter :: XC_MGGA_X_LTA                  = 201  ! Local tau approximation of Ernzerhof & Scuseria
  integer, parameter :: XC_MGGA_X_TPSS                 = 202  ! Perdew, Tao, Staroverov & Scuseria exchange
  integer, parameter :: XC_MGGA_X_M06_L                = 203  ! M06-Local functional of Minnesota
  integer, parameter :: XC_MGGA_X_GVT4                 = 204  ! GVT4 from Van Voorhis and Scuseria
  integer, parameter :: XC_MGGA_X_TAU_HCTH             = 205  ! tau-HCTH from Boese and Handy
  integer, parameter :: XC_MGGA_X_BR89                 = 206  ! Becke-Roussel 89
  integer, parameter :: XC_MGGA_X_BJ06                 = 207  ! Becke & Johnson correction to Becke-Roussel 89
  integer, parameter :: XC_MGGA_X_TB09                 = 208  ! Tran & Blaha correction to Becke & Johnson
  integer, parameter :: XC_MGGA_X_RPP09                = 209  ! Rasanen, Pittalis, and Proetto correction to Becke & Johnson
  integer, parameter :: XC_MGGA_X_2D_PRHG07            = 210  ! Pittalis, Rasanen, Helbig, Gross Exchange Functional
  integer, parameter :: XC_MGGA_X_2D_PRHG07_PRP10      = 211  ! PRGH07 with PRP10 correction
  integer, parameter :: XC_MGGA_X_REVTPSS              = 212  ! revised Perdew, Tao, Staroverov & Scuseria exchange
  integer, parameter :: XC_MGGA_X_PKZB                 = 213  ! Perdew, Kurth, Zupan, and Blaha
  integer, parameter :: XC_MGGA_X_M05                  = 214  ! Worker for M05 functional
  integer, parameter :: XC_MGGA_X_M05_2X               = 215  ! Worker for M05-2X functional
  integer, parameter :: XC_MGGA_X_M06_HF               = 216  ! Worker for M06-HF functional
  integer, parameter :: XC_MGGA_X_M06                  = 217  ! Worker for M06 functional
  integer, parameter :: XC_MGGA_X_M06_2X               = 218  ! Worker for M06-2X functional
  integer, parameter :: XC_MGGA_X_M08_HX               = 219  ! Worker for M08-HX functional
  integer, parameter :: XC_MGGA_X_M08_SO               = 220  ! Worker for M08-SO functional
  integer, parameter :: XC_MGGA_X_MS0                  = 221  ! MS exchange of Sun, Xiao, and Ruzsinszky
  integer, parameter :: XC_MGGA_X_MS1                  = 222  ! MS1 exchange of Sun, et al
  integer, parameter :: XC_MGGA_X_MS2                  = 223  ! MS2 exchange of Sun, et al
  integer, parameter :: XC_MGGA_X_M11                  = 225  ! Worker for M11 functional
  integer, parameter :: XC_MGGA_X_M11_L                = 226  ! M11-L functional from Minnesota
  integer, parameter :: XC_MGGA_X_MN12_L               = 227  ! MN12-L functional from Minnesota
  integer, parameter :: XC_MGGA_C_CC06                 = 229  ! Cancio and Chou 2006
  integer, parameter :: XC_MGGA_X_MK00                 = 230  ! Exchange for accurate virtual orbital energies
  integer, parameter :: XC_MGGA_C_TPSS                 = 231  ! Perdew, Tao, Staroverov & Scuseria correlation
  integer, parameter :: XC_MGGA_C_VSXC                 = 232  ! VSxc from Van Voorhis and Scuseria (correlation part)
  integer, parameter :: XC_MGGA_C_M06_L                = 233  ! M06-Local functional from Minnesota
  integer, parameter :: XC_MGGA_C_M06_HF               = 234  ! Worker for M06-HF functional
  integer, parameter :: XC_MGGA_C_M06                  = 235  ! Worker for M06 functional
  integer, parameter :: XC_MGGA_C_M06_2X               = 236  ! Worker for M06-2X functional
  integer, parameter :: XC_MGGA_C_M05                  = 237  ! Worker for M05 functional
  integer, parameter :: XC_MGGA_C_M05_2X               = 238  ! Worker for M05-2X functional
  integer, parameter :: XC_MGGA_C_PKZB                 = 239  ! Perdew, Kurth, Zupan, and Blaha
  integer, parameter :: XC_MGGA_C_BC95                 = 240  ! Becke correlation 95
  integer, parameter :: XC_MGGA_C_REVTPSS              = 241  ! revised TPSS correlation
  integer, parameter :: XC_MGGA_XC_TPSSLYP1W           = 242  ! Functionals fitted for water
  integer, parameter :: XC_MGGA_X_MK00B                = 243  ! Exchange for accurate virtual orbital energies (v. B)
  integer, parameter :: XC_MGGA_X_BLOC                 = 244  ! functional with balanced localization
  integer, parameter :: XC_MGGA_X_MODTPSS              = 245  ! Modified Perdew, Tao, Staroverov & Scuseria exchange
  integer, parameter :: XC_MGGA_C_TPSSLOC              = 247  ! Semilocal dynamical correlation
  integer, parameter :: XC_MGGA_X_MBEEF                = 249  ! mBEEF exchange
  integer, parameter :: XC_MGGA_X_MBEEFVDW             = 250  ! mBEEF-vdW exchange
  integer, parameter :: XC_MGGA_XC_B97M_V              = 254  ! Mardirossian and Head-Gordon
  integer, parameter :: XC_MGGA_X_MVS                  = 257  ! MVS exchange of Sun, Perdew, and Ruzsinszky
  integer, parameter :: XC_MGGA_X_MN15_L               = 260  ! MN15-L functional from Minnesota
  integer, parameter :: XC_MGGA_C_MN15_L               = 261  ! MN15-L functional from Minnesota
  integer, parameter :: XC_MGGA_X_SCAN                 = 263  ! SCAN exchange of Sun, Ruzsinszky, and Perdew
  integer, parameter :: XC_MGGA_C_SCAN                 = 267  ! SCAN correlation
  integer, parameter :: XC_MGGA_C_MN15                 = 269  ! MN15 functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_X_DLDF             =  36  ! Dispersionless Density Functional
  integer, parameter :: XC_HYB_MGGA_X_MS2H             = 224  ! MS2 hybrid exchange of Sun, et al
  integer, parameter :: XC_HYB_MGGA_X_MN12_SX          = 248  ! MN12-SX hybrid functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_X_SCAN0            = 264  ! SCAN hybrid
  integer, parameter :: XC_HYB_MGGA_X_MN15             = 268  ! MN15 functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_M05             = 438  ! M05 functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_M05_2X          = 439  ! M05-2X functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_B88B95          = 440  ! Mixture of B88 with BC95 (B1B95)
  integer, parameter :: XC_HYB_MGGA_XC_B86B95          = 441  ! Mixture of B86 with BC95
  integer, parameter :: XC_HYB_MGGA_XC_PW86B95         = 442  ! Mixture of PW86 with BC95
  integer, parameter :: XC_HYB_MGGA_XC_BB1K            = 443  ! Mixture of B88 with BC95 from Zhao and Truhlar
  integer, parameter :: XC_HYB_MGGA_XC_M06_HF          = 444  ! M06-HF functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_MPW1B95         = 445  ! Mixture of mPW91 with BC95 from Zhao and Truhlar
  integer, parameter :: XC_HYB_MGGA_XC_MPWB1K          = 446  ! Mixture of mPW91 with BC95 for kinetics
  integer, parameter :: XC_HYB_MGGA_XC_X1B95           = 447  ! Mixture of X with BC95
  integer, parameter :: XC_HYB_MGGA_XC_XB1K            = 448  ! Mixture of X with BC95 for kinetics
  integer, parameter :: XC_HYB_MGGA_XC_M06             = 449  ! M06 functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_M06_2X          = 450  ! M06-2X functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_PW6B95          = 451  ! Mixture of PW91 with BC95 from Zhao and Truhlar
  integer, parameter :: XC_HYB_MGGA_XC_PWB6K           = 452  ! Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics
  integer, parameter :: XC_HYB_MGGA_XC_TPSSH           = 457  ! TPSS hybrid
  integer, parameter :: XC_HYB_MGGA_XC_REVTPSSH        = 458  ! revTPSS hybrid
  integer, parameter :: XC_HYB_MGGA_XC_M08_HX          = 460  ! M08-HX functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_M08_SO          = 461  ! M08-SO functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_XC_M11             = 462  ! M11    functional from Minnesota
  integer, parameter :: XC_HYB_MGGA_X_MVSH             = 474  ! MVS hybrid
  integer, parameter :: XC_HYB_MGGA_XC_WB97M_V         = 531  ! Mardirossian and Head-Gordon

end module libxc_funcs_m

module xc_f90_lib_m

  use xc_f90_types_m
  use libxc_funcs_m

  implicit none

  public

  ! Families of xc functionals
  integer, parameter :: &
    XC_FAMILY_UNKNOWN = -1, &
    XC_FAMILY_NONE = 0, &
    XC_FAMILY_LDA = 1, &
    XC_FAMILY_GGA = 2, &
    XC_FAMILY_MGGA = 4, &
    XC_FAMILY_LCA = 8, &
    XC_FAMILY_OEP = 16, &
    XC_FAMILY_HYB_GGA = 32, &
    XC_FAMILY_HYB_MGGA = 64

  ! Spin
  integer, parameter :: &
    XC_UNPOLARIZED = 1, & ! Spin unpolarized
    XC_POLARIZED = 2 ! Spin polarized

  integer, parameter :: &
    XC_NON_RELATIVISTIC = 0, & ! Functional includes or not relativistic
    XC_RELATIVISTIC = 1 ! corrections. Only available in some functionals.

  ! Kinds
  integer, parameter :: &
    XC_EXCHANGE = 0, &
    XC_CORRELATION = 1, &
    XC_EXCHANGE_CORRELATION = 2, &
    XC_KINETIC = 3

  integer, parameter :: &
    XC_FLAGS_HAVE_EXC = 1, &
    XC_FLAGS_HAVE_VXC = 2, &
    XC_FLAGS_HAVE_FXC = 4, &
    XC_FLAGS_HAVE_KXC = 8, &
    XC_FLAGS_HAVE_LXC = 16, &
    XC_FLAGS_1D = 32, &
    XC_FLAGS_2D = 64, &
    XC_FLAGS_3D = 128, &
    XC_FLAGS_STABLE = 512, &
    XC_FLAGS_DEVELOPMENT = 1024

  ! These are old names kept for compatibility, and that should disappear soon
  integer, parameter :: XC_GGA_C_VPBE = 83
  integer, parameter :: XC_GGA_XC_LB = 160
  integer, parameter :: XC_GGA_K_ABSR1 = 506
  integer, parameter :: XC_GGA_K_ABSR2 = 507

  contains

  !----------------------------------------------------------------
    subroutine xc_f90_version(major, minor, micro)
      use localorb_io
      use mpi_tasks
      integer, intent(inout) :: major, minor, micro
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_version

    subroutine xc_f90_version_string(version)
      use localorb_io
      use mpi_tasks
      character(len=10), intent(inout) :: version
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_version_string

  !----------------------------------------------------------------
    integer function xc_f90_info_number(info) result(number)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      ! Output an error
      character*128              :: info_str
      number = 0
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end function xc_f90_info_number

    integer function xc_f90_info_kind(info) result(number)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      ! Output an error
      character*128              :: info_str
      number = 0
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      character(len=*), intent(inout) :: s
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info) result(number)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      ! Output an error
      character*128              :: info_str
      number = 0
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end function xc_f90_info_family

    integer function xc_f90_info_flags(info) result(number)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      ! Output an error
      character*128              :: info_str
      number = 0
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end function xc_f90_info_flags

    subroutine xc_f90_info_refs(info, number, s)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      integer, intent(inout) :: number ! number of the reference. Must be 0 in the first call
      character(len=*), intent(inout) :: s ! the string that is output
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_info_refs

    subroutine xc_f90_functional_get_name(func_number, func_string)
      use localorb_io
      use mpi_tasks
      integer, intent(in) :: func_number
      character(len=256), intent(inout) :: func_string
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_functional_get_name

    integer function xc_f90_functional_get_number(func_string) result(number)
      use localorb_io
      use mpi_tasks
      character(len=*), intent(in) :: func_string
      ! Output an error
      character*128              :: info_str
      number = 0
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end function xc_f90_functional_get_number

    integer function xc_f90_family_from_id(id) result(number)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      integer, intent(in) :: id
      ! Output an error
      character*128              :: info_str
      number = 0
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end function xc_f90_family_from_id


  !----------------------------------------------------------------
    subroutine xc_f90_func_init(p, info, functional, nspin)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      type(xc_f90_pointer_t), intent(inout) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_func_init

    subroutine xc_f90_func_end(p)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_func_end


  ! LDAs
  !----------------------------------------------------------------
    subroutine xc_f90_lda(p, np, rho, zk, vrho, fxc, kxc)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: zk ! the energy per unit particle
      real(xc_f90_kind), intent(inout) :: vrho ! v(nspin) the potential
      real(xc_f90_kind), intent(inout) :: fxc ! v(nspin,nspin) the xc kernel
      real(xc_f90_kind), intent(inout) :: kxc ! v(nspin,nspin,nspin) the derivative of xc kernel
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda

    subroutine xc_f90_lda_exc(p, np, rho, zk)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: zk ! the energy per unit particle
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_exc

    subroutine xc_f90_lda_exc_vxc(p, np, rho, e, v)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: e ! the energy per unit particle
      real(xc_f90_kind), intent(inout) :: v ! v(nspin) the potential
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_exc_vxc

    subroutine xc_f90_lda_vxc(p, np, rho, v)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: v ! v(nspin) the potential
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_vxc

    subroutine xc_f90_lda_vxc_fxc(p, np, rho, v, fxc)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: v ! v(nspin) the potential
      real(xc_f90_kind), intent(inout) :: fxc ! v(nspin,nspin) the xc kernel
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_vxc_fxc

    subroutine xc_f90_lda_fxc(p, np, rho, fxc)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: fxc ! v(nspin,nspin) the xc kernel
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_fxc

    subroutine xc_f90_lda_kxc(p, np, rho, kxc)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(inout) :: kxc
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_kxc


    subroutine xc_f90_lda_x_1d_set_par(p, interaction, bb)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      integer, intent(in) :: interaction
      real(xc_f90_kind), intent(in) :: bb
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_x_1d_set_par

    subroutine xc_f90_lda_c_xalpha_set_par(p, alpha)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: alpha
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_c_xalpha_set_par

    subroutine xc_f90_lda_x_set_par(p, alpha, relativistic, omega)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: alpha ! of Xalpha, set to 4/3 to obtain standard LDA
      integer, intent(in) :: relativistic
      real(xc_f90_kind), intent(in) :: omega
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_x_set_par

    subroutine xc_f90_lda_c_1d_csc_set_par(p, interaction, bb)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      integer, intent(in) :: interaction
      real(xc_f90_kind), intent(in) :: bb
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_c_1d_csc_set_par

    subroutine xc_f90_lda_c_2d_prm_set_par(p, N)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: N
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_lda_c_2d_prm_set_par

  ! GGAs
  !----------------------------------------------------------------
    subroutine xc_f90_gga(p, np, rho, sigma, zk, vrho, vsigma, &
        v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: zk
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      real(xc_f90_kind), intent(inout) :: v2rho2
      real(xc_f90_kind), intent(inout) :: v2rhosigma
      real(xc_f90_kind), intent(inout) :: v2sigma2
      real(xc_f90_kind), intent(inout) :: v3rho3
      real(xc_f90_kind), intent(inout) :: v3rho2sigma
      real(xc_f90_kind), intent(inout) :: v3rhosigma2
      real(xc_f90_kind), intent(inout) :: v3sigma3
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga

    subroutine xc_f90_gga_exc(p, np, rho, sigma, zk)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: zk
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_exc

    subroutine xc_f90_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: zk
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_exc_vxc

    subroutine xc_f90_gga_vxc(p, np, rho, sigma, vrho, vsigma)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_vxc

    subroutine xc_f90_gga_vxc_fxc(p, np, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      real(xc_f90_kind), intent(inout) :: v2rho2
      real(xc_f90_kind), intent(inout) :: v2rhosigma
      real(xc_f90_kind), intent(inout) :: v2sigma2
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_vxc_fxc

    subroutine xc_f90_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: v2rho2
      real(xc_f90_kind), intent(inout) :: v2rhosigma
      real(xc_f90_kind), intent(inout) :: v2sigma2
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_fxc

    subroutine xc_f90_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(inout) :: v3rho3
      real(xc_f90_kind), intent(inout) :: v3rho2sigma
      real(xc_f90_kind), intent(inout) :: v3rhosigma2
      real(xc_f90_kind), intent(inout) :: v3sigma3
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_kxc

  !----------------------------------------------------------------
    subroutine xc_f90_gga_lb_set_par(p, modified, threshold, ip, qtot)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: modified ! should we use the modified version
      real(xc_f90_kind), intent(in) :: threshold ! if so, the threshold to use the asymptotic version
      real(xc_f90_kind), intent(in) :: ip ! ionization potential
      real(xc_f90_kind), intent(in) :: qtot ! total charge
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_lb_set_par


  !----------------------------------------------------------------
    subroutine xc_f90_gga_lb_modified(p, np, rho, grho, r, dedd)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(in) :: grho ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind), intent(in) :: r ! distance from center of finite system
      real(xc_f90_kind), intent(inout) :: dedd
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_lb_modified

  !----------------------------------------------------------------
    subroutine xc_f90_gga_x_wpbeh_set_par(p, omega)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: omega ! range separation
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_x_wpbeh_set_par

  !----------------------------------------------------------------
    subroutine xc_f90_gga_x_hjs_set_par(p, omega)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: omega ! range separation
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_x_hjs_set_par

  !----------------------------------------------------------------
    subroutine xc_f90_gga_ak13_get_asymptotic(homo, asymp)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      real(xc_f90_kind), intent(in) :: homo
      real(xc_f90_kind), intent(inout) :: asymp
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_gga_ak13_get_asymptotic

  !----------------------------------------------------------------
    subroutine xc_f90_hyb_exx_coef(p, coef)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(inout) :: coef
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_hyb_exx_coef

    subroutine xc_f90_hyb_cam_coef(p, omega, alpha, beta)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(inout) :: omega, alpha, beta
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_hyb_cam_coef


  !----------------------------------------------------------------
    subroutine xc_f90_hyb_gga_xc_hse_set_par(p, beta, omega)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: beta ! mixing
      real(xc_f90_kind), intent(in) :: omega ! range separation
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_hyb_gga_xc_hse_set_par

    subroutine xc_f90_hyb_gga_xc_pbeh_set_par(p, alpha)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: alpha ! mixing
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_hyb_gga_xc_pbeh_set_par


  ! the meta-GGAs
  !----------------------------------------------------------------
    subroutine xc_f90_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(inout) :: zk
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      real(xc_f90_kind), intent(inout) :: vlapl
      real(xc_f90_kind), intent(inout) :: vtau
      real(xc_f90_kind), intent(inout) :: v2rho2
      real(xc_f90_kind), intent(inout) :: v2sigma2
      real(xc_f90_kind), intent(inout) :: v2lapl2
      real(xc_f90_kind), intent(inout) :: v2tau2
      real(xc_f90_kind), intent(inout) :: v2rhosigma
      real(xc_f90_kind), intent(inout) :: v2rholapl
      real(xc_f90_kind), intent(inout) :: v2rhotau
      real(xc_f90_kind), intent(inout) :: v2sigmalapl
      real(xc_f90_kind), intent(inout) :: v2sigmatau
      real(xc_f90_kind), intent(inout) :: v2lapltau
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga

    subroutine xc_f90_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(inout) :: zk
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga_exc

    subroutine xc_f90_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(inout) :: zk
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      real(xc_f90_kind), intent(inout) :: vlapl
      real(xc_f90_kind), intent(inout) :: vtau
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga_exc_vxc

    subroutine xc_f90_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      real(xc_f90_kind), intent(inout) :: vlapl
      real(xc_f90_kind), intent(inout) :: vtau
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga_vxc

    subroutine xc_f90_mgga_vxc_fxc(p, np, rho, sigma, lapl, tau, &
      vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(inout) :: vrho
      real(xc_f90_kind), intent(inout) :: vsigma
      real(xc_f90_kind), intent(inout) :: vlapl
      real(xc_f90_kind), intent(inout) :: vtau
      real(xc_f90_kind), intent(inout) :: v2rho2
      real(xc_f90_kind), intent(inout) :: v2sigma2
      real(xc_f90_kind), intent(inout) :: v2lapl2
      real(xc_f90_kind), intent(inout) :: v2tau2
      real(xc_f90_kind), intent(inout) :: v2rhosigma
      real(xc_f90_kind), intent(inout) :: v2rholapl
      real(xc_f90_kind), intent(inout) :: v2rhotau
      real(xc_f90_kind), intent(inout) :: v2sigmalapl
      real(xc_f90_kind), intent(inout) :: v2sigmatau
      real(xc_f90_kind), intent(inout) :: v2lapltau
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga_vxc_fxc

    subroutine xc_f90_mgga_fxc(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(inout) :: v2rho2
      real(xc_f90_kind), intent(inout) :: v2sigma2
      real(xc_f90_kind), intent(inout) :: v2lapl2
      real(xc_f90_kind), intent(inout) :: v2tau2
      real(xc_f90_kind), intent(inout) :: v2rhosigma
      real(xc_f90_kind), intent(inout) :: v2rholapl
      real(xc_f90_kind), intent(inout) :: v2rhotau
      real(xc_f90_kind), intent(inout) :: v2sigmalapl
      real(xc_f90_kind), intent(inout) :: v2sigmatau
      real(xc_f90_kind), intent(inout) :: v2lapltau
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga_fxc

    subroutine xc_f90_mgga_x_tb09_set_par(p, cc)
      use localorb_io
      use mpi_tasks
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: cc
      ! Output an error
      character*128              :: info_str
      write(info_str,'(a)') '***** This is a call to libxc_stubs. This means you have selected'
      call localorb_info(info_str)
      write(info_str,'(a)') '***** functionals from libxc but DID NOT link against it.'
      call localorb_info(info_str)
      call aims_stop(' General STOP. You have to compile with libxc to get functionals from it!')
    end subroutine xc_f90_mgga_x_tb09_set_par

end module xc_f90_lib_m
