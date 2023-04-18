module xc_library

implicit none

private

public :: get_xc_name, translate_to_sr, translate_from_sr, &
    is_xc_lda, is_xc_gga, is_xc_meta_gga, is_xc_hybrid, is_xc_nonlocal, &
    is_xc_double_hybrid, normalize_flag_xc, is_dfauto

integer, parameter, public :: &
    xc__none = -1, &
    xc__hf = 0, &
    xc__pbe0 = 1, &
    xc__pz_lda = 3, &
    xc__pw91 = 4, &
    xc__pbe_vdwdf = 5, &
    xc__pbe = 6, &
    xc__hse = 7, &
    xc__pw_lda = 8, &
    xc__blyp = 9, &
    xc__b3lyp = 10, &
    xc__rpbe = 11, &
    xc__revpbe = 12, &
    xc__pbesol0 = 13, &
    xc__revpbe_vdwdf = 14, &
    xc__vwn = 15, &  ! this is "vwn5" (the full LDA parametrization by vwn)
    xc__vwn_gauss = 16, &  ! this is the RPA paramerization by vwn
    xc__pbesol = 17, &
    xc__am05 = 20, &
    xc__pbeint = 21, &
    xc__lc_wpbeh = 23, &
    xc__m06_l = 25, &
    xc__m06_hf = 26, &
    xc__m06 = 27, &
    xc__m06_2x = 28, &
    xc__b1lyp = 29, &
    xc__xyg3 = 30, &
    xc__xdh_pbe0 = 31, &
    xc__xygj_os = 32, &
    xc__zrps = 33, &   ! what is this?
    xc__tpss = 51, &
    xc__revtpss = 52, &
    xc__tpssloc = 53, &
    xc__m08_hx = 55, &
    xc__m08_so = 56, &
    xc__m11 = 57, &
    xc__m11_l = 58, &
    xc__scan = 59, &
    xc__scan0 = 60 !, &
!    xc__libxc = 100

integer, parameter, public :: &
    xc_sr__pw_lda = 8, &
    xc_sr__pbe = 6, &
    xc_sr__vwn = 7, &
    xc_sr__blyp = 9, &
    xc_sr__rpbe = 14, &
    xc_sr__revpbe = 15

integer, parameter, public :: &
    xc_flag_max = 100, &
    xc_dfauto_offset = 1000

contains

logical function is_dfauto(flag_xc)
    integer, intent(in) :: flag_xc

    is_dfauto = flag_xc > xc_dfauto_offset .and. flag_xc < xc_dfauto_offset+xc_flag_max
end function

! Subtract potential offsets from flag_xc
integer function normalize_flag_xc(flag_xc) result(normed)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case (xc_dfauto_offset:xc_dfauto_offset+xc_flag_max)
        normed = flag_xc-xc_dfauto_offset
    case default
        normed = flag_xc
    end select
end function

logical function is_xc_lda(flag_xc)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case (xc__pw_lda, xc__pz_lda, xc__vwn, xc__vwn_gauss)
        is_xc_lda = .true.
    case default
        is_xc_lda = .false.
    end select
end function

logical function is_xc_gga(flag_xc)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case ( &
            xc__am05, xc__blyp, xc__pbe, xc__pbeint, xc__pbesol, xc__rpbe, &
            xc__revpbe, xc__pw91, xc__pbe_vdwdf, xc__revpbe_vdwdf &
    )
        is_xc_gga = .true.
    case default
        is_xc_gga = .false.
    end select
end function

logical function is_xc_meta_gga(flag_xc)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case ( &  ! non-hybrid meta-GGAs
            xc__m06_l, xc__m11_l, xc__revtpss, xc__tpss, xc__tpssloc, &
            xc__scan &
    )
        is_xc_meta_gga = .true.
    case ( &  ! hybrid meta-GGAs
            xc__m06, xc__m06_2x, xc__m06_hf, xc__m08_hx, xc__m08_so, &
            xc__m11, xc__scan0 &
    )
        is_xc_meta_gga = .true.
    case default
        is_xc_meta_gga = .false.
    end select
end function

logical function is_xc_hybrid(flag_xc)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case ( &  ! hybrid GGAs
            xc__hf, xc__pbe0, xc__hse, xc__b3lyp, xc__pbesol0, &
            xc__lc_wpbeh, xc__b1lyp &
    )
        is_xc_hybrid = .true.
    case ( &  ! hybrid meta-GGAs
            xc__m06, xc__m06_2x, xc__m06_hf, xc__m08_hx, xc__m08_so, &
            xc__m11, xc__scan0 &
    )
        is_xc_hybrid = .true.
    case default
        is_xc_hybrid = .false.
    end select
end function

logical function is_xc_double_hybrid(flag_xc)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case (xc__xdh_pbe0, xc__xyg3, xc__xygj_os, xc__zrps)
        is_xc_double_hybrid = .true.
    case default
        is_xc_double_hybrid = .false.
    end select
end function

logical function is_xc_nonlocal(flag_xc)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case (xc__pbe_vdwdf, xc__revpbe_vdwdf)
        is_xc_nonlocal = .true.
    case default
        is_xc_nonlocal = .false.
    end select
end function

integer function translate_to_sr(flag_xc) result(xc_sr_flag)
    integer, intent(in) :: flag_xc

    select case (flag_xc)
    case (xc__pw_lda); xc_sr_flag = xc_sr__pw_lda
    case (xc__pbe); xc_sr_flag = xc_sr__pbe
    case (xc__vwn); xc_sr_flag = xc_sr__vwn
    case (xc__blyp); xc_sr_flag = xc_sr__blyp
    case (xc__rpbe); xc_sr_flag = xc_sr__rpbe
    case (xc__revpbe); xc_sr_flag = xc_sr__revpbe
    case default; xc_sr_flag = -1
    end select
end function

integer function translate_from_sr(xc_sr_flag) result(flag_xc)
    integer, intent(in) :: xc_sr_flag

    select case (xc_sr_flag)
    case (xc_sr__pw_lda); flag_xc = xc__pw_lda
    case (xc_sr__pbe); flag_xc = xc__pbe
    case (xc_sr__vwn); flag_xc = xc__vwn
    case (xc_sr__blyp); flag_xc = xc__blyp
    case (xc_sr__rpbe); flag_xc = xc__rpbe
    case (xc_sr__revpbe); flag_xc = xc__revpbe
    case default; flag_xc = -1
    end select
end function

 function get_xc_name(flag_xc) result(xc_name)
    integer, intent(in) :: flag_xc
    character(len=:), allocatable :: xc_name

    select case (flag_xc)
    case (xc__hf); xc_name = 'HF'
    case (xc__pbe0); xc_name = 'PBE0'
    case (xc__pz_lda); xc_name = 'PZ-LDA'
    case (xc__pbe_vdwdf); xc_name = 'PBE+vdW-DF'
    case (xc__pbe); xc_name = 'PBE'
    case (xc__hse); xc_name = 'HSE'
    case (xc__pw_lda); xc_name = 'PW-LDA'
    case (xc__blyp); xc_name = 'BLYP'
    case (xc__b3lyp); xc_name = 'B3LYP'
    case (xc__rpbe); xc_name = 'RPBE'
    case (xc__revpbe); xc_name = 'revPBE'
    case (xc__pbesol0); xc_name = 'PBEsol'
    case (xc__b1lyp); xc_name = 'B1LYP'
    case (xc__revpbe_vdwdf); xc_name = 'revPBE+vdW-DF'
    case (xc__vwn); xc_name = 'VWN'
    case (xc__vwn_gauss); xc_name = 'VNW(RPA)'
    case (xc__pbesol); xc_name = 'PBEsol'
    case (xc__am05); xc_name = 'AM05'
    case (xc__pbeint); xc_name = 'PBEint'
    case (xc__lc_wpbeh); xc_name = 'LC-wPBEh'
    case (xc__m06_l); xc_name = 'M06-L'
    case (xc__m06_hf); xc_name = 'M06-HF'
    case (xc__m06); xc_name = 'M06'
    case (xc__m06_2x); xc_name = 'M06-2X'
    case (xc__xyg3); xc_name = 'XYG3'
    case (xc__xdh_pbe0); xc_name = 'xDH-PBE0'
    case (xc__xygj_os); xc_name = 'XYGJ-OS'
    case (xc__zrps); xc_name = '?ZRPS'
    case (xc__tpss); xc_name = 'TPSS'
    case (xc__revtpss); xc_name = 'revTPSS'
    case (xc__tpssloc); xc_name = 'TPSSloc'
    case (xc__m08_hx); xc_name = 'M08-HX'
    case (xc__m08_so); xc_name = 'M08-SO'
    case (xc__m11); xc_name = 'M11'
    case (xc__m11_l); xc_name = 'M11-L'
    case (xc__scan); xc_name = 'SCAN'
    case (xc__scan0); xc_name = 'SCAN0'
!    case (xc__libxc); xc_name = 'Libxc'
    case default; xc_name = ''
    end select
end function

end module
