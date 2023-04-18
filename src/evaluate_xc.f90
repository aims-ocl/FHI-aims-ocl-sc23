!***** FHI-aims/evaluate_xc
! NAME
!   evaluate_xc
! SYNOPSIS
subroutine evaluate_xc( &
        rho, rho_gradient, kinetic_density, &
        en_density_xc, en_density_x, en_density_c, &
        local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, &
        initial, coord_current)
    ! PURPOSE
    !   Subroutine evaluate_xc evaluates the exchange correlation potential
    !   and energy contribution for a given density at one integration point.
    !
    !   VB: The present version supports LDA and gradient functionals.
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! SEE ALSO
    !   Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !   Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !   "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !   Computer Physics Communications (2008), submitted.
    ! COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    ! HISTORY
    !   Release version, FHI-aims (2008).
    ! USES
    use runtime_choices
    use dimensions
    use xc
    use xc_library
    use constants
    use mpi_tasks, only: stderr

    implicit none

    ! INPUTS
    !   o rho : Density at current point, given as spin up and spin down if spin-polarized
    !   o rho_gradient : Density gradient at current point
    !   o kinetic_density : Kinetic density at current point (use_meta_gga)
    !   o initial : Is initial evaluation? (no EXX)
    !   o coord_current : ?
    real*8, intent(in) :: rho(n_spin)
    real*8, intent(in) :: rho_gradient(3, n_spin)
    real*8, intent(in) :: kinetic_density(n_spin)
    logical, intent(in) :: initial
    real*8, intent(in), optional :: coord_current(3)

    ! OUTPUT
    !   o en_density_xc : Exchange-correlation energy density at current point
    !   o en_density_x, en_density_y : Components of en_density_xc
    !   o local_xc_derivs : Local parts of the XC potential:
    !       Partial derivative of the exchange-correlation
    !       energy functional by the density plus spin-weighted
    !       partial derivative of the XC energy functional by spin
    !   o xc_gradient_deriv : Partial derivative of the exchange-correlation
    !       energy functional by the modulus square of the density gradient
    !   o xc_tau_deriv : Partial derivative of the exchange-correlation energy
    !       functional by tau, the kinetic density
    real*8, intent(out) :: en_density_xc
    real*8, intent(out) :: en_density_x(n_spin)
    real*8, intent(out) :: en_density_c
    real*8, intent(out) :: local_xc_derivs(n_spin)
    real*8, intent(out) :: xc_gradient_deriv(3, n_spin)
    real*8, intent(out) :: xc_tau_deriv(n_spin)

    ! SOURCE
    logical :: xc_undefined
    logical :: initial_was_evaluated
    ! counters
    integer :: i_spin
    ! LibXC
    integer :: flag_libxc_x_id
    integer :: flag_libxc_c_id
    ! for case(18) R48PBE functional mix
    real*8 :: en_density_xc2
    real*8 :: en_density_x2(n_spin)
    real*8 :: en_density_c2
    real*8 :: local_xc_derivs2(n_spin)
    real*8 :: xc_gradient_deriv2(3, n_spin)
    ! temporary variables
    integer :: temp

    flag_libxc_x_id = 0
    flag_libxc_c_id = 0

    en_density_xc = 0.d0
    en_density_x = 0.d0
    en_density_c = 0.d0
    local_xc_derivs = 0.d0
    xc_gradient_deriv = 0.d0
    xc_tau_deriv = 0.d0

    ! (Rundong) The fully-rel(X2C and 4C-DKS) code only treat closed shell 
    ! currently. Due to the existing code structure, I have to set n_spin=1 
    ! here (although it could be 2). And in the end of this subroutine, we set
    ! it back to the initial value.
    temp = n_spin
    if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) n_spin=1

    xc_undefined = .true.
    do i_spin = 1, n_spin, 1
        ! This is true if both densities are below zero
        xc_undefined = xc_undefined .and. (rho(i_spin) .le. 0.d0)
    enddo
    do i_spin = 1, n_spin, 1
        ! This is true if one of the densities is less than zero
        xc_undefined = xc_undefined .or. (rho(i_spin) .lt. 0.d0)
    enddo
    if (xc_undefined)then
      if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) n_spin=temp
      return
    endif

    ! When initial == .true., exact exchange is not yet evaluated and hybrid
    ! functionals must be initialized via substitutes. This is handled here. The
    ! substitutes are indicated in comments.
    if (initial) then
        initial_was_evaluated = .true.
        ! This will be set to false in "case default" below if the functional is
        ! not evaluated here. If so, the second select case construct is
        ! entered; otherwise, an early exit from this routine is followed.

        ! AJL, Feb2017: Jan - I picked up a bug for dfauto below when
        ! initialising a meta-GGA calculation. The fix is horrid and you
        ! might want to look and see if it can be tidied up.
        select case (flag_xc)
        case (0)  ! HF -> LDA
            call xcpot_pw91_lda(rho, local_xc_derivs, en_density_xc, &
                en_density_x, en_density_c)
        case (1)  ! PBE0 -> PBE
            call xc_partials_pbe( &
                rho, rho_gradient, en_density_xc, en_density_x, &
                en_density_c, local_xc_derivs, xc_gradient_deriv)
        case (7)  ! HSE -> PBE
            call xc_partials_pbe( &
                rho, rho_gradient, en_density_xc, en_density_x, &
                en_density_c, local_xc_derivs, xc_gradient_deriv)
        case (10)  ! B3LYP -> BLYP
            call xc_partials_blyp( &
                rho, rho_gradient, en_density_xc, en_density_x, &
                en_density_c, local_xc_derivs, xc_gradient_deriv)
        case (13)  ! PBEsol0 -> PBEsol
            call xc_partials_pbesol( &
                rho, rho_gradient, en_density_xc, en_density_x, &
                en_density_c, local_xc_derivs, xc_gradient_deriv)
        case (23)  ! LC-wPBEh -> PBE
            call xc_partials_pbe( &
                rho, rho_gradient, en_density_xc, en_density_x, &
                en_density_c, local_xc_derivs, xc_gradient_deriv)
        case (29)  ! B1LYP -> BLYP
            call xc_partials_blyp( &
                rho, rho_gradient, en_density_xc, en_density_x, &
                en_density_c, local_xc_derivs, xc_gradient_deriv)
        case (25:28, 51:53, 55:58, 59, :-10)
            ! M06, M08, M11, TPSS, SCAN, LibXC -> LDA
            call xcpot_pw91_lda( &
                rho, local_xc_derivs, en_density_xc, en_density_x, en_density_c)
        case default
            ! AJL: Added a horrible catch for dfauto meta-GGA calculations - this
            ! needs to be made tidier somehow but works for now
            if ((is_dfauto(flag_xc)) .and. &
                (normalize_flag_xc(flag_xc).eq.xc__tpss.or. &
                 normalize_flag_xc(flag_xc).eq.xc__scan.or. &
                 normalize_flag_xc(flag_xc).eq.xc__scan0)) then
                    call xcpot_pw91_lda( &
                    rho, local_xc_derivs, en_density_xc, en_density_x, en_density_c)
                    initial_was_evaluated = .true.
            else
                initial_was_evaluated = .false.
            endif
        end select
        if (initial_was_evaluated) return

    end if

    ! We get here if .not. initial or the functional was not evaluated above
    select case (flag_xc)
    case (0)  ! HF, nothing to be done here
    case (1)  ! PBE0: only 3/4 of the PBE exchange is included
        call xc_partials_pbe0( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (3)  ! Perdew-Zunger LDA
        call xcpot_pz_lda(rho, local_xc_derivs, en_density_xc)
    case (4)  ! PW91 GGA
        ! FIXME: The following routine does not use the Hesse matrix of
        ! the density. Therefore, it returns only the energy density and
        ! its partial derivatives. To compute the local potential explicitly,
        ! use another routine (which actually requires the Hessian of the density).
        ! Johan has already supplied that subroutine. Need to integrate it here.
        call xc_partials_pw91gga( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (5)  ! PBE with VDW
        call xc_partials_pbe_vdw( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv, coord_current)
    case (6)  ! Perdew-Burke-Ernzerhoff PRL 77, 3865 (1996)
        ! FIXME: The following routine does not use the Hesse matrix of
        ! the density. Therefore, it returns only the energy density and
        ! its partial derivatives. To compute the local potential explicitly,
        ! use another routine (which actually requires the Hessian of the density).
        ! Johan has already supplied that subroutine. Need to integrate it here.
        call xc_partials_pbe( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (7)
        call xc_partials_hse( &
            hse_omega_pbe, rho, rho_gradient, en_density_xc, en_density_x, &
            en_density_c, local_xc_derivs, xc_gradient_deriv)
    case (8)  ! Perdew-Wang 1991 LDA
        call xcpot_pw91_lda( &
            rho, local_xc_derivs, en_density_xc, en_density_x, en_density_c)
    case (9)  ! BLYP
        call xc_partials_blyp( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (10)  ! B3LYP - rpa parametrized
        call xc_partials_b3lyp( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (11)  ! RPBE
        call xc_partials_rpbe( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (12)  ! revPBE
        call xc_partials_revpbe( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (13)  ! PBEsol0
        call xc_partials_pbesol0( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (14)  ! revPBE_vdw
        call xc_partials_revpbe_vdw( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv, coord_current)
    case (15)  ! VWN5 (this is the full VWN LDA parametrization)
        call xc_vwn5_lda( &
            rho, local_xc_derivs, en_density_xc, en_density_x, en_density_c)
    case (16)  ! VWN-GAUSS (this is the VWN RPA parametrization)
        call xc_vwn_lda( &
            rho, local_xc_derivs, en_density_xc, en_density_x, en_density_c)
    case (17)  ! PBEsol
        call xc_partials_pbesol( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (18)  ! mixed functional 0.52*PBE+0.48*RPBE, Phys. Rev. Lett. 108, 236104
        en_density_xc2 = 0.d0
        en_density_x2 = 0.d0
        en_density_c2 = 0.d0
        local_xc_derivs2 = 0.d0
        xc_gradient_deriv2 = 0.d0
        call xc_partials_pbe( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
        call xc_partials_rpbe( &
            rho, rho_gradient, en_density_xc2, en_density_x2, en_density_c2, &
            local_xc_derivs2, xc_gradient_deriv2)
        en_density_xc = en_density_xc*0.52 + en_density_xc2*0.48
        en_density_x = en_density_x*0.52 + en_density_x2*0.48
        en_density_c = en_density_c*0.52 + en_density_c2*0.48
        local_xc_derivs = local_xc_derivs*0.52 + local_xc_derivs2*0.48
        xc_gradient_deriv = xc_gradient_deriv*0.52 + xc_gradient_deriv2*0.48
    case (20)  ! AM05
        call xc_partials_am05( &
            rho, rho_gradient, en_density_xc, local_xc_derivs, xc_gradient_deriv)
    case (21)  ! PBEint
        call xc_partials_pbeint( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (22)  ! xPBE
        call xc_partials_xpbe( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (23)
        call xc_partials_lc_wpbeh( &
            hse_omega_pbe, rho, rho_gradient, en_density_xc, en_density_x, &
            en_density_c, local_xc_derivs, xc_gradient_deriv)
    case (25:28)
        ! Options for M06 subroutines: 1 = M06-L
        !                              2 = M06-HF
        !                              3 = M06
        !                              4 = M06-2X
        ! Our flag_xc options are in this order, so can do flag_xc-24.
        call xc_partials_m06_family( &
            rho, rho_gradient, kinetic_density, en_density_xc, &
            en_density_x(1), en_density_c, local_xc_derivs, xc_gradient_deriv, &
            xc_tau_deriv, flag_xc-24, tau_threshold)
    case (29)  ! B1LYP
        call xc_partials_b1lyp( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (30)  ! DFT part in XYG3 functional
        call xc_partials_xyg3_dft_part( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (31)  ! DFT part in xDH-PBE0 functional
        call xc_partials_xdhpbe0_dft_part( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (32)  ! DFT part in XYGJOS functional
        call xc_partials_xygjos_dft_part( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv )
    case (33)  ! DFT part in ZRS1 functional
        call xc_partials_zrs1_dft_part( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (34)  ! DFT part in XYG5 functional
        call xc_partials_xyg5_dft_part( &
            rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
            local_xc_derivs, xc_gradient_deriv)
    case (51:53)
        ! Options for TPSS subroutine: 1 = TPSS
        !                              2 = revTPSS
        !                              3 = TPSSloc
        ! Our flag_xc options are in this order, so can do flag_xc-50.
        call xc_partials_tpss_family( &
            rho, rho_gradient, kinetic_density, en_density_xc, &
            en_density_x(1), en_density_c, local_xc_derivs, xc_gradient_deriv, &
            xc_tau_deriv, flag_xc-50, tau_threshold)
    case (55:58)
        ! Options for M08 subroutines: 1 = M08-HX
        !                              2 = M08-SO
        !                              3 = M11
        !                              4 = M11-L
        ! Our flag_xc options are in this order, so can do flag_xc-54.
        call xc_partials_m08m11_family( &
            rho, rho_gradient, kinetic_density, en_density_xc, en_density_x(1), &
            en_density_c, local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, &
            flag_xc-54, tau_threshold)
    case (59)  ! SCAN meta-GGA functional
        call xc_partials_meta_scan( &
            rho, rho_gradient, kinetic_density, en_density_xc, &
            en_density_x(1), en_density_c, local_xc_derivs, xc_gradient_deriv, &
            xc_tau_deriv, 1, tau_threshold)
        if (n_spin == 2) en_density_x(2) = en_density_x(1)
    case (:-10)
        ! Recompose the libXC flags
        flag_libxc_x_id = (flag_xc/(-10))/1000
        if (flag_libxc_x_id .eq. 0) then
            ! Catch the scenario where we've defined a combined functional e.g. PBE0
            flag_libxc_x_id = mod((flag_xc/(-10)), 1000)
        else
            flag_libxc_c_id = mod((flag_xc/(-10)), 1000)
        endif
        
        call xc_partials_libxc( &
            rho, rho_gradient, kinetic_density, en_density_xc, &
            en_density_x(1), en_density_c, local_xc_derivs, xc_gradient_deriv, &
            xc_tau_deriv, flag_libxc_x_id, flag_libxc_c_id, tau_threshold)
    case (xc_dfauto_offset:xc_dfauto_offset+xc_flag_max)
        call xc_partials_dfauto( &
            flag_xc-xc_dfauto_offset, rho, rho_gradient, kinetic_density, &
            en_density_xc, en_density_x(1), en_density_c, &
            local_xc_derivs, xc_gradient_deriv, xc_tau_deriv)
        if (n_spin == 2) en_density_x(2) = en_density_x(1)
    case default
        write (stderr, *) "Chosen type of XC is not yet implemented."
        stop
    end select

    if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then !(Rundong) See my comments at the beginning.
       n_spin=temp
       local_xc_derivs(n_spin)=local_xc_derivs(1)
       en_density_x(n_spin)=en_density_x(1)
       xc_gradient_deriv(:,n_spin)=xc_gradient_deriv(:,1)
    endif
end subroutine evaluate_xc
