!----------------------------------------------------------------------
!************
!****s* FHI-aims/evaluate_post_xc
!  NAME
!  evaluate_post_xc
!  SYNOPSIS
      subroutine evaluate_post_xc &
      ( rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
        coord_current, kinetic_density_pointwise  &
      )

!  PURPOSE
!  Subroutine evaluate_post_xc evaluates the exchange correlation potential
!  and energy contribution for a given density at one integration point.
!
!
!  USES
!
      use runtime_choices
      use dimensions
      use xc
      use xc_library
      use constants
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      real*8 rho(n_spin)
      real*8 rho_gradient(3,n_spin)
      real*8 en_density_xc, en_density_x, en_density_c
      real*8 en_density_x2(n_spin)
      real*8 :: coord_current(3)
      real*8, optional :: kinetic_density_pointwise(n_spin)

!  INPUTS
!   o rho     : Density at current point, given as spin up and spin down if spin-polarized
!   o rho_gradient : Density gradient at current point
!   o kinetic_density : Kinetic density at current point (use_meta_gga)
!  OUTPUT
!   o en_density_xc : Exchange-correlation energy density at current point
!   o en_density_x : Exchange energy density
!   o en_density_c : Correlation energy density
!  LOCAL OUTPUT
!   o local_xc_derivs : Local parts of the XC potential:
!               Partial derivative of the exchange-correlation
!               energy functional by the density plus spin-weighted
!               partial derivative of the XC energy functional by spin
!   o xc_gradient_deriv : Partial derivative of the exchange-correlation
!               energy functional by the modulus square of the density gradient
!   o xc_tau_deriv : Partial derivative of the exchange-correlation energy
!                      functional by tau, the kinetic density      
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!    Updated by Andrew Logsdail, Jan 2015, for SCF Meta-GGA calculations
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

!  local variables

      logical :: xc_undefined

      real*8 local_xc_derivs(n_spin)
      real*8 xc_gradient_deriv(3,n_spin)
      real*8 xc_tau_deriv(n_spin)
  
      integer :: m06_option

!     counters

      integer :: i_spin

!  begin work

        xc_undefined = .true.

        do i_spin = 1, n_spin, 1
          ! This is true if both densities are below zero
          xc_undefined = xc_undefined .and. (rho(i_spin).le.0.d0)
        enddo

        do i_spin = 1, n_spin, 1
          ! This is true if one of the densities is less than zero
          xc_undefined = xc_undefined .or. (rho(i_spin).lt.0.d0)
        enddo

        en_density_xc = 0.d0
        en_density_x =0.d0
        en_density_c=0.d0

        local_xc_derivs = 0.d0
        xc_gradient_deriv = 0.d0
        xc_tau_deriv = 0.d0

        if (xc_undefined) then

           continue

        else

!         calculate xc-potential, thus initializing new potential
          
          select case(flag_post_xc)
          case(0)
!        Hartree-Fock calculation, set the XC contribution to zero here.

          case(1, 2, 5, 9)
!           M06 Family of functionals - post processing

            if (flag_post_xc.eq.1) then
               m06_option = 1 ! M06-L
            elseif (flag_post_xc.eq.9) then
               m06_option = 2 ! M06-HF
            elseif (flag_post_xc.eq.2) then
               m06_option = 3 ! M06
            elseif (flag_post_xc.eq.5) then
               m06_option = 4 ! M06-2X
            endif

            call xc_partials_m06_family &
            ( rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c, &
              local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, m06_option, tau_threshold)

          case(3)
!          pbe_vdw XC - post processing
            call xc_partials_pbe_vdw &
                 ( rho, rho_gradient, en_density_xc,  en_density_x2, &
                 en_density_c, local_xc_derivs, xc_gradient_deriv, &
                 coord_current)
            en_density_x = en_density_x2(1)+ en_density_x2(2) !strictly, en_density_x should be an array...
            
            
          case(4)
!          revpbe_vdw XC - post processing
            call xc_partials_revpbe_vdw &
                 ( rho, rho_gradient, en_density_xc, &
                 en_density_x2, en_density_c, local_xc_derivs, &
                 xc_gradient_deriv, coord_current )
            en_density_x = en_density_x2(1)+en_density_x2(2)

          case(6:8)
!           TPSS functional - post processing

            call xc_partials_tpss_family &
            ( rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c, &
              local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, flag_post_xc-5, tau_threshold)

          case(10:13)
!           M08/M11 Family of functionals - post processing
              ! Options for M08/M11 subroutines: 1 = M08-HX
              !                                  2 = M08-SO
              !                                  3 = M11
              !                                  4 = M11-L

            call xc_partials_m08m11_family &
            ( rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c, &
              local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, flag_post_xc-9, tau_threshold)           
          case(14) ! SCAN functional - post-processing
! Removing duplicate functionality. AJL/Dec2016
!              call xc_meta_scan ( 1, &
              call xc_partials_meta_scan ( &
                  rho, rho_gradient, kinetic_density_pointwise, &
                  en_density_xc, en_density_x, en_density_c, &
                  local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, 1, tau_threshold )
          case (xc_dfauto_offset:xc_dfauto_offset+xc_flag_max)
              call xc_partials_dfauto( &
                  flag_post_xc-xc_dfauto_offset, &
                  rho, rho_gradient, kinetic_density_pointwise, &
                  en_density_xc, en_density_x, en_density_c, &
                  local_xc_derivs, xc_gradient_deriv, xc_tau_deriv)
!          case (15)  ! B3LYP - rpa parametrized
!              call xc_partials_b3lyp( &
!                  rho, rho_gradient, en_density_xc, en_density_x2, &
!                  en_density_c, local_xc_derivs, xc_gradient_deriv)
!
!              en_density_x = en_density_x2(1)+ en_density_x2(2)
          case default
            write(use_unit,*) "Chosen type of XC post-processing is not yet implemented."
            stop
          end select

        end if


      end subroutine evaluate_post_xc
!---------------------------------------------------------------------

