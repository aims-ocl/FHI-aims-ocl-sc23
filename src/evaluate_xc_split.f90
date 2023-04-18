!****s* FHI-aims/evaluate_xc_split
!  NAME
!  evaluate_xc_split
!  SYNOPSIS
      subroutine evaluate_xc_split &
      ( rho, rho_gradient, kinetic_density, &
        en_density_xc, &
        en_density_x, en_density_c, &
        local_x_derivs, local_c_derivs, &
        x_gradient_deriv, c_gradient_deriv, &
        x_tau_deriv, c_tau_deriv, &
        coord_current &
      )

!  PURPOSE
!  Subroutine evaluate_xc_split evaluates the exchange correlation potential
!  and energy contribution for a given density at one integration point.
!
!  VB: The present version supports LDA and gradient functionals.
!
!  USES
!
      use runtime_choices
      use dimensions
      use xc
      use constants
      use localorb_io, only: use_unit
      implicit none

!  ARGUMENTS

      real*8 rho(n_spin)
      real*8 rho_gradient(3,n_spin)
      real*8 kinetic_density(n_spin)
      real*8 en_density_xc
      real*8 en_density_x(n_spin)
      real*8 en_density_c
      real*8 local_x_derivs(n_spin)
      real*8 local_c_derivs(n_spin)
      real*8 x_gradient_deriv(3,n_spin)
      real*8 c_gradient_deriv(3,n_spin)
      real*8 x_tau_deriv(n_spin)
      real*8 c_tau_deriv(n_spin)
      real*8, optional :: coord_current(3)  !SAG

!  INPUTS
!   o flag_xc : Determines type of XC functional, see read_control.f
!   o rho     : Density at current point, given as spin up and spin down if spin-polarized
!   o rho_gradient : Density gradient at current point
!   o kinetic_density: Kinetic-energy density at the current point
!  OUTPUT
!   o en_density_xc : Exchange-correlation energy density at current point
!   o local_x_derivs : Local parts of the exchange potential:
!               Partial derivative of the exchange
!               energy functional by the density 
!   o local_c_derivs : Local parts of the correlation potential:
!               Partial derivative of the correlation
!               energy functional by the density plus spin-weighted
!               partial derivative of the correlation energy functional by spin
!   o x_gradient_deriv : Partial derivative of the exchange
!               energy functional by the modulus square of the density gradient
!   o c_gradient_deriv : Partial derivative of the correlation
!               energy functional by the modulus square of the density gradient
!   o x_tau_deriv: Partial derivative of the exchange energy functional
!                  by the kinetic-energy density
!   o c_tau_deriv: Partial derivative of the correlation energy functional
!                  by the kinetic-energy density
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


!  local variables

      real*8 local_xc_derivs(n_spin)
      logical :: xc_undefined

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

        x_tau_deriv = 0.d0
        c_tau_deriv = 0.d0

        if (xc_undefined) then

          en_density_xc = 0.d0
          en_density_x = 0.d0
          en_density_c = 0.d0

          local_x_derivs = 0.d0
          local_c_derivs = 0.d0
          x_gradient_deriv = 0.d0
          c_gradient_deriv = 0.d0

        else

!         calculate xc-potential, thus initializing new potential
          select case(flag_xc)
           case(0)
!        Hartree-Fock calculation, set the XC contribution to zero here.
            en_density_xc = 0.d0
            en_density_x = 0.d0
            en_density_c = 0.d0
            local_x_derivs = 0.d0
            local_c_derivs = 0.d0
            x_gradient_deriv = 0.d0
            c_gradient_deriv = 0.d0

!          case(1)
!!           Hybrid-PBE0 functional: here only 3/4 of the PBE exchange are included.
!            call xc_partials_pbe0 &
!            ( rho, rho_gradient, en_density_xc, &
!              en_density_x, en_density_c, local_xc_derivs, &
!              xc_gradient_deriv )
!
!          case(3)
!!           Perdew-Zunger LDA
!
!            call xcpot_pz_lda(rho, local_xc_derivs, &
!                         en_density_xc)
!
!          case(4)
!!           PW91_gga
!!           FIXME: The following routine does not use the Hesse matrix of
!!           the density. Therefore, it returns only the energy density and
!!           its partial derivatives. To compute the local potential explicitly,
!!           use another routine (which actually requires the Hessian of the density).
!!           Johan has already supplied that subroutine. Need to integrate it here.
!            call xc_partials_pw91gga &
!            ( rho, rho_gradient, en_density_xc,  en_density_x, &
!              en_density_c, local_xc_derivs, xc_gradient_deriv )
!
!         case(5)  
!            !PBE with VDW. SAG
!            
!            
!           call xc_partials_pbe_vdw &
!            ( rho, rho_gradient, en_density_xc,  en_density_x, &
!              en_density_c, local_xc_derivs, xc_gradient_deriv, &
!              coord_current)
!      
!
         case(1)
!           Burke-Ernzerhoff-Perdew CPL 265, 115 (1997)
            call xc_split_partials_pbe0 &
            ( rho, rho_gradient, en_density_xc,  en_density_x, &
              en_density_c, local_x_derivs,local_c_derivs, &
              x_gradient_deriv, c_gradient_deriv )

          case(6)
!           Perdew-Burke-Ernzerhoff PRL 77, 3865 (1996)
!           FIXME: The following routine does not use the Hesse matrix of
!           the density. Therefore, it returns only the energy density and
!           its partial derivatives. To compute the local potential explicitly,
!           use another routine (which actually requires the Hessian of the density).
!           Johan has already supplied that subroutine. Need to integrate it here.
            call xc_split_partials_pbe &
            ( rho, rho_gradient, en_density_xc,  en_density_x, &
              en_density_c, local_x_derivs,local_c_derivs, &
              x_gradient_deriv, c_gradient_deriv )

          case(8)
!           Perdew-Wang 1991 LDA
            call x_cpot_pw91_lda(rho, local_x_derivs, local_c_derivs, &
              local_xc_derivs, &
              en_density_x(1),  en_density_c,  en_density_xc)
          case(23)
            call xc_split_partials_lc_wpbeh &
            ( hse_omega_pbe, rho, rho_gradient, en_density_xc, &
              en_density_x, en_density_c, local_x_derivs, local_c_derivs, &
              x_gradient_deriv, c_gradient_deriv )
         case default
            write(use_unit,*) "Chosen type of XC is not yet implemented."
            stop
          end select

        end if

      end subroutine evaluate_xc_split
!---------------------------------------------------------------------
!****** 
