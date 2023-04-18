!****h* FHI-aims/xc
!  NAME
!    xc
!  SYNOPSIS

module xc

  !  PURPOSE
  !  xc.f90 provides functions for the calculation of the XC energy and/or potential
  !
  !  USES
  !  Note:  I (WPH) ordinarily make a point of eliminating all module-level use 
  !         statements, but this module is 7000 lines of code, the modules used are low 
  !         level, and the variables used are small in number and very common.  So I'm 
  !         letting it slide here.
  use dimensions,      only : n_spin, use_hartree_fock
  use runtime_choices, only : spin_treatment, hybrid_coeff, lc_dielectric_constant
  implicit none
  !  AUTHOR
  !    FHI-aims team.
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften 
  !  e.V. Please note that any use of the "FHI-aims-Software" is subject to 
  !  the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !****


contains
  !------------------------------------------------------------------
  !****s* xc/xcpot_pz_lda
  !  NAME
  !    xcpot_pz_lda
  !  SYNOPSIS

  subroutine xcpot_pz_lda &
       ( density, pot_xc, en_density_xc &

       )
    !  PURPOSE
    !    Calculates pz_lda xc energy and its derivatives with respect to the density.
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8 :: density(n_spin)    
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: pot_xc

    !  INPUTS
    !   o density --  Full density rho(r)
    !  OUTPUT
    !   o pot_xc --  Exchange-correlation potential at current grid point
    !   o en_density_xc -- Exchange-correlation energy density at current grid point
    !   
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    !  local variables:
    !       en_density_x : Exchange energy density at current grid point
    !       pot_x :        Exchange potential at current grid point


    real*8 :: dVxc_drho(n_spin)

    !       for spin-polarized version

    real*8 :: rho
    real*8 :: zeta
    real*8 :: vxc_av
    real*8 :: dvxc

    !  begin work

    if (spin_treatment.eq.0) then
       ! unpolarized version - use routine from fhipp code


       call pz_lda(density(1),en_density_xc,en_density_x, en_density_c, & 
                   pot_xc(1),dVxc_drho(1))

    else if (spin_treatment.eq.1) then
       ! Spin-polarized version.

       rho = density(1)+density(2)
       zeta = (density(1)-density(2)) / rho

       call stvxc_spin &
            (rho, zeta, vxc_av, dvxc, en_density_xc)

       pot_xc(1) = vxc_av + dvxc/2.d0
       pot_xc(2) = vxc_av - dvxc/2.d0

    end if

  end subroutine xcpot_pz_lda
  !******
  !---------------------------------------------------------------------
  !****s* xc/xcpot_pw91_lda
  !  NAME
  !    xcpot_pw91_lda
  !  SYNOPSIS

  subroutine xcpot_pw91_lda &
       ( density, pot_xc, en_density_xc, &
         en_density_x, en_density_c &
       )

    !  PURPOSE
    !  Subroutine xcpot_pw91_lda provides the exchange-correlation potential
    !  according to the Perdew-Wang 1991 parametrisation of Ceperley-Alder LDA
    !  USES

    use constants

    implicit none
    !  ARGUMENTS

    real*8, dimension(n_spin) :: density
    real*8, dimension(n_spin) :: pot_xc
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c

    !  INPUTS
    !   o density -- Full density rho(r)
    !   o en_density_c : Correlation energy density at current grid point
    !   o en_density_x : Exchange energy density at current grid point
    !                    (possibly different for spin up, spin down)
    !  OUTPUT
    !   o pot_xc -- XC potential (possibly different for spin up, spin down) 
    !   o en_density_xc -- Exchange-correlation energy density at current grid 
    !                      point
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE





    !  local variables
    !       aux_density : Extra variable to specifiy the actual input density
    !                     which enters each subroutine
    !       pot_c :        correlation potential at current grid point 
    !                      (spin-up and spin-down)
    !       inv_k_fermi :  Inverse of the Fermi wave vector of a homogeneous
    !                      electron gas of given density
    !       r_seitz:       Wigner-Seitz radius
    !       zeta :         would be spin polarisation
    !       ec_rs :        derivative of en_density_c w.r.t. r_seitz
    !       ec_zeta :      derivative of en_density_c w.r.t. zeta
    !       alpha_c :      "correlation contribution to spin stiffness"

    real*8 :: aux_density

    real*8, dimension(n_spin) :: pot_c
    real*8 :: inv_k_f
    real*8 :: r_seitz
    real*8 :: zeta = 0.d0
    real*8 :: ec_rs = 0.d0
    real*8 :: ec_zeta = 0.d0
    real*8 :: alpha_c = 0.d0

    !       counters

    integer :: i_spin

    !  begin work

    ! Exchange energy first, for both spin channels
    do i_spin = 1, n_spin, 1
       aux_density = dble(n_spin) * density(i_spin)
       call xlda( aux_density, pot_xc(i_spin), &
            en_density_x(i_spin) )
    enddo

    if (spin_treatment.eq.1) then

       aux_density = density(1)+density(2)
       zeta = ( density(1)-density(2) ) / aux_density

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*density(i_spin)
       enddo

       en_density_xc = en_density_xc / aux_density

    else

       en_density_xc = en_density_x(1)
       aux_density = density(1)
       zeta = 0.d0

    end if

    ! Correlation energy is next

    inv_k_f = (pisq3*aux_density)**third
    r_seitz = const_rs/inv_k_f

    if (spin_treatment.eq.0) then
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c(1), &
            pot_c(1), ec_rs, ec_zeta, alpha_c )
    else
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c(1), &
            pot_c(2), ec_rs, ec_zeta, alpha_c )
    end if

    en_density_xc = &
         en_density_xc + en_density_c

    do i_spin = 1, n_spin, 1
       pot_xc(i_spin) = pot_xc(i_spin) + pot_c(i_spin)
    enddo

  end subroutine xcpot_pw91_lda

  !******
  !---------------------------------------------------------------------
  !****s* xc/x_cpot_pw91_lda
  !  NAME
  !    x_cpot_pw91_lda
  !  SYNOPSIS

  subroutine x_cpot_pw91_lda &
       ( density, pot_x, pot_c, pot_xc, en_density_x, en_density_c, en_density_xc &
       )

    !  PURPOSE
    !  Subroutine x_cpot_pw91_lda provides the exchange and correlation potential
    !  according to the Perdew-Wang 1991 parametrisation of Ceperley-Alder LDA
    !  USES

    use constants

    implicit none
    !  ARGUMENTS

    real*8 density(n_spin)
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: pot_xc

    !  INPUTS
    !   o density -- Full density rho(r)
    !  OUTPUT
    !   o en_density_x -- Exchange energy density at current grid point
    !   o en_density_c -- Correlation energy density at current grid point
    !   o pot_x -- X potential (possibly different for spin up, spin down 
    !   o pot_c -- C potential (possibly different for spin up, spin down 
    !   o pot_xc -- XC potential (possibly different for spin up, spin down 
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE





    !  local variables
    !       aux_density : Extra variable to specifiy the actual input density which enters each subroutine
    !       en_density_c : Correlation energy density at current grid point
    !       pot_c :        correlation potential at current grid point (spin-up and spin-down)
    !       inv_k_fermi :  Inverse of the Fermi wave vector of a homogeneous electron gas of given density
    !       r_seitz:       Wigner-Seitz radius
    !       zeta :         would be spin polarisation
    !       ec_rs :        derivative of en_density_c w.r.t. r_seitz
    !       ec_zeta :      derivative of en_density_c w.r.t. zeta
    !       alpha_c :      "correlation contribution to spin stiffness"

    real*8 :: aux_density
    real*8 :: en_density_x
    real*8, dimension(n_spin) :: en_density_x0

    real*8 :: en_density_c
    real*8, dimension(n_spin) :: pot_c, pot_x
    real*8, dimension(n_spin) :: pot_c0, pot_x0
    real*8 :: inv_k_f
    real*8 :: r_seitz
    real*8 :: zeta = 0.d0
    real*8 :: ec_rs = 0.d0
    real*8 :: ec_zeta = 0.d0
    real*8 :: alpha_c = 0.d0

    !       counters

    integer :: i_spin

    !  begin work

    ! Exchange energy first, for both spin channels
    do i_spin = 1, n_spin, 1
       aux_density = dble(n_spin) * density(i_spin)
       call xlda( aux_density, pot_x0(i_spin), &
            en_density_x0(i_spin) )
    enddo

    if (spin_treatment.eq.1) then

       aux_density = density(1)+density(2)
       zeta = ( density(1)-density(2) ) / aux_density

       en_density_x = 0.d0
       do i_spin = 1,n_spin,1
          en_density_x = en_density_x &
               + en_density_x0(i_spin)*density(i_spin)
       enddo

       en_density_x = en_density_x / aux_density

    else

       en_density_x = en_density_x0(1)
       aux_density = density(1)
       zeta = 0.d0

    end if

    ! Correlation energy is next

    en_density_c =0.d0
    inv_k_f = (pisq3*aux_density)**third
    r_seitz = const_rs/inv_k_f

    if (spin_treatment.eq.0) then
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c0(1), &
            pot_c0(1), ec_rs, ec_zeta, alpha_c )
    else
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c0(1), &
            pot_c0(2), ec_rs, ec_zeta, alpha_c )
    end if

    en_density_xc = &
         en_density_c + en_density_x

    pot_xc=0.0
    pot_x=0.0
    pot_c=0.0
    do i_spin = 1, n_spin, 1
      pot_x(i_spin) = pot_x(i_spin) + pot_x0(i_spin)
      pot_c(i_spin) = pot_c(i_spin) + pot_c0(i_spin)
    enddo
    do i_spin = 1, n_spin, 1
      pot_xc(i_spin) = pot_xc(i_spin) + pot_c(i_spin) + pot_x(i_spin)
    enddo

  end subroutine x_cpot_pw91_lda

  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_vwn5_lda
  !  NAME
  !    xc_vwn5_lda
  !  SYNOPSIS
  subroutine xc_vwn5_lda &
       ( density, pot_xc, en_density_xc, &
         en_density_x, en_density_c &
       )
    !  PURPOSE
    !  Subroutine xc_vwn5_lda provides the exchange-correlation potential
    !  according to the VWN5 (Vosko-Wilk-Nusair (1980)) parametrisation (as they name it in the Gaussian code)
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8 density(n_spin)
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: pot_xc

    !  INPUTS
    !   o density -- Full density rho(r)
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density at current grid point
    !   o pot_xc --  XC potential (possibly different for spin up, spin down
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE




    !  local variables

    !       aux_density : Extra variable to specifiy the actual input density which enters each subroutine

    !       en_density_c : Correlation energy density at current grid point
    !       pot_c :        correlation potential at current grid point (spin-up and spin-down)
    !       inv_k_fermi :  Inverse of the Fermi wave vector of a homogeneous electron gas of given density
    !       r_seitz:       Wigner-Seitz radius
    !       zeta :         would be spin polarisation
    !       ec_rs :        derivative of en_density_c w.r.t. r_seitz
    !       ec_zeta :      derivative of en_density_c w.r.t. zeta
    !       alpha_c :      "correlation contribution to spin stiffness"

    real*8 :: aux_density
    real*8, dimension(n_spin) :: pot_c
    real*8 :: zeta = 0.d0


    real*8, dimension(3) :: dork
    !       counters

    integer :: i_spin

    !  begin work

    ! Exchange energy first, for both spin channels
    do i_spin = 1, n_spin, 1
       aux_density = dble(n_spin) * density(i_spin)
       call xlda( aux_density, pot_xc(i_spin), &
            en_density_x(i_spin) )
    enddo

    if (spin_treatment.eq.1) then

       aux_density = density(1)+density(2)
       zeta = ( density(1)-density(2) ) / aux_density

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*density(i_spin)
       enddo

       en_density_xc = en_density_xc / aux_density

    else

       en_density_xc = en_density_x(1)
       aux_density = density(1)
       zeta = 0.d0

    end if

    ! Correlation energy is next
    ! to be consistent with vexcor, that is how it should be called...
    if (spin_treatment.eq.0) then

       ! no spin polarization - only the full density (density(1)) is stored
       dork(1)=0.5d0*density(1)
       dork(2)=0.5d0*density(1)

       dork(3)=dork(1)+dork(2)

       call cepvwn &
            (dork, en_density_c, pot_c, 0)
    else

       ! with spin
       dork(1)=density(1)
       dork(2)=density(2)
       dork(3)=dork(1)+dork(2)

       call cepvwn &
            (dork, en_density_c, pot_c, 1)
    end if

    en_density_xc = &
         en_density_xc + en_density_c

    do i_spin = 1, n_spin, 1
       pot_xc(i_spin) = pot_xc(i_spin) + pot_c(i_spin)
    enddo

  end subroutine xc_vwn5_lda
  !******

  !---------------------------------------------------------------------

  !****s* xc/xc_vwn_lda
  !  NAME
  !    xc_vwn_lda
  !  SYNOPSIS

  subroutine xc_vwn_lda &
       ( density, pot_xc, en_density_xc, &
         en_density_x, en_density_c &
       )

    !  PURPOSE
    !  Subroutine xc_vwn_lda provides the exchange-correlation potential
    !  according to the VWN-RPA parametrisation (Vosko-Wilk-Nusair (1980)), as they name it in the GAUSSIAN code.
    !  
    !  USES

    use constants

    !  ARGUMENTS

    implicit none
    real*8, dimension(n_spin) :: density
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: pot_xc

    !  INPUTS
    !   o density -- Full density rho(r)
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density at current grid point
    !   o pot_xc -- XC potential (possibly different for spin up, spin down
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE





    !  local variables

    !       aux_density : Extra variable to specifiy the actual input density which enters each subroutine

    !       en_density_c : Correlation energy density at current grid point
    !       pot_c :        correlation potential at current grid point (spin-up and spin-down)
    !       inv_k_fermi :  Inverse of the Fermi wave vector of a homogeneous electron gas of given density
    !       r_seitz:       Wigner-Seitz radius
    !       zeta :         would be spin polarisation
    !       ec_rs :        derivative of en_density_c w.r.t. r_seitz
    !       ec_zeta :      derivative of en_density_c w.r.t. zeta
    !       alpha_c :      "correlation contribution to spin stiffness"

    real*8 :: aux_density
    real*8, dimension(n_spin) :: pot_c
    real*8 :: zeta = 0.d0


    !       counters

    integer :: i_spin

    !  begin work

    ! Exchange energy first, for both spin channels
    do i_spin = 1, n_spin, 1
       aux_density = dble(n_spin) * density(i_spin)
       call xlda( aux_density, pot_xc(i_spin), &
            en_density_x(i_spin) )
    enddo

    if (spin_treatment.eq.1) then

       aux_density = density(1)+density(2)
       zeta = ( density(1)-density(2) ) / aux_density

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*density(i_spin)
       enddo

       en_density_xc = en_density_xc / aux_density

    else

       en_density_xc = en_density_x(1)
       aux_density = density(1)
       zeta = 0.d0

    end if

    ! Correlation energy is next
    ! to be consistent with vexcor, that is how it should be called...
    !                dork(1)=density(1)
    !                dork(2)=density(2)
    !                dork(3)=dork(1)+dork(2)

    !                if (spin_treatment.eq.0) then
    call vwnrpa_c_derivs &
         (density, en_density_c, pot_c, n_spin)
    !                else
    !                  call cepvwn_gauss
    !     +            (dork, en_density_c, pot_c, 1)
    !                end if

    en_density_xc = &
         en_density_xc + en_density_c

    do i_spin = 1, n_spin, 1
       pot_xc(i_spin) = pot_xc(i_spin) + pot_c(i_spin)
    enddo

  end subroutine xc_vwn_lda
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_pbe
  !  NAME
  !    xc_partials_pbe
  !  SYNOPSIS

  subroutine xc_partials_pbe &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, &
         local_xc_derivs, xc_gradient_deriv )

    !  PURPOSE
    !  Subroutine xc_partials_pbe provides the exchange-correlation energy functional
    !  of Perdew, Burke, Ernzerhoff, PRL 77, 3865 (1996) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ??
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE






    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho=0.0

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin) .gt. NUM_ZERO) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms
          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.NUM_ZERO) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_pbe


  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_split_partials_pbe
  !  NAME
  !    xc_split_partials_pbe
  !  SYNOPSIS

  subroutine xc_split_partials_pbe &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, &
         local_x_derivs, local_c_derivs, &
         x_density_gradient_deriv, c_density_gradient_deriv )

    !  PURPOSE
    !  Subroutine xc_split_partials_pbe provides the exchange-correlation energy functional
    !  of Perdew, Burke, Ernzerhoff, PRL 77, 3865 (1996) and its partial
    !  derivatives by the density and the density gradient. The only difference to
    !  xc_partials_pbe is that the exchange and correlation components of the xc potential 
    !  the gradient are separated.
    !
    !  We need this feature in some GW-type analysis -- XR, 24.09.2011
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_x_derivs
    real*8, dimension(n_spin) :: local_c_derivs
    real*8, dimension(3,n_spin) :: x_density_gradient_deriv
    real*8, dimension(3,n_spin) :: c_density_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_x_derivs -- Partial derivative d/d(rho) of e_x(rho, |grad(rho)|^2)
    !   o local_c_derivs -- Partial derivative d/d(rho) of e_c(rho, |grad(rho)|^2)
    !   o x_density_gradient_deriv -- Partial derivative d/d(|grad(rho)|) of e_x(rho, |grad(rho)|^2)
    !   o c_density_gradient_deriv -- Partial derivative d/d(|grad(rho)|) of e_c(rho, |grad(rho)|^2)
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE






    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
!    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms
          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), local_x_derivs(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          local_x_derivs(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
!       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
!            c_density_deriv
       local_c_derivs(i_spin) = c_density_deriv 
       do i_coord = 1,3,1
!          xc_gradient_deriv(i_coord, i_spin) = &
!               dble(n_spin) * x_gradient_deriv(i_spin) * &
!               rho_gradient(i_coord,i_spin) + &
!               c_gradient_deriv*total_rho_gradient(i_coord)
          x_density_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) 

          c_density_gradient_deriv(i_coord, i_spin) = &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_c_derivs(1) = c_density_deriv - &
            (spin_pol-1) * c_spin_deriv
       local_c_derivs(2) = c_density_deriv - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_split_partials_pbe


  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_pbe_vdw
  !  NAME
  !    xc_partials_pbe_vdw
  !  SYNOPSIS

  subroutine xc_partials_pbe_vdw &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, &
         local_xc_derivs, xc_gradient_deriv, coord_current )

 



    !  PURPOSE
    !  Subroutine xc_partials_pbe_vdw is a copy of xc_partials_pbe, but uses pbe exchange, lda correlation and nonlocal correlation.  SAG

    use constants
    use dimensions,      only : use_nlcorr_in_xc
    use octree_routines, only : cube_ready
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8, optional :: coord_current(3)

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ??
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE






    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin
    
    real*8 :: epsilon, epsilon2, epsilon3
    
    !  begin work

    epsilon = 0d0
    epsilon2 = 0d0
    epsilon3 = 0d0
   

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      squared_grad_rho = 0d0   !for pure LDA. SAG
       call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if
    
    c_gradient_deriv = 0d0  !correlation terms now strictly LDA. 



!nonlocal correlation.  SAG
    
    
    if(Present(coord_current).and.cube_ready.and.use_nlcorr_in_xc)then
       
       call nlcorr(coord_current,full_rho,total_rho_gradient,epsilon,epsilon2, epsilon3) 
       
    endif
    
    en_density_c = en_density_c + epsilon   !true energy density 
    c_density_deriv = c_density_deriv + epsilon2
    c_gradient_deriv = c_gradient_deriv + epsilon3





    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_pbe_vdw





  !******
  !---------------------------------------------------------------------
  !****s* xc/x_c_partials_pbe
  !  NAME
  !    x_c_partials_pbe
  !  SYNOPSIS

  subroutine x_c_partials_pbe &
       ( rho, rho_gradient, en_density_x, en_density_c, en_density_xc, &
       local_xc_derivs, xc_gradient_deriv )

    !  PURPOSE
    !  Subroutine x_c_partials_pbe provides the exchange and correlation energy functional
    !  of Perdew, Burke, Ernzerhoff, PRL 77, 3865 (1996) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs, local_x_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ??
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE






    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x0
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x0(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x0(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_x = 0.d0
       do i_spin = 1,n_spin,1
          en_density_x = en_density_x &
               + en_density_x0(i_spin)*rho(i_spin)
       enddo

       en_density_x = en_density_x / full_rho

    else

       en_density_x = en_density_x0(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_x + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

  end subroutine x_c_partials_pbe
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_pbesol
  !  NAME
  !    xc_partials_pbesol
  !  SYNOPSIS

  subroutine xc_partials_pbesol &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_pbesol provides the exchange-correlation energy functional
    !  proposed in arXiv:0711.0156v2 (2008/02) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv


    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE





    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbesol_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbesol_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if
    !       add relevant terms together
    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_pbesol
  
!@@edu>
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_pbeint
  !  NAME
  !    xc_partials_pbeint
  !  SYNOPSIS
  subroutine xc_partials_pbeint &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_pbeint provides the exchange-correlation energy functional
    !  of PBEint proposed by  E. Fabiano, L. A. Constantin,and F. Della Sala
    !  in Phys. Rev. B 82, 113104 (2010),  and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    E. Fabiano.
    !  HISTORY
    !    February 2011.
    !  SEE ALSO
    !    Phys. Rev. B 82, 113104 (2010)
    !  SOURCE




    !  Local variables
    !
    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbeint_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbeint_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
         c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_pbeint


!@@edu< 
   
  !****s* xc/xc_partials_pw91gga
  !  NAME
  !    xc_partials_pw91gga
  !  SYNOPSIS

  subroutine xc_partials_pw91gga &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, &
         local_xc_derivs, xc_gradient_deriv )

    !  PURPOSE
    !  Subroutine xc_partials_pw91gga provides the exchange-correlation energy functional
    !  of Perdew and Wang Phys Rev. B 46 6671 (1992) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  This implementation is based on the version in external/vexcor.f 
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ??
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE






    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PW91 subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                (dble(n_spin)*rho_gradient(i_coord,i_spin))**2
          enddo

          !           get exchange terms
          call exchpw91_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
    
      call corpw91_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_pw91gga

  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_rpbe
  !  NAME
  !    xc_partials_rpbe
  !  SYNOPSIS

  subroutine xc_partials_rpbe &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_rpbe provides the exchange-correlation energy functional
    !  of RPBE proposed by Hammer and Norskov PRB 59, 7413 (1999) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    



    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchrpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_rpbe
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_revpbe
  !  NAME
  !    xc_partials_revpbe
  !  SYNOPSIS

  subroutine xc_partials_revpbe &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_revpbe provides the exchange-correlation energy functional
    !  of RevPBE proposed by Y. Zhang and W. Yang, Phys. Rew. Lett. 80 890 (1998) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE




    !  Local variables
    !
    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchrevpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_revpbe

  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_revpbe
  !  NAME
  !    xc_partials_revpbe
  !  SYNOPSIS

  subroutine xc_partials_xpbe &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_revpbe provides the exchange-correlation energy functional
    !  of xPBE proposed by X. Xu and W. A. Gordard III, JCP. 9 4068 (2004) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE




    !  Local variables
    !
    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchxpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corxpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_xpbe
 



!******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_revpbe_vdw_vdw
  !  NAME
  !    xc_partials_revpbe_vdw
  !  SYNOPSIS

  subroutine xc_partials_revpbe_vdw &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv, coord_current &
       )


    !  PURPOSE
    !  Subroutine xc_partials_revpbe_vdw is a copy of xc_partials_revpbe, but which uses revpbe exchange, lda correlation and nonlocal correlation. SAG

    !  USES

    use constants
    use dimensions,      only : use_nlcorr_in_xc
    use octree_routines, only : cube_ready
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
      real*8, optional :: coord_current(3)

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE




    !  Local variables
    !
    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin
   
    real*8 ::  epsilon, epsilon2, epsilon3
    
    !  begin work

    epsilon = 0.d0
    epsilon2 = 0d0
    epsilon3 = 0d0
   

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchrevpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      squared_grad_rho = 0d0   !for pure LDA. SAG
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if
    
    c_gradient_deriv = 0d0  !correlation terms now pure LDA. 
   

    
    !nonlocal correlation. SAG
    if(Present(coord_current).and.cube_ready.and.use_nlcorr_in_xc)then
       
       call nlcorr(coord_current,full_rho,total_rho_gradient,epsilon,epsilon2, epsilon3) 
       
    endif
    
    en_density_c = en_density_c + epsilon   !true energy density 
    c_density_deriv = c_density_deriv + epsilon2
    c_gradient_deriv = c_gradient_deriv + epsilon3




    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_revpbe_vdw





 !******
  !---------------------------------------------------------------------
  !****s* xc/x_c_partials_revpbe
  !  NAME
  !    x_c_partials_revpbe
  !  SYNOPSIS

  subroutine x_c_partials_revpbe &
       ( rho, rho_gradient, en_density_x, en_density_c, en_density_xc, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine x_c_partials_revpbe provides the exchange and correlation energy functional
    !  of RevPBE proposed by Y. Zhang and W. Yang, Phys. Rew. Lett. 80 890 (1998) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE




    !  Local variables
    !
    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x0
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c0
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchrevpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x0(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x0(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_x = 0.d0
       do i_spin = 1,n_spin,1
          en_density_x = en_density_x &
               + en_density_x0(i_spin)*rho(i_spin)
       enddo

       en_density_x = en_density_x / full_rho

    else

       en_density_x = en_density_x0(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_x + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if


  end subroutine x_c_partials_revpbe
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_pbesol0
  !  NAME
  !    xc_partials_pbesol0
  !  SYNOPSIS

  subroutine xc_partials_pbesol0 &
       ( rho, rho_gradient, en_density_xc, & 
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_pbesol0 provides the exchange-correlation energy functional
    !  of Perdew et al. proposed in arXiv:0711.0156v2 (2008/02) and its
    !  derivatives by the density and the density gradient + HF hybridization.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE





    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbesol_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo

       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
    call corpbesol_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together
    ! PBEsol0: here we add 3/4 of the PBE exchange to the PBE correlation.
    en_density_xc = (1.d0 - hybrid_coeff) * en_density_xc &
         + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = (1.d0 - hybrid_coeff) * &
            x_density_deriv(i_spin) +  c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               (1.d0 - hybrid_coeff) * dble(n_spin) * &
               x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if
    do i_spin = 1, n_spin, 1
      en_density_x(i_spin)=(1.d0 - hybrid_coeff) * en_density_x(i_spin)
    enddo

    !  that's all folks

  end subroutine xc_partials_pbesol0
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_pbe0
  !  NAME
  !    xc_partials_pbe0
  !  SYNOPSIS

  subroutine xc_partials_pbe0 &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
         xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_pbe0 provides the exchange-correlation energy functional
    !  of Burke, Ernzerhoff, Perdew, CPL 265, 115 (1997) and its partial
    !  derivatives by the density and the density gradient + HF hybridization.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE







    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together
    ! PBE0: here we add 3/4 of the PBE exchange to the PBE correlation.
    en_density_xc = (1.d0 - hybrid_coeff) * en_density_xc &
         + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = (1.d0 - hybrid_coeff) * &
            x_density_deriv(i_spin) +  c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               (1.d0 - hybrid_coeff) * dble(n_spin) * &
               x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    do i_spin = 1, n_spin, 1
     en_density_x(i_spin)=(1.d0 - hybrid_coeff) * en_density_x(i_spin)
    enddo

    !  that's all folks

  end subroutine xc_partials_pbe0
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_split_partials_pbe0
  !  NAME
  !    xc_split_partials_pbe0
  !  SYNOPSIS

  subroutine xc_split_partials_pbe0 &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, &
         local_x_derivs, local_c_derivs, &
         x_density_gradient_deriv, c_density_gradient_deriv )

    !  PURPOSE
    !  Subroutine xc_split_partials_pbe0 provides the hybrid exchange-correlation energy functional
    !  of Burke, Ernzerhoff, Perdew, CPL 265, 115 (1997) and its partial
    !  derivatives by the density and the density gradient. The only difference to
    !  xc_partials_pbe0 is that the exchange and correlation components of the xc potential 
    !  the gradient are separated.
    !
    !  We need this feature in some GW-type analysis -- XR, 26.09.2011
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_x_derivs
    real*8, dimension(n_spin) :: local_c_derivs
    real*8, dimension(3,n_spin) :: x_density_gradient_deriv
    real*8, dimension(3,n_spin) :: c_density_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_x_derivs -- Partial derivative d/d(rho) of e_x(rho, |grad(rho)|^2)
    !   o local_c_derivs -- Partial derivative d/d(rho) of e_c(rho, |grad(rho)|^2)
    !   o x_density_gradient_deriv -- Partial derivative d/d(|grad(rho)|) of e_x(rho, |grad(rho)|^2)
    !   o c_density_gradient_deriv -- Partial derivative d/d(|grad(rho)|) of e_c(rho, |grad(rho)|^2)
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
!    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms
          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), local_x_derivs(i_spin), &
               x_gradient_deriv(i_spin) &
               )

          local_x_derivs(i_spin)=(1.d0-hybrid_coeff)*local_x_derivs(i_spin)
          en_density_x(i_spin)=(1.d0-hybrid_coeff)*en_density_x(i_spin)
          x_gradient_deriv(i_spin)=(1.d0-hybrid_coeff)*x_gradient_deriv(i_spin)
       else

          en_density_x(i_spin) = 0.d0
          local_x_derivs(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
!       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
!            c_density_deriv
       local_c_derivs(i_spin) = c_density_deriv 
       do i_coord = 1,3,1
!          xc_gradient_deriv(i_coord, i_spin) = &
!               dble(n_spin) * x_gradient_deriv(i_spin) * &
!               rho_gradient(i_coord,i_spin) + &
!               c_gradient_deriv*total_rho_gradient(i_coord)
          x_density_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin)

          c_density_gradient_deriv(i_coord, i_spin) = &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_c_derivs(1) = c_density_deriv - &
            (spin_pol-1) * c_spin_deriv
       local_c_derivs(2) = c_density_deriv - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks
  end subroutine xc_split_partials_pbe0
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_hse
  !  NAME
  !    xc_partials_hse
  !  SYNOPSIS

  subroutine xc_partials_hse &
       (hse_omega_pbe, rho, rho_gradient, & 
        en_density_xc, en_density_x, en_density_c, &
        local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_hse provides the hybrid exchange-correlation energy functional
    !  of J. Heyd, G. E. Scuseria, and M. Ernzerhof, J. Chem. Phys. 118, 8207 (2003) and its partial
    !  derivatives by the density and the density gradient + HF hybridization.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8 :: hse_omega_pbe

    !  INPUTS
    !   o hse_omega_pbe -- ???
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE







    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)
    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    aux_squared_grad_rho=0.0

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       !if (rho(i_spin).gt.0.d0) then
       if (rho(i_spin).gt.1.d-20) then !Actual enforcing zero caused problem with restart of HSE calculations
 

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo
          !
          !           get exchange terms
          !
          !     NOTICE, that the output is allready :
          !
          !     en_density_pbe - hybrid_coef * en_density_pbe_shortrange,
          !
          !     x_density_derive_pbe - hybrid_coef * x_density_derive_pbe_shortrange
          !
          !         and
          !
          !     x_gradient_derive_pbe - hybrid_coef * x_gradient_derive_pbe_shortrange
          !
          call exch_hse_derivs &
               (hybrid_coeff,hse_omega_pbe,aux_density,aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_hse
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_lc_wpbeh
  !  NAME
  !    xc_partials_lc_wpbeh
  !  SYNOPSIS

  subroutine xc_partials_lc_wpbeh &
       (lc_wpbeh_omega, rho, rho_gradient, &
        en_density_xc, en_density_x, en_density_c, &
        local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_lc_wpbeh provides the long-range corrected hybrid exchange-correlation energy functional
    !  base on wPBE and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8 :: lc_wpbeh_omega


    !  INPUTS
    !   o lc_wpbeh_omega -- screening factor
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    Lukas Gallandi
    !  HISTORY
    !    ???
    !  SEE ALSO
    !    ???
    !  SOURCE
    !    FHI-aims team (subroutine 'xc_partials_hse') as well as other subroutines. Basically it's just
    !    a rearrangement of exisiting subroutines.








        !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)
    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho


    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    aux_squared_grad_rho=0.0

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       !if (rho(i_spin).gt.0.d0) then
       if (rho(i_spin).gt.1.d-20) then !Actual enforcing zero caused problem with restart of LC-wPBEh calculations


          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo
          !
          !           get exchange terms
          !
          !     NOTICE, that the output is allready :
          !
          !     (1/epsilon - hybrid_coef) * en_density_wpbe_shortrange + (1-1/epsilon)* en_density_pbe,
          !
          !     (1/epsilon - hybrid_coef)  * x_density_derive_wpbe_shortrange+ (1-1/epsilon)* x_density_derive_pbe
          !
          !         and
          !
          !     (1/epsilon - hybrid_coef)  * x_gradient_derive_wpbe_shortrange+ (1-1/epsilon)* x_gradient_derive_pbe,
          !
          call exch_lc_wpbeh_derivs &
               (hybrid_coeff,lc_wpbeh_omega,aux_density,aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin), lc_dielectric_constant &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then


       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_xc_derivs(1) = local_xc_derivs(1) - &
            (spin_pol-1) * c_spin_deriv
       local_xc_derivs(2) = local_xc_derivs(2) - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_partials_lc_wpbeh
  !******************************************
  !****s* xc/xc_split_partials_lc_wpbeh
  !  NAME
  !    xc_split_partials_lc_wpbeh
  !  SYNOPSIS

  subroutine xc_split_partials_lc_wpbeh &
       (lc_wpbeh_omega, rho, rho_gradient, &
        en_density_xc, en_density_x, en_density_c, &
        local_x_derivs, local_c_derivs, &
         x_density_gradient_deriv, c_density_gradient_deriv  )

    !  PURPOSE
    !  Subroutine xc_partials_lc_wpbeh provides the long-range corrected hybrid exchange-correlation energy functional
    !  base on wPBE and its partial
    !  derivatives by the density and the density gradient. Like the subroutine xc_partials_lc_wpbeh but with split output.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_x_derivs
    real*8, dimension(n_spin) :: local_c_derivs
    real*8, dimension(3,n_spin) :: x_density_gradient_deriv
    real*8, dimension(3,n_spin) :: c_density_gradient_deriv
    real*8 :: lc_wpbeh_omega


    !  INPUTS
    !   o lc_wpbeh_omega -- screening factor
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_x_derivs -- Partial derivative d/d(rho) of e_x
    !   o local_c_derivs -- Partial derivative d/d(rho) of e_c
    !   o x_density_gradient_deriv -- Partial derivative d/d(|grad(rho)|) of e_x
    !   o c_density_gradient_deriv -- Partial derivative d/d(|grad(rho)|) of e_c
    !
    !  AUTHOR
    !    Lukas Gallandi
    !  HISTORY
    !    ???
    !  SEE ALSO
    !    ???
    !  SOURCE
    !    FHI-aims team (subroutine 'xc_partials_hse') as well as other subroutines. Basically it's just
    !    a rearrangement of exisiting subroutines.








        !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)
    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho


    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    aux_squared_grad_rho=0.0

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       !if (rho(i_spin).gt.0.d0) then
       if (rho(i_spin).gt.1.d-20) then !Actual enforcing zero caused problem with restart of LC-wPBEh calculations


          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo
          !
          !           get exchange terms
          !
          !     NOTICE, that the output is allready :
          !
          !     (1/epsilon - hybrid_coef) * en_density_wpbe_shortrange + (1-1/epsilon)* en_density_pbe,
          !
          !     (1/epsilon - hybrid_coef)  * x_density_derive_wpbe_shortrange+ (1-1/epsilon)* x_density_derive_pbe
          !
          !         and
          !
          !     (1/epsilon - hybrid_coef)  * x_gradient_derive_wpbe_shortrange+ (1-1/epsilon)* x_gradient_derive_pbe,
          !
          call exch_lc_wpbeh_derivs &
               (hybrid_coeff,lc_wpbeh_omega,aux_density,aux_squared_grad_rho, &
               en_density_x(i_spin), local_x_derivs(i_spin), &
               x_gradient_deriv(i_spin),lc_dielectric_constant &
               )

       else

          en_density_x(i_spin) = 0.d0
          local_x_derivs(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then


       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    do i_spin = 1, n_spin, 1
       local_c_derivs(i_spin) = c_density_deriv 
       do i_coord = 1,3,1
          x_density_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin)
          c_density_gradient_deriv(i_coord, i_spin) = &
               c_gradient_deriv*total_rho_gradient(i_coord)
       enddo
    enddo

    if (spin_treatment.eq.1) then
       ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
       local_c_derivs(1) = c_density_deriv - &
            (spin_pol-1) * c_spin_deriv
       local_c_derivs(2) = c_density_deriv - &
            (spin_pol+1) * c_spin_deriv
    end if

    !  that's all folks

  end subroutine xc_split_partials_lc_wpbeh
  !******
  !----------------------------------------------------------------------
  !****s* xc/xc_partials_blyp
  !  NAME
  !    xc_partials_blyp
  !  SYNOPSIS

  subroutine xc_partials_blyp &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c,  local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_blyp provides the exchange-correlation energy functional
    !  of BLYP and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !  refs: A.D. Becke, J.Chem.Phys.96, 2155, 1992; Johnson, Gill and Pople - JCP 98. 5612 (1993),
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE





    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8, dimension(n_spin) :: c_density_deriv
    real*8, dimension(3,n_spin) :: c_gradient_deriv


    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo



          !           get exchange terms

          call xbecke_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )
          !            en_density_x(i_spin) = 0.d0
          !            x_density_deriv(i_spin) = 0.d0
          !            x_gradient_deriv(i_spin) = 0.d0
       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo


    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       !           spin_pol = ( rho(1) - rho(2) ) / full_rho
       !           squared_grad_rho = 0.d0
       !           do i_coord = 1, 3, 1
       !             total_rho_gradient(i_coord) =
       !     +         rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
       !             squared_grad_rho = squared_grad_rho +
       !     +       (total_rho_gradient(i_coord))**2.d0
       !           enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       !           spin_pol = 0.d0
       !           squared_grad_rho = aux_squared_grad_rho
       !           do i_coord = 1,3,1
       !             total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       !           enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if
    !       write(use_unit,*) "only exchange:", en_density_xc

    !       get correlation terms

    call lyp_part_derivs &
         (  rho, rho_gradient, en_density_c, &
         c_density_deriv, c_gradient_deriv, n_spin &
         )



    !       add relevant terms together

    en_density_xc = en_density_xc + en_density_c
    !      write(use_unit,*) "all terms:", en_density_xc

    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = x_density_deriv(i_spin) + &
            c_density_deriv(i_spin)
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin) * x_gradient_deriv(i_spin) * &
               rho_gradient(i_coord,i_spin) + &
               c_gradient_deriv(i_coord, i_spin)
       enddo
    enddo

    ! In this LYP formulation, the spin polarization is already taken implicitly into acount in the derivatives
    !        if (spin_treatment.eq.1) then
    !          ! Must add spin derivative terms here - check explicitly whether they are in corpbe, or not.
    !          local_xc_derivs(1) = local_xc_derivs(1) -
    !     +      (spin_pol-1) * c_spin_deriv
    !          local_xc_derivs(2) = local_xc_derivs(2) -
    !     +      (spin_pol+1) * c_spin_deriv
    !        end if

    !  that's all folks

  end subroutine xc_partials_blyp
  !******
  !---------------------------------------------------------------------------------------------
  !****s* xc/xc_partials_b3lyp
  !  NAME
  !    xc_partials_b3lyp
  !  SYNOPSIS

  subroutine xc_partials_b3lyp &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_b3lyp provides the exchange-correlation energy functional
    !  of B3LYP and its partial
    !  derivatives by the density and the density gradient. This uses the RPA parametrization
    !  as in the GAUSSIAN code.
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    ! refs: A.D. Becke, J.Chem.Phys.96, 2155, 1992; Johnson, Gill and Pople - JCP 98. 5612 (1993),
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho

    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x_becke
    real*8, dimension(n_spin) :: en_density_x_lda
    real*8, dimension(n_spin) :: x_density_deriv_becke
    real*8, dimension(n_spin) :: x_density_deriv_lda
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c_vwn, en_density_c_lyp
    real*8, dimension(n_spin) :: c_density_deriv_vwn
    real*8, dimension(n_spin) :: c_density_deriv_lyp
    real*8, dimension(3,n_spin) :: c_gradient_deriv




    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho
    !      real*8, dimension(3) :: dork


    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for B3LYP
    real*8, parameter :: ax=7.2d-1, ac=8.1d-1
    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo



          !           get exchange terms
          call xlda( aux_density, x_density_deriv_lda(i_spin), &
               en_density_x_lda(i_spin) )
          call xbecke_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x_becke(i_spin), &
               x_density_deriv_becke(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x_becke(i_spin) = 0.d0
          x_density_deriv_becke(i_spin) = 0.d0
          en_density_x_lda(i_spin) = 0.d0
          x_density_deriv_lda(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo



    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       !           spin_pol = ( rho(1) - rho(2) ) / full_rho

       do i_spin = 1,n_spin,1
       !   en_density_x(i_spin)= &
       !         ((1.d0-hybrid_coeff)*en_density_x_lda(i_spin)*rho(i_spin)+ &
       !        ax*(en_density_x_becke(i_spin)-en_density_x_lda(i_spin)) &
       !        *rho(i_spin)) / full_rho
          en_density_x(i_spin)= &
                (1.d0-hybrid_coeff)*en_density_x_lda(i_spin)+ &
               ax*(en_density_x_becke(i_spin)-en_density_x_lda(i_spin)) 
       enddo

       en_density_xc = (en_density_x(1)*rho(1)+ en_density_x(2)*rho(2)) / full_rho

    else
       en_density_x(1) =  &
            (1.d0-hybrid_coeff)*en_density_x_lda(1)+ &
            ax*(en_density_x_becke(1)-en_density_x_lda(1))
       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0

    end if


    !       get correlation terms

    call lyp_part_derivs &
         (  rho, rho_gradient, en_density_c_lyp, &
         c_density_deriv_lyp, c_gradient_deriv, n_spin &
         )

    call vwnrpa_c_derivs &
         (rho, en_density_c_vwn, c_density_deriv_vwn, n_spin)





    !       add relevant terms together

    en_density_c =  (1.d0-ac)*en_density_c_vwn + &
         ac*en_density_c_lyp
    en_density_xc = en_density_xc + en_density_c


    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin)=(1.d0-hybrid_coeff)* &
            x_density_deriv_lda(i_spin)+ &
            ax*(x_density_deriv_becke(i_spin)-x_density_deriv_lda(i_spin))+ &
            (1.d0-ac)*c_density_deriv_vwn(i_spin)+ &
            ac*c_density_deriv_lyp(i_spin)

       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin)*ax*x_gradient_deriv(i_spin)* &
               rho_gradient(i_coord,i_spin) + &
               ac*c_gradient_deriv(i_coord, i_spin)
       enddo
    enddo



    !  that's all folks

  end subroutine xc_partials_b3lyp
  !******
  !---------------------------------------------------------------------------------------------
  !****s* xc/xc_partials_b1lyp
  !  NAME
  !    xc_partials_b1lyp
  !  SYNOPSIS

  subroutine xc_partials_b1lyp &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_b1lyp provides the exchange-correlation energy functional
    !  of B1LYP and its partial
    !  derivatives by the density and the density gradient. 
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  Exc[B1LYP] = Exc[BLYP] + 0.25*(Ex[HF]-Ex[Becke88])
    !
    ! refs: A.D. Becke, J.Chem.Phys.104, 1040, 1996;
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho

    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x_becke
    real*8, dimension(n_spin) :: en_density_x_lda
    real*8, dimension(n_spin) :: x_density_deriv_becke
    real*8, dimension(n_spin) :: x_density_deriv_lda
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c_vwn, en_density_c_lyp
    real*8, dimension(n_spin) :: c_density_deriv_vwn
    real*8, dimension(n_spin) :: c_density_deriv_lyp
    real*8, dimension(3,n_spin) :: c_gradient_deriv




    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho
    !      real*8, dimension(3) :: dork


    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for B1LYP
    real*8, parameter :: ax=0.75, ac=1.0
    !real*8, parameter :: ax=7.2d-1, ac=8.1d-1
    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo



          !           get exchange terms
          call xlda( aux_density, x_density_deriv_lda(i_spin), &
               en_density_x_lda(i_spin) )
          call xbecke_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x_becke(i_spin), &
               x_density_deriv_becke(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x_becke(i_spin) = 0.d0
          x_density_deriv_becke(i_spin) = 0.d0
          en_density_x_lda(i_spin) = 0.d0
          x_density_deriv_lda(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo



    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       !           spin_pol = ( rho(1) - rho(2) ) / full_rho

       do i_spin = 1,n_spin,1
          en_density_x(i_spin)= &
                (1.d0-hybrid_coeff)*en_density_x_lda(i_spin)+ &
               ax*(en_density_x_becke(i_spin)-en_density_x_lda(i_spin)) 
          !en_density_x(i_spin)= &
          !      (1.d0-hybrid_coeff)*en_density_x_becke(i_spin)
       enddo

       en_density_xc = (en_density_x(1)*rho(1)+ en_density_x(2)*rho(2)) / full_rho

    else
       en_density_x(1) =  &
            (1.d0-hybrid_coeff)*en_density_x_becke(1)
       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0

    end if


    !       get correlation terms

    call lyp_part_derivs &
         (  rho, rho_gradient, en_density_c_lyp, &
         c_density_deriv_lyp, c_gradient_deriv, n_spin &
         )

    call vwnrpa_c_derivs &
         (rho, en_density_c_vwn, c_density_deriv_vwn, n_spin)





    !!       add relevant terms together

    !en_density_c =  en_density_c_lyp
    !en_density_xc = en_density_xc + en_density_c


    !do i_spin = 1, n_spin, 1
    !   local_xc_derivs(i_spin)=(1.d0-hybrid_coeff)* &
    !        x_density_deriv_becke(i_spin)+ &
    !        c_density_deriv_lyp(i_spin)

    !   do i_coord = 1,3,1
    !      xc_gradient_deriv(i_coord, i_spin) = &
    !           dble(n_spin)*x_gradient_deriv(i_spin)* &
    !           rho_gradient(i_coord,i_spin) + &
    !           c_gradient_deriv(i_coord, i_spin)
    !   enddo
    !enddo

    !       add relevant terms together

    en_density_c =  (1.d0-ac)*en_density_c_vwn + &
         ac*en_density_c_lyp
    en_density_xc = en_density_xc + en_density_c


    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin)=(1.d0-hybrid_coeff)* &
            x_density_deriv_lda(i_spin)+ &
            ax*(x_density_deriv_becke(i_spin)-x_density_deriv_lda(i_spin))+ &
            (1.d0-ac)*c_density_deriv_vwn(i_spin)+ &
            ac*c_density_deriv_lyp(i_spin)

       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               dble(n_spin)*ax*x_gradient_deriv(i_spin)* &
               rho_gradient(i_coord,i_spin) + &
               ac*c_gradient_deriv(i_coord, i_spin)
       enddo
    enddo




    !  that's all folks

  end subroutine xc_partials_b1lyp
  !******
  !---------------------------------------------------------------------------------------------
  !****s* xc/xc_en_pw91_lda
  !  NAME
  !    xc_en_pw91_lda
  !  SYNOPSIS

  subroutine xc_en_pw91_lda &
       ( density, en_density_ldax, en_density_ldac &
       )

    !  PURPOSE
    !  Subroutine xcpot_pw91_lda provides the exchange-correlation potential
    !  according to the Perdew-Wang 1991 parametrisation of Ceperley-Alder LDA
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8 density(n_spin)

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !
    !  OUTPUT
    !   o en_density_ldax -- ???
    !   o en_density_ldac -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE




    real*8 :: en_density_ldax
    real*8 :: en_density_ldac


    !  local variables

    !       aux_density : Extra variable to specifiy the actual input density which enters each subroutine

    !       en_density_c : Correlation energy density at current grid point
    !       pot_c :        correlation potential at current grid point (spin-up and spin-down)
    !       inv_k_fermi :  Inverse of the Fermi wave vector of a homogeneous electron gas of given density
    !       r_seitz:       Wigner-Seitz radius
    !       zeta :         would be spin polarisation
    !       ec_rs :        derivative of en_density_c w.r.t. r_seitz
    !       ec_zeta :      derivative of en_density_c w.r.t. zeta
    !       alpha_c :      "correlation contribution to spin stiffness"

    real*8 :: aux_density
    real*8, dimension(n_spin) :: en_density_x
    real*8, dimension(n_spin) :: pot_xc
    real*8, dimension(n_spin) :: pot_c

    real*8 :: en_density_c
    real*8 :: inv_k_f
    real*8 :: r_seitz
    real*8 :: zeta = 0.d0
    real*8 :: ec_rs = 0.d0
    real*8 :: ec_zeta = 0.d0
    real*8 :: alpha_c = 0.d0

    !       counters

    integer :: i_spin

    !  begin work

    ! Exchange energy first, for both spin channels
    do i_spin = 1, n_spin, 1
       aux_density = dble(n_spin) * density(i_spin)
       call xlda( aux_density, pot_xc(i_spin), &
            en_density_x(i_spin) )
    enddo

    if (spin_treatment.eq.1) then

       aux_density = density(1)+density(2)
       zeta = 0.d0

       en_density_ldax = 0.d0
       if(aux_density .gt. 1.e-16) then
          do i_spin = 1,n_spin,1
             en_density_ldax = en_density_ldax &
                  + en_density_x(i_spin)*density(i_spin)
          enddo
          en_density_ldax = en_density_ldax / aux_density
          zeta = ( density(1)-density(2) ) / aux_density
       endif

    else

       en_density_ldax = en_density_x(1)
       aux_density = density(1)
       zeta = 0.d0

    end if

    ! Correlation energy is next

    inv_k_f = (pisq3*aux_density)**third
    r_seitz = const_rs/inv_k_f

    if (spin_treatment.eq.0) then
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c(1), &
            pot_c(1), ec_rs, ec_zeta, alpha_c )
    else
       call corlsd &
            ( r_seitz, zeta, en_density_c, pot_c(1), &
            pot_c(2), ec_rs, ec_zeta, alpha_c )
    end if

    en_density_ldac = en_density_c


  end subroutine xc_en_pw91_lda
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_am05
  !  NAME
  !    xc_partials_am05
  !  SYNOPSIS

  subroutine xc_partials_am05 &
       ( rho, rho_gradient, en_density_xc, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_am05 provides the exchange-correlation energy functional
    !  of Armiento and Mattson, PRB 72, 085108 (2005) and its partial
    !  derivatives by the density and the density gradient.
    !
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    !
    !  This routine does not include spin polarization since, for this functional
    !  the spin polarized "case" was not yet released.
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !   
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    !  Local variables

    real*8, dimension(n_spin) :: density_deriv
    real*8, dimension(n_spin) :: gradient_deriv

    real*8, dimension(2) :: aux_rho
    real*8, dimension(2) :: aux_abs_rho_gradient
    real*8, dimension(2) :: aux_local_xc_derivs, aux_gradient_deriv 
    real*8 :: aux_gradient_squared
    real*8 :: aux_rho_gradient

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work

    !       line up variables for use in partial derivatives:

! Re-shape quantities to input in the subroutine

    if (n_spin.eq.1) then
        do i_spin = 1, 2, 1
        aux_rho(i_spin) = rho(1)/2.d0
        aux_gradient_squared = 0.0
          do i_coord = 1, 3, 1
             aux_rho_gradient = rho_gradient(i_coord, 1)
             aux_gradient_squared= aux_gradient_squared + &
                         aux_rho_gradient**2.d0
          enddo
        aux_abs_rho_gradient(i_spin) =  sqrt(aux_gradient_squared)
        enddo
    else
        do i_spin = 1, 2, 1
        aux_rho(i_spin) = rho(i_spin)
        aux_gradient_squared = 0.0
          do i_coord = 1, 3, 1
             aux_rho_gradient = 2.d0*rho_gradient(i_coord, i_spin)
             aux_gradient_squared= aux_gradient_squared + &
                         aux_rho_gradient**2.d0
          enddo
        aux_abs_rho_gradient(i_spin) =  sqrt(aux_gradient_squared) 
        enddo
    endif

    call am05_partial_derivs(aux_rho(1), aux_rho(2), aux_abs_rho_gradient(1), aux_abs_rho_gradient(2), &
                             en_density_xc, aux_local_xc_derivs(1), aux_local_xc_derivs(2), & 
                             aux_gradient_deriv(1), aux_gradient_deriv(2))

! Put quantities back in their right variables 
    do i_spin = 1, n_spin, 1
       local_xc_derivs(i_spin) = aux_local_xc_derivs(i_spin)
       gradient_deriv(i_spin) = aux_gradient_deriv(i_spin) 
       do i_coord = 1,3,1
          xc_gradient_deriv(i_coord, i_spin) = &
               gradient_deriv(i_spin) * &
               dble(n_spin)*rho_gradient(i_coord,i_spin)
       enddo
    enddo

  end subroutine xc_partials_am05
!--------------------------------------------------------------------------
subroutine xc_partials_m06_family &
           ( rho, rho_gradient, kinetic_density, en_density_xc, en_density_x, en_density_c, &
             local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, m06_option, &
             tau_threshold )

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8, dimension(n_spin) :: xc_tau_deriv
    real*8, dimension(n_spin) :: kinetic_density
    integer, intent(in) :: m06_option
    real*8, intent(in)  :: tau_threshold

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !   o kinetic_density -- Gradient of the orbitals
    !   o m06_option -- Choice for m06 functionals. The options are:
    !                   1 : M06-L
    !                   2 : M06-HF (100% HF Exchange)
    !                   3 : M06    (27% HF Exchange)
    !                   4 : M06-2X (54% HF Exchange)
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho, tau )
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, tau )
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_gradient_deriv -- Partial derivative d/d(|grad(rho)|^2) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_tau_deriv -- Partial derivative d/d(tau) of of e_xc(rho,|grad(rho)|^2,tau)
    !
    !  AUTHOR
    !    FHI-aims team.
    !    Updated by Andrew Logsdail, Jan 2015
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    !  Local variables

    real*8, dimension(7) :: derivatives_x
    real*8, dimension(7) :: derivatives_c
    real*8 :: total_rho
    real*8 :: total_tau
    real*8 :: threshold

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    en_density_xc     = 0.d0
    en_density_x      = 0.d0
    en_density_c      = 0.d0
    local_xc_derivs   = 0.d0
    xc_gradient_deriv = 0.d0
    xc_tau_deriv      = 0.0d0

    derivatives_x     = 0.d0
    derivatives_c     = 0.d0
    threshold         = 1.d-10
    total_rho         = 0.d0
    total_tau         = 0.d0

    do i_spin=1, n_spin, 1
     total_rho = total_rho + rho(i_spin)
     total_tau = total_tau + kinetic_density(i_spin)
    enddo

    ! exchange and correlation

    if (total_rho.gt.threshold.and.total_tau.gt.tau_threshold) then

      if (spin_treatment.eq.0) then

        ! Note for self, and anyone else debugging this.
        ! For spin-paired calculations, this implementation seems to vary very slightly in what it produces
        ! depending on the number of processors we are parallelising across.
        !
        ! However, convergence is always to the correct final number. I have tried to debug by zeroing
        ! everything, but cannot find the source. As the end of the calculation is correct, I'm not going to
        ! trouble myself any more with it.
        !
        ! The spin-polarised implementation returns identical results irrespective of CPUs, 
        ! according to my testing. Same goes for TPSS spin-paired and -polarised. AJL

        call M06x(en_density_x,derivatives_x(:),rho(1)*0.5d0,rho(1)*0.5d0,rho_gradient(:,1)*0.5d0,rho_gradient(:,1)*0.5d0, &
                  kinetic_density(1)*0.5d0,kinetic_density(1)*0.5d0,m06_option)

        call M06c(en_density_c,derivatives_c(:),rho(1)*0.5d0,rho(1)*0.5d0,rho_gradient(:,1)*0.5d0,rho_gradient(:,1)*0.5d0, & 
                  kinetic_density(1)*0.5d0,kinetic_density(1)*0.5d0,m06_option) 

      else ! (spin_treatment.eq.1)

        call M06x(en_density_x,derivatives_x(:),rho(1),rho(2),rho_gradient(:,1),rho_gradient(:,2),& 
                  kinetic_density(1),kinetic_density(2),m06_option)

        call M06c(en_density_c,derivatives_c(:),rho(1),rho(2),rho_gradient(:,1),rho_gradient(:,2),& 
                  kinetic_density(1),kinetic_density(2),m06_option) 

      endif        

! The parameters of the M06 functionals are already optimized and
! set WITH the hybrid coefficient included. Therefore, there is no
! need to multiply the exchange part bu the coefficient.
      !en_density_xc = ((1.d0-hybrid_coeff)*en_density_x + en_density_c)/total_rho
      !en_density_x  = ((1.d0-hybrid_coeff)*en_density_x)/total_rho

      en_density_xc = (en_density_x + en_density_c)/total_rho
      en_density_x  = en_density_x/total_rho
      en_density_c  = en_density_c/total_rho

! We'll use the potential arrays as they are. AJL
! Indices for derivative_x and derivative_c:
! 1 - dF/d(rho_alpha)
! 2 - dF/d(rho_beta)
! 3 - dF/d(gamma_alpha.alpha)
! 4 - dF/d(gamma_beta.beta)
! 5 - dF/d(gamma_alpha.beta)
! 6 - dF/d(tau_alpha)
! 7 - dF/d(tau_beta)

      ! Scale for Hybrid Calculations
      ! derivatives_x(:) = (1.d0-hybrid_coeff)*derivatives_x(:)

      do i_spin = 1, n_spin, 1

         ! Derivatives with respect to rho 
         local_xc_derivs(i_spin) = 1.0d0 * (derivatives_x(i_spin) + derivatives_c(i_spin))

         ! Derivatives with respect to gamma
         do i_coord = 1,3,1
            ! There is no contribution from gamma_alpha.beta
            xc_gradient_deriv(i_coord, i_spin) = &
                   0.5d0 * rho_gradient(i_coord,i_spin) * dble(n_spin) * &
                   (derivatives_x(i_spin+2) + derivatives_c(i_spin+2))
         enddo

         ! Derivatives with respect to tau
         xc_tau_deriv(i_spin) = 1.0d0 * (derivatives_x(i_spin+5) + derivatives_c(i_spin+5))

      enddo

    endif ! total_rho.gt.threshold
 
end subroutine xc_partials_m06_family

subroutine xc_partials_m08m11_family &
           ( rho, rho_gradient, kinetic_density, en_density_xc, en_density_x, en_density_c, &
             local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, m08m11_option, &
             tau_threshold )

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8, dimension(n_spin) :: xc_tau_deriv
    real*8, dimension(n_spin) :: kinetic_density
    integer, intent(in) :: m08m11_option
    real*8, intent(in) :: tau_threshold

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !   o kinetic_density -- Gradient of the orbitals
    !   o m08m11_option -- Choice for m08 and m11 functionals. The options are:
    !                   1 : M08-HX (52.23% Exact Exchange)
    !                   2 : M08-SO (56.79% Exact Exchange)
    !                   3 : M11    (42.80% Exact Exchange SR; 100% Exact Exchange LR)
    !                   4 : M11-L  (00.00% Exact Exchange)
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho, tau )
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, tau )
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_gradient_deriv -- Partial derivative d/d(|grad(rho)|^2) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_tau_deriv -- Partial derivative d/d(tau) of of e_xc(rho,|grad(rho)|^2,tau)
    !
    !  AUTHOR
    !    FHI-aims team.
    !    Updated by Andrew Logsdail, Jan 2015
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    !  Local variables

    real*8, dimension(7) :: derivatives_x
    real*8, dimension(7) :: derivatives_c
    real*8 :: total_rho
    real*8 :: total_tau
    real*8 :: threshold

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    en_density_xc     = 0.d0
    en_density_x      = 0.d0
    en_density_c      = 0.d0
    local_xc_derivs   = 0.d0
    xc_gradient_deriv = 0.d0
    xc_tau_deriv      = 0.0d0

    derivatives_x     = 0.d0
    derivatives_c     = 0.d0
    threshold         = 1.d-10
    total_rho         = 0.d0
    total_tau         = 0.d0

    do i_spin=1, n_spin, 1
     total_rho = total_rho + rho(i_spin)
     total_tau = total_tau + kinetic_density(i_spin)
    enddo

    ! exchange and correlation

!    if (total_rho.gt.threshold) then
    if (total_rho.gt.threshold.and.total_tau.gt.tau_threshold) then

      if (spin_treatment.eq.0) then

        ! Note for self, and anyone else debugging this.
        ! For spin-paired calculations, this implementation seems to vary very slightly in what it produces
        ! depending on the number of processors we are parallelising across.
        !
        ! However, convergence is always to the correct final number. I have tried to debug by zeroing
        ! everything, but cannot find the source. As the end of the calculation is correct, I'm not going to
        ! trouble myself any more with it.
        !
        ! The spin-polarised implementation has the same problem, according to my tests. However, as I wrote above,
        ! everything has been benchmarked against other codes (NWChem) and the literature, and the final results
        ! are always correct.

        call M08M11x(en_density_x,derivatives_x(:),rho(1)*0.5d0,rho(1)*0.5d0,rho_gradient(:,1)*0.5d0,rho_gradient(:,1)*0.5d0, &
                  kinetic_density(1)*0.5d0,kinetic_density(1)*0.5d0,m08m11_option)

        call M08M11c(en_density_c,derivatives_c(:),rho(1)*0.5d0,rho(1)*0.5d0,rho_gradient(:,1)*0.5d0,rho_gradient(:,1)*0.5d0, &
                  kinetic_density(1)*0.5d0,kinetic_density(1)*0.5d0,m08m11_option)

      else ! (spin_treatment.eq.1)

        call M08M11x(en_density_x,derivatives_x(:),rho(1),rho(2),rho_gradient(:,1),rho_gradient(:,2),&
                  kinetic_density(1),kinetic_density(2),m08m11_option)

        call M08M11c(en_density_c,derivatives_c(:),rho(1),rho(2),rho_gradient(:,1),rho_gradient(:,2),&
                  kinetic_density(1),kinetic_density(2),m08m11_option)

      endif

! The parameters of the M08/M11 functionals are already optimized and
! set WITH the hybrid coefficient included. Therefore, there is no
! need to multiply the exchange part bu the coefficient.
      !en_density_xc = ((1.d0-hybrid_coeff)*en_density_x + en_density_c)/total_rho
      !en_density_x  = ((1.d0-hybrid_coeff)*en_density_x)/total_rho

      en_density_xc = (en_density_x + en_density_c)/total_rho
      en_density_x  = en_density_x/total_rho
      en_density_c  = en_density_c/total_rho

! We'll use the potential arrays as they are. AJL
! Indices for derivative_x and derivative_c:
! 1 - dF/d(rho_alpha)
! 2 - dF/d(rho_beta)
! 3 - dF/d(gamma_alpha.alpha)
! 4 - dF/d(gamma_beta.beta)
! 5 - dF/d(gamma_alpha.beta)
! 6 - dF/d(tau_alpha)
! 7 - dF/d(tau_beta)

      ! Scale for Hybrid Calculations
      ! derivatives_x(:) = (1.d0-hybrid_coeff)*derivatives_x(:)

      do i_spin = 1, n_spin, 1

         ! Derivatives with respect to rho
         local_xc_derivs(i_spin) = 1.0d0 * (derivatives_x(i_spin) + derivatives_c(i_spin))

         ! Derivatives with respect to gamma
         do i_coord = 1,3,1
            ! There is a gamma_alpha.beta contribution
            if (spin_treatment.eq.0) then
               ! Remove dble(n_spin) from these sums as it is always 1 in this case.
               xc_gradient_deriv(i_coord, i_spin) = &
                    0.5d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(3) + derivatives_c(3)) + &
                    0.25d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(5) + derivatives_c(5))
            else
               ! Remove dble(n_spin) from these sums as it is always 2 in this case.
               xc_gradient_deriv(i_coord, i_spin) = &
                    1.0d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(i_spin+2) + derivatives_c(i_spin+2)) + &
                    0.5d0 * rho_gradient(i_coord,3-i_spin) * &
                    (derivatives_x(5) + derivatives_c(5))
            endif
         enddo

         ! Derivatives with respect to tau
         xc_tau_deriv(i_spin) = 1.0d0 * (derivatives_x(i_spin+5) + derivatives_c(i_spin+5))

      enddo

    endif ! total_rho.gt.threshold

end subroutine xc_partials_m08m11_family

  !---------------------------------------------------------------------------------------------
  !****s* xc/xc_partials_xyg3_dft_part
  !  NAME
  !    xc_partials_xyg3_dft_part
  !  SYNOPSIS

  subroutine xc_partials_xyg3_dft_part &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_xyg3_dft_part provides the DFT part in XYG3
    !  These are the XC components which are later needed to compute the
    !  XYG3 total energy
    ! refs: Y. Zhang, X. Xu, W. A. Gorddard, PNAS, 106, 4963 (2009).
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho

    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x_becke
    real*8, dimension(n_spin) :: en_density_x_lda
    real*8, dimension(n_spin) :: x_density_deriv_becke
    real*8, dimension(n_spin) :: x_density_deriv_lda
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c_vwn, en_density_c_lyp
    real*8, dimension(n_spin) :: c_density_deriv_vwn
    real*8, dimension(n_spin) :: c_density_deriv_lyp
    real*8, dimension(3,n_spin) :: c_gradient_deriv




    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho
    !      real*8, dimension(3) :: dork


    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for XYG3
    !  Exc[XYG3-DFT-Part] = ax * Ex[lda] + adx * (Ex[B88]-Ex[lda]) + ac * Ec[LYP]
    !  Y. Zhang, X. Xu, W. A. Gorddard, PNAS, 2009, 106, 4963
    real*8, parameter :: ax=1.967d-1, adx=2.107d-1, ac=6.789d-1
    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo



          !           get exchange terms
          call xlda( aux_density, x_density_deriv_lda(i_spin), &
               en_density_x_lda(i_spin) )
          call xbecke_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x_becke(i_spin), &
               x_density_deriv_becke(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x_becke(i_spin) = 0.d0
          x_density_deriv_becke(i_spin) = 0.d0
          en_density_x_lda(i_spin) = 0.d0
          x_density_deriv_lda(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo


    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       !           spin_pol = ( rho(1) - rho(2) ) / full_rho

       do i_spin = 1,n_spin,1
       !   en_density_x(i_spin)= &
       !         ((1.d0-hybrid_coeff)*en_density_x_lda(i_spin)*rho(i_spin)+ &
       !        ax*(en_density_x_becke(i_spin)-en_density_x_lda(i_spin)) &
       !        *rho(i_spin)) / full_rho
          en_density_x(i_spin)= &
                ax*en_density_x_lda(i_spin)+ &
               adx*(en_density_x_becke(i_spin)-en_density_x_lda(i_spin)) 
       enddo

       en_density_xc = (en_density_x(1)*rho(1)+ en_density_x(2)*rho(2)) / full_rho

    else
       en_density_x(1) =  &
            ax*en_density_x_lda(1)+ &
            adx*(en_density_x_becke(1)-en_density_x_lda(1))
       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0

    end if


    !       get correlation terms

    call lyp_part_derivs &
         (  rho, rho_gradient, en_density_c_lyp, &
         c_density_deriv_lyp, c_gradient_deriv, n_spin &
         )

    !call vwnrpa_c_derivs &
    !     (rho, en_density_c_vwn, c_density_deriv_vwn, n_spin)


    !  add relevant terms together
    en_density_c  = ac*en_density_c_lyp
    en_density_xc = en_density_xc + en_density_c

    !  that's all folks

  end subroutine xc_partials_xyg3_dft_part

  !---------------------------------------------------------------------------------------------
  !****s* xc/xc_partials_xygjos_dft_part
  !  NAME
  !    xc_partials_xygjos_dft_part
  !  SYNOPSIS

  subroutine xc_partials_xygjos_dft_part &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_xygjos_dft_part provides the DFT-part energy in XYGJOS
    !  These are the XC components which are later needed to compute the
    !  XYGJOS total energy.
    ! refs: Y. Zhang, X. Xu, Y. Jung, W. A. Gorddard, PNAS, 108, 19897 (2011).
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho

    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x_becke
    real*8, dimension(n_spin) :: en_density_x_lda
    real*8, dimension(n_spin) :: x_density_deriv_becke
    real*8, dimension(n_spin) :: x_density_deriv_lda
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c_vwn, en_density_c_lyp
    real*8, dimension(n_spin) :: c_density_deriv_vwn
    real*8, dimension(n_spin) :: c_density_deriv_lyp
    real*8, dimension(3,n_spin) :: c_gradient_deriv




    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho
    !      real*8, dimension(3) :: dork


    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for XYGJOS
    !  Exc[XYGJOS-DFT-Part] = ax * Ex[LDA] + ab * Ec[VWN] + ac * Ec[LYP]
    !  Y. Zhang, X. Xu, Y. Jung, W. A. Gorddard, PNAS, 108, 19897 (2011).
    real*8, parameter :: ax=2.269d-1, ab=2.309d-1, ac=2.754d-1
    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo



          !           get exchange terms
          call xlda( aux_density, x_density_deriv_lda(i_spin), &
               en_density_x_lda(i_spin) )
          !call xbecke_derivs &
          !     (  aux_density, aux_squared_grad_rho, &
          !     en_density_x_becke(i_spin), &
          !     x_density_deriv_becke(i_spin), &
          !     x_gradient_deriv(i_spin) &
          !     )

       else

          !en_density_x_becke(i_spin) = 0.d0
          !x_density_deriv_becke(i_spin) = 0.d0
          en_density_x_lda(i_spin) = 0.d0
          x_density_deriv_lda(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo



    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       !           spin_pol = ( rho(1) - rho(2) ) / full_rho

       do i_spin = 1,n_spin,1
       !   en_density_x(i_spin)= &
       !         ax*en_density_x_lda(i_spin)+ &
       !        adx*(en_density_x_becke(i_spin)-en_density_x_lda(i_spin)) 
          en_density_x(i_spin)= &
                ax*en_density_x_lda(i_spin)
       enddo

       en_density_xc = (en_density_x(1)*rho(1)+ en_density_x(2)*rho(2)) / full_rho

    else
       !en_density_x(1) =  &
       !     ax*en_density_x_lda(1)+ &
       !     adx*(en_density_x_becke(1)-en_density_x_lda(1))
       en_density_x(1) =  &
            ax*en_density_x_lda(1)
       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0

    end if


    !       get correlation terms

    call lyp_part_derivs &
         (  rho, rho_gradient, en_density_c_lyp, &
         c_density_deriv_lyp, c_gradient_deriv, n_spin &
         )

    call vwnrpa_c_derivs &
         (rho, en_density_c_vwn, c_density_deriv_vwn, n_spin)

    !       add relevant terms together

    en_density_c  = ab*en_density_c_vwn+ac*en_density_c_lyp
    en_density_xc = en_density_xc + en_density_c

    !  that's all folks

  end subroutine xc_partials_xygjos_dft_part

  !---------------------------------------------------------------------------------------------
  !****s* xc/xc_partials_xyg5_dft_part
  !  NAME
  !    xc_partials_xyg5_dft_part
  !  SYNOPSIS

  subroutine xc_partials_xyg5_dft_part &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
       xc_gradient_deriv &
       )

    !  PURPOSE
    !  Subroutine xc_partials_xygjos_dft_part provides the DFT-part energy in XYG5
    !  These are the XC components which are later needed to compute the
    !  XYG5 total energy.
    ! refs: Y. Zhang, X. Xu, in writting.
    !
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS    

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho

    real*8 :: spin_pol
    real*8, dimension(n_spin) :: en_density_x_becke
    real*8, dimension(n_spin) :: en_density_x_lda
    real*8, dimension(n_spin) :: x_density_deriv_becke
    real*8, dimension(n_spin) :: x_density_deriv_lda
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: en_density_c_vwn, en_density_c_lyp
    real*8, dimension(n_spin) :: c_density_deriv_vwn
    real*8, dimension(n_spin) :: c_density_deriv_lyp
    real*8, dimension(3,n_spin) :: c_gradient_deriv




    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho
    !      real*8, dimension(3) :: dork


    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for XYG5
    !  Exc[XYG5-DFT-Part] = ax*Ex[LDA] + ab*Ex[Beck88] + ac*Ec[VWN] + ad*Ec[LYP]
    real*8, parameter ::    ax=5.000d-2, ab=7.200d-2,    ac=1.100d-1, ad=4.000d-1
    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo



          !           get exchange terms
          call xlda( aux_density, x_density_deriv_lda(i_spin), &
               en_density_x_lda(i_spin) )
          call xbecke_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x_becke(i_spin), &
               x_density_deriv_becke(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x_becke(i_spin) = 0.d0
          x_density_deriv_becke(i_spin) = 0.d0
          en_density_x_lda(i_spin) = 0.d0
          x_density_deriv_lda(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo



    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       !           spin_pol = ( rho(1) - rho(2) ) / full_rho

       do i_spin = 1,n_spin,1
          en_density_x(i_spin)= &
                ax*en_density_x_lda(i_spin) + &
                ab*en_density_x_becke(i_spin)
       enddo

       en_density_xc = (en_density_x(1)*rho(1)+ en_density_x(2)*rho(2)) / full_rho

    else
       en_density_x(1) =  &
            ax*en_density_x_lda(1) + ab*en_density_x_becke(1)
       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0

    end if


    !       get correlation terms

    call lyp_part_derivs &
         (  rho, rho_gradient, en_density_c_lyp, &
         c_density_deriv_lyp, c_gradient_deriv, n_spin &
         )

    call vwnrpa_c_derivs &
         (rho, en_density_c_vwn, c_density_deriv_vwn, n_spin)

    !       add relevant terms together

    en_density_c  = ac*en_density_c_vwn+ad*en_density_c_lyp
    en_density_xc = en_density_xc + en_density_c

    !  that's all folks

  end subroutine xc_partials_xyg5_dft_part

  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_xdhpbe0_dft_part
  !  NAME
  !    xc_partials_pbe0
  !  SYNOPSIS

  subroutine xc_partials_xdhpbe0_dft_part &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
         xc_gradient_deriv &
       )

    !  PURPOSE
    !  PURPOSE
    !  Subroutine xc_partials_xdhpbe0_dft_part provides the DFT part in xDH-PBE0
    !  These are the XC components which are later needed to compute the
    !  Hamiltonian matrix elements by integration by parts; the actual XC potential
    !  is NOT calculated here as we would need the density Hessian (second
    !  derivatives) for that task.
    ! refs: I.Y. Zhang et. al. JCP, 136 174103 (2012)
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE


    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for xDH-PBE0
    !  Exc[xDH-PBE0-DFT-Part] = ax * Ex[PBE] + ac * Ec[PBE]
    !  I.Y. Zhang, et. al., JPC, 136, 174103 (2012)
    real*8, parameter :: ax=1.665d-1, ac=5.292d-1

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    ! add relevant terms together
    en_density_xc = ax * en_density_xc &
         + ac * en_density_c

    do i_spin = 1, n_spin, 1
     en_density_x(i_spin)= ax * en_density_x(i_spin)
    enddo

    !  that's all folks

  end subroutine xc_partials_xdhpbe0_dft_part
  !******
  !---------------------------------------------------------------------
  !****s* xc/xc_partials_zrs1_dft_part
  !  NAME
  !    xc_partials_pbe0
  !  SYNOPSIS

  subroutine xc_partials_zrs1_dft_part &
       ( rho, rho_gradient, en_density_xc, &
         en_density_x, en_density_c, local_xc_derivs, &
         xc_gradient_deriv &
       )

    !  PURPOSE
    !  PURPOSE
    !  Subroutine xc_partials_zrs1_dft_part provides the DFT part in ZRPS
    !  These are the XC components which are later needed to compute the
    !  ZRPS total energy
    ! refs: I.Y. Zhang et. al. (to be submitted, 2015) 
    !
    !  USES

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !   o  local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho, |grad(rho)|^2)
    !   o xc_gradient_deriv -- ???
    !
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    !  Local variables

    !       squared_grad_rho = |grad(rho)|^2 enters our PBE subroutines
    !       spin_pol is a placeholder for when we actually introduce
    !                spin-polarization; spin_pol = zeta := (rho_up-rho_d)/rho
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho )
    !       x_density_deriv  : partial derivative d(rho*ex)/d(rho)
    !       x_gradient_deriv : partial derivative d(rho*ex)/d(|grad(rho)|^2)
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, zeta )
    !       c_density_deriv  : partial derivative d(rho*ec)/d(rho)
    !       c_gradient_deriv : partial derivative d(rho*ec)/d(|grad(rho)|^2)
    !       c_spin_deriv     : partial derivative (1/rho)*d(rho*ec)/d(zeta)

    real*8 :: full_rho
    real*8 :: squared_grad_rho
    real*8 :: spin_pol
    real*8, dimension(n_spin) :: x_density_deriv
    real*8, dimension(n_spin) :: x_gradient_deriv
    real*8 :: c_density_deriv
    real*8 :: c_gradient_deriv
    real*8 :: c_spin_deriv

    real*8, dimension(3) :: total_rho_gradient

    real*8 :: aux_density
    real*8 :: aux_squared_grad_rho

    !       counters

    integer :: i_coord
    integer :: i_spin
    !  parameters for xDH-PBE0
    !  Exc[ZRS1-DFT-Part] = ax * Ex[PBE] + ac * Ec[PBE]
    !  I.Y. Zhang, et. al., (to be submitted, 2015)
    real*8, parameter :: ax=5.000d-1, ac=7.500d-1

    !  begin work

    !       line up variables for use in partial derivatives:

    ! exchange first
    do i_spin = 1, n_spin, 1

       if (rho(i_spin).gt.0.d0) then

          aux_density = dble(n_spin)*rho(i_spin)

          aux_squared_grad_rho = 0.d0
          do i_coord = 1, 3, 1
             aux_squared_grad_rho = aux_squared_grad_rho + &
                  (dble(n_spin)*rho_gradient(i_coord,i_spin))**2.d0
          enddo

          !           get exchange terms

          call exchpbe_derivs &
               (  aux_density, aux_squared_grad_rho, &
               en_density_x(i_spin), x_density_deriv(i_spin), &
               x_gradient_deriv(i_spin) &
               )

       else

          en_density_x(i_spin) = 0.d0
          x_density_deriv(i_spin) = 0.d0
          x_gradient_deriv(i_spin) = 0.d0

       end if

    enddo

    if (spin_treatment.eq.1) then

       full_rho = rho(1) + rho(2)
       spin_pol = ( rho(1) - rho(2) ) / full_rho
       squared_grad_rho = 0.d0
       do i_coord = 1, 3, 1
          total_rho_gradient(i_coord) = &
               rho_gradient(i_coord,1)+rho_gradient(i_coord,2)
          squared_grad_rho = squared_grad_rho + &
               (total_rho_gradient(i_coord))**2.d0
       enddo

       en_density_xc = 0.d0
       do i_spin = 1,n_spin,1
          en_density_xc = en_density_xc &
               + en_density_x(i_spin)*rho(i_spin)
       enddo

       en_density_xc = en_density_xc / full_rho

    else

       en_density_xc = en_density_x(1)
       full_rho = rho(1)
       spin_pol = 0.d0
       squared_grad_rho = aux_squared_grad_rho
       do i_coord = 1,3,1
          total_rho_gradient(i_coord) = rho_gradient(i_coord,1)
       enddo
       ! non-polarized case; squared_grad_rho is already known from exchange part

    end if

    !       get correlation terms

    if (dabs(full_rho).gt.1d-20) then
      call corpbe_derivs &
         (  full_rho, squared_grad_rho, spin_pol, &
         en_density_c, c_density_deriv, c_gradient_deriv, &
         c_spin_deriv &
         )
    else
      en_density_c = 0.d0
      c_density_deriv = 0.d0
      c_gradient_deriv = 0.d0
      c_spin_deriv = 0.d0
    end if

    ! add relevant terms together
    ! here we add "ax" of the PBE exchange to "ac" PBE correlation.
    en_density_xc = ax * en_density_xc &
         + ac * en_density_c

    do i_spin = 1, n_spin, 1
     en_density_x(i_spin)= ax * en_density_x(i_spin)
    enddo

    !  that's all folks

  end subroutine xc_partials_zrs1_dft_part
!******
!-------------------------------------------------------------------
subroutine xc_partials_tpss_family &
           ( rho, rho_gradient, kinetic_density, en_density_xc, en_density_x, en_density_c, &
             local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, tpss_option, &
             tau_threshold )

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8, dimension(n_spin) :: xc_tau_deriv
    real*8, dimension(n_spin) :: kinetic_density
    integer, intent(in) :: tpss_option
    real*8, intent(in) :: tau_threshold

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !   o kinetic_density -- Gradient of the orbitals
    !   o tpss_option -- Choice for TPSS functionals. The options are:
    !                   1 : TPSS 
    !                   2 : revTPSS
    !                   3 : TPSSloc
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !       en_density_x     : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                          Really ex = ex( rho, squared_grad_rho, tau )
    !       en_density_c     : Correlation energy density ec(x) = ec_lda(x) + H(x)
    !                          Really ec = ec( rho, squared_grad_rho, tau )
    !   o local_xc_derivs -- Partial derivative d/d(rho) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_gradient_deriv -- Partial derivative d/d(|grad(rho)|^2) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_tau_deriv -- Partial derivative d/d(tau) of of e_xc(rho,|grad(rho)|^2,tau)
    !
    !  AUTHOR
    !    FHI-aims team.
    !    Updated by Andrew Logsdail, Jan 2015
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    !  Local variables

    real*8, dimension(7) :: derivatives_x
    real*8, dimension(7) :: derivatives_c
    real*8 :: total_rho
    real*8 :: total_tau
    real*8 :: threshold

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    en_density_xc     = 0.d0
    en_density_x      = 0.d0
    en_density_c      = 0.d0
    local_xc_derivs   = 0.d0
    xc_gradient_deriv = 0.d0
    xc_tau_deriv      = 0.0d0

    derivatives_x     = 0.d0
    derivatives_c     = 0.d0
    threshold         = 1.d-10
    total_rho         = 0.d0
    total_tau         = 0.d0

    do i_spin=1, n_spin, 1
     total_rho = total_rho + rho(i_spin)
     total_tau = total_tau + kinetic_density(i_spin)
    enddo

    ! exchange and correlation

    if (total_rho.gt.threshold.and.total_tau.gt.tau_threshold) then
!    if (total_rho.gt.threshold) then

      if (spin_treatment.eq.0) then

        call tpssx(en_density_x,derivatives_x(:),rho(1)*0.5d0,rho(1)*0.5d0,rho_gradient(:,1)*0.5d0,rho_gradient(:,1)*0.5d0, &
                   kinetic_density(1)*0.25d0,kinetic_density(1)*0.25d0,tpss_option)

        call tpssc(en_density_c,derivatives_c(:),rho(1)*0.5d0,rho(1)*0.5d0,rho_gradient(:,1)*0.5d0, rho_gradient(:,1)*0.5d0, &
                   kinetic_density(1)*0.25d0,kinetic_density(1)*0.25d0,tpss_option)

      else ! (spin_treatment.eq.1)

        call tpssx(en_density_x,derivatives_x(:),rho(1),rho(2),rho_gradient(:,1),rho_gradient(:,2),&
                   kinetic_density(1)*0.5d0,kinetic_density(2)*0.5d0,tpss_option)

        call tpssc(en_density_c,derivatives_c(:),rho(1),rho(2),rho_gradient(:,1),rho_gradient(:,2),&
                   kinetic_density(1)*0.5d0,kinetic_density(2)*0.5d0,tpss_option)

      endif        

      en_density_xc = (en_density_x + en_density_c)/total_rho
      en_density_x  = en_density_x/total_rho
      en_density_c  = en_density_c/total_rho

! We'll use the potential arrays as they are. AJL
! Indices for derivative_x and derivative_c:
! 1 - dF/d(rho_alpha)
! 2 - dF/d(rho_beta)
! 3 - dF/d(gamma_alpha.alpha)
! 4 - dF/d(gamma_beta.beta)
! 5 - dF/d(gamma_alpha.beta)
! 6 - dF/d(tau_alpha)
! 7 - dF/d(tau_beta)

      do i_spin = 1, n_spin, 1

         ! Derivatives with respect to rho 
         ! local_xc_derivs(i_spin) = dble(n_spin) * 0.5d0 * (derivatives_x(i_spin) + derivatives_c(i_spin))
         ! if (spin_treatment.eq.0) then
         !    ! Remove dble(n_spin) from these sums as it is always 1 in this case.
         !    local_xc_derivs(i_spin) = local_xc_derivs(i_spin) + 0.5d0 * (derivatives_x(i_spin+1) + derivatives_c(i_spin+1))
         ! endif

         ! Simplifies to:
         local_xc_derivs(i_spin) = 1.0d0 * (derivatives_x(i_spin) + derivatives_c(i_spin))

         ! Derivatives with respect to gamma
         do i_coord = 1,3,1
            ! There is a gamma_alpha.beta contribution
            if (spin_treatment.eq.0) then
               ! Remove dble(n_spin) from these sums as it is always 1 in this case.
               xc_gradient_deriv(i_coord, i_spin) = & 
                    0.5d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(3) + derivatives_c(3)) + &
                    0.25d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(5) + derivatives_c(5))
            else
               !xc_gradient_deriv(i_coord, i_spin) = &
               !      dble(n_spin) * 0.25d0 * rho_gradient(i_coord,i_spin) * dble(n_spin) * &
               !      (derivatives_x(i_spin+2) + derivatives_c(i_spin+2)) + &
               !      0.25d0 * rho_gradient(i_coord,3-i_spin) * dble(n_spin) * &
               !      (derivatives_x(5) + derivatives_c(5))

               ! Remove dble(n_spin) from these sums as it is always 2 in this case.
               xc_gradient_deriv(i_coord, i_spin) = &
                    1.0d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(i_spin+2) + derivatives_c(i_spin+2)) + &
                    0.5d0 * rho_gradient(i_coord,3-i_spin) * &
                    (derivatives_x(5) + derivatives_c(5))
            endif
         enddo

         ! Derivatives with respect to tau
         ! xc_tau_deriv(i_spin) = dble(n_spin) * 0.25d0 * (derivatives_x(i_spin+5) + derivatives_c(i_spin+5))
         ! if (spin_treatment.eq.0) then
         !    ! Remove dble(n_spin) from these sums as it is always 1 in this case.
         !    xc_tau_deriv(i_spin) = xc_tau_deriv(i_spin) + 0.25d0 * (derivatives_x(i_spin+6) + derivatives_c(i_spin+6))
         ! endif

         ! Simplifies to:
         xc_tau_deriv(i_spin) = 0.5d0 * (derivatives_x(i_spin+5) + derivatives_c(i_spin+5))

      enddo

    endif ! total_rho.gt.threshold
 
end subroutine xc_partials_tpss_family
!**********************************************************************************************


subroutine xc_partials_dfauto( &
        func_code, rho, rho_gradient, kinetic_density, &
        en_density_xc, en_density_x, en_density_c, &
        local_xc_derivs, xc_gradient_deriv, xc_tau_deriv)
    use xc_library,  only : xc__pw_lda, xc__pbe, xc__pbe0, xc__tpss, xc__scan, xc__scan0
    use constants

    implicit none

    integer, intent(in) :: func_code
    real*8, intent(in) :: rho(n_spin)
    real*8, intent(in) :: rho_gradient(3, n_spin)
    real*8, intent(in) :: kinetic_density(n_spin)
    real*8, intent(out) :: en_density_xc, en_density_x, en_density_c
    real*8, intent(out) :: local_xc_derivs(n_spin)
    real*8, intent(out) :: xc_gradient_deriv(3, n_spin)
    real*8, intent(out) :: xc_tau_deriv(n_spin)

    ! AJL: Debug mode doesn't like that e.g. tauc and tauo can be uninitialised,
    ! so I'm just going to initialise everything to zero here

    real*8, dimension(1) :: &
        rhoc, rhoo, sigmacc, sigmaco, sigmaoo, tauc, tauo, upsilonc, upsilono, &
        zk, vrhoc, vrhoo, vsigmacc, vsigmaco, vsigmaoo, vtauc, vtauo, &
        vupsilonc, vupsilono
    real*8 :: sigmaaa, sigmabb, sigmaab, coeff
    character(len=20) :: c_func, x_func

    rhoc(:) = 0.d0
    rhoo(:) = 0.d0
    sigmacc(:) = 0.d0
    sigmaoo(:) = 0.d0
    tauc(:) = 0.d0
    tauo(:) = 0.d0
    upsilonc(:) = 0.d0
    upsilono(:) = 0.d0
    zk(:) = 0.d0
    sigmaaa = 0.d0
    sigmabb = 0.d0
    sigmaab = 0.d0
    coeff = 0.d0

    if (spin_treatment == 1) then
        rhoc = sum(rho)
        rhoo = rho(1)-rho(2)
        if (func_code > 0) then
            sigmaaa = sum(rho_gradient(:, 1)**2)
            sigmabb = sum(rho_gradient(:, 2)**2)
            sigmaab = sum(rho_gradient(:, 1)*rho_gradient(:, 2))
            sigmacc = sigmaaa+2*sigmaab+sigmabb
            sigmaco = sigmaaa-sigmabb
            sigmaoo = sigmaaa-2*sigmaab+sigmabb
        end if
        if (func_code > 1) then
            tauc = sum(kinetic_density)
            tauo = kinetic_density(1)-kinetic_density(2)
        end if
    else
        rhoc = rho
        if (func_code > 0) then
            sigmacc = sum(rho_gradient(:, 1)**2)
        end if
        if (func_code > 1) then
            tauc = kinetic_density
        end if
    end if

    vrhoc(:) = 0.d0
    vrhoo(:) = 0.d0
    vsigmacc(:) = 0.d0
    vsigmaco(:) = 0.d0
    vsigmaoo(:) = 0.d0
    vtauc(:) = 0.d0
    vtauo(:) = 0.d0
    vupsilonc(:) = 0.d0
    vupsilono(:) = 0.d0

    select case (func_code)
        case (xc__pw_lda)
            x_func = 'DIRACX'
            c_func = 'PW92C'
        case (xc__pbe, xc__pbe0)
            x_func = 'PBEX'
            c_func = 'PBEC'
        case (xc__tpss)
            x_func = 'TPSSX'
            c_func = 'TPSSC'
        case (xc__scan, xc__scan0)
            x_func = 'SCANX'
            c_func = 'SCANC'
    end select
    call xc_func(trim(x_func), spin_treatment == 1, &
        rhoc, rhoo, sigmacc, sigmaco, sigmaoo, tauc, tauo, &
        upsilonc, upsilono, &
        zk, vrhoc, vrhoo, &
        vsigmacc, vsigmaco, vsigmaoo, vtauc, vtauo, &
        vupsilonc, vupsilono)
    en_density_x = zk(1)/rhoc(1)
    if (use_hartree_fock) then
        coeff = 1.d0-hybrid_coeff
        en_density_x = coeff*en_density_x
        vrhoc(1) = coeff*vrhoc(1)
        vrhoo(1) = coeff*vrhoo(1)
        vsigmacc(1) = coeff*vsigmacc(1)
        vsigmaco(1) = coeff*vsigmaco(1)
        vsigmaoo(1) = coeff*vsigmaoo(1)
        vtauc(1) = coeff*vtauc(1)
        vtauo(1) = coeff*vtauo(1)
        vupsilonc(1) = coeff*vupsilonc(1)
        vupsilono(1) = coeff*vupsilono(1)
    end if
    call xc_func(trim(c_func), spin_treatment == 1, &
        rhoc, rhoo, sigmacc, sigmaco, sigmaoo, tauc, tauo, &
        upsilonc, upsilono, &
        zk, vrhoc, vrhoo, &
        vsigmacc, vsigmaco, vsigmaoo, vtauc, vtauo, &
        vupsilonc, vupsilono)
    en_density_c = zk(1)/rhoc(1)

    en_density_xc = en_density_x+en_density_c
    if (spin_treatment == 1) then
        local_xc_derivs(1) = vrhoc(1)+vrhoo(1)
        local_xc_derivs(2) = vrhoc(1)-vrhoo(1)
        if (func_code > 0) then
            xc_gradient_deriv(:, 1) = &
                (rho_gradient(:, 1)+rho_gradient(:, 2))*vsigmacc(1) &
                +rho_gradient(:, 1)*vsigmaco(1) &
                +(rho_gradient(:, 1)-rho_gradient(:, 2))*vsigmaoo(1)
            xc_gradient_deriv(:, 2) = &
                (rho_gradient(:, 1)+rho_gradient(:, 2))*vsigmacc(1) &
                -rho_gradient(:, 2)*vsigmaco(1) &
                -(rho_gradient(:, 1)-rho_gradient(:, 2))*vsigmaoo(1)
        end if
        if (func_code > 1) then
            xc_tau_deriv(1) = vtauc(1)+vtauo(1)
            xc_tau_deriv(2) = vtauc(1)-vtauo(1)
        end if
    else
        local_xc_derivs = vrhoc
        if (func_code > 0) then
            xc_gradient_deriv(:, 1) = rho_gradient(:, 1)*vsigmacc(1)
        end if
        if (func_code > 1) then
            xc_tau_deriv = vtauc
        end if
    end if
end subroutine xc_partials_dfauto


subroutine xc_partials_meta_scan &
           ( rho, rho_gradient, kinetic_density, en_density_xc, en_density_x, en_density_c, &
             local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, option, &
             tau_threshold )

    use constants
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin) :: rho
    real*8, dimension(3,n_spin) :: rho_gradient
    real*8 :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8, dimension(n_spin) :: xc_tau_deriv
    real*8, dimension(n_spin) :: kinetic_density
    integer :: option
    real*8, intent(in) :: tau_threshold

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !   o kinetic_density -- Gradient of the orbitals
    !   
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !       en_density_x  : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                       Really ex = ex( rho, squared_grad_rho, tau )
    !       en_density_c  : Correlation energy density ec(x) = ec_lda(x)+H(x)
    !                       Really ec = ec( rho, squared_grad_rho, tau )
    !   o local_xc_derivs -Partial d/d(rho) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_gradient_deriv -- Partial derivative d/d(|grad(rho)|^2) 
    !                          of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_tau_deriv -Partial d/d(tau) of e_xc(rho,|grad(rho)|^2,tau)
    !
    !  AUTHOR
    !    FHI-aims team.
    !    Updated by Adrienn Ruzsinszky, July 2015
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    !  Local variables

    real*8, dimension(7) :: derivatives_x
    real*8, dimension(7) :: derivatives_c
    real*8 :: total_rho
    real*8 :: total_tau
    real*8 :: threshold

    !       counters

    integer :: i_coord
    integer :: i_spin

    !  begin work
    en_density_xc     = 0.d0
    en_density_x      = 0.d0
    en_density_c      = 0.d0
    local_xc_derivs   = 0.d0
    xc_gradient_deriv = 0.d0
    xc_tau_deriv      = 0.0d0

    derivatives_x     = 0.d0
    derivatives_c     = 0.d0
    threshold         = 1.d-10
    total_rho         = 0.d0
    total_tau         = 0.d0

    do i_spin=1, n_spin, 1
     total_rho = total_rho + rho(i_spin)
     total_tau = total_tau + kinetic_density(i_spin)
    enddo

    ! exchange and correlation

!    if (total_rho.gt.threshold) then
    if (total_rho.gt.threshold.and.total_tau.gt.tau_threshold) then

      if (spin_treatment.eq.0) then

        call meta_scan_x(en_density_x,derivatives_x(:),rho(1)*0.5d0,rho(1)*0.5d0, &
   &               rho_gradient(:,1)*0.5d0,rho_gradient(:,1)*0.5d0,&
   &               kinetic_density(1)*0.25d0,kinetic_density(1)*0.25d0,option)

        call meta_scan_c(en_density_c,derivatives_c(:),rho(1)*0.5d0,rho(1)*0.5d0, &
   &               rho_gradient(:,1)*0.5d0, rho_gradient(:,1)*0.5d0, & 
   &               kinetic_density(1)*0.25d0,kinetic_density(1)*0.25d0,option)

      else ! (spin_treatment.eq.1)

        call meta_scan_x(en_density_x,derivatives_x(:),rho(1),rho(2), &
   &               rho_gradient(:,1),rho_gradient(:,2), &
                   kinetic_density(1)*0.5d0,kinetic_density(2)*0.5d0,option)

        call meta_scan_c(en_density_c,derivatives_c(:),rho(1),rho(2), &
   &               rho_gradient(:,1),rho_gradient(:,2), &
   &               kinetic_density(1)*0.5d0,kinetic_density(2)*0.5d0,option)

      endif        

      en_density_xc = (en_density_x + en_density_c)/total_rho
      en_density_x  = en_density_x/total_rho
      en_density_c  = en_density_c/total_rho

! We'll use the potential arrays as they are. AJL
! Indices for derivative_x and derivative_c:
!  1 - dF/d(rho_alpha)
!  2 - dF/d(rho_beta)
!  3 - dF/d(gamma_rho)
! 3 - dF/d(gamma_alpha.alpha) ! do not use
! 4 - dF/d(gamma_beta.beta)   ! do not use
! 5 - dF/d(gamma_alpha.beta)  ! do not use
! 6 - dF/d(tau_rho)
! 7 - dF/d(tau_beta)          ! do not use

      do i_spin = 1, n_spin, 1

! Derivatives with respect to rho 
! local_xc_derivs(i_spin) = dble(n_spin) * 0.5d0 * (derivatives_x(i_spin) + &
!     derivatives_c(i_spin))
! if (spin_treatment.eq.0) then
!    ! Remove dble(n_spin) from these sums as it is always 1 in this case.
!    local_xc_derivs(i_spin) = local_xc_derivs(i_spin) + 0.5d0 * (derivatives_x(i_spin+1) + &
!    derivatives_c(i_spin+1))
! endif

! Simplifies to:
         local_xc_derivs(i_spin) = 1.0d0 * (derivatives_x(i_spin) + &
   &         derivatives_c(i_spin))

! Derivatives with respect to gamma
         do i_coord = 1,3,1
            ! There is a gamma_alpha.beta contribution
            if (spin_treatment.eq.0) then
            ! Remove dble(n_spin) from these sums as it is always 1 in this case.
               xc_gradient_deriv(i_coord, i_spin) =  &
   &                0.5d0 * rho_gradient(i_coord,i_spin) * &
   &                (derivatives_x(3) + derivatives_c(3)) + &
   &                0.25d0 * rho_gradient(i_coord,i_spin) * &
   &                (derivatives_x(5) + derivatives_c(5))
            else
!xc_gradient_deriv(i_coord, i_spin) = &
!      dble(n_spin) * 0.25d0 * rho_gradient(i_coord,i_spin) * dble(n_spin) * &
!      (derivatives_x(i_spin+2) + derivatives_c(i_spin+2)) + &
!      0.25d0 * rho_gradient(i_coord,3-i_spin) * dble(n_spin) * &
!      (derivatives_x(5) + derivatives_c(5))

               ! Remove dble(n_spin) from these sums as it is always 2 in this case.
               xc_gradient_deriv(i_coord, i_spin) = &
                    1.0d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x(i_spin+2) + derivatives_c(i_spin+2)) + &
                    0.5d0 * rho_gradient(i_coord,3-i_spin) * &
                    (derivatives_x(5) + derivatives_c(5))
            endif
         enddo

! Derivatives with respect to tau
! xc_tau_deriv(i_spin) = dble(n_spin) * 0.25d0 * (derivatives_x(i_spin+5) + 
!    derivatives_c(i_spin+5))
! if (spin_treatment.eq.0) then
!    ! Remove dble(n_spin) from these sums as it is always 1 in this case.
!    xc_tau_deriv(i_spin) = xc_tau_deriv(i_spin) + 0.25d0 * (derivatives_x(i_spin+6) + 
!    derivatives_c(i_spin+6))
! endif

! Simplifies to:
         xc_tau_deriv(i_spin) = 0.5d0 * (derivatives_x(i_spin+5) + derivatives_c(i_spin+5))

      enddo

    endif ! total_rho.gt.threshold
 
end subroutine xc_partials_meta_scan
! Duplicate functionality disabled. AJL/Dec2016
!*******************************************************************************************
!****f* xc/xc_meta_scan
! PURPOSE
!   Wrapper around an external subroutine.
!******
!subroutine xc_meta_scan ( &
!        nfun, rho, rho_gradient, kinetic_density, en_density_xc, &
!        en_density_x, en_density_c)
!    use constants
!
!    implicit none
!
!    integer, intent(in) :: nfun
!    real*8, intent(in) :: rho(n_spin)
!    real*8, intent(in):: rho_gradient(3,n_spin)
!    real*8, intent(in) :: kinetic_density(n_spin)
!    real*8, intent(out) :: en_density_xc, en_density_x, en_density_c
!
!    real*8 :: derivatives_x(7)
!    real*8 :: derivatives_c(7)
!    real*8 :: total_rho
!    real*8 :: threshold
!
!    integer :: i_coord
!    integer :: i_spin
!
!    en_density_xc = 0.d0
!    en_density_x = 0.d0
!    en_density_c = 0.d0
!    derivatives_x(:) = 0.d0
!    derivatives_c(:) = 0.d0
!    threshold = 1.d-10
!    total_rho = 0.d0
!    do i_spin = 1, n_spin
!        total_rho = total_rho+rho(i_spin)
!    enddo
!
!    if (n_spin == 1) then
!        if (total_rho > threshold) then
!            call meta_scan_x ( &
!                en_density_x, derivatives_x(:), rho(1)/2.d0, rho(1)/2.d0, &
!                rho_gradient(:, 1)/2.d0, rho_gradient(:, 1)/2.d0, &
!                kinetic_density(1)/4.d0, kinetic_density(1)/4.d0, nfun)
!            call meta_scan_c ( &
!                en_density_c, derivatives_c(:), rho(1)/2.d0, rho(1)/2.d0, &
!                rho_gradient(:, 1)/2.d0, rho_gradient(:, 1)/2.d0, & 
!                kinetic_density(1)/4.d0, kinetic_density(1)/4.d0, nfun) 
!         else
!             en_density_x = 0.d0
!             derivatives_x = 0.d0
!             en_density_c = 0.d0
!             derivatives_c = 0.d0
!         endif
!     else if (n_spin == 2) then
!          if (total_rho > threshold) then
!             call meta_scan_x ( &
!                 en_density_x, derivatives_x(:), rho(1), rho(2), &
!                 rho_gradient(:, 1), rho_gradient(:, 2), & 
!                 kinetic_density(1)/2.d0, kinetic_density(2)/2.d0, nfun)
!             call meta_scan_c ( &
!                 en_density_c, derivatives_c(:), rho(1), rho(2), &
!                 rho_gradient(:, 1), rho_gradient(:, 2), & 
!                 kinetic_density(1)/2.d0, kinetic_density(2)/2.d0, nfun) 
!          else
!              en_density_x = 0.d0
!              derivatives_x = 0.d0
!              en_density_c = 0.d0
!              derivatives_c = 0.d0
!          endif        
!      end if
!      en_density_xc = (en_density_x + en_density_c)/total_rho
!      en_density_x = en_density_x/total_rho
!      en_density_c = en_density_c/total_rho
!    end subroutine xc_meta_scan
!
!******
!-------------------------------------------------------------------
subroutine xc_partials_libxc &
           ( rho, rho_gradient, kinetic_density, en_density_xc, en_density_x, en_density_c, &
             local_xc_derivs, xc_gradient_deriv, xc_tau_deriv, libxc_option1, libxc_option2, &
             tau_threshold)

    use constants
    use xc_f03_lib_m
    
    implicit none

    !  ARGUMENTS

    real*8, dimension(n_spin), intent(in) :: rho
    real*8, dimension(3,n_spin), intent(in) :: rho_gradient
    real*8, intent(out) :: en_density_xc, en_density_x, en_density_c
    real*8, dimension(n_spin), intent(out) :: local_xc_derivs
    real*8, dimension(3,n_spin), intent(out) :: xc_gradient_deriv
    real*8, dimension(n_spin), intent(out) :: xc_tau_deriv
    real*8, dimension(n_spin), intent(in) :: kinetic_density
    integer, intent(in) :: libxc_option1, libxc_option2
    real*8, intent(in) :: tau_threshold

    !  INPUTS
    !   o rho -- Density (no factor 4 pi)
    !   o rho_gradient -- Gradient of the density in cartesian coordinates
    !   o kinetic_density -- Gradient of the orbitals
    !
    !
    !  OUTPUT
    !   o en_density_xc -- Exchange-correlation energy density e_xc
    !       en_density_x  : Exchange energy density ex(x) = ex_lda(x)*Fx(x)
    !                       Really ex = ex( rho, squared_grad_rho, tau )
    !       en_density_c  : Correlation energy density ec(x) = ec_lda(x)+H(x)
    !                       Really ec = ec( rho, squared_grad_rho, tau )
    !   o local_xc_derivs -Partial d/d(rho) of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_gradient_deriv -- Partial derivative d/d(|grad(rho)|^2)
    !                          of e_xc(rho,|grad(rho)|^2,tau)
    !   o xc_tau_deriv -Partial d/d(tau) of e_xc(rho,|grad(rho)|^2,tau)
    !
    !  AUTHOR
    !    FHI-aims team, specifically Andrew Logsdail, August 2016
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  COMMENTS: AJL
    !  
    ! If you're reading this it's probably because you think 
    ! there is something dodgey in the implementation (at least that's the only
    ! reason I ever look at these bits of the code!). Included below are
    ! all my benchmarks so you can see the process I went through of testing,
    ! for which I'm happy of LDA, GGA, MGGA and hybrids therein that don't
    ! require fancy settings.
    !
    ! I've mainly benchmarked against MGGA as they are the most complex cases - 
    ! please remember the in-code implementations have been benchmarked separately,
    ! and here we are comparing different implementations of the same XC functionals
    ! so some variance is to be expected with code optimisation/variable truncation
    ! 
    ! Energies for H2O with geometry of: 
    ! atom 0.000000 0.000000 0.000000 O
    ! atom 0.000000 0.000000 0.956914 H
    ! atom 0.926363 0.000000 -0.239868 H
    !
    ! Basis is 6-31+G**, as per pbe.m11 regression test. Relativistic effects turned off.
    !
    ! Benchmark LibXC vs. in-code implementations, no spin (spin)
    ! xc libxc LDA_X+LDA_C_PZ  : -2064.426597577 eV (with spin: -2064.426597577 eV)
    ! xc pz-lda                : -2064.426597522 eV (with spin: -2064.426597522 eV)
    !
    ! xc libxc  GGA_C_PBE+GGA_X_PBE : -2077.538543951 eV (with spin: -2077.538543951 eV)
    ! xc pbe                        : -2077.538543819 eV (with spin: -2077.538543819 eV)
    !
    ! xc libxc  MGGA_X_TPSS+MGGA_C_TPSS : -2079.872165547 eV (with spin: -2079.872165547 eV) (Needed a preconverged PBE wavefunction)
    ! xc tpss                           : -2079.872209123 eV (with spin: -2079.872209123 eV)
    !
    ! xc libxc  MGGA_X_M06_L+MGGA_C_M06_L : -2079.498099090 eV (with spin: -2079.498099090 eV) 
    ! xc m06-l                            : -2079.497951750 eV (with spin: -2079.497951750 eV)
    !
    ! xc libxc  MGGA_X_M11_L+MGGA_C_M11_L : -2079.072594105 eV (with spin: -2079.072594105 eV)
    ! xc m11-l                            : -2079.072594102 eV (with spin: -2079.072594102 eV)
    !
    ! xc libxc  HYB_GGA_XC_B3LYP : -2079.869206090 eV (with spin: -2079.869206089 eV)
    ! xc b3lyp                   : -2079.869202642 eV (with spin: -2079.869202642 eV)
    !
    ! xc libxc  HYB_MGGA_XC_M06_HF : -2078.718023406 eV (with spin: -2078.718023408 eV)
    ! xc m06-hf                    : -2078.717913480 eV (with spin: -2078.717913485 eV)
    !
    ! Forces for CH3, spin-polarised, as in regression tests for FHI-aims (from SCF energy solutions, hence slight deviation)
    ! xc tpss:
    !   1         -0.857996883407696E-11          0.151924259454896E-01         -0.604514199858955E-11
    !   2          0.284715674822049E-11          0.273938218946396E+00          0.201504733286318E-11
    !   3          0.255092180109905E+00         -0.144565322445954E+00          0.201504733286318E-11
    !   4         -0.255092180104172E+00         -0.144565322445932E+00          0.201504733286318E-11
    ! xc libxc MGGA_X_TPSS+MGGA_C_TPSS: (Needed a preconverged PBE wavefunction)
    !   1          0.748299045115608E-11          0.151898002381401E-01          0.104916727076745E-10
    !   2         -0.247128935787725E-11          0.272056997993830E+00         -0.349722423589149E-11
    !   3          0.253463768152332E+00         -0.143623399115965E+00         -0.349722423589149E-11
    !   4         -0.253463768157344E+00         -0.143623399116005E+00         -0.349722423589149E-11
    ! xc m06-hf:
    !   1         -0.206423678706593E-08          0.134989976048364E-01          0.176710540760745E-11
    !   2          0.501743713153363E-09         -0.621634705510792E-01         -0.589035135869151E-12
    !   3         -0.365101833375305E-01          0.243322363116509E-01         -0.589035135869151E-12
    !   4          0.365101849000236E-01          0.243322366345920E-01         -0.589035135869151E-12
    ! xc libxc  HYB_MGGA_XC_M06_HF:
    !   1          0.168350096432201E-08          0.134990020967936E-01         -0.280798993522305E-11
    !   2         -0.446919287995635E-09         -0.621640937328692E-01          0.935996645074349E-12
    !   3         -0.365107222927703E-01          0.243325459170402E-01          0.935996645074349E-12
    !   4          0.365107210561887E-01          0.243325457190353E-01          0.935996645074350E-12
    !
    !  SOURCE

    !  Local variables

    type(xc_f03_func_t) :: xc_func_x, xc_func_c
    type(xc_f03_func_info_t) :: xc_info_x, xc_info_c

    real*8, dimension(3) :: sigma
    real*8, dimension(n_spin) :: lap_rho
    real*8, dimension(n_spin) :: half_kinetic_density
    integer :: libxc_family_x, libxc_family_c 

    real*8, dimension(n_spin) :: derivatives_x_local 
    real*8, dimension(3) :: derivatives_x_sigma      
    real*8, dimension(n_spin) :: derivatives_x_lap_rho 
    real*8, dimension(n_spin) :: derivatives_x_tau   
    real*8, dimension(n_spin) :: derivatives_c_local 
    real*8, dimension(3) :: derivatives_c_sigma      
    real*8, dimension(n_spin) :: derivatives_c_lap_rho 
    real*8, dimension(n_spin) :: derivatives_c_tau   
    real*8, dimension(1) :: temp_en_x, temp_en_c 
    real*8 :: total_rho
    real*8 :: total_tau
    real*8 :: threshold

    !       counters
    integer :: i_coord
    integer :: i_spin

    !  begin work
    en_density_xc = 0.d0
    en_density_x  = 0.d0
    en_density_c  = 0.d0
    temp_en_x(:)  = 0.d0
    temp_en_c(:)  = 0.d0

    derivatives_x_local(:) = 0.d0
    derivatives_x_sigma(:) = 0.d0
    derivatives_x_lap_rho(:) = 0.d0
    derivatives_x_tau(:)   = 0.d0
    derivatives_c_local(:) = 0.d0
    derivatives_c_sigma(:) = 0.d0
    derivatives_c_lap_rho(:) = 0.d0
    derivatives_c_tau(:)   = 0.d0

    threshold         = 1.d-10
    total_rho         = 0.d0
    total_tau         = 0.d0

    ! sigma(:) = 0.d0
    lap_rho(:) = 0.d0
    half_kinetic_density(:) = 0.5d0*kinetic_density(:)
    libxc_family_x = 0
    libxc_family_c = 0

    ! Annoying this needs to be initalised here - would be nice if we pass in a
    ! flag perhaps? As all we need is to know if we are doing a meta-GGA
    ! calculation

    ! initialise - last variable is spin or no spin
    call xc_f03_func_init(xc_func_x, libxc_option1, n_spin)
    xc_info_x = xc_f03_func_get_info(xc_func_x)
    libxc_family_x = xc_f03_func_info_get_family(xc_info_x)
    ! I think I need to close these and reopen them - I can't have multiple
    ! instances open at once for x and c
    call xc_f03_func_end(xc_func_x)

    if (libxc_option2.gt.0) then
      call xc_f03_func_init(xc_func_c, libxc_option2, n_spin)
      xc_info_c = xc_f03_func_get_info(xc_func_c)
      libxc_family_c = xc_f03_func_info_get_family(xc_info_c)
      call xc_f03_func_end(xc_func_c)
    endif

    do i_spin=1, n_spin, 1
      total_rho = total_rho + rho(i_spin)
      if (libxc_family_x.eq.XC_FAMILY_MGGA.or. &
          libxc_family_x.eq.XC_FAMILY_HYB_MGGA.or. &
          libxc_family_c.eq.XC_FAMILY_MGGA.or. &
          libxc_family_c.eq.XC_FAMILY_HYB_MGGA) then
        total_tau = total_tau + kinetic_density(i_spin)
      else
        ! arbitrary setting to make sure we pass logic
        total_tau = 2*tau_threshold
     endif 
    enddo

    ! write(use_unit,*) n_spin, "RHO     : ", rho

    !if (total_rho.gt.threshold.and.total_tau.gt.tau_threshold) then
    if (total_rho.gt.threshold) then
      ! initialise - last variable is spin or no spin
      call xc_f03_func_init(xc_func_x, libxc_option1, n_spin)
      ! xc_info_x = xc_f03_func_get_info(xc_func_x)
      ! libxc_family_x = xc_f03_func_info_get_family(xc_info_x)

      ! Organise the required sigma variable (gamma in my old terminology above)
      if (libxc_family_x.gt.XC_FAMILY_LDA .or. &
          libxc_family_c.gt.XC_FAMILY_LDA) then
        sigma(1) = sum(rho_gradient(:, 1)**2)
        if (spin_treatment.eq.1) then
          sigma(2) = sum(rho_gradient(:, 1)*rho_gradient(:, 2))
          sigma(3) = sum(rho_gradient(:, 2)**2)
        endif
        ! write(use_unit,*) n_spin, "SIGMA   : ", sigma
        ! write(use_unit,*) n_spin, "KED     : ", kinetic_density
      endif

      ! Calculate values of importance
      select case (libxc_family_x)
      case(XC_FAMILY_LDA)
        call xc_f03_lda_exc_vxc(xc_func_x, 1, rho, temp_en_x, derivatives_x_local)
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call xc_f03_gga_exc_vxc(xc_func_x, 1, rho, sigma, temp_en_x, &
          derivatives_x_local, derivatives_x_sigma)
      case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
        call xc_f03_mgga_exc_vxc(xc_func_x, 1, rho, sigma, lap_rho, &
          half_kinetic_density, temp_en_x, derivatives_x_local, derivatives_x_sigma, &
          derivatives_x_lap_rho, derivatives_x_tau)
      end select

      call xc_f03_func_end(xc_func_x)

      ! Check if we are using a second functional component
      if (libxc_option2.gt.0) then
        ! Setup new (correlation) functional choice
        call xc_f03_func_init(xc_func_c, libxc_option2, n_spin)
        ! xc_info_c = xc_f03_func_get_info(xc_func_c)
        ! libxc_family_c = xc_f03_func_info_get_family(xc_info_c)

        ! Calculate values of importance
        select case (libxc_family_c)
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc_vxc(xc_func_c, 1, rho, temp_en_c, derivatives_c_local)
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc(xc_func_c, 1, rho, sigma, temp_en_c, &
            derivatives_c_local, derivatives_c_sigma)
        case(XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc_vxc(xc_func_c, 1, rho, sigma, lap_rho, &
            half_kinetic_density, temp_en_c, derivatives_c_local, derivatives_c_sigma, &
            derivatives_c_lap_rho, derivatives_c_tau)
        end select

        call xc_f03_func_end(xc_func_c)
      endif

      ! Lets bring together the derivatives now
      ! Indices for derivative_x and derivative_c:
      ! 1 - dF/d(rho_alpha)
      ! 2 - dF/d(rho_beta)
      ! 3 - dF/d(gamma_alpha.alpha) - sigma_alpha.alpha in LibXC terminology
      ! 4 - dF/d(gamma_alpha.beta) - sigma_alpha.beta in LibXC terminology
      ! 5 - dF/d(gamma_beta.beta) - sigma_beta.beta in LibXC terminology
      ! 6 - dF/d(tau_alpha)
      ! 7 - dF/d(tau_beta)
 
      ! Derivative_c might be empty, but that shouldn't effect the outcome here as we always
      ! add it to derivatives_x and never multiply with it, and it is initialised to 0.0d0

      ! write(use_unit,*) n_spin, "X_LOCAL: ", derivatives_x_local
      ! write(use_unit,*) n_spin, "C_LOCAL: ", derivatives_c_local
      ! write(use_unit,*) n_spin, "X_SIGMA: ", derivatives_x_sigma
      ! write(use_unit,*) n_spin, "C_SIGMA: ", derivatives_c_sigma
      ! write(use_unit,*) n_spin, "X_TAU  : ", derivatives_x_tau
      ! write(use_unit,*) n_spin, "C_TAU  : ", derivatives_c_tau

      ! Derivatives with respect to rho
      do i_spin = 1, n_spin, 1
        local_xc_derivs(i_spin) = 1.0d0 * (derivatives_x_local(i_spin) + derivatives_c_local(i_spin))
        
        if (libxc_family_x.gt.XC_FAMILY_LDA .or. &
            libxc_family_c.gt.XC_FAMILY_LDA ) then
          ! Derivatives with respect to gamma (sigma in LibXC terminology)
          do i_coord = 1,3,1
            ! Remember: there is a gamma_alpha.beta contribution
            ! First we do spin-paired, then alternatively spin-unpaired
            if (spin_treatment.eq.0) then
               xc_gradient_deriv(i_coord, i_spin) = &
                    1.0d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x_sigma(1) + derivatives_c_sigma(1)) + &
                    0.5d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x_sigma(2) + derivatives_c_sigma(2))
            else
               xc_gradient_deriv(i_coord, i_spin) = &
                    1.0d0 * rho_gradient(i_coord,i_spin) * &
                    (derivatives_x_sigma((2*i_spin)-1) + derivatives_c_sigma((2*i_spin)-1)) + &
                    0.5d0 * rho_gradient(i_coord,3-i_spin) * &
                    (derivatives_x_sigma(2) + derivatives_c_sigma(2))
            endif
          enddo

          ! Derivatives with respect to tau
          if (libxc_family_x.eq.XC_FAMILY_MGGA .or. libxc_family_x.eq.XC_FAMILY_HYB_MGGA .or. &
              libxc_family_c.eq.XC_FAMILY_MGGA .or. libxc_family_c.eq.XC_FAMILY_HYB_MGGA) then
            xc_tau_deriv(i_spin) = 0.5d0 * (derivatives_x_tau(i_spin) + derivatives_c_tau(i_spin))
          endif
        endif
      enddo
 
      ! write(use_unit,*) n_spin, "LOCAL_XC_DERIVS ", local_xc_derivs
      ! write(use_unit,*) n_spin, "XC_GRADIENT_DERIV ", xc_gradient_deriv
      ! write(use_unit,*) n_spin, "XC_TAU_DERIV ", xc_tau_deriv

      ! Finally, organise the energies so they are correct on exit
      ! The energy is already per unit particle, so no need to divide by rho
      if (xc_f03_func_info_get_kind(xc_info_x).eq.XC_EXCHANGE_CORRELATION) then
        ! There is no distinction between exchange and correlation
        en_density_xc = temp_en_x(1) 
      else if (xc_f03_func_info_get_kind(xc_info_x).eq.XC_CORRELATION) then
        ! Exchange and correlation were calculated the wrong way round
        en_density_x = temp_en_c(1)
        en_density_c = temp_en_x(1)  
        en_density_xc = en_density_x + en_density_c 
      else
        ! Everything is in the right place
        en_density_x = temp_en_x(1)
        en_density_c = temp_en_c(1)
        en_density_xc = en_density_x + en_density_c
      endif
    endif

end subroutine xc_partials_libxc

end module xc
