!****s* FHI-aims/cutoff_pot
!  NAME
!   cutoff_pot
!  SYNOPSIS

      real*8 function cutoff_pot ( radius, type, r_cutoff, width, pot_scale )

! PURPOSE
!  cutoff potential:
!  returns cutoff pot'l in Hartree
! ARGUMENTS
      real*8  :: radius
      integer :: type
      real*8  :: r_cutoff
      real*8  :: width
      real*8  :: pot_scale
! INPUTS 
!   o radius -- input radius
!   o type -- cutoff potential type
!   o r_cutoff -- cutoff onset radius
!   o width -- width of potential rise
!   o pot_scale -- overall scaling factor
! OUTPUTS 
!   o cutoff_pot -- cutoff potential in Hartree
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

      real*8 :: poly_1, poly_2
      real*8 :: factor_1, factor_2

!  begin work

      if ( type .eq. 1 ) then
        !
        ! Cutoff shape: v_cut(r) = scale * 5.0 * (r-rcut)^2 / (rcut + width - r)^2
        !
        ! This should be default cutoff. It combines the advantage of a Junquera type cutoff
        ! (strictly infinite cutoff at r_cut + w_cut) with the advantage of an x^2 onset
        ! of the cutoff (curvature of u(r) = kinetic energy part sets on earlier and is thus
        ! more "smeared out" in the Junqera potential). It also improves on Junquera as the
        ! increase towards r_cut+w_cut is much more rapid, i.e. the basis function is pushed
        ! down towards zero more reliably.

        if (radius.le.r_cutoff) then

          cutoff_pot = 0.d0

        else if (radius.lt.(r_cutoff+width)) then
          ! 5.d0 Ha prefactor is hardwired here to make pot_scale=1.d0 a reasonable default.

          factor_1 = (radius - r_cutoff)**2.d0
          factor_2 = ((r_cutoff + width) - radius)**2.d0

          cutoff_pot = &
          ( 5.d0*pot_scale*factor_1/factor_2 )

        else
          ! this part should never be touched -> quasi-infinity here.

          cutoff_pot = 1.0e20

        end if

      else if ( type .eq. 2 ) then
        ! This is the cutoff potential given in Junquera et. al. PRB 64, 235111
        ! The form is V(r) = V0*( exp(-(r_c - r_i)/(r - r_i)) / (r_c - r) )
        ! In our notation V0 = V0, r_i = r_cutoff, r_c = r_cutoff + width

        ! standard parameters:
        ! pot_scale = 1.0d0
        ! width = 5.0d0

        if (radius.le.r_cutoff) then

          cutoff_pot = 0.d0

        else if (radius.lt.(r_cutoff+width)) then
          ! 125 is hardwired here to make pot_scale=1.d0 a reasonable default.

          factor_1 = radius - r_cutoff
          factor_2 = (r_cutoff + width) - radius

!  Ville's default was 125 instead of 25 here ...
!  for a width of 5 bohrs

          cutoff_pot = &
          ( 100.d0*pot_scale*exp(-width/factor_1))/factor_2

        else

          cutoff_pot = 1.0e20

        end if

      else if ( type .eq. 3 ) then
        !
        ! Cutoff shape: v_cut(r) = scale * 200.d0 * exp[ - width/(r-rcut) ] / (rcut + width - r)^2
        ! This is the present default in species_defaults

        if (radius.le.r_cutoff) then

          cutoff_pot = 0.d0

        else if (radius.lt.(r_cutoff+width)) then
          ! 200.d0 Ha prefactor is hardwired here to make pot_scale=1.d0 a reasonable default.

          factor_1 = exp ( -width / (radius - r_cutoff) )
          factor_2 = ((r_cutoff + width) - radius)**2.d0

          cutoff_pot = &
          ( 250.d0*pot_scale*factor_1/factor_2 )

          ! cap the cutoff potential at 10^5 Hartree
          if (cutoff_pot.gt.1.d5) then
            cutoff_pot = 1.d5
          end if

        else
          ! this part should never be touched -> quasi-infinity here.

          cutoff_pot = 1.0d5

        end if

      end if

      return
    end function cutoff_pot
!******
