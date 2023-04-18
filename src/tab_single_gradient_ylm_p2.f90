!****s* FHI-aims/tab_single_gradient_ylm_p2
!  NAME
!   tab_single_gradient_ylm_p2
!  SYNOPSIS

subroutine tab_single_gradient_ylm_p2 & 
  ( trigonom_tab, l_max, l_ylm_max, &
    ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab )


!  PURPOSE
!  Subroutine tab_gradient_ylm is a wrapper which computes the additional 
!  Y_lm functions needed for wave function gradients for single atom center.
!
!
!  USES

  implicit none
!  ARGUMENTS

  real*8 trigonom_tab ( 4 )
  integer l_max 
  integer l_ylm_max
  real*8 ylm_tab( (l_ylm_max+1)**2 )
  real*8 dylm_dtheta_tab( (l_ylm_max+1)**2 )
  real*8 scaled_dylm_dphi_tab( (l_ylm_max+1)**2 )


!  INPUTS
!   o trigonom_tab -- trigonometric functions from subroutine tab_trigonom_p0 
!   o basis_l_max -- maximum l index of basis functions
!   o l_ylm_max -- maximum l index of wave functions
!
!  OUTPUT
!   o ylm_tab -- Y_lm
!   o dylm_dtheta_tab -- dy_(lm)/d(theta)
!   o scaled_dylm_dphi_tab -- 1/sin(theta)*dy_(lm)/d(phi), the
!           scaled derivative which is needed to compute cartesian
!           gradients and avoids the apparent coordinate singularity
!           on the z axis ( sin(theta)=0 ).
!
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
!



!  counters

!  begin work

!       increment table of ylm's
!       add all ylm functions from l = lmax+1 to l = lmax+1
  call increment_ylm_deriv &
       ( trigonom_tab(1), trigonom_tab(2),  &
       trigonom_tab(3), trigonom_tab(4),  &
       0, l_max, &
       ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab )

!  that's all folks

  return
end subroutine tab_single_gradient_ylm_p2

!----------------------------------------------------------------------
!******
