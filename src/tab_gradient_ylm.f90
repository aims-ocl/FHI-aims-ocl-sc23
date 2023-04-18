!****s* FHI-aims/tab_gradient_ylm
!  NAME
!   tab_gradient_ylm
!  SYNOPSIS

subroutine tab_gradient_ylm &
     ( trigonom_tab, basis_l_max, l_ylm_max, &
     ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab &
     )

!  PURPOSE
!  Subroutine tab_gradient_ylm is a wrapper which computes the additional
!  Y_lm functions needed for wave function gradients
!
!  USES

  use dimensions
  use geometry
  implicit none

!  ARGUMENTS

  real*8 trigonom_tab ( 4, n_atoms )
  integer basis_l_max (n_species)
  integer l_ylm_max
  real*8 ylm_tab( (l_ylm_max+1)**2, n_atoms )
  real*8 dylm_dtheta_tab( (l_ylm_max+1)**2, n_atoms )
  real*8 scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_atoms )

!  INPUTS
!    o trigonom_tab -- trigonometric functions from subroutine tab_trigonom_p0 
!    o basis_l_max -- maximum l index of basis functions
!    o l_ylm_max -- maximum l index of wave functions
!
!  OUTPUT
!   o ylm_tab -- Y_lm
!   o  dylm_dtheta_tab -- dy_(lm)/d(theta)
!   o scaled_dylm_dphi_tab -- 1/sin(theta)*dy_(lm)/d(phi), the
!           scaled derivative which is needed to compute cartesian
!           gradients and avoids the apparent coordinate singularity
!           on the z axis ( sin(theta)=0 ).
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




  !  counters

  integer i_atom

  !  begin work

  do i_atom = 1, n_atoms, 1

     !       increment table of ylm's
     !       add all ylm functions from l = lmax+1 to l = lmax+1
     call increment_ylm_deriv &
          ( trigonom_tab(1,i_atom), trigonom_tab(2,i_atom), &
          trigonom_tab(3,i_atom), trigonom_tab(4,i_atom), &
          0, &
          basis_l_max(species(i_atom)), &
          ylm_tab(1,i_atom), &
          dylm_dtheta_tab(1,i_atom), &
          scaled_dylm_dphi_tab(1,i_atom) )

  enddo

  !  that's all folks

  return
end subroutine tab_gradient_ylm

!----------------------------------------------------------------------
!******
