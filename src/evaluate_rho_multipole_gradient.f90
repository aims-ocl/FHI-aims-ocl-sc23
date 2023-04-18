!****s* FHI-aims/evaluate_rho_multipole_gradient
!  NAME
!    evaluate_rho_multipole_gradient
!  SYNOPSIS

subroutine evaluate_rho_multipole_gradient & 
     ( dist_tab, dir_tab, trigonom_tab, ylm_tab, dylm_dtheta_tab, &
     scaled_dylm_dphi_tab, & 
     rho_multipole_component, &
     rho_multipole_deriv, d_rho_free_d_r, &
     l_h_dim, gradient_rho)

!  PURPOSE
!     Subroutine provides the gradient of the multipole expansion 
!     of electron density at a given integration point
!
!  USES

  use dimensions
  use basis
  use geometry
  use grids
  use spline
  use species_data
  use free_atoms
  use runtime_choices
  use constants
  implicit none

!  ARGUMENTS

  real*8, intent(in) :: dist_tab
  real*8, intent(in) :: dir_tab(3)
  real*8, intent(in) :: trigonom_tab(4)
  real*8, intent(in) :: ylm_tab((l_pot_max+1)**2)

  real*8, intent(in) :: dylm_dtheta_tab((l_pot_max+1)**2)
  real*8, intent(in) :: scaled_dylm_dphi_tab((l_pot_max+1)**2)

  integer, intent(in) :: l_h_dim
  real*8, intent(in) :: rho_multipole_component((l_pot_max+1)**2)
  real*8, intent(in) :: rho_multipole_deriv((l_pot_max+1)**2)
  real*8, intent(in) :: d_rho_free_d_r

  real*8, intent(out) :: gradient_rho(3)


! INPUTS
! o dist_tab -- distance to atoms
! o dir_tab -- direction to atoms
! o trigonom_tab -- values of trigonometric functions
! o ylm_tab -- Y_lm functions
! o dylm_dtheta_tab -- ??????
! o scaled_dylm_dphi_tab -- ??????
! o rho_multipole_component -- radial part of charge density
! o rho_multipole_deriv -- gradient of radial part of charge density
! o d_rho_free_d_r  -- gradient of radial part of free atoms change density
! o l_h_dim -- total sum of l and m pairs
!
!  OUTPUT
!  o gradient_rho -- gradient of multipole expansion of the electron density
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


  ! local variables

  real*8 :: metric_dpsi_dtheta_rho
  real*8 :: metric_dpsi_dphi_rho
  real*8 :: d_rho_multipole_dr

  real*8, external :: ddot


  ! begin work

  ! Notice that we multiply the necessary metric tensor elements
  ! 1/r (for dpsi/dtheta) and 1/(r sin theta) (for dpsi/dphi)
  ! into the derivatives already here (saves a division operation later I guess)

  d_rho_multipole_dr = ddot(l_h_dim,rho_multipole_deriv, 1,ylm_tab,1)

  metric_dpsi_dtheta_rho = &
      ddot(l_h_dim,rho_multipole_component,1, dylm_dtheta_tab,1) / dist_tab

  metric_dpsi_dphi_rho   = &
      ddot(l_h_dim,rho_multipole_component,1, scaled_dylm_dphi_tab,1) / dist_tab

  gradient_rho(1) = & 
      dir_tab(1) * (d_rho_multipole_dr + d_rho_free_d_r) + & 
      (trigonom_tab(4) * trigonom_tab(2)) * metric_dpsi_dtheta_rho - &
      trigonom_tab(3) * metric_dpsi_dphi_rho

  gradient_rho(2) = & 
      dir_tab(2) * (d_rho_multipole_dr + d_rho_free_d_r) + & 
      (trigonom_tab(3) * trigonom_tab(2)) * metric_dpsi_dtheta_rho + &
      trigonom_tab(4) * metric_dpsi_dphi_rho

  gradient_rho(3) = dir_tab(3) * (d_rho_multipole_dr + d_rho_free_d_r) - &
      trigonom_tab(1) * metric_dpsi_dtheta_rho

end subroutine evaluate_rho_multipole_gradient
!---------------------------------------------------------------------
!******	
