!****s* FHI-aims/evaluate_v_hartree_and_rho_multipole_gradient
!  NAME
!    evaluate_v_hartree_and_rho_multipole_gradient
!  SYNOPSIS

subroutine evaluate_v_hartree_and_rho_multipole_gradient & 
     ( dist_tab, dir_tab, trigonom_tab, ylm_tab, dylm_dtheta_tab, &
     scaled_dylm_dphi_tab, delta_v_hartree_multipole_component, & 
     delta_v_hartree_multipole_deriv, rho_multipole_component, &
     rho_multipole_deriv, d_v_hartree_free_d_r, d_rho_free_d_r, &
     l_h_dim, gradient_pot, gradient_rho)

!  PURPOSE
!     Subroutine provides the gradient of the partitioned hartree potential and 
!     the  gradient of multipole expansion of electron density in a given integration point
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
  real*8, intent(in) :: delta_v_hartree_multipole_component((l_pot_max+1)**2)
  real*8, intent(in) :: delta_v_hartree_multipole_deriv((l_pot_max+1)**2)

  integer, intent(in) :: l_h_dim
  real*8, intent(in) :: rho_multipole_component((l_pot_max+1)**2)
  real*8, intent(in) :: rho_multipole_deriv((l_pot_max+1)**2)
  real*8, intent(in) :: d_v_hartree_free_d_r
  real*8, intent(in) :: d_rho_free_d_r

  real*8, intent(out) :: gradient_pot(3)
  real*8, intent(out) :: gradient_rho(3)


! INPUTS
! o dist_tab -- distance to atoms
! o dir_tab -- direction to atoms
! o trigonom_tab -- values of trigonometric functions
! o ylm_tab -- Y_lm functions
! o dylm_dtheta_tab -- ??????
! o scaled_dylm_dphi_tab -- ??????
! o delta_v_hartree_multipole_component -- radial part of Hartree potential
! o delta_v_hartree_multipole_deriv -- derivative of radial part of Hartree potential
! o rho_multipole_component -- radial part of charge density
! o rho_multipole_deriv -- gradient of radial part of charge density
! o d_v_hartree_free_d_r -- gradient of radial part of free atoms Hartree potential 
! o d_rho_free_d_r  -- gradient of radial part of free atoms change density
! o l_h_dim -- total sum of l and m pairs
!
!  OUTPUT
!  o gradient_pot -- gradient of multipole expansion of the hartree potential
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


  ! begin work

  ! MS: divided the task into separate subroutine
  !     as there was no performance reason to have them glued together

  call evaluate_rho_multipole_gradient( &
     dist_tab, dir_tab, trigonom_tab, ylm_tab, dylm_dtheta_tab, &
     scaled_dylm_dphi_tab, &
     rho_multipole_component, &
     rho_multipole_deriv, &
     d_rho_free_d_r, &
     l_h_dim, gradient_rho )

  call evaluate_v_hartree_gradient( &
     dist_tab, dir_tab, trigonom_tab, ylm_tab, dylm_dtheta_tab, &
     scaled_dylm_dphi_tab, &
     delta_v_hartree_multipole_component, &
     delta_v_hartree_multipole_deriv,  &
     d_v_hartree_free_d_r, &
     l_h_dim, gradient_pot )

end subroutine evaluate_v_hartree_and_rho_multipole_gradient
!---------------------------------------------------------------------
!******	
