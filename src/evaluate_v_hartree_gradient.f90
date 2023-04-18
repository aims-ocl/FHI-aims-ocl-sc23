!****s* FHI-aims/evaluate_v_hartree_gradient
!  NAME
!    evaluate_v_hartree_gradient
!  SYNOPSIS

subroutine evaluate_v_hartree_gradient & 
  ( dist_tab, dir_tab, trigonom_tab, ylm_tab, dylm_dtheta_tab, &
    scaled_dylm_dphi_tab, delta_v_hartree_multipole_component, & 
    delta_v_hartree_multipole_deriv,  &
    d_v_hartree_free_d_r, &
    l_h_dim, gradient_pot )

!  PURPOSE
!     Subroutine evaluate_delta_v_hartree_part_gradient provides the gradient
!     of the partitioned hartree potential for the difference density for
!     a given integration point
!
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
  real*8, intent(in) :: d_v_hartree_free_d_r
  
  real*8, intent(out) :: gradient_pot(3)


! INPUTS
! o dist_tab -- distance to atoms
! o dir_tab -- direction to atoms
! o trigonom_tab -- values of trigonometric functions
! o ylm_tab -- Y_lm functions
! o dylm_dtheta_tab -- ??????
! o scaled_dylm_dphi_tab -- ??????
! o delta_v_hartree_multipole_component -- radial part of Hartree potential
! o delta_v_hartree_multipole_deriv -- derivative of radial part of Hartree potential
! o d_v_hartree_free_d_r -- gradient of radial part of free atoms Hartree potential 
! o l_h_dim -- total sum of l and m pairs
!
!  OUTPUT
!  o gradient_pot -- gradient of multipole expansion of the hartree potential
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
  
  real*8 :: d_dv_hartree_dr
  real*8 :: metric_dpsi_dtheta_pot
  real*8 :: metric_dpsi_dphi_pot

  real*8, external :: ddot


  ! begin work

  ! MS: Judging from Paulas commit (6c71b41ad4300efadad4c482251951fc3f71b716)
  !     this should result in a better gradient in the limit of r->0
  !     Thus, this improvement is kept and now also included in a call of
  !     evaluate_v_hartree_and_rho_multipole_gradient


  if(sum(abs(dir_tab))< 1e-10)then
    !SR: to clarify this:
    !lim_{x,y,z ->0} Y_1m = sqrt(3/pi)/2. m e_m
    !lim_{r->0} (sum_at R(r) * Y_lm(theta,phi)) = [lim_{r->0} dR(r)/dr] * sqrt(3/pi)/2.
    !angular derivatives of Y_lm vanish, since only the monopole of R(r) does not vanish
    !in the limit r->0 but grad dY_00 = 0
    
     gradient_pot(1) = & 
         -delta_v_hartree_multipole_deriv(4)/2*sqrt(3.0/pi) !* 4*pi4_inv

     gradient_pot(2) = & 
          delta_v_hartree_multipole_deriv(2)/2*sqrt(3.0/pi) !* 4*pi4_inv

     gradient_pot(3) = &
          delta_v_hartree_multipole_deriv(3)/2*sqrt(3.0/pi) !* 4*pi4_inv

  else
     ! Notice that we multiply the necessary metric tensor elements
     ! 1/r (for dpsi/dtheta) and 1/(r sin theta) (for dpsi/dphi)
     ! into the derivatives already here (saves a division operation later I guess)

     d_dv_hartree_dr = ddot(l_h_dim,delta_v_hartree_multipole_deriv,1, ylm_tab,1)

     metric_dpsi_dtheta_pot = & 
          ddot(l_h_dim,delta_v_hartree_multipole_component,1, &
                       dylm_dtheta_tab,1) / dist_tab

     metric_dpsi_dphi_pot = & 
          ddot(l_h_dim,delta_v_hartree_multipole_component,1, &
                       scaled_dylm_dphi_tab,1) / dist_tab

     gradient_pot(1) = & 
          dir_tab(1) * (d_dv_hartree_dr + d_v_hartree_free_d_r) + & 
          (trigonom_tab(4) * trigonom_tab(2)) * metric_dpsi_dtheta_pot - &
          trigonom_tab(3) * metric_dpsi_dphi_pot

     gradient_pot(2) = & 
          dir_tab(2) * (d_dv_hartree_dr + d_v_hartree_free_d_r) + & 
          (trigonom_tab(3) * trigonom_tab(2)) * metric_dpsi_dtheta_pot + &
          trigonom_tab(4) * metric_dpsi_dphi_pot

     gradient_pot(3) = & 
          dir_tab(3) * (d_dv_hartree_dr + d_v_hartree_free_d_r) &
          - trigonom_tab(1) * metric_dpsi_dtheta_pot

  end if


end subroutine evaluate_v_hartree_gradient
!---------------------------------------------------------------------
!******	
