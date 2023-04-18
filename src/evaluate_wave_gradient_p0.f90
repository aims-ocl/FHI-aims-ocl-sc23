!****s* FHI-aims/evaluate_wave_gradient_p0
!  NAME
!   evaluate_wave_gradient_p0
!  SYNOPSIS

subroutine evaluate_wave_gradient_p0 &
     ( dist_tab, dir_tab, trigonom_tab,  &
     l_ylm_max, ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab,  &
     index_lm, n_compute, i_basis, radial_wave, radial_wave_deriv, &
     gradient, &
     n_compute_atoms, atom_index_inv, &
     n_compute_fns, i_basis_fns_inv, n_basis_list &
     )

!  PURPOSE
!     Subroutine evaluate_wave_gradient provides the gradient of each basis
!     function at a given integration point; needed for instance for 
!     relativistic treatment, ...
!
!  USES

  use dimensions
  use basis
  use geometry
  use grids
  use spline
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer n_compute_atoms
  integer n_compute_fns
  integer n_basis_list

  real*8 dist_tab(n_compute_atoms)
  real*8 dir_tab(3, n_compute_atoms)
  real*8 trigonom_tab(4, n_compute_atoms)
  integer :: l_ylm_max
  real*8 ylm_tab ((l_ylm_max+1)**2, n_compute_atoms)
  real*8 dylm_dtheta_tab ((l_ylm_max+1)**2, n_compute_atoms)
  real*8 scaled_dylm_dphi_tab ( (l_ylm_max+1)**2, n_compute_atoms)

  integer atom_index_inv(n_centers)

  integer index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )
  integer n_compute

  integer i_basis(n_compute)
  integer i_basis_fns_inv(n_basis_fns,n_centers)
  real*8 radial_wave(n_basis_list)
  real*8 radial_wave_deriv(n_basis_list)

  real*8 gradient( n_compute, 3)


!  INPUTS
!  o n_compute -- number of relevant basis functions
!  o n_compute_atoms -- number of relevant atoms
!  o n_compute_fns -- number of relevant radial functions
!  o dist_tab -- distance to relevant atoms
!  o dir_tab -- direction to relevant atoms
!  o trigonom_tab -- values of trigonometric functions
!  o l_ylm_max -- maximum l component
!  o ylm_tab -- Y_lm functions
!  o tdylm_dtheta_tab -- ????????
!  o scaled_dylm_dphi_tab -- ?????????????
!  o i_basis -- list of relevant basis functions
!  o i_basis_fns_inv -- inverse citing list of radial basis functions
!  o radial_wave -- radial basis functions
!  o radial_wave_deriv -- derivative of radial basis functions
! 
!  OUTPUT
!   o gradient -- gradient stores the 3d cartesian gradient of each basis function at the 
!                 present integration point
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





  !  local variables

  real*8 dpsi_dr
  real*8 metric_dpsi_dtheta, metric_dpsi_dphi
  real*8 e_theta(n_compute_atoms, 3)

  !     counters

  integer :: i_atom
  integer :: current_basis
  integer :: i_compute
  integer :: i_compute_fn
  integer :: i_basis_fn_1
  integer :: atom_compute
  integer :: current_lm

  !     begin work

  !     tabulate unit vector e_theta in local spherical coordinate system
  !     ( e_r is already given through dir_tab, 
  !       e_phi is simply [-sin(phi),cos(phi),0] in trigonom_tab -
  !       no need to compute e_r, e_phi again.)

  !      do i_atom = 1, n_atoms      

  !     cos(phi) * cos(theta)
  e_theta(1:n_compute_atoms, 1) =  &
       trigonom_tab(4,1:n_compute_atoms) *  &
       trigonom_tab(2,1:n_compute_atoms)

  !     sin(phi) * cos(theta)
  e_theta(1:n_compute_atoms, 2) =  &
       trigonom_tab(3,1:n_compute_atoms) *  &
       trigonom_tab(2,1:n_compute_atoms)

  !       - sin(theta)
  e_theta(1:n_compute_atoms, 3) =  &
       - trigonom_tab(1,1:n_compute_atoms)

  !      enddo

  !     tabulate grad ( wave function ) for each basis function


  do i_compute = 1, n_compute, 1
     current_basis = i_basis(i_compute)
     i_basis_fn_1 = i_basis_fns_inv(basis_fn(Cbasis_to_basis(current_basis)), &
          Cbasis_to_center(current_basis))

     !         FIXME - this needs to be computed per radial function only, not per basis fn.

     !         wave function derivatives in spherical coordinates
     !         Notice that we multiply the necessary metric tensor elements
     !         1/r (for dpsi/dtheta) and 1/(r sin theta) (for dpsi/dphi)
     !         into the derivatives already here (saves a division operation later I guess)

     !         The following are evil expressions. Perhaps we should rewrite them in
     !         a more legible form. Will there be a performance gain?

     if  ( i_basis_fn_1 .gt. 0 ) then

        atom_compute = atom_index_inv(Cbasis_to_center(current_basis))
        current_lm = index_lm(basis_m(Cbasis_to_basis(current_basis)), &
             basis_l(Cbasis_to_basis(current_basis)))

        dpsi_dr =  &
             ( radial_wave_deriv(i_basis_fn_1) *  &
             dist_tab(atom_compute)  &
             - radial_wave(i_basis_fn_1) ) /  &
             ( dist_tab (atom_compute) )**2.d0 * &
             ylm_tab ( current_lm, atom_compute )

        metric_dpsi_dtheta = &
             radial_wave( i_basis_fn_1 ) /  &
             ( dist_tab(atom_compute) )**2.d0  * &
             dylm_dtheta_tab( current_lm, atom_compute )

        metric_dpsi_dphi =  &
             radial_wave( i_basis_fn_1 ) /  &
             ( dist_tab(atom_compute) )**2.d0 * &
             scaled_dylm_dphi_tab( current_lm, atom_compute )

        !     gradient in spherical coordinates
        !     [ This is effectively a 3x3 matrix-vector-product, written out.
        !     Yes I know it's ugly, and if it becomes a speed bottleneck, 
        !     we may streamline the notation. ]

        gradient( i_compute, 1) = &
             dir_tab( 1,atom_compute )*dpsi_dr + &
             e_theta( atom_compute,1 ) *  &
             metric_dpsi_dtheta - &
             trigonom_tab( 3,atom_compute ) *  &
             metric_dpsi_dphi

        gradient( i_compute, 2 ) =  &
             dir_tab( 2,atom_compute )* &
             dpsi_dr + &
             e_theta( atom_compute,2 ) *  &
             metric_dpsi_dtheta + &
             trigonom_tab( 4,atom_compute ) *  &
             metric_dpsi_dphi

        gradient(i_compute, 3) =  &
             dir_tab( 3,atom_compute )* &
             dpsi_dr + &
             e_theta( atom_compute,3 ) *  &
             metric_dpsi_dtheta

        !    write(use_unit,*)  'gradient( i_compute, 1)', gradient( i_compute, 1), current_basis

     else
        ! zero result if basis function need not be computed at current point.
        !

        gradient( i_compute,1 ) = 0.d0
        gradient( i_compute,2 ) = 0.d0
        gradient( i_compute,3 ) = 0.d0

     end if

  enddo

end subroutine evaluate_wave_gradient_p0
!---------------------------------------------------------------------
!******	
