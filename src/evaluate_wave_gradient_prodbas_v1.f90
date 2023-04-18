subroutine evaluate_wave_gradient_prodbas_v1 &
     ( dist_tab, i_r, dir_tab, trigonom_tab, &
     l_ylm_max, ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab, &
     index_lm, n_compute, i_basis, radial_wave, gradient &
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
  implicit none

!  ARGUMENTS

  real*8 dist_tab(n_atoms)
  real*8 i_r(n_atoms)
  real*8 dir_tab(3, n_atoms)
  real*8 trigonom_tab(4, n_atoms)

  integer :: l_ylm_max
  real*8 ylm_tab ( (l_ylm_max+1)**2, n_atoms )
  real*8 dylm_dtheta_tab ( (l_ylm_max+1)**2, n_atoms )
  real*8 scaled_dylm_dphi_tab ( (l_ylm_max+1)**2, n_atoms )

  integer index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )
  integer n_compute
  integer i_basis(n_compute)
  real*8 radial_wave(n_compute)

  real*8 gradient(n_basis, 3)

!  INPUTS
!  o dist_tab -- distance to relevant atoms
!  o i_r -- the distance from current integration point to all atoms
!           in units of the logarithmic grid, i(r)
!  o dir_tab -- direction to relevant atoms
!  o trigonom_tab -- values of trigonometric functions
!  o l_ylm_max -- maximum l component
!  o ylm_tab -- Y_lm functions
!  o index_lm --  order of l and m components
!  o tdylm_dtheta_tab -- ????????
!  o scaled_dylm_dphi_tab -- ?????????????
!  o i_basis -- list of relevant basis functions
!  o radial_wave -- radial basis functions
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

  real*8 radial_wave_deriv
  real*8 dpsi_dr
  real*8 metric_dpsi_dtheta, metric_dpsi_dphi
  real*8 e_theta(3, n_atoms)

  !     counters

  integer :: i_atom
  integer :: i_basis_1
  integer :: i_compute

  !     begin work

  !     tabulate unit vector e_theta in local spherical coordinate system
  !     ( e_r is already given through dir_tab,
  !       e_phi is simply [-sin(phi),cos(phi),0] in trigonom_tab -
  !       no need to compute e_r, e_phi again.)

  do i_atom = 1, n_atoms

     !       cos(phi) * cos(theta)
     e_theta(1, i_atom) = &
          trigonom_tab(4,i_atom) * trigonom_tab(2,i_atom)

     !       sin(phi) * cos(theta)
     e_theta(2, i_atom) = &
          trigonom_tab(3,i_atom) * trigonom_tab(2,i_atom)

     !       - sin(theta)
     e_theta(3, i_atom) = &
          - trigonom_tab(1,i_atom)

  enddo

  !     tabulate grad ( wave function ) for each basis function

  do i_compute = 1, n_compute, 1

     i_basis_1 = i_basis(i_compute)

     if (dist_tab(basis_atom(i_basis_1)).gt.0.d0) then

        !         FIXME - this needs to be computed per radial function only, not per basis fn.
        radial_wave_deriv = &
             val_spline &
             ( i_r(basis_atom(i_basis_1)), &
             basis_deriv_spl(1,1,basis_fn(i_basis_1)), &
             n_grid(species(basis_atom(i_basis_1))) &
             )

        !         wave function derivatives in spherical coordinates
        !         Notice that we multiply the necessary metric tensor elements
        !         1/r (for dpsi/dtheta) and 1/(r sin theta) (for dpsi/dphi)
        !         into the derivatives already here (saves a division operation later I guess)

        !         The following are evil expressions. Perhaps we should rewrite them in
        !         a more legible form. Will there be a performance gain?

        dpsi_dr = &
             ( radial_wave_deriv * dist_tab(basis_atom(i_basis_1)) &
             - radial_wave(i_compute) ) / &
             dist_tab(basis_atom(i_basis_1))**2.d0 * &
             ylm_tab ( index_lm(basis_m(i_basis_1),basis_l(i_basis_1)), &
             basis_atom(i_basis_1) )

        metric_dpsi_dtheta = &
             radial_wave(i_compute) / &
             ( dist_tab(basis_atom(i_basis_1))**2.d0 ) * &
             dylm_dtheta_tab( &
             index_lm(basis_m(i_basis_1),basis_l(i_basis_1)), &
             basis_atom(i_basis_1) )

        metric_dpsi_dphi = &
             radial_wave(i_compute) / &
             ( dist_tab(basis_atom(i_basis_1))**2.d0 ) * &
             scaled_dylm_dphi_tab( &
             index_lm(basis_m(i_basis_1),basis_l(i_basis_1)), &
             basis_atom(i_basis_1) )

        !         gradient in spherical coordinates
        !         [ This is effectively a 3x3 matrix-vector-product, written out.
        !           Yes I know it's ugly, and if it becomes a speed bottleneck,
        !           we may streamline the notation. ]

        gradient(i_compute, 1) = &
             dir_tab(1,basis_atom(i_basis_1)) * dpsi_dr + &
             e_theta(1,basis_atom(i_basis_1)) * metric_dpsi_dtheta - &
             trigonom_tab(3,basis_atom(i_basis_1)) * metric_dpsi_dphi

        !test
        !          if (i_basis_1.eq.1) then
        !            write(use_unit,*) "dpsi_dr ", dpsi_dr
        !            write(use_unit,*) "metric_dpsi_dtheta ", metric_dpsi_dtheta
        !            write(use_unit,*) "metric_dpsi_dphi ", metric_dpsi_dphi
        !            write(use_unit,*) "gradient(1) ", gradient(i_basis_1, 1)
        !          end if
        !test end

        gradient(i_compute, 2) = &
             dir_tab(2,basis_atom(i_basis_1)) * dpsi_dr + &
             e_theta(2,basis_atom(i_basis_1)) * metric_dpsi_dtheta + &
             trigonom_tab(4,basis_atom(i_basis_1)) * metric_dpsi_dphi

        !test
        !          if (i_basis_1.eq.1) then
        !            write(use_unit,*) "gradient(2) ", gradient(i_basis_1, 2)
        !          end if
        !test end

        gradient(i_compute, 3) = &
             dir_tab(3,basis_atom(i_basis_1)) * dpsi_dr + &
             e_theta(3,basis_atom(i_basis_1)) * metric_dpsi_dtheta

        !test
        !          if (i_basis_1.eq.1) then
        !            write(use_unit,*) "gradient(3) ", gradient(i_basis_1, 3)
        !          end if
        !test end

     else
        !          the current integration point is exactly on the atom which
        !          we're looking at - set wave function gradient to zero, it's not defined.

        gradient(i_compute, 1) = 0.d0
        gradient(i_compute, 2) = 0.d0
        gradient(i_compute, 3) = 0.d0

     end if

  enddo

end subroutine evaluate_wave_gradient_prodbas_v1
