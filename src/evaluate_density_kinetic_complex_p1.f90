!***** FHI-aims/evaluate_density_kinetic_complex_p1
!  NAME
!    evaluate_density_kinetic_complex_p1
!  SYNOPSIS

subroutine evaluate_density_kinetic_complex_p1   &
     (n_points, n_compute,   &
     KS_orbital_gradient, max_occ_number,  &
     rho_kinetic &
     )

!  PURPOSE
!  Subroutine evaluates the gradient of the
!  charge density at one particulat point in space from the (given)
!  gradient of KS orbitals. This supports complex eigenvectors
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: max_occ_number
  integer :: n_compute
  complex*16, dimension(n_states*n_k_points_group,n_max_batch_size,3) :: KS_orbital_gradient
  real*8, dimension( n_points) :: rho_kinetic

! INPUTS
! o n_points -- number of grid points
! o n_compute -- number of relevant basis functions
! o KS_orbital_gradient -- gradients of Kohn-Sham orbitals
! o max_occ_number -- maximum number of occupated states
!
!  OUTPUT
! o rho_kinetic -- kinetic electron density
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

!  local variables

  real*8, external :: ddot

!     counters

  integer :: i_coord
  integer :: i_point

!  begin work

!     We do a straightforward implementation of
!     grad(rho) = Sum_i [f_i * grad(psi_i) \cdot grad(psi_i) ]
!     where psi_i are the (occupied) Kohn-Sham orbitals.

     do i_point = 1, n_points, 1

        rho_kinetic(i_point) = &
             dot_product( &
             KS_orbital_gradient(1:max_occ_number,i_point,1), &
             KS_orbital_gradient(1:max_occ_number,i_point,1)) + &
             dot_product( &
             KS_orbital_gradient(1:max_occ_number,i_point,2), &
             KS_orbital_gradient(1:max_occ_number,i_point,2)) + &
             dot_product( &
             KS_orbital_gradient(1:max_occ_number,i_point,3), &
             KS_orbital_gradient(1:max_occ_number,i_point,3))

     end do

!  end work

end subroutine evaluate_density_kinetic_complex_p1

