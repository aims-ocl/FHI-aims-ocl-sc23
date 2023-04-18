!****s* FHI-aims/evaluate_density_gradient_complex_p1
!  NAME
!    evaluate_density_gradient_complex_p1
!  SYNOPSIS

subroutine evaluate_density_gradient_complex_p1   &
     (n_points, n_compute, KS_orbital,   &
     KS_orbital_gradient, max_occ_number,  &
     rho_gradient &
     )

!  PURPOSE
!  Subroutine evaluates the gradient of the
!  charge density at one particulat point in space from the (given)
!  KS orbitals and their gradients. This supports complex eigenvectors
!
!  USES

  use dimensions
  use runtime_choices  
  implicit none

!  ARGUMENTS

  integer :: n_points
  integer :: max_occ_number
  integer :: n_compute
  complex*16, dimension(max_occ_number, n_points) :: KS_orbital
  complex*16, dimension(n_states*n_k_points_group,n_max_batch_size,3) :: KS_orbital_gradient
  real*8, dimension(3, n_points) :: rho_gradient

! INPUTS
! o n_points -- number of grid points
! o n_compute -- number of relevant basis functions
! o KS_orbital -- Kohn-Sham orbitals
! o KS_orbital_gradient -- gradients of Kohn-Sham orbitals
! o max_occ_number -- maximum number of occupated states
!
!  OUTPUT
! o rho_gradient -- gradient of electron density
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
!     grad(rho) = 2 * Sum_i [f_i * psi_i * grad(psi_i) ]
!     where psi_i are the (occupied) Kohn-Sham orbitals.

  do i_coord = 1, 3, 1
     do i_point = 1, n_points, 1

        rho_gradient(i_coord,i_point) = &
             0.5*dot_product(&
             KS_orbital(1:max_occ_number,i_point), &
             KS_orbital_gradient(1:max_occ_number,i_point,i_coord)) + &
!
             0.5*dot_product( &
             KS_orbital_gradient(1:max_occ_number,i_point,i_coord), &
             KS_orbital(1:max_occ_number,i_point))

     end do
  end do
      
!     second, multiply by 2
  rho_gradient = 2.0d0 * rho_gradient

!  end work 

end subroutine evaluate_density_gradient_complex_p1
!---------------------------------------------------------------------
!******	
