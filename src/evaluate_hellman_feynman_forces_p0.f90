!****s* FHI-aims/evaluate_hellman_feynman_forces_p0
!  NAME
!    evaluate_hellman_feynman_forces_p0
!  SYNOPSIS

subroutine evaluate_hellman_feynman_forces_p0(hellman_feynman_forces, rho, partition_tab, &
     dist_tab_sq_global, dir_tab_global, n_points)

!  PURPOSE
!    Sums up hellman-feynman-forces for a given batch of integration points.
!
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices

!  ARGUMENTS

  implicit none
  real*8, dimension(n_points), intent(in) :: rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: dir_tab_global
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_sq_global
  real*8, dimension(n_points), intent(in) :: partition_tab
  integer :: n_points
  real*8, dimension(3, n_atoms), intent(inout) :: hellman_feynman_forces


!  INPUTS
!  o  rho -- electron density
!  o  dir_tab_global -- direction to atoms
!  o  dist_tab_sq_global -- (distance to atoms)**2
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  o  hellman_feynman_forces -- Hellman-Feyman force component is added here.
!
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


  integer :: i_point, i_atom, i_coords
  real*8 :: point_term
  real*8 :: atomic_term
  
  do i_point = 1, n_points, 1
     point_term = partition_tab(i_point) * rho(i_point) 
     do i_atom = 1, n_atoms, 1
        atomic_term = species_z(species(i_atom)) * point_term &
                      / (sqrt(dist_tab_sq_global(i_atom, i_point)) ** 3) 
        do i_coords = 1, 3, 1
           hellman_feynman_forces(i_coords, i_atom) = hellman_feynman_forces(i_coords, i_atom) + &
                 atomic_term * dir_tab_global(i_coords, i_atom, i_point)  
        end do
     end do

  end do
  
end subroutine evaluate_hellman_feynman_forces_p0
!******
