!****s* FHI-aims/evaluate_wave_gradient_cartesian_p1
!  NAME
!    evaluate_wave_gradient_cartesian_p1
!  SYNOPSIS

subroutine evaluate_wave_gradient_cartesian_p1 &
     (dist_tab, i_r, dir_tab, index_lm, l_ylm_max, n_compute, i_basis, atom_index_inv, i_basis_fns_inv, &
     radial_wave, radial_wave_deriv, ylm_tab, sum_gradient, gradient, n_compute_atoms)

!  PURPOSE
!     Subroutine evaluate_wave_gradient_cartesian provides the gradient of each basis
!     function for a given integration point based on the expansion of ylm-functions
!     into cartesian terms
!     
!     The subroutine assumes that cartesians and gradient_terms are already evaluated:
!     * subroutine evaluate_cartesians()
!     * subroutine evaluate_gradient_terms() or subroutine evaluate_gradient_and_hessian_terms()
!     (module cartesian_ylm)
!
!  USES

  use dimensions
  use basis
  use geometry
  use grids
  use spline
  use cartesian_ylm
  use runtime_choices
  use pbc_lists
  implicit none


!  ARGUMENTS


  integer, intent(in) :: n_compute, n_compute_atoms
  real*8, dimension(n_compute_atoms), intent(in) :: dist_tab
  real*8, dimension(n_compute_atoms), intent(in) :: i_r
  real*8, dimension(3, n_compute_atoms), intent(in) :: dir_tab
  integer :: l_ylm_max

  integer, dimension(n_centers_basis_T), intent(in) :: i_basis
  integer, dimension(n_centers), intent(in) :: atom_index_inv
  integer, dimension(n_basis_fns, n_centers), intent(in) :: i_basis_fns_inv
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave_deriv

  integer, intent(in) :: index_lm(-l_wave_max:l_wave_max, 0:l_wave_max )
  real*8, dimension((l_wave_max+1) ** 2, n_compute_atoms) :: ylm_tab
  real*8, dimension(3, (l_wave_max+1) ** 2, n_compute_atoms) :: sum_gradient

  real*8, dimension(n_compute, 3) :: gradient

!  INPUTS
!  o n_compute -- number of relevant basis functions
!  o n_compute_atoms -- number of relevant atoms
!  o dist_tab -- distance to relevant atoms
!  o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!  o dir_tab -- direction to relevant atoms
!  o l_ylm_max -- maximum l component
!  o i_basis  -- non-zero  basis functions
!  o atom_index_inv -- inverse citing list of relevant atoms
!  o i_basis_fns_inv -- inverse citing list of radial functions
!  o radial_wave -- radial basis functions
!  o radial_wave_deriv -- derivative of radial basis functions
!  o index_lm -- order of l and m indexis
!  o ylm_tab -- Y_lm functions
!  o sum_gradient -- ????????
!
!  OUTPUT
!  o gradient --  gradient of each basis function for a given integration point
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

	



  ! local variables
  real*8 :: radial_wave_scaled
  real*8 :: radial_wave_deriv_scaled
  real*8 :: prefactor_one
  real*8 :: prefactor_two
  integer :: lm_index
  real*8 :: one_over_dist
  integer :: current_basis, current_basis_atom
  integer :: compute_atom
  integer :: current_fn

  ! counters

  integer :: i_atom
  integer :: i_compute
  integer :: i_coords

  ! begin work





  do i_compute = 1, n_compute, 1
     current_basis = i_basis(i_compute)
     current_basis_atom = Cbasis_to_basis(current_basis)

     compute_atom = atom_index_inv(Cbasis_to_center(current_basis))

     current_fn = i_basis_fns_inv(basis_fn(Cbasis_to_basis(current_basis)), &
          Cbasis_to_center(current_basis))

!   current_fn = i_basis_fns_inv(basis_fn(current_basis_atom), Cbasis_to_atom(current_basis))
     lm_index = index_lm(basis_m(current_basis_atom), basis_l(current_basis_atom))



     if (current_fn .gt. 0 .and. compute_atom > 0) then

        if (dist_tab(compute_atom) .gt. 0.d0) then
           
           ! first scale radial terms in the appropriate way (u/r^2, u'/r)
           one_over_dist = 1. / dist_tab(compute_atom)
           radial_wave_deriv_scaled = radial_wave_deriv(current_fn) * one_over_dist
           radial_wave_scaled = radial_wave(current_fn) * one_over_dist * one_over_dist
           
           prefactor_one = radial_wave_deriv_scaled - (basis_l(current_basis_atom) + 1) * radial_wave_scaled
           prefactor_two = radial_wave_scaled
           do i_coords = 1, 3, 1
              gradient(i_compute, i_coords) = dir_tab(i_coords, compute_atom) * &
                   prefactor_one * ylm_tab(lm_index, compute_atom) + &
                   prefactor_two * sum_gradient(i_coords, lm_index, compute_atom)
           end do
           
        else
           gradient(i_compute, :) = 0.d0
        end if
        
     else
        ! wave function is zero anyway; set gradient to zero explicitly to be sure.
        gradient(i_compute, :) = 0.d0
     end if
     

  end do
end subroutine evaluate_wave_gradient_cartesian_p1

!---------------------------------------------------------------------
!******	
