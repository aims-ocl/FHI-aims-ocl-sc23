!****s* FHI-aims/evaluate_wave_hessian_cartesian_p0
!  NAME
!    evaluate_wave_hessian_hessian_p0
!  SYNOPSIS

subroutine evaluate_wave_hessian_cartesian_p0 & 
     ( dist_tab, i_r, dir_tab, index_lm, l_ylm_max, n_compute, i_basis, atom_index_inv, i_basis_fns_inv, &
       radial_wave, radial_wave_deriv, radial_wave_2nd_deriv, ylm_tab, sum_gradient, sum_hessian, hessian, n_atoms_compute)

!  PURPOSE
!     Subroutine evaluate_wave_hessian_cartesian provides the hessian of each basis
!     function at a given integration point based on the expansion of ylm-functions
!     into cartesian terms.
!
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

  integer:: n_atoms_compute

  real*8, dimension(n_atoms_compute), intent(in) :: dist_tab
  real*8, dimension(n_atoms_compute), intent(in) :: i_r
  real*8, dimension(3, n_atoms_compute), intent(in) :: dir_tab
  integer :: l_ylm_max
  integer, intent(in) :: n_compute
  integer, dimension(n_centers_basis_T), intent(in) :: i_basis
  integer, dimension(n_centers), intent(in) :: atom_index_inv
  integer, dimension(n_basis_fns, n_centers), intent(in) :: i_basis_fns_inv
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave_deriv
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave_2nd_deriv
  integer, intent(in) :: index_lm(-l_wave_max:l_wave_max, 0:l_wave_max )
  real*8, dimension((l_wave_max+1) ** 2, n_atoms_compute) :: ylm_tab
  real*8, dimension(3, (l_wave_max+1) ** 2, n_atoms_compute) :: sum_gradient
  real*8, dimension(6, (l_wave_max+1) ** 2, n_atoms_compute) :: sum_hessian

  real*8, dimension(n_max_compute_dens, 6) :: hessian


!  INPUTS
!   o n_atoms_compute -- number of relevant atoms
!   o dist_tab -- distance to atoms
!   o i_r --  Not used anymore
!   o dir_tab -- direction to relevant atoms
!   o l_ylm_max --  Not used anymore
!   o n_compute -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
!   o atom_index_inv -- inverse of atom index list
!   o i_basis_fns_inv -- inverse of basis fns list
!   o radial_wave -- radial basis functions
!   o radial_wave_deriv -- derivative of radial basis functions
!   o radial_wave_2nd_deriv -- 2nd derivative of radial basis functions
!   o index_lm -- order of l and m index pairs
!   o ylm_tab -- Y_lm functions
!   o sum_gradient -- ????????
!   o sum_hessian -- ?????????
!
!  OUTPUT
!   o hessian - ??????????
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
  real*8 :: radial_wave_2nd_deriv_scaled
  real*8 :: prefactor_one
  real*8 :: prefactor_two
  real*8 :: prefactor_two_diag
  real*8 :: prefactor_three
  integer :: lm_index
  real*8 :: temp_inv_dist
  real*8 :: one_over_dist
  integer :: current_basis, current_basis_atom
  integer :: compute_atom
  integer :: current_fn

  ! counters

  integer :: i_atom
  integer :: i_compute
  integer :: i_coords
  integer :: i_coords_2
  integer :: i_counter

  ! begin work

  ! tabulate hessian (wave function) for each basis function
  ! assuming that cartesians and gradient_terms are already evaluated !!!
  ! -> subroutine evaluate_cartesians()
  ! -> subroutine evaluate_gradient_and_hessian_terms()
  ! (module cartesian_ylm)

  do i_compute = 1, n_compute, 1
     current_basis = i_basis(i_compute)
     compute_atom = atom_index_inv(Cbasis_to_center(current_basis))
     current_basis_atom = Cbasis_to_basis(current_basis)

     current_fn = i_basis_fns_inv(basis_fn(current_basis_atom), Cbasis_to_center(current_basis))
     lm_index = index_lm(basis_m(current_basis_atom), basis_l(current_basis_atom))


     if (current_fn .gt. 0) then

           ! first scale radial terms in the appropriate way (u/r^3, u'/r^2, u''/r)
           one_over_dist = 1. / dist_tab(compute_atom)
           radial_wave_2nd_deriv_scaled = radial_wave_2nd_deriv(current_fn) * one_over_dist
           
           temp_inv_dist = one_over_dist * one_over_dist
           radial_wave_deriv_scaled = radial_wave_deriv(current_fn) * temp_inv_dist
           
           temp_inv_dist = temp_inv_dist * one_over_dist
           radial_wave_scaled = radial_wave(current_fn) * temp_inv_dist
           
           prefactor_one   = radial_wave_2nd_deriv_scaled - (2 * basis_l(current_basis_atom) + 3) * radial_wave_deriv_scaled + &
                (basis_l(current_basis_atom) + 1) * (basis_l(current_basis_atom) + 3) * radial_wave_scaled

           prefactor_one = prefactor_one * ylm_tab(lm_index, compute_atom)

           prefactor_two   = radial_wave_deriv_scaled - (basis_l(current_basis_atom) + 1) * radial_wave_scaled

           prefactor_two_diag = prefactor_two * ylm_tab(lm_index, compute_atom)

           prefactor_three = radial_wave_scaled
           
           i_counter = 0
           do i_coords = 1, 3, 1
              ! diagonal elements
              i_counter = i_counter + 1
              hessian(i_compute, i_counter) = dir_tab(i_coords, compute_atom) * (dir_tab(i_coords, compute_atom) * &
                   prefactor_one + 2 * prefactor_two * &
                   sum_gradient(i_coords, lm_index, compute_atom)) + & 
                   prefactor_two_diag + &
                   prefactor_three * sum_hessian(i_counter, lm_index, compute_atom)
              
              ! non-diagonal elements
              do i_coords_2 = i_coords + 1, 3, 1
                 i_counter = i_counter + 1
                 hessian(i_compute, i_counter) = dir_tab(i_coords, compute_atom) * (dir_tab(i_coords_2, compute_atom) * &
                      prefactor_one + & 
                      prefactor_two * sum_gradient(i_coords_2, lm_index, compute_atom))  &
                      + prefactor_two * dir_tab(i_coords_2, compute_atom) & 
                        * sum_gradient(i_coords, lm_index, compute_atom) + &
                      prefactor_three * sum_hessian(i_counter, lm_index, compute_atom)
              end do
           end do
     else
           hessian(i_compute, :) = 0.d0
     end if
        

  end do
  
end subroutine evaluate_wave_hessian_cartesian_p0
!---------------------------------------------------------------------
!******
