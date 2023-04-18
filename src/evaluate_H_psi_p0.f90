!****s* FHI-aims/evaluate_H_psi_p0
!  NAME
!   evaluate_H_psi_p0
!  SYNOPSIS

subroutine evaluate_H_psi_p0 &
     ( l_ylm_max, ylm_tab, dist_tab, &
     index_lm, H_times_psi, &
     radial_wave, local_potential_parts,  &
     n_compute, i_basis, &
     n_compute_atoms, &
     atom_index_inv, n_compute_fns,  &
     i_basis_fns_inv, &
     kinetic_wave, zora_operator, n_basis_list)


!  PURPOSE
!  Subroutine evaluates Hamiltonian times basis function.
!
!  USES

  use dimensions
  use basis
  use grids
  use geometry
  use spline
  use runtime_choices
  use pbc_lists
  implicit none
 
!  ARGUMENTS

  integer :: n_compute_atoms
  integer :: n_compute_fns
  integer :: n_basis_list

  integer :: l_ylm_max
  real*8 ylm_tab ( (l_ylm_max+1)**2, n_compute_atoms )
  real*8 dist_tab( n_compute_atoms )
  integer index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max )
  integer :: n_compute
  integer :: i_basis(n_compute)
  real*8 radial_wave(n_basis_list)
  real*8 local_potential_parts

  integer :: atom_index_inv(n_centers)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)

  real*8 :: kinetic_wave(n_basis_list)



  real*8 H_times_psi(n_compute)


!  INPUTS
!   o n_compute_atoms -- number of relevant atoms
!   o n_compute_fns -- number of relevant radial functions
!   o n_basis_list -- total number of basis functions
!   o l_ylm_max -- maximum l index 
!   o ylm_tab -- Y_lm functions
!   o dist_tab -- distance to relevant atoms
!   o index_lm -- order of l and m indexes
!   o n_compute -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
!   o radial_wave -- radial part of basis functions
!   o local_potential_parts -- total potentials
!   o atom_index_inv -- inverse of atom citing list
!   o i_basis_fns_inv -- inverse or basis_fns list
!   o kinetic_wave -- kinetic part of the radian basis functions
!   
!  OUTPUT
!   o H_times_psi -- Hamiltonian times basis function
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

  real*8 T_plus_V(n_basis_list)
  real*8 zora_operator

  !     counters

  integer :: i_compute
  integer :: i_compute_fn
  integer :: i_basis_1


  !     begin work

  !     tabulate T_plus_V for each basis function
  !     tabulate total wave function value for each basis function


  if( (flag_rel.eq.REL_none).or.(flag_rel.eq.REL_atomic_zora) &
       .or.(flag_rel.eq.REL_own)) then  

     do i_compute_fn = 1, n_compute_fns, 1

        !         add T = (e-v_basis)*u(r) to V*u(r)

        T_plus_V(i_compute_fn) =  &
             local_potential_parts *  &
             radial_wave( i_compute_fn  ) &
             + kinetic_wave( i_compute_fn )

     enddo

  else if(flag_rel.eq.REL_zora)then

     do i_compute_fn = 1, n_compute_fns, 1

        !         add T = ()*u(r) to V*u(r)

        T_plus_V(i_compute_fn) =  &
             local_potential_parts *  &
             radial_wave( i_compute_fn  ) &
             + zora_operator * kinetic_wave( i_compute_fn )

     enddo
  end if


  !         multiply (T+V)*u(dist) by (1/dist)*Y_lm(theta,phi)
  do i_compute = 1, n_compute, 1
     i_basis_1 = i_basis(i_compute)

     if (  &
          i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), Cbasis_to_center(i_basis_1)) &
          .gt. 0 ) then

        H_times_psi(i_compute) = T_plus_V( &
             i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), &
             Cbasis_to_center(i_basis_1))) * &
             ylm_tab(  &
             index_lm(basis_m(Cbasis_to_basis(i_basis_1)),basis_l(Cbasis_to_basis(i_basis_1))),  &
             atom_index_inv(Cbasis_to_center(i_basis_1)) ) / &
             dist_tab(  &
             atom_index_inv(Cbasis_to_center(i_basis_1)) )

     else
        H_times_psi(i_compute) = 0.d0
     end if

  enddo

end subroutine evaluate_H_psi_p0
!---------------------------------------------------------------------
!******
