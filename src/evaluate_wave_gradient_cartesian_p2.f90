!****s* FHI-aims/evaluate_wave_gradient_cartesian_p2
!  NAME
!    evaluate_wave_gradient_cartesian_p2
!  SYNOPSIS

subroutine evaluate_wave_gradient_cartesian_p2 &
  (  n_compute, n_compute_atoms, n_compute_fns, one_over_dist_tab, &
     dir_tab, l_ylm_max, ylm_tab, radial_wave, radial_wave_deriv, &
     sum_gradient, gradient, &
     rad_index, wave_index, l_index, l_count, fn_atom, &
     n_zero_compute, zero_index_point,n_max_sum_gradient &
   )

!  PURPOSE
!     Subroutine evaluate_wave_gradient_cartesian provides the gradient of each basis
!     function for a given integration point based on the expansion of ylm-functions
!     into cartesian terms
!
!
!
!     VB: Modified version to take advantage of the vectorized "p2" integration
!     infrastructure.
!
!
!     VB - Comment: For some reason, my restructuring of this routine yields hardly
!     any time improvement whatsoever. Presumably the routine did not use much time
!     before, either, but it is not clear to me just where the time used in the
!     GGA forces goes. 
!
!  USES

  use dimensions, only : l_wave_max, n_max_compute_fns_dens
  implicit none

!  ARGUMENTS

  integer, intent(in) :: n_compute, n_compute_atoms, n_compute_fns
  real*8, dimension(n_compute_atoms), intent(in) :: one_over_dist_tab
  real*8, dimension(3, n_compute_atoms), intent(in) :: dir_tab
  integer :: l_ylm_max
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave
  real*8, dimension(n_max_compute_fns_dens) :: radial_wave_deriv
  real*8, dimension((l_ylm_max+1) ** 2, n_compute_atoms) :: ylm_tab
  real*8, dimension((l_wave_max+1) ** 2, 3, n_compute_atoms) :: sum_gradient

  integer :: rad_index(n_compute_atoms)
  integer :: wave_index(n_compute_fns)
  integer :: l_index(n_compute_fns)
  integer :: l_count(n_compute_fns)
  integer :: fn_atom(n_compute_fns)

  integer :: n_zero_compute
  integer :: zero_index_point(n_compute)
  integer :: n_max_sum_gradient

  real*8, dimension(n_max_sum_gradient, 3) :: gradient


!  INPUTS
!  o n_compute -- number of relevant basis functions
!  o n_compute_atoms -- number of relevant atoms
!  o n_compute_fns -- number of relevant radial functions
!  o dir_tab -- direction to relevant atoms
!  o l_ylm_max -- maximum l component
!  o radial_wave -- radial basis functions
!  o radial_wave_deriv -- derivative of radial basis functions
!  o ylm_tab -- Y_lm functions
!  o sum_gradient -- ????????
!  o one_over_dist_tab -- 1/(distance to atoms)
!  o  rad_index -- in the list of n_compute_fns radial functions, 
!                    the END of the radial function list associated with atom i_compute_atom
!  o wave_index -- in the list of non-zero wave functions, the START of
!                     those wave functions associated with current radial function i_compute_fn
!  o l_index -- in the list of ylm functions, the START of
!                    those ylm functions associated with current radial function i_compute_fn
!  o l_count -- in the list of ylm functions, the NUMBER of
!                     ylm functions associated with current radial function i_compute_fn, MINUS 1
!  o i_atom_fns -- in the list of non-zero atoms at the current point, the index
!                     of the atom at which radial function i_compute_fn is centered
!  o fn_atom -- ??????
!  o n_zero_compute -- number of  known zero basis functions at current point
!  o zero_index_point -- list of  known zero basis functions at current point
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
  real*8 :: radial_wave_deriv_scaled(n_compute_fns)

  real*8 :: prefactor_one(n_compute_fns,3)
  real*8 :: prefactor_two(n_compute_fns)


  ! counters

  integer :: i_compute_point

  integer :: i_compute
  integer :: i_compute_fn
  integer :: i_compute_atom

  integer :: index_start
  integer :: index_end


  ! begin work

  ! tabulate grad (wave function) for each basis function
  ! assuming that cartesians and gradient_terms are already evaluated !!!
  ! -> subroutine evaluate_cartesians()
  ! -> subroutine evaluate_gradient_terms() or subroutine evaluate_gradient_and_hessian_terms()
  ! (module cartesian_ylm)

  ! first, do all radial-fn. dependent steps per atom.

      index_start = 1
      do i_compute_atom = 1, n_compute_atoms, 1

        index_end = rad_index(i_compute_atom)

        prefactor_two   ( index_start:index_end ) = &
          radial_wave( index_start:index_end ) * &
          one_over_dist_tab( i_compute_atom )**2.d0

        radial_wave_deriv_scaled   ( index_start:index_end ) = &
          radial_wave_deriv( index_start:index_end ) * &
          one_over_dist_tab( i_compute_atom ) - &
          prefactor_two( index_start:index_end ) * dble( l_count( index_start:index_end ) / 2 + 1 )

        prefactor_one( index_start:index_end , 1) = &
          dir_tab(1, i_compute_atom) * radial_wave_deriv_scaled   ( index_start:index_end )

        prefactor_one( index_start:index_end , 2) = &
          dir_tab(2, i_compute_atom) * radial_wave_deriv_scaled   ( index_start:index_end )

        prefactor_one( index_start:index_end , 3) = &
          dir_tab(3, i_compute_atom) * radial_wave_deriv_scaled   ( index_start:index_end )

        index_start = index_end+1

      enddo

  ! next, all wave function dependent steps for each non-zero radial function

      ! VB: The following loop does not vectorize because the coord index of sum_gradient is on the
      !     inside instead of in the middle, where it should be.

      do i_compute_fn = 1, n_compute_fns, 1

        ! explicit subroutine call for vectorization
        call add_mul_vec &
        ( gradient(wave_index(i_compute_fn), 1), l_count(i_compute_fn)+1, & 
          ylm_tab( l_index(i_compute_fn), fn_atom(i_compute_fn)), prefactor_one(i_compute_fn,1), &
          sum_gradient(l_index(i_compute_fn), 1, fn_atom(i_compute_fn)), prefactor_two(i_compute_fn) &
        )

        call add_mul_vec &
        ( gradient(wave_index(i_compute_fn), 2), l_count(i_compute_fn)+1, & 
          ylm_tab( l_index(i_compute_fn), fn_atom(i_compute_fn)), prefactor_one(i_compute_fn,2), &
          sum_gradient(l_index(i_compute_fn), 2, fn_atom(i_compute_fn)), prefactor_two(i_compute_fn) &
        )

        call add_mul_vec &
        ( gradient(wave_index(i_compute_fn), 3), l_count(i_compute_fn)+1, & 
          ylm_tab( l_index(i_compute_fn), fn_atom(i_compute_fn)), prefactor_one(i_compute_fn,3), &
          sum_gradient(l_index(i_compute_fn), 3, fn_atom(i_compute_fn)), prefactor_two(i_compute_fn) &
        )

      enddo

      ! then, the zero functions
      do i_compute_point = 1, n_zero_compute, 1
        i_compute = zero_index_point(i_compute_point)

              gradient( i_compute,: ) = 0.d0

      enddo

end subroutine evaluate_wave_gradient_cartesian_p2
!******	
!---------------------------------------------------------------------
!****s* FHI-aims/add_mul_vec
!  NAME
!   add_mul_vec
!  SYNOPSIS

subroutine add_mul_vec ( wave, n_mul, ylm, factor, gradient, factor2 )

!  PURPOSE
!  write an explicitly vectorizable multiplication to avoid an index mess
!
!  USES

  implicit none
!  ARGUMENTS

  integer :: n_mul
  real*8 :: wave(1:n_mul)
  real*8 :: ylm(1:n_mul)
  real*8 :: factor
  real*8 :: gradient(1:n_mul)
  real*8 :: factor2


!  INPUTS
!   o n_mul -- dimension of the vector
!   o ylm -- the first vector
!   o factor -- multiplication factor for the first vector
!   o gradient -- the second vector
!   o factor2  -- multiplication factor for the second vector
!  OUTPUT
!   o wave -- results of the multiplication
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



  wave(1:n_mul) = ylm(1:n_mul) * factor + gradient(1:n_mul) * factor2



end subroutine add_mul_vec

!---------------------------------------------------------------------
!******	
