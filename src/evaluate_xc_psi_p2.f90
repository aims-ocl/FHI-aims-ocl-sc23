!****s* FHI-aims/evaluate_xc_psi_p2
!  NAME
!   evaluate_xc_psi_p2
!  SYNOPSIS

subroutine evaluate_xc_psi_p2 &
     ( n_compute, n_compute_atoms, n_compute_fns, &
     l_ylm_max, ylm_tab, one_over_dist_tab, &
     radial_wave, xc_times_psi, &
     local_xc_derivs,  &
     rad_index, wave_index, l_index, l_count, fn_atom, &
     n_zero_compute, zero_index_point )


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

  integer :: l_ylm_max
  real*8  :: ylm_tab ( (l_ylm_max+1)**2, n_compute_atoms )
  real*8  :: one_over_dist_tab( n_compute_atoms )
  integer :: n_compute
  real*8  :: radial_wave(n_compute_fns)
  real*8  :: local_xc_derivs

  integer :: rad_index(n_compute_atoms)
  integer :: wave_index(n_compute_fns)
  integer :: l_index(n_compute_fns)
  integer :: l_count(n_compute_fns)
  integer :: fn_atom(n_compute_fns)

  integer :: n_zero_compute
  integer :: zero_index_point(n_compute)

  real*8 xc_times_psi(n_compute)

!  INPUTS
!   o n_compute_atoms -- number of relevant atoms
!   o n_compute_fns -- number of relevant radial functions
!   o l_ylm_max -- maximum l index 
!   o ylm_tab -- Y_lm functions
!   o one_over_dist_tab -- 1/(distance to atoms)
!   o n_compute -- number of relevant basis functions
!   o radial_wave -- radial part of basis functions
!   o local_xc_derivs -- exchange-correlation potentials
!   o rad_index -- in the list of n_compute_fns radial functions, 
!                     the END of the radial function list associated with atom i_compute_atom
!   o wave_index -- in the list of non-zero wave functions, the START of
!                     those wave functions associated with current radial function i_compute_fn
!   o l_index -- in the list of ylm functions, the START of
!                     those ylm functions associated with current radial function i_compute_fn
!   o l_count -- in the list of ylm functions, the NUMBER of
!                     ylm functions associated with current radial function i_compute_fn, MINUS 1
!   o fn_atom -- ??????
!   o n_zero_compute -- number of  known zero basis functions at current point
!   o zero_index_point -- list of  known zero basis functions at current point
!
!  OUTPUT
!   o xc_times_psi -- Hamiltonian times basis function
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

  real*8 V_xc_only(n_compute_fns)

  !     counters

  integer :: i_compute
  integer :: i_compute_point
  integer :: i_compute_fn
  integer :: i_compute_atom

  integer :: index_start
  integer :: index_end

  !     begin work

  !     tabulate V_xc_only for each basis function
  !     tabulate total wave function value for each basis function


  index_start = 1
  do i_compute_atom = 1, n_compute_atoms, 1

      index_end = rad_index(i_compute_atom)

      V_xc_only ( index_start:index_end ) = &
             local_xc_derivs * radial_wave( index_start:index_end ) * &
             one_over_dist_tab(i_compute_atom)

      index_start = index_end+1

  enddo

  ! Now tabulate full wave function kinetic energy for each radial function

  ! first, the nonzero functions
  do i_compute_fn = 1, n_compute_fns, 1

        call mul_vec_2 ( &
             xc_times_psi(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, &
             ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), &
             V_xc_only(i_compute_fn) &
             )
  enddo

  ! then, the zero functions
  do i_compute_point = 1, n_zero_compute, 1
     i_compute = zero_index_point(i_compute_point)

     xc_times_psi(i_compute) = 0.0d0

  enddo

end subroutine evaluate_xc_psi_p2

!--------------------------------------------------------------------
!******
