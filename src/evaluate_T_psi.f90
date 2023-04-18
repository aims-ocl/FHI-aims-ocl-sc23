!****s* FHI-aims/evaluate_T_psi
!  NAME
!    evaluate_T_psi
!  SYNOPSIS

subroutine evaluate_T_psi &
     ( n_compute, n_compute_atoms, n_compute_fns, l_ylm_max, ylm_tab, one_over_dist_tab, &
     T_times_psi, kinetic_wave, zora_operator, rad_index, wave_index, l_index, l_count, fn_atom, &
     n_zero_compute, zero_index_point )

  !  PURPOSE
  !    Evaluates kinetic operator times basis function:  T|phi_il>
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
  real*8  :: kinetic_wave(n_compute_fns)
  real*8  :: zora_operator
  integer :: rad_index(n_compute_atoms)
  integer :: wave_index(n_compute_fns)
  integer :: l_index(n_compute_fns)
  integer :: l_count(n_compute_fns)
  integer :: fn_atom(n_compute_fns)
  integer :: n_zero_compute
  integer :: zero_index_point(n_compute)
  real*8  :: T_times_psi(n_compute)


  !  INPUTS
  !   o n_compute_atoms -- number of relevant atoms
  !   o n_compute_fns -- number of relevant different radial parts of the basis functions
  !   o l_ylm_max -- maximum l-component of the basis functions
  !   o ylm_tab -- sperical harmonics functions
  !   o one_over_dist_tab -- 1/r
  !   o n_compute -- number of non-zero basis functions
  !   o kinetic_wave -- kinetic operator time  basis function
  !   o zora_operator -- ZORA factor
  !
  !   indices for basis functions that are nonzero at current point
  !    o rad_index -- radial index
  !    o wave_index -- wave index
  !    o l_index -- l index
  !    o l_count -- 
  !    o fn_atom -- which atom radial fn functions belong 
  !
  !   indices for known zero basis functions at current point
  !    o n_zero_compute
  !    o zero_index_point
  !
  !  OUTPUT
  !   o T_times_psi -- kinetic operator times wave
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




  real*8  ::T(n_compute_fns)

  !     counters

  integer :: i_compute
  integer :: i_compute_point
  integer :: i_compute_fn
  integer :: i_compute_atom

  integer :: index_start
  integer :: index_end

  !     begin work

  !     tabulate T for each basis function
  !     tabulate total wave function value for each basis function

  if( (flag_rel.eq.REL_none).or.(flag_rel.eq.REL_atomic_zora) &
       .or.(flag_rel.eq.REL_own)) then  

     index_start = 1
     do i_compute_atom = 1, n_compute_atoms, 1

        index_end = rad_index(i_compute_atom)

        T ( index_start:index_end ) = &
             kinetic_wave ( index_start:index_end )

        index_start = index_end+1

     enddo

  else if(flag_rel.eq.REL_zora.or. flag_rel==REL_KOLNING_HARMON)then

     index_start = 1
     do i_compute_atom = 1, n_compute_atoms, 1

        index_end = rad_index(i_compute_atom)

        T ( index_start:index_end ) = &
             zora_operator * kinetic_wave ( index_start:index_end )

        index_start = index_end+1

     enddo

  end if

  index_start = 1
  do i_compute_atom = 1, n_compute_atoms, 1

     index_end = rad_index(i_compute_atom)

     T ( index_start:index_end ) = &
          T( index_start:index_end ) * one_over_dist_tab(i_compute_atom)

     index_start = index_end+1

  enddo

  ! Now tabulate full wave function kinetic energy for each radial function

  ! first, the nonzero functions
  do i_compute_fn = 1, n_compute_fns, 1

     !        write(use_unit,*) i_compute_fn, l_count(i_compute_fn)+1, wave_index(i_compute_fn)
     !        write(use_unit,*) l_index(i_compute_fn), fn_atom(i_compute_fn)

     call mul_vec_2 ( &
          T_times_psi(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, &
          ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), &
          T(i_compute_fn) &
          )

  enddo

  ! then, the zero functions
  do i_compute_point = 1, n_zero_compute, 1
     i_compute = zero_index_point(i_compute_point)

     T_times_psi(i_compute) = 0.0d0

  enddo

end subroutine evaluate_T_psi
!******
!---------------------------------------------------------------------

