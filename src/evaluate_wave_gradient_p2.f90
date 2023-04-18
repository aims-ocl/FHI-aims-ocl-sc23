!****s* FHI-aims/evaluate_wave_gradient_p2
!  NAME
!    evaluate_wave_gradient_p2
!  SYNOPSIS

subroutine evaluate_wave_gradient_p2 &
     ( n_compute, n_compute_atoms, n_compute_fns, &
     one_over_dist_tab, dir_tab, trigonom_tab,  &
     l_ylm_max, ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab,  &
     radial_wave, radial_wave_deriv, &
     gradient,  rad_index, wave_index, l_index, l_count, fn_atom, &
     n_zero_compute, zero_index_point )

!  PURPOSE
!     Subroutine evaluate_wave_gradient provides the gradient of each basis
!     function at a given integration point; needed for instance for 
!     relativistic treatment, forces, ...
!
!  USES

  use dimensions, only : l_wave_max
  implicit none

!  ARGUMENTS

  integer n_compute
  integer n_compute_atoms
  integer n_compute_fns

  real*8 one_over_dist_tab(n_compute_atoms)
  real*8 dir_tab(3, n_compute_atoms)
  real*8 trigonom_tab(4, n_compute_atoms)
  integer :: l_ylm_max
  real*8 ylm_tab ((l_ylm_max+1)**2, n_compute_atoms)
  real*8 dylm_dtheta_tab ((l_ylm_max+1)**2, n_compute_atoms)
  real*8 scaled_dylm_dphi_tab ((l_ylm_max+1)**2, n_compute_atoms)

  real*8 radial_wave(n_compute_fns)
  real*8 radial_wave_deriv(n_compute_fns)

  integer :: rad_index(n_compute_atoms)
  integer :: wave_index(n_compute_fns)
  integer :: l_index(n_compute_fns)
  integer :: l_count(n_compute_fns)
  integer :: fn_atom(n_compute_fns)

  integer :: n_zero_compute
  integer :: zero_index_point(n_compute)

  real*8 gradient( n_compute, 3)

!  INPUTS
!  o n_compute -- number of relevant basis functions
!  o n_compute_atoms -- number of relevant atoms
!  o n_compute_fns -- number of relevant radial functions
!  o one_over_dist_tab -- 1/(distance to atoms)
!  o dir_tab -- direction to relevant atoms
!  o trigonom_tab -- values of trigonometric functions
!  o l_ylm_max -- maximum l component
!  o ylm_tab -- Y_lm functions
!  o tdylm_dtheta_tab -- ????????
!  o scaled_dylm_dphi_tab -- ?????????????
!  o radial_wave -- radial basis functions
!  o radial_wave_deriv -- derivative of radial basis functions
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

  real*8 dpsi_dr(2*l_ylm_max+1)
  real*8 metric_dpsi_dtheta(2*l_ylm_max+1)  
  real*8 metric_dpsi_dphi(2*l_ylm_max+1)

  real*8 e_theta(n_compute_atoms, 3)

  real*8  :: aux_radial_1 (n_compute_fns)
  real*8  :: aux_radial_2 (n_compute_fns)


  !     counters

  integer :: i_compute_point

  integer :: i_compute
  integer :: i_compute_fn
  integer :: i_compute_atom

  integer :: index_start
  integer :: index_end

  !     begin work

  !     tabulate unit vector e_theta in local spherical coordinate system
  !     ( e_r is already given through dir_tab, 
  !       e_phi is simply [-sin(phi),cos(phi),0] in trigonom_tab -
  !       no need to compute e_r, e_phi again.)

  !     cos(phi) * cos(theta)
  e_theta(1:n_compute_atoms, 1) =  &
       trigonom_tab(4,1:n_compute_atoms) *  &
       trigonom_tab(2,1:n_compute_atoms)

  !     sin(phi) * cos(theta)
  e_theta(1:n_compute_atoms, 2) =  &
       trigonom_tab(3,1:n_compute_atoms) *  &
       trigonom_tab(2,1:n_compute_atoms)

  !       - sin(theta)
  e_theta(1:n_compute_atoms, 3) =  &
       - trigonom_tab(1,1:n_compute_atoms)
  !      enddo

  !      now go over all atoms and perform tasks that depend only on the atom
  index_start = 1
  do i_compute_atom = 1, n_compute_atoms, 1

     index_end = rad_index(i_compute_atom)

     aux_radial_2   ( index_start:index_end ) = &
          radial_wave( index_start:index_end ) * &
          one_over_dist_tab( i_compute_atom )**2

     aux_radial_1   ( index_start:index_end ) = &
          radial_wave_deriv( index_start:index_end ) * &
          one_over_dist_tab( i_compute_atom ) - &
          aux_radial_2( index_start:index_end )

     index_start = index_end+1

  enddo

  !     tabulate wave function derivatives for each basis function

  ! first, the nonzero functions
  do i_compute_fn = 1, n_compute_fns, 1

        call mul_vec_3 ( &
             dpsi_dr, l_count(i_compute_fn)+1, &
             ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), &
             aux_radial_1(i_compute_fn) &
             )

        call mul_vec_3 ( &
             metric_dpsi_dtheta, l_count(i_compute_fn)+1, &
             dylm_dtheta_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), &
             aux_radial_2(i_compute_fn) &
             )

        call mul_vec_3 ( &
             metric_dpsi_dphi, l_count(i_compute_fn)+1, &
             scaled_dylm_dphi_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), &
             aux_radial_2(i_compute_fn) &
             )

        ! The following would vectorize but does not because I wrote it
        ! badly. It can be written as a 3x3 transformation of an array of vectors 1:l_count+1
        ! can someone rewrite it that way?
        gradient( wave_index(i_compute_fn):wave_index(i_compute_fn)+l_count(i_compute_fn), 1) = &
             dir_tab( 1,fn_atom(i_compute_fn) )*dpsi_dr(1:l_count(i_compute_fn)+1) + &
             e_theta( fn_atom(i_compute_fn),1 ) *  &
             metric_dpsi_dtheta(1:l_count(i_compute_fn)+1) - &
             trigonom_tab( 3,fn_atom(i_compute_fn) ) *  &
             metric_dpsi_dphi(1:l_count(i_compute_fn)+1)
        
        gradient( wave_index(i_compute_fn):wave_index(i_compute_fn)+l_count(i_compute_fn), 2 ) =  &
             dir_tab( 2,fn_atom(i_compute_fn) )* &
             dpsi_dr(1:l_count(i_compute_fn)+1) + &
             e_theta( fn_atom(i_compute_fn),2 ) *  &
             metric_dpsi_dtheta(1:l_count(i_compute_fn)+1) + &
             trigonom_tab( 4,fn_atom(i_compute_fn) ) *  &
             metric_dpsi_dphi(1:l_count(i_compute_fn)+1)

        gradient( wave_index(i_compute_fn):wave_index(i_compute_fn)+l_count(i_compute_fn), 3) =  &
             dir_tab( 3,fn_atom(i_compute_fn) )* &
             dpsi_dr(1:l_count(i_compute_fn)+1) + &
             e_theta( fn_atom(i_compute_fn),3 ) *  &
             metric_dpsi_dtheta(1:l_count(i_compute_fn)+1)

  enddo

  ! then, the zero functions
  do i_compute_point = 1, n_zero_compute, 1
     i_compute = zero_index_point(i_compute_point)

     gradient( i_compute,1:3) = 0.d0

  enddo

end subroutine evaluate_wave_gradient_p2
!******	
!---------------------------------------------------------------------
!****s* FHI-aims/mul_vec_3
!  NAME
!   mul_vec_3
!  SYNOPSIS

subroutine mul_vec_3 ( wave, n_mul, ylm, factor  )

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

!  INPUTS
!   o n_mul -- dimension of the vector
!   o ylm -- the vector
!   o factor -- multiplication factor for the vector
!
!  OUTPUT
!   o wave -- results of the multiplication
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


  ! begin work

  wave(1:n_mul) = ylm(1:n_mul) * factor

  !    write(use_unit,*) ylm
  !    write(use_unit,*) factor

end subroutine mul_vec_3

!---------------------------------------------------------------------
!******	
