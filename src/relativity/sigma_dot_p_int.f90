
 subroutine evaluate_sigma_p_waves &
      ( iop, n_compute, n_compute_atoms, n_compute_fns, &
        l_ylm_max, ylm_tab, one_over_dist_tab,  &
        radial_wave, radial_deriv, wave, i_basis, &
        rad_index, wave_index, l_index, l_count, fn_atom, &
        n_zero_compute, zero_index_point )

  use basis, only: basis_k, basis_small_k
  use localorb_io, only: use_unit
  use pbc_lists, only: Cbasis_to_basis, Cbasis_to_basis_s

  implicit none
  integer, intent(IN) :: iop ! =1, t1; =2, t2
  integer, intent(IN) :: n_compute  ! number of non-zero basis functions
 ! (Rundong) In the present strategy I consider n_compute_atoms = n_compute_atoms_s.
 ! So to n_compute_fns, rad_index, l_index, l_count, fn_atom.
 ! If not, the indices would be very complicated. New indices should be rebuilt.
  integer, intent(IN) :: n_compute_atoms ! number of relevant atoms
  integer, intent(IN) :: n_compute_fns   ! number of non-zero radial basis fns
  integer, intent(IN) :: i_basis(n_compute)
  integer, intent(IN) :: l_ylm_max !  maximum of l
  real*8, intent(IN) :: ylm_tab ((l_ylm_max+1)**2, n_compute_atoms ) ! spherical harmonic functions
  real*8, intent(IN) :: one_over_dist_tab(n_compute_atoms) ! 1/r
  real*8, intent(IN) :: radial_wave(n_compute_fns)
  real*8, intent(IN) :: radial_deriv(n_compute_fns)
 ! indices for basis functions that are nonzero at current point:
  integer, intent(IN) :: rad_index(n_compute_atoms)
  integer, intent(IN) :: wave_index(n_compute_fns)
  integer, intent(IN) :: l_index(n_compute_fns)
  integer, intent(IN) :: l_count(n_compute_fns)
  integer, intent(IN) :: fn_atom(n_compute_fns)
 ! indices for known zero basis functions at current point:
  integer, intent(IN) :: n_zero_compute
  integer, intent(IN) :: zero_index_point(n_compute)

  real*8, intent(OUT)  :: wave(n_compute) ! total basis functions.

 ! local variables
  integer :: l_aux, k_aux
  real*8  :: j_aux(n_compute_fns)
  real*8  :: aux_radial (n_compute_fns)
  real*8  :: aux_radial_deriv (n_compute_fns)
  real*8  :: aux_radial_r (n_compute_fns)
  real*8  :: wave1(n_compute)
  real*8  :: wave2(n_compute)
 ! counters
  integer :: i_compute_point
  integer :: i_compute
  integer :: i_compute_fn
  integer :: i_compute_atom
  integer :: index_start
  integer :: index_end

 ! do work on radial functions first
  index_start = 1
  do i_compute_atom = 1, n_compute_atoms, 1

    index_end = rad_index(i_compute_atom)

   ! (Rundong) The radial bases and their derivative output from dft_atom some times differs by a '-' sign from 
   ! the correct value. This is in a mess. Anyway, the following code is the correct way to use them.
    if(iop.eq.1)then ! use small comp. radial basis
       aux_radial (index_start:index_end)       = radial_wave ( index_start:index_end ) * one_over_dist_tab( i_compute_atom )
       aux_radial_r (index_start:index_end)     = aux_radial  ( index_start:index_end ) * one_over_dist_tab( i_compute_atom )
       aux_radial_deriv (index_start:index_end) = -radial_deriv( index_start:index_end ) * one_over_dist_tab( i_compute_atom )
    elseif(iop.eq.2)then ! use large comp. radial basis
       aux_radial (index_start:index_end)       = -radial_wave ( index_start:index_end ) * one_over_dist_tab( i_compute_atom )
       aux_radial_r (index_start:index_end)     = aux_radial  ( index_start:index_end ) * one_over_dist_tab( i_compute_atom )
       aux_radial_deriv (index_start:index_end) = radial_deriv( index_start:index_end ) * one_over_dist_tab( i_compute_atom )
    endif

    index_start = index_end+1

  enddo

 ! get k values and thus j values
  do i_compute_fn = 1, n_compute_fns ! radial wave index

     i_compute = wave_index(i_compute_fn)
     l_aux = l_count(i_compute_fn)/2

     if(iop.eq.1)then
        k_aux = basis_k(Cbasis_to_basis(i_basis(i_compute)))
        if(k_aux.eq.l_aux)then ! spin down
           j_aux(i_compute_fn) = l_aux - 0.5d0
        elseif(k_aux.eq.-l_aux-1)then ! spin up
           j_aux(i_compute_fn) = l_aux + 0.5d0
        else
           write(use_unit,"('INCORRECT k value! subroutine: evaluate_sigma_p_waves in integrate_hamiltonian_matrix_p2')")
        endif
     elseif(iop.eq.2)then
        k_aux = basis_small_k(Cbasis_to_basis_s(i_basis(i_compute)))
       ! Note, j_aux(i_compute_fn) is not l-0.5 or l+0.5 now!
       ! Reason: for small comp. (iop=2) the l_aux we saved is already degraded/upgraded, and 
       ! -- in the present code structure, it is used only for the angular part.
       ! Here, the j value is an index for the radial part. We should use the initial j value.
        if(k_aux.eq.l_aux)then ! spin down
          ! l_aux was degraded, now turn it back to l_aux+1
           j_aux(i_compute_fn) = (l_aux+1.d0) - 0.5d0
        elseif(k_aux.eq.-l_aux-1)then ! spin up
          ! l_aux was upgraded, now turn it back to l_aux-1
           j_aux(i_compute_fn) = (l_aux-1.d0) + 0.5d0
        else
           write(use_unit,"('INCORRECT k value! subroutine: evaluate_sigma_p_waves in integrate_hamiltonian_matrix_p2')")
        endif
     endif
       ! write(use_unit,"('i_compute_fn=',i5,5x,'i_compute=',i5,5x,'j_aux=',f8.4,5x,'l_aux=',i3)") i_compute_fn, i_compute, j_aux(i_compute_fn), l_aux

  enddo


 ! tabulate total wave function value for each basis function

  ! first, the nonzero functions
  do i_compute_fn = 1, n_compute_fns, 1

     call mul_vec ( wave1(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, &
          ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), aux_radial_deriv(i_compute_fn) )

     call mul_vec_sigma_p ( wave2(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, &
          ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), j_aux(i_compute_fn), aux_radial_r(i_compute_fn) )

     wave ( wave_index(i_compute_fn) : wave_index(i_compute_fn)+l_count(i_compute_fn)+1 ) = &
     wave1( wave_index(i_compute_fn) : wave_index(i_compute_fn)+l_count(i_compute_fn)+1 ) + &
     wave2( wave_index(i_compute_fn) : wave_index(i_compute_fn)+l_count(i_compute_fn)+1 )

  enddo

  ! then, the zero functions
  do i_compute_point = 1, n_zero_compute, 1
    i_compute = zero_index_point(i_compute_point)
    wave(i_compute) = 0.0d0
  enddo

 end subroutine evaluate_sigma_p_waves

 subroutine mul_vec_sigma_p ( wave, n_mul, ylm, j, factor )

    implicit none
    integer,intent(in) :: n_mul         ! dimensions of the vectors
    real*8,intent(in)  :: ylm(1:n_mul)  ! spherical harmonics
    real*8,intent(in)  :: j             ! j value of the precent spherical harmonics
    real*8,intent(in)  :: factor        ! multiplication factor
    real*8,intent(out) :: wave(1:n_mul) ! basis functions, results from the multiplication

    integer :: l ! l value of the current angular term

    integer :: i

    l = (n_mul-1)/2

    do i=1, n_mul
       wave(i) = ylm(i) * factor * ( l*(l+1) - j*(j+1) -0.25d0 )
    enddo

 end subroutine mul_vec_sigma_p


 subroutine evaluate_T12_psi_rel (n_compute, wave, Tx_times_psi )

!  PURPOSE
!  Subroutine evaluates V,T times large component basis function, and V times small comp. basis function.

  use constants, only: light_speed

  implicit none

  integer,intent(in) :: n_compute ! Number of relevant basis functions
  real*8,intent(in)  :: wave(n_compute) ! small/large comp. basis functions

!  OUTPUT
!   T1_times_psi -- c times sigma_dot_p times small comp. basis function
!   T2_times_psi -- c times sigma_dot_p times large comp. basis function
  real*8,intent(out) :: Tx_times_psi(n_compute)

  integer :: i_compute

  do i_compute=1, n_compute
     Tx_times_psi(i_compute) = wave( i_compute ) * light_speed
  enddo

 end subroutine evaluate_T12_psi_rel


 subroutine get_small_comp_sigmadotp(i_species, n_grid, i_fn, r_grid, large, large_deriv, eigenval, pot, small, small_kinetic )
  use constants, only : light_speed,light_speed_sq
  implicit none
  integer,intent(in) :: i_species, n_grid, i_fn
  real*8,intent(in)  :: r_grid(n_grid), eigenval, pot(n_grid)
  real*8,intent(in)  :: large(n_grid), large_deriv(n_grid)
  real*8,intent(out) :: small(n_grid), small_kinetic(n_grid)

  integer :: i,k,l,m,n
  real*8 :: j, factor

  ! For Neon:
  if(i_fn.eq.1)then
    l=0; j=0.5d0
  elseif(i_fn.eq.2)then
    l=0; j=0.5d0
  elseif(i_fn.eq.3)then
    l=1; j=0.5d0
  elseif(i_fn.eq.4)then
    l=1; j=1.5d0
  endif

  ! For NaCl:
! if(i_species.eq.1)then ! Cl
!   if(i_fn.eq.1)then
!     l=0; j=0.5d0
!   elseif(i_fn.eq.2)then
!     l=0; j=0.5d0
!   elseif(i_fn.eq.3)then
!     l=0; j=0.5d0
!   elseif(i_fn.eq.4)then
!     l=1; j=0.5d0
!   elseif(i_fn.eq.5)then
!     l=1; j=1.5d0
!   elseif(i_fn.eq.6)then
!     l=1; j=0.5d0
!   elseif(i_fn.eq.7)then
!     l=1; j=1.5d0
!   endif
! else ! Na
!   if(i_fn.eq.1)then
!     l=0; j=0.5d0
!   elseif(i_fn.eq.2)then
!     l=0; j=0.5d0
!   elseif(i_fn.eq.3)then
!     l=0; j=0.5d0
!   elseif(i_fn.eq.4)then
!     l=1; j=0.5d0
!   elseif(i_fn.eq.5)then
!     l=1; j=1.5d0
!   endif
! endif

  factor = l*(l+1) - j*(j+1) - 0.25d0
  do n=1, n_grid
     small(n) = ( large_deriv(n)  +  large(n)/r_grid(n) * factor )/2.d0/light_speed
     small_kinetic(n) = (eigenval-pot(n))/light_speed_sq + 2.d0
     small_kinetic(n) = small_kinetic(n)*small(n)
  enddo

 end subroutine


