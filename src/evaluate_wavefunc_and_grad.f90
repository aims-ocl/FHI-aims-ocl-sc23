!****s* FHI-aims/evaluate_wavefunc_and_grad
!  NAME
!    evaluate_wavefunc
!  SYNOPSIS

subroutine evaluate_wavefunc_and_grad(n_points, points, &
&                            n_centers, center2coords, center2species, &
&                            n_species, max_bas_L, max_n_bas_fnLsp, &
&                            Lsp2n_bas_fnLsp, Lsp2bas_fn, Lsp2bas_sp, &
&                            n_wave_size, n_grid, r_grid_min, r_grid_inc, &
&                            n_bas_fns, basfn2radius, basfn2wave, use_bwave, &
&                            center2bas_off, n_bas, n_func, &
&                            is_real, center2phase_real, center2phase_cmplx, &
&                                     coeffs_real, coeffs_cmplx, &
&                                     wavefunc_real, wavefunc_cmplx, grad_real, grad_cmplx)

  !  PURPOSE
  !
  !     Evaluate a wave function and its gradients given by coefficients for given grid points.
  !
  !  USES

  use mpi_tasks
  use constants
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  integer, intent(IN) :: n_centers
  real*8, intent(IN) :: center2coords(3, n_centers)
  integer, intent(IN) :: center2species(n_centers)
  integer, intent(IN) :: n_species, max_bas_L, max_n_bas_fnLsp
  integer, intent(IN) :: Lsp2n_bas_fnLsp(0:max_bas_L, n_species)
  integer, intent(IN) :: Lsp2bas_fn(max_n_bas_fnLsp, 0:max_bas_L, n_species)
  integer, intent(IN) :: Lsp2bas_sp(max_n_bas_fnLsp, 0:max_bas_L, n_species)
  integer, intent(IN) :: n_wave_size
  integer, intent(IN) :: n_grid(n_species)
  real*8, intent(IN) :: r_grid_min(n_species), r_grid_inc(n_species)
  integer, intent(IN) :: n_bas_fns
  real*8, intent(IN) :: basfn2radius(n_bas_fns)
  real*8, intent(IN) :: basfn2wave(n_wave_size, n_bas_fns)
  logical, intent(IN) :: use_bwave
  integer, intent(IN) :: center2bas_off(n_centers)
  integer, intent(IN) :: n_bas
  integer, intent(IN) :: n_func
  logical, intent(IN) :: is_real
  real*8, intent(IN) :: center2phase_real(n_centers)
  real*8, intent(IN) :: coeffs_real(n_bas, n_func)
  real*8, intent(OUT) :: wavefunc_real(n_points, n_func)
  complex*16, intent(IN) :: center2phase_cmplx(n_centers)
  complex*16, intent(IN) :: coeffs_cmplx(n_bas, n_func)
  complex*16, intent(OUT) :: wavefunc_cmplx(n_points, n_func)
  real*8, intent(OUT) :: grad_real(3, n_points, n_func)
  complex*16, intent(OUT) :: grad_cmplx(3, n_points, n_func)

  !  INPUTS
  !   grid points
  !    o n_points -- number of grid points
  !    o points -- point coordinates
  !   center properties
  !    o n_centers -- Number of atomic centers
  !    o center2coords -- Where are they
  !    o center2species -- What are they
  !    o n_species -- Number of kinds of centers
  !    o center2phase -- Phase factor for all functions at that center
  !   basis properties
  !    o max_bas_L -- Highest angular momentum
  !    o max_n_bas_fnLsp -- Maximum number of radial parts in L-channel
  !    o Lsp2n_bas_fnLsp -- Number of radial parts i_fn in L-channel
  !    o Lsp2bas_fn -- (i_fn, L, i_species) -> i_bas_fn (see basis.f90)
  !    o Lsp2bas_sp -- (i_fn, L, i_species) -> i_bas_sp (see basis.f90)
  !   basfn properties
  !    o n_bas_fns -- Total number of radial parts
  !    o basfn2radius -- Radii of radial parts (when are they zero)
  !    o basfn2wave -- Actual radial part (u(r) = r*f(r))
  !    o use_bwave -- Is it B-spline?
  !   global function properties
  !    o is_real -- .true. if phases and coeffs are real-valued
  !    o n_bas -- Number of basis functions in coeffs.
  !    o center2bas_off -- Offset of center in coeffs
  !                        [Different centers can map to the same position,
  !                         possibly with differing center2phase.]
  !    o n_func -- Number of wave functions.
  !    o coeffs -- Expansion coefficients for wave functions in basis functions
  !  OUTPUTS
  !    o wave_func -- Values of wave function at points.
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_max_comp, n_comp
  integer, allocatable :: comp2center(:), comp2bas_sp(:)
  integer, allocatable :: comp2fn(:), comp2l(:), comp2m(:)
  real*8, allocatable :: wave(:,:), compcoeff_real(:,:), dwave(:,:,:)
  real*8, allocatable :: wf_real(:,:), wf_imag(:,:)
  complex*16, allocatable :: compcoeff_imag(:,:)
  integer :: i_comp, i_center, i_bas, i_point, i_func, bas_off
  integer :: info
  real*8 :: dummy(1)
  character(*), parameter :: func = 'evaluate_wavefunc_and_grad'
  real*8, allocatable :: dummy_radius(:)
  integer :: outer_radius_dimension

  !Create a dummy array for outer radius with large entries
  !This will prevent kinky cube files, without
  !explicitely changing outer_radius
  outer_radius_dimension=size(basfn2radius)
  allocate(dummy_radius(outer_radius_dimension), stat=info)
  call check_allocation(info, 'dummy_radius', func)
  dummy_radius=maxval(basfn2radius)


  ! --- Prune (get comp2...)

  ! Worst case: All functions on all centers have maximum angular momentum.
  ! Improbable.  But we are talking only about some indexing array sizes.
  n_max_comp = n_centers * max_n_bas_fnLsp * (max_bas_L+1) * (2*max_bas_L+1)

  allocate(comp2center(n_max_comp), comp2bas_sp(n_max_comp), stat=info)
  call check_allocation(info, 'comp2center, comp2bas_sp', func)
  allocate(comp2fn(n_max_comp), stat=info)
  call check_allocation(info, 'comp2fn', func)
  allocate(comp2l(n_max_comp), comp2m(n_max_comp), stat=info)
  call check_allocation(info, 'comp2l, comp2m', func)

  call prune_general_basis(n_points, points, &
  &                        n_centers, center2coords, center2species, &
  &                        n_species, max_bas_L, max_n_bas_fnLsp, n_max_comp, &
  &                        Lsp2n_bas_fnLsp, Lsp2bas_fn, Lsp2bas_sp, &
  &                        n_bas_fns, dummy_radius, &
  &                        n_comp, comp2center, comp2bas_sp, &
  &                        comp2fn, comp2l, comp2m)

  if (n_comp == 0) then
     if (is_real) then
        wavefunc_real = 0.d0
     else
        wavefunc_cmplx = 0.d0
     end if
     return
  end if
  

  ! --- Evaluate atomic orbitals

  allocate(wave(n_points, n_comp), stat=info)
  call check_allocation(info, 'wave', func)
  allocate(dwave(3, n_points, n_comp), stat=info)
  call check_allocation(info, 'dwave', func)

  call evaluate_waves_mult_point_center_fn(&
  &                  n_points, points, &
  &                  n_centers, center2coords, &
  &                  n_comp, comp2center, comp2fn, comp2l, comp2m, &
  &                  n_species, center2species, &
  &                  n_wave_size, n_grid, r_grid_min, r_grid_inc, &
  &                  n_bas_fns, basfn2wave, dummy_radius, &
  &                  use_bwave, .false., wave, .true., dwave)


  ! --- Get phased coefficients

  allocate(compcoeff_real(n_comp, n_func), stat=info)
  call check_allocation(info, 'compcoeff_real', func)
  compcoeff_real=0d0
  if (.not. is_real) then
     allocate(compcoeff_imag(n_comp, n_func), stat=info)
     call check_allocation(info, 'compcoeff_imag', func)
     compcoeff_imag=cmplx(0d0,0d0)
  end if


  do i_comp = 1, n_comp
     i_center = comp2center(i_comp)
     bas_off = center2bas_off(i_center)
     if (bas_off < 0) then
        call aims_stop('Assertion error: bas_off < 0', func)
        ! JW: If this happens, I'd bet on one of the points not being in the
        ! central parallel epipede.  That is, it has not been mapped back by
        ! map_to_center_cell().  Another cause would be that the basis
        ! functions' radius is larger than the outer_radius.  Why? The code in
        ! evaluate_wavefunc_*.f90 makes use of the fact that all possibly
        ! relevant centers (for points in the central parallel epipede and
        ! radii around centers not larger than outer_radius) are already
        ! defined in the center_to_cell list in pbc_lists.f90.  Solution: If
        ! the point has not been mapped back, do so.  If your radius is too
        ! large, you will have to rewrite the corresponding copy&paste clone
        ! of evaluate_wavefunc_*.f90. :-( If neither of this is the cause, I
        ! do not know either.
     end if
     i_bas = bas_off + comp2bas_sp(i_comp)
     if (is_real) then
        compcoeff_real(i_comp,:) &
        & = coeffs_real(i_bas, :)  * center2phase_real(i_center)
     else
        compcoeff_real(i_comp,:) &
        & = dble(coeffs_cmplx(i_bas, :)  * center2phase_cmplx(i_center))
        compcoeff_imag(i_comp,:) &
        & = aimag(coeffs_cmplx(i_bas, :) * center2phase_cmplx(i_center))
     end if
  end do

  ! --- Actually calculate wave functions

  ! wavefunc(i_point, i_func)
  ! := \sum_{i_comp} wave(i_point, i_comp) * compcoeff(i_comp, i_func)
  if (is_real) then
     ! wavefunc_real = matmul(wave, compcoeff_real)
     call dgemm('N', 'N', n_points, n_func, n_comp, &
     &          1.d0, wave, n_points, &
     &                compcoeff_real, n_comp, &
     &          0.d0, wavefunc_real, n_points)
     do i_comp = 1, 3
        call dgemm('N', 'N', n_points, n_func, n_comp, &
             &          1.d0, dwave(i_comp,:,:), n_points, &
             &                compcoeff_real, n_comp, &
             &          0.d0, grad_real(i_comp,:,:), n_points)
     end do
  else
     allocate(wf_real(n_points, n_func), wf_imag(n_points, n_func), stat=info)
     call check_allocation(info, 'wf_real, wf_imag', func)
     call dgemm('N', 'N', n_points, n_func, n_comp, &
     &          1.d0, wave, n_points, &
     &                compcoeff_real, n_comp, &
     &          0.d0, wf_real, n_points)
     call dgemm('N', 'N', n_points, n_func, n_comp, &
     &          1.d0, wave, n_points, &
     &                compcoeff_imag, n_comp, &
     &          0.d0, wf_imag, n_points)
     wavefunc_cmplx = cmplx(wf_real, wf_imag, kind(0.d0))
     do i_comp = 1, 3
        call dgemm('N', 'N', n_points, n_func, n_comp, &
             &          1.d0, dwave(i_comp,:,:), n_points, &
             &                compcoeff_real, n_comp, &
             &          0.d0, wf_real, n_points)
        call dgemm('N', 'N', n_points, n_func, n_comp, &
             &          1.d0, dwave(i_comp,:,:), n_points, &
             &                compcoeff_imag, n_comp, &
             &          0.d0, wf_imag, n_points)
        grad_cmplx(i_comp,:,:) = cmplx(wf_real, wf_imag, kind(0.d0))
     end do
     deallocate(wf_real, wf_imag)
  end if
  
  deallocate(comp2center, comp2bas_sp, comp2fn, comp2l, comp2m)
  deallocate(wave)
  deallocate(dwave)
  deallocate(compcoeff_real)
  !deallocate(dummy_radius)
  if (allocated(compcoeff_imag)) deallocate(compcoeff_imag)

end subroutine evaluate_wavefunc_and_grad
!******
