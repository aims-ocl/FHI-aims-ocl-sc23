!****s* FHI-aims/evaluate_wavefunc_and_grad_real
!  NAME
!    evaluate_wavefunc_and_grad_real
!  SYNOPSIS

subroutine evaluate_wavefunc_and_grad_real(n_points, points, i_k_point, &
&                                 use_basbas, n_bas, &
&                                 n_func, coeffs, wavefunc, grad)

  !  PURPOSE
  !
  !    Evaluate real-valued wave function and its gradient defined by coeffs
  !    w.r.t. basis.f90 or prodbas.f90.
  !
  !
  !  USES

  use dimensions
  use pbc_lists
  use basis
  use prodbas
  use grids
  use mpi_tasks
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  integer, intent(IN) :: i_k_point
  integer, intent(IN) :: n_func
  logical, intent(IN) :: use_basbas
  integer, intent(IN) :: n_bas
  real*8, intent(IN) :: coeffs(n_bas, n_func)
  real*8, intent(OUT) :: wavefunc(n_points, n_func)
  real*8, intent(OUT) :: grad(3, n_points, n_func)

  !  INPUTS
  !    o n_points -- number of grid points
  !    o points -- grid points
  !    o i_k_point -- index of k-point (choose =1 for cluster)
  !    o use_basbas -- if (.true.), use prodbas.f90 instead of basis.f90.
  !    o n_bas -- either n_basis (use_basbas) or n_basbas (.not. use_basbas)
  !    o n_func -- number of functions
  !    o coeffs -- expansion coefficients
  !  OUTPUTS
  !    o wavefunc -- wave function values at points
  !    o grad -- wave function gradient values at points
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

  real*8, allocatable :: center2phase(:)
  integer, allocatable :: center2bas_off(:)
  integer :: i_center, i_cell, i_atom
  complex*16 :: dummy_cmplx
  integer :: info
  character(*), parameter :: func = 'evaluate_wavefunc_real'

  allocate(center2phase(n_centers), center2bas_off(n_centers), stat=info)
  call check_allocation(info, 'center2phase, center2bas_off', func)

  if (use_basbas) then
     if (n_bas /= n_basbas) call aims_stop('n_bas /= n_basbas', func)
  else
     if (n_bas /= n_basis) call aims_stop('n_bas /= n_basis', func)
  end if

  center2bas_off = -1  ! invalid
  do i_center = 1, n_centers
     i_cell = center_to_cell(i_center)
     if (i_cell == 0) cycle
     center2phase(i_center) = dble(k_phase(i_cell, i_k_point))
     i_atom = center_to_atom(i_center)
     if (use_basbas) then
        center2bas_off(i_center) = atom2basbas_off(i_atom)
     else
        center2bas_off(i_center) = atom2basis_off(i_atom)
     end if
  end do

  if (use_basbas) then
     call evaluate_wavefunc_and_grad(n_points, points, &
     &                      n_centers, coords_center, species_center, &
     &                      n_species, max_basbas_L, max_n_basbas_fnLsp, &
     &                      Lsp2n_basbas_fnLsp, Lsp2basbas_fn, Lsp2basbas_sp, &
     &                      4*n_max_grid, n_grid, r_grid_min, r_grid_inc, &
     &                      n_basbas_fns, charge_radius_basbas_fn, &
     &                                                basbas_wave_spl, .false.,&
     &                      center2bas_off, n_bas, n_func, &
     &                      .true., center2phase, dummy_cmplx, &
     &                              coeffs, dummy_cmplx, &
     &                              wavefunc, dummy_cmplx, grad, dummy_cmplx)
  else
     call evaluate_wavefunc_and_grad(n_points, points, &
     &                      n_centers, coords_center, species_center, &
     &                      n_species, max_basis_L, max_n_basis_fnLsp, &
     &                      Lsp2n_basis_fnLsp, Lsp2basis_fn, Lsp2basis_sp, &
     &                      4*n_max_grid, n_grid, r_grid_min, r_grid_inc, &
     &                      n_basis_fns, outer_radius, basis_wave_spl, .false.,&
     &                      center2bas_off, n_bas, n_func, &
     &                      .true., center2phase, dummy_cmplx, &
     &                              coeffs, dummy_cmplx, &
     &                              wavefunc, dummy_cmplx, grad, dummy_cmplx)
  end if

  deallocate(center2phase, center2bas_off)

end subroutine evaluate_wavefunc_and_grad_real
!******
