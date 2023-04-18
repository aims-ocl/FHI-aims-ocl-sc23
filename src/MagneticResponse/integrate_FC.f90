!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Calculates the Fermi contact (FC) matrix elements in the
!!  nonrelativistic regime. Since only a single integration point is
!!  required, no parallelism is employed here.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
subroutine integrate_FC(i_atom, matrix_BC_out, timer)

  use aims_memory_tracking, only: aims_allocate
  use basis,                only: basis_atom, basis_fn, basis_wave_ordered
  use dimensions,           only: l_wave_max, n_basis, n_centers
  use geometry,             only: coords
  use MR_global
  use pbc_lists,            only: centers_basis_integrals, &
       & inv_centers_basis_integrals
  use scalapack_wrapper,    only: l_col, l_row, mxld, mxcol
  use species_data,         only: l_shell_max
  use timing,               only: start_timer, stop_timer
  use tools,                only: safe_deallocate
  use types,                only: dp

  use psi_at_nucleus_mod,   only: psi_at_nucleus

  implicit none

  ! atom at which the matrix elements are calculated
  integer, intent(in) :: i_atom
  ! block cyclic array that contains result of integration.
  real(dp), intent(out) :: matrix_BC_out(mxld,mxcol)
  real(dp), intent(in out), optional :: timer(4)

  character(*), parameter :: THIS_SUB = 'integrate_FC::'

  integer :: n_compute_atoms, n_zero_compute, n_batch_centers, n_compute, &
       & n_compute_fns

  integer, dimension(:), allocatable :: i_basis, i_basis_fns, i_atom_fns, &
       & atom_index, atom_index_inv, spline_array_start, spline_array_end, &
       & rad_index, wave_index, l_index, l_count, fn_atom, zero_index_point, &
       & batch_center

  integer, allocatable :: i_basis_fns_inv(:,:)

  real(dp), dimension(:), allocatable :: i_r, one_over_dist_tab, radial_wave
  real(dp), dimension(:,:), allocatable :: trigonom_tab, ylm_tab, dist_tab, &
       & dist_tab_sq, wave
  real(dp), dimension(:,:,:), allocatable :: dir_tab

  ! Wavefunction value at the nucleus
  real(dp) :: wave_at_p0(n_basis)
  ! Row and column iterators
  integer :: mb, nb
  ! Other
  integer :: nC

  ! STEP 1 - Allocations
  nC = n_centers
  nb = n_basis
  call aims_allocate(i_r, nC, THIS_SUB//'i_r')
  call aims_allocate(trigonom_tab, 4, nC, THIS_SUB//'trigonom_tab')
  call aims_allocate(one_over_dist_tab, nC, THIS_SUB//'one_over_dist_tab')
  call aims_allocate(dist_tab, n_centers, 1, THIS_SUB//'dist_tab')
  call aims_allocate(dist_tab_sq, nC, 1, THIS_SUB//'dist_tab_sq')
  call aims_allocate(dir_tab, 3, nC, 1, THIS_SUB//'dir_tab')
  call aims_allocate(i_basis_fns, nb*nC, THIS_SUB//'i_basis_fns')
  call aims_allocate(i_basis_fns_inv, nb, nC, THIS_SUB//'i_basis_fns_inv')
  call aims_allocate(i_atom_fns, nb*nC, THIS_SUB//'i_atom_fns')
  call aims_allocate(ylm_tab, (l_wave_max+1)**2, nC, THIS_SUB//'ylm_tab')
  call aims_allocate(radial_wave, nb, THIS_SUB//'radial_wave')
  call aims_allocate(wave, nb, 1, THIS_SUB//'wave')
  call aims_allocate(i_basis, nb, THIS_SUB//'i_basis')
  call aims_allocate(atom_index, nC, THIS_SUB//'atom_index')
  call aims_allocate(atom_index_inv, nC, THIS_SUB//'atom_index_inv')
  call aims_allocate(spline_array_start, nC, THIS_SUB//'spline_array_start')
  call aims_allocate(spline_array_end, nC, THIS_SUB//'spline_array_end')
  call aims_allocate(rad_index, nC, THIS_SUB//'rad_index')
  call aims_allocate(wave_index, nb, THIS_SUB//'wave_index')
  call aims_allocate(l_index, nb, THIS_SUB//'l_index')
  call aims_allocate(l_count, nb, THIS_SUB//'l_count')
  call aims_allocate(fn_atom, nb, THIS_SUB//'fn_atom')
  call aims_allocate(zero_index_point, nb, THIS_SUB//'zero_index_point')
  call aims_allocate(batch_center, nC, THIS_SUB//'batch_center')

  ! STEP 2 - Extract wavefunction values at the nucleus.
  if (present(timer)) call start_timer(timer)
  i_basis = 0
  n_compute = 0
  call tab_atom_centered_coords_p0(coords(:,i_atom), dist_tab_sq, dir_tab, &
       & n_centers, centers_basis_integrals)
  call prune_basis_p2(dist_tab_sq, n_compute, i_basis, n_basis, n_centers, &
       & inv_centers_basis_integrals)
  ! Here we remove zeros from dist_tab_sq in order to avoid numerical
  ! errors. These values are not used anyway.
  where (dist_tab_sq < tiny(1d0)) dist_tab_sq = 1d-1
  call collect_batch_centers_p2(n_compute, i_basis, n_basis, n_centers, &
       & inv_centers_basis_integrals, n_batch_centers, batch_center)
  i_basis_fns_inv = 0
  n_compute_atoms = 0
  n_compute_fns = 0
  call prune_radial_basis_p2(n_centers, n_basis, dist_tab_sq, dist_tab, &
       & dir_tab, n_compute_atoms, atom_index, atom_index_inv, n_compute_fns, &
       & i_basis_fns, i_basis_fns_inv, i_atom_fns, spline_array_start, &
       & spline_array_end, n_centers, centers_basis_integrals, n_compute, &
       & i_basis, n_batch_centers, batch_center, one_over_dist_tab, rad_index, &
       & wave_index, l_index, l_count, fn_atom, n_zero_compute, &
       & zero_index_point)
  call tab_local_geometry_p2(n_compute_atoms, atom_index, dist_tab, i_r)
  call tab_trigonom_p0(n_compute_atoms, dir_tab, trigonom_tab)
  call tab_wave_ylm_p0(n_compute_atoms, atom_index, trigonom_tab, l_shell_max, &
       & l_wave_max, ylm_tab)
  call evaluate_radial_functions_p0(spline_array_start, spline_array_end, &
       & n_compute_atoms, n_compute_fns, dist_tab, i_r, atom_index, &
       & i_basis_fns_inv, basis_wave_ordered, radial_wave, .false. , &
       & n_compute, n_basis)
  call evaluate_waves_p2(n_compute, n_compute_atoms, n_compute_fns, &
       & l_wave_max, ylm_tab, one_over_dist_tab, radial_wave, wave, rad_index, &
       & wave_index, l_index, l_count, fn_atom, n_zero_compute, &
       & zero_index_point)

  where (basis_atom(i_basis(:n_compute)) == i_atom)
     ! For orbitals originating from the current atom, extract the
     ! wavefunction value from psi_at_nucleus.
     wave_at_p0(:n_compute) = psi_at_nucleus(basis_fn(i_basis(:n_compute)))
  elsewhere
     ! For orbitals originating from other atoms, use the value of
     ! wave as constructed above.
     wave_at_p0(:n_compute) = wave(:n_compute,1)
  end where

  ! STEP 3 - Compute the matrix elements at the nucleus. Here each
  !          task has all the matrix elements. We can thus skip the
  !          step of writing the data into a packed matrix_packed and
  !          write it directly to a block cyclic matrix.
  matrix_BC_out = 0d0
  do mb = 1, n_compute
     do nb = 1, n_compute
        if (l_row(i_basis(mb)) /= 0 .and. l_col(i_basis(nb)) /= 0) &
             & matrix_BC_out(l_row(i_basis(mb)),l_col(i_basis(nb))) = &
             & wave_at_p0(mb)*wave_at_p0(nb)
     end do
  end do
  if (present(timer)) call stop_timer(timer, .true.)

  ! STEP 4 - Deallocations
  call safe_deallocate(i_r)
  call safe_deallocate(trigonom_tab)
  call safe_deallocate(one_over_dist_tab)
  call safe_deallocate(dist_tab)
  call safe_deallocate(dist_tab_sq)
  call safe_deallocate(dir_tab)
  call safe_deallocate(i_basis_fns)
  call safe_deallocate(i_basis_fns_inv)
  call safe_deallocate(i_atom_fns)
  call safe_deallocate(ylm_tab)
  call safe_deallocate(radial_wave)
  call safe_deallocate(wave)
  call safe_deallocate(i_basis)
  call safe_deallocate(atom_index)
  call safe_deallocate(atom_index_inv)
  call safe_deallocate(spline_array_start)
  call safe_deallocate(spline_array_end)
  call safe_deallocate(rad_index)
  call safe_deallocate(wave_index)
  call safe_deallocate(l_index)
  call safe_deallocate(l_count)
  call safe_deallocate(fn_atom)
  call safe_deallocate(zero_index_point)
  call safe_deallocate(batch_center)
end subroutine integrate_FC
