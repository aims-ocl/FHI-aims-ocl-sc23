!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  C Fortran interfaces for the magnetic response part of GPU.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module gpuMR

  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int

  implicit none

  public

  interface
     subroutine mr_initialize_gpu(wave_size, matrix_wave_size, &
          & matrix_batch_size, matrix_batch_GIAO_size, r_mn_basis_size, &
          & matrix_packed_size, i_basis_size, basis_glb_to_loc_size, &
          & basis_glb_to_loc, matrix_wave_size_f) &
          & bind(c, name='mr_initialize_gpu')
       import c_int
       integer(c_int), intent(in) :: wave_size, matrix_wave_size, &
            & matrix_batch_size, matrix_batch_GIAO_size, r_mn_basis_size, &
            & matrix_packed_size, i_basis_size, basis_glb_to_loc_size, &
            & basis_glb_to_loc(*), matrix_wave_size_f
     end subroutine mr_initialize_gpu

     subroutine evaluate_mr_batch(n_points, n_compute, wave, matrix_wave, &
          & i_symmetrization, i_dir, max_dims) &
          & bind(c, name='evaluate_mr_batch')
       import c_double, c_int
       integer(c_int), intent(in) :: n_points, n_compute, i_symmetrization, &
            & i_dir, max_dims
       real(c_double), intent(in) :: wave(*), matrix_wave(*)
     end subroutine evaluate_mr_batch

     ! Like evaluate_mr_batch but without the
     ! symmetrization/antisymmetrization part
     subroutine evaluate_mr_batch_no_symm(n_points, n_compute, wave, &
          & matrix_wave, i_symmetrization, i_dim, max_dims) &
          & bind(c, name='evaluate_mr_batch_no_symm')
       import c_double, c_int
       integer(c_int), intent(in) :: n_points, n_compute, i_symmetrization, &
            & i_dim, max_dims
       real(c_double), intent(in) :: wave(*), matrix_wave(*)
     end subroutine evaluate_mr_batch_no_symm

     subroutine get_matrix_packed_gpu(matrix_packed, matrix_packed_size) &
          & bind(c, name='get_matrix_packed_gpu')
       import c_double, c_int
       real(c_double) :: matrix_packed(*)
       integer(c_int) :: matrix_packed_size
     end subroutine get_matrix_packed_gpu

     subroutine update_mr_batch_gpu(n_compute, i_shift, i_basis) &
          & bind(c, name='update_mr_batch_gpu')
       import c_int
       integer(c_int), intent(in) :: n_compute, i_shift, i_basis(*)
     end subroutine update_mr_batch_gpu

     ! Store both the lower and upper triangle in the packed array
     subroutine update_mr_batch_gpu_full(n_compute, i_shift, i_basis, &
          & matrix_packed_ld) &
          & bind(c, name='update_mr_batch_gpu_full')
       import c_int
       integer(c_int), intent(in) :: n_compute, i_shift, i_basis(*), &
            & matrix_packed_ld
     end subroutine update_mr_batch_gpu_full

     subroutine symm_antisymm_gpu(i_symmetrization, n_compute) &
          & bind(c, name='symm_antisymm_gpu')
       import c_int
       integer(c_int), intent(in) :: i_symmetrization, n_compute
     end subroutine symm_antisymm_gpu

     subroutine giao_mult_psi_gpu(c_compact_indexing, c_do_transpose, &
          & c_mult_R_mn_both, i_dim, n_dims, n_compute, n_points, wave, &
          & matrix_wave, r_mn_basis, matrix_batch_aux, matrix_batch, n_spins) &
          & bind(c, name='giao_mult_psi_gpu')
       import c_bool, c_double, c_int
       logical(c_bool), intent(in) :: c_compact_indexing, c_do_transpose, &
            & c_mult_R_mn_both
       integer(c_int), intent(in) :: i_dim, n_dims, n_compute, n_points
       real(c_double), intent(in) :: wave(*), matrix_wave(*), r_mn_basis(*)
       real(c_double), intent(in out) :: matrix_batch_aux(*), matrix_batch(*)
       integer(c_int), intent(in) :: n_spins
     end subroutine giao_mult_psi_gpu

     subroutine mr_destroy_gpu() bind(c, name='mr_destroy_gpu')
     end subroutine mr_destroy_gpu
  end interface
end module gpuMR
