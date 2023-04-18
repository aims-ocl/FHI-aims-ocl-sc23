#ifndef CUDA_MR_H
#define CUDA_MR_H

#include "gpuMacro.h"
#include <stdio.h>
#include <stdlib.h>

extern "C"
void mr_initialize_gpu(const int &wave_size,
                       const int &matrix_wave_size,
                       const int &matrix_batch_size,
                       const int &matrix_batch_GIAO_size,
                       const int &r_mn_basis_size,
                       const int &matrix_packed_size,
                       const int &i_basis_size,
                       const int &basis_glb_to_loc_size,
                       const int *basis_glb_to_loc,
                       const int &matrix_wave_size_f);

__global__
void d_zero_matrix(double *d_matrix, const int matrix_size);

extern "C"
void evaluate_mr_batch(const int &n_points,
                       const int &n_compute,
                       const double *wave,
                       const double *matrix_wave,
                       const int &i_symmetrization,
                       const int &i_dir,
                       const int &max_dims);

extern "C"
void evaluate_mr_batch_no_symm(const int &n_points,
                               const int &n_compute,
                               const double *wave,
                               const double *matrix_wave,
                               const int &i_symmetrization,
                               const int &i_dir,
                               const int &max_dims);

__global__
void symmetrize(double *d_matrix_batch, const int n_compute);

__global__
void antisymmetrize(double *d_matrix_batch, const int n_compute);

extern "C"
void get_matrix_packed_gpu(double *matrix_packed,
                           const int &matrix_packed_size);

extern "C"
void update_mr_batch_gpu(const int &n_compute,
                         const int &starting_point,
                         const int &i_basis);

__global__
void d_batch_to_packed(double *d_matrix_packed,
                       const double *d_matrix_batch,
                       const int n_compute,
                       const int starting_point,
                       const int *d_i_basis,
                       const int *d_basis_glb_to_loc);

extern "C"
void update_mr_batch_gpu_full(const int &n_compute,
                              const int &starting_point,
                              const int *i_basis,
                              const int &matrix_packed_ld);

__global__
void d_batch_to_packed_full(double *d_matrix_packed,
                            const double *d_matrix_batch,
                            const int n_compute,
                            const int starting_point,
                            const int *d_i_basis,
                            const int *d_basis_glb_to_loc,
                            const int matrix_packed_ld);

extern "C"
void symm_antisymm_gpu(const int &i_symmetrization, const int &n_compute);

// This function essentially the one in GIAO_mult_psi.f90, but
// translated to C++ and CUDA. See that source file for more comments
// about the structure of this function.
extern "C"
void giao_mult_psi_gpu(const bool &compact_indexing,
                       const bool &do_transpose,
                       const bool &mult_R_mn_both,
                       const int &i_dim,
                       const int &n_dims,
                       const int &n_compute,
                       const int &n_points,
                       const double *wave,
                       const double *matrix_wave,
                       const double *r_mn_basis,
                       double *matrix_batch_GIAO,
                       double *matrix_batch,
                       const int &n_spins);

// Perform operations on a single matrix
__global__
void d_operate1(const int n_compute,
                const double *d_r_mn_basis,
                const int r_GIAO_1,
                const int r_mn_1,
                const int r_GIAO_2,
                const int r_mn_2,
                double *matrix_batch,
                double *matrix_batch_GIAO);

// Perform operations on three matrices
__global__
void d_operate3(const int n_compute,
                const double *d_r_mn_basis,
                const int l_1,
                const int r_GIAO_11,
                const int r_mn_11,
                const int r_GIAO_12,
                const int r_mn_12,
                const int l_2,
                const int r_GIAO_21,
                const int r_mn_21,
                const int r_GIAO_22,
                const int r_mn_22,
                const int r_GIAO_31,
                const int r_mn_31,
                const int r_GIAO_32,
                const int r_mn_32,
                double *matrix_batch,
                double *matrix_batch_GIAO);

// Perform operations on four matrices
__global__
void d_operate4(const int n_compute,
                const double *d_r_mn_basis,
                const int l_1,
                const int r_GIAO_11,
                const int r_mn_11,
                const int r_GIAO_12,
                const int r_mn_12,
                const int l_2,
                const int r_GIAO_21,
                const int r_mn_21,
                const int r_GIAO_22,
                const int r_mn_22,
                const int l_3,
                const int r_GIAO_31,
                const int r_mn_31,
                const int r_GIAO_32,
                const int r_mn_32,
                const int r_GIAO_41,
                const int r_mn_41,
                const int r_GIAO_42,
                const int r_mn_42,
                double *matrix_batch,
                double *matrix_batch_GIAO);

__global__
void d_do_transpose(const int n_batches, const int n_compute, double *matrix);

__global__
void d_cross_product(const int n_compute,
                     const double *d_r_mn_basis,
                     const double *d_matrix_batch_aux,
                     const int shift1,
                     const int shift2,
                     double *d_matrix_batch);

__global__
void d_cross_product_modified(const int n_compute,
                              const double *d_r_mn_basis,
                              const double *d_matrix_batch_aux,
                              const int shift1,
                              const int shift2,
                              const bool negative,
                              double *d_matrix_batch);

extern "C"
void mr_destroy_gpu();

#endif
