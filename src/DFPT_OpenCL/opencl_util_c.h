#pragma once
#include <stdlib.h>

void loadProgramSource(const char **files, size_t length, char **buffer, size_t *sizes);

void opencl_init_();

void opencl_finish();

void opencl_common_buffer_init_();

void opencl_common_buffer_free_();

void sum_up_first_begin();

void sum_up_begin();

void sum_up_final_end();

void rho_pass_vars_(
  int* l_ylm_max,
  int* n_local_matrix_size,
  int* n_basis_local,
  int* perm_n_full_points,
  int* first_order_density_matrix_size,
  int* basis_l_max,
  int* n_points_all_batches,
  int* n_batch_centers_all_batches,
  int* batch_center_all_batches,
  int* batch_point_to_i_full_point,
  int* ins_idx_all_batches,
  double* first_order_rho,
  double* first_order_density_matrix,
  double* partition_tab
);

void rho_first_begin();
void H_first_begin();

void rho_begin();
void H_begin();

void rho_end();

void h_pass_vars_(
  int* j_coord_,
  int* n_spin_,
  int* l_ylm_max_,
  int* n_basis_local_,
  int* n_matrix_size_,
  int* basis_l_max,
  int* n_points_all_batches,
  int* n_batch_centers_all_batches,
  int* batch_center_all_batches,
  int* ins_idx_all_batches,
  int* batches_batch_i_basis_h,
  double* partition_all_batches,
  double* first_order_H,
  double* local_potential_parts_all_points,
  double* local_first_order_rho_all_batches,
  double* local_first_order_potential_all_batches,
  double* local_dVxc_drho_all_batches,
  double* local_rho_gradient,
  double* first_order_gradient_rho
);