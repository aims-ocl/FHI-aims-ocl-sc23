#ifndef __SAVE_VAR_H__
#define __SAVE_VAR_H__

void m_save_load(const char *file_name, int is_load, int print);
void m_save_load_sumup(const char *file_name, int is_load, int print);
void m_save_load_rho(const char *file_name, int is_load, int print);
void m_save_load_H(const char *file_name, int is_load, int print);
void set_sum_up_param(int forces_on, double *partition_tab_std, double *delta_v_hartree, double *rho_multipole,
                      // double *centers_rho_multipole_spl, double *centers_delta_v_hart_part_spl,
                      double *adap_outer_radius_sq, double *multipole_radius_sq, int *l_hartree_max_far_distance,
                      double *outer_potential_radius, double *multipole_c);
void m_save_load_sumup_();
void m_save_check_sumup_(double *delta_v_hartree, double *rho_multipole);
void m_save_check_rho_(double *first_order_rho);
void m_save_check_h_(double *first_order_H, int* n_spin_, int* n_matrix_size_);

void debug_check_output_04_(double* array, int* dims, int* dims_num_);
void debug_check_output_04_file_(double* array, int* dims, int* dims_num_, const char* file_name);

void print_mem_info_();
void m_load_free();
void m_load_sumup_free();
#endif