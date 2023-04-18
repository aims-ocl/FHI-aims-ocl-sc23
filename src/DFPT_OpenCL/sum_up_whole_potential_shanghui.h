#pragma once

void tab_single_atom_centered_coords_p0_c_(int *current_center, double coord_current[3], double *dist_tab_sq,
                                           double dir_tab[3]);
void tab_single_atom_centered_coords_radial_log_p0_c_(int *current_center, double *dist_tab_sq, double dir_tab[3],
                                                      double *dist_tab_in, double *i_r, double *i_r_log,
                                                      double dir_tab_in[3]);
double invert_log_grid_c_(double r_current, double r_min, double scale);
double invert_log_grid_p2_c_(double r_current, int i_species);
double invert_radial_grid_c_(double r_current, int n_scale, double r_scale);
void tab_single_trigonom_p0_c_(double dir_tab[3], double trigonom_tab[4]);
void tab_single_wave_ylm_p2_c_(double trigonom_tab[4], int *l_max, int *l_ylm_max, double *ylm_tab);
// void increment_ylm_c_(double SINTH, double COSTH, double SINPH, double COSPH, double LMIN, double LMAX, double *Y);
void SHEval_c_(int lmax, double sintheta, double costheta, double sinphi, double cosphi, double *pSH);
void spline_vector_v2_c_(double *r_output, double *spl_param, int *n_l_dim, int *n_coeff, int *n_grid_dim,
                         int *n_points, int *n_vector, double *out_result);
void spline_vector_c_(double r_output, double *spl_param, int n_grid_dim, int n_l_dim, int n_points, int n_vector,
                      double *out_result);
void F_erf_c_(double *F, double r, int p, int c);
void F_erf_table_original_c_(double *F_erf_table, double r, int p_max);
void F_erfc_table_original_c_(double *F_table, double r, int p_max);
void far_distance_hartree_fp_periodic_single_atom_c_(int *current_atom, int *i_center, double *dist,
                                                     int *l_hartree_max_far_distance, // (n_atoms)
                                                     int *inside, int *forces_on, double *multipole_radius_sq,
                                                     double *adap_outer_radius
                                                     //  , int *non_peri_extd
);
void far_distance_hartree_fp_cluster_single_atom_p2_c_(double *dist_tab, int *l_max, int *forces_on);
void far_distance_real_hartree_potential_single_atom_p2_c_(int *i_center, double *potential, int *l_max,
                                                           double coord_current[3]);
void far_distance_real_hartree_potential_single_atom_c_(int *current_center, int *i_center, double *potential,
                                                        int *l_hartree_max_far_distance, double coord_current[3]);

void init_sum_up_c_(int *forces_on, double *partition_tab_std, double *delta_v_hartree, double *rho_multipole,
                    double *adap_outer_radius_sq, double *multipole_radius_sq, int *l_hartree_max_far_distance,
                    double *outer_potential_radius, double *multipole_c);
void release_sum_up_c_();
void sum_up_whole_potential_shanghui_sub_t_(int *forces_on, double *partition_tab_std, double *delta_v_hartree,
                                            double *rho_multipole, double *centers_rho_multipole_spl,
                                            // double *centers_delta_v_hart_part_spl, double *adap_outer_radius_sq,
                                            double *multipole_radius_sq, int *l_hartree_max_far_distance,
                                            double *outer_potential_radius, double *multipole_c);

void cubic_spline_v2_c_(double *spl_param, int *n_l_dim_, int *n_coeff_, int *n_grid_dim_, int *n_points_,
                        int *n_vector_);

void get_rho_multipole_spl_c_(double *rho_multipole_spl, int *spl_atom_);
double compensating_density_c_(double radius, double r_outer, int l);
void spline_angular_integral_log_c_(double *angular_part_spl, double *angular_integral_log, int *i_atom_);
void integrate_delta_v_hartree_internal_c_(double *angular_integral_log, double *delta_v_hartree, int *n_coeff_hartree_,
                                           int *i_atom_);
void integrate_delta_v_hartree_c_(double *angular_part_spl, double *delta_v_hartree, int *n_coeff_hartree_,
                                  int *i_atom_);