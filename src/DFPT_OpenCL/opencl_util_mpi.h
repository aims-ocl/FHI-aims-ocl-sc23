typedef struct OCL_UTIL_VARS_T {
  // dimensions
  int n_centers_hartree_potential;
  int n_periodic;
  int n_max_radial;
  int l_pot_max;
  int n_max_spline;
  int n_hartree_grid;
  int n_species;
  int n_atoms;
  int n_centers;
  int n_centers_basis_integrals;
  int n_centers_integrals;
  int n_max_compute_fns_ham;
  int n_basis_fns;
  int n_basis;
  int n_centers_basis_T;
  int n_centers_basis_I;
  int n_max_grid;
  int n_max_compute_atoms;
  int n_max_compute_ham;
  int n_max_compute_dens;
  int n_max_batch_size;
  // runtime_choices
  int use_hartree_non_periodic_ewald;
  int hartree_fp_function_splines;
  // int fast_ylm;
  // int new_ylm;
  int flag_rel;
  int Adams_Moulton_integrator;
  int compensate_multipole_errors;
  // pbc_lists
  int index_hamiltonian_dim2;
  int position_in_hamiltonian_dim1;
  int position_in_hamiltonian_dim2;
  int column_index_hamiltonian_size;
  // analytic_multipole_coefficients
  int l_max_analytic_multipole;
  // hartree_potential_real_p0
  int n_hartree_atoms;
  int hartree_force_l_add;
  // hartree_f_p_functions
  int Fp_max_grid;
  int lmax_Fp;
  // hartree_potential_storage
  int n_rho_multipole_atoms;
  // sumup batch
  int n_my_batches_work_sumup;
  int n_full_points_work_sumup;
  // rho batch
  int n_my_batches_work_rho;
  int n_full_points_work_rho;
  // h batch
  int n_my_batches_work_h;
  int n_full_points_work_h;
  // hartree_f_p_functions double
  double Fp_grid_min;
  double Fp_grid_inc;
  double Fp_grid_max;
  // sumup
  int forces_on;
  // help
  int n_cc_lm_ijk_l_pot_max;
  // ======= arrays ===========
  int *species;
  int *empty;
  int *centers_hartree_potential;
  int *center_to_atom;
  int *species_center;
  int *center_to_cell;
  int *cbasis_to_basis;
  int *cbasis_to_center;
  int *centers_basis_integrals;
  int *index_hamiltonian;
  int *position_in_hamiltonian;
  int *column_index_hamiltonian;
  double *coords_center;
  int *l_hartree;
  double *multipole_radius_free;
  int *n_grid;
  int *n_radial;
  double *r_grid_min;
  double *r_grid_inc;
  double *log_r_grid_inc;
  double *scale_radial;
  double *r_radial;
  double *r_grid;
  int *n_cc_lm_ijk;
  int *index_cc;
  int *index_ijk_max_cc;
  double *b0;
  double *b2;
  double *b4;
  double *b6;
  double *a_save;
  double *fp_function_spline;
  double *fpc_function_spline;
  double *ewald_radius_to;
  double *inv_ewald_radius_to;
  double *p_erfc_4;
  double *p_erfc_5;
  double *p_erfc_6;
  int *rho_multipole_index;
  int *compensation_norm;
  int *compensation_radius;
  // double *rho_multipole;
  int *perm_basis_fns_spl;
  double *outer_radius_sq;
  int *basis_fn;
  int *basis_l;
  double *atom_radius_sq;
  int *basis_fn_start_spl;
  int *basis_fn_atom;
  double *basis_wave_ordered;
  double *basis_kinetic_ordered;
  int *batches_size_sumup;
  double *batches_points_coords_sumup;
  // sumup arrays
  double *partition_tab;
  double *delta_v_hartree;
  double *rho_multipole;
  double *adap_outer_radius_sq;
  double *multipole_radius_sq;
  int *l_hartree_max_far_distance;
  double *outer_potential_radius;
  double *multipole_c;
  // sumup generate
  int* point_to_i_batch;
  int* point_to_i_index;
  int* valid_point_to_i_full_point;
  // ============ array_size =============
  int species___sizen;
  int empty___sizen;
  int centers_hartree_potential___sizen;
  int center_to_atom___sizen;
  int species_center___sizen;
  int center_to_cell___sizen;
  int cbasis_to_basis___sizen;
  int cbasis_to_center___sizen;
  int centers_basis_integrals___sizen;
  int index_hamiltonian___sizen;
  int position_in_hamiltonian___sizen;
  int column_index_hamiltonian___sizen;
  int coords_center___sizen;
  int l_hartree___sizen;
  int multipole_radius_free___sizen;
  int n_grid___sizen;
  int n_radial___sizen;
  int r_grid_min___sizen;
  int r_grid_inc___sizen;
  int log_r_grid_inc___sizen;
  int scale_radial___sizen;
  int r_radial___sizen;
  int r_grid___sizen;
  int n_cc_lm_ijk___sizen;
  int index_cc___sizen;
  int index_ijk_max_cc___sizen;
  int b0___sizen;
  int b2___sizen;
  int b4___sizen;
  int b6___sizen;
  int a_save___sizen;
  int fp_function_spline___sizen;
  int fpc_function_spline___sizen;
  int ewald_radius_to___sizen;
  int inv_ewald_radius_to___sizen;
  int p_erfc_4___sizen;
  int p_erfc_5___sizen;
  int p_erfc_6___sizen;
  int rho_multipole_index___sizen;
  int compensation_norm___sizen;
  int compensation_radius___sizen;
  // int rho_multipole___sizen;
  int perm_basis_fns_spl___sizen;
  int outer_radius_sq___sizen;
  int basis_fn___sizen;
  int basis_l___sizen;
  int atom_radius_sq___sizen;
  int basis_fn_start_spl___sizen;
  int basis_fn_atom___sizen;
  int basis_wave_ordered___sizen;
  int basis_kinetic_ordered___sizen;
  int batches_size_sumup___sizen;
  int batches_points_coords_sumup___sizen;
  // sumup arrays
  int partition_tab___sizen;
  int delta_v_hartree___sizen;
  int rho_multipole___sizen;
  int adap_outer_radius_sq___sizen;
  int multipole_radius_sq___sizen;
  int l_hartree_max_far_distance___sizen;
  int outer_potential_radius___sizen;
  int multipole_c___sizen;
  // sumup generate
  int point_to_i_batch___sizen;
  int point_to_i_index___sizen;
  int valid_point_to_i_full_point___sizen;
} OCL_UTIL_VARS;

extern OCL_UTIL_VARS ocl_util_vars_all[8];

void set_ocl_util_vars_single(
  OCL_UTIL_VARS* ocl_util_vars,
  // dimensions
  int n_centers_hartree_potential,
  int n_periodic,
  int n_max_radial,
  int l_pot_max,
  int n_max_spline,
  int n_hartree_grid,
  int n_species,
  int n_atoms,
  int n_centers,
  int n_centers_basis_integrals,
  int n_centers_integrals,
  int n_max_compute_fns_ham,
  int n_basis_fns,
  int n_basis,
  int n_centers_basis_T,
  int n_centers_basis_I,
  int n_max_grid,
  int n_max_compute_atoms,
  int n_max_compute_ham,
  int n_max_compute_dens,
  int n_max_batch_size,
  // runtime_choices
  int use_hartree_non_periodic_ewald,
  int hartree_fp_function_splines,
  // int fast_ylm,
  // int new_ylm,
  int flag_rel,
  int Adams_Moulton_integrator,
  int compensate_multipole_errors,
  // pbc_lists
  int index_hamiltonian_dim2,
  int position_in_hamiltonian_dim1,
  int position_in_hamiltonian_dim2,
  int column_index_hamiltonian_size,
  // analytic_multipole_coefficients
  int l_max_analytic_multipole,
  // hartree_potential_real_p0
  int n_hartree_atoms,
  int hartree_force_l_add,
  // hartree_f_p_functions
  int Fp_max_grid,
  int lmax_Fp,
  // hartree_potential_storage
  int n_rho_multipole_atoms,
  // sumup batch
  int n_my_batches_work_sumup,
  int n_full_points_work_sumup,
  // rho batch
  int n_my_batches_work_rho,
  int n_full_points_work_rho,
  // h batch
  int n_my_batches_work_h,
  int n_full_points_work_h,
  // hartree_f_p_functions double
  double Fp_grid_min,
  double Fp_grid_inc,
  double Fp_grid_max,
  // sumup
  int forces_on,
  // help
  int n_cc_lm_ijk_l_pot_max
);

void set_ocl_util_vars_arrays(
  OCL_UTIL_VARS* ocl_util_vars,
  int *batches_size_sumup,
  double *batches_points_coords_sumup,
  double *partition_tab,
  double *delta_v_hartree,
  double *rho_multipole,
  double *adap_outer_radius_sq,
  double *multipole_radius_sq,
  int *l_hartree_max_far_distance,
  double *outer_potential_radius,
  double *multipole_c
  // // generate
  // int* point_to_i_batch,
  // int* point_to_i_index,
  // int* valid_point_to_i_full_point
);

// void opencl_util_mpi_vars_default_(int* is_send, int* send_or_recv_rank);
void opencl_util_mpi_vars(OCL_UTIL_VARS* ocl_util_vars, int is_send, int send_or_recv_rank);
void opencl_util_mpi_arrays_(OCL_UTIL_VARS* ocl_util_vars, int is_send, int send_or_recv_rank);
void opencl_util_mpi_arrays_results_(OCL_UTIL_VARS* ocl_util_vars, int is_send, int send_or_recv_rank);
void init_opencl_util_mpi_();