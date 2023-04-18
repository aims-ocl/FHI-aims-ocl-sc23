#include "save_load_var.h"
#include "pass_mod_var.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/sysinfo.h>

extern double *Fp_function_spline_slice;
extern double *Fpc_function_spline_slice;
extern double *Fp;

#define _CHK_(size1, size2)                                                                                            \
  if ((size1) != (size2)) {                                                                                            \
    printf("Error! size %lu != %lu (or %ld != %ld)\n", (size_t)(size1), (size_t)(size2), (long)(size1),                \
           (long)(size2));                                                                                             \
    printf("Error! rank%d, %s:%d, %s\n", myid, __FILE__, __LINE__, __FUNCTION__);                                     \
    fflush(stdout);                                                                                                    \
    exit(-1);                                                                                                          \
  }
#define _NL_ status = func(&endl, sizeof(char), 1, file_p);
#define _FWD_(type, var)                                                                                               \
  status = func(&var, sizeof(type), 1, file_p);                                                                        \
  if (print) {                                                                                                         \
    printf("rank%d, %s = %d\n", myid, #var, var);                                                                      \
  }
#define _FW_(type, var, size)                                                                                          \
  if (is_load) {                                                                                                       \
    var = (type *)malloc(sizeof(type) * (size));                                                                       \
  }                                                                                                                    \
  status = func(var, sizeof(type), size, file_p);                                                                      \
  _CHK_(size, status);
#define _FWS_(type, var, size)                                                                                         \
  status = func(var, sizeof(type), size, file_p);                                                                      \
  _CHK_(size, status);

void m_save_load(const char *file_name, int is_load, int print) {
  size_t (*func)(void *restrict, size_t, size_t, FILE *restrict);
  FILE *file_p;
  if (print)
    printf("rank%d, %s file <%s>\n", myid, is_load ? "read" : "write", file_name);

  if (is_load) {
    func = fread;
    file_p = fopen(file_name, "rb");
  } else {
    func = (size_t (*)(void *restrict, size_t, size_t, FILE *restrict))fwrite;
    file_p = fopen(file_name, "wb");
  }

  size_t status;
  // char endl = '\n';
  // dimensions
  _FWD_(int, n_centers_hartree_potential);
  _FWD_(int, n_periodic);
  _FWD_(int, n_max_radial);
  _FWD_(int, l_pot_max);
  _FWD_(int, n_max_spline);
  _FWD_(int, n_hartree_grid);
  _FWD_(int, n_species);
  _FWD_(int, n_atoms);
  _FWD_(int, n_centers);
  _FWD_(int, n_centers_basis_integrals);
  _FWD_(int, n_centers_integrals);
  _FWD_(int, n_max_compute_fns_ham);
  _FWD_(int, n_basis_fns);
  _FWD_(int, n_basis);
  _FWD_(int, n_centers_basis_T);
  _FWD_(int, n_centers_basis_I);
  _FWD_(int, n_max_grid);
  _FWD_(int, n_max_compute_atoms);
  _FWD_(int, n_max_compute_ham);
  _FWD_(int, n_max_compute_dens);
  _FWD_(int, n_max_batch_size);
  // runtime_choices
  _FWD_(int, use_hartree_non_periodic_ewald);
  _FWD_(int, hartree_fp_function_splines);
  // _FWD_(int, fast_ylm);
  // _FWD_(int, new_ylm);
  _FWD_(int, flag_rel);
  _FWD_(int, Adams_Moulton_integrator);
  _FWD_(int, compensate_multipole_errors);
  // pbc_lists
  _FWD_(int, index_hamiltonian_dim2);
  _FWD_(int, position_in_hamiltonian_dim1);
  _FWD_(int, position_in_hamiltonian_dim2);
  _FWD_(int, column_index_hamiltonian_size);
  // analytic_multipole_coefficients
  _FWD_(int, l_max_analytic_multipole);
  // hartree_potential_real_p0
  _FWD_(int, n_hartree_atoms);
  _FWD_(int, hartree_force_l_add);
  // hartree_f_p_functions
  _FWD_(int, Fp_max_grid);
  _FWD_(int, lmax_Fp);
  // hartree_potential_storage
  _FWD_(int, n_rho_multipole_atoms);
  // sumup batch
  _FWD_(int, n_my_batches_work_sumup);
  _FWD_(int, n_full_points_work_sumup);
  // rho batch
  _FWD_(int, n_my_batches_work_rho);
  _FWD_(int, n_full_points_work_rho);
  // h batch
  _FWD_(int, n_my_batches_work_h);
  _FWD_(int, n_full_points_work_h);

  status = func(&Fp_grid_min, sizeof(double), 1, file_p);
  status = func(&Fp_grid_inc, sizeof(double), 1, file_p);
  status = func(&Fp_grid_max, sizeof(double), 1, file_p);

  if (print)
    printf("rank%d, <save> 1 <constant var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(geometry, species), n_atoms);
  _FW_(int, MV(geometry, empty), n_atoms);
  if (print)
    printf("rank%d, <save> 2 <geometry var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(pbc_lists, centers_hartree_potential), n_centers_hartree_potential);
  _FW_(int, MV(pbc_lists, center_to_atom), n_centers);
  _FW_(int, MV(pbc_lists, species_center), n_centers);
  _FW_(int, MV(pbc_lists, center_to_cell), n_centers);
  _FW_(int, MV(pbc_lists, cbasis_to_basis), n_centers_basis_T);
  _FW_(int, MV(pbc_lists, cbasis_to_center), n_centers_basis_T);
  _FW_(int, MV(pbc_lists, centers_basis_integrals), n_centers_basis_integrals);
  _FW_(int, MV(pbc_lists, index_hamiltonian), 2 * index_hamiltonian_dim2 * n_basis);
  _FW_(int, MV(pbc_lists, position_in_hamiltonian), position_in_hamiltonian_dim1 *position_in_hamiltonian_dim2);
  _FW_(int, MV(pbc_lists, column_index_hamiltonian), column_index_hamiltonian_size);
  _FW_(double, MV(pbc_lists, coords_center), 3 * n_centers);
  if (print)
    printf("rank%d, <save> 3 <pbc_lists var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(species_data, l_hartree), n_species);
  _FW_(double, MV(species_data, multipole_radius_free), n_species);
  if (print)
    printf("rank%d, <save> 4 <species_data var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(grids, n_grid), n_species);
  _FW_(int, MV(grids, n_radial), n_species);
  _FW_(double, MV(grids, r_grid_min), n_species);
  _FW_(double, MV(grids, r_grid_inc), n_species);
  _FW_(double, MV(grids, log_r_grid_inc), n_species);
  _FW_(double, MV(grids, scale_radial), n_species);
  _FW_(double, MV(grids, r_radial), n_max_radial *n_species);
  _FW_(double, MV(grids, r_grid), n_max_grid *n_species);
  if (print)
    printf("rank%d, <save> 5 <grids var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(analytic_multipole_coefficients, n_cc_lm_ijk), (l_max_analytic_multipole + 1));
  _FW_(int, MV(analytic_multipole_coefficients, index_cc), n_cc_lm_ijk(l_max_analytic_multipole) * 6);
  _FW_(int, MV(analytic_multipole_coefficients, index_ijk_max_cc), 3 * (l_max_analytic_multipole + 1));
  if (print)
    printf("rank%d, <save> 6 <analytic_multipole_coefficients var> %ld bytes\n", myid, ftell(file_p));

  // _FW_(double, MV(hartree_potential_real_p0, multipole_c), n_cc_lm_ijk(l_pot_max) * n_atoms); // 迭代中会改变
  _FWS_(double, MV(hartree_potential_real_p0, b0), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, b2), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, b4), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, b6), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, a_save), pmaxab + 1);
  if (print)
    printf("rank%d, <save> 7 <hartree_potential_real_p0 var> %ld bytes\n", myid, ftell(file_p));

  _FW_(double, MV(hartree_f_p_functions, fp_function_spline), (lmax_Fp + 1) * n_max_spline * Fp_max_grid);
  _FW_(double, MV(hartree_f_p_functions, fpc_function_spline), (lmax_Fp + 1) * n_max_spline * Fp_max_grid);
  _FWS_(double, MV(hartree_f_p_functions, ewald_radius_to), 11);
  _FWS_(double, MV(hartree_f_p_functions, inv_ewald_radius_to), 2);
  _FWS_(double, MV(hartree_f_p_functions, p_erfc_4), 6);
  _FWS_(double, MV(hartree_f_p_functions, p_erfc_5), 7);
  _FWS_(double, MV(hartree_f_p_functions, p_erfc_6), 8);
  if (print)
    printf("rank%d, <save> 8 <hartree_f_p_functions var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(hartree_potential_storage, rho_multipole_index), n_atoms);
  if (compensate_multipole_errors) {
    _FW_(int, MV(hartree_potential_storage, compensation_norm), n_atoms);
    _FW_(int, MV(hartree_potential_storage, compensation_radius), n_atoms);
  }
  if (MV(hartree_potential_storage, use_rho_multipole_shmem)){
    _FW_(double, MV(hartree_potential_storage, rho_multipole_shmem_ptr),
         (l_pot_max + 1) * (l_pot_max + 1) * (n_max_radial + 2) * n_rho_multipole_atoms);
  } else {
    _FW_(double, MV(hartree_potential_storage, rho_multipole),
        (l_pot_max + 1) * (l_pot_max + 1) * (n_max_radial + 2) * n_rho_multipole_atoms);
  }
  if (print)
    printf("rank%d, <save> 8 <hartree_potential_storage var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(basis, perm_basis_fns_spl), n_basis_fns);
  _FW_(double, MV(basis, outer_radius_sq), n_basis_fns);
  _FW_(int, MV(basis, basis_fn), n_basis);
  _FW_(int, MV(basis, basis_l), n_basis);
  _FW_(double, MV(basis, atom_radius_sq), n_species);
  _FW_(int, MV(basis, basis_fn_start_spl), n_species);
  _FW_(int, MV(basis, basis_fn_atom), n_basis_fns *n_atoms);
  _FW_(double, MV(basis, basis_wave_ordered), n_basis_fns *n_max_spline *n_max_grid);
  _FW_(double, MV(basis, basis_kinetic_ordered), n_basis_fns *n_max_spline *n_max_grid);
  if (print)
    printf("rank%d, <save> 8 <basis var> %ld bytes\n", myid, ftell(file_p));

  // sumup batch
  _FW_(int, MV(opencl_util, batches_size_sumup), n_my_batches_work_sumup); // 进程间有差别
  _FW_(double, MV(opencl_util, batches_points_coords_sumup),
       3 * n_max_batch_size * n_my_batches_work_sumup); // 进程间有差别
  if (print)
    printf("rank%d, <save> 8 <opencl_util sumup var> %ld bytes\n", myid, ftell(file_p));

  if (print)
    printf("rank%d, <save> all %ld bytes\n", myid, ftell(file_p));
  fclose(file_p);
}

void m_save_load_check(const char *file_name, int is_load, int print) {
  size_t (*func)(void *restrict, size_t, size_t, FILE *restrict);
  FILE *file_p;
  if (print)
    printf("rank%d, %s file <%s>\n", myid, is_load ? "read" : "write", file_name);

  if (is_load) {
    func = fread;
    file_p = fopen(file_name, "rb");
  } else {
    func = (size_t (*)(void *restrict, size_t, size_t, FILE *restrict))fwrite;
    file_p = fopen(file_name, "wb");
  }

  size_t status;
  // char endl = '\n';
  // dimensions
  _FWD_(int, n_centers_hartree_potential);
  _FWD_(int, n_periodic);
  _FWD_(int, n_max_radial);
  _FWD_(int, l_pot_max);
  _FWD_(int, n_max_spline);
  _FWD_(int, n_hartree_grid);
  _FWD_(int, n_species);
  _FWD_(int, n_atoms);
  _FWD_(int, n_centers);
  _FWD_(int, n_centers_basis_integrals);
  _FWD_(int, n_centers_integrals);
  _FWD_(int, n_max_compute_fns_ham);
  _FWD_(int, n_basis_fns);
  _FWD_(int, n_basis);
  _FWD_(int, n_centers_basis_T);
  _FWD_(int, n_centers_basis_I);
  _FWD_(int, n_max_grid);
  _FWD_(int, n_max_compute_atoms);
  _FWD_(int, n_max_compute_ham);
  _FWD_(int, n_max_compute_dens);
  _FWD_(int, n_max_batch_size);
  // runtime_choices
  _FWD_(int, use_hartree_non_periodic_ewald);
  _FWD_(int, hartree_fp_function_splines);
  // _FWD_(int, fast_ylm);
  // _FWD_(int, new_ylm);
  _FWD_(int, flag_rel);
  _FWD_(int, Adams_Moulton_integrator);
  _FWD_(int, compensate_multipole_errors);
  // pbc_lists
  // _FWD_(int, index_hamiltonian_dim2);
  // _FWD_(int, position_in_hamiltonian_dim1);
  // _FWD_(int, position_in_hamiltonian_dim2);
  // -_FWD_(int, column_index_hamiltonian_size);
  // analytic_multipole_coefficients
  _FWD_(int, l_max_analytic_multipole);
  // hartree_potential_real_p0
  _FWD_(int, n_hartree_atoms);
  _FWD_(int, hartree_force_l_add);
  // // hartree_f_p_functions
  // _FWD_(int, Fp_max_grid);
  // _FWD_(int, lmax_Fp);
  // // hartree_potential_storage
  // _FWD_(int, n_rho_multipole_atoms);
  // // sumup batch
  // _FWD_(int, n_my_batches_work_sumup);
  // _FWD_(int, n_full_points_work_sumup);
  // // rho batch
  // _FWD_(int, n_my_batches_work_rho);
  // _FWD_(int, n_full_points_work_rho);
  // // h batch
  // _FWD_(int, n_my_batches_work_h);
  // _FWD_(int, n_full_points_work_h);

  // status = func(&Fp_grid_min, sizeof(double), 1, file_p);
  // status = func(&Fp_grid_inc, sizeof(double), 1, file_p);
  // status = func(&Fp_grid_max, sizeof(double), 1, file_p);

  if (print)
    printf("rank%d, <save> 1 <constant var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(geometry, species), n_atoms);
  _FW_(int, MV(geometry, empty), n_atoms);
  if (print)
    printf("rank%d, <save> 2 <geometry var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(pbc_lists, centers_hartree_potential), n_centers_hartree_potential);
  _FW_(int, MV(pbc_lists, center_to_atom), n_centers);
  _FW_(int, MV(pbc_lists, species_center), n_centers);
  _FW_(int, MV(pbc_lists, center_to_cell), n_centers);
  _FW_(int, MV(pbc_lists, cbasis_to_basis), n_centers_basis_T);
  _FW_(int, MV(pbc_lists, cbasis_to_center), n_centers_basis_T);
  _FW_(int, MV(pbc_lists, centers_basis_integrals), n_centers_basis_integrals);
  // _FW_(int, MV(pbc_lists, index_hamiltonian), 2 * index_hamiltonian_dim2 * n_basis);
  // _FW_(int, MV(pbc_lists, position_in_hamiltonian), position_in_hamiltonian_dim1 *position_in_hamiltonian_dim2);
  // _FW_(int, MV(pbc_lists, column_index_hamiltonian), column_index_hamiltonian_size);
  _FW_(double, MV(pbc_lists, coords_center), 3 * n_centers);
  if (print)
    printf("rank%d, <save> 3 <pbc_lists var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(species_data, l_hartree), n_species);
  _FW_(double, MV(species_data, multipole_radius_free), n_species);
  if (print)
    printf("rank%d, <save> 4 <species_data var> %ld bytes\n", myid, ftell(file_p));


  for(int i=0; i<n_species; i++){
    printf("n_radial[%d] = %d\n", i, MV(grids, n_radial)[i]);
    for(int j=MV(grids, n_radial)[i]; j<n_max_radial; j++)
      MV(grids, r_radial)[i*n_max_radial+j] = 0;
  }
  for(int i=0; i<n_species; i++){
    printf("n_grid[%d] = %d\n", i, MV(grids, n_grid)[i]);
    for(int j=MV(grids, n_grid)[i]; j<n_max_grid; j++){
      MV(grids, r_grid)[i*n_max_grid+j] = 0;
    }
  }
  _FW_(int, MV(grids, n_grid), n_species);
  _FW_(int, MV(grids, n_radial), n_species);
  _FW_(double, MV(grids, r_grid_min), n_species);
  _FW_(double, MV(grids, r_grid_inc), n_species);
  _FW_(double, MV(grids, log_r_grid_inc), n_species);
  _FW_(double, MV(grids, scale_radial), n_species);
  _FW_(double, MV(grids, r_radial), n_max_radial * n_species);
  _FW_(double, MV(grids, r_grid), n_max_grid *n_species);
  if (print)
    printf("rank%d, <save> 5 <grids var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(analytic_multipole_coefficients, n_cc_lm_ijk), (l_max_analytic_multipole + 1));
  _FW_(int, MV(analytic_multipole_coefficients, index_cc), n_cc_lm_ijk(l_max_analytic_multipole) * 6);
  _FW_(int, MV(analytic_multipole_coefficients, index_ijk_max_cc), 3 * (l_max_analytic_multipole + 1));
  if (print)
    printf("rank%d, <save> 6 <analytic_multipole_coefficients var> %ld bytes\n", myid, ftell(file_p));

  // _FW_(double, MV(hartree_potential_real_p0, multipole_c), n_cc_lm_ijk(l_pot_max) * n_atoms); // 迭代中会改变
  _FWS_(double, MV(hartree_potential_real_p0, b0), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, b2), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, b4), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, b6), pmaxab + 1);
  _FWS_(double, MV(hartree_potential_real_p0, a_save), pmaxab + 1);
  if (print)
    printf("rank%d, <save> 7 <hartree_potential_real_p0 var> %ld bytes\n", myid, ftell(file_p));

  // _FW_(double, MV(hartree_f_p_functions, fp_function_spline), (lmax_Fp + 1) * n_max_spline * Fp_max_grid);
  // _FW_(double, MV(hartree_f_p_functions, fpc_function_spline), (lmax_Fp + 1) * n_max_spline * Fp_max_grid);
  // _FWS_(double, MV(hartree_f_p_functions, ewald_radius_to), 11);
  // _FWS_(double, MV(hartree_f_p_functions, inv_ewald_radius_to), 2);
  // _FWS_(double, MV(hartree_f_p_functions, p_erfc_4), 6);
  // _FWS_(double, MV(hartree_f_p_functions, p_erfc_5), 7);
  // _FWS_(double, MV(hartree_f_p_functions, p_erfc_6), 8);
  if (print)
    printf("rank%d, <save> 8 <hartree_f_p_functions var> %ld bytes\n", myid, ftell(file_p));

  // _FW_(int, MV(hartree_potential_storage, rho_multipole_index), n_atoms);
  if (compensate_multipole_errors) {
    _FW_(int, MV(hartree_potential_storage, compensation_norm), n_atoms);
    _FW_(int, MV(hartree_potential_storage, compensation_radius), n_atoms);
  }
  // _FW_(double, MV(hartree_potential_storage, rho_multipole_shmem_ptr),
  //      (l_pot_max + 1) * (l_pot_max + 1) * (n_max_radial + 2) * n_rho_multipole_atoms);
  if (print)
    printf("rank%d, <save> 8 <hartree_potential_storage var> %ld bytes\n", myid, ftell(file_p));

  _FW_(int, MV(basis, perm_basis_fns_spl), n_basis_fns);
  _FW_(double, MV(basis, outer_radius_sq), n_basis_fns);
  _FW_(int, MV(basis, basis_fn), n_basis);
  _FW_(int, MV(basis, basis_l), n_basis);
  _FW_(double, MV(basis, atom_radius_sq), n_species);
  _FW_(int, MV(basis, basis_fn_start_spl), n_species);
  _FW_(int, MV(basis, basis_fn_atom), n_basis_fns *n_atoms);
  // _FW_(double, MV(basis, basis_wave_ordered), n_basis_fns *n_max_spline *n_max_grid);
  // _FW_(double, MV(basis, basis_kinetic_ordered), n_basis_fns *n_max_spline *n_max_grid);
  if (print)
    printf("rank%d, <save> 8 <basis var> %ld bytes\n", myid, ftell(file_p));

  // // sumup batch
  // _FW_(int, MV(opencl_util, batches_size_sumup), n_my_batches_work_sumup); // 进程间有差别
  // _FW_(double, MV(opencl_util, batches_points_coords_sumup),
  //      3 * n_max_batch_size * n_my_batches_work_sumup); // 进程间有差别
  if (print)
    printf("rank%d, <save> 8 <opencl_util sumup var> %ld bytes\n", myid, ftell(file_p));

  if (print)
    printf("rank%d, <save> all %ld bytes\n", myid, ftell(file_p));
  fclose(file_p);
}

void m_load_free() {
  free(MV(geometry, species));
  free(MV(geometry, empty));

  free(MV(pbc_lists, centers_hartree_potential));
  free(MV(pbc_lists, center_to_atom));
  free(MV(pbc_lists, species_center));
  free(MV(pbc_lists, center_to_cell));
  free(MV(pbc_lists, cbasis_to_basis));
  free(MV(pbc_lists, cbasis_to_center));
  free(MV(pbc_lists, centers_basis_integrals));
  free(MV(pbc_lists, index_hamiltonian));
  free(MV(pbc_lists, position_in_hamiltonian));
  free(MV(pbc_lists, column_index_hamiltonian));
  free(MV(pbc_lists, coords_center));

  free(MV(species_data, l_hartree));
  free(MV(species_data, multipole_radius_free));

  free(MV(grids, n_grid));
  free(MV(grids, n_radial));
  free(MV(grids, r_grid_min));
  free(MV(grids, r_grid_inc));
  free(MV(grids, log_r_grid_inc));
  free(MV(grids, scale_radial));
  free(MV(grids, r_radial));
  free(MV(grids, r_grid));

  free(MV(analytic_multipole_coefficients, n_cc_lm_ijk));
  free(MV(analytic_multipole_coefficients, index_cc));
  free(MV(analytic_multipole_coefficients, index_ijk_max_cc));

  free(MV(hartree_f_p_functions, fp_function_spline));
  free(MV(hartree_f_p_functions, fpc_function_spline));

  free(MV(hartree_potential_storage, rho_multipole_index));
  if (compensate_multipole_errors) {
    free(MV(hartree_potential_storage, compensation_norm));
    free(MV(hartree_potential_storage, compensation_radius));
  }
  free(MV(hartree_potential_storage, rho_multipole_shmem_ptr));

  free(MV(basis, perm_basis_fns_spl));
  free(MV(basis, outer_radius_sq));
  free(MV(basis, basis_fn));
  free(MV(basis, basis_l));
  free(MV(basis, atom_radius_sq));
  free(MV(basis, basis_fn_start_spl));
  free(MV(basis, basis_fn_atom));
  free(MV(basis, basis_wave_ordered));
  free(MV(basis, basis_kinetic_ordered));

  // sumup batch
  free(MV(opencl_util, batches_size_sumup));
  free(MV(opencl_util, batches_points_coords_sumup));
  // // rho batch
  // free(MV(opencl_util, batches_size_rho));
  // free(MV(opencl_util, batches_batch_n_compute_rho));
  // free(MV(opencl_util, batches_batch_i_basis_rho));
  // free(MV(opencl_util, batches_points_coords_rho));
  // // H batch
  // free(MV(opencl_util, batches_size_h));
  // free(MV(opencl_util, batches_batch_n_compute_h));
  // free(MV(opencl_util, batches_batch_i_basis_h));
  // free(MV(opencl_util, batches_points_coords_h));
}

static int m_save_load_count = 0;
void m_save_load_() {
  int count = m_save_load_count;
  m_save_load_count++;
  // if (count == 1 || myid == 0) {
  if ((myid == 0 || myid == 1 || myid == 3) && count <= 1) {
    char save_file_name[64];
    sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid, count);
    m_save_load(save_file_name, 0, 1);
  }
}

void m_save_load_not_count_() {
  int count = m_save_load_count;
  if ((myid == 0 || myid == 1 || myid == 3) && count <= 1) {
    char save_file_name[64];
    sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid, count);
    m_save_load(save_file_name, 0, 1);
  }
}

static int m_save_load_sumup_count = 0;
void m_save_load_sumup_() {
  int count = m_save_load_sumup_count;
  m_save_load_sumup_count++;
  // if (count == 1 || myid == 0) {
  if ((myid == 0 || myid == 3) && count <= 1) {
    char save_file_name[64];
    sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid, count);
    m_save_load_sumup(save_file_name, 0, 0);
  }
}

void m_save_load_sumup(const char *file_name, int is_load, int print) {
  size_t (*func)(void *restrict, size_t, size_t, FILE *restrict);
  FILE *file_p;
  if (is_load) {
    func = fread;
  } else {
    func = (size_t (*)(void *restrict, size_t, size_t, FILE *restrict))fwrite;
  }

  size_t status;

  char more_file_name[128];
  const char *tmp = file_name;
  const char *last_begin = NULL;
  while (*tmp != '\0') {
    if (*tmp == '/')
      last_begin = tmp;
    ++tmp;
  }
  if (last_begin != NULL) {
    last_begin++;
    memcpy(more_file_name, file_name, last_begin - file_name);
    sprintf(more_file_name + (last_begin - file_name), "sumup_%s", last_begin);
  } else
    sprintf(more_file_name, "sumup_%s", file_name);
  if (print)
    printf("%s\n", more_file_name);
  if (is_load) {
    file_p = fopen(more_file_name, "rb");
  } else {
    file_p = fopen(more_file_name, "wb");
  }
  _FWD_(int, sum_up_param.forces_on);
  _FW_(double, sum_up_param.partition_tab, n_full_points_work_sumup);
  _FW_(double, sum_up_param.delta_v_hartree, n_full_points_work_sumup);
  _FW_(double, sum_up_param.rho_multipole, n_full_points_work_sumup);
  if (print)
    printf("rank%d, %ld bytes, mdim=%d\n", myid, ftell(file_p), n_cc_lm_ijk(l_pot_max) * n_atoms);
  // _FW_(double, sum_up_param.centers_rho_multipole_spl,
  //      (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2) * n_atoms);
  // _FW_(double, sum_up_param.centers_delta_v_hart_part_spl,
  //      (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * n_atoms);
  // if (print)
  //   printf("rank%d, %ld bytes, mdim=%d\n", myid, ftell(file_p), n_cc_lm_ijk(l_pot_max) * n_atoms);
  _FW_(double, sum_up_param.adap_outer_radius_sq, n_atoms);
  _FW_(double, sum_up_param.multipole_radius_sq, n_atoms);
  _FW_(int, sum_up_param.l_hartree_max_far_distance, n_atoms);
  _FW_(double, sum_up_param.outer_potential_radius, (l_pot_max + 1) * n_atoms);
  if (print)
    printf("rank%d, %ld bytes, mdim=%d\n", myid, ftell(file_p), n_cc_lm_ijk(l_pot_max) * n_atoms);
  // _FW_(double, sum_up_param.multipole_c, n_cc_lm_ijk(l_pot_max) * n_atoms);
  if (is_load) {
    sum_up_param.multipole_c = (double *)malloc(
        sizeof(double) *
        (MV(analytic_multipole_coefficients, n_cc_lm_ijk)[MV(dimensions, l_pot_max)] * MV(dimensions, n_atoms)));
  }
  status =
      func(sum_up_param.multipole_c, sizeof(double),
           MV(analytic_multipole_coefficients, n_cc_lm_ijk)[MV(dimensions, l_pot_max)] * MV(dimensions, n_atoms),
           file_p);
  _CHK_(MV(analytic_multipole_coefficients, n_cc_lm_ijk)[MV(dimensions, l_pot_max)] * MV(dimensions, n_atoms),
        status);
  fclose(file_p);
}

void m_load_sumup_free() {
  free(sum_up_param.partition_tab);
  free(sum_up_param.delta_v_hartree);
  free(sum_up_param.rho_multipole);

  free(sum_up_param.adap_outer_radius_sq);
  free(sum_up_param.multipole_radius_sq);
  free(sum_up_param.l_hartree_max_far_distance);
  free(sum_up_param.outer_potential_radius);

  free(sum_up_param.multipole_c);
}

static int m_save_check_sumup_count = -1;
void m_save_check_sumup_(double *delta_v_hartree, double *rho_multipole) {
  char save_file_name[64];
  sprintf(save_file_name, "sumup_check_rank%d_%d.bin", myid, m_save_check_sumup_count++);
  if(!(m_save_check_sumup_count <= 1)){
    return;
  }
  FILE *file_p = fopen(save_file_name, "w");
  for (int i = 0; i < n_full_points_work_sumup; i++) {
    fprintf(file_p, "%6d, %.13lf, %.13lf\n", i, delta_v_hartree[i], rho_multipole[i]);
  }
  fclose(file_p);
}

static int m_save_check_rho_test_count = -1;
void m_save_check_rho_test_(int *j_cells, 
  int* i_place_begins,
  int* i_place_ends,
  int* i_basis_uc,
  int* j_basis_uc,
  int* poss,
  int* n_compute_c) {
  char save_file_name[64];
  sprintf(save_file_name, "rho_check_test_rank%d_%d.bin", myid, m_save_check_rho_test_count++);
  FILE *file_p = fopen(save_file_name, "w");
  for (int i = 0; i < n_max_batch_size; i++) {
    for (int j = 0; j < *n_compute_c; j++) {
      fprintf(file_p, "(%6d, %6d), %6d, %6d, %6d, %6d, %6d, %6d\n", i, j, 
        j_cells[i * *n_compute_c + j], i_place_begins[i * *n_compute_c + j], i_place_ends[i * *n_compute_c + j],
              i_basis_uc[i * *n_compute_c + j], j_basis_uc[i * *n_compute_c + j], poss[i * *n_compute_c + j]);
    }
  }
  fclose(file_p);
}

static int m_save_check_rho_wave_count = -1;
void m_save_check_rho_wave_(double *wave, int* n_compute_c) {
  char save_file_name[64];
  sprintf(save_file_name, "rho_check_wave_rank%d_%d.bin", myid, m_save_check_rho_wave_count++);
  FILE *file_p = fopen(save_file_name, "w");
  // for (int i = 0; i < *n_compute_c; i++) {
  //   for (int j = 0; j < *n_compute_c; j++) {
  //     fprintf(file_p, "%6d, %6d, %.13lf\n", i, j, wave[i * *n_compute_c + j]);
  //   }
  // }
  for (int i = 0; i < n_max_batch_size; i++) {
    for (int j = 0; j < *n_compute_c; j++) {
      fprintf(file_p, "%6d, %6d, %.13lf\n", i, j, wave[i * n_max_compute_ham + j]);
    }
  }
  fclose(file_p);
}

static int m_save_check_rho_count = -1;
void m_save_check_rho_(double *first_order_rho) {
  char save_file_name[64];
  sprintf(save_file_name, "rho_check_rank%d_%d.bin", myid, m_save_check_rho_count++);
  if(!(m_save_check_rho_count <= 1)){
    return;
  }
  FILE *file_p = fopen(save_file_name, "w");
  int size = -1;
  if(rho_param.n_basis_local > 0){
    size = rho_param.perm_n_full_points;
  } else {
    size = n_full_points_work_rho;
  }
  for (int i = 0; i < size; i++) {
    fprintf(file_p, "%6d, %.13lf\n", i, first_order_rho[i]);
  }
  fclose(file_p);
}

static int m_save_check_H_count = -1;
void m_save_check_h_(double *first_order_H, int* n_spin_, int* n_matrix_size_) {
  char save_file_name[64];
  sprintf(save_file_name, "H_check_rank%d_%d.bin", myid, m_save_check_H_count++);
  if(!(m_save_check_H_count <= 4)){
    return;
  }
  FILE *file_p = fopen(save_file_name, "w");
  for(int i=0; i<*n_spin_; i++){
    for(int j=0; j<*n_matrix_size_; j++){
      fprintf(file_p, "%6d, %6d, %.13lf\n", i, j, first_order_H[i * (*n_matrix_size_) + j]);
    }
  }
  fclose(file_p);
}

void m_save_load_rho(const char *file_name, int is_load, int print) {
  size_t (*func)(void *restrict, size_t, size_t, FILE *restrict);
  FILE *file_p;
  if (is_load) {
    func = fread;
  } else {
    func = (size_t (*)(void *restrict, size_t, size_t, FILE *restrict))fwrite;
  }
  size_t status;
  char more_file_name[128];
  const char *tmp = file_name;
  const char *last_begin = NULL;
  while (*tmp != '\0') {
    if (*tmp == '/')
      last_begin = tmp;
    ++tmp;
  }
  if (last_begin != NULL) {
    last_begin++;
    memcpy(more_file_name, file_name, last_begin - file_name);
    sprintf(more_file_name + (last_begin - file_name), "rho_%s", last_begin);
  } else
    sprintf(more_file_name, "rho_%s", file_name);
  printf("%s\n", more_file_name);
  if (is_load) {
    file_p = fopen(more_file_name, "rb");
  } else {
    file_p = fopen(more_file_name, "wb");
  }
  _FWD_(int, rho_param.l_ylm_max);
  _FWD_(int, rho_param.n_local_matrix_size);
  _FWD_(int, rho_param.first_order_density_matrix_size);
  _FWD_(int, rho_param.n_basis_local);
  _FWD_(int, rho_param.perm_n_full_points);
  _FWD_(int, max_n_batch_centers);
  _FW_(int, rho_param.basis_l_max, n_species);
  _FW_(int, rho_param.n_points_all_batches, n_my_batches_work_rho);
  _FW_(int, rho_param.n_batch_centers_all_batches, n_my_batches_work_rho);
  _FW_(int, rho_param.batch_center_all_batches, n_my_batches_work_rho * max_n_batch_centers);
  _FW_(int, rho_param.batch_point_to_i_full_point, n_my_batches_work_rho *n_max_batch_size);
  if(rho_param.n_basis_local > 0){
    // 注意这不是一行是多行 !!! 去掉大括号就有问题了
    _FW_(int, rho_param.ins_idx_all_batches, rho_param.n_basis_local *n_my_batches_work_rho);
    _FW_(double, rho_param.first_order_rho, rho_param.perm_n_full_points);
  } else {
    _FW_(double, rho_param.first_order_rho, n_full_points_work_rho);
  }
  _FW_(double, rho_param.first_order_density_matrix, rho_param.first_order_density_matrix_size);
  _FW_(double, rho_param.partition_tab, n_full_points_work_rho);
  // rho batch
  _FW_(int, MV(opencl_util, batches_size_rho), n_my_batches_work_rho);                             // 进程间有差别
  _FW_(int, MV(opencl_util, batches_batch_n_compute_rho), n_my_batches_work_rho);                  // 进程间有差别
  // _FW_(int, MV(opencl_util, batches_batch_i_basis_rho), n_centers_basis_I *n_my_batches_work_rho); // 进程间有差别
  _FW_(int, MV(opencl_util, batches_batch_i_basis_rho), n_max_compute_dens * n_my_batches_work_rho); // 进程间有差别
  _FW_(double, MV(opencl_util, batches_points_coords_rho),
       3 * n_max_batch_size * n_my_batches_work_rho); // 进程间有差别
  if (print)
    printf("rank%d, <save> 8 <opencl_util rho var> %ld bytes\n", myid, ftell(file_p));
  printf("rank%d, %ld bytes, mdim=%d\n", myid, ftell(file_p), n_cc_lm_ijk(l_pot_max) * n_atoms);
  fclose(file_p);
}

void m_save_load_H(const char *file_name, int is_load, int print) {
  size_t (*func)(void *restrict, size_t, size_t, FILE *restrict);
  FILE *file_p;
  if (is_load) {
    func = fread;
  } else {
    func = (size_t (*)(void *restrict, size_t, size_t, FILE *restrict))fwrite;
  }
  size_t status;
  char more_file_name[128];
  const char *tmp = file_name;
  const char *last_begin = NULL;
  while (*tmp != '\0') {
    if (*tmp == '/')
      last_begin = tmp;
    ++tmp;
  }
  if (last_begin != NULL) {
    last_begin++;
    memcpy(more_file_name, file_name, last_begin - file_name);
    sprintf(more_file_name + (last_begin - file_name), "H_%s", last_begin);
  } else
    sprintf(more_file_name, "H_%s", file_name);
  printf("%s\n", more_file_name);
  if (is_load) {
    file_p = fopen(more_file_name, "rb");
  } else {
    file_p = fopen(more_file_name, "wb");
  }
  _FWD_(int, H_param.j_coord);
  _FWD_(int, H_param.n_spin);
  _FWD_(int, H_param.l_ylm_max);
  _FWD_(int, H_param.n_basis_local);
  _FWD_(int, H_param.n_matrix_size);
  _FWD_(int, max_n_batch_centers);
  _FW_(int, H_param.basis_l_max, n_species);
  _FW_(int, H_param.n_points_all_batches, n_my_batches_work_h);
  _FW_(int, H_param.n_batch_centers_all_batches, n_my_batches_work_h);
  _FW_(int, H_param.batch_center_all_batches, n_my_batches_work_h * max_n_batch_centers);
  // _FW_(int, H_param.batch_point_to_i_full_point, n_my_batches_work_h *n_max_batch_size);
  if(H_param.n_basis_local > 0){
    _FW_(int, H_param.ins_idx_all_batches, H_param.n_basis_local *n_my_batches_work_h);
  }
  // _FW_(int, H_param.batches_batch_i_basis_h, n_centers_basis_I * n_my_batches_work_h);
  _FW_(double, H_param.partition_all_batches, n_max_batch_size * n_my_batches_work_h);
  _FW_(double, H_param.first_order_H, H_param.n_matrix_size * H_param.n_spin);
  _FW_(double, H_param.local_potential_parts_all_points, H_param.n_spin * n_full_points_work_h);
  _FW_(double, H_param.local_first_order_rho_all_batches, H_param.n_spin * n_max_batch_size * n_my_batches_work_h);
  _FW_(double, H_param.local_first_order_potential_all_batches, n_max_batch_size * n_my_batches_work_h);
  _FW_(double, H_param.local_dVxc_drho_all_batches, 3 * n_max_batch_size * n_my_batches_work_h);
  _FW_(double, H_param.local_rho_gradient, 3 * H_param.n_spin * n_max_batch_size);
  _FW_(double, H_param.first_order_gradient_rho, 3 * H_param.n_spin * n_max_batch_size);
  // H batch
  _FW_(int, MV(opencl_util, batches_size_h), n_my_batches_work_h);                             // 进程间有差别
  _FW_(int, MV(opencl_util, batches_batch_n_compute_h), n_my_batches_work_h);                  // 进程间有差别
  _FW_(int, MV(opencl_util, batches_batch_i_basis_h), n_max_compute_dens * n_my_batches_work_h); // 进程间有差别
  _FW_(double, MV(opencl_util, batches_points_coords_h),
       3 * n_max_batch_size * n_my_batches_work_h); // 进程间有差别
  if (print)
    printf("rank%d, <save> 8 <opencl_util h var> %ld bytes\n", myid, ftell(file_p));
  printf("rank%d, %ld bytes\n", myid, ftell(file_p));
  fclose(file_p);
}

void set_sum_up_param(int forces_on, double *partition_tab_std, double *delta_v_hartree, double *rho_multipole,
                      // double *centers_rho_multipole_spl, double *centers_delta_v_hart_part_spl,
                      double *adap_outer_radius_sq, double *multipole_radius_sq, int *l_hartree_max_far_distance,
                      double *outer_potential_radius, double *multipole_c) {
  sum_up_param.forces_on = forces_on;
  sum_up_param.partition_tab = partition_tab_std;
  sum_up_param.delta_v_hartree = delta_v_hartree;
  sum_up_param.rho_multipole = rho_multipole;
  // sum_up_param.centers_rho_multipole_spl = centers_rho_multipole_spl;
  // sum_up_param.centers_delta_v_hart_part_spl = centers_delta_v_hart_part_spl;
  sum_up_param.adap_outer_radius_sq = adap_outer_radius_sq;
  sum_up_param.multipole_radius_sq = multipole_radius_sq;
  sum_up_param.l_hartree_max_far_distance = l_hartree_max_far_distance;
  sum_up_param.outer_potential_radius = outer_potential_radius;
  sum_up_param.multipole_c = multipole_c;
}

void debug_check_output_04_(double* array, int* dims, int* dims_num_){
  int dims_num = *dims_num_;

  int* dims_count = (int*)malloc(sizeof(int) * dims_num);
  memset(dims_count, 0, sizeof(int) * dims_num);

  int offset = 0;
  int status_continue = 1;

  while(status_continue){
    for(int i=0; i<dims[dims_num-1]; i++){
      for(int j=0; j<dims_num-1; j++){
        printf("%d, ", dims_count[j]);
      }
      printf("%d, %.14f\n", i, array[offset]);
      offset++;
    }
    if(dims_num == 1)
      break;
    for(int i=dims_num-2; i>=0; i--){
      dims_count[i]++;
      if(dims_count[i] >= dims[i]){
        if(i == 0)
          status_continue = 0;
        dims_count[i] = 0;
      } else{
        break;
      }
    }
  }

  free(dims_count);
}

void debug_check_output_04_file_(double* array, int* dims, int* dims_num_, const char* file_name){
  int dims_num = *dims_num_;
  
  FILE *file_p = fopen(file_name, "w");

  int* dims_count = (int*)malloc(sizeof(int) * dims_num);
  memset(dims_count, 0, sizeof(int) * dims_num);

  int offset = 0;
  int status_continue = 1;

  while(status_continue){
    for(int i=0; i<dims[dims_num-1]; i++){
      for(int j=0; j<dims_num-1; j++){
        fprintf(file_p, "%d, ", dims_count[j]);
      }
      fprintf(file_p, "%d, %.14f\n", i, array[offset]);
      offset++;
    }
    if(dims_num == 1)
      break;
    for(int i=dims_num-2; i>=0; i--){
      dims_count[i]++;
      if(dims_count[i] >= dims[i]){
        if(i == 0)
          status_continue = 0;
        dims_count[i] = 0;
      } else{
        break;
      }
    }
  }

  free(dims_count);
  fclose(file_p);
}

void print_mem_info_(){
  struct sysinfo s_info;
  int error = sysinfo(&s_info);
  printf("rank%d, error0: %d, totalram: %.3f, freeram: %.3f\n", myid, error, s_info.totalram / (1024.0 * 1024 * 1024), s_info.freeram / (1024.0 * 1024 * 1024));
}

#undef _CHK_
#undef _FWD_
#undef _FW_
#undef _FWS_
#undef _NL_