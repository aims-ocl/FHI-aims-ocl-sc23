#ifndef __PASS_MOD_VAR_H__
#define __PASS_MOD_VAR_H__

// 需要注意，INTEL 编译器定义了 GNUC 这个宏，不能用来判断
// ifort 和 gfortran 都将 Fortran 中模块和变量名小写化了
// 因此后面的模块和变量名得改
#ifdef __INTEL_COMPILER
#define MV(mod_name, var_name) mod_name##_mp_##var_name##_
#else
#define MV(mod_name, var_name) __##mod_name##_MOD_##var_name
#endif

#define pmaxab 30
extern int n_coeff_hartree; // = 2;
extern double *Fp_function_spline_slice;
extern double *Fpc_function_spline_slice;
extern int sumup_c_count;
// ------

extern int MV(mpi_tasks, n_tasks);
extern int MV(mpi_tasks, myid);

extern int MV(dimensions, n_centers_hartree_potential);
extern int MV(dimensions, n_periodic);
extern int MV(dimensions, n_max_radial);
extern int MV(dimensions, l_pot_max);
extern int MV(dimensions, n_max_spline);
extern int MV(dimensions, n_hartree_grid);
extern int MV(dimensions, n_species);
extern int MV(dimensions, n_atoms);
extern int MV(dimensions, n_centers);
extern int MV(dimensions, n_centers_basis_integrals); // TODO
extern int MV(dimensions, n_centers_integrals);       // TODO
extern int MV(dimensions, n_max_compute_fns_ham);     // TODO
extern int MV(dimensions, n_basis_fns);               // TODO
extern int MV(dimensions, n_basis);                   // TODO
extern int MV(dimensions, n_centers_basis_t);         // TODO
extern int MV(dimensions, n_centers_basis_i);         // TODO
extern int MV(dimensions, n_max_grid);                // TODO
extern int MV(dimensions, n_max_compute_atoms);       // TODO
extern int MV(dimensions, n_max_compute_ham);         // TODO
extern int MV(dimensions, n_max_compute_dens);        // TODO
extern int MV(dimensions, n_max_batch_size);
extern int MV(dimensions, n_my_batches);  // 进程间不同
extern int MV(dimensions, n_full_points); // 进程间不同

extern int MV(runtime_choices, use_hartree_non_periodic_ewald);
extern int MV(runtime_choices, hartree_fp_function_splines);
extern int MV(runtime_choices, fast_ylm);
extern int MV(runtime_choices, new_ylm);
extern int MV(runtime_choices, flag_rel); // TODO
extern int MV(runtime_choices, adams_moulton_integrator);
extern int MV(runtime_choices, compensate_multipole_errors);

extern int *MV(geometry, species); // (n_atoms)
extern int *MV(geometry, empty);   // (n_atoms)

// extern int MV(pbc_lists, n_cells_in_hamiltonian);
extern int MV(pbc_lists, index_hamiltonian_dim2);
extern int MV(pbc_lists, position_in_hamiltonian_dim1);  // TODO
extern int MV(pbc_lists, position_in_hamiltonian_dim2);  // TODO
extern int MV(pbc_lists, column_index_hamiltonian_size); // TODO
extern int *MV(pbc_lists, centers_hartree_potential);    // (n_centers_hartree_potential)
extern int *MV(pbc_lists, center_to_atom);               // (n_centers)
extern int *MV(pbc_lists, species_center);               // (n_centers)
extern int *MV(pbc_lists, center_to_cell);               // (n_centers)  // TODO
extern int *MV(pbc_lists, cbasis_to_basis);              // (n_centers_basis_T) // TODO
extern int *MV(pbc_lists, cbasis_to_center);             // (n_centers_basis_T) // TODO
extern int *MV(pbc_lists, centers_basis_integrals);      // (n_centers_basis_integrals) // TODO
extern int *MV(pbc_lists, index_hamiltonian);            // (2, n_cells_in_hamiltonian, n_basis)  // TODO
extern int *MV(pbc_lists, position_in_hamiltonian);      // (.._dim1, .._dim2)  // TODO
extern int *MV(pbc_lists, column_index_hamiltonian);     // (..size)  // TODO
extern double *MV(pbc_lists, coords_center);             // (3,n_centers)  // 实际中涉及到重名的情况

extern int *MV(species_data, l_hartree); // (n_species)
extern double *MV(species_data, multipole_radius_free); // (n_species)

extern int *MV(grids, n_grid);            // (n_species)
extern int *MV(grids, n_radial);          // (n_species)
extern double *MV(grids, r_grid_min);     // (n_species)
extern double *MV(grids, r_grid_inc);     // (n_species)
extern double *MV(grids, log_r_grid_inc); // (n_species)
extern double *MV(grids, scale_radial);   // (n_species)
extern double *MV(grids, r_radial);       // (n_max_radial, n_species)
extern double *MV(grids, r_grid);         // (n_max_grid, n_species)

extern int MV(analytic_multipole_coefficients, l_max_analytic_multipole);
extern int *MV(analytic_multipole_coefficients, n_cc_lm_ijk);      // (0:l_max_analytic_multipole)
extern int *MV(analytic_multipole_coefficients, index_cc);         // (n_cc_lm_ijk(l_max_analytic_multipole),6)
extern int *MV(analytic_multipole_coefficients, index_ijk_max_cc); // (3,0:l_max_analytic_multipole)

extern int MV(hartree_potential_real_p0, n_hartree_atoms);
extern int MV(hartree_potential_real_p0, hartree_force_l_add);
// extern double *MV(hartree_potential_real_p0, multipole_c); // ( n_cc_lm_ijk(l_pot_max), n_atoms) // 每次 sumup 初始化
extern double MV(hartree_potential_real_p0, b0)[pmaxab + 1];     // (0:pmaxab)
extern double MV(hartree_potential_real_p0, b2)[pmaxab + 1];     // (0:pmaxab)
extern double MV(hartree_potential_real_p0, b4)[pmaxab + 1];     // (0:pmaxab)
extern double MV(hartree_potential_real_p0, b6)[pmaxab + 1];     // (0:pmaxab)
extern double MV(hartree_potential_real_p0, a_save)[pmaxab + 1]; // (0:pmaxab) // 改 fortran，新建一个变量放出来

extern int MV(hartree_f_p_functions, fp_max_grid);
extern int MV(hartree_f_p_functions, lmax_fp);
extern double MV(hartree_f_p_functions, fp_grid_min);
extern double MV(hartree_f_p_functions, fp_grid_inc);
extern double MV(hartree_f_p_functions, fp_grid_max);
extern double *MV(hartree_f_p_functions, fp_function_spline);    // (0:lmax_Fp,n_max_spline,Fp_max_grid)
extern double *MV(hartree_f_p_functions, fpc_function_spline);   // (0:lmax_Fp,n_max_spline,Fp_max_grid)
extern double MV(hartree_f_p_functions, ewald_radius_to)[11];    // 11
extern double MV(hartree_f_p_functions, inv_ewald_radius_to)[2]; // 2
extern double MV(hartree_f_p_functions, p_erfc_4)[6];            // 6
extern double MV(hartree_f_p_functions, p_erfc_5)[7];            // 7
extern double MV(hartree_f_p_functions, p_erfc_6)[8];            // 8

extern int MV(hartree_potential_storage, n_rho_multipole_atoms);
extern int *MV(hartree_potential_storage, rho_multipole_index); // (n_atoms)
extern int *MV(hartree_potential_storage, compensation_norm); // (n_atoms)
extern int *MV(hartree_potential_storage, compensation_radius); // (n_atoms)
extern double *MV(hartree_potential_storage,
                  rho_multipole); // ((l_pot_max+1)**2, n_max_radial+2, n_rho_multipole_atoms)
extern double *MV(hartree_potential_storage, rho_multipole_shmem_ptr);
extern int MV(hartree_potential_storage, use_rho_multipole_shmem);

extern int *MV(basis, perm_basis_fns_spl);       // (n_basis_fns)  // TOOD
extern double *MV(basis, outer_radius_sq);       // (n_basis_fns)  // TOOD
extern int *MV(basis, basis_fn);                 // (n_basis)  // TOOD
extern int *MV(basis, basis_l);                  // (n_basis)  // TOOD
extern double *MV(basis, atom_radius_sq);        // (n_species)  // TOOD
extern int *MV(basis, basis_fn_start_spl);       // (n_species)  // TOOD
extern int *MV(basis, basis_fn_atom);            // (n_basis_fns,n_atoms)  // TOOD
extern double *MV(basis, basis_wave_ordered);    // (n_basis_fns,n_max_spline, n_max_grid)  // TOOD
extern double *MV(basis, basis_kinetic_ordered); // (n_basis_fns,n_max_spline, n_max_grid)  // TOOD

// 千万注意，rho, sumup, H 的 batch 可能因为 batch_permutation 而出现不同
// 为保证兼容性，暂时每次都复制
// extern int MV(opencl_util, n_my_batches_work_stable);     // TODO
// extern int MV(opencl_util, n_full_points_work_stable);    // TODO
// extern int *MV(opencl_util, batches_size_s);                // (n_my_batches_work)  // 进程间不同
// extern int *MV(opencl_util, batches_batch_n_compute_s);     // (n_my_batches_work)  // 进程间不同
// extern int *MV(opencl_util, batches_batch_i_basis_s);       // (n_centers_basis_I, n_my_batches_work) // 进程间不同
// extern double *MV(opencl_util, batches_points_coords_s);    // (3, n_max_batch_size, n_my_batches_work) // 进程间不同
// // TODO extern double *MV(opencl_util, batches_points_coords_mc_s); // (3, n_max_batch_size, n_my_batches_work) //
// 进程间不同  // TODO

extern int MV(opencl_util, use_opencl_version);
extern int MV(opencl_util, use_sumup_pre_c_cl_version);

// sumup batch
extern int MV(opencl_util, n_my_batches_work_sumup);
extern int MV(opencl_util, n_full_points_work_sumup);
extern int *MV(opencl_util, batches_size_sumup); // (n_my_batches_work)  // 进程间不同
extern double *MV(opencl_util,
                  batches_points_coords_sumup); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO

// rho batch
extern int MV(opencl_util, n_my_batches_work_rho);
extern int MV(opencl_util, n_full_points_work_rho);
extern int *MV(opencl_util, batches_size_rho);            // (n_my_batches_work)  // 进程间不同
extern int *MV(opencl_util, batches_batch_n_compute_rho); // (n_my_batches_work)  // 进程间不同
extern int *MV(opencl_util, batches_batch_i_basis_rho);   // (n_centers_basis_I, n_my_batches_work) // 进程间不同
extern double *MV(opencl_util,
                  batches_points_coords_rho); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO

// H batch
extern int MV(opencl_util, n_my_batches_work_h);
extern int MV(opencl_util, n_full_points_work_h);
extern int *MV(opencl_util, batches_size_h);            // (n_my_batches_work)  // 进程间不同
extern int *MV(opencl_util, batches_batch_n_compute_h); // (n_my_batches_work)  // 进程间不同
extern int *MV(opencl_util, batches_batch_i_basis_h);   // (n_centers_basis_I, n_my_batches_work) // 进程间不同
extern double *MV(opencl_util,
                  batches_points_coords_h); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO

extern int MV(mpi_tasks, mpi_platform_relative_id);
extern int MV(mpi_tasks, mpi_per_node);
extern int MV(mpi_tasks, mpi_task_per_gpu);
extern int MV(opencl_util, max_n_batch_centers);
extern int MV(opencl_util, debug_io);

// ------

#define n_tasks MV(mpi_tasks, n_tasks)
#define myid MV(mpi_tasks, myid)

#define n_centers_hartree_potential MV(dimensions, n_centers_hartree_potential)
#define n_periodic MV(dimensions, n_periodic)
// #define n_my_batches_work MV(dimensions, n_my_batches)
#define n_max_radial MV(dimensions, n_max_radial)
#define l_pot_max MV(dimensions, l_pot_max)
#define n_max_spline MV(dimensions, n_max_spline)
#define n_hartree_grid MV(dimensions, n_hartree_grid)
#define n_species MV(dimensions, n_species)
#define n_atoms MV(dimensions, n_atoms)
#define n_centers MV(dimensions, n_centers)
#define n_centers_basis_integrals MV(dimensions, n_centers_basis_integrals)
#define n_centers_integrals MV(dimensions, n_centers_integrals)
#define n_max_compute_fns_ham MV(dimensions, n_max_compute_fns_ham)
#define n_basis_fns MV(dimensions, n_basis_fns)
#define n_basis MV(dimensions, n_basis)
#define n_centers_basis_T MV(dimensions, n_centers_basis_t)
#define n_centers_basis_I MV(dimensions, n_centers_basis_i)
#define n_max_grid MV(dimensions, n_max_grid)
#define n_max_compute_atoms MV(dimensions, n_max_compute_atoms)
#define n_max_compute_ham MV(dimensions, n_max_compute_ham)
#define n_max_compute_dens MV(dimensions, n_max_compute_dens)
#define n_max_batch_size MV(dimensions, n_max_batch_size)
// #define n_my_batches MV(dimensions, n_my_batches)
// #define n_full_points MV(dimensions, n_full_points)

#define use_hartree_non_periodic_ewald MV(runtime_choices, use_hartree_non_periodic_ewald)
#define hartree_fp_function_splines MV(runtime_choices, hartree_fp_function_splines)
#define fast_ylm MV(runtime_choices, fast_ylm)
#define new_ylm MV(runtime_choices, new_ylm)
#define flag_rel MV(runtime_choices, flag_rel)
#define Adams_Moulton_integrator MV(runtime_choices, adams_moulton_integrator)
#define compensate_multipole_errors MV(runtime_choices, compensate_multipole_errors)

#define species(i) MV(geometry, species)[(i)-1]
#define empty(i) MV(geometry, empty)[(i)-1]

// #define n_cells_in_hamiltonian MV(pbc_lists, n_cells_in_hamiltonian)
#define index_hamiltonian_dim2 MV(pbc_lists, index_hamiltonian_dim2)
#define position_in_hamiltonian_dim1 MV(pbc_lists, position_in_hamiltonian_dim1)
#define position_in_hamiltonian_dim2 MV(pbc_lists, position_in_hamiltonian_dim2)
#define column_index_hamiltonian_size MV(pbc_lists, column_index_hamiltonian_size)
#define centers_hartree_potential(i) MV(pbc_lists, centers_hartree_potential)[(i)-1]
#define center_to_atom(i) MV(pbc_lists, center_to_atom)[(i)-1]
#define species_center(i) MV(pbc_lists, species_center)[(i)-1]
#define center_to_cell(i) MV(pbc_lists, center_to_cell)[(i)-1]
#define centers_basis_integrals(i) MV(pbc_lists, centers_basis_integrals)[(i)-1]
#define Cbasis_to_basis(i) MV(pbc_lists, cbasis_to_basis)[(i)-1]
#define Cbasis_to_center(i) MV(pbc_lists, cbasis_to_center)[(i)-1]
#define coords_center(i, j) MV(pbc_lists, coords_center)[((j)-1) * 3 + (i)-1]
#define column_index_hamiltonian(i) MV(pbc_lists, column_index_hamiltonian)[(i)-1]
#define index_hamiltonian(i, j, k)                                                                                     \
  MV(pbc_lists, index_hamiltonian)[(((k)-1) * index_hamiltonian_dim2 + (j)-1) * 2 + (i)-1]
#define position_in_hamiltonian(i, j)                                                                                  \
  MV(pbc_lists, position_in_hamiltonian)[((i)-1) + ((j)-1) * position_in_hamiltonian_dim1]

#define l_hartree(i) MV(species_data, l_hartree)[(i)-1]
#define multipole_radius_free(i) MV(species_data, multipole_radius_free)[(i)-1]

#define n_grid(i) MV(grids, n_grid)[(i)-1]
#define n_radial(i) MV(grids, n_radial)[(i)-1]
#define r_grid_min(i) MV(grids, r_grid_min)[(i)-1]
#define r_grid_inc(i) MV(grids, r_grid_inc)[(i)-1]
#define log_r_grid_inc(i) MV(grids, log_r_grid_inc)[(i)-1]
#define scale_radial(i) MV(grids, scale_radial)[(i)-1]
#define r_radial(i, j) MV(grids, r_radial)[(i)-1+n_max_radial*((j)-1)]
#define r_grid(i, j) MV(grids, r_grid)[(i)-1+n_max_grid*((j)-1)]

#define l_max_analytic_multipole MV(analytic_multipole_coefficients, l_max_analytic_multipole)
#define n_cc_lm_ijk(i) MV(analytic_multipole_coefficients, n_cc_lm_ijk)[(i)]
#define index_cc(i, j, i_dim) MV(analytic_multipole_coefficients, index_cc)[((j)-1) * i_dim + (i)-1]
#define index_ijk_max_cc(i, j) MV(analytic_multipole_coefficients, index_ijk_max_cc)[(j)*3 + (i)-1]

#define n_hartree_atoms MV(hartree_potential_real_p0, n_hartree_atoms)
#define hartree_force_l_add MV(hartree_potential_real_p0, hartree_force_l_add)
// #define multipole_c(i, j) MV(hartree_potential_real_p0, multipole_c)[((j)-1) * n_cc_lm_ijk(l_pot_max) + (i)-1]
#define b0 MV(hartree_potential_real_p0, b0)         // info 没有减 1，因为从 0 开始
#define b2 MV(hartree_potential_real_p0, b2)         // info 没有减 1，因为从 0 开始
#define b4 MV(hartree_potential_real_p0, b4)         // info 没有减 1，因为从 0 开始
#define b6 MV(hartree_potential_real_p0, b6)         // info 没有减 1，因为从 0 开始
#define a_save MV(hartree_potential_real_p0, a_save) // info 没有减 1，因为从 0 开始

#define Fp_max_grid MV(hartree_f_p_functions, fp_max_grid)
#define lmax_Fp MV(hartree_f_p_functions, lmax_fp)
#define Fp_grid_min MV(hartree_f_p_functions, fp_grid_min)
#define Fp_grid_inc MV(hartree_f_p_functions, fp_grid_inc)
#define Fp_grid_max MV(hartree_f_p_functions, fp_grid_max)
#define Fp_function_spline MV(hartree_f_p_functions, fp_function_spline)
#define Fpc_function_spline MV(hartree_f_p_functions, fpc_function_spline)
#define Ewald_radius_to(i) MV(hartree_f_p_functions, ewald_radius_to)[(i)-1]
#define inv_Ewald_radius_to(i) MV(hartree_f_p_functions, inv_ewald_radius_to)[(i)-1]
#define P_erfc_4(i) MV(hartree_f_p_functions, p_erfc_4)[(i)-1]
#define P_erfc_5(i) MV(hartree_f_p_functions, p_erfc_5)[(i)-1]
#define P_erfc_6(i) MV(hartree_f_p_functions, p_erfc_6)[(i)-1]

#define n_rho_multipole_atoms MV(hartree_potential_storage, n_rho_multipole_atoms)
#define rho_multipole_index(i) MV(hartree_potential_storage, rho_multipole_index)[(i)-1]
#define compensation_norm(i) MV(hartree_potential_storage, compensation_norm)[(i)-1]
#define compensation_radius(i) MV(hartree_potential_storage, compensation_radius)[(i)-1]

#define perm_basis_fns_spl(i) MV(basis, perm_basis_fns_spl)[(i)-1]
#define outer_radius_sq(i) MV(basis, outer_radius_sq)[(i)-1]
#define basis_fn(i) MV(basis, basis_fn)[(i)-1]
#define basis_l(i) MV(basis, basis_l)[(i)-1]
#define atom_radius_sq(i) MV(basis, atom_radius_sq)[(i)-1]
#define basis_fn_start_spl(i) MV(basis, basis_fn_start_spl)[(i)-1]
#define basis_fn_atom(i, j) MV(basis, basis_fn_atom)[(i)-1 + ((j)-1) * n_basis_fns]
#define basis_wave_ordered MV(basis, basis_wave_ordered)
#define basis_kinetic_ordered MV(basis, basis_kinetic_ordered)

// #define n_my_batches_work_stable MV(opencl_util, n_my_batches_work_stable)
// #define n_full_points_work_stable MV(opencl_util, n_full_points_work_stable)
// #define batches_size_s(i) MV(opencl_util, batches_size_s)[(i)-1]
// #define batches_batch_n_compute_s(i) MV(opencl_util, batches_batch_n_compute_s)[(i)-1]
// #define batches_batch_i_basis_s(i, j) MV(opencl_util, batches_batch_i_basis_s)[(i)-1 + n_centers_basis_I * ((j)-1)]
// #define batches_points_coords_s(i, j, k) MV(opencl_util, batches_points_coords_s)[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]
// #define batches_points_coords_mc_s(i, j, k) MV(opencl_util, batches_points_coords_mc_s)[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]

#define ctrl_use_opencl_version MV(opencl_util, use_opencl_version)
#define ctrl_use_sumup_pre_c_cl_version MV(opencl_util, use_sumup_pre_c_cl_version)

#define max_n_batch_centers MV(opencl_util, max_n_batch_centers)

#define n_my_batches_work_h MV(opencl_util, n_my_batches_work_h)
#define n_full_points_work_h MV(opencl_util, n_full_points_work_h)
#define batches_size_h(i) MV(opencl_util, batches_size_h)[(i)-1]
#define batches_batch_n_compute_h(i) MV(opencl_util, batches_batch_n_compute_h)[(i)-1]
#define batches_batch_i_basis_h(i, j) MV(opencl_util, batches_batch_i_basis_h)[(i)-1 + n_max_compute_dens * ((j)-1)]
#define batches_points_coords_h(i, j, k)                                                                               \
  MV(opencl_util, batches_points_coords_h)[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]

#define n_my_batches_work_rho MV(opencl_util, n_my_batches_work_rho)
#define n_full_points_work_rho MV(opencl_util, n_full_points_work_rho)
#define batches_size_rho(i) MV(opencl_util, batches_size_rho)[(i)-1]
#define batches_batch_n_compute_rho(i) MV(opencl_util, batches_batch_n_compute_rho)[(i)-1]
#define batches_batch_i_basis_rho(i, j) MV(opencl_util, batches_batch_i_basis_rho)[(i)-1 + n_max_compute_dens * ((j)-1)]
#define batches_points_coords_rho(i, j, k)                                                                             \
  MV(opencl_util, batches_points_coords_rho)[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]

#define n_my_batches_work_sumup MV(opencl_util, n_my_batches_work_sumup)
#define n_full_points_work_sumup MV(opencl_util, n_full_points_work_sumup)
#define batches_size_sumup(i) MV(opencl_util, batches_size_sumup)[(i)-1]
#define batches_points_coords_sumup(i, j, k)                                                                           \
  MV(opencl_util, batches_points_coords_sumup)[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]

#define partition_tab(i) partition_tab_std[(i)-1]
#define multipole_radius_sq(i) multipole_radius_sq[(i)-1]
#define outer_potential_radius(i, j) outer_potential_radius[(((j)-1) * (l_pot_max + 1)) + (i)]

#define Fp(i, j) Fp[((j)-1) * (l_pot_max + 2) + (i)] // info 没有减 1，因为从 0 开始数 // TODO 验证大小

#define opencl_util_debug_io MV(opencl_util, debug_io)
// ------

typedef struct SUM_UP_PARAM_T {
  int forces_on;                         // int
  double *partition_tab;                 // (n_full_points_work)
  double *delta_v_hartree;               // (n_full_points_work)
  double *rho_multipole;                 // (n_full_points_work)
  double *centers_rho_multipole_spl;     // (l_pot_max+1)**2, n_max_spline, n_max_radial+2, n_atoms)
  double *centers_delta_v_hart_part_spl; // (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, n_atoms)
  double *adap_outer_radius_sq;          // (n_atoms)
  double *multipole_radius_sq;           // (n_atoms)
  int *l_hartree_max_far_distance;       // (n_atoms)
  double *outer_potential_radius;        // (0:l_pot_max, n_atoms)
  double *multipole_c;                   // (n_cc_lm_ijk(l_pot_max), n_atoms)
} SUM_UP_PARAM;

extern SUM_UP_PARAM sum_up_param;

typedef struct RHO_PARAM_T {
  int l_ylm_max;
  int n_local_matrix_size; // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
  int n_basis_local;       // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
  int perm_n_full_points;  // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
  int first_order_density_matrix_size;
  int *basis_l_max;
  int *n_points_all_batches;
  int *n_batch_centers_all_batches;
  int *batch_center_all_batches;
  int *batch_point_to_i_full_point;
  int *ins_idx_all_batches;
  double *first_order_rho;            // (n_full_points_work)
  double *first_order_density_matrix; // n_hamiltonian_matrix_size
  double *partition_tab;
} RHO_PARAM;

extern RHO_PARAM rho_param;

typedef struct H_PARAM_T {
  int j_coord;
  int n_spin;
  int l_ylm_max;
  int n_basis_local;
  int n_matrix_size;
  int *basis_l_max;                                // (n_species)
  int *n_points_all_batches;                       // (n_my_batches_work)
  int *n_batch_centers_all_batches;                // (n_my_batches_work)
  int *batch_center_all_batches;                   // (n_centers_integrals, n_my_batches_work)
  int *ins_idx_all_batches;                        // (n_basis_local, n_my_batches_work)
  int *batches_batch_i_basis_h;                    // (n_centers_basis_I, n_my_batches_work)
  double *partition_all_batches;                   // (n_max_batch_size, n_my_batches_work)
  double *first_order_H;                           // (n_matrix_size, n_spin)
  double *local_potential_parts_all_points;        // (n_spin, n_full_points)
  double *local_first_order_rho_all_batches;       // (n_spin, n_max_batch_size, n_my_batches_work)
  double *local_first_order_potential_all_batches; // (n_max_batch_size, n_my_batches_work)
  double *local_dVxc_drho_all_batches;             // (3, n_max_batch_size, n_my_batches_work)
  double *local_rho_gradient;                      // (3, n_spin, n_max_batch_size)
  double *first_order_gradient_rho;                // (3, n_spin, n_max_batch_size)
} H_PARAM;

extern H_PARAM H_param;

#endif