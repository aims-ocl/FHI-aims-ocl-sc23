#include "opencl_util_mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

OCL_UTIL_VARS ocl_util_vars_all[8];

void clear_ocl_util_vars(OCL_UTIL_VARS* ocl_util_vars){
    memset(ocl_util_vars, 0, sizeof(OCL_UTIL_VARS));
}

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
){
  ocl_util_vars->n_centers_hartree_potential = n_centers_hartree_potential;
  ocl_util_vars->n_periodic = n_periodic;
  ocl_util_vars->n_max_radial = n_max_radial;
  ocl_util_vars->l_pot_max = l_pot_max;
  ocl_util_vars->n_max_spline = n_max_spline;
  ocl_util_vars->n_hartree_grid = n_hartree_grid;
  ocl_util_vars->n_species = n_species;
  ocl_util_vars->n_atoms = n_atoms;
  ocl_util_vars->n_centers = n_centers;
  ocl_util_vars->n_centers_basis_integrals = n_centers_basis_integrals;
  ocl_util_vars->n_centers_integrals = n_centers_integrals;
  ocl_util_vars->n_max_compute_fns_ham = n_max_compute_fns_ham;
  ocl_util_vars->n_basis_fns = n_basis_fns;
  ocl_util_vars->n_basis = n_basis;
  ocl_util_vars->n_centers_basis_T = n_centers_basis_T;
  ocl_util_vars->n_centers_basis_I = n_centers_basis_I;
  ocl_util_vars->n_max_grid = n_max_grid;
  ocl_util_vars->n_max_compute_atoms = n_max_compute_atoms;
  ocl_util_vars->n_max_compute_ham = n_max_compute_ham;
  ocl_util_vars->n_max_compute_dens = n_max_compute_dens;
  ocl_util_vars->n_max_batch_size = n_max_batch_size;
  // runtime_choices
  ocl_util_vars->use_hartree_non_periodic_ewald = use_hartree_non_periodic_ewald;
  ocl_util_vars->hartree_fp_function_splines = hartree_fp_function_splines;
  // int fast_ylm,
  // int new_ylm,
  ocl_util_vars->flag_rel = flag_rel;
  ocl_util_vars->Adams_Moulton_integrator = Adams_Moulton_integrator;
  ocl_util_vars->compensate_multipole_errors = compensate_multipole_errors;
  // pbc_lists
  ocl_util_vars->index_hamiltonian_dim2 = index_hamiltonian_dim2;
  ocl_util_vars->position_in_hamiltonian_dim1 = position_in_hamiltonian_dim1;
  ocl_util_vars->position_in_hamiltonian_dim2 = position_in_hamiltonian_dim2;
  ocl_util_vars->column_index_hamiltonian_size = column_index_hamiltonian_size;
  // analytic_multipole_coefficients
  ocl_util_vars->l_max_analytic_multipole = l_max_analytic_multipole;
  // hartree_potential_real_p0
  ocl_util_vars->n_hartree_atoms = n_hartree_atoms;
  ocl_util_vars->hartree_force_l_add = hartree_force_l_add;
  // hartree_f_p_functions
  ocl_util_vars->Fp_max_grid = Fp_max_grid;
  ocl_util_vars->lmax_Fp = lmax_Fp;
  // hartree_potential_storage
  ocl_util_vars->n_rho_multipole_atoms = n_rho_multipole_atoms;
  // sumup batch
  ocl_util_vars->n_my_batches_work_sumup = n_my_batches_work_sumup;
  ocl_util_vars->n_full_points_work_sumup = n_full_points_work_sumup;
  // rho batch
  ocl_util_vars->n_my_batches_work_rho = n_my_batches_work_rho;
  ocl_util_vars->n_full_points_work_rho = n_full_points_work_rho;
  // h batch
  ocl_util_vars->n_my_batches_work_h = n_my_batches_work_h;
  ocl_util_vars->n_full_points_work_h = n_full_points_work_h;
  // hartree_f_p_functions double
  ocl_util_vars->Fp_grid_min = Fp_grid_min;
  ocl_util_vars->Fp_grid_inc = Fp_grid_inc;
  ocl_util_vars->Fp_grid_max = Fp_grid_max;
  // sumup
  ocl_util_vars->forces_on = forces_on;
  // help
  ocl_util_vars->n_cc_lm_ijk_l_pot_max = n_cc_lm_ijk_l_pot_max;
}

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
){
  ocl_util_vars->batches_size_sumup = batches_size_sumup;
  ocl_util_vars->batches_points_coords_sumup = batches_points_coords_sumup;
  ocl_util_vars->partition_tab = partition_tab;
  ocl_util_vars->delta_v_hartree = delta_v_hartree;
  ocl_util_vars->rho_multipole = rho_multipole;
  ocl_util_vars->adap_outer_radius_sq = adap_outer_radius_sq;
  ocl_util_vars->multipole_radius_sq = multipole_radius_sq;
  ocl_util_vars->l_hartree_max_far_distance = l_hartree_max_far_distance;
  ocl_util_vars->outer_potential_radius = outer_potential_radius;
  ocl_util_vars->multipole_c = multipole_c;

  // ocl_util_vars->point_to_i_batch = ocl_util_vars;
  // ocl_util_vars->point_to_i_index = ocl_util_vars;
  // ocl_util_vars->valid_point_to_i_full_point = ocl_util_vars;
}

// ---------- 宏污染区 -----------------

#include "pass_mod_var.h"

static int init_opencl_util_mpi_first_call = 1;

void init_opencl_util_mpi_(){
  if(init_opencl_util_mpi_first_call){
    for(int i=0; i<8; i++){
        clear_ocl_util_vars(&ocl_util_vars_all[i]);
    }
  }
  init_opencl_util_mpi_first_call = 0;

  set_ocl_util_vars_single(
    &ocl_util_vars_all[0],
    n_centers_hartree_potential,
    n_periodic,
    n_max_radial,
    l_pot_max,
    n_max_spline,
    n_hartree_grid,
    n_species,
    n_atoms,
    n_centers,
    n_centers_basis_integrals,
    n_centers_integrals,
    n_max_compute_fns_ham,
    n_basis_fns,
    n_basis,
    n_centers_basis_T,
    n_centers_basis_I,
    n_max_grid,
    n_max_compute_atoms,
    n_max_compute_ham,
    n_max_compute_dens,
    n_max_batch_size,
    // runtime_choices
    use_hartree_non_periodic_ewald,
    hartree_fp_function_splines,
    // fast_ylm,
    // new_ylm,
    flag_rel,
    Adams_Moulton_integrator,
    compensate_multipole_errors,
    // pbc_lists
    index_hamiltonian_dim2,
    position_in_hamiltonian_dim1,
    position_in_hamiltonian_dim2,
    column_index_hamiltonian_size,
    // analytic_multipole_coefficients
    l_max_analytic_multipole,
    // hartree_potential_real_p0
    n_hartree_atoms,
    hartree_force_l_add,
    // hartree_f_p_functions
    Fp_max_grid,
    lmax_Fp,
    // hartree_potential_storage
    n_rho_multipole_atoms,
    // sumup batch
    n_my_batches_work_sumup,
    n_full_points_work_sumup,
    // rho batch
    n_my_batches_work_rho,
    n_full_points_work_rho,
    // h batch
    n_my_batches_work_h,
    n_full_points_work_h,
    // hartree_f_p_functions double
    Fp_grid_min,
    Fp_grid_inc,
    Fp_grid_max,
    // sumup
    sum_up_param.forces_on,
    n_cc_lm_ijk(l_pot_max)
  );

  set_ocl_util_vars_arrays(
    &ocl_util_vars_all[0],
    MV(opencl_util, batches_size_sumup),
    MV(opencl_util, batches_points_coords_sumup),
    sum_up_param.partition_tab,
    sum_up_param.delta_v_hartree,
    sum_up_param.rho_multipole,
    sum_up_param.adap_outer_radius_sq,
    sum_up_param.multipole_radius_sq,
    sum_up_param.l_hartree_max_far_distance,
    sum_up_param.outer_potential_radius,
    sum_up_param.multipole_c
  );
}

#ifndef DEBUG_NO_MPI

void init_opencl_util_mpi3_(){
  if(init_opencl_util_mpi_first_call){
    for(int i=0; i<8; i++){
        clear_ocl_util_vars(&ocl_util_vars_all[i]);
    }
  }
  init_opencl_util_mpi_first_call = 0;

  set_ocl_util_vars_single(
    &ocl_util_vars_all[0],
    n_centers_hartree_potential,
    n_periodic,
    n_max_radial,
    l_pot_max,
    n_max_spline,
    n_hartree_grid,
    n_species,
    n_atoms,
    n_centers,
    n_centers_basis_integrals,
    n_centers_integrals,
    n_max_compute_fns_ham,
    n_basis_fns,
    n_basis,
    n_centers_basis_T,
    n_centers_basis_I,
    n_max_grid,
    n_max_compute_atoms,
    n_max_compute_ham,
    n_max_compute_dens,
    n_max_batch_size,
    // runtime_choices
    use_hartree_non_periodic_ewald,
    hartree_fp_function_splines,
    // fast_ylm,
    // new_ylm,
    flag_rel,
    Adams_Moulton_integrator,
    compensate_multipole_errors,
    // pbc_lists
    index_hamiltonian_dim2,
    position_in_hamiltonian_dim1,
    position_in_hamiltonian_dim2,
    column_index_hamiltonian_size,
    // analytic_multipole_coefficients
    l_max_analytic_multipole,
    // hartree_potential_real_p0
    n_hartree_atoms,
    hartree_force_l_add,
    // hartree_f_p_functions
    Fp_max_grid,
    lmax_Fp,
    // hartree_potential_storage
    n_rho_multipole_atoms,
    // sumup batch
    n_my_batches_work_sumup,
    n_full_points_work_sumup,
    // rho batch
    n_my_batches_work_rho,
    n_full_points_work_rho,
    // h batch
    n_my_batches_work_h,
    n_full_points_work_h,
    // hartree_f_p_functions double
    Fp_grid_min,
    Fp_grid_inc,
    Fp_grid_max,
    // sumup
    sum_up_param.forces_on,
    n_cc_lm_ijk(l_pot_max)
  );

  set_ocl_util_vars_arrays(
    &ocl_util_vars_all[0],
    MV(opencl_util, batches_size_sumup),
    MV(opencl_util, batches_points_coords_sumup),
    sum_up_param.partition_tab,
    sum_up_param.delta_v_hartree,
    sum_up_param.rho_multipole,
    sum_up_param.adap_outer_radius_sq,
    sum_up_param.multipole_radius_sq,
    sum_up_param.l_hartree_max_far_distance,
    sum_up_param.outer_potential_radius,
    sum_up_param.multipole_c
  );


    // 这一批还没处理呢
    // fp_function_spline
    // fpc_function_spline
    // ewald_radius_to
    // inv_ewald_radius_to
    // p_erfc_4
    // p_erfc_5
    // p_erfc_6


  if(MV(opencl_util, mpi_platform_relative_id) % MV(opencl_util, mpi_task_per_gpu) == 0){
    // recv
    for(int i=1; i<MV(opencl_util, mpi_task_per_gpu); i++){
      // printf("rank %d, call opencl_util_mpi_vars, prepare to recv from %d\n", myid, myid+i);
      opencl_util_mpi_vars(&ocl_util_vars_all[i], 0, myid+i);
      // printf("rank %d, call opencl_util_mpi_arrays_, prepare to recv from %d\n", myid, myid+i);
      opencl_util_mpi_arrays_(&ocl_util_vars_all[i], 0, myid+i);
      // printf("rank %d, call opencl_util_mpi_arrays_, recv from %d\n", myid, myid+i);
    }
  } else {
    // printf("rank %d, call opencl_util_mpi_vars, prepare to send to %d\n", myid, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
    // send
    opencl_util_mpi_vars(&ocl_util_vars_all[0], 1, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
    // printf("rank %d, call opencl_util_mpi_arrays_, prepare to send to %d\n", myid, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
    opencl_util_mpi_arrays_(&ocl_util_vars_all[0], 1, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
    // printf("rank %d, call opencl_util_mpi_arrays_, send to %d\n", myid, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
  }
  // printf("rank %d, finish init_opencl_util_mpi_\n", myid);
}

void finish_opencl_util_mpi_(){
  if(MV(opencl_util, mpi_platform_relative_id) % MV(opencl_util, mpi_task_per_gpu) == 0){
    // send
    for(int i=1; i<MV(opencl_util, mpi_task_per_gpu); i++){
      // printf("rank %d, call opencl_util_mpi_arrays_results_, prepare to send to %d\n", myid, myid+i);
      opencl_util_mpi_arrays_results_(&ocl_util_vars_all[i], 1, myid+i);
      // printf("rank %d, finish opencl_util_mpi_arrays_results_, send to %d\n", myid, myid+i);
    }
  } else {
    // printf("rank %d, call opencl_util_mpi_arrays_results_, prepare to recv from %d\n", myid, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
    // recv
    opencl_util_mpi_arrays_results_(&ocl_util_vars_all[0], 0, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
    // printf("rank %d, finish opencl_util_mpi_arrays_results_, recv from %d\n", myid, (myid/MV(opencl_util, mpi_task_per_gpu))*MV(opencl_util, mpi_task_per_gpu));
  }
  // printf("rank %d, finish init_opencl_util_mpi_\n", myid);
}

#endif