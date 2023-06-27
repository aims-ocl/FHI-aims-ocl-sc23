#ifdef fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

// #define WAVE_SIZE 64

// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #define kernel
// #define global
// #define local
// #define constant
#define L_POT_MAX 8
#define mconstant global const
#define LOCALSIZE_SUM_UP_PRE_PROC 1

#define MV(mod_name, var_name) var_name
#define pmaxab 30

#define n_centers_hartree_potential MV(dimensions, n_centers_hartree_potential)
#define n_periodic MV(dimensions, n_periodic)
#define n_my_batches_work MV(dimensions, n_my_batches)
#define n_max_radial MV(dimensions, n_max_radial)
#define l_pot_max MV(dimensions, l_pot_max)
#define n_max_spline MV(dimensions, n_max_spline)
#define n_hartree_grid MV(dimensions, n_hartree_grid)
#define n_species MV(dimensions, n_species)
#define n_atoms MV(dimensions, n_atoms)
#define n_centers MV(dimensions, n_centers)
#define n_max_batch_size MV(dimensions, n_max_batch_size)
#define n_my_batches MV(dimensions, n_my_batches)
#define n_full_points MV(dimensions, n_full_points)

#define use_hartree_non_periodic_ewald MV(runtime_choices, use_hartree_non_periodic_ewald)
#define hartree_fp_function_splines MV(runtime_choices, hartree_fp_function_splines)
#define fast_ylm MV(runtime_choices, fast_ylm)
#define new_ylm MV(runtime_choices, new_ylm)

#define species(i) MV(geometry, species)[(i)-1]
#define empty(i) MV(geometry, empty)[(i)-1]

#define centers_hartree_potential(i) MV(pbc_lists, centers_hartree_potential)[(i)-1]
#define center_to_atom(i) MV(pbc_lists, center_to_atom)[(i)-1]
#define species_center(i) MV(pbc_lists, species_center)[(i)-1]
#define coords_center(i, j) MV(pbc_lists, coords_center)[((j)-1) * 3 + (i)-1]

#define l_hartree(i) MV(species_data, l_hartree)[(i)-1]

#define n_grid(i) MV(grids, n_grid)[(i)-1]
#define n_radial(i) MV(grids, n_radial)[(i)-1]
#define batches_size_s(i) MV(grids, batches_size_s)[(i)-1]
#define batches_batch_n_compute_s(i) MV(grids, batches_batch_n_compute_s)[(i)-1]
#define batches_points_coords_s(i, j, k)                                                                               \
  MV(grids, batches_points_coords_s)[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]
#define r_grid_min(i) MV(grids, r_grid_min)[(i)-1]
#define log_r_grid_inc(i) MV(grids, log_r_grid_inc)[(i)-1]
#define scale_radial(i) MV(grids, scale_radial)[(i)-1]

#define l_max_analytic_multipole MV(analytic_multipole_coefficients, l_max_analytic_multipole)
#define n_cc_lm_ijk(i) MV(analytic_multipole_coefficients, n_cc_lm_ijk)[(i)]
#define index_cc(i, j, i_dim) MV(analytic_multipole_coefficients, index_cc)[((j)-1) * i_dim + (i)-1]
#define index_ijk_max_cc(i, j) MV(analytic_multipole_coefficients, index_ijk_max_cc)[(j)*3 + (i)-1]

#define n_hartree_atoms MV(hartree_potential_real_p0, n_hartree_atoms)
#define hartree_force_l_add MV(hartree_potential_real_p0, hartree_force_l_add)
#define multipole_c(i, j) MV(hartree_potential_real_p0, multipole_c)[((j)-1) * n_cc_lm_ijk(l_pot_max) + (i)-1]
// #define b0 MV(hartree_potential_real_p0, b0)         // info 没有减 1，因为从 0 开始
// #define b2 MV(hartree_potential_real_p0, b2)         // info 没有减 1，因为从 0 开始
// #define b4 MV(hartree_potential_real_p0, b4)         // info 没有减 1，因为从 0 开始
// #define b6 MV(hartree_potential_real_p0, b6)         // info 没有减 1，因为从 0 开始
// #define a_save MV(hartree_potential_real_p0, a_save) // info 没有减 1，因为从 0 开始

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

#define partition_tab(i) partition_tab_std[(i)-1]
#define multipole_radius_sq(i) multipole_radius_sq[(i)-1]
#define outer_potential_radius(i, j) outer_potential_radius[(((j)-1) * (l_pot_max + 1)) + (i)]

// #define Fp(i, j) Fp[(i) * lsize + lid] // info 没有减 1，因为从 0 开始数 // TODO 验证大小
#define Fp(i, j) Fp[(i)] // info 没有减 1，因为从 0 开始数 // TODO 验证大小

#define CEIL_DIV(a, b) (((a) + (b) - 1) / (b))

void SHEval_c_(int lmax, double sintheta, double costheta, double sinphi, double cosphi, double *pSH);

double invert_log_grid_c_(double r_current, double r_min, double scale) {
  return 1.0 + log(r_current / r_min) / log(scale);
}

// n_coeff -- number of spline coefficients, must be 2 or 4
void spline_vector_v2_c_(double r_output, global double *spl_param, int n_l_dim, int n_coeff, int n_grid_dim,
                         int n_points, int n_vector, double *out_result) {
#define spl_param(i, j, k)                                                                                             \
  spl_param[((((k)-1) * n_coeff) + ((j)-1)) * n_l_dim + (i)] // TODO 检查 spl_param 是否会导致切片问题
  int i_spl;
  double t, t2, t3, ta, tb, tc, td;
  i_spl = (int)r_output;
  i_spl = 1 > i_spl ? 1 : i_spl;
  i_spl = (n_points - 1) < i_spl ? (n_points - 1) : i_spl;
  t = r_output - (double)i_spl;
  if (n_coeff == 4) {
    t2 = t * t;
    t3 = t * t2;
    for (int i = 0; i < n_vector; i++)
      out_result[i] = spl_param(i, 1, i_spl) + spl_param(i, 2, i_spl) * t + spl_param(i, 3, i_spl) * t2 +
                      spl_param(i, 4, i_spl) * t3;
  } else {
    ta = (t - 1) * (t - 1) * (1 + 2 * t);
    tb = (t - 1) * (t - 1) * t;
    tc = t * t * (3 - 2 * t);
    td = t * t * (t - 1);
    for (int i = 0; i < n_vector; i++)
      out_result[i] = spl_param(i, 1, i_spl) * ta + spl_param(i, 2, i_spl) * tb + spl_param(i, 1, i_spl + 1) * tc +
                      spl_param(i, 2, i_spl + 1) * td;
  }
#undef spl_param
}

// n_coeff -- number of spline coefficients, must be 2 or 4
void spline_vector_v2_c_step(double r_output, global double *spl_param, int n_l_dim, int n_coeff, int n_grid_dim,
                         int n_points, int n_vector, global double *out_result, int step) {
#define spl_param(i, j, k)                                                                                             \
  spl_param[((((k)-1) * n_coeff) + ((j)-1)) * n_l_dim + (i)] // TODO 检查 spl_param 是否会导致切片问题
  int i_spl;
  double t, t2, t3, ta, tb, tc, td;
  i_spl = (int)r_output;
  i_spl = 1 > i_spl ? 1 : i_spl;
  i_spl = (n_points - 1) < i_spl ? (n_points - 1) : i_spl;
  t = r_output - (double)i_spl;
  if (n_coeff == 4) {
    t2 = t * t;
    t3 = t * t2;
    for (int i = 0; i < n_vector; i++)
      out_result[i * step] = spl_param(i, 1, i_spl) + spl_param(i, 2, i_spl) * t + spl_param(i, 3, i_spl) * t2 +
                      spl_param(i, 4, i_spl) * t3;
  } else {
    ta = (t - 1) * (t - 1) * (1 + 2 * t);
    tb = (t - 1) * (t - 1) * t;
    tc = t * t * (3 - 2 * t);
    td = t * t * (t - 1);
    for (int i = 0; i < n_vector; i++)
      out_result[i * step] = spl_param(i, 1, i_spl) * ta + spl_param(i, 2, i_spl) * tb + spl_param(i, 1, i_spl + 1) * tc +
                      spl_param(i, 2, i_spl + 1) * td;
  }
#undef spl_param
}

void spline_vector_v2_c_g_(double r_output, global double *spl_param, int n_l_dim, int n_coeff, int n_grid_dim,
                         int n_points, int n_vector, global double *out_result) {
#define spl_param(i, j, k)                                                                                             \
  spl_param[((((k)-1) * n_coeff) + ((j)-1)) * n_l_dim + (i)] // TODO 检查 spl_param 是否会导致切片问题
  int i_spl;
  double t, t2, t3, ta, tb, tc, td;
  i_spl = (int)r_output;
  i_spl = 1 > i_spl ? 1 : i_spl;
  i_spl = (n_points - 1) < i_spl ? (n_points - 1) : i_spl;
  t = r_output - (double)i_spl;
  if (n_coeff == 4) {
    t2 = t * t;
    t3 = t * t2;
    for (int i = 0; i < n_vector; i++)
      out_result[i] = spl_param(i, 1, i_spl) + spl_param(i, 2, i_spl) * t + spl_param(i, 3, i_spl) * t2 +
                      spl_param(i, 4, i_spl) * t3;
  } else {
    ta = (t - 1) * (t - 1) * (1 + 2 * t);
    tb = (t - 1) * (t - 1) * t;
    tc = t * t * (3 - 2 * t);
    td = t * t * (t - 1);
    for (int i = 0; i < n_vector; i++)
      out_result[i] = spl_param(i, 1, i_spl) * ta + spl_param(i, 2, i_spl) * tb + spl_param(i, 1, i_spl + 1) * tc +
                      spl_param(i, 2, i_spl + 1) * td;
  }
#undef spl_param
}

double spline_vector_v2_c_reduce(double r_output, global double *spl_param, int n_l_dim, int n_coeff, int n_grid_dim,
                         int n_points, int n_vector, double *ylm_tab) {
#define spl_param(i, j, k)                                                                                             \
  spl_param[((((k)-1) * n_coeff) + ((j)-1)) * n_l_dim + (i)] // TODO 检查 spl_param 是否会导致切片问题
  int i_spl;
  double t, t2, t3, ta, tb, tc, td;
  i_spl = (int)r_output;
  i_spl = 1 > i_spl ? 1 : i_spl;
  i_spl = (n_points - 1) < i_spl ? (n_points - 1) : i_spl;
  t = r_output - (double)i_spl;
  double out_result = 0.0;
  if (n_coeff == 4) {
    t2 = t * t;
    t3 = t * t2;
    // #pragma unroll 2
    for (int i = 0; i < n_vector; i++)
      out_result += (spl_param(i, 1, i_spl) + spl_param(i, 2, i_spl) * t + spl_param(i, 3, i_spl) * t2 +
                      spl_param(i, 4, i_spl) * t3) * ylm_tab[i];
  } else {
    ta = (t - 1) * (t - 1) * (1 + 2 * t);
    tb = (t - 1) * (t - 1) * t;
    tc = t * t * (3 - 2 * t);
    td = t * t * (t - 1);
    // #pragma unroll 2
    for (int i = 0; i < n_vector; i++)
      out_result += (spl_param(i, 1, i_spl) * ta + spl_param(i, 2, i_spl) * tb + spl_param(i, 1, i_spl + 1) * tc +
                      spl_param(i, 2, i_spl + 1) * td) * ylm_tab[i];
  }
  return out_result;
#undef spl_param
}

// INFO 小心 spl_param 的大小，c/c++ 可没有 fortran 的自动切片，应准确地为 spl_param(n_l_dim,4,n_grid_dim)
void spline_vector_c_(double r_output, global double *spl_param, int n_grid_dim, int n_l_dim, int n_points,
                      int n_vector, double *out_result) {
#define spl_param(i, j, k) spl_param[((((k)-1) * 4) + ((j)-1)) * n_l_dim + (i)]
  double t, term;
  int i_spl = (int)r_output;
  i_spl = 1 > i_spl ? 1 : i_spl;
  i_spl = (n_points - 1) < i_spl ? (n_points - 1) : i_spl;
  t = r_output - (double)i_spl;
  for (int i = 0; i < n_vector; i++)
    out_result[i] = spl_param(i, 1, i_spl);
  term = 1.0;
  for (int i_term = 2; i_term <= 4; i_term++) {
    term = term * t;
    for (int i = 0; i < n_vector; i++)
      out_result[i] += term * spl_param(i, i_term, i_spl);
  }
#undef spl_param
}

void F_erf_table_original_c_(double *F_erf_table, double r, int p_max) {
  printf("%s, not finished\n", __func__); // TODO
  // exit(-19);
}
void F_erfc_table_original_c_(double *F_table, double r, int p_max) {
  printf("%s, not finished\n", __func__); // TODO
  // exit(-19);
}

// F_erf + F_erfc
void F_erf_c_(double *F, double r, int p, int c,
              // outer
              int hartree_fp_function_splines, int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc,
              global double *Fp_function_spline_slice, global double *Fpc_function_spline_slice) {
  if (hartree_fp_function_splines) {
    // F_erf_spline + F_erfc_spline
    double rlog = invert_log_grid_c_(r, Fp_grid_min, Fp_grid_inc);
    global double *spl_param = (c != 0) ? Fpc_function_spline_slice : Fp_function_spline_slice;
    spline_vector_c_(rlog, spl_param, Fp_max_grid, lmax_Fp + 1, Fp_max_grid, p + 1, F);
  } else {
    if (c)
      F_erfc_table_original_c_(F, r, p);
    else
      F_erf_table_original_c_(F, r, p);
  }
}

// info: non_peri_extd 是可选参数，由于 C 没有可选参数，不用请置空
void far_distance_hartree_fp_periodic_single_atom_c_(
    int current_atom, int i_center, double dist,
    global int *l_hartree_max_far_distance, // (n_atoms)
    int inside, int forces_on, double multipole_radius_sq,
    double adap_outer_radius, // int *non_peri_extd
    // outer
    int l_pot_max, double *Fp, global double *b0, global double *b2, global double *b4, global double *b6,
    global double *a_save, int hartree_force_l_add, int use_hartree_non_periodic_ewald, int hartree_fp_function_splines,
    int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc, global double *Fp_function_spline_slice,
    global double *Fpc_function_spline_slice) {
  double drel;
  int lid = get_local_id(0);
  int lsize = get_local_size(0);

  int lmax = l_hartree_max_far_distance[current_atom - 1] + hartree_force_l_add;
  if (!use_hartree_non_periodic_ewald) {
    // if (inside) {
    //   F_erf_c_(&Fp(0, i_center), dist, lmax, 0, hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min,
    //            Fp_grid_inc, Fp_function_spline_slice, Fpc_function_spline_slice);
    //   for (int i = 0; i <= lmax; i++)
    //     Fp(i, i_center) = -Fp(i, i_center);

    // } else {
    //   F_erf_c_(&Fp(0, i_center), dist, lmax, 1, hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min,
    //            Fp_grid_inc, Fp_function_spline_slice, Fpc_function_spline_slice);
    // }
    F_erf_c_(&Fp(0, i_center), dist, lmax, !inside, hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min,
               Fp_grid_inc, Fp_function_spline_slice, Fpc_function_spline_slice);
    if (inside) {
      for (int i = 0; i <= lmax; i++)
        Fp(i, i_center) = -Fp(i, i_center);
    }
  } else { // WARNING 这部分暂时没测到
    for (int i = 0; i <= lmax; i++)
      Fp(i, i_center) = 0.0;
    if (dist < adap_outer_radius) {
      drel = dist / adap_outer_radius;
      double drel_2 = drel * drel;
      double drel_4 = drel_2 * drel_2;
      double drel_6 = drel_2 * drel_4;
      double adap_outer_radius_power = adap_outer_radius;
      for (int p = 0; p <= lmax; p++) {
        Fp(p, i_center) = 1 / adap_outer_radius_power * (b0[p] + b2[p] * drel_2 + b4[p] * drel_4 + b6[p] * drel_6);
        adap_outer_radius_power *= adap_outer_radius * adap_outer_radius;
      }
      // if (non_peri_extd != NULL)
      //   for (int i = 0; i <= lmax; i++)
      //     Fp(i, i_center) = -Fp(i, i_center);
    }
    // if ((dist * dist) >= multipole_radius_sq && !(dist < adap_outer_radius && non_peri_extd != NULL)) {
    if ((dist * dist) >= multipole_radius_sq) {
      double dist_power = dist;
      for (int p = 0; p <= lmax; p++) {
        Fp(p, i_center) += a_save[p] * dist_power;
        dist_power *= dist * dist;
      }
    }
  }
}

kernel void sum_up_whole_potential_shanghui_sub_t_(
    int forces_on, global double *partition_tab_std, global double *delta_v_hartree, global double *rho_multipole,
    global double *centers_rho_multipole_spl, global double *centers_delta_v_hart_part_spl,
    global double *adap_outer_radius_sq, global double *multipole_radius_sq, global int *l_hartree_max_far_distance,
    // global double *outer_potential_radius, global double *multipole_c,
    global double *outer_potential_radius, mconstant double *multipole_c,
    // outer
    // dimensions
    int n_centers_hartree_potential, int n_periodic, int n_max_radial, int l_pot_max, int n_max_spline,
    int n_hartree_grid, int n_species, int n_atoms, int n_centers, int n_max_batch_size, int n_my_batches,
    int n_full_points,
    // runtime_choices
    int use_hartree_non_periodic_ewald, int hartree_fp_function_splines, int fast_ylm, int new_ylm,
    // analytic_multipole_coefficients
    int l_max_analytic_multipole,
    // hartree_potential_real_p0
    int n_hartree_atoms, int hartree_force_l_add,
    // hartree_f_p_functions
    int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc, double Fp_grid_max,
    // outer arrays
    // geometry
    global int *species, // 从0开始数，第35个
    // pbc_lists
    global int *centers_hartree_potential, global int *center_to_atom, global int *species_center,
    global double *coords_center,
    // species_data
    global int *l_hartree,
    // grids
    global int *n_grid, global int *n_radial, global int *batches_size_s, global double *batches_points_coords_s,
    global double *r_grid_min, global double *log_r_grid_inc, global double *scale_radial,
    // analytic_multipole_coefficients
    global int *n_cc_lm_ijk, global int *index_cc, global int *index_ijk_max_cc,
    // hartree_potential_real_p0
    global double *b0, global double *b2, global double *b4, global double *b6, global double *a_save,
    // hartree_f_p_functions
    global double *Fp_function_spline_slice, global double *Fpc_function_spline_slice,
    // ------ loop helper ------
    int valid_max_point, global int *point_to_i_batch, global int *point_to_i_index,
    global int* valid_point_to_i_full_point,
    mconstant int *index_cc_aos,
    // global int *index_cc_aos,
    // ------ intermediate ------
    global double *Fp_all, global double *coord_c_all, global double *coord_mat_all, global double *rest_mat_all,
    global double *vector_all, global double *delta_v_hartree_multipole_component_all,
    global double *rho_multipole_component_all, global double *ylm_tab_all,
    int i_center_begin, int i_center_end, global int* i_center_to_centers_index

) {
#define centers_rho_multipole_spl(i, j, k, l)                                                                          \
  centers_rho_multipole_spl[((((l)-1) * (n_max_radial + 2) + (k)-1) * n_max_spline + (j)-1) * n_l_dim + (i)-1]
#define centers_delta_v_hart_part_spl(i, j, k, l)                                                                      \
  centers_delta_v_hart_part_spl[((((l)-1) * n_hartree_grid + (k)-1) * n_coeff_hartree + (j)-1) * n_l_dim + (i)-1]
  int gid = get_global_id(0);
  int gsize = get_global_size(0);
  int lid = get_local_id(0);
  int lsize = get_local_size(0);
  double Fp[(L_POT_MAX+2)];
  // global double *Fp = &Fp_all[get_group_id(0) * lsize * (l_pot_max + 2)];
  // global double *Fp = &Fp_all[gid * (l_pot_max + 2)];
  // global double *coord_c = coord_c_all + gid * 3 * (l_pot_max + 1);
  // global double *coord_mat = coord_mat_all + gid * (l_pot_max + 1) * (l_pot_max + 1);
  // global double *rest_mat = rest_mat_all + gid * (l_pot_max + 1) * (l_pot_max + 1);
  global double *vector = vector_all + gid * n_cc_lm_ijk(l_pot_max);
  // global double *delta_v_hartree_multipole_component =
  //     delta_v_hartree_multipole_component_all + gid * (l_pot_max + 1) * (l_pot_max + 1);
  // global double *rho_multipole_component = rho_multipole_component_all + gid * (l_pot_max + 1) * (l_pot_max + 1);
  // global double *ylm_tab = ylm_tab_all + gid * (l_pot_max + 1) * (l_pot_max + 1);

  // if(get_group_id(0) != 0)
  //   return;

  for (int i_center = i_center_begin+1; i_center <= i_center_end; i_center++) {
  // for (int i_center = 1; i_center <= n_centers_hartree_potential; i_center++) {
  // for (int i_center = 1; i_center <= 1; i_center++) {
    // int i_full_points = 0;
    // for (int i_batch = 1; i_batch <= n_my_batches; i_batch++)
    {
      // for (int i_index = 1; i_index <= batches_size_s(i_batch); i_index++)
      // for (int i_valid_points_ = gid; i_valid_points_ < valid_max_point; i_valid_points_ += gsize) {
      for (int i_valid_points_ = 0; i_valid_points_ < (valid_max_point+gsize-1) / gsize; i_valid_points_ += 1) {
        int i_valid_points = i_valid_points_ + gid * ((valid_max_point+gsize-1) / gsize);
        if(i_valid_points >= valid_max_point)
          break;
        // int i_full_points_ = i_valid_points_; // 没有 partition_tab 预判
        int i_full_points_ = valid_point_to_i_full_point[i_valid_points]; // 有 partition_tab 预判
        int i_batch = point_to_i_batch[i_full_points_];
        int i_index = point_to_i_index[i_full_points_];
        int i_full_points = i_full_points_ + 1;
        // i_full_points++;
        if (partition_tab(i_full_points) > 0.0) {
          int n_l_dim = (l_pot_max + 1) * (l_pot_max + 1);
          int n_coeff_hartree = 2;
          int current_center = centers_hartree_potential(i_center);
          int current_spl_atom = center_to_atom(current_center);
          int centers_xxx_index = i_center_to_centers_index[(i_center-1) - i_center_begin];
          double coord_current[3];
          for (int i = 0; i < 3; i++)
            coord_current[i] = batches_points_coords_s(i + 1, i_index, i_batch);
          double dist_tab_sq;
          double dir_tab[3];
          // tab_single_atom_centered_coords_p0_c_(&current_center, coord_current, &dist_tab_sq, dir_tab);
          {
            dist_tab_sq = 0.0;
            // printf("%ld\n", &coords_center(1, *current_center));
            for (int i_coord = 0; i_coord < 3; i_coord++) {
              dir_tab[i_coord] = coord_current[i_coord] - coords_center(i_coord + 1, current_center);
              dist_tab_sq += dir_tab[i_coord] * dir_tab[i_coord];
            }
          }
          int l_atom_max = l_hartree(species(current_spl_atom));
          while (outer_potential_radius(l_atom_max, current_spl_atom) < dist_tab_sq && l_atom_max > 0) {
            l_atom_max--;
          }
          int use_far_distance_hartree_fp_periodic_single_atom_c_ = 0;
          double dist_tab_in;
          double dist_tab_out;
          if (dist_tab_sq < multipole_radius_sq[current_spl_atom - 1]) {
            double i_r, i_r_log;
            double dir_tab_in[3];
            double trigonom_tab[4];
            // tab_single_atom_centered_coords_radial_log_p0_c_(&current_center, &dist_tab_sq, dir_tab, &dist_tab_in,
            // &i_r, &i_r_log, dir_tab_in);
            {
              dist_tab_in = sqrt(dist_tab_sq);
              dir_tab_in[0] = dir_tab[0] / dist_tab_in;
              dir_tab_in[1] = dir_tab[1] / dist_tab_in;
              dir_tab_in[2] = dir_tab[2] / dist_tab_in;
              // i_r_log = invert_log_grid_p2_c_(dist_tab_in, species_center(current_center));
              int i_species = species_center(current_center);
              i_r_log = 1.0 + log(dist_tab_in / r_grid_min(i_species)) / log_r_grid_inc(i_species);
              // i_r = invert_radial_grid_c_(dist_tab_in, n_radial(species_center(current_center)),
              //                             scale_radial(species_center(current_center)));
              i_r = (double)(n_radial(species_center(current_center)) + 1) *
                    sqrt(1.0 - exp(-dist_tab_in / scale_radial(species_center(current_center))));
            }
            // tab_single_trigonom_p0_c_(dir_tab_in, trigonom_tab);
            {
              double abmax, abcmax, ab, abc;
              abmax = fmax(fabs(dir_tab_in[0]), fabs(dir_tab_in[1]));
              if (abmax > 1.0e-36) {
                ab = sqrt(dir_tab_in[0] * dir_tab_in[0] + dir_tab_in[1] * dir_tab_in[1]);
                trigonom_tab[3] = dir_tab_in[0] / ab;
                trigonom_tab[2] = dir_tab_in[1] / ab;
              } else {
                trigonom_tab[3] = 1.0;
                trigonom_tab[2] = 0.0;
                ab = 0.0;
              }
              abcmax = fmax(abmax, fabs(dir_tab_in[2]));
              if (abcmax > 1.0e-36) {
                abc = sqrt(ab * ab + dir_tab_in[2] * dir_tab_in[2]);
                trigonom_tab[1] = dir_tab_in[2] / abc;
                trigonom_tab[0] = ab / abc;
              } else {
                trigonom_tab[1] = 0.0;
                trigonom_tab[0] = 1.0;
              }
            }
            double ylm_tab[(L_POT_MAX + 1) * (L_POT_MAX + 1)];
            // tab_single_wave_ylm_p2_c_(trigonom_tab, &l_atom_max, &l_pot_max, ylm_tab);
            {
              // if (fast_ylm) {
              //   SHEval_c_(l_atom_max, trigonom_tab[0], trigonom_tab[1], trigonom_tab[2], trigonom_tab[3], ylm_tab);
              // } else {
              //   printf("%s, not finished\n", __func__); // TODO
              //   // exit(-19);
              // }
              /* variables for tabulate ylm */
              double YLLI, YLL1I, YL1L1I, YLMI;
              double YLLR, YLL1R, YL1L1R, YLMR;
              int I2L, I4L2, INDEX, INDEX2, L, M, MSIGN;
              /* VB */
              int I22L, I24L2;
              double TEMP1, TEMP2, TEMP3;

              double D4LL1C, D2L13;
              const double PI = 3.14159265358979323846;

              #define trigonom_tab(i1) trigonom_tab[(i1) - 1]
              #define ylm_tab(i) ylm_tab[(i)-1]
                    if (0 <= 0)
                    {
                        YLLR = 1.0/sqrt(4.0*PI);
                        YLLI = 0.0;
                        ylm_tab(1) = YLLR;
                    }

                    if ( (0 <= 1) && (l_atom_max >= 1))
                    {
                        ylm_tab(3) = sqrt(3.00)*YLLR*trigonom_tab(2);
                        TEMP1 = -sqrt(3.00)*YLLR*trigonom_tab(1);
                        ylm_tab(4) = TEMP1*trigonom_tab(4);
                        ylm_tab(2) = -TEMP1*trigonom_tab(3);
                    }

                    // L = max(2,0)
                    for (L = 2; L <= l_atom_max; L++)
                    {
                        INDEX  = L*L+1;
                        INDEX2 = INDEX + 2*L;
                        MSIGN  = 1 - 2*(L%2);

                        YL1L1R = ylm_tab(INDEX-1);
                        YL1L1I = - MSIGN * ylm_tab(INDEX-2*L+1);
                        TEMP1 = -sqrt((double)(2*L+1)/(double)(2*L))*trigonom_tab(1);
                        YLLR = TEMP1*(trigonom_tab(4)*YL1L1R - trigonom_tab(3)*YL1L1I);
                        YLLI = TEMP1*(trigonom_tab(4)*YL1L1I + trigonom_tab(3)*YL1L1R);
                        ylm_tab(INDEX2) = YLLR;
                        ylm_tab(INDEX)  = MSIGN * YLLI;
                        INDEX2 = INDEX2 - 1;
                        INDEX  = INDEX  + 1;

                        TEMP2 = sqrt((double)(2*L+1))*trigonom_tab(2);
                        YLL1R = TEMP2*YL1L1R;
                        YLL1I = TEMP2*YL1L1I;
                        ylm_tab(INDEX2) = YLL1R;
                        ylm_tab(INDEX)  = - MSIGN * YLL1I;
                        INDEX2 = INDEX2 - 1;
                        INDEX  = INDEX  + 1;

                        I4L2 = INDEX - 4*L + 2;
                        I2L  = INDEX - 2*L;
                        I24L2 = INDEX2 - 4*L + 2;
                        I22L  = INDEX2 - 2*L;
                        D4LL1C = trigonom_tab(2)*sqrt((double)(4*L*L-1));
                        D2L13  = -sqrt((double)(2*L+1)/(double)(2*L-3));

                        for (M = L - 2; M >= 0; M--)
                        {
                            TEMP1 = 1.00/sqrt((double)((L+M)*(L-M)));
                            TEMP2 = D4LL1C*TEMP1;
                            TEMP3 = D2L13*sqrt((double)((L+M-1)*(L-M-1)))*TEMP1;
                            YLMR = TEMP2*ylm_tab(I22L) + TEMP3*ylm_tab(I24L2);
                            YLMI = TEMP2*ylm_tab(I2L) + TEMP3*ylm_tab(I4L2);
                            ylm_tab(INDEX2) = YLMR;
                            ylm_tab(INDEX)  = YLMI;

                            INDEX2 = INDEX2 - 1;
                            INDEX  = INDEX  + 1;
                            I24L2   = I24L2   - 1;
                            I22L    = I22L    - 1;
                            I4L2   = I4L2   + 1;
                            I2L    = I2L    + 1;
                        }
                    }
              #undef trigonom_tab
              #undef ylm_tab
            }
            int l_h_dim = (l_atom_max + 1) * (l_atom_max + 1);
            // // double delta_v_hartree_multipole_component[(l_pot_max + 1) * (l_pot_max + 1)];
            // // double rho_multipole_component[(l_pot_max + 1) * (l_pot_max + 1)];
            // double tmp_multipole_component[(L_POT_MAX + 1) * (L_POT_MAX + 1)];
            double delta_v_hartree_aux = spline_vector_v2_c_reduce(i_r_log, &centers_delta_v_hart_part_spl(1, 1, 1, centers_xxx_index+1), n_l_dim,
                                n_coeff_hartree, n_hartree_grid, n_grid(species_center(current_center)), l_h_dim,
                                ylm_tab);
            delta_v_hartree[i_full_points - 1] += delta_v_hartree_aux;
            double rho_multipole_aux = spline_vector_v2_c_reduce(i_r + 1, &centers_rho_multipole_spl(1, 1, 1, centers_xxx_index+1), n_l_dim, n_max_spline,
                                n_max_radial + 2, n_radial(species_center(current_center)) + 2, l_h_dim,
                                ylm_tab);
            rho_multipole[i_full_points - 1] += rho_multipole_aux;
            if (n_periodic > 0 || use_hartree_non_periodic_ewald) {
              use_far_distance_hartree_fp_periodic_single_atom_c_ = 1;
              //   // TODO WARNING far_distance_hartree_Fp_periodic_single_atom_c_ 里面有个
              //   // firstcall，要Fortran里预先调用一遍处理一下
              double tmp2 = sqrt(adap_outer_radius_sq[current_spl_atom - 1]);
              far_distance_hartree_fp_periodic_single_atom_c_(
                  current_spl_atom, i_center, dist_tab_in, l_hartree_max_far_distance, 1, forces_on,
                  multipole_radius_sq(current_spl_atom), tmp2,
                  // outer
                  l_pot_max, Fp, b0, b2, b4, b6, a_save, hartree_fp_function_splines, use_hartree_non_periodic_ewald,
                  hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min, Fp_grid_inc, Fp_function_spline_slice,
                  Fpc_function_spline_slice);
            }
          }
          else if (dist_tab_sq < adap_outer_radius_sq[current_spl_atom - 1]) {
            dist_tab_out = sqrt(dist_tab_sq);
            if (n_periodic == 0 && !use_hartree_non_periodic_ewald) {
              // if(!empty(*current_spl_atom)){}
              // far_distance_hartree_fp_cluster_single_atom_p2_c_(&dist_tab_out, &l_atom_max, forces_on);
              {
                double dist_tab = dist_tab_out;
                int l_max = l_atom_max;
                double dist_sq = dist_tab * dist_tab;
                int one_minus_2l = 1;
                Fp(0, 1) = 1.0 / dist_tab;
                for (int i_l = 1; i_l <= l_max + hartree_force_l_add; i_l++) {
                  one_minus_2l -= 2;
                  Fp(i_l, 1) = Fp(i_l - 1, 1) * (double)one_minus_2l / dist_sq;
                }
              }
              // far_distance_real_hartree_potential_single_atom_p2_c_(&i_center,
              // &delta_v_hartree[i_full_points - 1], &l_atom_max, coord_current);
              {
                int l_max = l_atom_max;
                double dpot = 0.0;
double coord_c[3][(L_POT_MAX + 1)];
#define coord_c(i, j) coord_c[i][j]
// #define coord_c(i, j) coord_c[(i) * (l_pot_max + 1) + (j)]
// #define coord_c(i, j) coord_c_all[((i) * (l_pot_max + 1) + (j)) * gsize + gid]
                double dir[3];
                coord_c(0, 0) = 1.0;
                coord_c(1, 0) = 1.0;
                coord_c(2, 0) = 1.0;
                dir[0] = coord_current[0] - coords_center(1, i_center);
                dir[1] = coord_current[1] - coords_center(2, i_center);
                dir[2] = coord_current[2] - coords_center(3, i_center);
                int maxval = -1;
                for (int i = 1; i <= 3; i++)
                  maxval = maxval > index_ijk_max_cc(i, l_max) ? maxval : index_ijk_max_cc(i, l_max);
                for (int i_l = 1; i_l <= maxval; i_l++) {
                  coord_c(0, i_l) = dir[0] * coord_c(0, i_l - 1);
                }
                for (int i_l = 1; i_l <= maxval; i_l++) {
                  coord_c(1, i_l) = dir[1] * coord_c(1, i_l - 1);
                }
                for (int i_l = 1; i_l <= maxval; i_l++) {
                  coord_c(2, i_l) = dir[2] * coord_c(2, i_l - 1);
                }
                int index_cc_i_dim = n_cc_lm_ijk(l_max_analytic_multipole);
                int multipole_c_size1 = n_cc_lm_ijk(l_pot_max);
                mconstant double * multipole_c_tmp = &multipole_c[1 - 1 + (center_to_atom(i_center) - 1) * multipole_c_size1];
                mconstant int* index_tmp;
                int ii, jj, kk, nn;

                int nmax = n_cc_lm_ijk(l_max) & 0xFFFFFFFC;
                for (int n = 1; n <= nmax; n+=4) {
                  index_tmp = index_cc_aos + (n-1)*4;
                  ii = index_tmp[0];  jj = index_tmp[1];  kk = index_tmp[2];  nn = index_tmp[3];
                  dpot = dpot + coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1) * (*(multipole_c_tmp++));
                  ii = index_tmp[4];  jj = index_tmp[5];  kk = index_tmp[6];  nn = index_tmp[7];
                  dpot = dpot + coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1) * (*(multipole_c_tmp++));
                  ii = index_tmp[8];  jj = index_tmp[9];  kk = index_tmp[10];  nn = index_tmp[11];
                  dpot = dpot + coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1) * (*(multipole_c_tmp++));
                  ii = index_tmp[12];  jj = index_tmp[13];  kk = index_tmp[14];  nn = index_tmp[15];
                  dpot = dpot + coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1) * (*(multipole_c_tmp++));
                }
                for (int n = nmax + 1; n <= n_cc_lm_ijk(l_max); n++) {
                  index_tmp = index_cc_aos + (n-1)*4;
                  ii = index_tmp[0];  jj = index_tmp[1];  kk = index_tmp[2];  nn = index_tmp[3];
                  dpot = dpot + coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1) * (*(multipole_c_tmp++));
                }
                
                // int nmax_ = n_cc_lm_ijk(l_max);
                // #pragma unroll 4
                // for (int n = 1; n <= nmax_; n++) {
                //   index_tmp = index_cc_aos + (n-1)*4;
                //   ii = index_tmp[0];  jj = index_tmp[1];  kk = index_tmp[2];  nn = index_tmp[3];
                //   dpot = dpot + coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1) * (*(multipole_c_tmp++));
                // }
                delta_v_hartree[i_full_points - 1] += dpot;
              }
            } else {
              use_far_distance_hartree_fp_periodic_single_atom_c_ = 2;
              double tmp2 = sqrt(adap_outer_radius_sq[current_spl_atom - 1]);
              far_distance_hartree_fp_periodic_single_atom_c_(
                  current_spl_atom, i_center, dist_tab_out, l_hartree_max_far_distance, 0, forces_on,
                  multipole_radius_sq(current_spl_atom), tmp2,
                  // outer
                  l_pot_max, Fp, b0, b2, b4, b6, a_save, hartree_fp_function_splines, use_hartree_non_periodic_ewald,
                  hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min, Fp_grid_inc, Fp_function_spline_slice,
                  Fpc_function_spline_slice);
            }
          }
// #ifdef N_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD1
          // if (n_periodic > 0 || use_hartree_non_periodic_ewald) {

          //   if(use_far_distance_hartree_fp_periodic_single_atom_c_ > 0){
          //     double arg2;
          //     int choose;
          //     if(use_far_distance_hartree_fp_periodic_single_atom_c_ == 1){
          //       arg2 = dist_tab_in;
          //       choose = 1;
          //     }
          //     else{
          //       arg2 = dist_tab_out;
          //       choose = 0;
          //     }
          //     double tmp2 = sqrt(adap_outer_radius_sq[current_spl_atom - 1]);
          //     far_distance_hartree_fp_periodic_single_atom_c_(
          //         current_spl_atom, i_center, arg2, l_hartree_max_far_distance, choose, forces_on,
          //         multipole_radius_sq(current_spl_atom), tmp2,
          //         // outer
          //         l_pot_max, Fp, b0, b2, b4, b6, a_save, hartree_fp_function_splines, use_hartree_non_periodic_ewald,
          //         hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min, Fp_grid_inc, Fp_function_spline_slice,
          //         Fpc_function_spline_slice);
          //   }
          // }

          if (n_periodic > 0 || use_hartree_non_periodic_ewald) {
            double tmp1 = adap_outer_radius_sq[current_spl_atom - 1];
            double tmp2 = multipole_radius_sq[current_spl_atom - 1];
            double tmp_max = tmp1 > tmp2 ? tmp1 : tmp2;
            if (dist_tab_sq < tmp_max) {
              // far_distance_real_hartree_potential_single_atom_c_(&current_center, &i_center,
              //    &delta_v_hartree[i_full_points - 1], l_hartree_max_far_distance,
              // coord_current);
              {
                double c_pot = 0.0;
                double dpot = 0.0;
                int l_max = l_hartree_max_far_distance[center_to_atom(current_center) - 1];
// double coord_c[3][l_pot_max + 1];
// double coord_mat[l_pot_max + 1][l_pot_max + 1];
// double rest_mat[l_pot_max + 1][l_pot_max + 1];
// double vector[n_cc_lm_ijk(l_pot_max)];
double coord_c[3][(L_POT_MAX + 1)];
#define coord_c(i, j) coord_c[i][j]
// #define coord_c(i, j) coord_c[(i) * (l_pot_max + 1) + (j)]
double coord_mat[(L_POT_MAX + 1)][(L_POT_MAX + 1)];
#define coord_mat(i, j) coord_mat[i][j]
// #define coord_mat(i, j) coord_mat[i * (l_pot_max + 1) + j]
double rest_mat[(L_POT_MAX + 1)][(L_POT_MAX + 1)];
#define rest_mat(i, j) rest_mat[i][j]
// #define rest_mat(i, j) rest_mat[i * (l_pot_max + 1) + j]
                double dir[3];
                coord_c(0, 0) = 1.0;
                coord_c(1, 0) = 1.0;
                coord_c(2, 0) = 1.0;
                dir[0] = coord_current[0] - coords_center(1, current_center);
                dir[1] = coord_current[1] - coords_center(2, current_center);
                dir[2] = coord_current[2] - coords_center(3, current_center);
                for (int i_coord = 0; i_coord < 3; i_coord++)
                  for (int i_l = 1; i_l <= index_ijk_max_cc(i_coord + 1, l_max); i_l++)
                    coord_c(i_coord, i_l) = dir[i_coord] * coord_c(i_coord, i_l - 1);
                int index_cc_i_dim = n_cc_lm_ijk(l_max_analytic_multipole);
                mconstant int* index_tmp;
                int multipole_c_size1 = n_cc_lm_ijk(l_pot_max);
                // mconstant double * multipole_c_tmp = &multipole_c[1 - 1 + (center_to_atom(current_center) - 1) * multipole_c_size1];
                int ii, jj, kk, nn;
                int nmax_ = n_cc_lm_ijk(l_max);
                #pragma unroll 4
                for (int n = 1; n <= nmax_; n++) {
                  index_tmp = index_cc_aos + (n-1)*4;
                  ii = index_tmp[0];  jj = index_tmp[1];  kk = index_tmp[2];  nn = index_tmp[3];
                  dpot += coord_c(0, ii) * coord_c(1, jj) * coord_c(2, kk) * Fp(nn, 1)
                                  * multipole_c[n - 1 + (center_to_atom(current_center) - 1) * multipole_c_size1];
                }
                if (fabs(dpot) > 1e-30)
                  c_pot += dpot;
                delta_v_hartree[i_full_points - 1] += c_pot;
              }
            }
          }
// #endif  // N_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD1
        }
      }
    }
  }

#undef centers_rho_multipole_spl
#undef centers_delta_v_hart_part_spl
}

kernel void sum_up_count_centers_hit(
    // int forces_on, 
    global double *partition_tab_std, 
    // global double *delta_v_hartree, global double *rho_multipole,
    // global double *centers_rho_multipole_spl, global double *centers_delta_v_hart_part_spl,
    // global double *adap_outer_radius_sq, 
    global double *multipole_radius_sq, 
    // global int *l_hartree_max_far_distance,
    // global double *outer_potential_radius, global double *multipole_c,
    // global double *outer_potential_radius, mconstant double *multipole_c,
    // outer
    // dimensions
    // int n_centers_hartree_potential, int n_periodic, int n_max_radial, 
    int l_pot_max, 
    // int n_max_spline,
    // int n_hartree_grid, int n_species, int n_atoms, int n_centers, 
    int n_max_batch_size, 
    // int n_my_batches,
    // int n_full_points,
    // runtime_choices
    // int use_hartree_non_periodic_ewald, int hartree_fp_function_splines, int fast_ylm, int new_ylm,
    // analytic_multipole_coefficients
    // int l_max_analytic_multipole,
    // hartree_potential_real_p0
    // int n_hartree_atoms, int hartree_force_l_add,
    // hartree_f_p_functions
    // int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc, double Fp_grid_max,
    // outer arrays
    // geometry
    // global int *species, // 从0开始数，第35个
    // pbc_lists
    global int *centers_hartree_potential, global int *center_to_atom, global int *species_center,
    global double *coords_center,
    // species_data
    // global int *l_hartree,
    // grids
    // global int *n_grid, global int *n_radial, global int *batches_size_s, 
    global double *batches_points_coords_s,
    // global double *r_grid_min, global double *log_r_grid_inc, global double *scale_radial,
    // analytic_multipole_coefficients
    // global int *n_cc_lm_ijk, global int *index_cc, global int *index_ijk_max_cc,
    // hartree_potential_real_p0
    // global double *b0, global double *b2, global double *b4, global double *b6, global double *a_save,
    // hartree_f_p_functions
    // global double *Fp_function_spline_slice, global double *Fpc_function_spline_slice,
    // ------ loop helper ------
    int valid_max_point, global int *point_to_i_batch, global int *point_to_i_index,
    global int* valid_point_to_i_full_point,
    // mconstant int *index_cc_aos,

    // global int *index_cc_aos,
    // ------ intermediate ------
    // global double *Fp_all, global double *coord_c_all, global double *coord_mat_all, global double *rest_mat_all,
    // global double *vector_all, global double *delta_v_hartree_multipole_component_all,
    // global double *rho_multipole_component_all, global double *ylm_tab_all,
    int i_center_begin, int i_center_end, global int* i_center_to_centers_index, global int* spl_atom_needs

) {
  int gid = get_global_id(0);
  int gsize = get_global_size(0);
  int lid = get_local_id(0);
  int lsize = get_local_size(0);


  for (int i_center = i_center_begin+1; i_center <= i_center_end; i_center++) {
    {
      for (int i_valid_points_ = 0; i_valid_points_ < (valid_max_point+gsize-1) / gsize; i_valid_points_ += 1) {
        int i_valid_points = i_valid_points_ + gid * ((valid_max_point+gsize-1) / gsize);
        if(i_valid_points >= valid_max_point)
          break;
        // int i_full_points_ = i_valid_points_; // 没有 partition_tab 预判
        int i_full_points_ = valid_point_to_i_full_point[i_valid_points]; // 有 partition_tab 预判
        int i_batch = point_to_i_batch[i_full_points_];
        int i_index = point_to_i_index[i_full_points_];
        int i_full_points = i_full_points_ + 1;
        // i_full_points++;
        if (partition_tab(i_full_points) > 0.0) {
          int n_l_dim = (l_pot_max + 1) * (l_pot_max + 1);
          int n_coeff_hartree = 2;
          int current_center = centers_hartree_potential(i_center);
          int current_spl_atom = center_to_atom(current_center);
          int centers_xxx_index = i_center_to_centers_index[(i_center-1) - i_center_begin];
          double coord_current[3];
          for (int i = 0; i < 3; i++)
            coord_current[i] = batches_points_coords_s(i + 1, i_index, i_batch);
          double dist_tab_sq;
          double dir_tab[3];
          {
            dist_tab_sq = 0.0;
            for (int i_coord = 0; i_coord < 3; i_coord++) {
              dir_tab[i_coord] = coord_current[i_coord] - coords_center(i_coord + 1, current_center);
              dist_tab_sq += dir_tab[i_coord] * dir_tab[i_coord];
            }
          }
          if (dist_tab_sq < multipole_radius_sq[current_spl_atom - 1]) {
            spl_atom_needs[current_spl_atom] = 1;
          }
        }
      }
    }
  }
}

double invert_radial_grid_c_(double r_current, int n_scale, double r_scale) {
  return (double)(n_scale + 1) * sqrt(1.0 - exp(-r_current / r_scale));
}

void cubic_spline_v2_c_(global double *spl_param, int *n_l_dim_, int *n_coeff_, int *n_grid_dim_, int n_points,
                        int *n_vector_) {
#define spl_param(i, j, k) spl_param[(i)-1 + n_l_dim * ((j)-1 + n_coeff * ((k)-1))]
#define d_inv(i) d_inv[(i)-1]
#define _MIN(i, j) ((i) < (j) ? (i) : (j))
  int n_l_dim = *n_l_dim_;
  int n_coeff = *n_coeff_;
  int n_grid_dim = *n_grid_dim_;
  // int n_points = *n_points_;
  int n_vector = *n_vector_;
  double d_inv[20];
  if (n_points == 1) {
    for (int i = 1; i <= n_vector; i++)
      spl_param(i, 2, 1) = 0;
    if (n_coeff == 4) {
      for (int i = 1; i <= n_vector; i++)
        spl_param(i, 3, 1) = 0;
      for (int i = 1; i <= n_vector; i++)
        spl_param(i, 4, 1) = 0;
    }
  }
  d_inv(1) = 0.5;
  for (int i = 2; i <= 20; i++)
    d_inv(i) = 1.0 / (4 - d_inv(i - 1));

  for (int i = 1; i <= n_vector; i++)
    spl_param(i, 2, 1) = 3 * (spl_param(i, 1, 2) - spl_param(i, 1, 1));

  for (int x = 2; x <= n_points - 1; x++) {
    int d_inv_id = _MIN(x - 1, 20);
    for (int i = 1; i <= n_vector; i++)
      spl_param(i, 2, x) =
          3 * (spl_param(i, 1, x + 1) - spl_param(i, 1, x - 1)) - d_inv(d_inv_id) * spl_param(i, 2, x - 1);
  }

  int d_inv_id = _MIN(n_points - 1, 20);
  for (int i = 1; i <= n_vector; i++)
    spl_param(i, 2, n_points) = 3 * (spl_param(i, 1, n_points) - spl_param(i, 1, n_points - 1)) -
                                d_inv(d_inv_id) * spl_param(i, 2, n_points - 1);
  for (int i = 1; i <= n_vector; i++)
    spl_param(i, 2, n_points) = spl_param(i, 2, n_points) / (2 - d_inv(d_inv_id)); // TODO 合并

  for (int x = n_points - 1; x >= 1; x--) {
    int d_inv_id = _MIN(x, 20);
    for (int i = 1; i <= n_vector; i++)
      spl_param(i, 2, x) = (spl_param(i, 2, x) - spl_param(i, 2, x + 1)) * d_inv(d_inv_id);
    if (n_coeff == 4) {
      for (int i = 1; i <= n_vector; i++)
        spl_param(i, 3, x) =
            3 * (spl_param(i, 1, x + 1) - spl_param(i, 1, x)) - 2 * spl_param(i, 2, x) - spl_param(i, 2, x + 1);
      for (int i = 1; i <= n_vector; i++)
        spl_param(i, 4, x) =
            2 * (spl_param(i, 1, x) - spl_param(i, 1, x + 1)) + spl_param(i, 2, x) + spl_param(i, 2, x + 1);
    }
  }

#undef _MIN
#undef d_inv
#undef spl_param
}

void get_rho_multipole_spl_c_(global double *rho_multipole_spl, int spl_atom, 
  // outer
  int l_pot_max, int n_max_radial, int n_max_spline,
  global int *species, global int *l_hartree, global int* n_radial, global double* scale_radial,
  global double* r_radial, global double* multipole_radius_free,
  global double* rho_multipole, global int* rho_multipole_index) {
#define r_radial(i, j) r_radial[(i)-1+n_max_radial*((j)-1)]
#define rho_multipole(i, j, k)                                                                                         \
  MV(hartree_potential_storage, rho_multipole)[(i)-1 + l_pot_max_help * ((j)-1 + n_max_radial_help * ((k)-1))]
#define rho_multipole_spl(i, j, k) rho_multipole_spl[(i)-1 + l_pot_max_help * ((j)-1 + n_max_spline * ((k)-1))]
  // int spl_atom = *spl_atom_;
  int l_pot_max_help = (l_pot_max + 1) * (l_pot_max + 1);
  int n_max_radial_help = n_max_radial + 2;

  int i_atom_index = rho_multipole_index[spl_atom-1];
  int species_tmp = species(spl_atom);
  int l_h_dim = (l_hartree(species_tmp) + 1) * (l_hartree(species_tmp) + 1);  // HIV 中此值可参考为 25
  int n_rad = n_radial(species_tmp);
  int n_rad_help = n_rad + 2;

  const int lid = get_local_id(0);
  const int lsize = get_local_size(0);

  for (int k = 1; k <= (n_rad + 2); k++) {
    for (int i = 1; i <= l_h_dim; i++) {
      rho_multipole_spl(i, 1, k) = rho_multipole(i, k, i_atom_index);
    }
  }

  // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

  // TODO first check this func
  cubic_spline_v2_c_(rho_multipole_spl, &l_pot_max_help, &n_max_spline, &n_max_radial_help, n_rad_help, &l_h_dim);

  // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

  // if(get_local_id(0) == 0){
    int i_radial = n_rad;
    while ((r_radial(i_radial, species_tmp) >= multipole_radius_free[species_tmp-1]) && i_radial > 1) {
      for (int j = 1; j <= n_max_spline; j++)
        for (int i = 1; i <= l_h_dim; i++)
          rho_multipole_spl(i, j, i_radial + 1) = 0;
      i_radial--;
    }

    double i_r_outer = invert_radial_grid_c_(multipole_radius_free[species_tmp-1], n_rad, scale_radial(species_tmp));
    double delta = (double)(i_r_outer - i_radial);
    double delta_2 = delta * delta;
    double delta_3 = delta_2 * delta;
    i_radial = i_radial + 1;
    for (int j = 1; j <= l_h_dim; j++) {
      rho_multipole_spl(j, 3, i_radial) =
          -3.0 / delta_2 * rho_multipole_spl(j, 1, i_radial) - 2.0 / delta * rho_multipole_spl(j, 2, i_radial);
      rho_multipole_spl(j, 4, i_radial) =
          2.0 / delta_3 * rho_multipole_spl(j, 1, i_radial) + 1.0 / delta_2 * rho_multipole_spl(j, 2, i_radial);
    }
  // }
#undef r_radial
#undef rho_multipole
#undef rho_multipole_spl
}

double compensating_density_c_(double radius, double r_outer, int l) {
  if (radius >= r_outer) {
    return 0;
  } else {
    double rl = 1.0;
    for (int i_l = 1; i_l <= l; i_l++) {
      rl *= radius;
    }
    double rfrac = radius / r_outer;
    double rfrac2 = rfrac * rfrac;
    double rfrac3 = rfrac2 * rfrac;
    return rl * (2.0 * rfrac3 - 3.0 * rfrac2 + 1.0);
  }
}

void spline_angular_integral_log_c_(global double *angular_part_spl, global double *angular_integral_log, int i_atom,
  // outer parameters
  int n_max_radial, int l_pot_max, int n_max_spline, int n_max_grid,
  int compensate_multipole_errors,
  // outer arrays
  global int* species, global int* l_hartree, global double* multipole_radius_free,
  global int* n_grid, global int* n_radial, global double* scale_radial, global double* r_grid,
  global int* compensation_norm, global int* compensation_radius
) {
#define r_grid(i, j) MV(grids, r_grid)[(i)-1+n_max_grid*((j)-1)]
#define multipole_radius_free(i) multipole_radius_free[(i)-1]
#define angular_part_spl(i, j, k) angular_part_spl[(i)-1 + l_pot_max_help * ((j)-1 + n_max_spline * ((k)-1))]
#define angular_integral_log(i, j) angular_integral_log[(i)-1 + l_pot_max_help * ((j)-1)]
// #define angular_integral_log(i, j) angular_integral_log[((i)-1) * n_max_grid + (j)-1]
  int l_pot_max_help = (l_pot_max + 1) * (l_pot_max + 1);
  int n_max_radial_help = n_max_radial + 2;
  // int i_atom = *i_atom_;

  const int lid = get_local_id(0);
  const int lsize = get_local_size(0);

  int l_h_dim = (l_hartree(species(i_atom)) + 1) * (l_hartree(species(i_atom)) + 1);

  int i_grid_max = n_grid(species(i_atom));
  for (int i_grid = 1; i_grid <= i_grid_max; i_grid++) {
    if (r_grid(i_grid, species(i_atom)) < multipole_radius_free(species(i_atom))) {
      double i_r_radial = 1 + invert_radial_grid_c_(r_grid(i_grid, species(i_atom)), n_radial(species(i_atom)),
                                                    scale_radial(species(i_atom)));
      spline_vector_v2_c_g_(i_r_radial, angular_part_spl, l_pot_max_help, n_max_spline, n_max_radial_help,
                          n_radial(species(i_atom)) + 2, l_h_dim, &angular_integral_log(1, i_grid));
    } else {
      for (int i = 1; i <= l_h_dim; i++)
        angular_integral_log(i, i_grid) = 0.0;
    }
  }
  // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
  if (compensate_multipole_errors) {
    double c_n_tmp = compensation_norm[i_atom-1];
    double c_r_tmp = compensation_radius[i_atom-1];
    int s_tmp = species(i_atom);
    int i_l = 0;
    for (int i_grid = 1; i_grid <= i_grid_max; i_grid++) {
      angular_integral_log(1, i_grid) += c_n_tmp * compensating_density_c_(r_grid(i_grid, s_tmp), c_r_tmp, i_l);
    }
  }

#undef r_grid
#undef multipole_radius_free
#undef angular_part_spl
#undef angular_integral_log
}

void integrate_delta_v_hartree_internal_c_(global double *angular_integral_log, 
  global double *delta_v_hartree, int *n_coeff_hartree_, int i_atom,
  // outer parameters
  int l_pot_max, int n_hartree_grid, int n_max_grid,
  int Adams_Moulton_integrator,
  // outer arrays
  global int* species, global int* l_hartree, global int* n_grid,
  global double* r_grid_inc, global double* r_grid,
  // local array
  local double* integral_zero_r_deprecate, local double* integral_r_infty_deprecate,
  local double (*local_grids_AM)[L_POT_MAX + 1],
  // local double (*local_iri)[(L_POT_MAX + 1) * (L_POT_MAX + 1)],
  local double (*local_iri)[LOCALSIZE_SUM_UP_PRE_PROC+1],
  local int* local_i_index_to_i_l
) {
#define r_grid(i, j) MV(grids, r_grid)[(i)-1+n_max_grid*((j)-1)]
#define angular_integral_log(i, j) angular_integral_log[(i)-1 + l_pot_max_help * ((j)-1)]
// #define angular_integral_log(i, j) angular_integral_log[((i)-1) * n_max_grid + (j)-1]
#define delta_v_hartree(i, j, k) delta_v_hartree[(i)-1 + l_pot_max_help * ((j)-1 + n_coeff_hartree * ((k)-1))]
#define integral_zero_r(i) integral_zero_r[(i)-1]
#define integral_r_infty(i) integral_r_infty[(i)-1]
#define d_1_24 (1.0 / 24.0)
  int l_pot_max_help = (l_pot_max + 1) * (l_pot_max + 1);
  int n_coeff_hartree = *n_coeff_hartree_;
  // int i_atom = *i_atom_;

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int lsize = get_local_size(0);

  double prefactor[L_POT_MAX + 1];  // generate once
  double integral_zero_r[(L_POT_MAX + 1) * (L_POT_MAX + 1)];
  double integral_r_infty[(L_POT_MAX + 1) * (L_POT_MAX + 1)];
  double r_l_AM[4];
  double r_neg_l1_AM[4];
  double r_inv_AM[4];
  double dr_coef_AM[4];

  // for(int i_index=lid; i_index<(l_pot_max + 1) * (l_pot_max + 1); i_index+=lsize){
  //   int i_l = (int)sqrt((float)i_index);  // the sqrt may be replaced by integer sqrt
  //   local_i_index_to_i_l[i_index] = i_l;
  // }

  for (int i_l = 0; i_l <= l_pot_max; i_l++)
    prefactor[i_l] = 12.56637061435917295376 / (2.0 * (double)i_l + 1.0); // pi4 = 12.56637061435917295376

  int l_h_dim = (l_hartree(species(i_atom)) + 1) * (l_hartree(species(i_atom)) + 1);
  double alpha = log(r_grid_inc[species(i_atom)-1]);

  for (int i = 0; i < (l_pot_max + 1) * (l_pot_max + 1); i++)
    integral_zero_r[i] = 0;
  for (int i = 0; i < (l_pot_max + 1) * (l_pot_max + 1); i++)
    integral_r_infty[i] = 0;
  // barrier(CLK_LOCAL_MEM_FENCE);

  // TODO print Adams_Moulton_integrator
  // WARNING the other path may mot be checked
  int s_tmp = species(i_atom);
  int n_tmp = n_grid(s_tmp);
  if (Adams_Moulton_integrator) {
    int i_l_max = l_hartree(species(i_atom));
    int i_index = 0; // TODO: change to static form in the loop
    // if(get_local_id(0) == 0){
      double r_grid1 = r_grid(1, species(i_atom));
      double r_grid2 = r_grid(2, species(i_atom));
      double r_grid3 = r_grid(3, species(i_atom));

      // len[-i_l, i_l] = 2 * i_l + 1
      // let a_i = 2i + 1
      // then \sum_{i=0}^{n} a_i = (1 + 2n+1) * (n+1) / 2 = (n+1)^2
      // so the loop1 can be convert to loop2
      // loop1
      // for (int i_l = 0; i_l <= i_l_max; i_l++) {
      //   for (int i_m = -i_l; i_m <= i_l; i_m++) {
      //     i_index++;
      // loop2
      // for(int i_index_=0; i_index_<(n+1)*(n+1); i_index_++){
      //   int i_l = (int)sqrt((float)i_index);
      //   int i_m = i_index - i_l * i_l - i_l;
      //   int i_index = i_index_ + 1;
      // }
      for (int i_l = 0; i_l <= i_l_max; i_l++) {
        for (int i_m = -i_l; i_m <= i_l; i_m++) {
          ++i_index;
          double tmp1 = angular_integral_log(i_index, 1) * alpha * pow(r_grid1, i_l + 3);
          double tmp2 = angular_integral_log(i_index, 2) * alpha * pow(r_grid2, i_l + 3);
          double tmp3 = angular_integral_log(i_index, 3) * alpha * pow(r_grid3, i_l + 3);
          integral_zero_r(i_index) += tmp1;
          delta_v_hartree(i_index, 1, 1) = integral_zero_r(i_index) / pow(r_grid1, i_l + 1);
          integral_zero_r(i_index) += (tmp1 + tmp2) * 0.5;
          delta_v_hartree(i_index, 1, 2) = integral_zero_r(i_index) / pow(r_grid2, i_l + 1);
          integral_zero_r(i_index) += (5 * tmp3 + 8 * tmp2 - tmp1) / 12.0;
          delta_v_hartree(i_index, 1, 2) = integral_zero_r(i_index) / pow(r_grid3, i_l + 1);
        }
      }
    // }
    
    // 0.016
    // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
    // if(get_local_id(0) == 0)
    // local double (*local_r_neg_l1_AM)[L_POT_MAX + 1] = local_grids_AM;
    for (int i_grid = 4; i_grid <= n_tmp; i_grid++) {
      for (int i_g = 0; i_g <= 3; i_g++) {
        r_inv_AM[i_g] = 1.0 / r_grid(i_grid - i_g, s_tmp);
        r_l_AM[i_g] = r_inv_AM[i_g];
        r_neg_l1_AM[i_g] = 1.0;
        dr_coef_AM[i_g] =
            r_grid(i_grid - i_g, s_tmp) * r_grid(i_grid - i_g, s_tmp) * alpha * r_grid(i_grid - i_g, s_tmp);
      }

      int i_index = 0;
      for (int i_l = 0; i_l <= i_l_max; i_l++) {
        for (int i_g = 0; i_g <= 3; i_g++)
          r_l_AM[i_g] *= r_grid(i_grid - i_g, s_tmp);
        r_neg_l1_AM[0] *= r_inv_AM[0];
        for (int i_m = 0; i_m <= i_l+i_l; i_m++) {
          ++i_index;
          integral_zero_r(i_index) += (9 * angular_integral_log(i_index, i_grid) * r_l_AM[0] * dr_coef_AM[0] +
                                       19 * angular_integral_log(i_index, i_grid - 1) * r_l_AM[1] * dr_coef_AM[1] -
                                       5 * angular_integral_log(i_index, i_grid - 2) * r_l_AM[2] * dr_coef_AM[2] +
                                       angular_integral_log(i_index, i_grid - 3) * r_l_AM[3] * dr_coef_AM[3]) *
                                      d_1_24;
          delta_v_hartree(i_index, 1, i_grid) = integral_zero_r(i_index) * r_neg_l1_AM[0];
        }
      }
    }
    // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
    // 0.034
    // if(get_local_id(0) == 0){
      for (int i = 0; i < (l_pot_max + 1) * (l_pot_max + 1); i++)
        integral_r_infty[i] = 0.0;
      // int i_index = 0;
      i_index = 0;
      double r_tmp = r_grid(n_tmp, s_tmp);
      double r_tmp1 = r_grid(n_tmp - 1, s_tmp);
      double r_tmp2 = r_grid(n_tmp - 2, s_tmp);
      for (int i_l = 0; i_l < i_l_max; i_l++) {
        for (int i_m = -i_l; i_m <= i_l; i_m++) {
          ++i_index;
          // TERM 1 : Integral_N = h*f_N; but h = 1
          // TODO: pow 改随迭代静态算
          delta_v_hartree(i_index, 1, n_tmp) += integral_r_infty(i_index) * pow(r_tmp, i_l);
          integral_r_infty(i_index) +=
              angular_integral_log(i_index, n_tmp) / pow(r_tmp, i_l + 1) * r_tmp * r_tmp * alpha * r_tmp;
          delta_v_hartree(i_index, 1, n_tmp) *= prefactor[i_l];
          // TERM 2 : Integral_(N-1) = Integral_N + h(f_(N-1)+f_N)/2
          delta_v_hartree(i_index, 1, n_tmp - 1) += integral_r_infty(i_index) * pow(r_tmp1, i_l);
          integral_r_infty(i_index) +=
              (angular_integral_log(i_index, n_tmp) / pow(r_tmp, (i_l + 1)) * r_tmp * r_tmp * alpha * r_tmp +
               angular_integral_log(i_index, n_tmp - 1) / pow(r_tmp1, (i_l + 1)) * r_tmp1 * r_tmp1 * alpha * r_tmp1) *
              0.5;
          delta_v_hartree(i_index, 1, n_tmp - 1) *= prefactor[i_l];
          // TERM 3 : Integral_(N-2) = Integral_(N-1) + h(5f_(N-2) + 8f_(N-1) - f_N)/12
          delta_v_hartree(i_index, 1, n_tmp - 2) += integral_r_infty(i_index) * pow(r_tmp2, i_l);
          integral_r_infty(i_index) +=
              (-1 * angular_integral_log(i_index, n_tmp) / pow(r_tmp, i_l + 1) * r_tmp * r_tmp * alpha * r_tmp +
               8 * angular_integral_log(i_index, n_tmp - 1) / pow(r_tmp1, i_l + 1) * r_tmp1 * r_tmp1 * alpha * r_tmp1 +
               5 * angular_integral_log(i_index, n_tmp - 2) / pow(r_tmp2, i_l + 1) * r_tmp2 * r_tmp2 * alpha * r_tmp2) /
              12.0;
          delta_v_hartree(i_index, 1, n_tmp - 2) *= prefactor[i_l];
        }
      }
    // }
    // all remaining terms
    // Integral_i = Integral_(i+1) + h[9 f_i + 19 f_(i+1) - 5 f_(i+2) + f_(i+3)]/24
    // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
    local double (*local_r_l_AM_0)[L_POT_MAX + 1] = local_grids_AM;
    for (int i_grid = n_tmp - 3; i_grid >= 1; i_grid--) {
      for (int i_g = 0; i_g <= 3; i_g++) {
        r_inv_AM[i_g] = 1.0 / r_grid(i_grid + i_g, s_tmp);
        r_l_AM[i_g] = r_inv_AM[i_g];
        r_neg_l1_AM[i_g] = 1.0;
        dr_coef_AM[i_g] =
            r_grid(i_grid + i_g, s_tmp) * r_grid(i_grid + i_g, s_tmp) * alpha * r_grid(i_grid + i_g, s_tmp);
      }
      int i_index = 0;
      for (int i_l = 0; i_l <= i_l_max; i_l++) {
        for (int i_g = 0; i_g <= 3; i_g++)
          r_neg_l1_AM[i_g] *= r_inv_AM[i_g];
        r_l_AM[0] *= r_grid(i_grid, s_tmp);
        for (int i_m = -i_l; i_m <= i_l; i_m++) {
          ++i_index;
          integral_r_infty(i_index) +=
              (9 * angular_integral_log(i_index, i_grid) * r_neg_l1_AM[0] * dr_coef_AM[0] +
               19 * angular_integral_log(i_index, i_grid + 1) * r_neg_l1_AM[1] * dr_coef_AM[1] -
               5 * angular_integral_log(i_index, i_grid + 2) * r_neg_l1_AM[2] * dr_coef_AM[2] +
               angular_integral_log(i_index, i_grid + 3) * r_neg_l1_AM[3] * dr_coef_AM[3]) *
              d_1_24;
          delta_v_hartree(i_index, 1, i_grid) += integral_r_infty(i_index) * r_l_AM[0];
          delta_v_hartree(i_index, 1, i_grid) *= prefactor[i_l];
        }
      }
    }
    // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
  } else { // TODO NO CHECKED (without cases which Adams_Moulton_integrator=False)
    // if(get_local_id(0) == 0){
    // Now to the integrations
    // First part of the integral 0 -> r
    for (int i_grid = 1; i_grid <= n_tmp; i_grid++) {
      double r_tmp = r_grid(i_grid, species(i_atom));
      double r_inv = 1.0 / r_tmp;
      double r_l = r_inv;
      double r_neg_l1 = 1.0;

      double dr_coef = r_tmp * r_tmp * alpha * r_tmp;
      int i_index = 0;
      int i_l_max = l_hartree(species(i_atom));
      for (int i_l = 0; i_l <= i_l_max; i_l++) {
        r_l *= r_tmp;
        r_neg_l1 *= r_inv;
        for (int i_m = -i_l; i_m <= i_l; i_m++) {
          ++i_index;
          integral_zero_r(i_index) += angular_integral_log(i_index, i_grid) * r_l * dr_coef;
          delta_v_hartree(i_index, 1, i_grid) = integral_zero_r(i_index) * r_neg_l1;
        }
      }
    }
    // run a second time through the radial grid from outward to inward
    // (not the whole integration grid!)
    // to evaluate integral_r_infty via tabulated angular_integral
    for (int i_grid = n_tmp; i_grid >= 1; i_grid--) {
      double r_tmp = r_grid(i_grid, species(i_atom));
      double r_inv = 1.0 / r_tmp;
      double r_l = r_inv;
      double r_neg_l1 = 1.0;
      double dr_coef = r_tmp * r_tmp * alpha * r_tmp;

      int i_index = 0;
      int i_l_max = l_hartree(species(i_atom));
      for (int i_l = 0; i_l <= i_l_max; i_l++) {
        r_l *= r_tmp;
        r_neg_l1 *= r_inv;
        for (int i_m = -i_l; i_m <= i_l; i_m++) {
          ++i_index;
          delta_v_hartree(i_index, 1, i_grid) += integral_r_infty(i_index) * r_l;
          integral_r_infty(i_index) += angular_integral_log(i_index, i_grid) * r_neg_l1 * dr_coef;
          delta_v_hartree(i_index, 1, i_grid) *= prefactor[i_l];
        }
      }
    }
    // }
  } // end if !(Adams_Moulton_integrator) and else
  // Calculate spline coefficients
  cubic_spline_v2_c_(delta_v_hartree, &l_pot_max_help, &n_coeff_hartree, &n_hartree_grid, n_grid(species(i_atom)),
                     &l_h_dim);
                    
#undef r_grid
#undef angular_integral_log
#undef delta_v_hartree
#undef integral_zero_r
#undef integral_r_infty
#undef d_1_24
}

// void integrate_delta_v_hartree_c_(double *angular_part_spl, double *delta_v_hartree, int *n_coeff_hartree_,
//                                   int *i_atom_) {
//   double angular_integral_log[(l_pot_max + 1) * (l_pot_max + 1) * n_max_grid];
//   spline_angular_integral_log_c_(angular_part_spl, angular_integral_log, i_atom_);
//   integrate_delta_v_hartree_internal_c_(angular_integral_log, delta_v_hartree, n_coeff_hartree_, i_atom_);
// }

kernel void sum_up_whole_potential_shanghui_pre_proc_(
  // outer parameters
  int n_max_radial, int l_pot_max, int n_max_spline, 
  int n_hartree_grid, int n_atoms, int n_max_grid,
  int Adams_Moulton_integrator, int compensate_multipole_errors,
  // outer arrays
  global int *species, global int *l_hartree, global double* multipole_radius_free,
  global int* n_grid, global int* n_radial, global double* r_grid_inc,
  global double* scale_radial, global double* r_grid, global double* r_radial,
  global double* rho_multipole, global int* rho_multipole_index,
  global int* compensation_norm, global int* compensation_radius,
  global int* centers_hartree_potential, global int* center_to_atom,
  // tmp
  global double *angular_integral_log_,
  // func spec
  global double *rho_multipole_spl, global double* delta_v_hartree, int n_coeff_hartree,
  int i_center_begin, int i_center_end, global int* i_center_to_centers_index,
  global int* spl_atom_needs
){
  const int gid = get_global_id(0);
  const int bid = get_group_id(0);
  const int bnum = get_num_groups(0);
  const int lid = get_local_id(0);
  const int lsize = get_local_size(0);
  // double angular_integral_log[(l_pot_max + 1) * (l_pot_max + 1) * n_max_grid];
  global double* angular_integral_log = &angular_integral_log_[bid * (l_pot_max + 1) * (l_pot_max + 1) * n_max_grid];

  local double integral_zero_r[(L_POT_MAX + 1) * (L_POT_MAX + 1)];
  local double integral_r_infty[(L_POT_MAX + 1) * (L_POT_MAX + 1)];
  local double local_iri[(L_POT_MAX + 1) * (L_POT_MAX + 1)][LOCALSIZE_SUM_UP_PRE_PROC+1];
  local double local_grids_AM[LOCALSIZE_SUM_UP_PRE_PROC][L_POT_MAX + 1];
  local int local_i_index_to_i_l[(L_POT_MAX + 1) * (L_POT_MAX + 1)];

  for(int i_center_=i_center_begin+bid; i_center_<i_center_end; i_center_+=bnum){
    // int i_center = i_center_ + 1;
    int current_center = MV(pbc_lists, centers_hartree_potential)[i_center_];
    int current_spl_atom = MV(pbc_lists, center_to_atom)[current_center - 1];
    if(!(spl_atom_needs[current_spl_atom] > 0))
      continue;
    int index = i_center_to_centers_index[i_center_];
    // if(gid <= 1)
    //   printf("%d: i_center_=%d, index=%d, i_center_begin=%d, current_center=%d, current_spl_atom=%d\n", gid, i_center_, index, 
    //               i_center_begin, current_center, current_spl_atom);
    if(index == (i_center_ - i_center_begin)){
      int rho_multipole_spl_offset = (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial+2) * index;
      int delta_v_hartree_offset = (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * index;
      // for(int i=lid; i<(l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial+2); i+=lsize)
      //   rho_multipole_spl[rho_multipole_spl_offset + i] = 0.0;
      // for(int i=lid; i<(l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid; i+=lsize)
      //   delta_v_hartree[delta_v_hartree_offset + i] = 0.0;
      
      get_rho_multipole_spl_c_(&rho_multipole_spl[rho_multipole_spl_offset], current_spl_atom, l_pot_max, n_max_radial, n_max_spline, 
        species, l_hartree, n_radial, scale_radial, r_radial, multipole_radius_free,
        rho_multipole, rho_multipole_index);

      // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

      spline_angular_integral_log_c_(&rho_multipole_spl[rho_multipole_spl_offset], angular_integral_log, current_spl_atom,
        n_max_radial, l_pot_max, n_max_spline, n_max_grid, compensate_multipole_errors,
        species, l_hartree, multipole_radius_free,
        n_grid, n_radial, scale_radial, r_grid,
        compensation_norm, compensation_radius);

      // barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

      integrate_delta_v_hartree_internal_c_(angular_integral_log, &delta_v_hartree[delta_v_hartree_offset], &n_coeff_hartree, current_spl_atom, 
        l_pot_max, n_hartree_grid, n_max_grid, Adams_Moulton_integrator,
        species, l_hartree, n_grid, r_grid_inc, r_grid, integral_zero_r, integral_r_infty,
        local_grids_AM, local_iri, local_i_index_to_i_l);
    }
  }
}

// void SHEval_c_(int lmax, double sintheta, double costheta, double sinphi, double cosphi, double *pSH) {
// // intent(inout) :: pSH((lmax+1)*(lmax+1))
// #define pSH(i) pSH[(i)-1]

//   double fX, fY, fZ;
//   double fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
//   // double fC0_1, fC1_1, fS0_1, fS1_1, fTmpA_1, fTmpB_1, fTmpC_1;
//   // double fTmpA_2, fTmpB_2, fTmpC_2;
//   double fZ2;
//   fX = sintheta * cosphi;
//   fY = sintheta * sinphi;
//   fZ = costheta;
//   fZ2 = fZ * fZ;

// #undef pSH
// }


#undef MV
#undef pmaxab

#undef n_centers_hartree_potential
#undef n_periodic
#undef n_my_batches_work
#undef n_max_radial
#undef l_pot_max
#undef n_max_spline
#undef n_hartree_grid
#undef n_species
#undef n_atoms
#undef n_centers
#undef n_max_batch_size
#undef n_my_batches
#undef n_full_points

#undef use_hartree_non_periodic_ewald
#undef hartree_fp_function_splines
#undef fast_ylm
#undef new_ylm

#undef species
#undef empty

#undef centers_hartree_potential
#undef center_to_atom
#undef species_center
#undef coords_center

#undef l_hartree

#undef n_grid
#undef n_radial
#undef batches_size_s
#undef batches_batch_n_compute_s
#undef batches_points_coords_s
#undef r_grid_min
#undef log_r_grid_inc
#undef scale_radial

#undef l_max_analytic_multipole
#undef n_cc_lm_ijk
#undef index_cc
#undef index_ijk_max_cc

#undef n_hartree_atoms
#undef hartree_force_l_add
#undef multipole_c
// #undef b0 MV(hartree_potential_real_p0, b0)         // info 没有减 1，因为从 0 开始
// #undef b2 MV(hartree_potential_real_p0, b2)         // info 没有减 1，因为从 0 开始
// #undef b4 MV(hartree_potential_real_p0, b4)         // info 没有减 1，因为从 0 开始
// #undef b6 MV(hartree_potential_real_p0, b6)         // info 没有减 1，因为从 0 开始
// #undef a_save MV(hartree_potential_real_p0, a_save) // info 没有减 1，因为从 0 开始

#undef Fp_max_grid
#undef lmax_Fp
#undef Fp_grid_min
#undef Fp_grid_inc
#undef Fp_grid_max
#undef Fp_function_spline
#undef Fpc_function_spline
#undef Ewald_radius_to
#undef inv_Ewald_radius_to
#undef P_erfc_4
#undef P_erfc_5
#undef P_erfc_6

#undef partition_tab
#undef multipole_radius_sq
#undef outer_potential_radius

#undef Fp