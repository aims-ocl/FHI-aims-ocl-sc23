#include "pass_mod_var.h"
#include "save_load_var.h"
#include "opencl_util_c.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int init_sum_up = 0;
int sumup_c_count = -1;
double *Fp_function_spline_slice = NULL;
double *Fpc_function_spline_slice = NULL;
double *Fp = NULL;
SUM_UP_PARAM sum_up_param;
int n_coeff_hartree = 2;

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
void spline_vector_v2_c_(double r_output, double *spl_param, int n_l_dim, int n_coeff, int n_grid_dim, int n_points,
                         int n_vector, double *out_result);
void spline_vector_c_(double r_output, double *spl_param, int n_grid_dim, int n_l_dim, int n_points, int n_vector,
                      double *out_result);
void F_erf_c_(double *F, double r, int p, int c,
              // outer
              int hartree_fp_function_splines, int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc,
              double *Fp_function_spline_slice, double *Fpc_function_spline_slice);
void F_erf_table_original_c_(double *F_erf_table, double r, int p_max);
void F_erfc_table_original_c_(double *F_table, double r, int p_max);
void far_distance_hartree_fp_periodic_single_atom_c_(
    int current_atom, int i_center, double dist,
    int *l_hartree_max_far_distance, // (n_atoms)
    int inside, int forces_on, double multipole_radius_sq,
    double adap_outer_radius, // int *non_peri_extd
    // outer
    int l_pot_max, double *Fp, double b0[pmaxab + 1], double b2[pmaxab + 1], double b4[pmaxab + 1],
    double b6[pmaxab + 1], double a_save[pmaxab + 1], int hartree_force_l_add, int use_hartree_non_periodic_ewald,
    int hartree_fp_function_splines, int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc,
    double *Fp_function_spline_slice, double *Fpc_function_spline_slice);
void far_distance_hartree_fp_cluster_single_atom_p2_c_(double *dist_tab, int *l_max, int *forces_on);
void far_distance_real_hartree_potential_single_atom_p2_c_(int *i_center, double *potential, int *l_max,
                                                           double coord_current[3]);
void far_distance_real_hartree_potential_single_atom_c_(int *current_center, int *i_center, double *potential,
                                                        int *l_hartree_max_far_distance, double coord_current[3]);

void init_sum_up_c_(int *forces_on, double *partition_tab_std, double *delta_v_hartree, double *rho_multipole,
                    double *adap_outer_radius_sq, double *multipole_radius_sq, int *l_hartree_max_far_distance,
                    double *outer_potential_radius, double *multipole_c) {
  // sumup_c_count++;
  // char save_file_name[64];
  set_sum_up_param(*forces_on, partition_tab_std, delta_v_hartree, rho_multipole,
                   adap_outer_radius_sq, multipole_radius_sq, l_hartree_max_far_distance,
                   outer_potential_radius, multipole_c);
  // sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid, sumup_c_count);
  // if (myid == 0 && sumup_c_count <= 1)
  // if (myid == 0)
  //   m_save_load_sumup(save_file_name, 0, 1);
  // if(sumup_c_count == 1 || myid == 0)
  //   m_save_load_sumup(save_file_name, 0, 0);

  m_save_load_sumup_();

  if (init_sum_up) {
    return;
  }
  init_sum_up = 1;
  // if (Fp_function_spline_slice == NULL)
  //   Fp_function_spline_slice = (double *)malloc(sizeof(double) * (Fp_max_grid+1) * 4 * (lmax_Fp + 1));
  // if (Fpc_function_spline_slice == NULL)
  //   Fpc_function_spline_slice = (double *)malloc(sizeof(double) * (Fp_max_grid+1) * 4 * (lmax_Fp + 1));
  // for (int i = 0; i < Fp_max_grid; i++) {
  //   memcpy(&Fp_function_spline_slice[i * 4 * (lmax_Fp + 1)], &Fp_function_spline[i * n_max_spline * (lmax_Fp + 1)],
  //          sizeof(double) * 4 * (lmax_Fp + 1));
  //   memcpy(&Fpc_function_spline_slice[i * 4 * (lmax_Fp + 1)], &Fpc_function_spline[i * n_max_spline * (lmax_Fp + 1)],
  //          sizeof(double) * 4 * (lmax_Fp + 1));
  // }
  if (Fp == NULL)
    Fp = (double *)malloc(sizeof(double) * (l_pot_max + 2) * n_centers_hartree_potential);
  // FILE *file_p = fopen("tmp.bin", "w");
  // for (int i = 0; i < Fp_max_grid; i++) {
  //   for (int j = 0; j < n_max_spline; j++) {
  //     for (int k = 0; k < (lmax_Fp + 1); k++) {
  //       fprintf(file_p, "%6d,%6d,%6d, %.13lf, %.13lf\n", i, j, k,
  //               Fp_function_spline[((i * n_max_spline) + j) * (lmax_Fp + 1) + k],
  //               Fpc_function_spline[((i * n_max_spline) + j) * (lmax_Fp + 1) + k]);
  //     }
  //   }
  // }
  // fclose(file_p);

  // opencl_init_();
  // opencl_common_buffer_init_();
}

void release_sum_up_c_(){
#define NOT_NULL_THEN_FREE(ptr) {if((ptr) != NULL) free(ptr);}
  NOT_NULL_THEN_FREE(Fp_function_spline_slice);
  NOT_NULL_THEN_FREE(Fpc_function_spline_slice);
  NOT_NULL_THEN_FREE(Fp);
#undef NOT_NULL_THEN_FREE
}

void sum_up_whole_potential_shanghui_sub_t_(int *forces_on, double *partition_tab_std, double *delta_v_hartree,
                                            double *rho_multipole, double *adap_outer_radius_sq,
                                            double *multipole_radius_sq, int *l_hartree_max_far_distance,
                                            double *outer_potential_radius, double *multipole_c) {
  if(ctrl_use_opencl_version){
    sum_up_begin();
    return;
  } else {
    printf("No implemented sum_up c version !!! ( %s:%d %s )\n", __FILE__, __LINE__, __func__);
    printf("Since the way to generate centers_rho_multipole_spl and centers_delta_v_hart_part_spl has been changed\n");
    fflush(stdout);
    exit(-1);
  }

  // TODO generate every i_center like fortran
  double* centers_rho_multipole_spl = (double*)malloc(sizeof(double) * (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2));
  double* centers_delta_v_hart_part_spl = (double*)malloc(sizeof(double) * (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid);

#define centers_rho_multipole_spl(i, j, k, l)                                                                          \
  centers_rho_multipole_spl[((((l)-1) * (n_max_radial + 2) + (k)-1) * n_max_spline + (j)-1) * n_l_dim + (i)-1]
#define centers_delta_v_hart_part_spl(i, j, k, l)                                                                      \
  centers_delta_v_hart_part_spl[((((l)-1) * n_hartree_grid + (k)-1) * n_coeff_hartree + (j)-1) * n_l_dim + (i)-1]
  int n_l_dim = (l_pot_max + 1) * (l_pot_max + 1);
  double coord_current[3];
  for (int i_center = 1; i_center <= n_centers_hartree_potential; i_center++) {
    int current_center = centers_hartree_potential(i_center);
    int current_spl_atom = center_to_atom(current_center);
    int i_full_points = 0;
    for (int i_batch = 1; i_batch <= n_my_batches_work_sumup; i_batch++) {
      for (int i_index = 1; i_index <= batches_size_sumup(i_batch); i_index++) {
        i_full_points++;
        if (partition_tab(i_full_points) > 0.0) {
          for (int i = 0; i < 3; i++)
            coord_current[i] = batches_points_coords_sumup(i + 1, i_index, i_batch);
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
          if (dist_tab_sq < multipole_radius_sq[current_spl_atom - 1]) {
            double dist_tab_in;
            double i_r, i_r_log;
            double dir_tab_in[3];
            double trigonom_tab[4];
            double ylm_tab[(l_pot_max + 1) * (l_pot_max + 1)];
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
            // tab_single_wave_ylm_p2_c_(trigonom_tab, &l_atom_max, &l_pot_max, ylm_tab);
            {
              // if (fast_ylm) {
              SHEval_c_(l_atom_max, trigonom_tab[0], trigonom_tab[1], trigonom_tab[2], trigonom_tab[3], ylm_tab);
              // } else {
              //   printf("%s, not finished\n", __func__); // TODO
              //   exit(-19);
              // }
            }
            int l_h_dim = (l_atom_max + 1) * (l_atom_max + 1);
            double delta_v_hartree_multipole_component[(l_pot_max + 1) * (l_pot_max + 1)];
            double rho_multipole_component[(l_pot_max + 1) * (l_pot_max + 1)];
            spline_vector_v2_c_(i_r_log, &centers_delta_v_hart_part_spl(1, 1, 1, current_spl_atom), n_l_dim,
                                n_coeff_hartree, n_hartree_grid, n_grid(species_center(current_center)), l_h_dim,
                                delta_v_hartree_multipole_component);
            double delta_v_hartree_aux = 0.0;
            for (int i = 0; i < l_h_dim; i++)
              delta_v_hartree_aux += delta_v_hartree_multipole_component[i] * ylm_tab[i];
            delta_v_hartree[i_full_points - 1] += delta_v_hartree_aux;
            spline_vector_v2_c_(i_r + 1, &centers_rho_multipole_spl(1, 1, 1, current_spl_atom), n_l_dim, n_max_spline,
                                n_max_radial + 2, n_radial(species_center(current_center)) + 2, l_h_dim,
                                rho_multipole_component);
            double rho_multipole_aux = 0.0;
            for (int i = 0; i < l_h_dim; i++) {
              rho_multipole_aux += rho_multipole_component[i] * ylm_tab[i];
            }
            rho_multipole[i_full_points - 1] += rho_multipole_aux;
            if (n_periodic > 0 || use_hartree_non_periodic_ewald) {
              // TODO WARNING far_distance_hartree_Fp_periodic_single_atom_c_ 里面有个
              // firstcall，要Fortran里预先调用一遍处理一下
              double tmp2 = sqrt(adap_outer_radius_sq[current_spl_atom - 1]);
              far_distance_hartree_fp_periodic_single_atom_c_(
                  current_spl_atom, i_center, dist_tab_in, l_hartree_max_far_distance, 1, *forces_on,
                  multipole_radius_sq(current_spl_atom), tmp2,
                  // outer
                  l_pot_max, Fp, b0, b2, b4, b6, a_save, hartree_fp_function_splines, use_hartree_non_periodic_ewald,
                  hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min, Fp_grid_inc, Fp_function_spline_slice,
                  Fpc_function_spline_slice);
            }
          } 
          else if (dist_tab_sq < adap_outer_radius_sq[current_spl_atom - 1]) {
            double dist_tab_out = sqrt(dist_tab_sq);
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
              // far_distance_real_hartree_potential_single_atom_p2_c_(&i_center, &delta_v_hartree[i_full_points - 1],
              //                                                       &l_atom_max, coord_current);
              {
                int l_max = l_atom_max;

                double dpot = 0.0;
                double coord_c[3][(l_pot_max + 1)];
                double dir[3];

                coord_c[0][0] = 1.0;
                coord_c[1][0] = 1.0;
                coord_c[2][0] = 1.0;

                dir[0] = coord_current[0] - coords_center(1, i_center);
                dir[1] = coord_current[1] - coords_center(2, i_center);
                dir[2] = coord_current[2] - coords_center(3, i_center);

                int maxval = -1;
                for (int i = 1; i <= 3; i++)
                  maxval = maxval > index_ijk_max_cc(i, l_max) ? maxval : index_ijk_max_cc(i, l_max);
                for (int i_l = 1; i_l <= maxval; i_l++) {
                  coord_c[0][i_l] = dir[0] * coord_c[0][i_l - 1];
                  coord_c[1][i_l] = dir[1] * coord_c[1][i_l - 1];
                  coord_c[2][i_l] = dir[2] * coord_c[2][i_l - 1];
                }
                int index_cc_i_dim = n_cc_lm_ijk(l_max_analytic_multipole);
                int multipole_c_size1 = n_cc_lm_ijk(l_pot_max);
                for (int n = 1; n <= n_cc_lm_ijk(l_max); n++) {
                  int ii = index_cc(n, 3, index_cc_i_dim);
                  int jj = index_cc(n, 4, index_cc_i_dim);
                  int kk = index_cc(n, 5, index_cc_i_dim);
                  int nn = index_cc(n, 6, index_cc_i_dim);
                  dpot = dpot + coord_c[0][ii] * coord_c[1][jj] * coord_c[2][kk] * Fp(nn, 1) *
                                    multipole_c[n - 1 + (center_to_atom(i_center) - 1) * multipole_c_size1];
                }
                delta_v_hartree[i_full_points - 1] += dpot;
              }
            } else {
              double tmp2 = sqrt(adap_outer_radius_sq[current_spl_atom - 1]);
              far_distance_hartree_fp_periodic_single_atom_c_(
                  current_spl_atom, i_center, dist_tab_out, l_hartree_max_far_distance, 0, *forces_on,
                  multipole_radius_sq(current_spl_atom), tmp2,
                  // outer
                  l_pot_max, Fp, b0, b2, b4, b6, a_save, hartree_fp_function_splines, use_hartree_non_periodic_ewald,
                  hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min, Fp_grid_inc, Fp_function_spline_slice,
                  Fpc_function_spline_slice);
            }
          }
          if (n_periodic > 0 || use_hartree_non_periodic_ewald) {
            double tmp1 = adap_outer_radius_sq[current_spl_atom - 1];
            double tmp2 = multipole_radius_sq[current_spl_atom - 1];
            double tmp_max = tmp1 > tmp2 ? tmp1 : tmp2;
            if (dist_tab_sq < tmp_max) {
              // far_distance_real_hartree_potential_single_atom_c_(&current_center, &i_center,
              //                                                    &delta_v_hartree[i_full_points - 1],
              //                                                    l_hartree_max_far_distance, coord_current);
              {
                double c_pot = 0.0;
                double dpot = 0.0;
                int l_max = l_hartree_max_far_distance[center_to_atom(current_center) - 1];
                double coord_c[3][l_pot_max + 1];
                double coord_mat[l_pot_max + 1][l_pot_max + 1];
                double rest_mat[l_pot_max + 1][l_pot_max + 1];
                double vector[n_cc_lm_ijk(l_pot_max)];
                double dir[3];
                coord_c[0][0] = 1.0;
                coord_c[1][0] = 1.0;
                coord_c[2][0] = 1.0;
                dir[0] = coord_current[0] - coords_center(1, current_center);
                dir[1] = coord_current[1] - coords_center(2, current_center);
                dir[2] = coord_current[2] - coords_center(3, current_center);
                for (int i_coord = 0; i_coord < 3; i_coord++)
                  for (int i_l = 1; i_l <= index_ijk_max_cc(i_coord + 1, l_max); i_l++)
                    coord_c[i_coord][i_l] = dir[i_coord] * coord_c[i_coord][i_l - 1];
                for (int i_l = 0; i_l <= index_ijk_max_cc(1, l_max); i_l++)
                  for (int i_l2 = 0; i_l2 <= index_ijk_max_cc(2, l_max); i_l2++)
                    coord_mat[i_l2][i_l] = coord_c[0][i_l] * coord_c[1][i_l2];
                for (int i_l = 0; i_l <= index_ijk_max_cc(3, l_max); i_l++)
                  for (int i_l2 = 0; i_l2 <= l_max; i_l2++)
                    rest_mat[i_l2][i_l] = coord_c[2][i_l] * Fp(i_l2, i_center);
                int index_cc_i_dim = n_cc_lm_ijk(l_max_analytic_multipole);
                for (int n = 1; n <= n_cc_lm_ijk(l_max); n++) {
                  int ii = index_cc(n, 3, index_cc_i_dim);
                  int jj = index_cc(n, 4, index_cc_i_dim);
                  int kk = index_cc(n, 5, index_cc_i_dim);
                  int nn = index_cc(n, 6, index_cc_i_dim);
                  vector[n - 1] = coord_mat[jj][ii] * rest_mat[nn][kk];
                }
                int multipole_c_size1 = n_cc_lm_ijk(l_pot_max);
                for (int n = 0; n < n_cc_lm_ijk(l_max); n++)
                  dpot +=
                      vector[n] * multipole_c[(n + 1) - 1 + (center_to_atom(current_center) - 1) * multipole_c_size1];
                if (fabs(dpot) > 1e-30)
                  c_pot += dpot;
                delta_v_hartree[i_full_points - 1] += c_pot;
              }
            }
          }
        }
      }
    }
  }
  if(myid == 0) m_save_check_sumup_(delta_v_hartree, rho_multipole);
#undef centers_rho_multipole_spl
#undef centers_delta_v_hart_part_spl
}

void sum_up_help_find_max_ir(double *partition_tab_std, double *multipole_radius_sq, double* i_r_log_max_, double* i_r_max_) {
  double i_r_log_max = 0;
  double i_r_max = 0;
  int n_l_dim = (l_pot_max + 1) * (l_pot_max + 1);
  double coord_current[3];
  for (int i_center = 1; i_center <= n_centers_hartree_potential; i_center++) {
    int current_center = centers_hartree_potential(i_center);
    int current_spl_atom = center_to_atom(current_center);
    int i_full_points = 0;
    for (int i_batch = 1; i_batch <= n_my_batches_work_sumup; i_batch++) {
      for (int i_index = 1; i_index <= batches_size_sumup(i_batch); i_index++) {
        i_full_points++;
        if (partition_tab(i_full_points) > 0.0) {
          for (int i = 0; i < 3; i++)
            coord_current[i] = batches_points_coords_sumup(i + 1, i_index, i_batch);
          double dist_tab_sq;
          double dir_tab[3];
          {
            dist_tab_sq = 0.0;
            // printf("%ld\n", &coords_center(1, *current_center));
            for (int i_coord = 0; i_coord < 3; i_coord++) {
              dir_tab[i_coord] = coord_current[i_coord] - coords_center(i_coord + 1, current_center);
              dist_tab_sq += dir_tab[i_coord] * dir_tab[i_coord];
            }
          }
          if (dist_tab_sq < multipole_radius_sq[current_spl_atom - 1]) {
            double dist_tab_in;
            double i_r, i_r_log;
            double dir_tab_in[3];
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
              // printf("%f, %f\n", i_r_log, i_r);
              if(i_r_log > i_r_log_max)
                i_r_log_max = i_r_log;
              if(i_r > i_r_max)
                i_r_max = i_r;
            }
          }
        }
      }
    }
  }
  *i_r_log_max_ = i_r_log_max;
  *i_r_max_ = i_r_max;
}

// in:  current_center, coord_current[3]
// out: dist_tab_sq, dir_tab[3]
void tab_single_atom_centered_coords_p0_c_(int *current_center, double coord_current[3], double *dist_tab_sq,
                                           double dir_tab[3]) {
  *dist_tab_sq = 0.0;
  // printf("%ld\n", &coords_center(1, *current_center));
  for (int i_coord = 0; i_coord < 3; i_coord++) {
    dir_tab[i_coord] = coord_current[i_coord] - coords_center(i_coord + 1, *current_center);
    *dist_tab_sq += dir_tab[i_coord] * dir_tab[i_coord];
  }
  return;
}

void tab_single_atom_centered_coords_radial_log_p0_c_(int *current_center, double *dist_tab_sq, double dir_tab[3],
                                                      double *dist_tab_in, double *i_r, double *i_r_log,
                                                      double dir_tab_in[3]) {
  *dist_tab_in = sqrt(*dist_tab_sq);
  dir_tab_in[0] = dir_tab[0] / *dist_tab_in;
  dir_tab_in[1] = dir_tab[1] / *dist_tab_in;
  dir_tab_in[2] = dir_tab[2] / *dist_tab_in;
  *i_r_log = invert_log_grid_p2_c_(*dist_tab_in, species_center(*current_center));
  *i_r = invert_radial_grid_c_(*dist_tab_in, n_radial(species_center(*current_center)),
                               scale_radial(species_center(*current_center)));
}

double invert_log_grid_c_(double r_current, double r_min, double scale) {
  return 1.0 + log(r_current / r_min) / log(scale);
}
double invert_log_grid_p2_c_(double r_current, int i_species) {
  return 1.0 + log(r_current / r_grid_min(i_species)) / log_r_grid_inc(i_species);
}
double invert_radial_grid_c_(double r_current, int n_scale, double r_scale) {
  return (double)(n_scale + 1) * sqrt(1.0 - exp(-r_current / r_scale));
}

void tab_single_trigonom_p0_c_(double dir_tab[3], double trigonom_tab[4]) {
  double abmax, abcmax, ab, abc;
  abmax = fmax(fabs(dir_tab[0]), fabs(dir_tab[1]));
  if (abmax > 1.0e-36) {
    ab = sqrt(dir_tab[0] * dir_tab[0] + dir_tab[1] * dir_tab[1]);
    trigonom_tab[3] = dir_tab[0] / ab;
    trigonom_tab[2] = dir_tab[1] / ab;
  } else {
    trigonom_tab[3] = 1.0;
    trigonom_tab[2] = 0.0;
    ab = 0.0;
  }
  abcmax = fmax(abmax, fabs(dir_tab[2]));
  if (abcmax > 1.0e-36) {
    abc = sqrt(ab * ab + dir_tab[2] * dir_tab[2]);
    trigonom_tab[1] = dir_tab[2] / abc;
    trigonom_tab[0] = ab / abc;
  } else {
    trigonom_tab[1] = 0.0;
    trigonom_tab[0] = 1.0;
  }
}
void tab_single_wave_ylm_p2_c_(double trigonom_tab[4], int *l_max_, int *l_ylm_max, double *ylm_tab) {
  int l_max = *l_max_;
  // if (fast_ylm) {
  SHEval_c_(l_max, trigonom_tab[0], trigonom_tab[1], trigonom_tab[2], trigonom_tab[3], ylm_tab);
  // } else {
  //   printf("%s, not finished\n", __func__); // TODO
  //   exit(-19);
  // }
}

void cubic_spline_v2_c_(double *spl_param, int *n_l_dim_, int *n_coeff_, int *n_grid_dim_, int *n_points_,
                        int *n_vector_) {
#define spl_param(i, j, k) spl_param[(i)-1 + n_l_dim * ((j)-1 + n_coeff * ((k)-1))]
#define d_inv(i) d_inv[(i)-1]
#define _MIN(i, j) ((i) < (j) ? (i) : (j))
  int n_l_dim = *n_l_dim_;
  int n_coeff = *n_coeff_;
  int n_grid_dim = *n_grid_dim_;
  int n_points = *n_points_;
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

void get_rho_multipole_spl_c_(double *rho_multipole_spl, int *spl_atom_) {
#define rho_multipole(i, j, k)                                                                                         \
  MV(hartree_potential_storage, rho_multipole_shmem_ptr)[(i)-1 + l_pot_max_help * ((j)-1 + n_max_radial_help * ((k)-1))]
#define rho_multipole_spl(i, j, k) rho_multipole_spl[(i)-1 + l_pot_max_help * ((j)-1 + n_max_spline * ((k)-1))]
  int spl_atom = *spl_atom_;
  int l_pot_max_help = (l_pot_max + 1) * (l_pot_max + 1);
  int n_max_radial_help = n_max_radial + 2;

  int i_atom_index = rho_multipole_index(spl_atom);
  int species_tmp = species(spl_atom);
  int l_h_dim = (l_hartree(species_tmp) + 1) * (l_hartree(species_tmp) + 1);
  int n_rad = n_radial(species_tmp);
  int n_rad_help = n_rad + 2;

  for (int k = 1; k <= (n_rad + 2); k++) {
    for (int i = 1; i <= l_h_dim; i++) {
      rho_multipole_spl(i, 1, k) = rho_multipole(i, k, i_atom_index);
    }
  }

  cubic_spline_v2_c_(rho_multipole_spl, &l_pot_max_help, &n_max_spline, &n_max_radial_help, &n_rad_help, &l_h_dim);

  int i_radial = n_rad;
  while ((r_radial(i_radial, species_tmp) >= multipole_radius_free(species_tmp)) && i_radial > 1) {
    for (int j = 1; j <= n_max_spline; j++)
      for (int i = 1; i <= l_h_dim; i++)
        rho_multipole_spl(i, j, i_radial + 1) = 0;
    i_radial--;
  }

  double i_r_outer = invert_radial_grid_c_(multipole_radius_free(species_tmp), n_rad, scale_radial(species_tmp));
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

void spline_angular_integral_log_c_(double *angular_part_spl, double *angular_integral_log, int *i_atom_) {
#define angular_part_spl(i, j, k) angular_part_spl[(i)-1 + l_pot_max_help * ((j)-1 + n_max_spline * ((k)-1))]
#define angular_integral_log(i, j) angular_integral_log[(i)-1 + l_pot_max_help * ((j)-1)]
  int l_pot_max_help = (l_pot_max + 1) * (l_pot_max + 1);
  int n_max_radial_help = n_max_radial + 2;
  int i_atom = *i_atom_;

  int l_h_dim = (l_hartree(species(i_atom)) + 1) * (l_hartree(species(i_atom)) + 1);

  int i_grid_max = n_grid(species(i_atom));
  for (int i_grid = 1; i_grid <= i_grid_max; i_grid++) {
    if (r_grid(i_grid, species(i_atom)) < multipole_radius_free(species(i_atom))) {
      double i_r_radial = 1 + invert_radial_grid_c_(r_grid(i_grid, species(i_atom)), n_radial(species(i_atom)),
                                                    scale_radial(species(i_atom)));
      spline_vector_v2_c_(i_r_radial, angular_part_spl, l_pot_max_help, n_max_spline, n_max_radial_help,
                          n_radial(species(i_atom)) + 2, l_h_dim, &angular_integral_log(1, i_grid));
    } else {
      for (int i = 1; i <= l_h_dim; i++)
        angular_integral_log(i, i_grid) = 0.0;
    }
  }

  if (compensate_multipole_errors) {
    double c_n_tmp = compensation_norm(i_atom);
    double c_r_tmp = compensation_radius(i_atom);
    int s_tmp = species(i_atom);
    int i_l = 0;
    for (int i_grid = 1; i_grid <= i_grid_max; i_grid++) {
      angular_integral_log(1, i_grid) += c_n_tmp * compensating_density_c_(r_grid(i_grid, s_tmp), c_r_tmp, i_l);
    }
  }

#undef angular_part_spl
#undef angular_integral_log
}

void integrate_delta_v_hartree_internal_c_(double *angular_integral_log, double *delta_v_hartree, int *n_coeff_hartree_,
                                           int *i_atom_) {
#define angular_integral_log(i, j) angular_integral_log[(i)-1 + l_pot_max_help * ((j)-1)]
#define delta_v_hartree(i, j, k) delta_v_hartree[(i)-1 + l_pot_max_help * ((j)-1 + n_coeff_hartree * ((k)-1))]
#define integral_zero_r(i) integral_zero_r[(i)-1]
#define integral_r_infty(i) integral_r_infty[(i)-1]
#define d_1_24 (1.0 / 24.0)
  int l_pot_max_help = (l_pot_max + 1) * (l_pot_max + 1);
  int n_coeff_hartree = *n_coeff_hartree_;
  int i_atom = *i_atom_;

  double prefactor[l_pot_max + 1];
  double integral_zero_r[(l_pot_max + 1) * (l_pot_max + 1)];
  double integral_r_infty[(l_pot_max + 1) * (l_pot_max + 1)];
  double r_l_AM[4];
  double r_neg_l1_AM[4];
  double r_inv_AM[4];
  double dr_coef_AM[4];

  for (int i_l = 0; i_l <= l_pot_max; i_l++)
    prefactor[i_l] = 12.56637061435917295376 / (2.0 * (double)i_l + 1.0); // pi4 = 12.56637061435917295376

  int l_h_dim = (l_hartree(species(i_atom)) + 1) * (l_hartree(species(i_atom)) + 1);
  double alpha = log(r_grid_inc(species(i_atom)));

  for (int i = 0; i < (l_pot_max + 1) * (l_pot_max + 1); i++)
    integral_zero_r[i] = 0;
  for (int i = 0; i < (l_pot_max + 1) * (l_pot_max + 1); i++)
    integral_r_infty[i] = 0;

  // TODO print Adams_Moulton_integrator
  // WARNING the other path may mot be checked
  int s_tmp = species(i_atom);
  int n_tmp = n_grid(s_tmp);
  if (Adams_Moulton_integrator) {
    int i_l_max = l_hartree(species(i_atom));
    int i_index = 0; // TODO: change to static form in the loop
    double r_grid1 = r_grid(1, species(i_atom));
    double r_grid2 = r_grid(2, species(i_atom));
    double r_grid3 = r_grid(3, species(i_atom));
    for (int i_l = 0; i_l <= i_l_max; i_l++) {
      for (int i_m = -i_l; i_m <= i_l; i_m++) {
        ++i_index;
        // TODO WARNING: 此处应为 pow 浮点数的整数次方，不知会不会使用到浮点数的浮点数次方版本
        // TODO: 配合这个循环来迭代次方理应更优
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
        for (int i_m = -i_l; i_m <= i_l; i_m++) {
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
    for (int i = 0; i < (l_pot_max + 1) * (l_pot_max + 1); i++)
      integral_r_infty[i] = 0.0;
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
    // all remaining terms
    // Integral_i = Integral_(i+1) + h[9 f_i + 19 f_(i+1) - 5 f_(i+2) + f_(i+3)]/24
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
  } else { // TODO NO CHECKED (without cases which Adams_Moulton_integrator=False)
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
  } // end if !(Adams_Moulton_integrator) and else
  // Calculate spline coefficients
  cubic_spline_v2_c_(delta_v_hartree, &l_pot_max_help, &n_coeff_hartree, &n_hartree_grid, &n_grid(species(i_atom)),
                     &l_h_dim);
#undef angular_integral_log
#undef delta_v_hartree
#undef integral_zero_r
#undef integral_r_infty
#undef d_1_24
}

void integrate_delta_v_hartree_c_(double *angular_part_spl, double *delta_v_hartree, int *n_coeff_hartree_,
                                  int *i_atom_) {
  double angular_integral_log[(l_pot_max + 1) * (l_pot_max + 1) * n_max_grid];
  spline_angular_integral_log_c_(angular_part_spl, angular_integral_log, i_atom_);
  // if(*i_atom_ == 1){
  //   m_save_check_angular_integral_log(angular_integral_log, (l_pot_max + 1) * (l_pot_max + 1), n_max_grid);
  // }
  integrate_delta_v_hartree_internal_c_(angular_integral_log, delta_v_hartree, n_coeff_hartree_, i_atom_);
}

void SHEval_c_(int lmax, double sintheta, double costheta, double sinphi, double cosphi, double *pSH) {
// intent(inout) :: pSH((lmax+1)*(lmax+1))
#define pSH(i) pSH[(i)-1]

  double fX, fY, fZ;
  double fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  // double fC0_1, fC1_1, fS0_1, fS1_1, fTmpA_1, fTmpB_1, fTmpC_1;
  // double fTmpA_2, fTmpB_2, fTmpC_2;
  double fZ2;
  fX = sintheta * cosphi;
  fY = sintheta * sinphi;
  fZ = costheta;
  fZ2 = fZ * fZ;
  switch (lmax) {
  case 0: {
    pSH(1) = 0.28209479177387814347;
  } break;
  case 1: {
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    fC0 = fX;
    fS0 = fY;
    fTmpB = -0.48860251190291992159;
    pSH(4) = fTmpB * fC0;
    pSH(2) = -fTmpB * fS0;
  } break;
  case 2: {
    // !m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    // !m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    // !m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.54627421529603953527;
    pSH(9) = fTmpC * fC1;
    pSH(5) = fTmpC * fS1;
  } break;
  case 3: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.59004358992664351035;
    pSH(16) = fTmpC * fC0;
    pSH(10) = -fTmpC * fS0;
  } break;
  case 4: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.62583573544917613459;
    pSH(25) = fTmpC * fC1;
    pSH(17) = fTmpC * fS1;
  } break;
  case 5: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.65638205684017010281;
    pSH(36) = fTmpC * fC0;
    pSH(26) = -fTmpC * fS0;
  } break;
  case 6: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.68318410519191432198;
    pSH(49) = fTmpC * fC1;
    pSH(37) = fTmpC * fS1;
  } break;
  case 7: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.70716273252459617823;
    pSH(64) = fTmpC * fC0;
    pSH(50) = -fTmpC * fS0;
  } break;
  case 8: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.72892666017482986887;
    pSH(81) = fTmpC * fC1;
    pSH(65) = fTmpC * fS1;
  } break;
  case 9: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.74890095185318829655;
    pSH(100) = fTmpC * fC0;
    pSH(82) = -fTmpC * fS0;
  } break;
  case 10: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.76739511822199001256;
    pSH(121) = fTmpC * fC1;
    pSH(101) = fTmpC * fS1;
  } break;
  case 11: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.78464210578719688375;
    pSH(144) = fTmpC * fC0;
    pSH(122) = -fTmpC * fS0;
  } break;
  case 12: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.80082199578397171410;
    pSH(169) = fTmpC * fC1;
    pSH(145) = fTmpC * fS1;
  } break;
  case 13: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.81607711883762830929;
    pSH(196) = fTmpC * fC0;
    pSH(170) = -fTmpC * fS0;
  } break;
  case 14: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.83052208306452400305;
    pSH(225) = fTmpC * fC1;
    pSH(197) = fTmpC * fS1;
  } break;
  case 15: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    pSH(241) = 1.9988885800753266487 * fZ * pSH(211) - 1.0000823011400101477 * pSH(183);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    fTmpC = 2.0033454163331038362 * fZ * fTmpB - 0.99975195336341716696 * fTmpA;
    pSH(242) = fTmpC * fC0;
    pSH(240) = -fTmpC * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    fTmpB = 2.0168969490698876097 * fZ * fTmpA - 0.99874229606879180776 * fTmpC;
    pSH(243) = fTmpB * fC1;
    pSH(239) = fTmpB * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    fTmpA = 2.0401071141087266541 * fZ * fTmpC - 0.99699479851094886098 * fTmpB;
    pSH(244) = fTmpA * fC0;
    pSH(238) = -fTmpA * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    fTmpC = 2.0739902137422357166 * fZ * fTmpB - 0.99440219512922986155 * fTmpA;
    pSH(245) = fTmpC * fC1;
    pSH(237) = fTmpC * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    fTmpB = 2.1201415047114190203 * fZ * fTmpA - 0.99079092984678996135 * fTmpC;
    pSH(246) = fTmpB * fC0;
    pSH(236) = -fTmpB * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    fTmpA = 2.1809662438042814911 * fZ * fTmpC - 0.98588907503509976255 * fTmpB;
    pSH(247) = fTmpA * fC1;
    pSH(235) = fTmpA * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    fTmpC = 2.2600784378986817489 * fZ * fTmpB - 0.97926740294193723593 * fTmpA;
    pSH(248) = fTmpC * fC0;
    pSH(234) = -fTmpC * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    fTmpB = 2.3630173363047971159 * fZ * fTmpA - 0.97022618722766530118 * fTmpC;
    pSH(249) = fTmpB * fC1;
    pSH(233) = fTmpB * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    fTmpA = 2.4986107250941583108 * fZ * fTmpC - 0.95756141751469769713 * fTmpB;
    pSH(250) = fTmpA * fC0;
    pSH(232) = -fTmpA * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    fTmpC = 2.6817904466978772501 * fZ * fTmpB - 0.93903023262181382111 * fTmpA;
    pSH(251) = fTmpC * fC1;
    pSH(231) = fTmpC * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    fTmpB = 2.9401072717216916592 * fZ * fTmpA - 0.90994035683194807483 * fTmpC;
    pSH(252) = fTmpB * fC0;
    pSH(230) = -fTmpB * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    fTmpA = fZ * (53.794072123058085929 * fZ2 - 5.5649040127301468202);
    pSH(253) = fTmpA * fC1;
    pSH(229) = fTmpA * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    fTmpC = -17.608243388844820556 * fZ2 + 0.60718080651189036399;
    pSH(254) = fTmpC * fC0;
    pSH(228) = -fTmpC * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.83052208306452400305;
    pSH(225) = fTmpA * fC1;
    pSH(197) = fTmpA * fS1;
    fTmpB = 4.6241512566300114788 * fZ;
    pSH(255) = fTmpB * fC1;
    pSH(227) = fTmpB * fS1;
    //! m = 15
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.84425065085737263613;
    pSH(256) = fTmpC * fC0;
    pSH(226) = -fTmpC * fS0;
  } break;
  case 16: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    pSH(241) = 1.9988885800753266487 * fZ * pSH(211) - 1.0000823011400101477 * pSH(183);
    pSH(273) = 1.9990231989649344737 * fZ * pSH(241) - 1.0000673468701305763 * pSH(211);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    fTmpC = 2.0033454163331038362 * fZ * fTmpB - 0.99975195336341716696 * fTmpA;
    pSH(242) = fTmpC * fC0;
    pSH(240) = -fTmpC * fS0;
    fTmpA = 2.0029390170153341294 * fZ * fTmpC - 0.99979713966725040627 * fTmpB;
    pSH(274) = fTmpA * fC0;
    pSH(272) = -fTmpA * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    fTmpB = 2.0168969490698876097 * fZ * fTmpA - 0.99874229606879180776 * fTmpC;
    pSH(243) = fTmpB * fC1;
    pSH(239) = fTmpB * fS1;
    fTmpC = 2.0148259998133361203 * fZ * fTmpB - 0.99897320026315349004 * fTmpA;
    pSH(275) = fTmpC * fC1;
    pSH(271) = fTmpC * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    fTmpA = 2.0401071141087266541 * fZ * fTmpC - 0.99699479851094886098 * fTmpB;
    pSH(244) = fTmpA * fC0;
    pSH(238) = -fTmpA * fS0;
    fTmpB = 2.0351168037383750115 * fZ * fTmpA - 0.99755389786357772286 * fTmpC;
    pSH(276) = fTmpB * fC0;
    pSH(270) = -fTmpB * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    fTmpC = 2.0739902137422357166 * fZ * fTmpB - 0.99440219512922986155 * fTmpA;
    pSH(245) = fTmpC * fC1;
    pSH(237) = fTmpC * fS1;
    fTmpA = 2.0645822822062578187 * fZ * fTmpC - 0.99546384960081245811 * fTmpB;
    pSH(277) = fTmpA * fC1;
    pSH(269) = fTmpA * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    fTmpB = 2.1201415047114190203 * fZ * fTmpA - 0.99079092984678996135 * fTmpC;
    pSH(246) = fTmpB * fC0;
    pSH(236) = -fTmpB * fS0;
    fTmpC = 2.1044171232366050512 * fZ * fTmpB - 0.99258333397093026682 * fTmpA;
    pSH(278) = fTmpC * fC0;
    pSH(268) = -fTmpC * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    fTmpA = 2.1809662438042814911 * fZ * fTmpC - 0.98588907503509976255 * fTmpB;
    pSH(247) = fTmpA * fC1;
    pSH(235) = fTmpA * fS1;
    fTmpB = 2.1563858652847824675 * fZ * fTmpA - 0.98872959240459255313 * fTmpC;
    pSH(279) = fTmpB * fC1;
    pSH(267) = fTmpB * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    fTmpC = 2.2600784378986817489 * fZ * fTmpB - 0.97926740294193723593 * fTmpA;
    pSH(248) = fTmpC * fC0;
    pSH(234) = -fTmpC * fS0;
    fTmpA = 2.2230674720995866294 * fZ * fTmpC - 0.98362403482177094968 * fTmpB;
    pSH(280) = fTmpA * fC0;
    pSH(266) = -fTmpA * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    fTmpB = 2.3630173363047971159 * fZ * fTmpA - 0.97022618722766530118 * fTmpC;
    pSH(249) = fTmpB * fC1;
    pSH(233) = fTmpB * fS1;
    fTmpC = 2.3082731640774234847 * fZ * fTmpB - 0.97683293669229670981 * fTmpA;
    pSH(281) = fTmpC * fC1;
    pSH(265) = fTmpC * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    fTmpA = 2.4986107250941583108 * fZ * fTmpC - 0.95756141751469769713 * fTmpB;
    pSH(250) = fTmpA * fC0;
    pSH(232) = -fTmpA * fS0;
    fTmpB = 2.4177911997760033443 * fZ * fTmpA - 0.96765421499777267545 * fTmpC;
    pSH(282) = fTmpB * fC0;
    pSH(264) = -fTmpB * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    fTmpC = 2.6817904466978772501 * fZ * fTmpB - 0.93903023262181382111 * fTmpA;
    pSH(251) = fTmpC * fC1;
    pSH(231) = fTmpC * fS1;
    fTmpA = 2.5607991541103546085 * fZ * fTmpC - 0.95488413617980451759 * fTmpB;
    pSH(283) = fTmpA * fC1;
    pSH(263) = fTmpA * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    fTmpB = 2.9401072717216916592 * fZ * fTmpA - 0.90994035683194807483 * fTmpC;
    pSH(252) = fTmpB * fC0;
    pSH(230) = -fTmpB * fS0;
    fTmpC = 2.7527763762750104249 * fZ * fTmpB - 0.93628433314374189412 * fTmpA;
    pSH(284) = fTmpC * fC0;
    pSH(262) = -fTmpC * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    fTmpA = fZ * (53.794072123058085929 * fZ2 - 5.5649040127301468202);
    pSH(253) = fTmpA * fC1;
    pSH(229) = fTmpA * fS1;
    fTmpB = 3.0222389997200041804 * fZ * fTmpA - 0.90717582656041188324 * fTmpC;
    pSH(285) = fTmpB * fC1;
    pSH(261) = fTmpB * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    fTmpC = -17.608243388844820556 * fZ2 + 0.60718080651189036399;
    pSH(254) = fTmpC * fC0;
    pSH(228) = -fTmpC * fS0;
    fTmpA = fZ * (-60.380154942952015813 * fZ2 + 5.8432408009308402400);
    pSH(286) = fTmpA * fC0;
    pSH(260) = -fTmpA * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.83052208306452400305;
    pSH(225) = fTmpA * fC1;
    pSH(197) = fTmpA * fS1;
    fTmpB = 4.6241512566300114788 * fZ;
    pSH(255) = fTmpB * fC1;
    pSH(227) = fTmpB * fS1;
    fTmpC = 19.093881509360250421 * fZ2 - 0.61593166159226614262;
    pSH(287) = fTmpC * fC1;
    pSH(259) = fTmpC * fS1;
    //! m = 15
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.84425065085737263613;
    pSH(256) = fTmpA * fC0;
    pSH(226) = -fTmpA * fS0;
    fTmpB = -4.8498507532306824834 * fZ;
    pSH(288) = fTmpB * fC0;
    pSH(258) = -fTmpB * fS0;
    //! m = 16
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.85734058883802509634;
    pSH(289) = fTmpC * fC1;
    pSH(257) = fTmpC * fS1;
  } break;
  case 17: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    pSH(241) = 1.9988885800753266487 * fZ * pSH(211) - 1.0000823011400101477 * pSH(183);
    pSH(273) = 1.9990231989649344737 * fZ * pSH(241) - 1.0000673468701305763 * pSH(211);
    pSH(307) = 1.9991347609372268760 * fZ * pSH(273) - 1.0000558082429209263 * pSH(241);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    fTmpC = 2.0033454163331038362 * fZ * fTmpB - 0.99975195336341716696 * fTmpA;
    pSH(242) = fTmpC * fC0;
    pSH(240) = -fTmpC * fS0;
    fTmpA = 2.0029390170153341294 * fZ * fTmpC - 0.99979713966725040627 * fTmpB;
    pSH(274) = fTmpA * fC0;
    pSH(272) = -fTmpA * fS0;
    fTmpB = 2.0026024734496526300 * fZ * fTmpA - 0.99983197513113354916 * fTmpC;
    pSH(308) = fTmpB * fC0;
    pSH(306) = -fTmpB * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    fTmpB = 2.0168969490698876097 * fZ * fTmpA - 0.99874229606879180776 * fTmpC;
    pSH(243) = fTmpB * fC1;
    pSH(239) = fTmpB * fS1;
    fTmpC = 2.0148259998133361203 * fZ * fTmpB - 0.99897320026315349004 * fTmpA;
    pSH(275) = fTmpC * fC1;
    pSH(271) = fTmpC * fS1;
    fTmpA = 2.0131148946216081338 * fZ * fTmpC - 0.99915074294659364525 * fTmpB;
    pSH(309) = fTmpA * fC1;
    pSH(305) = fTmpA * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    fTmpA = 2.0401071141087266541 * fZ * fTmpC - 0.99699479851094886098 * fTmpB;
    pSH(244) = fTmpA * fC0;
    pSH(238) = -fTmpA * fS0;
    fTmpB = 2.0351168037383750115 * fZ * fTmpA - 0.99755389786357772286 * fTmpC;
    pSH(276) = fTmpB * fC0;
    pSH(270) = -fTmpB * fS0;
    fTmpC = 2.0310096011589900901 * fZ * fTmpB - 0.99798183447169211035 * fTmpA;
    pSH(310) = fTmpC * fC0;
    pSH(304) = -fTmpC * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    fTmpC = 2.0739902137422357166 * fZ * fTmpB - 0.99440219512922986155 * fTmpA;
    pSH(245) = fTmpC * fC1;
    pSH(237) = fTmpC * fS1;
    fTmpA = 2.0645822822062578187 * fZ * fTmpC - 0.99546384960081245811 * fTmpB;
    pSH(277) = fTmpA * fC1;
    pSH(269) = fTmpA * fS1;
    fTmpB = 2.0568833780186057912 * fZ * fTmpA - 0.99627096277343579159 * fTmpC;
    pSH(311) = fTmpB * fC1;
    pSH(303) = fTmpB * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    fTmpB = 2.1201415047114190203 * fZ * fTmpA - 0.99079092984678996135 * fTmpC;
    pSH(246) = fTmpB * fC0;
    pSH(236) = -fTmpB * fS0;
    fTmpC = 2.1044171232366050512 * fZ * fTmpB - 0.99258333397093026682 * fTmpA;
    pSH(278) = fTmpC * fC0;
    pSH(268) = -fTmpC * fS0;
    fTmpA = 2.0916500663351888699 * fZ * fTmpC - 0.99393320993236341952 * fTmpB;
    pSH(312) = fTmpA * fC0;
    pSH(302) = -fTmpA * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    fTmpA = 2.1809662438042814911 * fZ * fTmpC - 0.98588907503509976255 * fTmpB;
    pSH(247) = fTmpA * fC1;
    pSH(235) = fTmpA * fS1;
    fTmpB = 2.1563858652847824675 * fZ * fTmpA - 0.98872959240459255313 * fTmpC;
    pSH(279) = fTmpB * fC1;
    pSH(267) = fTmpB * fS1;
    fTmpC = 2.1366369348357590877 * fZ * fTmpB - 0.99084165280112553369 * fTmpA;
    pSH(313) = fTmpC * fC1;
    pSH(301) = fTmpC * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    fTmpC = 2.2600784378986817489 * fZ * fTmpB - 0.97926740294193723593 * fTmpA;
    pSH(248) = fTmpC * fC0;
    pSH(234) = -fTmpC * fS0;
    fTmpA = 2.2230674720995866294 * fZ * fTmpC - 0.98362403482177094968 * fTmpB;
    pSH(280) = fTmpA * fC0;
    pSH(266) = -fTmpA * fS0;
    fTmpB = 2.1937410968480305151 * fZ * fTmpA - 0.98680814882156560011 * fTmpC;
    pSH(314) = fTmpB * fC0;
    pSH(300) = -fTmpB * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    fTmpB = 2.3630173363047971159 * fZ * fTmpA - 0.97022618722766530118 * fTmpC;
    pSH(249) = fTmpB * fC1;
    pSH(233) = fTmpB * fS1;
    fTmpC = 2.3082731640774234847 * fZ * fTmpB - 0.97683293669229670981 * fTmpA;
    pSH(281) = fTmpC * fC1;
    pSH(265) = fTmpC * fS1;
    fTmpA = 2.2656860623955237928 * fZ * fTmpC - 0.98155023315928857429 * fTmpB;
    pSH(315) = fTmpA * fC1;
    pSH(299) = fTmpA * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    fTmpA = 2.4986107250941583108 * fZ * fTmpC - 0.95756141751469769713 * fTmpB;
    pSH(250) = fTmpA * fC0;
    pSH(232) = -fTmpA * fS0;
    fTmpB = 2.4177911997760033443 * fZ * fTmpA - 0.96765421499777267545 * fTmpC;
    pSH(282) = fTmpB * fC0;
    pSH(264) = -fTmpB * fS0;
    fTmpC = 2.3564559438666820500 * fZ * fTmpB - 0.97463169858712211832 * fTmpA;
    pSH(316) = fTmpC * fC0;
    pSH(298) = -fTmpC * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    fTmpC = 2.6817904466978772501 * fZ * fTmpB - 0.93903023262181382111 * fTmpA;
    pSH(251) = fTmpC * fC1;
    pSH(231) = fTmpC * fS1;
    fTmpA = 2.5607991541103546085 * fZ * fTmpC - 0.95488413617980451759 * fTmpB;
    pSH(283) = fTmpA * fC1;
    pSH(263) = fTmpA * fS1;
    fTmpB = 2.4720661623652209829 * fZ * fTmpA - 0.96534949193391143148 * fTmpC;
    pSH(317) = fTmpB * fC1;
    pSH(297) = fTmpB * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    fTmpB = 2.9401072717216916592 * fZ * fTmpA - 0.90994035683194807483 * fTmpC;
    pSH(252) = fTmpB * fC0;
    pSH(230) = -fTmpB * fS0;
    fTmpC = 2.7527763762750104249 * fZ * fTmpB - 0.93628433314374189412 * fTmpA;
    pSH(284) = fTmpC * fC0;
    pSH(262) = -fTmpC * fS0;
    fTmpA = 2.6220221204253788675 * fZ * fTmpC - 0.95250095250142875238 * fTmpB;
    pSH(318) = fTmpA * fC0;
    pSH(296) = -fTmpA * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    fTmpA = fZ * (53.794072123058085929 * fZ2 - 5.5649040127301468202);
    pSH(253) = fTmpA * fC1;
    pSH(229) = fTmpA * fS1;
    fTmpB = 3.0222389997200041804 * fZ * fTmpA - 0.90717582656041188324 * fTmpC;
    pSH(285) = fTmpB * fC1;
    pSH(261) = fTmpB * fS1;
    fTmpC = 2.8223247937435036367 * fZ * fTmpB - 0.93385228435109810060 * fTmpA;
    pSH(319) = fTmpC * fC1;
    pSH(295) = fTmpC * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    fTmpC = -17.608243388844820556 * fZ2 + 0.60718080651189036399;
    pSH(254) = fTmpC * fC0;
    pSH(228) = -fTmpC * fS0;
    fTmpA = fZ * (-60.380154942952015813 * fZ2 + 5.8432408009308402400);
    pSH(286) = fTmpA * fC0;
    pSH(260) = -fTmpA * fS0;
    fTmpB = 3.1024184114977141490 * fZ * fTmpA - 0.90473663963430495830 * fTmpC;
    pSH(320) = fTmpB * fC0;
    pSH(294) = -fTmpB * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.83052208306452400305;
    pSH(225) = fTmpA * fC1;
    pSH(197) = fTmpA * fS1;
    fTmpB = 4.6241512566300114788 * fZ;
    pSH(255) = fTmpB * fC1;
    pSH(227) = fTmpB * fS1;
    fTmpC = 19.093881509360250421 * fZ2 - 0.61593166159226614262;
    pSH(287) = fTmpC * fC1;
    pSH(259) = fTmpC * fS1;
    fTmpA = fZ * (67.288948373844056316 * fZ2 - 6.1171771248949142105);
    pSH(321) = fTmpA * fC1;
    pSH(293) = fTmpA * fS1;
    //! m = 15
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.84425065085737263613;
    pSH(256) = fTmpA * fC0;
    pSH(226) = -fTmpA * fS0;
    fTmpB = -4.8498507532306824834 * fZ;
    pSH(288) = fTmpB * fC0;
    pSH(258) = -fTmpB * fS0;
    fTmpC = -20.602948605549728459 * fZ2 + 0.62433177592574934724;
    pSH(322) = fTmpC * fC0;
    pSH(292) = -fTmpC * fS0;
    //! m = 16
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.85734058883802509634;
    pSH(289) = fTmpA * fC1;
    pSH(257) = fTmpA * fS1;
    fTmpB = 5.0720953248553606107 * fZ;
    pSH(323) = fTmpB * fC1;
    pSH(291) = fTmpB * fS1;
    //! m = 17
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.86985717192062808222;
    pSH(324) = fTmpC * fC0;
    pSH(290) = -fTmpC * fS0;
  } break;
  case 18: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    pSH(241) = 1.9988885800753266487 * fZ * pSH(211) - 1.0000823011400101477 * pSH(183);
    pSH(273) = 1.9990231989649344737 * fZ * pSH(241) - 1.0000673468701305763 * pSH(211);
    pSH(307) = 1.9991347609372268760 * fZ * pSH(273) - 1.0000558082429209263 * pSH(241);
    pSH(343) = 1.9992282461607312886 * fZ * pSH(307) - 1.0000467628422711159 * pSH(273);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    fTmpC = 2.0033454163331038362 * fZ * fTmpB - 0.99975195336341716696 * fTmpA;
    pSH(242) = fTmpC * fC0;
    pSH(240) = -fTmpC * fS0;
    fTmpA = 2.0029390170153341294 * fZ * fTmpC - 0.99979713966725040627 * fTmpB;
    pSH(274) = fTmpA * fC0;
    pSH(272) = -fTmpA * fS0;
    fTmpB = 2.0026024734496526300 * fZ * fTmpA - 0.99983197513113354916 * fTmpC;
    pSH(308) = fTmpB * fC0;
    pSH(306) = -fTmpB * fS0;
    fTmpC = 2.0023206350873464509 * fZ * fTmpB - 0.99985926394976398460 * fTmpA;
    pSH(344) = fTmpC * fC0;
    pSH(342) = -fTmpC * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    fTmpB = 2.0168969490698876097 * fZ * fTmpA - 0.99874229606879180776 * fTmpC;
    pSH(243) = fTmpB * fC1;
    pSH(239) = fTmpB * fS1;
    fTmpC = 2.0148259998133361203 * fZ * fTmpB - 0.99897320026315349004 * fTmpA;
    pSH(275) = fTmpC * fC1;
    pSH(271) = fTmpC * fS1;
    fTmpA = 2.0131148946216081338 * fZ * fTmpC - 0.99915074294659364525 * fTmpB;
    pSH(309) = fTmpA * fC1;
    pSH(305) = fTmpA * fS1;
    fTmpB = 2.0116846174288851484 * fZ * fTmpA - 0.99928952033659667242 * fTmpC;
    pSH(345) = fTmpB * fC1;
    pSH(341) = fTmpB * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    fTmpA = 2.0401071141087266541 * fZ * fTmpC - 0.99699479851094886098 * fTmpB;
    pSH(244) = fTmpA * fC0;
    pSH(238) = -fTmpA * fS0;
    fTmpB = 2.0351168037383750115 * fZ * fTmpA - 0.99755389786357772286 * fTmpC;
    pSH(276) = fTmpB * fC0;
    pSH(270) = -fTmpB * fS0;
    fTmpC = 2.0310096011589900901 * fZ * fTmpB - 0.99798183447169211035 * fTmpA;
    pSH(310) = fTmpC * fC0;
    pSH(304) = -fTmpC * fS0;
    fTmpA = 2.0275875100994065630 * fZ * fTmpC - 0.99831507883683527632 * fTmpB;
    pSH(346) = fTmpA * fC0;
    pSH(340) = -fTmpA * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    fTmpC = 2.0739902137422357166 * fZ * fTmpB - 0.99440219512922986155 * fTmpA;
    pSH(245) = fTmpC * fC1;
    pSH(237) = fTmpC * fS1;
    fTmpA = 2.0645822822062578187 * fZ * fTmpC - 0.99546384960081245811 * fTmpB;
    pSH(277) = fTmpA * fC1;
    pSH(269) = fTmpA * fS1;
    fTmpB = 2.0568833780186057912 * fZ * fTmpA - 0.99627096277343579159 * fTmpC;
    pSH(311) = fTmpB * fC1;
    pSH(303) = fTmpB * fS1;
    fTmpC = 2.0504988306618110689 * fZ * fTmpB - 0.99689600906642312802 * fTmpA;
    pSH(347) = fTmpC * fC1;
    pSH(339) = fTmpC * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    fTmpB = 2.1201415047114190203 * fZ * fTmpA - 0.99079092984678996135 * fTmpC;
    pSH(246) = fTmpB * fC0;
    pSH(236) = -fTmpB * fS0;
    fTmpC = 2.1044171232366050512 * fZ * fTmpB - 0.99258333397093026682 * fTmpA;
    pSH(278) = fTmpC * fC0;
    pSH(268) = -fTmpC * fS0;
    fTmpA = 2.0916500663351888699 * fZ * fTmpC - 0.99393320993236341952 * fTmpB;
    pSH(312) = fTmpA * fC0;
    pSH(302) = -fTmpA * fS0;
    fTmpB = 2.0811303848941723349 * fZ * fTmpA - 0.99497063031224519069 * fTmpC;
    pSH(348) = fTmpB * fC0;
    pSH(338) = -fTmpB * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    fTmpA = 2.1809662438042814911 * fZ * fTmpC - 0.98588907503509976255 * fTmpB;
    pSH(247) = fTmpA * fC1;
    pSH(235) = fTmpA * fS1;
    fTmpB = 2.1563858652847824675 * fZ * fTmpA - 0.98872959240459255313 * fTmpC;
    pSH(279) = fTmpB * fC1;
    pSH(267) = fTmpB * fS1;
    fTmpC = 2.1366369348357590877 * fZ * fTmpB - 0.99084165280112553369 * fTmpA;
    pSH(313) = fTmpC * fC1;
    pSH(301) = fTmpC * fS1;
    fTmpA = 2.1205017749999120851 * fZ * fTmpC - 0.99244833805276922520 * fTmpB;
    pSH(349) = fTmpA * fC1;
    pSH(337) = fTmpA * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    fTmpC = 2.2600784378986817489 * fZ * fTmpB - 0.97926740294193723593 * fTmpA;
    pSH(248) = fTmpC * fC0;
    pSH(234) = -fTmpC * fS0;
    fTmpA = 2.2230674720995866294 * fZ * fTmpC - 0.98362403482177094968 * fTmpB;
    pSH(280) = fTmpA * fC0;
    pSH(266) = -fTmpA * fS0;
    fTmpB = 2.1937410968480305151 * fZ * fTmpA - 0.98680814882156560011 * fTmpC;
    pSH(314) = fTmpB * fC0;
    pSH(300) = -fTmpB * fS0;
    fTmpC = 2.1700439878239586254 * fZ * fTmpB - 0.98919785518075951600 * fTmpA;
    pSH(350) = fTmpC * fC0;
    pSH(336) = -fTmpC * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    fTmpB = 2.3630173363047971159 * fZ * fTmpA - 0.97022618722766530118 * fTmpC;
    pSH(249) = fTmpB * fC1;
    pSH(233) = fTmpB * fS1;
    fTmpC = 2.3082731640774234847 * fZ * fTmpB - 0.97683293669229670981 * fTmpA;
    pSH(281) = fTmpC * fC1;
    pSH(265) = fTmpC * fS1;
    fTmpA = 2.2656860623955237928 * fZ * fTmpC - 0.98155023315928857429 * fTmpB;
    pSH(315) = fTmpA * fC1;
    pSH(299) = fTmpA * fS1;
    fTmpB = 2.2317637040621551360 * fZ * fTmpA - 0.98502777640009739967 * fTmpC;
    pSH(351) = fTmpB * fC1;
    pSH(335) = fTmpB * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    fTmpA = 2.4986107250941583108 * fZ * fTmpC - 0.95756141751469769713 * fTmpB;
    pSH(250) = fTmpA * fC0;
    pSH(232) = -fTmpA * fS0;
    fTmpB = 2.4177911997760033443 * fZ * fTmpA - 0.96765421499777267545 * fTmpC;
    pSH(282) = fTmpB * fC0;
    pSH(264) = -fTmpB * fS0;
    fTmpC = 2.3564559438666820500 * fZ * fTmpB - 0.97463169858712211832 * fTmpA;
    pSH(316) = fTmpC * fC0;
    pSH(298) = -fTmpC * fS0;
    fTmpA = 2.3085099321848032226 * fZ * fTmpC - 0.97965333839290678348 * fTmpB;
    pSH(352) = fTmpA * fC0;
    pSH(334) = -fTmpA * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    fTmpC = 2.6817904466978772501 * fZ * fTmpB - 0.93903023262181382111 * fTmpA;
    pSH(251) = fTmpC * fC1;
    pSH(231) = fTmpC * fS1;
    fTmpA = 2.5607991541103546085 * fZ * fTmpC - 0.95488413617980451759 * fTmpB;
    pSH(283) = fTmpA * fC1;
    pSH(263) = fTmpA * fS1;
    fTmpB = 2.4720661623652209829 * fZ * fTmpA - 0.96534949193391143148 * fTmpC;
    pSH(317) = fTmpB * fC1;
    pSH(297) = fTmpB * fS1;
    fTmpC = 2.4044230077089180940 * fZ * fTmpB - 0.97263699666048446727 * fTmpA;
    pSH(353) = fTmpC * fC1;
    pSH(333) = fTmpC * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    fTmpB = 2.9401072717216916592 * fZ * fTmpA - 0.90994035683194807483 * fTmpC;
    pSH(252) = fTmpB * fC0;
    pSH(230) = -fTmpB * fS0;
    fTmpC = 2.7527763762750104249 * fZ * fTmpB - 0.93628433314374189412 * fTmpA;
    pSH(284) = fTmpC * fC0;
    pSH(262) = -fTmpC * fS0;
    fTmpA = 2.6220221204253788675 * fZ * fTmpC - 0.95250095250142875238 * fTmpB;
    pSH(318) = fTmpA * fC0;
    pSH(296) = -fTmpA * fS0;
    fTmpB = 2.5257296658248257996 * fZ * fTmpA - 0.96327549876469720969 * fTmpC;
    pSH(354) = fTmpB * fC0;
    pSH(332) = -fTmpB * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    fTmpA = fZ * (53.794072123058085929 * fZ2 - 5.5649040127301468202);
    pSH(253) = fTmpA * fC1;
    pSH(229) = fTmpA * fS1;
    fTmpB = 3.0222389997200041804 * fZ * fTmpA - 0.90717582656041188324 * fTmpC;
    pSH(285) = fTmpB * fC1;
    pSH(261) = fTmpB * fS1;
    fTmpC = 2.8223247937435036367 * fZ * fTmpB - 0.93385228435109810060 * fTmpA;
    pSH(319) = fTmpC * fC1;
    pSH(295) = fTmpC * fS1;
    fTmpA = 2.6822461565718468645 * fZ * fTmpC - 0.95036764107299718112 * fTmpB;
    pSH(355) = fTmpA * fC1;
    pSH(331) = fTmpA * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    fTmpC = -17.608243388844820556 * fZ2 + 0.60718080651189036399;
    pSH(254) = fTmpC * fC0;
    pSH(228) = -fTmpC * fS0;
    fTmpA = fZ * (-60.380154942952015813 * fZ2 + 5.8432408009308402400);
    pSH(286) = fTmpA * fC0;
    pSH(260) = -fTmpA * fS0;
    fTmpB = 3.1024184114977141490 * fZ * fTmpA - 0.90473663963430495830 * fTmpC;
    pSH(320) = fTmpB * fC0;
    pSH(294) = -fTmpB * fS0;
    fTmpC = 2.8904737863674562919 * fZ * fTmpB - 0.93168406158731500387 * fTmpA;
    pSH(356) = fTmpC * fC0;
    pSH(330) = -fTmpC * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.83052208306452400305;
    pSH(225) = fTmpA * fC1;
    pSH(197) = fTmpA * fS1;
    fTmpB = 4.6241512566300114788 * fZ;
    pSH(255) = fTmpB * fC1;
    pSH(227) = fTmpB * fS1;
    fTmpC = 19.093881509360250421 * fZ2 - 0.61593166159226614262;
    pSH(287) = fTmpC * fC1;
    pSH(259) = fTmpC * fS1;
    fTmpA = fZ * (67.288948373844056316 * fZ2 - 6.1171771248949142105);
    pSH(321) = fTmpA * fC1;
    pSH(293) = fTmpA * fS1;
    fTmpB = 3.1807526624998681277 * fZ * fTmpA - 0.90256893466271141000 * fTmpC;
    pSH(357) = fTmpB * fC1;
    pSH(329) = fTmpB * fS1;
    //! m = 15
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.84425065085737263613;
    pSH(256) = fTmpA * fC0;
    pSH(226) = -fTmpA * fS0;
    fTmpB = -4.8498507532306824834 * fZ;
    pSH(288) = fTmpB * fC0;
    pSH(258) = -fTmpB * fS0;
    fTmpC = -20.602948605549728459 * fZ2 + 0.62433177592574934724;
    pSH(322) = fTmpC * fC0;
    pSH(292) = -fTmpC * fS0;
    fTmpA = fZ * (-74.515507921532000473 * fZ2 + 6.3870435361313143263);
    pSH(358) = fTmpA * fC0;
    pSH(328) = -fTmpA * fS0;
    //! m = 16
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.85734058883802509634;
    pSH(289) = fTmpA * fC1;
    pSH(257) = fTmpA * fS1;
    fTmpB = 5.0720953248553606107 * fZ;
    pSH(323) = fTmpB * fC1;
    pSH(291) = fTmpB * fS1;
    fTmpC = 22.134404124649146409 * fZ2 - 0.63241154641854704026;
    pSH(359) = fTmpC * fC1;
    pSH(327) = fTmpC * fS1;
    //! m = 17
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.86985717192062808222;
    pSH(324) = fTmpA * fC0;
    pSH(290) = -fTmpA * fS0;
    fTmpB = -5.2911346120699725413 * fZ;
    pSH(360) = fTmpB * fC0;
    pSH(326) = -fTmpB * fS0;
    //! m = 18
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.88185576867832886131;
    pSH(361) = fTmpC * fC1;
    pSH(325) = fTmpC * fS1;
  } break;
  case 19: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    pSH(241) = 1.9988885800753266487 * fZ * pSH(211) - 1.0000823011400101477 * pSH(183);
    pSH(273) = 1.9990231989649344737 * fZ * pSH(241) - 1.0000673468701305763 * pSH(211);
    pSH(307) = 1.9991347609372268760 * fZ * pSH(273) - 1.0000558082429209263 * pSH(241);
    pSH(343) = 1.9992282461607312886 * fZ * pSH(307) - 1.0000467628422711159 * pSH(273);
    pSH(381) = 1.9993073592865872621 * fZ * pSH(343) - 1.0000395718327849261 * pSH(307);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    fTmpC = 2.0033454163331038362 * fZ * fTmpB - 0.99975195336341716696 * fTmpA;
    pSH(242) = fTmpC * fC0;
    pSH(240) = -fTmpC * fS0;
    fTmpA = 2.0029390170153341294 * fZ * fTmpC - 0.99979713966725040627 * fTmpB;
    pSH(274) = fTmpA * fC0;
    pSH(272) = -fTmpA * fS0;
    fTmpB = 2.0026024734496526300 * fZ * fTmpA - 0.99983197513113354916 * fTmpC;
    pSH(308) = fTmpB * fC0;
    pSH(306) = -fTmpB * fS0;
    fTmpC = 2.0023206350873464509 * fZ * fTmpB - 0.99985926394976398460 * fTmpA;
    pSH(344) = fTmpC * fC0;
    pSH(342) = -fTmpC * fS0;
    fTmpA = 2.0020822493926999835 * fZ * fTmpC - 0.99988094529394086354 * fTmpB;
    pSH(382) = fTmpA * fC0;
    pSH(380) = -fTmpA * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    fTmpB = 2.0168969490698876097 * fZ * fTmpA - 0.99874229606879180776 * fTmpC;
    pSH(243) = fTmpB * fC1;
    pSH(239) = fTmpB * fS1;
    fTmpC = 2.0148259998133361203 * fZ * fTmpB - 0.99897320026315349004 * fTmpA;
    pSH(275) = fTmpC * fC1;
    pSH(271) = fTmpC * fS1;
    fTmpA = 2.0131148946216081338 * fZ * fTmpC - 0.99915074294659364525 * fTmpB;
    pSH(309) = fTmpA * fC1;
    pSH(305) = fTmpA * fS1;
    fTmpB = 2.0116846174288851484 * fZ * fTmpA - 0.99928952033659667242 * fTmpC;
    pSH(345) = fTmpB * fC1;
    pSH(341) = fTmpB * fS1;
    fTmpC = 2.0104767610501468006 * fZ * fTmpB - 0.99939957965166423681 * fTmpA;
    pSH(383) = fTmpC * fC1;
    pSH(379) = fTmpC * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    fTmpA = 2.0401071141087266541 * fZ * fTmpC - 0.99699479851094886098 * fTmpB;
    pSH(244) = fTmpA * fC0;
    pSH(238) = -fTmpA * fS0;
    fTmpB = 2.0351168037383750115 * fZ * fTmpA - 0.99755389786357772286 * fTmpC;
    pSH(276) = fTmpB * fC0;
    pSH(270) = -fTmpB * fS0;
    fTmpC = 2.0310096011589900901 * fZ * fTmpB - 0.99798183447169211035 * fTmpA;
    pSH(310) = fTmpC * fC0;
    pSH(304) = -fTmpC * fS0;
    fTmpA = 2.0275875100994065630 * fZ * fTmpC - 0.99831507883683527632 * fTmpB;
    pSH(346) = fTmpA * fC0;
    pSH(340) = -fTmpA * fS0;
    fTmpB = 2.0247053657709850061 * fZ * fTmpA - 0.99857853517341885086 * fTmpC;
    pSH(384) = fTmpB * fC0;
    pSH(378) = -fTmpB * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    fTmpC = 2.0739902137422357166 * fZ * fTmpB - 0.99440219512922986155 * fTmpA;
    pSH(245) = fTmpC * fC1;
    pSH(237) = fTmpC * fS1;
    fTmpA = 2.0645822822062578187 * fZ * fTmpC - 0.99546384960081245811 * fTmpB;
    pSH(277) = fTmpA * fC1;
    pSH(269) = fTmpA * fS1;
    fTmpB = 2.0568833780186057912 * fZ * fTmpA - 0.99627096277343579159 * fTmpC;
    pSH(311) = fTmpB * fC1;
    pSH(303) = fTmpB * fS1;
    fTmpC = 2.0504988306618110689 * fZ * fTmpB - 0.99689600906642312802 * fTmpA;
    pSH(347) = fTmpC * fC1;
    pSH(339) = fTmpC * fS1;
    fTmpA = 2.0451427078940417827 * fZ * fTmpC - 0.99738789279580297796 * fTmpB;
    pSH(385) = fTmpA * fC1;
    pSH(377) = fTmpA * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    fTmpB = 2.1201415047114190203 * fZ * fTmpA - 0.99079092984678996135 * fTmpC;
    pSH(246) = fTmpB * fC0;
    pSH(236) = -fTmpB * fS0;
    fTmpC = 2.1044171232366050512 * fZ * fTmpB - 0.99258333397093026682 * fTmpA;
    pSH(278) = fTmpC * fC0;
    pSH(268) = -fTmpC * fS0;
    fTmpA = 2.0916500663351888699 * fZ * fTmpC - 0.99393320993236341952 * fTmpB;
    pSH(312) = fTmpA * fC0;
    pSH(302) = -fTmpA * fS0;
    fTmpB = 2.0811303848941723349 * fZ * fTmpA - 0.99497063031224519069 * fTmpC;
    pSH(348) = fTmpB * fC0;
    pSH(338) = -fTmpB * fS0;
    fTmpC = 2.0723520109148583406 * fZ * fTmpB - 0.99578192022804934259 * fTmpA;
    pSH(386) = fTmpC * fC0;
    pSH(376) = -fTmpC * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    fTmpA = 2.1809662438042814911 * fZ * fTmpC - 0.98588907503509976255 * fTmpB;
    pSH(247) = fTmpA * fC1;
    pSH(235) = fTmpA * fS1;
    fTmpB = 2.1563858652847824675 * fZ * fTmpA - 0.98872959240459255313 * fTmpC;
    pSH(279) = fTmpB * fC1;
    pSH(267) = fTmpB * fS1;
    fTmpC = 2.1366369348357590877 * fZ * fTmpB - 0.99084165280112553369 * fTmpA;
    pSH(313) = fTmpC * fC1;
    pSH(301) = fTmpC * fS1;
    fTmpA = 2.1205017749999120851 * fZ * fTmpC - 0.99244833805276922520 * fTmpB;
    pSH(349) = fTmpA * fC1;
    pSH(337) = fTmpA * fS1;
    fTmpB = 2.1071307505705477697 * fZ * fTmpA - 0.99369440545299007363 * fTmpC;
    pSH(387) = fTmpB * fC1;
    pSH(375) = fTmpB * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    fTmpC = 2.2600784378986817489 * fZ * fTmpB - 0.97926740294193723593 * fTmpA;
    pSH(248) = fTmpC * fC0;
    pSH(234) = -fTmpC * fS0;
    fTmpA = 2.2230674720995866294 * fZ * fTmpC - 0.98362403482177094968 * fTmpB;
    pSH(280) = fTmpA * fC0;
    pSH(266) = -fTmpA * fS0;
    fTmpB = 2.1937410968480305151 * fZ * fTmpA - 0.98680814882156560011 * fTmpC;
    pSH(314) = fTmpB * fC0;
    pSH(300) = -fTmpB * fS0;
    fTmpC = 2.1700439878239586254 * fZ * fTmpB - 0.98919785518075951600 * fTmpA;
    pSH(350) = fTmpC * fC0;
    pSH(336) = -fTmpC * fS0;
    fTmpA = 2.1505813167606566929 * fZ * fTmpC - 0.99103120896511485334 * fTmpB;
    pSH(388) = fTmpA * fC0;
    pSH(374) = -fTmpA * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    fTmpB = 2.3630173363047971159 * fZ * fTmpA - 0.97022618722766530118 * fTmpC;
    pSH(249) = fTmpB * fC1;
    pSH(233) = fTmpB * fS1;
    fTmpC = 2.3082731640774234847 * fZ * fTmpB - 0.97683293669229670981 * fTmpA;
    pSH(281) = fTmpC * fC1;
    pSH(265) = fTmpC * fS1;
    fTmpA = 2.2656860623955237928 * fZ * fTmpC - 0.98155023315928857429 * fTmpB;
    pSH(315) = fTmpA * fC1;
    pSH(299) = fTmpA * fS1;
    fTmpB = 2.2317637040621551360 * fZ * fTmpA - 0.98502777640009739967 * fTmpC;
    pSH(351) = fTmpB * fC1;
    pSH(335) = fTmpB * fS1;
    fTmpC = 2.2042200113840402628 * fZ * fTmpB - 0.98765832931686222325 * fTmpA;
    pSH(389) = fTmpC * fC1;
    pSH(373) = fTmpC * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    fTmpA = 2.4986107250941583108 * fZ * fTmpC - 0.95756141751469769713 * fTmpB;
    pSH(250) = fTmpA * fC0;
    pSH(232) = -fTmpA * fS0;
    fTmpB = 2.4177911997760033443 * fZ * fTmpA - 0.96765421499777267545 * fTmpC;
    pSH(282) = fTmpB * fC0;
    pSH(264) = -fTmpB * fS0;
    fTmpC = 2.3564559438666820500 * fZ * fTmpB - 0.97463169858712211832 * fTmpA;
    pSH(316) = fTmpC * fC0;
    pSH(298) = -fTmpC * fS0;
    fTmpA = 2.3085099321848032226 * fZ * fTmpC - 0.97965333839290678348 * fTmpB;
    pSH(352) = fTmpA * fC0;
    pSH(334) = -fTmpA * fS0;
    fTmpB = 2.2701478869385202614 * fZ * fTmpA - 0.98338233476432278865 * fTmpC;
    pSH(390) = fTmpB * fC0;
    pSH(372) = -fTmpB * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    fTmpC = 2.6817904466978772501 * fZ * fTmpB - 0.93903023262181382111 * fTmpA;
    pSH(251) = fTmpC * fC1;
    pSH(231) = fTmpC * fS1;
    fTmpA = 2.5607991541103546085 * fZ * fTmpC - 0.95488413617980451759 * fTmpB;
    pSH(283) = fTmpA * fC1;
    pSH(263) = fTmpA * fS1;
    fTmpB = 2.4720661623652209829 * fZ * fTmpA - 0.96534949193391143148 * fTmpC;
    pSH(317) = fTmpB * fC1;
    pSH(297) = fTmpB * fS1;
    fTmpC = 2.4044230077089180940 * fZ * fTmpB - 0.97263699666048446727 * fTmpA;
    pSH(353) = fTmpC * fC1;
    pSH(333) = fTmpC * fS1;
    fTmpA = 2.3513263559497452386 * fZ * fTmpC - 0.97791709213023767975 * fTmpB;
    pSH(391) = fTmpA * fC1;
    pSH(371) = fTmpA * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    fTmpB = 2.9401072717216916592 * fZ * fTmpA - 0.90994035683194807483 * fTmpC;
    pSH(252) = fTmpB * fC0;
    pSH(230) = -fTmpB * fS0;
    fTmpC = 2.7527763762750104249 * fZ * fTmpB - 0.93628433314374189412 * fTmpA;
    pSH(284) = fTmpC * fC0;
    pSH(262) = -fTmpC * fS0;
    fTmpA = 2.6220221204253788675 * fZ * fTmpC - 0.95250095250142875238 * fTmpB;
    pSH(318) = fTmpA * fC0;
    pSH(296) = -fTmpA * fS0;
    fTmpB = 2.5257296658248257996 * fZ * fTmpA - 0.96327549876469720969 * fTmpC;
    pSH(354) = fTmpB * fC0;
    pSH(332) = -fTmpB * fS0;
    fTmpC = 2.4520399670478456539 * fZ * fTmpB - 0.97082439194737994593 * fTmpA;
    pSH(392) = fTmpC * fC0;
    pSH(370) = -fTmpC * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    fTmpA = fZ * (53.794072123058085929 * fZ2 - 5.5649040127301468202);
    pSH(253) = fTmpA * fC1;
    pSH(229) = fTmpA * fS1;
    fTmpB = 3.0222389997200041804 * fZ * fTmpA - 0.90717582656041188324 * fTmpC;
    pSH(285) = fTmpB * fC1;
    pSH(261) = fTmpB * fS1;
    fTmpC = 2.8223247937435036367 * fZ * fTmpB - 0.93385228435109810060 * fTmpA;
    pSH(319) = fTmpC * fC1;
    pSH(295) = fTmpC * fS1;
    fTmpA = 2.6822461565718468645 * fZ * fTmpC - 0.95036764107299718112 * fTmpB;
    pSH(355) = fTmpA * fC1;
    pSH(331) = fTmpA * fS1;
    fTmpB = 2.5787147157554005437 * fZ * fTmpA - 0.96140121570767059495 * fTmpC;
    pSH(393) = fTmpB * fC1;
    pSH(369) = fTmpB * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    fTmpC = -17.608243388844820556 * fZ2 + 0.60718080651189036399;
    pSH(254) = fTmpC * fC0;
    pSH(228) = -fTmpC * fS0;
    fTmpA = fZ * (-60.380154942952015813 * fZ2 + 5.8432408009308402400);
    pSH(286) = fTmpA * fC0;
    pSH(260) = -fTmpA * fS0;
    fTmpB = 3.1024184114977141490 * fZ * fTmpA - 0.90473663963430495830 * fTmpC;
    pSH(320) = fTmpB * fC0;
    pSH(294) = -fTmpB * fS0;
    fTmpC = 2.8904737863674562919 * fZ * fTmpB - 0.93168406158731500387 * fTmpA;
    pSH(356) = fTmpC * fC0;
    pSH(330) = -fTmpC * fS0;
    fTmpA = 2.7414640249326636021 * fZ * fTmpC - 0.94844798034925006032 * fTmpB;
    pSH(394) = fTmpA * fC0;
    pSH(368) = -fTmpA * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.83052208306452400305;
    pSH(225) = fTmpA * fC1;
    pSH(197) = fTmpA * fS1;
    fTmpB = 4.6241512566300114788 * fZ;
    pSH(255) = fTmpB * fC1;
    pSH(227) = fTmpB * fS1;
    fTmpC = 19.093881509360250421 * fZ2 - 0.61593166159226614262;
    pSH(287) = fTmpC * fC1;
    pSH(259) = fTmpC * fS1;
    fTmpA = fZ * (67.288948373844056316 * fZ2 - 6.1171771248949142105);
    pSH(321) = fTmpA * fC1;
    pSH(293) = fTmpA * fS1;
    fTmpB = 3.1807526624998681277 * fZ * fTmpA - 0.90256893466271141000 * fTmpC;
    pSH(357) = fTmpB * fC1;
    pSH(329) = fTmpB * fS1;
    fTmpC = 2.9572714696920445985 * fZ * fTmpB - 0.92973952503676233437 * fTmpA;
    pSH(395) = fTmpC * fC1;
    pSH(367) = fTmpC * fS1;
    //! m = 15
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.84425065085737263613;
    pSH(256) = fTmpA * fC0;
    pSH(226) = -fTmpA * fS0;
    fTmpB = -4.8498507532306824834 * fZ;
    pSH(288) = fTmpB * fC0;
    pSH(258) = -fTmpB * fS0;
    fTmpC = -20.602948605549728459 * fZ2 + 0.62433177592574934724;
    pSH(322) = fTmpC * fC0;
    pSH(292) = -fTmpC * fS0;
    fTmpA = fZ * (-74.515507921532000473 * fZ2 + 6.3870435361313143263);
    pSH(358) = fTmpA * fC0;
    pSH(328) = -fTmpA * fS0;
    fTmpB = 3.2573446421352252765 * fZ * fTmpA - 0.90063003157873466786 * fTmpC;
    pSH(396) = fTmpB * fC0;
    pSH(366) = -fTmpB * fS0;
    //! m = 16
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.85734058883802509634;
    pSH(289) = fTmpA * fC1;
    pSH(257) = fTmpA * fS1;
    fTmpB = 5.0720953248553606107 * fZ;
    pSH(323) = fTmpB * fC1;
    pSH(291) = fTmpB * fS1;
    fTmpC = 22.134404124649146409 * fZ2 - 0.63241154641854704026;
    pSH(359) = fTmpC * fC1;
    pSH(327) = fTmpC * fS1;
    fTmpA = fZ * (82.055245832745453664 * fZ2 - 6.6531280404928746214);
    pSH(397) = fTmpA * fC1;
    pSH(365) = fTmpA * fS1;
    //! m = 17
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.86985717192062808222;
    pSH(324) = fTmpA * fC0;
    pSH(290) = -fTmpA * fS0;
    fTmpB = -5.2911346120699725413 * fZ;
    pSH(360) = fTmpB * fC0;
    pSH(326) = -fTmpB * fS0;
    fTmpC = -23.687309134978252702 * fZ2 + 0.64019754418860142438;
    pSH(398) = fTmpC * fC0;
    pSH(364) = -fTmpC * fS0;
    //! m = 18
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.88185576867832886131;
    pSH(361) = fTmpA * fC1;
    pSH(325) = fTmpA * fS1;
    fTmpB = 5.5071875102722439487 * fZ;
    pSH(399) = fTmpB * fC1;
    pSH(363) = fTmpB * fS1;
    //! m = 19
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpC = -0.89338378434994943301;
    pSH(400) = fTmpC * fC0;
    pSH(362) = -fTmpC * fS0;
  } break;
  case 20: {
    //! m = 0
    pSH(1) = 0.28209479177387814347;
    pSH(3) = 0.48860251190291992159 * fZ;
    pSH(7) = 0.94617469575756001809 * fZ2 - 0.31539156525252000603;
    pSH(13) = fZ * (1.8658816629505769571 * fZ2 - 1.1195289977703461742);
    pSH(21) = 1.9843134832984429429 * fZ * pSH(13) - 1.0062305898749053634 * pSH(7);
    pSH(31) = 1.9899748742132399095 * fZ * pSH(21) - 1.0028530728448139498 * pSH(13);
    pSH(43) = 1.9930434571835663369 * fZ * pSH(31) - 1.0015420209622192481 * pSH(21);
    pSH(57) = 1.9948914348241344528 * fZ * pSH(43) - 1.0009272139219581055 * pSH(31);
    pSH(73) = 1.9960899278339139999 * fZ * pSH(57) - 1.0006007810695147948 * pSH(43);
    pSH(91) = 1.9969111950679364953 * fZ * pSH(73) - 1.0004114379931337590 * pSH(57);
    pSH(111) = 1.9974984355438178916 * fZ * pSH(91) - 1.0002940744071803443 * pSH(73);
    pSH(133) = 1.9979328159850827788 * fZ * pSH(111) - 1.0002174622185106380 * pSH(91);
    pSH(157) = 1.9982631347136331423 * fZ * pSH(133) - 1.0001653302482984141 * pSH(111);
    pSH(183) = 1.9985201625794738002 * fZ * pSH(157) - 1.0001286256356210525 * pSH(133);
    pSH(211) = 1.9987240828047460812 * fZ * pSH(183) - 1.0001020356106936058 * pSH(157);
    pSH(241) = 1.9988885800753266487 * fZ * pSH(211) - 1.0000823011400101477 * pSH(183);
    pSH(273) = 1.9990231989649344737 * fZ * pSH(241) - 1.0000673468701305763 * pSH(211);
    pSH(307) = 1.9991347609372268760 * fZ * pSH(273) - 1.0000558082429209263 * pSH(241);
    pSH(343) = 1.9992282461607312886 * fZ * pSH(307) - 1.0000467628422711159 * pSH(273);
    pSH(381) = 1.9993073592865872621 * fZ * pSH(343) - 1.0000395718327849261 * pSH(307);
    pSH(421) = 1.9993749023132204957 * fZ * pSH(381) - 1.0000337832131310391 * pSH(343);
    //! m = 1
    fC0 = fX;
    fS0 = fY;
    fTmpA = -0.48860251190291992159;
    pSH(4) = fTmpA * fC0;
    pSH(2) = -fTmpA * fS0;
    fTmpB = -1.0925484305920790705 * fZ;
    pSH(8) = fTmpB * fC0;
    pSH(6) = -fTmpB * fS0;
    fTmpC = -2.2852289973223286808 * fZ2 + 0.45704579946446573616;
    pSH(14) = fTmpC * fC0;
    pSH(12) = -fTmpC * fS0;
    fTmpA = fZ * (-4.6833258049010241757 * fZ2 + 2.0071396306718675039);
    pSH(22) = fTmpA * fC0;
    pSH(20) = -fTmpA * fS0;
    fTmpB = 2.0310096011589900901 * fZ * fTmpA - 0.99103120896511485334 * fTmpC;
    pSH(32) = fTmpB * fC0;
    pSH(30) = -fTmpB * fS0;
    fTmpC = 2.0213149892370277761 * fZ * fTmpB - 0.99522670305623857702 * fTmpA;
    pSH(44) = fTmpC * fC0;
    pSH(42) = -fTmpC * fS0;
    fTmpA = 2.0155644370746374131 * fZ * fTmpC - 0.99715504402183205232 * fTmpB;
    pSH(58) = fTmpA * fC0;
    pSH(56) = -fTmpA * fS0;
    fTmpB = 2.0118695404073912315 * fZ * fTmpA - 0.99816681789017427595 * fTmpC;
    pSH(74) = fTmpB * fC0;
    pSH(72) = -fTmpB * fS0;
    fTmpC = 2.0093531297410119494 * fZ * fTmpB - 0.99874921777190894579 * fTmpA;
    pSH(92) = fTmpC * fC0;
    pSH(90) = -fTmpC * fS0;
    fTmpA = 2.0075614636426527858 * fZ * fTmpC - 0.99910833687128449455 * fTmpB;
    pSH(112) = fTmpA * fC0;
    pSH(110) = -fTmpA * fS0;
    fTmpB = 2.0062402647738879433 * fZ * fTmpA - 0.99934188870792151413 * fTmpC;
    pSH(134) = fTmpB * fC0;
    pSH(132) = -fTmpB * fS0;
    fTmpC = 2.0052378963551982949 * fZ * fTmpB - 0.99950037468777319163 * fTmpA;
    pSH(158) = fTmpC * fC0;
    pSH(156) = -fTmpC * fS0;
    fTmpA = 2.0044593143431828851 * fZ * fTmpC - 0.99961172586383361697 * fTmpB;
    pSH(184) = fTmpA * fC0;
    pSH(182) = -fTmpA * fS0;
    fTmpB = 2.0038424627162224562 * fZ * fTmpA - 0.99969226034045866500 * fTmpC;
    pSH(212) = fTmpB * fC0;
    pSH(210) = -fTmpB * fS0;
    fTmpC = 2.0033454163331038362 * fZ * fTmpB - 0.99975195336341716696 * fTmpA;
    pSH(242) = fTmpC * fC0;
    pSH(240) = -fTmpC * fS0;
    fTmpA = 2.0029390170153341294 * fZ * fTmpC - 0.99979713966725040627 * fTmpB;
    pSH(274) = fTmpA * fC0;
    pSH(272) = -fTmpA * fS0;
    fTmpB = 2.0026024734496526300 * fZ * fTmpA - 0.99983197513113354916 * fTmpC;
    pSH(308) = fTmpB * fC0;
    pSH(306) = -fTmpB * fS0;
    fTmpC = 2.0023206350873464509 * fZ * fTmpB - 0.99985926394976398460 * fTmpA;
    pSH(344) = fTmpC * fC0;
    pSH(342) = -fTmpC * fS0;
    fTmpA = 2.0020822493926999835 * fZ * fTmpC - 0.99988094529394086354 * fTmpB;
    pSH(382) = fTmpA * fC0;
    pSH(380) = -fTmpA * fS0;
    fTmpB = 2.0018788167600158716 * fZ * fTmpA - 0.99989838947288713042 * fTmpC;
    pSH(422) = fTmpB * fC0;
    pSH(420) = -fTmpB * fS0;
    //! m = 2
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.54627421529603953527;
    pSH(9) = fTmpA * fC1;
    pSH(5) = fTmpA * fS1;
    fTmpB = 1.4453057213202770277 * fZ;
    pSH(15) = fTmpB * fC1;
    pSH(11) = fTmpB * fS1;
    fTmpC = 3.3116114351514600633 * fZ2 - 0.47308734787878000905;
    pSH(23) = fTmpC * fC1;
    pSH(19) = fTmpC * fS1;
    fTmpA = fZ * (7.1903051774599856325 * fZ2 - 2.3967683924866618775);
    pSH(33) = fTmpA * fC1;
    pSH(29) = fTmpA * fS1;
    fTmpB = 2.1139418156609703623 * fZ * fTmpA - 0.97361012046232688422 * fTmpC;
    pSH(45) = fTmpB * fC1;
    pSH(41) = fTmpB * fS1;
    fTmpC = 2.0816659994661327353 * fZ * fTmpB - 0.98473192783466186187 * fTmpA;
    pSH(59) = fTmpC * fC1;
    pSH(55) = fTmpC * fS1;
    fTmpA = 2.0615528128088302749 * fZ * fTmpC - 0.99033793766028713580 * fTmpB;
    pSH(75) = fTmpA * fC1;
    pSH(71) = fTmpA * fS1;
    fTmpB = 2.0481223583578191106 * fZ * fTmpA - 0.99348527267040401407 * fTmpC;
    pSH(93) = fTmpB * fC1;
    pSH(89) = fTmpB * fS1;
    fTmpC = 2.0386883037875113095 * fZ * fTmpB - 0.99539380324041186222 * fTmpA;
    pSH(113) = fTmpC * fC1;
    pSH(109) = fTmpC * fS1;
    fTmpA = 2.0317984959648750082 * fZ * fTmpC - 0.99662047022596034223 * fTmpB;
    pSH(135) = fTmpA * fC1;
    pSH(131) = fTmpA * fS1;
    fTmpB = 2.0266087084444439303 * fZ * fTmpA - 0.99744571741206722642 * fTmpC;
    pSH(159) = fTmpB * fC1;
    pSH(155) = fTmpB * fS1;
    fTmpC = 2.0225995873897262587 * fZ * fTmpB - 0.99802175869569072522 * fTmpA;
    pSH(185) = fTmpC * fC1;
    pSH(181) = fTmpC * fS1;
    fTmpA = 2.0194368026754390117 * fZ * fTmpC - 0.99843627738579290867 * fTmpB;
    pSH(213) = fTmpA * fC1;
    pSH(209) = fTmpA * fS1;
    fTmpB = 2.0168969490698876097 * fZ * fTmpA - 0.99874229606879180776 * fTmpC;
    pSH(243) = fTmpB * fC1;
    pSH(239) = fTmpB * fS1;
    fTmpC = 2.0148259998133361203 * fZ * fTmpB - 0.99897320026315349004 * fTmpA;
    pSH(275) = fTmpC * fC1;
    pSH(271) = fTmpC * fS1;
    fTmpA = 2.0131148946216081338 * fZ * fTmpC - 0.99915074294659364525 * fTmpB;
    pSH(309) = fTmpA * fC1;
    pSH(305) = fTmpA * fS1;
    fTmpB = 2.0116846174288851484 * fZ * fTmpA - 0.99928952033659667242 * fTmpC;
    pSH(345) = fTmpB * fC1;
    pSH(341) = fTmpB * fS1;
    fTmpC = 2.0104767610501468006 * fZ * fTmpB - 0.99939957965166423681 * fTmpA;
    pSH(383) = fTmpC * fC1;
    pSH(379) = fTmpC * fS1;
    fTmpA = 2.0094473837049796908 * fZ * fTmpC - 0.99948799341275179535 * fTmpB;
    pSH(423) = fTmpA * fC1;
    pSH(419) = fTmpA * fS1;
    //! m = 3
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.59004358992664351035;
    pSH(16) = fTmpA * fC0;
    pSH(10) = -fTmpA * fS0;
    fTmpB = -1.7701307697799305310 * fZ;
    pSH(24) = fTmpB * fC0;
    pSH(18) = -fTmpB * fS0;
    fTmpC = -4.4031446949172534892 * fZ2 + 0.48923829943525038768;
    pSH(34) = fTmpC * fC0;
    pSH(28) = -fTmpC * fS0;
    fTmpA = fZ * (-10.133257854664158491 * fZ2 + 2.7636157785447704974);
    pSH(46) = fTmpA * fC0;
    pSH(40) = -fTmpA * fS0;
    fTmpB = 2.2079402165819617137 * fZ * fTmpA - 0.95940322360024695434 * fTmpC;
    pSH(60) = fTmpB * fC0;
    pSH(54) = -fTmpB * fS0;
    fTmpC = 2.1532216876958202242 * fZ * fTmpB - 0.97521738656001772954 * fTmpA;
    pSH(76) = fTmpC * fC0;
    pSH(70) = -fTmpC * fS0;
    fTmpA = 2.1180441711898057371 * fZ * fTmpC - 0.98366284497920962827 * fTmpB;
    pSH(94) = fTmpA * fC0;
    pSH(88) = -fTmpA * fS0;
    fTmpB = 2.0939473213563383757 * fZ * fTmpA - 0.98862306548596150408 * fTmpC;
    pSH(114) = fTmpB * fC0;
    pSH(108) = -fTmpB * fS0;
    fTmpC = 2.0766559657295187131 * fZ * fTmpB - 0.99174222032690902698 * fTmpA;
    pSH(136) = fTmpC * fC0;
    pSH(130) = -fTmpC * fS0;
    fTmpA = 2.0637972912229677746 * fZ * fTmpC - 0.99380798999990653174 * fTmpB;
    pSH(160) = fTmpA * fC0;
    pSH(154) = -fTmpA * fS0;
    fTmpB = 2.0539595906443729255 * fZ * fTmpA - 0.99523320404555565084 * fTmpC;
    pSH(186) = fTmpB * fC0;
    pSH(180) = -fTmpB * fS0;
    fTmpC = 2.0462565272714634688 * fZ * fTmpB - 0.99624965193668058779 * fTmpA;
    pSH(214) = fTmpC * fC0;
    pSH(208) = -fTmpC * fS0;
    fTmpA = 2.0401071141087266541 * fZ * fTmpC - 0.99699479851094886098 * fTmpB;
    pSH(244) = fTmpA * fC0;
    pSH(238) = -fTmpA * fS0;
    fTmpB = 2.0351168037383750115 * fZ * fTmpA - 0.99755389786357772286 * fTmpC;
    pSH(276) = fTmpB * fC0;
    pSH(270) = -fTmpB * fS0;
    fTmpC = 2.0310096011589900901 * fZ * fTmpB - 0.99798183447169211035 * fTmpA;
    pSH(310) = fTmpC * fC0;
    pSH(304) = -fTmpC * fS0;
    fTmpA = 2.0275875100994065630 * fZ * fTmpC - 0.99831507883683527632 * fTmpB;
    pSH(346) = fTmpA * fC0;
    pSH(340) = -fTmpA * fS0;
    fTmpB = 2.0247053657709850061 * fZ * fTmpA - 0.99857853517341885086 * fTmpC;
    pSH(384) = fTmpB * fC0;
    pSH(378) = -fTmpB * fS0;
    fTmpC = 2.0222546987202585513 * fZ * fTmpB - 0.99878961794038943127 * fTmpA;
    pSH(424) = fTmpC * fC0;
    pSH(418) = -fTmpC * fS0;
    //! m = 4
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.62583573544917613459;
    pSH(25) = fTmpA * fC1;
    pSH(17) = fTmpA * fS1;
    fTmpB = 2.0756623148810412790 * fZ;
    pSH(35) = fTmpB * fC1;
    pSH(27) = fTmpB * fS1;
    fTmpC = 5.5502139080159657518 * fZ2 - 0.50456490072872415925;
    pSH(47) = fTmpC * fC1;
    pSH(39) = fTmpC * fS1;
    fTmpA = fZ * (13.491805046726768313 * fZ2 - 3.1134934723215619185);
    pSH(61) = fTmpA * fC1;
    pSH(53) = fTmpA * fS1;
    fTmpB = 2.3048861143232218275 * fZ * fTmpA - 0.94817638735546538523 * fTmpC;
    pSH(77) = fTmpB * fC1;
    pSH(69) = fTmpB * fS1;
    fTmpC = 2.2291771507062351977 * fZ * fTmpB - 0.96715283972318221417 * fTmpA;
    pSH(95) = fTmpC * fC1;
    pSH(87) = fTmpC * fS1;
    fTmpA = 2.1794494717703367761 * fZ * fTmpC - 0.97769236109380361190 * fTmpB;
    pSH(115) = fTmpA * fC1;
    pSH(107) = fTmpA * fS1;
    fTmpB = 2.1447610589527216610 * fZ * fTmpA - 0.98408386463328365425 * fTmpC;
    pSH(137) = fTmpB * fC1;
    pSH(129) = fTmpB * fS1;
    fTmpC = 2.1194781197266462935 * fZ * fTmpB - 0.98821176880261854125 * fTmpA;
    pSH(161) = fTmpC * fC1;
    pSH(153) = fTmpC * fS1;
    fTmpA = 2.1004201260420147053 * fZ * fTmpC - 0.99100816681840077731 * fTmpB;
    pSH(187) = fTmpA * fC1;
    pSH(179) = fTmpA * fS1;
    fTmpB = 2.0856653614614210205 * fZ * fTmpA - 0.99297532698451274258 * fTmpC;
    pSH(215) = fTmpB * fC1;
    pSH(207) = fTmpB * fS1;
    fTmpC = 2.0739902137422357166 * fZ * fTmpB - 0.99440219512922986155 * fTmpA;
    pSH(245) = fTmpC * fC1;
    pSH(237) = fTmpC * fS1;
    fTmpA = 2.0645822822062578187 * fZ * fTmpC - 0.99546384960081245811 * fTmpB;
    pSH(277) = fTmpA * fC1;
    pSH(269) = fTmpA * fS1;
    fTmpB = 2.0568833780186057912 * fZ * fTmpA - 0.99627096277343579159 * fTmpC;
    pSH(311) = fTmpB * fC1;
    pSH(303) = fTmpB * fS1;
    fTmpC = 2.0504988306618110689 * fZ * fTmpB - 0.99689600906642312802 * fTmpA;
    pSH(347) = fTmpC * fC1;
    pSH(339) = fTmpC * fS1;
    fTmpA = 2.0451427078940417827 * fZ * fTmpC - 0.99738789279580297796 * fTmpB;
    pSH(385) = fTmpA * fC1;
    pSH(377) = fTmpA * fS1;
    fTmpB = 2.0406034646643134620 * fZ * fTmpA - 0.99778047604589777240 * fTmpC;
    pSH(425) = fTmpB * fC1;
    pSH(417) = fTmpB * fS1;
    //! m = 5
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.65638205684017010281;
    pSH(36) = fTmpA * fC0;
    pSH(26) = -fTmpA * fS0;
    fTmpB = -2.3666191622317520320 * fZ;
    pSH(48) = fTmpB * fC0;
    pSH(38) = -fTmpB * fS0;
    fTmpC = -6.7459025233633841567 * fZ2 + 0.51891557872026031975;
    pSH(62) = fTmpC * fC0;
    pSH(52) = -fTmpC * fS0;
    fTmpA = fZ * (-17.249553110490540088 * fZ2 + 3.4499106220981080175);
    pSH(78) = fTmpA * fC0;
    pSH(68) = -fTmpA * fS0;
    fTmpB = 2.4016363469220611496 * fZ * fTmpA - 0.93922460420437088487 * fTmpC;
    pSH(96) = fTmpB * fC0;
    pSH(86) = -fTmpB * fS0;
    fTmpC = 2.3065125189341591779 * fZ * fTmpB - 0.96039207679804948932 * fTmpA;
    pSH(116) = fTmpC * fC0;
    pSH(106) = -fTmpC * fS0;
    fTmpA = 2.2430448056157950943 * fZ * fTmpC - 0.97248325651937386751 * fTmpB;
    pSH(138) = fTmpA * fC0;
    pSH(128) = -fTmpA * fS0;
    fTmpB = 2.1981657747106435415 * fZ * fTmpA - 0.97999191510005049931 * fTmpC;
    pSH(162) = fTmpB * fC0;
    pSH(152) = -fTmpB * fS0;
    fTmpC = 2.1650635094610966169 * fZ * fTmpB - 0.98494096049061433733 * fTmpA;
    pSH(188) = fTmpC * fC0;
    pSH(178) = -fTmpC * fS0;
    fTmpA = 2.1398475105532759975 * fZ * fTmpC - 0.98835322899414756490 * fTmpB;
    pSH(216) = fTmpA * fC0;
    pSH(206) = -fTmpA * fS0;
    fTmpB = 2.1201415047114190203 * fZ * fTmpA - 0.99079092984678996135 * fTmpC;
    pSH(246) = fTmpB * fC0;
    pSH(236) = -fTmpB * fS0;
    fTmpC = 2.1044171232366050512 * fZ * fTmpB - 0.99258333397093026682 * fTmpA;
    pSH(278) = fTmpC * fC0;
    pSH(268) = -fTmpC * fS0;
    fTmpA = 2.0916500663351888699 * fZ * fTmpC - 0.99393320993236341952 * fTmpB;
    pSH(312) = fTmpA * fC0;
    pSH(302) = -fTmpA * fS0;
    fTmpB = 2.0811303848941723349 * fZ * fTmpA - 0.99497063031224519069 * fTmpC;
    pSH(348) = fTmpB * fC0;
    pSH(338) = -fTmpB * fS0;
    fTmpC = 2.0723520109148583406 * fZ * fTmpB - 0.99578192022804934259 * fTmpA;
    pSH(386) = fTmpC * fC0;
    pSH(376) = -fTmpC * fS0;
    fTmpA = 2.0649455198624490587 * fZ * fTmpC - 0.99642604585832904896 * fTmpB;
    pSH(426) = fTmpA * fC0;
    pSH(416) = -fTmpA * fS0;
    //! m = 6
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.68318410519191432198;
    pSH(49) = fTmpA * fC1;
    pSH(37) = fTmpA * fS1;
    fTmpB = 2.6459606618019002220 * fZ;
    pSH(63) = fTmpB * fC1;
    pSH(51) = fTmpB * fS1;
    fTmpC = 7.9849914908931386147 * fZ2 - 0.53233276605954257431;
    pSH(79) = fTmpC * fC1;
    pSH(67) = fTmpC * fS1;
    fTmpA = fZ * (21.392890190908636255 * fZ2 - 3.7752159160427005155);
    pSH(97) = fTmpA * fC1;
    pSH(85) = fTmpA * fS1;
    fTmpB = 2.4968730444297723645 * fZ * fTmpA - 0.93196897827695329104 * fTmpC;
    pSH(117) = fTmpB * fC1;
    pSH(105) = fTmpB * fS1;
    fTmpC = 2.3837686425440851889 * fZ * fTmpB - 0.95470158078801415952 * fTmpA;
    pSH(139) = fTmpC * fC1;
    pSH(127) = fTmpC * fS1;
    fTmpA = 2.3073955174772430146 * fZ * fTmpC - 0.96796118394051333064 * fTmpB;
    pSH(163) = fTmpA * fC1;
    pSH(151) = fTmpA * fS1;
    fTmpB = 2.2528177844479149153 * fZ * fTmpA - 0.97634660697919710715 * fTmpC;
    pSH(189) = fTmpB * fC1;
    pSH(177) = fTmpB * fS1;
    fTmpC = 2.2121821805628938752 * fZ * fTmpB - 0.98196232106939826309 * fTmpA;
    pSH(217) = fTmpC * fC1;
    pSH(205) = fTmpC * fS1;
    fTmpA = 2.1809662438042814911 * fZ * fTmpC - 0.98588907503509976255 * fTmpB;
    pSH(247) = fTmpA * fC1;
    pSH(235) = fTmpA * fS1;
    fTmpB = 2.1563858652847824675 * fZ * fTmpA - 0.98872959240459255313 * fTmpC;
    pSH(279) = fTmpB * fC1;
    pSH(267) = fTmpB * fS1;
    fTmpC = 2.1366369348357590877 * fZ * fTmpB - 0.99084165280112553369 * fTmpA;
    pSH(313) = fTmpC * fC1;
    pSH(301) = fTmpC * fS1;
    fTmpA = 2.1205017749999120851 * fZ * fTmpC - 0.99244833805276922520 * fTmpB;
    pSH(349) = fTmpA * fC1;
    pSH(337) = fTmpA * fS1;
    fTmpB = 2.1071307505705477697 * fZ * fTmpA - 0.99369440545299007363 * fTmpC;
    pSH(387) = fTmpB * fC1;
    pSH(375) = fTmpB * fS1;
    fTmpC = 2.0959143930173156901 * fZ * fTmpB - 0.99467695227256541789 * fTmpA;
    pSH(427) = fTmpC * fC1;
    pSH(415) = fTmpC * fS1;
    //! m = 7
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.70716273252459617823;
    pSH(64) = fTmpA * fC0;
    pSH(50) = -fTmpA * fS0;
    fTmpB = -2.9157066406993194755 * fZ;
    pSH(80) = fTmpB * fC0;
    pSH(66) = -fTmpB * fS0;
    fTmpC = -9.2633931828489042401 * fZ2 + 0.54490548134405319060;
    pSH(98) = fTmpC * fC0;
    pSH(84) = -fTmpC * fS0;
    fTmpA = fZ * (-25.910241313366302025 * fZ2 + 4.0910907336894161093);
    pSH(118) = fTmpA * fC0;
    pSH(104) = -fTmpA * fS0;
    fTmpB = 2.5900450446533421889 * fZ * fTmpA - 0.92598927658525138721 * fTmpC;
    pSH(140) = fTmpB * fC0;
    pSH(126) = -fTmpB * fS0;
    fTmpC = 2.4602096615832091356 * fZ * fTmpB - 0.94987138029195529652 * fTmpA;
    pSH(164) = fTmpC * fC0;
    pSH(150) = -fTmpC * fS0;
    fTmpA = 2.3717082451262844990 * fZ * fTmpC - 0.96402688037572713823 * fTmpB;
    pSH(190) = fTmpA * fC0;
    pSH(176) = -fTmpA * fS0;
    fTmpB = 2.3079277744862160134 * fZ * fTmpA - 0.97310779233865149521 * fTmpC;
    pSH(218) = fTmpB * fC0;
    pSH(204) = -fTmpB * fS0;
    fTmpC = 2.2600784378986817489 * fZ * fTmpB - 0.97926740294193723593 * fTmpA;
    pSH(248) = fTmpC * fC0;
    pSH(234) = -fTmpC * fS0;
    fTmpA = 2.2230674720995866294 * fZ * fTmpC - 0.98362403482177094968 * fTmpB;
    pSH(280) = fTmpA * fC0;
    pSH(266) = -fTmpA * fS0;
    fTmpB = 2.1937410968480305151 * fZ * fTmpA - 0.98680814882156560011 * fTmpC;
    pSH(314) = fTmpB * fC0;
    pSH(300) = -fTmpB * fS0;
    fTmpC = 2.1700439878239586254 * fZ * fTmpB - 0.98919785518075951600 * fTmpA;
    pSH(350) = fTmpC * fC0;
    pSH(336) = -fTmpC * fS0;
    fTmpA = 2.1505813167606566929 * fZ * fTmpC - 0.99103120896511485334 * fTmpB;
    pSH(388) = fTmpA * fC0;
    pSH(374) = -fTmpA * fS0;
    fTmpB = 2.1343747458109495622 * fZ * fTmpA - 0.99246409757984947724 * fTmpC;
    pSH(428) = fTmpB * fC0;
    pSH(414) = -fTmpB * fS0;
    //! m = 8
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.72892666017482986887;
    pSH(81) = fTmpA * fC1;
    pSH(65) = fTmpA * fS1;
    fTmpB = 3.1773176489546974773 * fZ;
    pSH(99) = fTmpB * fC1;
    pSH(83) = fTmpB * fS1;
    fTmpC = 10.577811721687949636 * fZ2 - 0.55672693272041840189;
    pSH(119) = fTmpC * fC1;
    pSH(103) = fTmpC * fS1;
    fTmpA = fZ * (30.791579703357485663 * fZ2 - 4.3987971004796408090);
    pSH(141) = fTmpA * fC1;
    pSH(125) = fTmpA * fS1;
    fTmpB = 2.6809513236909020762 * fZ * fTmpA - 0.92098549701625905369 * fTmpC;
    pSH(165) = fTmpB * fC1;
    pSH(149) = fTmpB * fS1;
    fTmpC = 2.5354627641855497325 * fZ * fTmpB - 0.94573248748692077368 * fTmpA;
    pSH(191) = fTmpC * fC1;
    pSH(175) = fTmpC * fS1;
    fTmpA = 2.4355324226579661401 * fZ * fTmpC - 0.96058694178469484714 * fTmpB;
    pSH(219) = fTmpA * fC1;
    pSH(203) = fTmpA * fS1;
    fTmpB = 2.3630173363047971159 * fZ * fTmpA - 0.97022618722766530118 * fTmpC;
    pSH(249) = fTmpB * fC1;
    pSH(233) = fTmpB * fS1;
    fTmpC = 2.3082731640774234847 * fZ * fTmpB - 0.97683293669229670981 * fTmpA;
    pSH(281) = fTmpC * fC1;
    pSH(265) = fTmpC * fS1;
    fTmpA = 2.2656860623955237928 * fZ * fTmpC - 0.98155023315928857429 * fTmpB;
    pSH(315) = fTmpA * fC1;
    pSH(299) = fTmpA * fS1;
    fTmpB = 2.2317637040621551360 * fZ * fTmpA - 0.98502777640009739967 * fTmpC;
    pSH(351) = fTmpB * fC1;
    pSH(335) = fTmpB * fS1;
    fTmpC = 2.2042200113840402628 * fZ * fTmpB - 0.98765832931686222325 * fTmpA;
    pSH(389) = fTmpC * fC1;
    pSH(373) = fTmpC * fS1;
    fTmpA = 2.1814968648679216996 * fZ * fTmpC - 0.98969107149070360218 * fTmpB;
    pSH(429) = fTmpA * fC1;
    pSH(413) = fTmpA * fS1;
    //! m = 9
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.74890095185318829655;
    pSH(100) = fTmpA * fC0;
    pSH(82) = -fTmpA * fS0;
    fTmpB = -3.4318952998917144349 * fZ;
    pSH(120) = fTmpB * fC0;
    pSH(102) = -fTmpB * fS0;
    fTmpC = -11.925527539452185581 * fZ2 + 0.56788226378343740862;
    pSH(142) = fTmpC * fC0;
    pSH(124) = -fTmpC * fS0;
    fTmpA = fZ * (-36.028090689310769890 * fZ2 + 4.6993161768666221596);
    pSH(166) = fTmpA * fC0;
    pSH(148) = -fTmpA * fS0;
    fTmpB = 2.7695585470349864865 * fZ * fTmpA - 0.91674152287482094273 * fTmpC;
    pSH(192) = fTmpB * fC0;
    pSH(174) = -fTmpB * fS0;
    fTmpC = 2.6093477445855914594 * fZ * fTmpB - 0.94215294613615865838 * fTmpA;
    pSH(220) = fTmpC * fC0;
    pSH(202) = -fTmpC * fS0;
    fTmpA = 2.4986107250941583108 * fZ * fTmpC - 0.95756141751469769713 * fTmpB;
    pSH(250) = fTmpA * fC0;
    pSH(232) = -fTmpA * fS0;
    fTmpB = 2.4177911997760033443 * fZ * fTmpA - 0.96765421499777267545 * fTmpC;
    pSH(282) = fTmpB * fC0;
    pSH(264) = -fTmpB * fS0;
    fTmpC = 2.3564559438666820500 * fZ * fTmpB - 0.97463169858712211832 * fTmpA;
    pSH(316) = fTmpC * fC0;
    pSH(298) = -fTmpC * fS0;
    fTmpA = 2.3085099321848032226 * fZ * fTmpC - 0.97965333839290678348 * fTmpB;
    pSH(352) = fTmpA * fC0;
    pSH(334) = -fTmpA * fS0;
    fTmpB = 2.2701478869385202614 * fZ * fTmpA - 0.98338233476432278865 * fTmpC;
    pSH(390) = fTmpB * fC0;
    pSH(372) = -fTmpB * fS0;
    fTmpC = 2.2388700687965297948 * fZ * fTmpB - 0.98622212309517369923 * fTmpA;
    pSH(430) = fTmpC * fC0;
    pSH(412) = -fTmpC * fS0;
    //! m = 10
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.76739511822199001256;
    pSH(121) = fTmpA * fC1;
    pSH(101) = fTmpA * fS1;
    fTmpB = 3.6802976988053108636 * fZ;
    pSH(143) = fTmpB * fC1;
    pSH(123) = fTmpB * fS1;
    fTmpC = 13.304254200257634746 * fZ2 - 0.57844583479381020637;
    pSH(167) = fTmpC * fC1;
    pSH(147) = fTmpC * fS1;
    fTmpA = fZ * (41.611931535496447639 * fZ2 - 4.9934317842595737167);
    pSH(193) = fTmpA * fC1;
    pSH(173) = fTmpA * fS1;
    fTmpB = 2.8559149146989656071 * fZ * fTmpA - 0.91309911838748371447 * fTmpC;
    pSH(221) = fTmpB * fC1;
    pSH(201) = fTmpB * fS1;
    fTmpC = 2.6817904466978772501 * fZ * fTmpB - 0.93903023262181382111 * fTmpA;
    pSH(251) = fTmpC * fC1;
    pSH(231) = fTmpC * fS1;
    fTmpA = 2.5607991541103546085 * fZ * fTmpC - 0.95488413617980451759 * fTmpB;
    pSH(283) = fTmpA * fC1;
    pSH(263) = fTmpA * fS1;
    fTmpB = 2.4720661623652209829 * fZ * fTmpA - 0.96534949193391143148 * fTmpC;
    pSH(317) = fTmpB * fC1;
    pSH(297) = fTmpB * fS1;
    fTmpC = 2.4044230077089180940 * fZ * fTmpB - 0.97263699666048446727 * fTmpA;
    pSH(353) = fTmpC * fC1;
    pSH(333) = fTmpC * fS1;
    fTmpA = 2.3513263559497452386 * fZ * fTmpC - 0.97791709213023767975 * fTmpB;
    pSH(391) = fTmpA * fC1;
    pSH(371) = fTmpA * fS1;
    fTmpB = 2.3086792761230391397 * fZ * fTmpA - 0.98186254336034942466 * fTmpC;
    pSH(431) = fTmpB * fC1;
    pSH(411) = fTmpB * fS1;
    //! m = 11
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.78464210578719688375;
    pSH(144) = fTmpA * fC0;
    pSH(122) = -fTmpA * fS0;
    fTmpB = -3.9232105289359851156 * fZ;
    pSH(168) = fTmpB * fC0;
    pSH(146) = -fTmpB * fS0;
    fTmpC = -14.712039483509941570 * fZ2 + 0.58848157934039766281;
    pSH(194) = fTmpC * fC0;
    pSH(172) = -fTmpC * fS0;
    fTmpA = fZ * (-47.536054360662613678 * fZ2 + 5.2817838178514015198);
    pSH(222) = fTmpA * fC0;
    pSH(200) = -fTmpA * fS0;
    fTmpB = 2.9401072717216916592 * fZ * fTmpA - 0.90994035683194807483 * fTmpC;
    pSH(252) = fTmpB * fC0;
    pSH(230) = -fTmpB * fS0;
    fTmpC = 2.7527763762750104249 * fZ * fTmpB - 0.93628433314374189412 * fTmpA;
    pSH(284) = fTmpC * fC0;
    pSH(262) = -fTmpC * fS0;
    fTmpA = 2.6220221204253788675 * fZ * fTmpC - 0.95250095250142875238 * fTmpB;
    pSH(318) = fTmpA * fC0;
    pSH(296) = -fTmpA * fS0;
    fTmpB = 2.5257296658248257996 * fZ * fTmpA - 0.96327549876469720969 * fTmpC;
    pSH(354) = fTmpB * fC0;
    pSH(332) = -fTmpB * fS0;
    fTmpC = 2.4520399670478456539 * fZ * fTmpB - 0.97082439194737994593 * fTmpA;
    pSH(392) = fTmpC * fC0;
    pSH(370) = -fTmpC * fS0;
    fTmpA = 2.3939888879647968682 * fZ * fTmpC - 0.97632539442130713904 * fTmpB;
    pSH(432) = fTmpA * fC0;
    pSH(410) = -fTmpA * fS0;
    //! m = 12
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.80082199578397171410;
    pSH(169) = fTmpA * fC1;
    pSH(145) = fTmpA * fS1;
    fTmpB = 4.1611931535496447639 * fZ;
    pSH(195) = fTmpB * fC1;
    pSH(171) = fTmpB * fS1;
    fTmpC = 16.147194793928202586 * fZ2 - 0.59804425162697046613;
    pSH(223) = fTmpC * fC1;
    pSH(199) = fTmpC * fS1;
    fTmpA = fZ * (53.794072123058085929 * fZ2 - 5.5649040127301468202);
    pSH(253) = fTmpA * fC1;
    pSH(229) = fTmpA * fS1;
    fTmpB = 3.0222389997200041804 * fZ * fTmpA - 0.90717582656041188324 * fTmpC;
    pSH(285) = fTmpB * fC1;
    pSH(261) = fTmpB * fS1;
    fTmpC = 2.8223247937435036367 * fZ * fTmpB - 0.93385228435109810060 * fTmpA;
    pSH(319) = fTmpC * fC1;
    pSH(295) = fTmpC * fS1;
    fTmpA = 2.6822461565718468645 * fZ * fTmpC - 0.95036764107299718112 * fTmpB;
    pSH(355) = fTmpA * fC1;
    pSH(331) = fTmpA * fS1;
    fTmpB = 2.5787147157554005437 * fZ * fTmpA - 0.96140121570767059495 * fTmpC;
    pSH(393) = fTmpB * fC1;
    pSH(369) = fTmpB * fS1;
    fTmpC = 2.4992186278915256197 * fZ * fTmpB - 0.96917220529352445597 * fTmpA;
    pSH(433) = fTmpC * fC1;
    pSH(409) = fTmpC * fS1;
    //! m = 13
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.81607711883762830929;
    pSH(196) = fTmpA * fC0;
    pSH(170) = -fTmpA * fS0;
    fTmpB = -4.3947097802721178604 * fZ;
    pSH(224) = fTmpB * fC0;
    pSH(198) = -fTmpB * fS0;
    fTmpC = -17.608243388844820556 * fZ2 + 0.60718080651189036399;
    pSH(254) = fTmpC * fC0;
    pSH(228) = -fTmpC * fS0;
    fTmpA = fZ * (-60.380154942952015813 * fZ2 + 5.8432408009308402400);
    pSH(286) = fTmpA * fC0;
    pSH(260) = -fTmpA * fS0;
    fTmpB = 3.1024184114977141490 * fZ * fTmpA - 0.90473663963430495830 * fTmpC;
    pSH(320) = fTmpB * fC0;
    pSH(294) = -fTmpB * fS0;
    fTmpC = 2.8904737863674562919 * fZ * fTmpB - 0.93168406158731500387 * fTmpA;
    pSH(356) = fTmpC * fC0;
    pSH(330) = -fTmpC * fS0;
    fTmpA = 2.7414640249326636021 * fZ * fTmpC - 0.94844798034925006032 * fTmpB;
    pSH(394) = fTmpA * fC0;
    pSH(368) = -fTmpA * fS0;
    fTmpB = 2.6309842116740119423 * fZ * fTmpA - 0.95970043296068228091 * fTmpC;
    pSH(434) = fTmpB * fC0;
    pSH(408) = -fTmpB * fS0;
    //! m = 14
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.83052208306452400305;
    pSH(225) = fTmpA * fC1;
    pSH(197) = fTmpA * fS1;
    fTmpB = 4.6241512566300114788 * fZ;
    pSH(255) = fTmpB * fC1;
    pSH(227) = fTmpB * fS1;
    fTmpC = 19.093881509360250421 * fZ2 - 0.61593166159226614262;
    pSH(287) = fTmpC * fC1;
    pSH(259) = fTmpC * fS1;
    fTmpA = fZ * (67.288948373844056316 * fZ2 - 6.1171771248949142105);
    pSH(321) = fTmpA * fC1;
    pSH(293) = fTmpA * fS1;
    fTmpB = 3.1807526624998681277 * fZ * fTmpA - 0.90256893466271141000 * fTmpC;
    pSH(357) = fTmpB * fC1;
    pSH(329) = fTmpB * fS1;
    fTmpC = 2.9572714696920445985 * fZ * fTmpB - 0.92973952503676233437 * fTmpA;
    pSH(395) = fTmpC * fC1;
    pSH(367) = fTmpC * fS1;
    fTmpA = 2.7996848562146502883 * fZ * fTmpC - 0.94671215845672612907 * fTmpB;
    pSH(435) = fTmpA * fC1;
    pSH(407) = fTmpA * fS1;
    //! m = 15
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.84425065085737263613;
    pSH(256) = fTmpA * fC0;
    pSH(226) = -fTmpA * fS0;
    fTmpB = -4.8498507532306824834 * fZ;
    pSH(288) = fTmpB * fC0;
    pSH(258) = -fTmpB * fS0;
    fTmpC = -20.602948605549728459 * fZ2 + 0.62433177592574934724;
    pSH(322) = fTmpC * fC0;
    pSH(292) = -fTmpC * fS0;
    fTmpA = fZ * (-74.515507921532000473 * fZ2 + 6.3870435361313143263);
    pSH(358) = fTmpA * fC0;
    pSH(328) = -fTmpA * fS0;
    fTmpB = 3.2573446421352252765 * fZ * fTmpA - 0.90063003157873466786 * fTmpC;
    pSH(396) = fTmpB * fC0;
    pSH(366) = -fTmpB * fS0;
    fTmpC = 3.0227707252027662084 * fZ * fTmpB - 0.92798615353802621180 * fTmpA;
    pSH(436) = fTmpC * fC0;
    pSH(406) = -fTmpC * fS0;
    //! m = 16
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.85734058883802509634;
    pSH(289) = fTmpA * fC1;
    pSH(257) = fTmpA * fS1;
    fTmpB = 5.0720953248553606107 * fZ;
    pSH(323) = fTmpB * fC1;
    pSH(291) = fTmpB * fS1;
    fTmpC = 22.134404124649146409 * fZ2 - 0.63241154641854704026;
    pSH(359) = fTmpC * fC1;
    pSH(327) = fTmpC * fS1;
    fTmpA = fZ * (82.055245832745453664 * fZ2 - 6.6531280404928746214);
    pSH(397) = fTmpA * fC1;
    pSH(365) = fTmpA * fS1;
    fTmpB = 3.3322915038553674929 * fZ * fTmpA - 0.89888569656853228599 * fTmpC;
    pSH(437) = fTmpB * fC1;
    pSH(405) = fTmpB * fS1;
    //! m = 17
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.86985717192062808222;
    pSH(324) = fTmpA * fC0;
    pSH(290) = -fTmpA * fS0;
    fTmpB = -5.2911346120699725413 * fZ;
    pSH(360) = fTmpB * fC0;
    pSH(326) = -fTmpB * fS0;
    fTmpC = -23.687309134978252702 * fZ2 + 0.64019754418860142438;
    pSH(398) = fTmpC * fC0;
    pSH(364) = -fTmpC * fS0;
    fTmpA = fZ * (-89.903887312140613733 * fZ2 + 6.9156836393954318256);
    pSH(438) = fTmpA * fC0;
    pSH(404) = -fTmpA * fS0;
    //! m = 18
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpA = 0.88185576867832886131;
    pSH(361) = fTmpA * fC1;
    pSH(325) = fTmpA * fS1;
    fTmpB = 5.5071875102722439487 * fZ;
    pSH(399) = fTmpB * fC1;
    pSH(363) = fTmpB * fS1;
    fTmpC = 25.260811578777112685 * fZ2 - 0.64771311740454135089;
    pSH(439) = fTmpC * fC1;
    pSH(403) = fTmpC * fS1;
    //! m = 19
    fC0 = fX * fC1 - fY * fS1;
    fS0 = fX * fS1 + fY * fC1;
    fTmpA = -0.89338378434994943301;
    pSH(400) = fTmpA * fC0;
    pSH(362) = -fTmpA * fS0;
    fTmpB = -5.7204473629006425018 * fZ;
    pSH(440) = fTmpB * fC0;
    pSH(402) = -fTmpB * fS0;
    //! m = 20
    fC1 = fX * fC0 - fY * fS0;
    fS1 = fX * fS0 + fY * fC0;
    fTmpC = 0.90448214509349098445;
    pSH(441) = fTmpC * fC1;
    pSH(401) = fTmpC * fS1;
  } break;
  default:
    break;
  }
#undef pSH
}

// n_coeff -- number of spline coefficients, must be 2 or 4
void spline_vector_v2_c_(double r_output, double *spl_param, int n_l_dim, int n_coeff, int n_grid_dim, int n_points,
                         int n_vector, double *out_result) {
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
void spline_vector_v2_c2_(double *r_output_, double *spl_param, int *n_l_dim_, int *n_coeff_, int *n_grid_dim_,
                          int *n_points_, int *n_vector_, double *out_result) {
  double r_output = *r_output_;
  int n_l_dim = *n_l_dim_;
  int n_coeff = *n_coeff_;
  int n_grid_dim = *n_grid_dim_;
  int n_points = *n_points_;
  int n_vector = *n_vector_;

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

// INFO 小心 spl_param 的大小，c/c++ 可没有 fortran 的自动切片，应准确地为 spl_param(n_l_dim,4,n_grid_dim)
void spline_vector_c_(double r_output, double *spl_param, int n_grid_dim, int n_l_dim, int n_points, int n_vector,
                      double *out_result) {
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

// info: non_peri_extd 是可选参数，由于 C 没有可选参数，不用请置空
void far_distance_hartree_fp_periodic_single_atom_c_(
    int current_atom, int i_center, double dist,
    int *l_hartree_max_far_distance, // (n_atoms)
    int inside, int forces_on, double multipole_radius_sq,
    double adap_outer_radius, // int *non_peri_extd
    // outer
    int l_pot_max, double *Fp, double b0[pmaxab + 1], double b2[pmaxab + 1], double b4[pmaxab + 1],
    double b6[pmaxab + 1], double a_save[pmaxab + 1], int hartree_force_l_add, int use_hartree_non_periodic_ewald,
    int hartree_fp_function_splines, int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc,
    double *Fp_function_spline_slice, double *Fpc_function_spline_slice) {
  double drel;

  int lmax = l_hartree_max_far_distance[current_atom - 1] + hartree_force_l_add;
  if (!use_hartree_non_periodic_ewald) {
    if (inside) {
      F_erf_c_(&Fp(0, i_center), dist, lmax, 0, hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min,
               Fp_grid_inc, Fp_function_spline_slice, Fpc_function_spline_slice);
      for (int i = 0; i <= lmax; i++)
        Fp(i, i_center) = -Fp(i, i_center);
    } else {
      F_erf_c_(&Fp(0, i_center), dist, lmax, 1, hartree_fp_function_splines, Fp_max_grid, lmax_Fp, Fp_grid_min,
               Fp_grid_inc, Fp_function_spline_slice, Fpc_function_spline_slice);
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

// F_erf + F_erfc
void F_erf_c_(double *F, double r, int p, int c,
              // outer
              int hartree_fp_function_splines, int Fp_max_grid, int lmax_Fp, double Fp_grid_min, double Fp_grid_inc,
              double *Fp_function_spline_slice, double *Fpc_function_spline_slice) {
  if (hartree_fp_function_splines) {
    // F_erf_spline + F_erfc_spline
    double rlog = invert_log_grid_c_(r, Fp_grid_min, Fp_grid_inc);
    double *spl_param = (c != 0) ? Fpc_function_spline_slice : Fp_function_spline_slice;
    spline_vector_c_(rlog, spl_param, Fp_max_grid, lmax_Fp + 1, Fp_max_grid, p + 1, F);
  } else {
    if (c)
      F_erfc_table_original_c_(F, r, p);
    else
      F_erf_table_original_c_(F, r, p);
  }
}

void F_erf_table_original_c_(double *F_erf_table, double r, int p_max) {
  printf("%s, not finished\n", __func__); // TODO
  exit(-19);
}
void F_erfc_table_original_c_(double *F_table, double r, int p_max) {
  printf("%s, not finished\n", __func__); // TODO
  exit(-19);
}

void far_distance_hartree_fp_cluster_single_atom_p2_c_(double *dist_tab_, int *l_max_, int *forces_on) {
  double dist_tab = *dist_tab_;
  int l_max = *l_max_;
  double dist_sq = dist_tab * dist_tab;
  int one_minus_2l = 1;
  Fp(0, 1) = 1.0 / dist_tab;
  for (int i_l = 1; i_l <= l_max + hartree_force_l_add; i_l++) {
    one_minus_2l -= 2;
    Fp(i_l, 1) = Fp(i_l - 1, 1) * (double)one_minus_2l / dist_sq;
  }
}

// void far_distance_real_hartree_potential_single_atom_p2_c_(int *i_center_, double *potential_, int *l_max_,
//                                                            double coord_current[3]) {
//   int i_center = *i_center_;
//   int l_max = *l_max_;

//   double dpot = 0.0;
//   double coord_c[3][(l_pot_max + 1)];
//   double dir[3];

//   coord_c[0][0] = 1.0;
//   coord_c[1][0] = 1.0;
//   coord_c[2][0] = 1.0;

//   dir[0] = coord_current[0] - coords_center(1, i_center);
//   dir[1] = coord_current[1] - coords_center(2, i_center);
//   dir[2] = coord_current[2] - coords_center(3, i_center);

//   int maxval = -1;
//   for (int i = 1; i <= 3; i++)
//     maxval = maxval > index_ijk_max_cc(i, l_max) ? maxval : index_ijk_max_cc(i, l_max);
//   for (int i_l = 1; i_l <= maxval; i_l++) {
//     coord_c[0][i_l] = dir[0] * coord_c[0][i_l - 1];
//     coord_c[1][i_l] = dir[1] * coord_c[1][i_l - 1];
//     coord_c[2][i_l] = dir[2] * coord_c[2][i_l - 1];
//   }
//   int index_cc_i_dim = n_cc_lm_ijk(l_max_analytic_multipole);
//   for (int n = 1; n <= n_cc_lm_ijk(l_max); n++) {
//     int ii = index_cc(n, 3, index_cc_i_dim);
//     int jj = index_cc(n, 4, index_cc_i_dim);
//     int kk = index_cc(n, 5, index_cc_i_dim);
//     int nn = index_cc(n, 6, index_cc_i_dim);
//     dpot =
//         dpot + coord_c[0][ii] * coord_c[1][jj] * coord_c[2][kk] * Fp(nn, 1) * multipole_c(n,
//         center_to_atom(i_center));
//   }
//   *potential_ += dpot;
// }

// void far_distance_real_hartree_potential_single_atom_c_(int *current_center_, int *i_center_, double *potential,
//                                                         int *l_hartree_max_far_distance, double coord_current[3]) {
//   int current_center = *current_center_;
//   int i_center = *i_center_;
//   double c_pot = 0.0;
//   double dpot = 0.0;
//   int l_max = l_hartree_max_far_distance[center_to_atom(current_center) - 1];
//   double coord_c[3][l_pot_max + 1];
//   double coord_mat[l_pot_max + 1][l_pot_max + 1];
//   double rest_mat[l_pot_max + 1][l_pot_max + 1];
//   double vector[n_cc_lm_ijk(l_pot_max)];
//   double dir[3];
//   coord_c[0][0] = 1.0;
//   coord_c[1][0] = 1.0;
//   coord_c[2][0] = 1.0;
//   dir[0] = coord_current[0] - coords_center(1, current_center);
//   dir[1] = coord_current[1] - coords_center(2, current_center);
//   dir[2] = coord_current[2] - coords_center(3, current_center);
//   for (int i_coord = 0; i_coord < 3; i_coord++)
//     for (int i_l = 1; i_l <= index_ijk_max_cc(i_coord + 1, l_max); i_l++)
//       coord_c[i_coord][i_l] = dir[i_coord] * coord_c[i_coord][i_l - 1];
//   for (int i_l = 0; i_l <= index_ijk_max_cc(1, l_max); i_l++)
//     for (int i_l2 = 0; i_l2 <= index_ijk_max_cc(2, l_max); i_l2++)
//       coord_mat[i_l2][i_l] = coord_c[0][i_l] * coord_c[1][i_l2];
//   for (int i_l = 0; i_l <= index_ijk_max_cc(3, l_max); i_l++)
//     for (int i_l2 = 0; i_l2 <= l_max; i_l2++)
//       rest_mat[i_l2][i_l] = coord_c[2][i_l] * Fp(i_l2, i_center);
//   int index_cc_i_dim = n_cc_lm_ijk(l_max_analytic_multipole);
//   for (int n = 1; n <= n_cc_lm_ijk(l_max); n++) {
//     int ii = index_cc(n, 3, index_cc_i_dim);
//     int jj = index_cc(n, 4, index_cc_i_dim);
//     int kk = index_cc(n, 5, index_cc_i_dim);
//     int nn = index_cc(n, 6, index_cc_i_dim);
//     vector[n - 1] = coord_mat[jj][ii] * rest_mat[nn][kk];
//   }
//   for (int n = 0; n < n_cc_lm_ijk(l_max); n++)
//     dpot += vector[n] * multipole_c(n + 1, center_to_atom(current_center));
//   if (fabs(dpot) > 1e-30)
//     c_pot += dpot;
//   *potential += c_pot;
// }