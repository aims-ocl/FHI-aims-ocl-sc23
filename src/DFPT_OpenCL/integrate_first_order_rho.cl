
#ifdef fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

// #include <math.h>
// #define kernel
// #define global
// #define local
#define mconstant global const
#define atomic_cmpxchg atom_cmpxchg
#define NULL 0

#define H_IC_TILE 16
#define H_JC_TILE 16
#define H_PT_TILE 16
#define H_IC_WORK 3
#define H_JC_WORK 3

typedef double m_float_type;
#define LOCAL_SIZE 256
#define TILE_N 16 // get_local_size(0) 必须为 16 的倍数
#define TILE_K 16 // TILE_K == TILE_N
#define WORK_M 4
#define WORK_N 8
#define TILE_M (LOCAL_SIZE / TILE_N)

void SHEvalderiv_c_(int lmax, double sintheta, double costheta, double sinphi, double cosphi, global double *pSH,
                    global double *pdSHdtheta, global double *pdSHdphi_sintheta);

void prune_radial_basis_p2_c_(int *dim_atoms_, int *dim_fns_, global double *dist_tab_sq, global double *dist_tab,
                              global double *dir_tab, // (3, n_atom_list)
                              int *n_compute_atoms_, global int *atom_index, global int *atom_index_inv, int *n_compute_fns_,
                              global int *i_basis_fns, global int *i_basis_fns_inv, // (n_basis_fns,n_centers)
                              global int *i_atom_fns, global int *spline_array_start, global int *spline_array_end, int *n_atom_list_,
                              global int *atom_list, int *n_compute_, mconstant int *i_basis, global int *n_batch_centers_, global int *batch_center,
                              private double *one_over_dist_tab, private int *rad_index, global int *wave_index, global int *l_index, global int *l_count,
                              global int *fn_atom, int *n_zero_compute_, global int *zero_index_point,
                              // outer
                              int n_basis_fns, // const
                              mconstant int *center_to_atom, mconstant int *species_center, mconstant int *Cbasis_to_basis, mconstant int *Cbasis_to_center,
                              mconstant int *perm_basis_fns_spl, mconstant double *outer_radius_sq, mconstant int *basis_fn, mconstant int *basis_l,
                              mconstant double *atom_radius_sq, mconstant int *basis_fn_start_spl, mconstant int *basis_fn_atom,
                              global double* pbc_lists_coords_center, 
                              double coords_center0, double coords_center1, double coords_center2) {
// #define dir_tab(i, j) dir_tab[(i)-1 + 3 * ((j)-1)]
#define dir_tab(i, j) dir_tab[lid + lsize * ((i)-1 + 3 * ((j)-1))]
// #define i_basis_fns_inv(i, j) i_basis_fns_inv[(i)-1 + n_basis_fns * ((j)-1)]
#define i_basis_fns_inv(i, j) i_basis_fns_inv[((i)-1 + n_basis_fns * ((j)-1)) * gsize + gid]
// outer
#define center_to_atom(i) center_to_atom[(i)-1]
#define species_center(i) species_center[(i)-1]
#define Cbasis_to_basis(i) Cbasis_to_basis[(i)-1]
#define Cbasis_to_center(i) Cbasis_to_center[(i)-1]
#define perm_basis_fns_spl(i) perm_basis_fns_spl[(i)-1]
#define outer_radius_sq(i) outer_radius_sq[(i)-1]
#define basis_fn(i) basis_fn[(i)-1]
#define basis_l(i) basis_l[(i)-1]
#define atom_radius_sq(i) atom_radius_sq[(i)-1]
#define basis_fn_start_spl(i) basis_fn_start_spl[(i)-1]
#define basis_fn_atom(i, j) basis_fn_atom[(i)-1 + ((j)-1) * n_basis_fns]
#define pbc_lists_coords_center(i, j) pbc_lists_coords_center[((j)-1) * 3 + (i)-1]
#define atom_index_inv(i) atom_index_inv[((i)-1) * lsize + lid]
  int gid = get_global_id(0);
  int gsize = get_global_size(0);
  int lid = get_local_id(0);
  int lsize = get_local_size(0);

  int dim_atoms = *dim_atoms_;
  int dim_fns = *dim_fns_;
  int n_atom_list = *n_atom_list_;
  int n_batch_centers = *n_batch_centers_;
  int n_compute = *n_compute_;
  int n_compute_atoms = *n_compute_atoms_;
  int n_compute_fns = *n_compute_fns_;
  int n_zero_compute = *n_zero_compute_;

  // local variables
  int i_offset_spl;
  int l_aux;

  // counters
  int i_batch_center;
  int i_basis_1;
  int i_compute;
  int i_fn;
  int i_lm;
  int i_center, i_center_L;
  int i_atom_compute;
  int i_spline;

  double dir_tab_local[3];
  for (i_batch_center = 1; i_batch_center <= n_batch_centers; i_batch_center++) {
    i_center_L = batch_center[i_batch_center - 1];
    i_center = atom_list[i_center_L - 1];

    dir_tab_local[0] = coords_center0 - pbc_lists_coords_center(1, atom_list[i_center_L - 1]);
    dir_tab_local[1] = coords_center1 - pbc_lists_coords_center(2, atom_list[i_center_L - 1]);
    dir_tab_local[2] = coords_center2 - pbc_lists_coords_center(3, atom_list[i_center_L - 1]);

    double dist_tab_sq_now = dir_tab_local[0] * dir_tab_local[0]
                           + dir_tab_local[1] * dir_tab_local[1]
                           + dir_tab_local[2] * dir_tab_local[2];

    if (dist_tab_sq_now <= atom_radius_sq(species_center(i_center)) && dist_tab_sq_now > 0.0) {
      n_compute_atoms = n_compute_atoms + 1;
      atom_index[(n_compute_atoms - 1)] = i_center;
      // atom_index_inv(i_center) = n_compute_atoms;
      dist_tab_sq[n_compute_atoms - 1] = dist_tab_sq_now;
      dist_tab[n_compute_atoms - 1] = sqrt(dist_tab_sq[n_compute_atoms - 1]);
      double tmp = 1.0 / dist_tab[n_compute_atoms - 1];
      one_over_dist_tab[n_compute_atoms - 1] = 1.0 / dist_tab[n_compute_atoms - 1];
      for (int i = 1; i <= 3; i++) {
        // dir_tab(i, n_compute_atoms) = dir_tab(i, i_center_L) * tmp;
        dir_tab(i, n_compute_atoms) = dir_tab_local[i-1] * tmp;
      }
    } else {
      // atom_index_inv[i_center - 1] = 0;
      // atom_index_inv[i_center - 1] = dim_atoms + 1; // dim_atoms == n_max_compute_atoms
    }
  }
  // private int private_center_to_atom[MACRO_n_centers+1];
  // // 确保非法的部分将被映射到从不会被使用的 MACRO_n_max_compute_atoms + 1 上，由此保证其值为 0
  // for (int i=0; i<MACRO_n_centers+1; i++){
  //   private_center_to_atom[i] = MACRO_n_max_compute_atoms + 1;
  // }
  // next, check for radial basis functions
  for (i_atom_compute = 1; i_atom_compute <= n_compute_atoms; i_atom_compute++) {
    i_center = atom_index[(i_atom_compute - 1)];
    // private_center_to_atom[i_center] = i_atom_compute;
    i_offset_spl = basis_fn_start_spl(species_center(i_center));
    double dist_tab_sq_reg = dist_tab_sq[i_atom_compute - 1];
    for (i_spline = 1; i_spline <= n_basis_fns; i_spline++) {
      i_basis_1 = perm_basis_fns_spl(i_spline);
      if (basis_fn_atom(i_basis_1, center_to_atom(i_center))) {
        if (dist_tab_sq_reg <= outer_radius_sq(i_basis_1)) {
          spline_array_start[i_atom_compute - 1] = (i_spline - 1) + 1;
          break;
        }
      }
    }

    for (i_basis_1 = 1; i_basis_1 <= n_basis_fns; i_basis_1++) {
      if (basis_fn_atom(i_basis_1, center_to_atom(i_center))) {
        if (dist_tab_sq_reg <= outer_radius_sq(i_basis_1)) {
          n_compute_fns = n_compute_fns + 1;
          // i_basis_fns[n_compute_fns - 1] = i_basis_1;
          // i_atom_fns[n_compute_fns - 1] = i_center;
          // i_basis_fns_inv(i_basis_1, i_center) = n_compute_fns;
          i_basis_fns_inv(i_basis_1, i_atom_compute) = n_compute_fns;
          fn_atom[n_compute_fns - 1] = i_atom_compute;
        }
        i_offset_spl = i_offset_spl + 1;
      }
    }
    rad_index[i_atom_compute - 1] = n_compute_fns;
    spline_array_end[i_atom_compute - 1] = i_offset_spl - 1;
  }
  // athread 版本不考虑下面这种情况，因此保持一致 WARNING
  // // 担忧出现 i_center = atom_index[(i_atom_compute - 1)]; 将不同的 i_atom_compute 映射到相同的 i_center
  // // 为了确保后文其他函数调用是正常，额外传一个 i_basis_fns_inv_remapping 出去
  // // 暂时不选用重新给 i_basis_fns_inv 赋值的方式
  // for (i_atom_compute = 1; i_atom_compute <= n_compute_atoms; i_atom_compute++) {
  //   i_basis_fns_inv_remapping[i_atom_compute-1] = private_center_to_atom[atom_index[(i_atom_compute - 1)]];
  // }

  n_zero_compute = 0;
  for (int i = 1; i <= n_compute_fns; i++) {
    wave_index[i - 1] = 0;
  }

  i_compute = 1;
  while (i_compute <= n_compute) {
    i_basis_1 = i_basis[i_compute - 1];
    // i_fn = i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), Cbasis_to_center(i_basis_1));
    // i_fn = i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), private_center_to_atom[Cbasis_to_center(i_basis_1)]);
    // i_fn = i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), atom_index_inv(Cbasis_to_center(i_basis_1)));
    int cc = dim_atoms + 1;
    int i_center = Cbasis_to_center(i_basis_1);
    for(int i=0; i<dim_atoms; i++){
      if(atom_index[i] == i_center){
        cc = i+1;
        break;
      }
    }
    i_fn = i_basis_fns_inv(basis_fn(Cbasis_to_basis(i_basis_1)), cc);
    l_aux = basis_l(Cbasis_to_basis(i_basis[i_compute - 1]));
    if (i_fn == 0) {
      for (i_lm = 0; i_lm <= 2 * l_aux; i_lm++) {
        n_zero_compute = n_zero_compute + 1;
        zero_index_point[n_zero_compute - 1] = i_compute + i_lm;
      }
    } else if (wave_index[i_fn - 1] == 0) {
      wave_index[i_fn - 1] = i_compute;
      // l_index[i_fn - 1] = l_aux * l_aux + 1;
      // l_count[i_fn - 1] = 2 * l_aux;
      l_count[i_fn - 1] = l_aux;
    }
    i_compute = i_compute + 2 * l_aux + 1;
  }
  *n_compute_atoms_ = n_compute_atoms;
  *n_compute_fns_ = n_compute_fns;
  *n_zero_compute_ = n_zero_compute;

#undef center_to_atom
#undef species_center
#undef Cbasis_to_basis
#undef Cbasis_to_center
#undef perm_basis_fns_spl
#undef outer_radius_sq
#undef basis_fn
#undef basis_l
#undef atom_radius_sq
#undef basis_fn_start_spl
#undef basis_fn_atom
#undef pbc_lists_coords_center
#undef atom_index_inv

#undef dir_tab
#undef i_basis_fns_inv
}

void tab_local_geometry_p2_c_(int *n_compute_atoms_, global int *atom_index, global double *dist_tab, global double *i_r,
                              // outer
                              mconstant int *species_center, global double *r_grid_min, global double *log_r_grid_inc) {
#define species_center(i) species_center[(i)-1]
#define r_grid_min(i) r_grid_min[(i)-1]
#define log_r_grid_inc(i) log_r_grid_inc[(i)-1]
  int gsize = get_global_size(0);
  int gid = get_global_id(0);
  int n_compute_atoms = *n_compute_atoms_;
  //  counters
  int i_compute_atom;
  for (i_compute_atom = 1; i_compute_atom <= n_compute_atoms; i_compute_atom++) {
    double r_current = dist_tab[i_compute_atom - 1];
    int i_species = species_center(atom_index[(i_compute_atom - 1)]);
    i_r[i_compute_atom - 1] = 1.0 + log(r_current / r_grid_min(i_species)) / log_r_grid_inc(i_species);
  }
#undef species_center
#undef r_grid_min
#undef log_r_grid_inc
}

void tab_trigonom_p0_c_(int *n_compute_atoms_, global double *dir_tab, global double *trigonom_tab) {
#define dir_tab(i, j) dir_tab[(i)-1 + ((j)-1) * 3]
#define trigonom_tab(i, j) trigonom_tab[(i)-1 + ((j)-1) * 4]
  int n_compute_atoms = *n_compute_atoms_;
  //  local variables
  double abmax, abcmax, ab, abc;

  //  counters
  int i_atom;
  int i_coord;

  for (i_atom = 1; i_atom <= n_compute_atoms; i_atom++) {
    abmax = fmax(fabs(dir_tab(1, i_atom)), fabs(dir_tab(2, i_atom)));
    if (abmax > 0.0) {
      ab = sqrt(pow(dir_tab(1, i_atom), 2.0) + pow(dir_tab(2, i_atom), 2.0));
      trigonom_tab(4, i_atom) = dir_tab(1, i_atom) / ab;
      trigonom_tab(3, i_atom) = dir_tab(2, i_atom) / ab;
    } else {
      trigonom_tab(4, i_atom) = 1.0;
      trigonom_tab(3, i_atom) = 0.0;
      ab = 0.0;
    }
    abcmax = fmax(abmax, fabs(dir_tab(3, i_atom)));
    if (abcmax > 0.0) {
      abc = sqrt(pow(ab, 2.0) + pow(dir_tab(3, i_atom), 2.0));
      trigonom_tab(2, i_atom) = dir_tab(3, i_atom) / abc;
      trigonom_tab(1, i_atom) = ab / abc;
    } else {
      trigonom_tab(2, i_atom) = 1.0;
      trigonom_tab(1, i_atom) = 0.0;
    }
  }
#undef dir_tab
#undef trigonom_tab
}

void tab_gradient_ylm_p0_c_(global double *trigonom_tab, // ( 4, n_compute_atoms )
                            global int *basis_l_max, int *l_ylm_max_, int *n_compute_atoms_, global int *atom_index,
                            global double *ylm_tab,             // ( (l_ylm_max+1)**2, n_compute_atoms )
                            global double *dylm_dtheta_tab,     // ( (l_ylm_max+1)**2, n_compute_atoms )
                            global double *scaled_dylm_dphi_tab // ( (l_ylm_max+1)**2, n_compute_atoms )
                            // outer
                            ,
                            global double *dir_tab,
                            mconstant int *species_center) {
#define species_center(i) species_center[(i)-1]
// #define dir_tab(i, j) dir_tab[(i)-1 + ((j)-1) * 3]
#define dir_tab(i, j) dir_tab[lid + lsize * ((i)-1 + ((j)-1) * 3)]
  int n_compute_atoms = *n_compute_atoms_;
  int l_ylm_max = *l_ylm_max_;
  int l_ylm_max_1pow2 = (l_ylm_max + 1) * (l_ylm_max + 1);
  int gid = get_global_id(0);
  int gsize = get_global_size(0);
  int lid = get_local_id(0);
  int lsize = get_local_size(0);
#define trigonom_tab(i, j) trigonom_tab[(i)-1 + ((j)-1) * 4]
#define ylm_tab(i, j) ylm_tab[i - 1 + l_ylm_max_1pow2 * (j - 1)]
#define dylm_dtheta_tab(i, j) dylm_dtheta_tab[i - 1 + l_ylm_max_1pow2 * (j - 1)]
#define scaled_dylm_dphi_tab(i, j) scaled_dylm_dphi_tab[i - 1 + l_ylm_max_1pow2 * (j - 1)]
  for (int i_atom = 1; i_atom <= n_compute_atoms; i_atom++) {
    double trigonom_tab_reg[4];
    {
      //  local variables
      double abmax, abcmax, ab, abc;
      abmax = fmax(fabs(dir_tab(1, i_atom)), fabs(dir_tab(2, i_atom)));
      if (abmax > 0.0) {
        ab = sqrt(pow(dir_tab(1, i_atom), 2.0) + pow(dir_tab(2, i_atom), 2.0));
        trigonom_tab_reg[3] = dir_tab(1, i_atom) / ab;
        trigonom_tab_reg[2] = dir_tab(2, i_atom) / ab;
      } else {
        trigonom_tab_reg[3] = 1.0;
        trigonom_tab_reg[2] = 0.0;
        ab = 0.0;
      }
      abcmax = fmax(abmax, fabs(dir_tab(3, i_atom)));
      if (abcmax > 0.0) {
        abc = sqrt(pow(ab, 2.0) + pow(dir_tab(3, i_atom), 2.0));
        trigonom_tab_reg[1] = dir_tab(3, i_atom) / abc;
        trigonom_tab_reg[0] = ab / abc;
      } else {
        trigonom_tab_reg[1] = 1.0;
        trigonom_tab_reg[0] = 0.0;
      }
    }
    //     increment_ylm_deriv(trigonom_tab(1, i_atom), trigonom_tab(2, i_atom), trigonom_tab(3, i_atom),
    //                         trigonom_tab(4, i_atom), 0, basis_l_max(species_center(atom_index(i_atom))),
    //                         &ylm_tab(1, i_atom), &dylm_dtheta_tab(1, i_atom), &scaled_dylm_dphi_tab(1, i_atom));
    if (1) {
      // SHEvalderiv_c_(basis_l_max[species_center(atom_index[i_atom - 1]) - 1], trigonom_tab(1, i_atom),
      //                trigonom_tab(2, i_atom), trigonom_tab(3, i_atom), trigonom_tab(4, i_atom), &ylm_tab(1, i_atom),
      //                &dylm_dtheta_tab(1, i_atom), &scaled_dylm_dphi_tab(1, i_atom));
      {
              int l_atom_max = basis_l_max[species_center(atom_index[(i_atom - 1)]) - 1];

              /* variables for tabulate ylm */
              double YLLI, YLL1I, YL1L1I, YLMI;
              double YLLR, YLL1R, YL1L1R, YLMR;
              int I2L, I4L2, INDEX, INDEX2, L, M, MSIGN;
              /* VB */
              int I22L, I24L2;
              double TEMP1, TEMP2, TEMP3;

              double D4LL1C, D2L13;
              const double PI = 3.14159265358979323846;

              #define trigonom_tab_(i1) trigonom_tab_reg[i1-1]
              #define ylm_tab_(i) ylm_tab(i, i_atom)
                    if (0 <= 0)
                    {
                        YLLR = 1.0/sqrt(4.0*PI);
                        YLLI = 0.0;
                        ylm_tab_(1) = YLLR;
                    }

                    if ( (0 <= 1) && (l_atom_max >= 1))
                    {
                        ylm_tab_(3) = sqrt(3.00)*YLLR*trigonom_tab_(2);
                        TEMP1 = -sqrt(3.00)*YLLR*trigonom_tab_(1);
                        ylm_tab_(4) = TEMP1*trigonom_tab_(4);
                        ylm_tab_(2) = -TEMP1*trigonom_tab_(3);
                    }

                    // L = max(2,0)
                    for (L = 2; L <= l_atom_max; L++)
                    {
                        INDEX  = L*L+1;
                        INDEX2 = INDEX + 2*L;
                        MSIGN  = 1 - 2*(L%2);

                        YL1L1R = ylm_tab_(INDEX-1);
                        YL1L1I = - MSIGN * ylm_tab_(INDEX-2*L+1);
                        TEMP1 = -sqrt((double)(2*L+1)/(double)(2*L))*trigonom_tab_(1);
                        YLLR = TEMP1*(trigonom_tab_(4)*YL1L1R - trigonom_tab_(3)*YL1L1I);
                        YLLI = TEMP1*(trigonom_tab_(4)*YL1L1I + trigonom_tab_(3)*YL1L1R);
                        ylm_tab_(INDEX2) = YLLR;
                        ylm_tab_(INDEX)  = MSIGN * YLLI;
                        INDEX2 = INDEX2 - 1;
                        INDEX  = INDEX  + 1;

                        TEMP2 = sqrt((double)(2*L+1))*trigonom_tab_(2);
                        YLL1R = TEMP2*YL1L1R;
                        YLL1I = TEMP2*YL1L1I;
                        ylm_tab_(INDEX2) = YLL1R;
                        ylm_tab_(INDEX)  = - MSIGN * YLL1I;
                        INDEX2 = INDEX2 - 1;
                        INDEX  = INDEX  + 1;

                        I4L2 = INDEX - 4*L + 2;
                        I2L  = INDEX - 2*L;
                        I24L2 = INDEX2 - 4*L + 2;
                        I22L  = INDEX2 - 2*L;
                        D4LL1C = trigonom_tab_(2)*sqrt((double)(4*L*L-1));
                        D2L13  = -sqrt((double)(2*L+1)/(double)(2*L-3));

                        for (M = L - 2; M >= 0; M--)
                        {
                            TEMP1 = 1.00/sqrt((double)((L+M)*(L-M)));
                            TEMP2 = D4LL1C*TEMP1;
                            TEMP3 = D2L13*sqrt((double)((L+M-1)*(L-M-1)))*TEMP1;
                            YLMR = TEMP2*ylm_tab_(I22L) + TEMP3*ylm_tab_(I24L2);
                            YLMI = TEMP2*ylm_tab_(I2L) + TEMP3*ylm_tab_(I4L2);
                            ylm_tab_(INDEX2) = YLMR;
                            ylm_tab_(INDEX)  = YLMI;

                            INDEX2 = INDEX2 - 1;
                            INDEX  = INDEX  + 1;
                            I24L2   = I24L2   - 1;
                            I22L    = I22L    - 1;
                            I4L2   = I4L2   + 1;
                            I2L    = I2L    + 1;
                        }
                    }
              #undef trigonom_tab_
              #undef ylm_tab_
      }
    } else {
      // printf("%s, not finished\n", __func__); // TODO
      // exit(-19);
    }
  }
  const int REL_x2c = 7;
  const int REL_q4c = 8;
  // if (flag_rel == REL_q4c || flag_rel == REL_x2c) {
  if (0) {
    for (int i_atom = 1; i_atom <= n_compute_atoms; i_atom++) {
      if (l_ylm_max >= 1) {
        ylm_tab(4, i_atom) = -ylm_tab(4, i_atom);
        dylm_dtheta_tab(4, i_atom) = -dylm_dtheta_tab(4, i_atom);
        scaled_dylm_dphi_tab(4, i_atom) = -scaled_dylm_dphi_tab(4, i_atom);
      }
      if (l_ylm_max >= 2) {
        ylm_tab(8, i_atom) = -ylm_tab(8, i_atom);
        dylm_dtheta_tab(8, i_atom) = -dylm_dtheta_tab(8, i_atom);
        scaled_dylm_dphi_tab(8, i_atom) = -scaled_dylm_dphi_tab(8, i_atom);
      }
      if (l_ylm_max >= 3) {
        ylm_tab(14, i_atom) = -ylm_tab(14, i_atom);
        dylm_dtheta_tab(14, i_atom) = -dylm_dtheta_tab(14, i_atom);
        scaled_dylm_dphi_tab(14, i_atom) = -scaled_dylm_dphi_tab(14, i_atom);
        ylm_tab(16, i_atom) = -ylm_tab(16, i_atom);
        dylm_dtheta_tab(16, i_atom) = -dylm_dtheta_tab(16, i_atom);
        scaled_dylm_dphi_tab(16, i_atom) = -scaled_dylm_dphi_tab(16, i_atom);
      }
      if (l_ylm_max >= 4) {
        ylm_tab(22, i_atom) = -ylm_tab(22, i_atom);
        dylm_dtheta_tab(22, i_atom) = -dylm_dtheta_tab(22, i_atom);
        scaled_dylm_dphi_tab(22, i_atom) = -scaled_dylm_dphi_tab(22, i_atom);
        ylm_tab(24, i_atom) = -ylm_tab(24, i_atom);
        dylm_dtheta_tab(24, i_atom) = -dylm_dtheta_tab(24, i_atom);
        scaled_dylm_dphi_tab(24, i_atom) = -scaled_dylm_dphi_tab(24, i_atom);
      }
      if (l_ylm_max >= 5) {
        ylm_tab(32, i_atom) = -ylm_tab(32, i_atom);
        dylm_dtheta_tab(32, i_atom) = -dylm_dtheta_tab(32, i_atom);
        scaled_dylm_dphi_tab(32, i_atom) = -scaled_dylm_dphi_tab(32, i_atom);
        ylm_tab(34, i_atom) = -ylm_tab(34, i_atom);
        dylm_dtheta_tab(34, i_atom) = -dylm_dtheta_tab(34, i_atom);
        scaled_dylm_dphi_tab(34, i_atom) = -scaled_dylm_dphi_tab(34, i_atom);
        ylm_tab(36, i_atom) = -ylm_tab(36, i_atom);
        dylm_dtheta_tab(36, i_atom) = -dylm_dtheta_tab(36, i_atom);
        scaled_dylm_dphi_tab(36, i_atom) = -scaled_dylm_dphi_tab(36, i_atom);
      }
      if (l_ylm_max >= 6) {
        ylm_tab(44, i_atom) = -ylm_tab(44, i_atom);
        dylm_dtheta_tab(44, i_atom) = -dylm_dtheta_tab(44, i_atom);
        scaled_dylm_dphi_tab(44, i_atom) = -scaled_dylm_dphi_tab(44, i_atom);
        ylm_tab(46, i_atom) = -ylm_tab(46, i_atom);
        dylm_dtheta_tab(46, i_atom) = -dylm_dtheta_tab(46, i_atom);
        scaled_dylm_dphi_tab(46, i_atom) = -scaled_dylm_dphi_tab(46, i_atom);
        ylm_tab(48, i_atom) = -ylm_tab(48, i_atom);
        dylm_dtheta_tab(48, i_atom) = -dylm_dtheta_tab(48, i_atom);
        scaled_dylm_dphi_tab(48, i_atom) = -scaled_dylm_dphi_tab(48, i_atom);
      }
    }
  }
#undef trigonom_tab
#undef ylm_tab
#undef dylm_dtheta_tab
#undef scaled_dylm_dphi_tab

#undef species_center
#undef dir_tab
}

void tab_gradient_ylm_p0_c_2(global double *trigonom_tab, // ( 4, n_compute_atoms )
                            global int *basis_l_max, int *l_ylm_max_, int *n_compute_atoms_, global int *atom_index,
                            global double *ylm_tab,             // ( (l_ylm_max+1)**2, n_compute_atoms )
                            global double *dylm_dtheta_tab,     // ( (l_ylm_max+1)**2, n_compute_atoms )
                            global double *scaled_dylm_dphi_tab // ( (l_ylm_max+1)**2, n_compute_atoms )
                            // outer
                            ,
                            global double *dir_tab,
                            mconstant int *species_center) {
#define species_center(i) species_center[(i)-1]
// #define dir_tab(i, j) dir_tab[(i)-1 + ((j)-1) * 3]
#define dir_tab(i, j) dir_tab[lid + lsize * ((i)-1 + ((j)-1) * 3)]
  int n_compute_atoms = *n_compute_atoms_;
  int l_ylm_max = *l_ylm_max_;
  int l_ylm_max_1pow2 = (l_ylm_max + 1) * (l_ylm_max + 1);
  int gid = get_global_id(0);
  int gsize = get_global_size(0);
  int lid = get_local_id(0);
  int lsize = get_local_size(0);
#define trigonom_tab(i, j) trigonom_tab[(i)-1 + ((j)-1) * 4]
#define ylm_tab(i, j) ylm_tab[(i - 1 + l_ylm_max_1pow2 * (j - 1)) * gsize + gid]
#define dylm_dtheta_tab(i, j) dylm_dtheta_tab[i - 1 + l_ylm_max_1pow2 * (j - 1)]
#define scaled_dylm_dphi_tab(i, j) scaled_dylm_dphi_tab[i - 1 + l_ylm_max_1pow2 * (j - 1)]
  for (int i_atom = 1; i_atom <= n_compute_atoms; i_atom++) {
    double trigonom_tab_reg[4];
    {
      //  local variables
      double abmax, abcmax, ab, abc;
      abmax = fmax(fabs(dir_tab(1, i_atom)), fabs(dir_tab(2, i_atom)));
      if (abmax > 0.0) {
        ab = sqrt(pow(dir_tab(1, i_atom), 2.0) + pow(dir_tab(2, i_atom), 2.0));
        trigonom_tab_reg[3] = dir_tab(1, i_atom) / ab;
        trigonom_tab_reg[2] = dir_tab(2, i_atom) / ab;
      } else {
        trigonom_tab_reg[3] = 1.0;
        trigonom_tab_reg[2] = 0.0;
        ab = 0.0;
      }
      abcmax = fmax(abmax, fabs(dir_tab(3, i_atom)));
      if (abcmax > 0.0) {
        abc = sqrt(pow(ab, 2.0) + pow(dir_tab(3, i_atom), 2.0));
        trigonom_tab_reg[1] = dir_tab(3, i_atom) / abc;
        trigonom_tab_reg[0] = ab / abc;
      } else {
        trigonom_tab_reg[1] = 1.0;
        trigonom_tab_reg[0] = 0.0;
      }
    }
    //     increment_ylm_deriv(trigonom_tab(1, i_atom), trigonom_tab(2, i_atom), trigonom_tab(3, i_atom),
    //                         trigonom_tab(4, i_atom), 0, basis_l_max(species_center(atom_index(i_atom))),
    //                         &ylm_tab(1, i_atom), &dylm_dtheta_tab(1, i_atom), &scaled_dylm_dphi_tab(1, i_atom));
    if (1) {
      // SHEvalderiv_c_(basis_l_max[species_center(atom_index[i_atom - 1]) - 1], trigonom_tab(1, i_atom),
      //                trigonom_tab(2, i_atom), trigonom_tab(3, i_atom), trigonom_tab(4, i_atom), &ylm_tab(1, i_atom),
      //                &dylm_dtheta_tab(1, i_atom), &scaled_dylm_dphi_tab(1, i_atom));
      {
              int l_atom_max = basis_l_max[species_center(atom_index[(i_atom - 1)]) - 1];

              /* variables for tabulate ylm */
              double YLLI, YLL1I, YL1L1I, YLMI;
              double YLLR, YLL1R, YL1L1R, YLMR;
              int I2L, I4L2, INDEX, INDEX2, L, M, MSIGN;
              /* VB */
              int I22L, I24L2;
              double TEMP1, TEMP2, TEMP3;

              double D4LL1C, D2L13;
              const double PI = 3.14159265358979323846;

              #define trigonom_tab_(i1) trigonom_tab_reg[i1-1]
              #define ylm_tab_(i) ylm_tab(i, i_atom)
                    if (0 <= 0)
                    {
                        YLLR = 1.0/sqrt(4.0*PI);
                        YLLI = 0.0;
                        ylm_tab_(1) = YLLR;
                    }

                    if ( (0 <= 1) && (l_atom_max >= 1))
                    {
                        ylm_tab_(3) = sqrt(3.00)*YLLR*trigonom_tab_(2);
                        TEMP1 = -sqrt(3.00)*YLLR*trigonom_tab_(1);
                        ylm_tab_(4) = TEMP1*trigonom_tab_(4);
                        ylm_tab_(2) = -TEMP1*trigonom_tab_(3);
                    }

                    // L = max(2,0)
                    for (L = 2; L <= l_atom_max; L++)
                    {
                        INDEX  = L*L+1;
                        INDEX2 = INDEX + 2*L;
                        MSIGN  = 1 - 2*(L%2);

                        YL1L1R = ylm_tab_(INDEX-1);
                        YL1L1I = - MSIGN * ylm_tab_(INDEX-2*L+1);
                        TEMP1 = -sqrt((double)(2*L+1)/(double)(2*L))*trigonom_tab_(1);
                        YLLR = TEMP1*(trigonom_tab_(4)*YL1L1R - trigonom_tab_(3)*YL1L1I);
                        YLLI = TEMP1*(trigonom_tab_(4)*YL1L1I + trigonom_tab_(3)*YL1L1R);
                        ylm_tab_(INDEX2) = YLLR;
                        ylm_tab_(INDEX)  = MSIGN * YLLI;
                        INDEX2 = INDEX2 - 1;
                        INDEX  = INDEX  + 1;

                        TEMP2 = sqrt((double)(2*L+1))*trigonom_tab_(2);
                        YLL1R = TEMP2*YL1L1R;
                        YLL1I = TEMP2*YL1L1I;
                        ylm_tab_(INDEX2) = YLL1R;
                        ylm_tab_(INDEX)  = - MSIGN * YLL1I;
                        INDEX2 = INDEX2 - 1;
                        INDEX  = INDEX  + 1;

                        I4L2 = INDEX - 4*L + 2;
                        I2L  = INDEX - 2*L;
                        I24L2 = INDEX2 - 4*L + 2;
                        I22L  = INDEX2 - 2*L;
                        D4LL1C = trigonom_tab_(2)*sqrt((double)(4*L*L-1));
                        D2L13  = -sqrt((double)(2*L+1)/(double)(2*L-3));

                        for (M = L - 2; M >= 0; M--)
                        {
                            TEMP1 = 1.00/sqrt((double)((L+M)*(L-M)));
                            TEMP2 = D4LL1C*TEMP1;
                            TEMP3 = D2L13*sqrt((double)((L+M-1)*(L-M-1)))*TEMP1;
                            YLMR = TEMP2*ylm_tab_(I22L) + TEMP3*ylm_tab_(I24L2);
                            YLMI = TEMP2*ylm_tab_(I2L) + TEMP3*ylm_tab_(I4L2);
                            ylm_tab_(INDEX2) = YLMR;
                            ylm_tab_(INDEX)  = YLMI;

                            INDEX2 = INDEX2 - 1;
                            INDEX  = INDEX  + 1;
                            I24L2   = I24L2   - 1;
                            I22L    = I22L    - 1;
                            I4L2   = I4L2   + 1;
                            I2L    = I2L    + 1;
                        }
                    }
              #undef trigonom_tab_
              #undef ylm_tab_
      }
    } else {
      // printf("%s, not finished\n", __func__); // TODO
      // exit(-19);
    }
  }
  const int REL_x2c = 7;
  const int REL_q4c = 8;
  // if (flag_rel == REL_q4c || flag_rel == REL_x2c) {
  if (0) {
    for (int i_atom = 1; i_atom <= n_compute_atoms; i_atom++) {
      if (l_ylm_max >= 1) {
        ylm_tab(4, i_atom) = -ylm_tab(4, i_atom);
        dylm_dtheta_tab(4, i_atom) = -dylm_dtheta_tab(4, i_atom);
        scaled_dylm_dphi_tab(4, i_atom) = -scaled_dylm_dphi_tab(4, i_atom);
      }
      if (l_ylm_max >= 2) {
        ylm_tab(8, i_atom) = -ylm_tab(8, i_atom);
        dylm_dtheta_tab(8, i_atom) = -dylm_dtheta_tab(8, i_atom);
        scaled_dylm_dphi_tab(8, i_atom) = -scaled_dylm_dphi_tab(8, i_atom);
      }
      if (l_ylm_max >= 3) {
        ylm_tab(14, i_atom) = -ylm_tab(14, i_atom);
        dylm_dtheta_tab(14, i_atom) = -dylm_dtheta_tab(14, i_atom);
        scaled_dylm_dphi_tab(14, i_atom) = -scaled_dylm_dphi_tab(14, i_atom);
        ylm_tab(16, i_atom) = -ylm_tab(16, i_atom);
        dylm_dtheta_tab(16, i_atom) = -dylm_dtheta_tab(16, i_atom);
        scaled_dylm_dphi_tab(16, i_atom) = -scaled_dylm_dphi_tab(16, i_atom);
      }
      if (l_ylm_max >= 4) {
        ylm_tab(22, i_atom) = -ylm_tab(22, i_atom);
        dylm_dtheta_tab(22, i_atom) = -dylm_dtheta_tab(22, i_atom);
        scaled_dylm_dphi_tab(22, i_atom) = -scaled_dylm_dphi_tab(22, i_atom);
        ylm_tab(24, i_atom) = -ylm_tab(24, i_atom);
        dylm_dtheta_tab(24, i_atom) = -dylm_dtheta_tab(24, i_atom);
        scaled_dylm_dphi_tab(24, i_atom) = -scaled_dylm_dphi_tab(24, i_atom);
      }
      if (l_ylm_max >= 5) {
        ylm_tab(32, i_atom) = -ylm_tab(32, i_atom);
        dylm_dtheta_tab(32, i_atom) = -dylm_dtheta_tab(32, i_atom);
        scaled_dylm_dphi_tab(32, i_atom) = -scaled_dylm_dphi_tab(32, i_atom);
        ylm_tab(34, i_atom) = -ylm_tab(34, i_atom);
        dylm_dtheta_tab(34, i_atom) = -dylm_dtheta_tab(34, i_atom);
        scaled_dylm_dphi_tab(34, i_atom) = -scaled_dylm_dphi_tab(34, i_atom);
        ylm_tab(36, i_atom) = -ylm_tab(36, i_atom);
        dylm_dtheta_tab(36, i_atom) = -dylm_dtheta_tab(36, i_atom);
        scaled_dylm_dphi_tab(36, i_atom) = -scaled_dylm_dphi_tab(36, i_atom);
      }
      if (l_ylm_max >= 6) {
        ylm_tab(44, i_atom) = -ylm_tab(44, i_atom);
        dylm_dtheta_tab(44, i_atom) = -dylm_dtheta_tab(44, i_atom);
        scaled_dylm_dphi_tab(44, i_atom) = -scaled_dylm_dphi_tab(44, i_atom);
        ylm_tab(46, i_atom) = -ylm_tab(46, i_atom);
        dylm_dtheta_tab(46, i_atom) = -dylm_dtheta_tab(46, i_atom);
        scaled_dylm_dphi_tab(46, i_atom) = -scaled_dylm_dphi_tab(46, i_atom);
        ylm_tab(48, i_atom) = -ylm_tab(48, i_atom);
        dylm_dtheta_tab(48, i_atom) = -dylm_dtheta_tab(48, i_atom);
        scaled_dylm_dphi_tab(48, i_atom) = -scaled_dylm_dphi_tab(48, i_atom);
      }
    }
  }
#undef trigonom_tab
#undef ylm_tab
#undef dylm_dtheta_tab
#undef scaled_dylm_dphi_tab

#undef species_center
#undef dir_tab
}

double spline_vector_waves_c_(double r_output, global double *spl_param, int n_grid_dim, int n_compute_fns, int spline_start,
                            int spline_end, int n_spl_points, int n_spline, global double *out_wave, int index) {
#define spl_param(i, j, k) spl_param[i - 1 + n_compute_fns * (j - 1 + 4 * (k - 1))]
#define out_wave(i) out_wave[i - 1]
  int i_spl;
  double t, term;
  int i_term;
  i_spl = (int)(r_output);
  i_spl = i_spl > 1 ? i_spl : 1;
  i_spl = i_spl < (n_spl_points - 1) ? i_spl : (n_spl_points - 1);
  // i_spl = fmax(1, i_spl);
  // i_spl = fmin(n_spl_points - 1, i_spl);
  t = r_output - (double)(i_spl);
  double ans = spl_param(index + 1 + spline_start - 1, 1, i_spl);
  // for (int i = 1; i <= n_spline; i++) {
  //   out_wave(i) = spl_param(i + spline_start - 1, 1, i_spl);
  // }
  term = 1.0;
  for (i_term = 2; i_term <= 4; i_term++) {
    term = term * t;
    // for (int i = 1; i <= n_spline; i++) {
    //   out_wave(i) = out_wave(i) + term * spl_param(i + spline_start - 1, i_term, i_spl);
    // }
    ans += term * spl_param(index + 1 + spline_start - 1, i_term, i_spl);
  }
  return ans;
#undef spl_param
#undef out_wave
}

// TODO 注意 spline_data 切片，n_max_spline == 4 ?
void evaluate_radial_functions_p0_c_(global int *spline_array_start, global int *spline_array_end, int *n_compute_atoms_,
                                     int *n_compute_fns_, global double *dist_tab, global double *i_r, global int *atom_index,
                                     global int *i_basis_fns_inv, // (n_basis_fns,n_centers)
                                     global double *spline_data,  // (n_basis_fns,n_max_spline, n_max_grid)
                                     global double *wave_aux, int *derivative_, int *n_compute_,
                                     int *n_basis_list_
                                     // outer
                                     ,
                                     int n_basis_fns, int n_max_grid, mconstant int *species_center,global  int *n_grid,
                                     mconstant int *perm_basis_fns_spl
                                     // tmp
                                     , global double* spline_array_aux
                                     ) {
#define species_center(i) species_center[(i)-1]
#define n_grid(i) n_grid[(i)-1]
#define perm_basis_fns_spl(i) perm_basis_fns_spl[(i)-1]

// #define i_basis_fns_inv(i, j) i_basis_fns_inv[(i)-1 + n_basis_fns * ((j)-1)]
#define i_basis_fns_inv(i, j) i_basis_fns_inv[((i)-1 + n_basis_fns * ((j)-1)) * gsize + gid]
  int gid = get_global_id(0);
  int gsize = get_global_size(0);

  int n_compute = *n_compute_;
  int n_compute_atoms = *n_compute_atoms_;
  int n_compute_fns = *n_compute_fns_;
  int n_basis_list = *n_basis_list_;
  int derivative = *derivative_;
  // double spline_array_aux[n_basis_fns];
  for (int i_atom_1 = 1; i_atom_1 <= n_compute_atoms; i_atom_1++) {
    int current_atom = atom_index[(i_atom_1 - 1)];
    int current_species = species_center(current_atom);
    int spline_start = spline_array_start[i_atom_1 - 1];
    int spline_end = spline_array_end[i_atom_1 - 1];
    int n_spline = spline_end - spline_start + 1;
    double r_point = i_r[i_atom_1 - 1];
    double distance_from_atom = dist_tab[i_atom_1 - 1];
    // spline_vector_waves
    int i_rad = (spline_start - 1);
    for (int i_spline = spline_start; i_spline <= spline_end; i_spline++) {
      i_rad = i_rad + 1;
      int spline_counter = i_spline - (spline_start - 1);
      int current_basis_fn = perm_basis_fns_spl(i_rad);
      // int current_basis_fn_comp = i_basis_fns_inv(current_basis_fn, current_atom);
      int current_basis_fn_comp = i_basis_fns_inv(current_basis_fn, i_atom_1);
      if (current_basis_fn_comp == 0)
        continue;
      // TODO 注意 spline_data 切片
      double tmp = spline_vector_waves_c_(r_point, spline_data, n_max_grid, n_basis_fns, spline_start, spline_end,
                             n_grid(current_species), n_spline, spline_array_aux, spline_counter - 1);
      // if (derivative) {
      //   if (distance_from_atom > 0)
      //     wave_aux[current_basis_fn_comp - 1] = spline_array_aux[spline_counter - 1];
      //   else
      //     wave_aux[current_basis_fn_comp - 1] = 0.0;
      // } else {
      //   wave_aux[current_basis_fn_comp - 1] = spline_array_aux[spline_counter - 1];
      // }
      // if (derivative) {
      //   if (distance_from_atom > 0)
      //     wave_aux[current_basis_fn_comp - 1] = tmp;
      //   else
      //     wave_aux[current_basis_fn_comp - 1] = 0.0;
      // } else {
      //   wave_aux[current_basis_fn_comp - 1] = tmp;
      // }
      wave_aux[current_basis_fn_comp - 1] = (derivative && (distance_from_atom <= 0)) ? 0.0 : tmp;
    }
  }
#undef i_basis_fns_inv

#undef species_center
#undef n_grid
#undef perm_basis_fns_spl
}

// void mul_vec_2_c_; 等价于 mul_vec_c_
void mul_vec_c_(global double *wave, int n_mul, global double *ylm, double factor) {
  for (int i = 0; i < n_mul; i++)
    wave[i] = ylm[i] * factor;
}
inline void mul_vec_c_2(global double *wave, int n_mul, global double *ylm, double factor, int array_factor, int ylm_factor) {
  for (int i = 0; i < n_mul; i++)
    wave[i*array_factor] = ylm[i*ylm_factor] * factor;
}

void evaluate_waves_p2_c_(int *n_compute_, int *n_compute_atoms_, int *n_compute_fns_, int *l_ylm_max_,
                          global double *ylm_tab, // ((l_ylm_max+1)**2, n_compute_atoms )
                          private double *one_over_dist_tab, global double *radial_wave, global double *wave, private int *rad_index, global int *wave_index,
                          global int *l_index, global int *l_count, global int *fn_atom, int *n_zero_compute_, global int *zero_index_point
                          // tmp
                          , global double *aux_radial
                          ) {
  int n_compute = *n_compute_;
  int n_compute_atoms = *n_compute_atoms_;
  int n_compute_fns = *n_compute_fns_;
  int l_ylm_max = *l_ylm_max_;
  int n_zero_compute = *n_zero_compute_;
  // double aux_radial[n_compute_fns];
  int index_start = 1;
  int index_end;
  int ylm_tab_dim1 = (l_ylm_max + 1) * (l_ylm_max + 1);
  for (int i_compute_atom = 1; i_compute_atom <= n_compute_atoms; i_compute_atom++) {
    index_end = rad_index[i_compute_atom - 1];
    for (int i = index_start; i <= index_end; i++)
      aux_radial[i - 1] = radial_wave[i - 1] * one_over_dist_tab[i_compute_atom - 1];
    index_start = index_end + 1;
  }
  for (int i_compute_fn = 1; i_compute_fn <= n_compute_fns; i_compute_fn++) {
    int l_aux = l_count[i_compute_fn - 1];
    int l_index_val = l_aux * l_aux + 1;
    int l_count_val = 2 * l_aux;
    mul_vec_c_(&wave[wave_index[i_compute_fn - 1] - 1], l_count_val + 1,
               &ylm_tab[l_index_val - 1 + (fn_atom[i_compute_fn - 1] - 1) * ylm_tab_dim1],
               aux_radial[i_compute_fn - 1]);
  }
  for (int i_compute_point = 1; i_compute_point <= n_zero_compute; i_compute_point++) {
    int i_compute = zero_index_point[i_compute_point - 1];
    wave[i_compute - 1] = 0.0;
  }
}
void evaluate_waves_p2_c_2(int *n_compute_, int *n_compute_atoms_, int *n_compute_fns_, int *l_ylm_max_,
                          global double *ylm_tab, // ((l_ylm_max+1)**2, n_compute_atoms )
                          private double *one_over_dist_tab, global double *radial_wave, global double *wave, private int *rad_index, global int *wave_index,
                          global int *l_index, global int *l_count, global int *fn_atom, int *n_zero_compute_, global int *zero_index_point
                          // tmp
                          , global double *aux_radial, int array_factor
                          ) {
  int gid = get_global_id(0);
  int gsize = get_global_size(0);

  int n_compute = *n_compute_;
  int n_compute_atoms = *n_compute_atoms_;
  int n_compute_fns = *n_compute_fns_;
  int l_ylm_max = *l_ylm_max_;
  int n_zero_compute = *n_zero_compute_;
  // double aux_radial[n_compute_fns];
  int index_start = 1;
  int index_end;
  int ylm_tab_dim1 = (l_ylm_max + 1) * (l_ylm_max + 1);
  for (int i_compute_atom = 1; i_compute_atom <= n_compute_atoms; i_compute_atom++) {
    index_end = rad_index[i_compute_atom - 1];
    for (int i = index_start; i <= index_end; i++)
      aux_radial[i - 1] = radial_wave[i - 1] * one_over_dist_tab[i_compute_atom - 1];
    index_start = index_end + 1;
  }
  for (int i_compute_fn = 1; i_compute_fn <= n_compute_fns; i_compute_fn++) {
    int l_aux = l_count[i_compute_fn - 1];
    int l_index_val = l_aux * l_aux + 1;
    int l_count_val = 2 * l_aux;
    mul_vec_c_2(&wave[(wave_index[i_compute_fn - 1] - 1) * array_factor], l_count_val + 1,
               &ylm_tab[(l_index_val - 1 + (fn_atom[i_compute_fn - 1] - 1) * ylm_tab_dim1) * gsize + gid],
               aux_radial[i_compute_fn - 1], array_factor, gsize);
  }
  for (int i_compute_point = 1; i_compute_point <= n_zero_compute; i_compute_point++) {
    int i_compute = zero_index_point[i_compute_point - 1];
    wave[(i_compute - 1) * array_factor] = 0.0;
  }
}

void evaluate_h_psi_p2_c_(int *n_compute_, int *n_compute_atoms_, int *n_compute_fns_, int *l_ylm_max_, 
                      global double *ylm_tab, global double *one_over_dist_tab, global double *radial_wave, 
                      global double *H_times_psi, double *local_potential_parts_, global double *kinetic_wave, 
                      double *zora_operator_, global int *rad_index, global int *wave_index, 
                      global int *l_index, global int *l_count, global int *fn_atom, 
                      int *n_zero_compute_, global int *zero_index_point,
                      // tmp
                      global double *T_plus_V) {
  int n_compute = *n_compute_;
  int n_compute_atoms = *n_compute_atoms_;
  int n_compute_fns = *n_compute_fns_;
  int l_ylm_max = *l_ylm_max_;
  int n_zero_compute = *n_zero_compute_;
  double local_potential_parts = *local_potential_parts_;
  double zora_operator = *zora_operator_;

  int index_start, index_end;
  // double T_plus_V[n_compute_fns];
  if (1) {
    // if( (flag_rel.eq.REL_none).or.(flag_rel.eq.REL_atomic_zora) &
    //      .or.(flag_rel.eq.REL_own) .or.(flag_rel.eq.REL_q4c) ) then
    index_start = 1;
    for (int i_compute_atom = 1; i_compute_atom <= n_compute_atoms; i_compute_atom++) {
      index_end = rad_index[i_compute_atom - 1];
      for (int i = index_start; i <= index_end; i++) {
        T_plus_V[i - 1] = local_potential_parts * radial_wave[i - 1] + kinetic_wave[i - 1];
      }
      index_start = index_end + 1;
    }
  } else {
    // if(flag_rel.eq.REL_zora.or. flag_rel==REL_KOLNING_HARMON)then
  }

  index_start = 1;
  for (int i_compute_atom = 1; i_compute_atom <= n_compute_atoms; i_compute_atom++) {
    index_end = rad_index[i_compute_atom - 1];
    for (int i = index_start; i <= index_end; i++) {
      T_plus_V[i - 1] += one_over_dist_tab[i_compute_atom - 1];
    }
    index_start = index_end + 1;
  }

  int ylm_tab_dim1 = (l_ylm_max + 1) * (l_ylm_max + 1);
  for (int i_compute_fn = 1; i_compute_fn <= n_compute_fns; i_compute_fn++) {
    mul_vec_c_(&H_times_psi[wave_index[i_compute_fn - 1] - 1], l_count[i_compute_fn - 1] + 1,
               &ylm_tab[l_index[i_compute_fn - 1] - 1 + (fn_atom[i_compute_fn - 1] - 1) * ylm_tab_dim1],
               T_plus_V[i_compute_fn - 1]);
  }
  for (int i_compute_point = 1; i_compute_point <= n_zero_compute; i_compute_point++) {
    H_times_psi[zero_index_point[i_compute_point - 1] - 1] = 0.0;
  }
}

void prune_density_matrix_sparse_dielectric_c_(global double *density_matrix_sparse, global double *density_matrix_con,
                                               int *n_compute_,
                                               global int *i_basis_index
                                               // outer
                                               ,
                                               int index_hamiltonian_dim2, int position_in_hamiltonian_dim1,
                                               global int *center_to_cell, global int *Cbasis_to_basis, global int *Cbasis_to_center,
                                               global int *index_hamiltonian, global int *position_in_hamiltonian,
                                               global int *column_index_hamiltonian) {
#define center_to_cell(i) center_to_cell[(i)-1]
#define Cbasis_to_basis(i) Cbasis_to_basis[(i)-1]
#define Cbasis_to_center(i) Cbasis_to_center[(i)-1]
#define index_hamiltonian(i, j, k) index_hamiltonian[(((k)-1) * index_hamiltonian_dim2 + (j)-1) * 2 + (i)-1]
#define position_in_hamiltonian(i, j) position_in_hamiltonian[((i)-1) + ((j)-1) * position_in_hamiltonian_dim1]
#define column_index_hamiltonian(i) column_index_hamiltonian[(i)-1]

  int lid = get_local_id(0);
  int lsize = get_local_size(0);
  barrier(CLK_GLOBAL_MEM_FENCE);
  // density_matrix_con = 0.0d0 // 不做全局赋0，靠后面的条件分支来赋0
  int n_compute = *n_compute_;
  for (int i_compute = 0; i_compute < n_compute; i_compute++) {
    for (int j_compute = lid; j_compute < n_compute; j_compute+=lsize) {
      density_matrix_con[j_compute + i_compute * n_compute] = 0.0;
    }
  }
  for (int i_compute = 0; i_compute < n_compute; i_compute++) {
    int i_basis = i_basis_index[i_compute];
    int i_basis_uc = Cbasis_to_basis(i_basis);
    int i_cell = center_to_cell(Cbasis_to_center(i_basis));

    for (int j_compute = lid; j_compute < n_compute; j_compute+=lsize) {
      int j_basis = i_basis_index[j_compute];
      int j_basis_uc = Cbasis_to_basis(j_basis);
      int j_cell = -1;
      int valid_place = -1;
      if (j_basis_uc <= i_basis_uc) {
        j_cell = center_to_cell(Cbasis_to_center(j_basis));
        for (int i_place = index_hamiltonian(1, position_in_hamiltonian(i_cell, j_cell), i_basis_uc);
             i_place <= index_hamiltonian(2, position_in_hamiltonian(i_cell, j_cell), i_basis_uc); i_place++) {
          valid_place = column_index_hamiltonian(i_place) != j_basis_uc ? valid_place : i_place;
          // if (column_index_hamiltonian(i_place) == j_basis_uc) {
          //   // tmp = density_matrix_sparse[i_place - 1];
          //   density_matrix_con[i_compute + j_compute * n_compute] = density_matrix_sparse[i_place - 1];
          //   density_matrix_con[j_compute + i_compute * n_compute] = density_matrix_sparse[i_place - 1];
          // }
        }
      }
      // 优化版
      if (valid_place != -1) {
        density_matrix_con[i_compute + j_compute * n_compute] = density_matrix_sparse[valid_place - 1];
        density_matrix_con[j_compute + i_compute * n_compute] = density_matrix_sparse[valid_place - 1];
      }
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
  }
#undef center_to_cell
#undef Cbasis_to_basis
#undef Cbasis_to_center
#undef index_hamiltonian
#undef position_in_hamiltonian
#undef column_index_hamiltonian
}

void prune_density_matrix_sparse_polar_reduce_memory_local_index(global double *density_matrix_sparse, global double *density_matrix_con,
                                               int* n_compute_,
                                               global int *ins_idx,
                                               int n_local_matrix_size
                                               ) {
  int n_compute_c = *n_compute_;
  int lid = get_local_id(0);
  int lsize = get_local_size(0);
  // if(lid == 0){
  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  for(int i=lid; i<n_compute_c*n_compute_c; i+=lsize){
    density_matrix_con[i] = 0.0;
  }
  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  for(int i=0; i<n_compute_c; i++){
    int i_off = (ins_idx[i] * (ins_idx[i] -1)) / 2;
    for(int j=lid; j<=i; j+=lsize){
      if(ins_idx[j] + i_off > n_local_matrix_size){
        break;
      } else {
        double tmp = density_matrix_sparse[ins_idx[j] + i_off - 1];
        density_matrix_con[j + i*n_compute_c] = i == j ? tmp : tmp * 2;
        // density_matrix_con[i + j*n_compute_c] = density_matrix_sparse[ins_idx[j] + i_off - 1];
      }
    }
  }
  // }
  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
}

void prune_density_matrix_sparse_polar_reduce_memory(global double *density_matrix_sparse, global double *density_matrix_con,
                                               int* n_compute_, mconstant int* i_basis_index,
                                               // outer
                                               int index_hamiltonian_dim2, global int *index_hamiltonian,
                                               global int *column_index_hamiltonian
                                               ) {
#define index_hamiltonian(i, j, k) index_hamiltonian[(((k)-1) * index_hamiltonian_dim2 + (j)-1) * 2 + (i)-1]
  int n_compute_c = *n_compute_;
  int lid = get_local_id(0);
  int lsize = get_local_size(0);
  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  for(int i=lid; i<n_compute_c*n_compute_c; i+=lsize){
    density_matrix_con[i] = 0.0;
  }
  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  for(int i_compute = 1+lid; i_compute <= n_compute_c; i_compute+=lsize){
    int i_basis = i_basis_index[i_compute-1];
    int i_start = index_hamiltonian(1, 1, i_basis);
    int i_end = index_hamiltonian(2, 1, i_basis);

    if(i_start <= i_end){
      for(int j_compute = 1; j_compute <= i_compute; j_compute++){
        int j_basis = i_basis_index[j_compute-1];
        int i_index_real;
        for(int i_place = i_start; i_place <= i_end; i_place++){
          if(column_index_hamiltonian[i_place-1] == j_basis){
            double tmp = density_matrix_sparse[i_place-1];
            density_matrix_con[(j_compute-1) + (i_compute-1) * n_compute_c] = i_compute == j_compute ? tmp : tmp * 2;
            // density_matrix_con[(i_compute-1) + (j_compute-1) * n_compute_c] = density_matrix_sparse[i_place-1];
            i_index_real = i_place;
            break;
          } else if(column_index_hamiltonian[i_place-1] > j_basis){
            i_index_real = i_place;
            break;
          }
        }
        i_start = i_index_real;
      }
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
#undef index_hamiltonian
}


inline void atomicAdd_g_f(volatile __global double *addr, double val){
  union {
  unsigned long u64;
  double f64;
  } next, expected, current;
  current.f64 = *addr;
  do {
  expected.f64 = current.f64;
  next.f64 = expected.f64 + val;
  current.u64 = atomic_cmpxchg( (volatile __global unsigned long *)addr,
  expected.u64, next.u64);
  } while( current.u64 != expected.u64 );
}

void prune_density_matrix_sparse_polar_reduce_memory_reverse(global double *density_matrix_sparse, global double *density_matrix_con,
                                               int* n_compute_, mconstant int* i_basis_index,
                                               // outer
                                               int index_hamiltonian_dim2, global int *index_hamiltonian,
                                               global int *column_index_hamiltonian
                                               ) {
#define index_hamiltonian(i, j, k) index_hamiltonian[(((k)-1) * index_hamiltonian_dim2 + (j)-1) * 2 + (i)-1]
  int n_compute_c = *n_compute_;
  int lid = get_local_id(0);
  int lsize = get_local_size(0);
  // barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  // for(int i=lid; i<n_compute_c*n_compute_c; i+=lsize){
  //   density_matrix_con[i] = 0.0;
  // }
  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  for(int i_compute = 1+lid; i_compute <= n_compute_c; i_compute+=lsize){
    int i_basis = i_basis_index[i_compute-1];
    int i_start = index_hamiltonian(1, 1, i_basis);
    int i_end = index_hamiltonian(2, 1, i_basis);

    if(i_start <= i_end){
      for(int j_compute = 1; j_compute <= i_compute; j_compute++){
        int j_basis = i_basis_index[j_compute-1];
        int i_index_real;
        for(int i_place = i_start; i_place <= i_end; i_place++){
          if(column_index_hamiltonian[i_place-1] == j_basis){
            // double tmp = density_matrix_sparse[i_place-1];
            // density_matrix_con[(j_compute-1) + (i_compute-1) * n_compute_c] = i_compute == j_compute ? tmp : tmp * 2;
            // density_matrix_con[(i_compute-1) + (j_compute-1) * n_compute_c] = density_matrix_sparse[i_place-1];
            if(j_compute<=i_compute){
              atomicAdd_g_f(&density_matrix_sparse[i_place-1], density_matrix_con[(j_compute-1) + (i_compute-1) * n_compute_c]);
            } else {
              atomicAdd_g_f(&density_matrix_sparse[i_place-1], density_matrix_con[(i_compute-1) + (j_compute-1) * n_compute_c]);
            }
            i_index_real = i_place;
            break;
          } else if(column_index_hamiltonian[i_place-1] > j_basis){
            i_index_real = i_place;
            break;
          }
        }
        i_start = i_index_real;
      }
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
#undef index_hamiltonian
}

void evaluate_first_order_h_polar_reduce_memory_c_(
  global double* first_order_H, 
  global int* n_points_,
  global double* partition_tab,  // (n_points)
  global double* grid_coord,     // (n_points)
  global double* H_times_psi,    // (n_max_compute_ham,n_points,n_spin)
  int* n_compute_c_, 
  mconstant int* i_basis_index,         // (n_compute_c)
  global const double* wave,           // (n_max_compute_ham, n_points)
  global double* gradient_basis_wave,  // (n_max_compute_ham, 3, n_points)
  global double* first_order_rho,      // (n_spin, n_points)
  global double* v_hartree_gradient,   // (n_points)
  // LDA
  global double* dVxc_drho,  // (3,n_points)
  // GGA
  global double* vsigma, global double* v2rho2, global double* v2rhosigma, global double* v2sigma2, 
  global double* gradient_rho, global double* first_order_gradient_rho, 
  // ---
  int* n_matrix_size_,
  // use
  global int *ins_idx,
  int *n_basis_local_, // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
  // outer
  int *n_spin_,
  int *n_max_compute_ham_,
  // tmp
  global double *contract,
  global double *wave_t,
  global double *first_order_H_dense,
  // local double local_contract[H_PT_TILE][H_JC_TILE * H_JC_WORK],
  // local double local_wave[H_PT_TILE][H_IC_TILE * H_IC_WORK]
  local double A_local[TILE_M * WORK_M][TILE_K],
  local double B_local[TILE_K][TILE_N * WORK_N]
  ){
  int n_compute_c = *n_compute_c_;
  int n_points = *n_points_;

  int n_spin = *n_spin_;
  int n_max_compute_ham = *n_max_compute_ham_;

  int lid = get_local_id(0);
  int lsize = get_local_size(0);
  // double contract[n_points * n_compute_c];
  // double wave_t[n_points * n_max_compute_ham];
  // double first_order_H_dense[n_compute_c * n_compute_c * n_spin];

  // local double local_contract[H_PT_TILE][H_JC_TILE * H_JC_WORK];
  // local double local_wave[H_PT_TILE][H_IC_TILE * H_IC_WORK];

// #define wave(i,j) wave[(i) - 1 + n_max_compute_ham * ((j) - 1)]
#define wave(i,j) wave[((i) - 1) * n_points + ((j) - 1)]

#define dVxc_drho(i,j) dVxc_drho[(i) - 1 + 3 * ((j) - 1)]
#define first_order_rho(i,j) first_order_rho[(i) - 1 + n_spin * ((j) - 1)]

#define contract(i,j) contract[(i) - 1 + n_points * ((j) - 1)]
#define wave_t(i,j) wave_t[(i) - 1 + n_points * ((j) - 1)]
#define first_order_H_dense(i,j,k) first_order_H_dense[(i) - 1 + n_compute_c * ((j) - 1 + ((k) - 1) * n_compute_c)]
  barrier(CLK_GLOBAL_MEM_FENCE);
  // not use_gga

  if(n_spin == 1){
    int i_spin = 1;
    // for(int i=lid; i<n_compute_c * n_compute_c * n_spin; i+=lsize)
    //   first_order_H_dense[i] = 0;
    for(int i_point=lid+1; i_point<=n_points; i_point+=lsize){
      grid_coord[i_point-1] = partition_tab[i_point-1] * (-grid_coord[i_point-1]+v_hartree_gradient[i_point-1]
           +dVxc_drho(i_spin,i_point)*first_order_rho(i_spin,i_point));
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
    // double one = 1.0;
    // double zero = 0.0;
    // dgemm("T", "N", &n_compute_c, &n_compute_c, &n_points, &one, contract, &n_points,
    //   wave_t, &n_points, &zero, first_order_H_dense, &n_compute_c);
{
  const int M = n_compute_c;
  const int N = n_compute_c;
  const int K = n_points;
  int lid = get_local_id(0);
  int lsize = get_local_size(0);

  // local m_float_type A_local[TILE_M * WORK_M][TILE_K];
  // local m_float_type B_local[TILE_K][TILE_N * WORK_N];

  for (int m_out = 0; m_out < M; m_out += TILE_M * WORK_M) {
    for (int n_out = 0; n_out < N; n_out += TILE_N * WORK_N) {
      int m_in = lid / TILE_N;
      int n_in = lid % TILE_N;
      m_float_type sum[WORK_M][WORK_N];
      for (int i = 0; i < WORK_M; i++) {
        for (int j = 0; j < WORK_N; j++) {
          sum[i][j] = 0;
        }
      }
      m_float_type B_regs[WORK_N];
      for (int k_out = 0; k_out < K; k_out += TILE_K) {
        {
          {
            int m_in = lid / TILE_K;
            int k_in = lid % TILE_K;
#pragma unroll WORK_M
            for (int i = 0; i < WORK_M; i++) {
              m_float_type val = wave((m_out + m_in + TILE_M * i) + 1, (k_out + k_in) + 1);
              long cond = (m_out + m_in + TILE_M * i) >= M || (k_out + k_in) >= K;
              val = select(val, 0.0, cond);
              A_local[m_in + TILE_M * i][k_in] = val;
            }
          }
          {
            int k_in = lid / TILE_K;
            int n_in = lid % TILE_K;
#pragma unroll WORK_N
            for (int i = 0; i < WORK_N; i++) {
              m_float_type val = grid_coord[(k_out + k_in)] * wave((n_out + n_in + TILE_N * i) + 1, (k_out + k_in) + 1);
              long cond = (n_out + n_in + TILE_N * i) >= N || (k_out + k_in) >= K;
              val = select(val, 0.0, cond);
              B_local[k_in][n_in + TILE_N * i] = val;
            }
          }
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        for (int k = 0; k < TILE_K; k++) {
          for (int j = 0; j < WORK_N; j++) {
            B_regs[j] = B_local[k][n_in + TILE_N * j];
          }
          for (int i = 0; i < WORK_M; i++) {
            m_float_type A_reg = A_local[m_in + TILE_M * i][k];
            for (int j = 0; j < WORK_N; j++) {
              sum[i][j] += A_reg * B_regs[j];
            }
          }
        }

        barrier(CLK_LOCAL_MEM_FENCE);
      }

      for (int i = 0; i < WORK_M; i++) {
        for (int j = 0; j < WORK_N; j++) {
          if ((m_out + m_in + TILE_M * i) < M && (n_out + n_in + TILE_N * j) < N)
            // C_group[(m_out + m_in + TILE_M * i) * N + (n_out + n_in + TILE_N * j)] = sum[i][j];
            first_order_H_dense((n_out + n_in + TILE_N * j) + 1, (m_out + m_in + TILE_M * i) + 1, i_spin) = sum[i][j];
        }
      }
    }
  }
}
    barrier(CLK_GLOBAL_MEM_FENCE);
    if(*n_basis_local_ > 0){
      for(int i=0; i<n_compute_c; i++){
        int i_off = (ins_idx[i] * (ins_idx[i] -1)) / 2;
        for(int j=lid; j<=i; j+=lsize){
          atomicAdd_g_f(&first_order_H[ins_idx[j] + i_off - 1], first_order_H_dense(j+1, i+1, i_spin));
        }
      }
    }
  }
  barrier(CLK_GLOBAL_MEM_FENCE);


#undef wave
#undef dVxc_drho
#undef first_order_rho

#undef contract
#undef wave_t
#undef first_order_H_dense
}

// void evaluate_first_order_rho_dielectric_c_(
//     int *restrict n_points_, int *restrict n_compute_c_,
//     int *restrict i_basis_index,                         // (n_compute_c)
//     double *restrict wave,                               // (n_max_compute_ham, n_points)
//     double *restrict first_order_density_matrix_compute, // (n_compute_c,n_compute_c)
//     double *restrict first_order_rho                     // (n_points)
// ) {
//   int n_points = *n_points_;
//   int n_compute_c = *n_compute_c_;

//   if (n_compute_c == 0) {
//     for (int i = 0; i < n_points; i++) {
//       first_order_rho[i] = 0.0;
//     }
//     return;
//   }

//   // (v1) this is the oringal method: nested loops
//   for (int i_point = 0; i_point < n_points; i_point++) {
//     double tmp = 0.0;
//     for (int j_compute = 0; j_compute < n_compute_c; j_compute++) {
//       for (int i_compute = 0; i_compute < n_compute_c; i_compute++) {
//         tmp += first_order_density_matrix_compute[j_compute * n_compute_c + i_compute] *
//                wave[i_point * n_max_compute_ham + i_compute] * wave[i_point * n_max_compute_ham + j_compute];
//       }
//     }
//     // 虽然 v1 写法上好像是 +=，但是从矩阵乘的版本来看应该为 =
//     first_order_rho[i_point] = tmp;
//   }
// }

#define centers_hartree_potential(i) centers_hartree_potential[(i)-1]
#define center_to_atom(i) center_to_atom[(i)-1]
#define species_center(i) species_center[(i)-1]
#define center_to_cell(i) center_to_cell[(i)-1]
#define centers_basis_integrals centers_basis_integrals
#define Cbasis_to_basis(i) cbasis_to_basis[(i)-1]
#define Cbasis_to_center(i) cbasis_to_center[(i)-1]
#define pbc_lists_coords_center(i, j) pbc_lists_coords_center[((j)-1) * 3 + (i)-1]
#define column_index_hamiltonian(i) column_index_hamiltonian[(i)-1]
#define index_hamiltonian(i, j, k) index_hamiltonian[(((k)-1) * index_hamiltonian_dim2 + (j)-1) * 2 + (i)-1]
#define position_in_hamiltonian(i, j) position_in_hamiltonian[((i)-1) + ((j)-1) * position_in_hamiltonian_dim1]

#define n_grid(i) n_grid[(i)-1]
#define r_grid_min(i) r_grid_min[(i)-1]
#define log_r_grid_inc(i) log_r_grid_inc[(i)-1]

#define perm_basis_fns_spl(i) perm_basis_fns_spl[(i)-1]
#define outer_radius_sq(i) outer_radius_sq[(i)-1]
#define basis_fn(i) basis_fn[(i)-1]
#define basis_l(i) basis_l[(i)-1]
#define atom_radius_sq(i) atom_radius_sq[(i)-1]
#define basis_fn_start_spl(i) basis_fn_start_spl[(i)-1]
#define basis_fn_atom(i, j) basis_fn_atom[(i)-1 + ((j)-1) * n_basis_fns]

#define batches_size_rho(i) batches_size_rho[(i)-1]
#define batches_batch_n_compute_rho(i) batches_batch_n_compute_rho[(i)-1]
// #define batches_batch_i_basis_rho(i, j) batches_batch_i_basis_rho[(i)-1 + n_centers_basis_I * ((j)-1)]
#define batches_batch_i_basis_rho(i, j) batches_batch_i_basis_rho[(i)-1 + n_max_compute_dens * ((j)-1)]
#define batches_points_coords_rho(i, j, k) batches_points_coords_rho[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]

kernel void integrate_first_order_rho_sub_tmp2_(
    int l_ylm_max_,
    int n_local_matrix_size_, // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
    int n_basis_local_,       // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
    int first_order_density_matrix_size_, global int *basis_l_max, global int *n_points_all_batches,
    global int *n_batch_centers_all_batches, global int *batch_center_all_batches,
    global int *batch_point_to_i_full_point,
    global int *ins_idx_all_batches, // 仅适用于使用了 local_index 且实际使用 ins_idx 转换矩阵的版本
    global double *first_order_rho,
    global double
        *first_order_density_matrix_sparse, // first_order_density_matrix 等价于 first_order_density_matrix_sparse
    global double *partition_tab,
    // outer nums
    // dimensions num 13
    int n_centers, int n_centers_integrals, int n_max_compute_fns_ham, int n_basis_fns, int n_centers_basis_I,
    int n_max_grid, int n_max_compute_atoms, int n_max_compute_ham, int n_max_compute_dens, int n_max_batch_size,
    // pbc_lists num
    int index_hamiltonian_dim2, int position_in_hamiltonian_dim1, int position_in_hamiltonian_dim2,
    int column_index_hamiltonian_size,
    // rho batch num 27
    int n_my_batches_work_rho, int n_full_points_work_rho, // !!!!!! 记得给这几个值赋值
    // outer arrays 29
    // pbc_lists
    mconstant int *center_to_atom, mconstant int *species_center, mconstant int *center_to_cell, mconstant int *cbasis_to_basis,
    mconstant int *cbasis_to_center, global int *centers_basis_integrals, global int *index_hamiltonian,
    global int *position_in_hamiltonian, global int *column_index_hamiltonian, global double *pbc_lists_coords_center,
    // grids
    global int *n_grid, global double *r_grid_min, global double *log_r_grid_inc,
    // basis
    mconstant int *perm_basis_fns_spl, mconstant double *outer_radius_sq, mconstant int *basis_fn, mconstant int *basis_l,
    mconstant double *atom_radius_sq, mconstant int *basis_fn_start_spl, mconstant int *basis_fn_atom,
    global double *basis_wave_ordered,
    // rho batch 50
    global int *batches_size_rho, global int *batches_batch_n_compute_rho, mconstant int *batches_batch_i_basis_rho,
    global double *batches_points_coords_rho,
    // tmp 54
    global double *dist_tab_sq__, global double *dist_tab__, global double *dir_tab__, global int *atom_index__, global int *atom_index_inv__,
    global int *i_basis_fns__, global int *i_basis_fns_inv__, global int *i_atom_fns__, global int *spline_array_start__, global int *spline_array_end__,
    global double *one_over_dist_tab__, global int *rad_index__, global int *wave_index__, global int *l_index__, global int *l_count__, global int *fn_atom__,
    global int *zero_index_point__, global double *wave__, global double *first_order_density_matrix_con__, global double *i_r__,
    global double *trigonom_tab__, global double *radial_wave__,
    global double *spline_array_aux__, global double *aux_radial__,
    global double *ylm_tab__, global double* dylm_dtheta_tab__, global double* scaled_dylm_dphi_tab__, int max_n_batch_centers
    // int *i_full_points__,                       // test
    // int *i_my_batch_,                         // test
    // int *n_points,                            // test
    // double *first_order_density_matrix_con__, // test
    // int *i_full_points_DM_rho_                // test
) {
  int gid = get_global_id(0);
  int lid = get_local_id(0);
  int gsize = get_global_size(0);
  int lsize = get_local_size(0);
  // int block_id = gid / lsize;

#define IC_TILE 16
#define JC_TILE 16
#define PT_TILE 16
#define IC_WORK 3
#define PT_WORK 3

  local double local_density[IC_TILE * IC_WORK][JC_TILE];
  local double local_wave[JC_TILE][PT_TILE*PT_WORK];
  local double local_tmp_rho[IC_TILE][PT_TILE*PT_WORK];

  // double dist_tab_sq[n_centers_integrals];
  // double dist_tab[n_centers_integrals];
  // double dir_tab[3 * n_centers_integrals];
#define dist_tab_sq(i) dist_tab_sq[(i)-1]
#define dist_tab(i) dist_tab[(i)-1]
// #define dir_tab(i, j) dir_tab[(i)-1 + 3 * ((j)-1)]
#define wave(i, j) wave[(i)-1 + n_max_compute_ham * ((j)-1)]
#define batch_center_all_batches(i, j) batch_center_all_batches[(i)-1 + max_n_batch_centers * ((j)-1)]
#define batch_point_to_i_full_point(i, j) batch_point_to_i_full_point[(i)-1 + n_max_batch_size * ((j)-1)]
  int l_ylm_max = l_ylm_max_;

  // int atom_index[n_centers_integrals];
  // int atom_index_inv[n_centers];                      // 没用？
  // int i_basis_fns[n_basis_fns * n_centers_integrals]; // 没用?
  // int i_basis_fns_inv[n_basis_fns * n_centers];
  // int i_atom_fns[n_basis_fns * n_centers_integrals]; // 没用？
  // int spline_array_start[n_centers_integrals];
  // int spline_array_end[n_centers_integrals];
  // double one_over_dist_tab[n_max_compute_atoms];
  // int rad_index[n_max_compute_atoms];
  // int wave_index[n_max_compute_fns_ham];
  // int l_index[n_max_compute_fns_ham];
  // int l_count[n_max_compute_fns_ham];
  // int fn_atom[n_max_compute_fns_ham];
  // int zero_index_point[n_max_compute_ham];
  // double wave[n_max_compute_ham * n_max_batch_size];
  // double first_order_density_matrix_con[n_max_compute_dens * n_max_compute_dens];
  global double *dist_tab_sq = dist_tab_sq__ + gid * n_max_compute_atoms;
  global double *dist_tab = dist_tab__ + gid * n_max_compute_atoms;
  // global double *dir_tab = dir_tab__ + gid * 3 * n_centers_integrals;
  global double *dir_tab = dir_tab__ + get_group_id(0) * lsize * 3 * n_max_compute_atoms;
  global int *atom_index = atom_index__ + gid * n_max_compute_atoms;                   // use private instead
  // global int *atom_index_inv = atom_index_inv__ + gid * n_centers;
  // global int *atom_index_inv = atom_index_inv__ + get_group_id(0) * get_local_size(0) * n_centers;
  // global int *i_basis_fns = i_basis_fns__ + gid * n_basis_fns * n_centers_integrals;   // NULL removed
  // global int *i_basis_fns_inv = i_basis_fns_inv__ + gid * n_basis_fns * n_centers;
  // global int *i_atom_fns = i_atom_fns__ + gid * n_basis_fns * n_centers_integrals;     // NULL removed
  global int *spline_array_start = spline_array_start__ + gid * n_max_compute_atoms;   // use private instead
  global int *spline_array_end = spline_array_end__ + gid * n_max_compute_atoms;       // use private instead
  // global double *one_over_dist_tab = one_over_dist_tab__ + gid * n_max_compute_atoms;  // use private instead
  // global int *rad_index = rad_index__ + gid * n_max_compute_atoms;
  // private int atom_index[MACRO_n_centers_integrals];
  // private int spline_array_start[MACRO_n_centers_integrals];
  // private int spline_array_end[MACRO_n_centers_integrals];
  private double one_over_dist_tab[MACRO_n_max_compute_atoms];
  private int rad_index[MACRO_n_max_compute_atoms];

  global int *wave_index = wave_index__ + gid * n_max_compute_fns_ham;
  // global int *l_index = l_index__ + gid * n_max_compute_fns_ham;  // val[i] = l_aux * l_aux + 1, store in l_count[i]=val  // NULL removed
  global int *l_count = l_count__ + gid * n_max_compute_fns_ham;  // val[i] = 2 * l_aux, store in l_count[i]=val
  global int *fn_atom = fn_atom__ + gid * n_max_compute_fns_ham;
  global int *zero_index_point = zero_index_point__ + gid * n_max_compute_ham;
  // global double *wave = wave__ + gid * n_max_compute_ham;
  global double *wave_group = wave__ + get_group_id(0) * ((n_max_batch_size+127)/128*128) * ((n_max_compute_ham+127)/128*128) + 128;
  global double *i_r = i_r__ + gid * n_max_compute_atoms;
  // global double *trigonom_tab = trigonom_tab__ + gid * 4 * n_max_compute_atoms;
  global double *radial_wave = radial_wave__ + gid * n_max_compute_fns_ham;

  global double *spline_array_aux = spline_array_aux__ + gid * n_basis_fns;
  global double *aux_radial = aux_radial__ + gid * n_max_compute_atoms * n_basis_fns; // 有风险, n_max_compute_atoms 是猜的

  global double *ylm_tab = ylm_tab__ + gid * ((l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms);
  global double *dylm_dtheta_tab = dylm_dtheta_tab__ + gid * ((l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms);
  // 暂时这两个没用的用一样的空间
  global double *scaled_dylm_dphi_tab = dylm_dtheta_tab__ + gid * ((l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms);

// #define i_basis_fns_inv(i, j) i_basis_fns_inv[(i)-1 + n_basis_fns * ((j)-1)]

  // int i_my_batch = *i_my_batch_;
  // for (int i_my_batch = 1; i_my_batch <= n_my_batches_work_rho; i_my_batch++) {
  for (int i_my_batch = get_group_id(0)+1; i_my_batch <= n_my_batches_work_rho; i_my_batch+=(get_global_size(0) / get_local_size(0))) {
  // for (int i_my_batch = get_group_id(0)+1; i_my_batch <= 2; i_my_batch+=(get_global_size(0) / get_local_size(0))) {
    // int i_my_batch_max = (*i_my_batch_ + 127) < n_my_batches_work_rho ? (*i_my_batch_ + 127) : n_my_batches_work_rho;
    // for (int i_my_batch = *i_my_batch_; i_my_batch <= i_my_batch_max; i_my_batch++) {
    int n_compute_c = batches_batch_n_compute_rho(i_my_batch);
    global double *first_order_density_matrix_con =
              first_order_density_matrix_con__ + get_group_id(0) * n_max_compute_dens * n_max_compute_dens;
    // prune_density_matrix_sparse_dielectric_c_(first_order_density_matrix_sparse, first_order_density_matrix_con,
    //           &n_compute_c, &batches_batch_i_basis_rho(1, i_my_batch), index_hamiltonian_dim2, position_in_hamiltonian_dim1, 
    //           &center_to_cell(1), &Cbasis_to_basis(1), &Cbasis_to_center(1), &index_hamiltonian(1, 1, 1),
    //           &position_in_hamiltonian(1, 1), &column_index_hamiltonian(1));
    if(n_basis_local_ > 0){
      prune_density_matrix_sparse_polar_reduce_memory_local_index(first_order_density_matrix_sparse, first_order_density_matrix_con,
                &n_compute_c, &ins_idx_all_batches[(i_my_batch-1) * (n_basis_local_)], n_local_matrix_size_);
    } else {
      prune_density_matrix_sparse_polar_reduce_memory(first_order_density_matrix_sparse, first_order_density_matrix_con,
                                               &n_compute_c, &batches_batch_i_basis_rho(1, i_my_batch),
                                               // outer
                                               index_hamiltonian_dim2, index_hamiltonian,
                                               column_index_hamiltonian);
    }
    // if(lid == 0 && i_my_batch <= 2){
    //   printf("%d: %d, %d\n", i_my_batch, n_points_all_batches[i_my_batch - 1], batch_point_to_i_full_point(n_points_all_batches[i_my_batch - 1], i_my_batch));
    // }
    if (n_compute_c > 0) {
      // int i_point = 0;
      // for (int i_index = 1; i_index <= batches_size_rho(i_my_batch); i_index++) {
      for (int i_point_div = 0; i_point_div < ((n_points_all_batches[i_my_batch - 1] + lsize - 1) / lsize); i_point_div++) {
      // for (int i_point = lid+1; i_point <= n_points_all_batches[i_my_batch - 1]; i_point+=lsize) {
        int i_point = i_point_div * lsize + lid + 1;
      if(i_point <= n_points_all_batches[i_my_batch - 1]){
        // *i_full_points = *i_full_points + 1;
        // if (partition_tab[*i_full_points - 1] > 0.0) {
        // i_point = i_point + 1;
        // for (int i = 0; i < n_centers; i++)
        //   for (int j = 0; j < n_basis_fns; j++)
        //     i_basis_fns_inv(j+1, i+1) = 0;
        for(int i=0; i < n_basis_fns * (n_max_compute_atoms+1); i++)
          i_basis_fns_inv__[i * gsize + gid] = 0.0;
        // for(int i=0; i < n_centers; i++){
        //   atom_index_inv[i * lsize + lid] = n_max_compute_atoms + 1;
        // }
        // memset(i_basis_fns_inv, 0, sizeof(int) * n_basis_fns * n_centers);
        // double coords_center[3];
        // coords_center[0] = batches_points_coords_rho(1, i_point, i_my_batch);
        // coords_center[1] = batches_points_coords_rho(2, i_point, i_my_batch);
        // coords_center[2] = batches_points_coords_rho(3, i_point, i_my_batch);
        // // tab_atom_centered_coords_p0
        // for (int i_center_L = 1; i_center_L <= n_centers_integrals; i_center_L++) {
        //   dir_tab(1, i_center_L) =
        //       coords_center[0] - pbc_lists_coords_center(1, centers_basis_integrals[i_center_L - 1]);
        //   dir_tab(2, i_center_L) =
        //       coords_center[1] - pbc_lists_coords_center(2, centers_basis_integrals[i_center_L - 1]);
        //   dir_tab(3, i_center_L) =
        //       coords_center[2] - pbc_lists_coords_center(3, centers_basis_integrals[i_center_L - 1]);
        //   dist_tab_sq(i_center_L) = dir_tab(1, i_center_L) * dir_tab(1, i_center_L) +
        //                             dir_tab(2, i_center_L) * dir_tab(2, i_center_L) +
        //                             dir_tab(3, i_center_L) * dir_tab(3, i_center_L);
        // }
        int n_compute_atoms = 0;
        int n_compute_fns = 0;
        int n_zero_compute;
        prune_radial_basis_p2_c_(&n_max_compute_atoms, &n_max_compute_fns_ham, &dist_tab_sq(1), &dist_tab(1),
                                //  &dir_tab(1, 1), // (3, n_atom_list)
                                 dir_tab, // (3, n_atom_list)
                                 &n_compute_atoms, atom_index, NULL, &n_compute_fns, NULL,
                                 i_basis_fns_inv__, // (n_basis_fns,n_centers)
                                 NULL, spline_array_start, spline_array_end, &n_centers_integrals,
                                 centers_basis_integrals, &n_compute_c, &batches_batch_i_basis_rho(1, i_my_batch),
                                 &n_batch_centers_all_batches[i_my_batch - 1], &batch_center_all_batches(1, i_my_batch),
                                 one_over_dist_tab, rad_index, wave_index, NULL, l_count, fn_atom, &n_zero_compute,
                                 zero_index_point
                                 // outer
                                 ,
                                 n_basis_fns, &center_to_atom(1), &species_center(1), &Cbasis_to_basis(1),
                                 &Cbasis_to_center(1), &perm_basis_fns_spl(1), &outer_radius_sq(1), &basis_fn(1),
                                 &basis_l(1), &atom_radius_sq(1), &basis_fn_start_spl(1), &basis_fn_atom(1, 1),
                                 pbc_lists_coords_center,
                                 batches_points_coords_rho(1, i_point, i_my_batch), 
                                 batches_points_coords_rho(2, i_point, i_my_batch),
                                 batches_points_coords_rho(3, i_point, i_my_batch));
        // double i_r[n_max_compute_atoms];
        tab_local_geometry_p2_c_(&n_compute_atoms, atom_index, &dist_tab(1),
                                 i_r
                                 // outer
                                 ,
                                 &species_center(1), &r_grid_min(1), &log_r_grid_inc(1));
        // double trigonom_tab[4 * n_max_compute_atoms];
        // tab_trigonom_p0_c_(&n_compute_atoms, &dir_tab(1, 1), trigonom_tab);
        // double ylm_tab[(l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms];              //
        // double dylm_dtheta_tab[(l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms];      // 没用
        // double scaled_dylm_dphi_tab[(l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms]; // 没用
        tab_gradient_ylm_p0_c_(NULL, basis_l_max, &l_ylm_max, &n_compute_atoms, atom_index, ylm_tab,
                               dylm_dtheta_tab, scaled_dylm_dphi_tab, dir_tab, &species_center(1));
                              //  dylm_dtheta_tab, scaled_dylm_dphi_tab, &dir_tab(1, 1), &species_center(1));
        int mfalse = 0;
        // double radial_wave[n_max_compute_fns_ham];
        evaluate_radial_functions_p0_c_(
            spline_array_start, spline_array_end, &n_compute_atoms, &n_compute_fns, &dist_tab(1), i_r, atom_index,
            i_basis_fns_inv__, basis_wave_ordered, radial_wave, &mfalse, &n_compute_c,
            &n_max_compute_fns_ham
            // outer
            ,
            n_basis_fns, n_max_grid, &species_center(1), &n_grid(1), &perm_basis_fns_spl(1),
            spline_array_aux);
        evaluate_waves_p2_c_(&n_compute_c, &n_compute_atoms, &n_compute_fns, &l_ylm_max, ylm_tab, one_over_dist_tab,
                             radial_wave, wave_group + n_max_compute_ham * (i_point-1), rad_index, wave_index, NULL, l_count, fn_atom,
                             &n_zero_compute, zero_index_point, aux_radial);
      } // if(i_point <= n_points_all_batches[i_my_batch - 1])
        int point_valid = (i_point <= n_points_all_batches[i_my_batch - 1]);
        int tmp_point = i_point - 1;
        double i_point_rho = 0.0;
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
        // global double* wave_group = &wave__[get_group_id(0) * get_local_size(0) * n_max_compute_ham];
        int x_point_off_max = min(lsize, n_points_all_batches[i_my_batch - 1] - i_point_div * lsize);
        for(int i_compute_off = 0; i_compute_off < n_compute_c; i_compute_off+=(IC_TILE*IC_WORK)){  // 1,2 是对的，怀疑因为波阵列是 32
          int i_compute_max = min(n_compute_c-i_compute_off, (IC_TILE*IC_WORK));
          for(int x_point_off = 0; x_point_off < x_point_off_max; x_point_off+=(PT_TILE*PT_WORK)){
            int x_point_max = min(x_point_off_max-x_point_off, (PT_TILE*PT_WORK));
            double private_rho[IC_WORK][PT_WORK];
            for(int i=0; i<IC_WORK; i++)
              for(int j=0; j<PT_WORK; j++)
                private_rho[i][j] = 0.0;
            for(int x_point_work = 0; x_point_work < PT_WORK; x_point_work++){
              local_tmp_rho[lid / JC_TILE][lid % JC_TILE + PT_TILE * x_point_work] = 0.0;
            }
            for(int j_compute_off = 0; j_compute_off < (i_compute_off+i_compute_max); j_compute_off+=JC_TILE){
              int j_compute_max = min((i_compute_off+i_compute_max)-j_compute_off, JC_TILE);
              // int j_compute_max = min((n_compute_c+i_compute_max)-j_compute_off, JC_TILE);
              int id = lid;
              int i_compute = id / JC_TILE;
              int x_point = id % JC_TILE;
              {
                int i_compute = id / JC_TILE;
                int x_point = id / JC_TILE;
                int j_compute = id % JC_TILE;
                for(int x_point_work = 0; x_point_work < PT_WORK; x_point_work++){
                  // 有一点点越界风险
                  local_wave[j_compute][x_point + PT_TILE * x_point_work] = j_compute < j_compute_max && (x_point + PT_TILE * x_point_work) < x_point_max
                    ? wave_group[(x_point_off + x_point + PT_TILE * x_point_work) * n_max_compute_ham + (j_compute_off + j_compute)]
                    : 0.0; 
                }
                for(int i_compute_work = 0; i_compute_work < IC_WORK; i_compute_work++){
                  if(i_compute + IC_TILE*i_compute_work < i_compute_max){
                    // double tmp;
                    // if((j_compute_off + j_compute) < (i_compute_off + i_compute + IC_TILE*i_compute_work))
                    //   tmp = first_order_density_matrix_con[(i_compute_off + i_compute + IC_TILE*i_compute_work) * n_compute_c + (j_compute_off + j_compute)];
                    // else
                    //   tmp = 0.0;
                    // local_density[i_compute + IC_TILE*i_compute_work][j_compute] = tmp;
                    local_density[i_compute + IC_TILE*i_compute_work][j_compute] = 
                      first_order_density_matrix_con[(i_compute_off + i_compute + IC_TILE*i_compute_work) * n_compute_c + (j_compute_off + j_compute)];
                  }
                }
              }
              barrier(CLK_LOCAL_MEM_FENCE);
              for(int i_compute_work = 0; i_compute_work < IC_WORK; i_compute_work++){
                if((i_compute+IC_TILE * i_compute_work) < i_compute_max && x_point < x_point_max){
                  for(int x_point_work = 0; x_point_work < PT_WORK; x_point_work++){
                    if((x_point + PT_TILE * x_point_work) < x_point_max){
                      for(int j_compute = 0; j_compute < JC_TILE; j_compute++){
                        private_rho[i_compute_work][x_point_work] += local_wave[j_compute][x_point+PT_TILE*x_point_work]
                                     * local_density[i_compute+IC_TILE*i_compute_work][j_compute];
                      }
                    }
                  }
                }
              }
              barrier(CLK_LOCAL_MEM_FENCE);
            }
            int id = lid;
            int i_compute = id / JC_TILE;
            int x_point = id % JC_TILE;
            for(int i_compute_work = 0; i_compute_work < IC_WORK; i_compute_work++){
              if((i_compute+IC_TILE * i_compute_work) < i_compute_max && x_point < x_point_max){
                  for(int x_point_work = 0; x_point_work < PT_WORK; x_point_work++){
                    if((x_point + PT_TILE * x_point_work) < x_point_max){
                      private_rho[i_compute_work][x_point_work] *= wave_group[(x_point_off + x_point + PT_TILE * x_point_work) * n_max_compute_ham + i_compute_off + i_compute + IC_TILE * i_compute_work];
                      local_tmp_rho[i_compute][x_point + PT_TILE * x_point_work] += private_rho[i_compute_work][x_point_work];
                    }
                  }
              }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if(x_point_off <= lid && lid < (x_point_off + x_point_max)){
              // i_point_rho += private_rho;
              int x_point = lid % (PT_TILE*PT_WORK);
              int maxi = min(IC_TILE, i_compute_max);
              for(int i=0; i<maxi; i++){
                i_point_rho += local_tmp_rho[i][x_point];
              }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
          }
        }
#undef IC_TILE
#undef JC_TILE
#undef PT_TILE
        if(point_valid){
          first_order_rho[batch_point_to_i_full_point(i_point, i_my_batch) - 1] = i_point_rho;
          }
        // if(511 <= batch_point_to_i_full_point(i_point, i_my_batch) - 1 && batch_point_to_i_full_point(i_point, i_my_batch) - 1 <= 530){
        //   printf("gid=%d, lid=%d, i_point=%d, i_point_div=%d, ifullpoint=%d, npoint=%d, i_my_batch=%d\n", 
        //     gid, lid, i_point, i_point_div, batch_point_to_i_full_point(i_point, i_my_batch) - 1, n_points_all_batches[i_my_batch - 1], i_my_batch);
        // }
#undef WAVEJ_TILE_SIZE
#undef WAVEI_TILE_SIZE
      }
      // evaluate_first_order_rho_dielectric_c_(n_points, &n_compute_c, &batches_batch_i_basis_s(1, i_my_batch), wave,
      //                                        first_order_density_matrix_con, local_first_order_rho);
      // i_point = 0;
      // for (int i_index = 1; i_index <= batches_size_rho(i_my_batch); i_index++) {
      //   *i_full_points_DM_rho_ = *i_full_points_DM_rho_ + 1;
      //   if (partition_tab[*i_full_points_DM_rho_ - 1] > 0.0) {
      //     i_point = i_point + 1;
      //     first_order_rho[batch_point_to_i_full_point(i_point, i_my_batch) - 1] = local_first_order_rho[i_point - 1];
      //   }
      // }
    } else {
      // 既然 first_order_rho 初始化为全 0, 那么理应无需赋 0
      // i_full_points = i_full_points + batches_size_rho(i_my_batch);
      // for (int i = *i_full_points_DM_rho_ + 1; i <= *i_full_points_DM_rho_ + batches_size_rho(i_my_batch); i++) {
      //   first_order_rho[i - 1] = 0.0;
      // }
      // *i_full_points_DM_rho_ = *i_full_points_DM_rho_ + batches_size_rho(i_my_batch);
    }
    // if (i_my_batch == 1) {
    //   m_save_check_rho_wave_(wave, &n_compute_c);
    // }
  }
#undef i_basis_fns_inv

  // m_save_check_rho_(first_order_rho);
#undef dist_tab_sq
#undef dist_tab
#undef dir_tab
#undef wave
#undef batch_center_all_batches
#undef batch_point_to_i_full_point
}


#define batches_size_h(i) batches_size_h[(i)-1]
#define batches_batch_n_compute_h(i) batches_batch_n_compute_h[(i)-1]
#define batches_batch_i_basis_h(i, j) batches_batch_i_basis_h[(i)-1 + n_max_compute_dens * ((j)-1)]
#define batches_points_coords_h(i, j, k) batches_points_coords_h[(((k)-1) * n_max_batch_size + (j)-1) * 3 + (i)-1]

kernel void integrate_first_order_h_sub_tmp2_(
  int j_coord,
  int n_spin,
  int l_ylm_max,
  int n_basis_local,
  int n_matrix_size,
  global int* basis_l_max,
  global int* n_points_all_batches,
  global int* n_batch_centers_all_batches,
  global int* batch_center_all_batches,
  global int* ins_idx_all_batches,
  global int* batches_batch_i_basis_h_not_use__,
  global double* partition_all_batches,
  global double* first_order_H,
  global double* local_potential_parts_all_points,
  global double* local_first_order_rho_all_batches,
  global double* local_first_order_potential_all_batches,
  global double* local_dVxc_drho_all_batches,
  global double* local_rho_gradient,
  global double* first_order_gradient_rho,
  // outer nums
  // dimensions num 19
  int n_centers, int n_centers_integrals, int n_max_compute_fns_ham, int n_basis_fns, int n_centers_basis_I,
  int n_max_grid, int n_max_compute_atoms, int n_max_compute_ham, int n_max_compute_dens, int n_max_batch_size,
  // pbc_lists num
  int index_hamiltonian_dim2, int position_in_hamiltonian_dim1, int position_in_hamiltonian_dim2,
  int column_index_hamiltonian_size,
  // H batch num
  int n_my_batches_work_h, int n_full_points_work_h,
  // outer arrays 35
  // pbc_lists
  mconstant int *center_to_atom, mconstant int *species_center, mconstant int *center_to_cell, mconstant int *cbasis_to_basis,
  mconstant int *cbasis_to_center, global int *centers_basis_integrals, global int *index_hamiltonian,
  global int *position_in_hamiltonian, global int *column_index_hamiltonian, global double *pbc_lists_coords_center,
  // grids
  global int *n_grid, global double *r_grid_min, global double *log_r_grid_inc,
  // basis
  mconstant int *perm_basis_fns_spl, mconstant double *outer_radius_sq, mconstant int *basis_fn, mconstant int *basis_l,
  mconstant double *atom_radius_sq, mconstant int *basis_fn_start_spl, mconstant int *basis_fn_atom,
  global double *basis_wave_ordered,
  global double *basis_kinetic_ordered, // new !!!!
  // H batch
  global int *batches_batch_n_compute_h, mconstant int *batches_batch_i_basis_h,
  global double *batches_points_coords_h,
  // tmp 60
  global double *dist_tab_sq__, global double *dist_tab__, global double *dir_tab__, global int *atom_index__, global int *atom_index_inv__,
  global int *i_basis_fns__, global int *i_basis_fns_inv__, global int *i_atom_fns__, global int *spline_array_start__, global int *spline_array_end__,
  global double *one_over_dist_tab__, global int *rad_index__, global int *wave_index__, global int *l_index__, global int *l_count__, global int *fn_atom__,
  global int *zero_index_point__, global double *wave__, global double *first_order_density_matrix_con__, global double *i_r__,
  global double *trigonom_tab__, global double *radial_wave__,
  global double *spline_array_aux__, global double *aux_radial__,
  global double *ylm_tab__, global double* dylm_dtheta_tab__, global double* scaled_dylm_dphi_tab__,
  // tmp more
  global double *kinetic_wave__, global double *grid_coord__, global double *H_times_psi__, global double *T_plus_V__,
  global double *contract__, global double *wave_t__, global double *first_order_H_dense__, int max_n_batch_centers
  // test
  // int *i_my_batch_                         // test
){

  int gid = get_global_id(0);
  int lid = get_local_id(0);
  int gsize = get_global_size(0);
  int lsize = get_local_size(0);
  // int block_id = gid / lsize;
  // local double local_contract[H_PT_TILE][H_JC_TILE * H_JC_WORK];
  // local double local_wave[H_PT_TILE][H_IC_TILE * H_IC_WORK];
  local double A_local[TILE_M * WORK_M][TILE_K];
  local double B_local[TILE_K][TILE_N * WORK_N];

#define dist_tab_sq(i) dist_tab_sq[(i)-1]
#define dist_tab(i) dist_tab[(i)-1]
// #define dir_tab(i, j) dir_tab[(i)-1 + 3 * ((j)-1)]
// #define dir_tab(i, j) dir_tab[lid + lsize * ((i)-1 + 3 * ((j)-1))]
#define wave(i, j) wave[(i)-1 + n_max_compute_ham * ((j)-1)]
#define batch_center_all_batches(i, j) batch_center_all_batches[(i)-1 + max_n_batch_centers * ((j)-1)]
#define batch_point_to_i_full_point(i, j) batch_point_to_i_full_point[(i)-1 + n_max_batch_size * ((j)-1)]
  // int l_ylm_max = l_ylm_max_;

  // global double *wave = wave__ + gid * n_max_compute_ham;
  global double *wave_group = wave__ + get_group_id(0) * n_max_batch_size * n_max_compute_ham;
  global double *i_r = i_r__ + gid * n_max_compute_atoms;
  global double *trigonom_tab = trigonom_tab__ + gid * 4 * n_max_compute_atoms;
  global double *radial_wave = radial_wave__ + gid * n_max_compute_fns_ham;

  // global double *spline_array_aux = spline_array_aux__ + gid * n_basis_fns;
  global double *aux_radial = aux_radial__ + gid * n_max_compute_atoms * n_basis_fns; // 有风险, n_max_compute_atoms 是猜的

  // global double *ylm_tab = ylm_tab__ + gid * ((l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms);
  global double *dylm_dtheta_tab = dylm_dtheta_tab__ + gid * ((l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms);
  // 暂时这两个没用的用一样的空间
  global double *scaled_dylm_dphi_tab = dylm_dtheta_tab__ + gid * ((l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms);

  // global double *kinetic_wave = kinetic_wave__ + gid * n_max_compute_fns_ham;
  global double *grid_coord_group = grid_coord__ + get_group_id(0) * n_max_batch_size; // 共用性类似 wave
  // global double *H_times_psi_group = H_times_psi__ + get_group_id(0) * (n_max_compute_ham * lsize * n_spin);   // 共用性类似 wave
  // global double *T_plus_V = T_plus_V__ + gid * n_max_compute_atoms * n_basis_fns; // 有风险, n_max_compute_atoms 是猜的

  // global double *contract_group = contract__ + get_group_id(0) * n_max_batch_size * n_max_compute_ham;
  // global double *wave_t_group = wave_t__ + get_group_id(0) * n_max_batch_size * n_max_compute_ham;
  global double *first_order_H_dense_group = first_order_H_dense__ + get_group_id(0) * n_max_compute_ham * n_max_compute_ham * n_spin;

// #define i_basis_fns_inv(i, j) i_basis_fns_inv[(i)-1 + n_basis_fns * ((j)-1))]
  // return;
  // int i_my_batch = *i_my_batch_;
  // for (int i_my_batch = 1; i_my_batch <= n_my_batches_work_h; i_my_batch++) {
  for (int i_my_batch = get_group_id(0)+1; i_my_batch <= n_my_batches_work_h; i_my_batch+=(get_global_size(0) / get_local_size(0))) {
    int n_compute_c = batches_batch_n_compute_h(i_my_batch);
    if (n_compute_c > 0) {
      for (int i_point_div = 0; i_point_div < ((n_points_all_batches[i_my_batch - 1] + lsize - 1) / lsize); i_point_div++) {
        int i_point = i_point_div * lsize + lid + 1;
      if(i_point <= n_points_all_batches[i_my_batch - 1]){

        global double *dist_tab_sq = dist_tab_sq__ + gid * n_max_compute_atoms;
        global double *dist_tab = dist_tab__ + gid * n_max_compute_atoms;
        global double *dir_tab = dir_tab__ + get_group_id(0) * lsize * 3 * n_max_compute_atoms;
        global int *atom_index = atom_index__ + gid * n_max_compute_atoms;                   // use private instead
        // global int *atom_index_inv = atom_index_inv__ + gid * n_centers;
        // global int *atom_index_inv = atom_index_inv__ + get_group_id(0) * get_local_size(0) * n_centers;
        // global int *i_basis_fns = i_basis_fns__ + gid * n_basis_fns * n_centers_integrals;   // NULL removed
        // global int *i_basis_fns_inv = i_basis_fns_inv__ + gid * n_basis_fns * n_centers;
        // global int *i_atom_fns = i_atom_fns__ + gid * n_basis_fns * n_centers_integrals;     // NULL removed
        global int *spline_array_start = spline_array_start__ + gid * n_max_compute_atoms;   // use private instead
        global int *spline_array_end = spline_array_end__ + gid * n_max_compute_atoms;       // use private instead
        // global double *one_over_dist_tab = one_over_dist_tab__ + gid * n_max_compute_atoms;  // use private instead
        // global int *rad_index = rad_index__ + gid * n_max_compute_atoms;                     // use private instead
        global int *wave_index = wave_index__ + gid * n_max_compute_fns_ham;
        // global int *l_index = l_index__ + gid * n_max_compute_fns_ham;  // val[i] = l_aux * l_aux + 1, store in l_count[i]=val  // NULL removed
        global int *l_count = l_count__ + gid * n_max_compute_fns_ham;  // val[i] = 2 * l_aux, store in l_count[i]=val
        global int *fn_atom = fn_atom__ + gid * n_max_compute_fns_ham;
        global int *zero_index_point = zero_index_point__ + gid * n_max_compute_ham;
        
        // private int atom_index[MACRO_n_centers_integrals];
        // private int spline_array_start[MACRO_n_centers_integrals];
        // private int spline_array_end[MACRO_n_centers_integrals];
        private int rad_index[MACRO_n_max_compute_atoms];
        private double one_over_dist_tab[MACRO_n_max_compute_atoms];

        for(int i=0; i < n_basis_fns * (n_max_compute_atoms+1); i++)
          i_basis_fns_inv__[i * gsize + gid] = 0.0;

        // for(int i=0; i < n_centers; i++){
        //   atom_index_inv[i * lsize + lid] = n_max_compute_atoms + 1;
        // }

        double zora_operator[2]; // double zora_operator[n_spin]; n_spin <= 2
        // double coords_center[3];
        // coords_center[0] = batches_points_coords_h(1, i_point, i_my_batch);
        // coords_center[1] = batches_points_coords_h(2, i_point, i_my_batch);
        // coords_center[2] = batches_points_coords_h(3, i_point, i_my_batch);
        // tab_atom_centered_coords_p0
        // for (int i_center_L = 1; i_center_L <= n_centers_integrals; i_center_L++) {
        //   dir_tab(1, i_center_L) =
        //       coords_center[0] - pbc_lists_coords_center(1, centers_basis_integrals[i_center_L - 1]);
        //   dir_tab(2, i_center_L) =
        //       coords_center[1] - pbc_lists_coords_center(2, centers_basis_integrals[i_center_L - 1]);
        //   dir_tab(3, i_center_L) =
        //       coords_center[2] - pbc_lists_coords_center(3, centers_basis_integrals[i_center_L - 1]);
        //   dist_tab_sq(i_center_L) = dir_tab(1, i_center_L) * dir_tab(1, i_center_L) +
        //                             dir_tab(2, i_center_L) * dir_tab(2, i_center_L) +
        //                             dir_tab(3, i_center_L) * dir_tab(3, i_center_L);
        // }
        int n_compute_atoms = 0;
        int n_compute_fns = 0;
        int n_zero_compute;

        // grid_coord__[gid] = batches_points_coords_h(j_coord, i_point, i_my_batch);
        grid_coord_group[i_point - 1] = batches_points_coords_h(j_coord, i_point, i_my_batch);

        prune_radial_basis_p2_c_(&n_max_compute_atoms, &n_max_compute_fns_ham, &dist_tab_sq(1), &dist_tab(1),
                                 dir_tab, // (3, n_atom_list)
                                 &n_compute_atoms, atom_index, NULL, &n_compute_fns, NULL,
                                //  &n_compute_atoms, atom_index, atom_index_inv, &n_compute_fns, NULL,
                                 i_basis_fns_inv__, // (n_basis_fns,n_centers)
                                 NULL, spline_array_start, spline_array_end, &n_centers_integrals,
                                 centers_basis_integrals, &n_compute_c, &batches_batch_i_basis_h(1, i_my_batch),
                                 &n_batch_centers_all_batches[i_my_batch - 1], &batch_center_all_batches(1, i_my_batch),
                                 one_over_dist_tab, rad_index, wave_index, NULL, l_count, fn_atom, &n_zero_compute,
                                 zero_index_point
                                 // outer
                                 ,
                                 n_basis_fns, &center_to_atom(1), &species_center(1), &Cbasis_to_basis(1),
                                 &Cbasis_to_center(1), &perm_basis_fns_spl(1), &outer_radius_sq(1), &basis_fn(1),
                                 &basis_l(1), &atom_radius_sq(1), &basis_fn_start_spl(1), &basis_fn_atom(1, 1),
                                 pbc_lists_coords_center,
                                 batches_points_coords_h(1, i_point, i_my_batch), 
                                 batches_points_coords_h(2, i_point, i_my_batch),
                                 batches_points_coords_h(3, i_point, i_my_batch));
        // double i_r[n_max_compute_atoms];
        tab_local_geometry_p2_c_(&n_compute_atoms, atom_index, &dist_tab(1),
                                 i_r
                                 // outer
                                 ,
                                 &species_center(1), &r_grid_min(1), &log_r_grid_inc(1));
        // double trigonom_tab[4 * n_max_compute_atoms];
        // tab_trigonom_p0_c_(&n_compute_atoms, &dir_tab(1, 1), trigonom_tab);
        // double ylm_tab[(l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms];              //
        // double dylm_dtheta_tab[(l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms];      // 没用
        // double scaled_dylm_dphi_tab[(l_ylm_max + 1) * (l_ylm_max + 1) * n_max_compute_atoms]; // 没用
        tab_gradient_ylm_p0_c_2(NULL, basis_l_max, &l_ylm_max, &n_compute_atoms, atom_index, ylm_tab__,
                               dylm_dtheta_tab, scaled_dylm_dphi_tab, dir_tab, &species_center(1));
        int mfalse = 0;
        // double radial_wave[n_max_compute_fns_ham];
        evaluate_radial_functions_p0_c_(
            spline_array_start, spline_array_end, &n_compute_atoms, &n_compute_fns, &dist_tab(1), i_r, atom_index,
            i_basis_fns_inv__, basis_wave_ordered, radial_wave, &mfalse, &n_compute_c,
            &n_max_compute_fns_ham
            // outer
            ,
            n_basis_fns, n_max_grid, &species_center(1), &n_grid(1), &perm_basis_fns_spl(1),
            NULL);

        evaluate_waves_p2_c_2(&n_compute_c, &n_compute_atoms, &n_compute_fns, &l_ylm_max, ylm_tab__, one_over_dist_tab,
                             radial_wave, &wave_group[(i_point-1)], rad_index, wave_index, NULL, l_count, fn_atom,
                             &n_zero_compute, zero_index_point, aux_radial, n_points_all_batches[i_my_batch - 1]);
        // printf("gid=%6d lid=%3d i_my_batch=%4d i_point_div=%3d i_point=%3d n_point=%4d\n", 
        //   gid, lid, i_my_batch, i_point_div, i_point, n_points_all_batches[i_my_batch - 1]);
        // evaluate_radial_functions_p0_c_(
        //     spline_array_start, spline_array_end, &n_compute_atoms, &n_compute_fns, &dist_tab(1), i_r, atom_index,
        //     i_basis_fns_inv, basis_kinetic_ordered, kinetic_wave, &mfalse, &n_compute_c,
        //     &n_max_compute_fns_ham
        //     // outer
        //     ,
        //     n_basis_fns, n_max_grid, &species_center(1), &n_grid(1), &perm_basis_fns_spl(1),
        //     spline_array_aux);
        // for(int i_spin = 1; i_spin <= n_spin; i_spin++){
        //   evaluate_h_psi_p2_c_(&n_compute_c, &n_compute_atoms, &n_compute_fns, &l_ylm_max, ylm_tab, one_over_dist_tab,
        //                       radial_wave, &H_times_psi_group[n_max_compute_ham * (lid + (i_spin-1) * lsize)],
        //                       &local_potential_parts_all_points[i_spin-1 + n_spin * (i_point - 1)],
        //                       kinetic_wave, &zora_operator[i_spin],
        //                       rad_index, wave_index, l_index, l_count, fn_atom, &n_zero_compute, zero_index_point,
        //                       T_plus_V);
        // }
      } // if(i_point <= n_points_all_batches[i_my_batch - 1])
      }
      global double* H_times_psi_group = NULL;
      evaluate_first_order_h_polar_reduce_memory_c_(first_order_H, &n_points_all_batches[i_my_batch-1],
              &partition_all_batches[(i_my_batch-1) * n_max_batch_size], grid_coord_group,
              H_times_psi_group, &n_compute_c, &batches_batch_i_basis_h(1, i_my_batch), 
              wave_group, NULL,
              &local_first_order_rho_all_batches[n_spin * n_max_batch_size * (i_my_batch-1)],
              &local_first_order_potential_all_batches[n_max_batch_size * (i_my_batch-1)],
              &local_dVxc_drho_all_batches[3 * n_max_batch_size * (i_my_batch-1)],
              NULL, NULL, NULL, NULL,
              local_rho_gradient,
              first_order_gradient_rho, &n_matrix_size,
              &ins_idx_all_batches[n_basis_local * (i_my_batch-1)], &n_basis_local, &n_spin, &n_max_compute_ham,
              NULL, NULL, first_order_H_dense_group, A_local, B_local);
      if(n_basis_local <= 0){
        prune_density_matrix_sparse_polar_reduce_memory_reverse(first_order_H, first_order_H_dense_group,
                                                 &n_compute_c, &batches_batch_i_basis_h(1, i_my_batch),
                                                 // outer
                                                 index_hamiltonian_dim2, index_hamiltonian,
                                                 column_index_hamiltonian);
      }
    }
  }
#undef i_basis_fns_inv

  // m_save_check_h_(first_order_h);
#undef dist_tab_sq
#undef dist_tab
#undef dir_tab
#undef wave
#undef batch_center_all_batches
#undef batch_point_to_i_full_point
}

#undef H_IC_TILE
#undef H_JC_TILE
#undef H_PT_TILE
#undef H_IC_WORK
#undef H_PT_WORK