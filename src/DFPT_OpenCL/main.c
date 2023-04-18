#include "opencl_util_c.h"
#include "pass_mod_var.h"
#include "save_load_var.h"
#include "sum_up_whole_potential_shanghui.h"
#include <stdio.h>
#include <sys/time.h>

// #define _MY_MPI_DEBUG

#ifdef _MY_MPI_DEBUG
#include <mpi.h>
#endif

int MV(mpi_tasks, n_tasks) = 1;
int MV(mpi_tasks, myid) = 0;
int MV(opencl_util, mpi_platform_relative_id) = 0;

#define n_tasks MV(mpi_tasks, n_tasks)
#define myid MV(mpi_tasks, myid)
#define mpi_platform_relative_id MV(opencl_util, mpi_platform_relative_id)

int MV(dimensions, n_centers_hartree_potential);
int MV(dimensions, n_periodic);
int MV(dimensions, n_max_radial);
int MV(dimensions, l_pot_max);
int MV(dimensions, n_max_spline);
int MV(dimensions, n_hartree_grid);
int MV(dimensions, n_species);
int MV(dimensions, n_atoms);
int MV(dimensions, n_centers);
int MV(dimensions, n_centers_basis_integrals); // rho
int MV(dimensions, n_centers_integrals);       // rho
int MV(dimensions, n_max_compute_fns_ham);     // rho
int MV(dimensions, n_basis_fns);               // rho
int MV(dimensions, n_basis);                   // rho
int MV(dimensions, n_centers_basis_t);         // rho
int MV(dimensions, n_centers_basis_i);         // rho
int MV(dimensions, n_max_grid);                // rho
int MV(dimensions, n_max_compute_atoms);       // rho
int MV(dimensions, n_max_compute_ham);         // rho
int MV(dimensions, n_max_compute_dens);        // rho
int MV(dimensions, n_max_batch_size);

int MV(runtime_choices, use_hartree_non_periodic_ewald);
int MV(runtime_choices, hartree_fp_function_splines);
int MV(runtime_choices, fast_ylm);
int MV(runtime_choices, new_ylm);
int MV(runtime_choices, flag_rel); // rho
int MV(runtime_choices, adams_moulton_integrator);
int MV(runtime_choices, compensate_multipole_errors);

int *MV(geometry, species); // (n_atoms)
int *MV(geometry, empty);   // (n_atoms)

// int MV(pbc_lists, n_cells_in_hamiltonian);
int MV(pbc_lists, index_hamiltonian_dim2);
int MV(pbc_lists, position_in_hamiltonian_dim1);
int MV(pbc_lists, position_in_hamiltonian_dim2);
int MV(pbc_lists, column_index_hamiltonian_size);
int *MV(pbc_lists, centers_hartree_potential); // (n_centers_hartree_potential)
int *MV(pbc_lists, center_to_atom);            // (n_centers)
int *MV(pbc_lists, species_center);            // (n_centers)
int *MV(pbc_lists, center_to_cell);            // (n_centers)
int *MV(pbc_lists, cbasis_to_basis);           // (n_centers_basis_T) // rho
int *MV(pbc_lists, cbasis_to_center);          // (n_centers_basis_T) // rho
int *MV(pbc_lists, centers_basis_integrals);   // (n_centers_basis_integrals) // rho
int *MV(pbc_lists, index_hamiltonian);         // (2, n_cells_in_hamiltonian, n_basis)
int *MV(pbc_lists, position_in_hamiltonian);   // (.._dim1, .._dim2)
int *MV(pbc_lists, column_index_hamiltonian);  // (..size)
double *MV(pbc_lists, coords_center);          // (3,n_centers)

int *MV(species_data, l_hartree); // (n_species)
double *MV(species_data, multipole_radius_free); // (n_species)

int *MV(grids, n_grid);            // (n_species) // rho
int *MV(grids, n_radial);          // (n_species)
double *MV(grids, r_grid_min);     // (n_species)
double *MV(grids, r_grid_inc);     // (n_species)
double *MV(grids, log_r_grid_inc); // (n_species)
double *MV(grids, scale_radial);   // (n_species)
double *MV(grids, r_radial);       // (n_max_radial, n_species)
double *MV(grids, r_grid);         // (n_max_grid, n_species)

int MV(analytic_multipole_coefficients, l_max_analytic_multipole);
int *MV(analytic_multipole_coefficients, n_cc_lm_ijk);      // (0:l_max_analytic_multipole)
int *MV(analytic_multipole_coefficients, index_cc);         // (n_cc_lm_ijk(l_max_analytic_multipole),6)
int *MV(analytic_multipole_coefficients, index_ijk_max_cc); // (3,0:l_max_analytic_multipole)

int MV(hartree_potential_real_p0, n_hartree_atoms);
int MV(hartree_potential_real_p0, hartree_force_l_add);
// double *MV(hartree_potential_real_p0, multipole_c); // ( n_cc_lm_ijk(l_pot_max), n_atoms) // 每次 sumup 初始化
double MV(hartree_potential_real_p0, b0)[pmaxab + 1];     // (0:pmaxab)
double MV(hartree_potential_real_p0, b2)[pmaxab + 1];     // (0:pmaxab)
double MV(hartree_potential_real_p0, b4)[pmaxab + 1];     // (0:pmaxab)
double MV(hartree_potential_real_p0, b6)[pmaxab + 1];     // (0:pmaxab)
double MV(hartree_potential_real_p0, a_save)[pmaxab + 1]; // (0:pmaxab) // 改 fortran，新建一个变量放出来

int MV(hartree_f_p_functions, fp_max_grid);
int MV(hartree_f_p_functions, lmax_fp);
double MV(hartree_f_p_functions, fp_grid_min);
double MV(hartree_f_p_functions, fp_grid_inc);
double MV(hartree_f_p_functions, fp_grid_max);
double *MV(hartree_f_p_functions, fp_function_spline);    // (0:lmax_Fp,n_max_spline,Fp_max_grid)
double *MV(hartree_f_p_functions, fpc_function_spline);   // (0:lmax_Fp,n_max_spline,Fp_max_grid)
double MV(hartree_f_p_functions, ewald_radius_to)[11];    // 11
double MV(hartree_f_p_functions, inv_ewald_radius_to)[2]; // 2
double MV(hartree_f_p_functions, p_erfc_4)[6];            // 6
double MV(hartree_f_p_functions, p_erfc_5)[7];            // 7
double MV(hartree_f_p_functions, p_erfc_6)[8];            // 8

int MV(hartree_potential_storage, n_rho_multipole_atoms);
int MV(hartree_potential_storage, use_rho_multipole_shmem) = 1;
int *MV(hartree_potential_storage, rho_multipole_index); // (n_atoms)
int *MV(hartree_potential_storage, compensation_norm); // (n_atoms)
int *MV(hartree_potential_storage, compensation_radius); // (n_atoms)
double *MV(hartree_potential_storage,
                  rho_multipole); // ((l_pot_max+1)**2, n_max_radial+2, n_rho_multipole_atoms)
double *MV(hartree_potential_storage, rho_multipole_shmem_ptr);

int *MV(basis, perm_basis_fns_spl);    // (n_basis_fns)  // TOOD
double *MV(basis, outer_radius_sq);    // (n_basis_fns)  // TOOD
int *MV(basis, basis_fn);              // (n_basis)  // TOOD
int *MV(basis, basis_l);               // (n_basis)  // TOOD
double *MV(basis, atom_radius_sq);     // (n_species)  // TOOD
int *MV(basis, basis_fn_start_spl);    // (n_species)  // TOOD
int *MV(basis, basis_fn_atom);         // (n_basis_fns,n_atoms)  // TOOD
double *MV(basis, basis_wave_ordered); // (n_basis_fns,n_max_spline, n_max_grid)  // TOOD
double *MV(basis, basis_kinetic_ordered); // (n_basis_fns,n_max_spline, n_max_grid)  // TOOD

// sumup batch
int MV(opencl_util, n_my_batches_work_sumup);
int MV(opencl_util, n_full_points_work_sumup);
int *MV(opencl_util, batches_size_sumup); // (n_my_batches_work)  // 进程间不同
double *MV(opencl_util,
           batches_points_coords_sumup); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO

// rho batch
int MV(opencl_util, n_my_batches_work_rho);
int MV(opencl_util, n_full_points_work_rho);
int *MV(opencl_util, batches_size_rho);            // (n_my_batches_work)  // 进程间不同
int *MV(opencl_util, batches_batch_n_compute_rho); // (n_my_batches_work)  // 进程间不同
int *MV(opencl_util, batches_batch_i_basis_rho);   // (n_centers_basis_I, n_my_batches_work) // 进程间不同
double *MV(opencl_util,
           batches_points_coords_rho); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO
// ------

// H batch
int MV(opencl_util, n_my_batches_work_h);
int MV(opencl_util, n_full_points_work_h);
int *MV(opencl_util, batches_size_h);            // (n_my_batches_work)  // 进程间不同
int *MV(opencl_util, batches_batch_n_compute_h); // (n_my_batches_work)  // 进程间不同
int *MV(opencl_util, batches_batch_i_basis_h);   // (n_centers_basis_I, n_my_batches_work) // 进程间不同
double *MV(opencl_util,
           batches_points_coords_h); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO
// ------
int MV(opencl_util, max_n_batch_centers);


// use: 1, not use: 0
int MV(opencl_util, use_opencl_version) = 1;
int MV(opencl_util, use_sumup_pre_c_cl_version) = 0;

extern long dgemm_time;
extern long dgemm_flop;

int main(int argc,char *argv[]) {
#ifdef _MY_MPI_DEBUG
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  printf("Process %d of %d\n", myid, n_tasks);
  mpi_platform_relative_id = myid % 32;
  int mpi_sr_number = 2000;
  // char tmpstr[128];
  // sprintf(tmpstr, "%d", (mpi_platform_relative_id/8));
  // setenv("GPU_DEVICE_ORDINAL", tmpstr, 1);
  // setenv("HIP_VISIBLE_DEVICES", tmpstr, 1);
#endif
  printf("---\n");

  // char *file_com_name = "/public/home/autopar/FHI-aims-test/Benchmarks/bigAtoms-cp/strong_weak_test/1w_atoms_mini_local_index/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp/strong_weak_test/3w_atoms_mini_local_index/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp/strong_weak_test/11w_atoms_mini_local_index_energy/sbatch_node512_ntask4096_gpu2048/mdata_outer_rank0_0.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp/strong_weak_test/6w_atoms_mini_local_index/sbatch_node512_ntask4096_gpu2048/mdata_outer_rank0_0.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/c2h4-3k/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/c2h4-3k/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/HIV-1/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/6w_atoms_mini_local_index/sbatch_node256_ntask8192_gpu1024/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_no_all2all/3w_atoms_mini_local_index/sbatch_node64_ntask2048_gpu256/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1w_atoms_mini_local_index/sbatch_node8_ntask256_gpu32/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1w_atoms_mini_local_index/sbatch_node16_ntask512_gpu64/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1w_atoms_mini_local_index/sbatch_node32_ntask1024_gpu128/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node256_ntask8192_gpu1024/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node512_ntask16384_gpu2048/mdata_outer_rank3_0.bin";

  char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1.5w/sbatch_node4_ntask128_gpu16/mdata_outer_rank0_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1.5w/sbatch_node8_ntask256_gpu32/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1.5w/sbatch_node16_ntask512_gpu64/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1.5w/sbatch_node32_ntask1024_gpu128/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/1.5w/sbatch_node64_ntask2048_gpu256/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/3w_atoms_mini_local_index/sbatch_node8_ntask256_gpu32/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/3w_atoms_mini_local_index/sbatch_node16_ntask512_gpu64/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/3w_atoms_mini_local_index/sbatch_node32_ntask1024_gpu128/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/3w_atoms_mini_local_index/sbatch_node64_ntask2048_gpu256/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/3w_atoms_mini_local_index/sbatch_node128_ntask4096_gpu512/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/6w_atoms_mini_local_index/sbatch_node32_ntask1024_gpu128/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/6w_atoms_mini_local_index/sbatch_node64_ntask2048_gpu256/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/6w_atoms_mini_local_index/sbatch_node128_ntask4096_gpu512/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/6w_atoms_mini_local_index/sbatch_node256_ntask8192_gpu1024/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node128_ntask4096_gpu512/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node256_ntask8192_gpu1024/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node512_ntask16384_gpu2048/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node1024_ntask32768_gpu4096/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/11w_atoms_mini_local_index_energy/sbatch_node2048_ntask65536_gpu8192/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/20w_atoms_mini_local_index/sbatch_node512_ntask16384_gpu2048/mdata_outer_rank3_1.bin";
  // char *file_com_name = "/work1/aicao/wzk/Benchmarks/bigAtoms-cp2/strong_weak_test_all2all/20w_atoms_mini_local_index/sbatch_node1024_ntask32768_gpu4096/mdata_outer_rank3_1.bin";

  // char *file_com_name = "/mdata_outer_rank3_1.bin";

  m_save_load(file_com_name, 1, !myid);
  struct timeval start, end;
  // -------- test sumup ----------
  m_save_load_sumup(file_com_name, 1, !myid);

  // printf("-------------------\n\n\n\n\n");
  // m_save_load_check("./check0.bin", 0, 1);
  // return;
  unsigned int numOfPlatforms;

  gettimeofday(&start, NULL);

  // if(myid < 16 && myid >= 8){
  //   // 在同一节点上，只要执行这行的进程数大于 8，那么就会导致 kernel 时间变慢
  //   // 但是对小小 kernel (正常 0.01 以内) 起效情况随机，不定时翻 10 倍
  //   // 对于正常 kernel (0.4 - 1.5 左右)，以典型的 1.0s 为例 (myid < 1)
  //   // (myid < 9) 9 个进程执行下面这行是，慢至 1.5s
  //   // (myid < 16) 递增至 16 个进程执行下面这行时，慢至 2.4s
  //   int error = clGetPlatformIDs(0, NULL, &numOfPlatforms);
  // }


  // if(myid == 0){
  // if(myid == 10){
    init_sum_up_c_(&sum_up_param.forces_on, sum_up_param.partition_tab, sum_up_param.delta_v_hartree,
                   sum_up_param.rho_multipole, sum_up_param.adap_outer_radius_sq,
                   sum_up_param.multipole_radius_sq, sum_up_param.l_hartree_max_far_distance,
                   sum_up_param.outer_potential_radius, sum_up_param.multipole_c);

// MPI_Barrier(MPI_COMM_WORLD);

// //   if(myid < 8){
// #ifdef _MY_MPI_DEBUG
//   // if((mpi_platform_relative_id % 8) != 0){
//   if((mpi_platform_relative_id % 4) != 0){
//     // call MPI_Recv(mpi_buffer, count, MPI_Integer, myid-mpi_platform_num, tag, MPI_COMM_WORLD, status, ierr)
//     MPI_Recv(&mpi_sr_number, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     printf("myid=%d recv %d from id=%d\n", myid, mpi_sr_number, myid-1);
//   }
//   // if((mpi_platform_relative_id / 4) != 0){
//   //   // call MPI_Recv(mpi_buffer, count, MPI_Integer, myid-mpi_platform_num, tag, MPI_COMM_WORLD, status, ierr)
//   //   MPI_Recv(&mpi_sr_number, 1, MPI_INT, myid-4, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//   //   printf("myid=%d recv %d from id=%d\n", myid, mpi_sr_number, myid-4);
//   // }
// #endif
    // if(myid == 13)
    sum_up_whole_potential_shanghui_sub_t_(
        &sum_up_param.forces_on, sum_up_param.partition_tab, sum_up_param.delta_v_hartree, sum_up_param.rho_multipole,
        sum_up_param.adap_outer_radius_sq, sum_up_param.multipole_radius_sq, sum_up_param.l_hartree_max_far_distance,
        sum_up_param.outer_potential_radius, sum_up_param.multipole_c);

// #ifdef _MY_MPI_DEBUG
//   // if((mpi_platform_relative_id % 8) != 7 && myid != (n_tasks - 1)) {
//   if((mpi_platform_relative_id % 4) != 3 && myid != (n_tasks - 1)) {
//     mpi_sr_number = 2000 + myid;
//     MPI_Send(&mpi_sr_number, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
//     printf("myid=%d send %d from id=%d\n", myid, mpi_sr_number, myid+1);
//   }
//   // if((mpi_platform_relative_id / 4) != 7 && (myid+4 < n_tasks)) {
//   //   mpi_sr_number = 2000 + myid;
//   //   MPI_Send(&mpi_sr_number, 1, MPI_INT, myid+4, 0, MPI_COMM_WORLD);
//   //   printf("myid=%d send %d to id=%d\n", myid, mpi_sr_number, myid+4);
//   // }
// #endif
  // }
    // sum_up_final_end();
    // m_load_sumup_free();
    // release_sum_up_c_();
  // }

  // --------- test rho -----------
  // m_save_load_rho(file_com_name, 1, 0);

  // // // struct timeval start, end;
  // gettimeofday(&start, NULL);

  // integrate_first_order_rho_sub_t_(
  //     &(rho_param.l_ylm_max), &(rho_param.n_local_matrix_size), &(rho_param.n_basis_local), &(rho_param.perm_n_full_points),
  //     &(rho_param.first_order_density_matrix_size), rho_param.basis_l_max, rho_param.n_points_all_batches,
  //     rho_param.n_batch_centers_all_batches, rho_param.batch_center_all_batches,
  //     rho_param.batch_point_to_i_full_point, rho_param.ins_idx_all_batches, rho_param.first_order_rho,
  //     rho_param.first_order_density_matrix, rho_param.partition_tab);
  // ------- test rho end ---------

//   m_save_load_H(file_com_name, 1, !myid);
  
// //   struct timeval start, end;
//   gettimeofday(&start, NULL);

//   h_pass_vars_(
//       &(H_param.j_coord), &(H_param.n_spin), &(H_param.l_ylm_max), &(H_param.n_basis_local), &(H_param.n_matrix_size),
//       H_param.basis_l_max, H_param.n_points_all_batches, H_param.n_batch_centers_all_batches,
//       H_param.batch_center_all_batches, H_param.ins_idx_all_batches, H_param.batches_batch_i_basis_h,
//       H_param.partition_all_batches, H_param.first_order_H, H_param.local_potential_parts_all_points,
//       H_param.local_first_order_rho_all_batches, H_param.local_first_order_potential_all_batches,
//       H_param.local_dVxc_drho_all_batches, H_param.local_rho_gradient, H_param.first_order_gradient_rho);


//   integrate_first_order_h_sub_t_(
//       &(H_param.j_coord), &(H_param.n_spin), &(H_param.l_ylm_max), &(H_param.n_basis_local), &(H_param.n_matrix_size),
//       H_param.basis_l_max, H_param.n_points_all_batches, H_param.n_batch_centers_all_batches,
//       H_param.batch_center_all_batches, H_param.ins_idx_all_batches, H_param.batches_batch_i_basis_h,
//       H_param.partition_all_batches, H_param.first_order_H, H_param.local_potential_parts_all_points,
//       H_param.local_first_order_rho_all_batches, H_param.local_first_order_potential_all_batches,
//       H_param.local_dVxc_drho_all_batches, H_param.local_rho_gradient, H_param.first_order_gradient_rho);

  gettimeofday(&end, NULL);
  long timeuse = 1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;
  m_load_free();
  opencl_common_buffer_free_();
  opencl_finish();
  printf("=== time = %8.3f s ===\n",timeuse /1000000.0);

#ifdef _MY_MPI_DEBUG
  MPI_Finalize();
#endif
}

