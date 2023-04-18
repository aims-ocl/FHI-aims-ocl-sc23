#include <mpi.h>
#include "opencl_util_mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include "pass_mod_var.h"

#ifdef __INTEL_COMPILER
#define MV(mod_name, var_name) mod_name##_mp_##var_name##_
#else
#define MV(mod_name, var_name) __##mod_name##_MOD_##var_name
#endif

// sumup batch
// extern int MV(opencl_util, n_my_batches_work_sumup);
// extern int MV(opencl_util, n_full_points_work_sumup);
extern int *MV(opencl_util, batches_size_sumup); // (n_my_batches_work)  // 进程间不同
extern double *MV(opencl_util,
                  batches_points_coords_sumup); // (3, n_max_batch_size, n_my_batches_work) // 进程间不同 // TODO

extern int MV(mpi_tasks, n_tasks);
extern int MV(mpi_tasks, myid);
#define n_tasks MV(mpi_tasks, n_tasks)
#define myid MV(mpi_tasks, myid)

#define MY_SET(var1, var2, is_var2_to_var1)    \
  if (is_var2_to_var1) {                       \
    (var1) = (var2);                           \
  } else {                                     \
    (var2) = (var1);                           \
  }

// 使用 buf2##_size < count 的目的应该是，如果有的值有变化，那么分配新的
#define ISEND_RECV(is_send, buf, _count, mpi_datatype, rank, tag, comm, request, buf2, datatype)                       \
  if (is_send) {                                                                                                       \
    MPI_Send((buf), (_count), (mpi_datatype), (rank), (tag), (comm));                                                  \
  } else {                                                                                                             \
    if (_count <= 0) {                                                                                                 \
      printf("rank %d, " #buf2 " count=0\n", myid);                                                                    \
      fflush(stdout);                                                                                                  \
      exit(-1);                                                                                                        \
    }                                                                                                                  \
    if (buf2##___sizen < _count) {                                                                                     \
      if (buf2 != NULL) {                                                                                              \
        free(buf2);                                                                                                    \
      }                                                                                                                \
      buf2##___sizen = _count;                                                                                         \
      buf2 = malloc(sizeof(datatype) * (_count));                                                                      \
    }                                                                                                                  \
    MPI_Status status;                                                                                                 \
    MPI_Recv((buf2), (_count), (mpi_datatype), (rank), (tag), (comm), &status);                                        \
    int count;                                                                                                         \
    MPI_Get_count(&status, mpi_datatype, &count);                                                                      \
    if (_count != count) {                                                                                             \
      printf("rank %d, " #buf " sizeof(datatype)=%d, count=%d, status.count=%d, recv from %d\n", myid,                 \
             sizeof(datatype), _count, count, rank);                                                                   \
      fflush(stdout);                                                                                                  \
      exit(-1);                                                                                                        \
    }                                                                                                                  \
  }

#define ISEND_RECV2(is_send, buf, _count, mpi_datatype, rank, tag, comm, request, buf2, datatype)                      \
  if (is_send) {                                                                                                       \
    MPI_Send((buf), (_count), (mpi_datatype), (rank), (tag), (comm));                                                  \
  } else {                                                                                                             \
    if (_count <= 0) {                                                                                                 \
      printf("rank %d, " #buf2 " count=0\n", myid);                                                                    \
      fflush(stdout);                                                                                                  \
      exit(-1);                                                                                                        \
    }                                                                                                                  \
    MPI_Status status;                                                                                                 \
    MPI_Recv((buf2), (_count), (mpi_datatype), (rank), (tag), (comm), &status);                                        \
    int count;                                                                                                         \
    MPI_Get_count(&status, mpi_datatype, &count);                                                                      \
    if (_count != count) {                                                                                             \
      printf("rank %d, " #buf " sizeof(datatype)=%d, count=%d, status.count=%d, recv from %d\n", myid,                 \
             sizeof(datatype), _count, count, rank);                                                                   \
      fflush(stdout);                                                                                                  \
      exit(-1);                                                                                                        \
    }                                                                                                                  \
  }

// version with print

// // 使用 buf2##_size < count 的目的应该是，如果有的值有变化，那么分配新的
// #define ISEND_RECV(is_send, buf, _count, mpi_datatype, rank, tag, comm, request, buf2, datatype)                       \
//   if (is_send) {                                                                                                       \
//     printf("rank %d, send " #buf "\n", myid);                                                                          \
//     MPI_Send((buf), (_count), (mpi_datatype), (rank), (tag), (comm));                                                  \
//     printf("rank %d, finish send " #buf "\n", myid);                                                                   \
//   } else {                                                                                                             \
//     if (_count <= 0) {                                                                                                 \
//       printf("rank %d, " #buf2 " count=0\n", myid);                                                                    \
//       fflush(stdout);                                                                                                  \
//       exit(-1);                                                                                                        \
//     }                                                                                                                  \
//     printf(#buf2 "___sizen=%d, count=%d\n", buf2##___sizen, _count);                                                   \
//     if (buf2##___sizen < _count) {                                                                                     \
//       if (buf2 != NULL) {                                                                                              \
//         free(buf2);                                                                                                    \
//         if (myid == 0)                                                                                                 \
//           printf("rank %d, resize " #buf2 " from %d to %d\n", myid, buf2##___sizen, _count);                           \
//       } else {                                                                                                         \
//         if (myid == 0)                                                                                                 \
//           printf("rank %d, new malloc " #buf2 " from %d to %d\n", myid, buf2##___sizen, _count);                       \
//       }                                                                                                                \
//       buf2##___sizen = _count;                                                                                         \
//       buf2 = malloc(sizeof(datatype) * (_count));                                                                      \
//     }                                                                                                                  \
//     MPI_Status status;                                                                                                 \
//     printf("rank %d, recv " #buf2 "\n", myid);                                                                         \
//     MPI_Recv((buf2), (_count), (mpi_datatype), (rank), (tag), (comm), &status);                                        \
//     printf("rank %d, finish recv " #buf2 "\n", myid);                                                                  \
//     int count;                                                                                                         \
//     MPI_Get_count(&status, mpi_datatype, &count);                                                                      \
//     printf("rank %d, sizeof(datatype)=%d, count=%d, status.count=%d\n", myid, sizeof(datatype), _count, count);        \
//   }

// #define ISEND_RECV2(is_send, buf, _count, mpi_datatype, rank, tag, comm, request, buf2, datatype)                       \
//   if (is_send) {                                                                                                       \
//     printf("rank %d, send " #buf "\n", myid);                                                                          \
//     MPI_Send((buf), (_count), (mpi_datatype), (rank), (tag), (comm));                                                  \
//     printf("rank %d, finish send " #buf "\n", myid);                                                                   \
//   } else {                                                                                                             \
//     if (_count <= 0) {                                                                                                 \
//       printf("rank %d, " #buf2 " count=0\n", myid);                                                                    \
//       fflush(stdout);                                                                                                  \
//       exit(-1);                                                                                                        \
//     }                                                                                                                  \
//     printf(#buf2 "___sizen=%d, count=%d\n", buf2##___sizen, _count);                                                   \
//     MPI_Status status;                                                                                                 \
//     printf("rank %d, recv " #buf2 "\n", myid);                                                                         \
//     MPI_Recv((buf2), (_count), (mpi_datatype), (rank), (tag), (comm), &status);                                        \
//     printf("rank %d, finish recv " #buf2 "\n", myid);                                                                  \
//     int count;                                                                                                         \
//     MPI_Get_count(&status, mpi_datatype, &count);                                                                      \
//     printf("rank %d, sizeof(datatype)=%d, count=%d, status.count=%d\n", myid, sizeof(datatype), _count, count);        \
//   }

// #define ISEND_RECV2(is_send, buf, _count, mpi_datatype, rank, tag, comm, request, buf2, datatype)                      \
//   if (is_send) {                                                                                                       \
//     MPI_Send((buf), (_count), (mpi_datatype), (rank), (tag), (comm));                                                  \
//   } else {                                                                                                             \
//     if (_count <= 0) {                                                                                                 \
//       printf("rank %d, " #buf " count=0\n", myid);                                                                     \
//       fflush(stdout);                                                                                                  \
//       exit(-1);                                                                                                        \
//     }                                                                                                                  \
//     MPI_Status status;                                                                                                 \
//     MPI_Recv((buf2), (_count), (mpi_datatype), (rank), (tag), (comm), &status);                                        \
//     int count;                                                                                                         \
//     MPI_Get_count(&status, mpi_datatype, &count);                                                                      \
//     printf("rank %d, sizeof(datatype)=%d, count=%d, status.count=%d\n", myid, sizeof(datatype), _count, count);        \
//   }

// #define ISEND_RECV(is_send, buf, count, datatype, rank, tag, comm, request, buf2)                                      \
//   if (is_send) {                                                                                                       \
//     MPI_Isend((buf), (count), (datatype), (rank), (tag), (comm), (request));                                           \
//   } else {                                                                                                             \
//     if (buf2##_size < count) {                                                                                         \
//       if (buf2 != NULL) {                                                                                              \
//         free(buf2);                                                                                                    \
//       }                                                                                                                \
//       buf2 = malloc(sizeof(datatype) * (count));                                                                       \
//     }                                                                                                                  \
//     MPI_Irecv((buf2), (count), (datatype), (rank), (tag), (comm), (request));                                          \
//   }


// void opencl_util_mpi_vars_default_(int* is_send, int* send_or_recv_rank) {
//   opencl_util_mpi_vars(&ocl_util_vars_local, *is_send, *send_or_recv_rank) {
// }


void opencl_util_mpi_vars(OCL_UTIL_VARS* ocl_util_vars, int is_send, int send_or_recv_rank) {
  MPI_Request *request = (MPI_Request *)malloc(sizeof(MPI_Request) * 1024);
  int request_count = 0;
  int int_vars[1024]; // WARNING 1024 随便写的
  int int_vars_count = 0;
  double double_vars[1024]; // WARNING 1024 随便写的
  int double_vars_count = 0;
  // dimensions

  int mpi_tag_num = 3544;

  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_hartree_potential, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_periodic, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_radial, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->l_pot_max, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_spline, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_hartree_grid, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_species, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_atoms, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_basis_integrals, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_integrals, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_fns_ham, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_basis_fns, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_basis, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_basis_T, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_basis_I, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_grid, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_atoms, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_ham, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_dens, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_batch_size, is_send)
  // runtime_choices
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->use_hartree_non_periodic_ewald, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->hartree_fp_function_splines, is_send)
  // MY_SET(int_vars[int_vars_count++],  fast_ylm, is_send)
  // MY_SET(int_vars[int_vars_count++],  new_ylm, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->flag_rel, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->Adams_Moulton_integrator, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->compensate_multipole_errors, is_send)
  // pbc_lists
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->index_hamiltonian_dim2, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->position_in_hamiltonian_dim1, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->position_in_hamiltonian_dim2, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->column_index_hamiltonian_size, is_send)
  // analytic_multipole_coefficients
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->l_max_analytic_multipole, is_send)
  // hartree_potential_real_p0
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_hartree_atoms, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->hartree_force_l_add, is_send)
  // hartree_f_p_functions
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->Fp_max_grid, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->lmax_Fp, is_send)
  // ---------- vars below may change -----------
  // hartree_potential_storage
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_rho_multipole_atoms, is_send)
  // sumup batch
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_my_batches_work_sumup, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_full_points_work_sumup, is_send)
  // rho batch
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_my_batches_work_rho, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_full_points_work_rho, is_send)
  // h batch
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_my_batches_work_h, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_full_points_work_h, is_send)
  // others
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->forces_on, is_send)
  MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_cc_lm_ijk_l_pot_max, is_send)
  // ======================= double vars ========================
  MY_SET(double_vars[double_vars_count++], ocl_util_vars->Fp_grid_min, is_send)
  MY_SET(double_vars[double_vars_count++], ocl_util_vars->Fp_grid_inc, is_send)
  MY_SET(double_vars[double_vars_count++], ocl_util_vars->Fp_grid_max, is_send)

  int mpi_error;
  if (is_send) {
    mpi_error = MPI_Send(int_vars, int_vars_count, MPI_INT, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD);
    if(mpi_error != MPI_SUCCESS){
      printf("rank %d, %s:%d, %s, is_send=%d, send_or_recv_rank=%d\n", myid, __FILE__, __LINE__, __func__, is_send, send_or_recv_rank);
      fflush(stdout);
      exit(-1);
    }
    mpi_error = MPI_Send(double_vars, double_vars_count, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD);
    if(mpi_error != MPI_SUCCESS){
      printf("rank %d, %s:%d, %s, is_send=%d, send_or_recv_rank=%d\n", myid, __FILE__, __LINE__, __func__, is_send, send_or_recv_rank);
      fflush(stdout);
      exit(-1);
    }
  } else {
    mpi_error = MPI_Recv(int_vars, int_vars_count, MPI_INT, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(mpi_error != MPI_SUCCESS){
      printf("rank %d, %s:%d, %s, is_send=%d, send_or_recv_rank=%d\n", myid, __FILE__, __LINE__, __func__, is_send, send_or_recv_rank);
      fflush(stdout);
      exit(-1);
    }
    mpi_error = MPI_Recv(double_vars, double_vars_count, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(mpi_error != MPI_SUCCESS){
      printf("rank %d, %s:%d, %s, is_send=%d, send_or_recv_rank=%d\n", myid, __FILE__, __LINE__, __func__, is_send, send_or_recv_rank);
      fflush(stdout);
      exit(-1);
    }
    int int_vars_count = 0;
    int double_vars_count = 0;
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_hartree_potential, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_periodic, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_radial, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->l_pot_max, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_spline, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_hartree_grid, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_species, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_atoms, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_basis_integrals, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_integrals, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_fns_ham, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_basis_fns, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_basis, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_basis_T, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_centers_basis_I, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_grid, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_atoms, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_ham, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_compute_dens, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_max_batch_size, is_send)
    // runtime_choices
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->use_hartree_non_periodic_ewald, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->hartree_fp_function_splines, is_send)
    // MY_SET(int_vars[int_vars_count++],  fast_ylm, is_send)
    // MY_SET(int_vars[int_vars_count++],  new_ylm, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->flag_rel, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->Adams_Moulton_integrator, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->compensate_multipole_errors, is_send)
    // pbc_lists
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->index_hamiltonian_dim2, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->position_in_hamiltonian_dim1, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->position_in_hamiltonian_dim2, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->column_index_hamiltonian_size, is_send)
    // analytic_multipole_coefficients
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->l_max_analytic_multipole, is_send)
    // hartree_potential_real_p0
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_hartree_atoms, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->hartree_force_l_add, is_send)
    // hartree_f_p_functions
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->Fp_max_grid, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->lmax_Fp, is_send)
    // ---------- vars below may change -----------
    // hartree_potential_storage
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_rho_multipole_atoms, is_send)
    // sumup batch
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_my_batches_work_sumup, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_full_points_work_sumup, is_send)
    // rho batch
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_my_batches_work_rho, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_full_points_work_rho, is_send)
    // h batch
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_my_batches_work_h, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_full_points_work_h, is_send)
    // others
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->forces_on, is_send)
    MY_SET(int_vars[int_vars_count++], ocl_util_vars->n_cc_lm_ijk_l_pot_max, is_send)
    // ======================= double vars ========================
    MY_SET(double_vars[double_vars_count++], ocl_util_vars->Fp_grid_min, is_send)
    MY_SET(double_vars[double_vars_count++], ocl_util_vars->Fp_grid_inc, is_send)
    MY_SET(double_vars[double_vars_count++], ocl_util_vars->Fp_grid_max, is_send)
  }

  free(request);
}

void opencl_util_mpi_arrays_(OCL_UTIL_VARS* ocl_util_vars, int is_send, int send_or_recv_rank) {
  MPI_Request *request = (MPI_Request *)malloc(sizeof(MPI_Request) * 1024);
  // MPI_Status *status = (MPI_Status *)malloc(sizeof(MPI_Status) * 1024);

  int request_count = 0;
  int mpi_tag_num = 5544;

  // sumup batch
  ISEND_RECV(is_send, ocl_util_vars->batches_size_sumup, ocl_util_vars->n_my_batches_work_sumup, MPI_INT, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->batches_size_sumup, int);
  ISEND_RECV(is_send, ocl_util_vars->batches_points_coords_sumup, 3 * ocl_util_vars->n_max_batch_size * ocl_util_vars->n_my_batches_work_sumup, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->batches_points_coords_sumup, double);

  // == MPI_Isend(&sum_up_param.forces_on, 1, MPI_INT, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++]);
  ISEND_RECV(is_send, ocl_util_vars->partition_tab, ocl_util_vars->n_full_points_work_sumup, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->partition_tab, double);
  ISEND_RECV(is_send, ocl_util_vars->delta_v_hartree, ocl_util_vars->n_full_points_work_sumup, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->delta_v_hartree, double);
  ISEND_RECV(is_send, ocl_util_vars->rho_multipole, ocl_util_vars->n_full_points_work_sumup, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->rho_multipole, double);
  ISEND_RECV(is_send, ocl_util_vars->adap_outer_radius_sq, ocl_util_vars->n_atoms, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->adap_outer_radius_sq, double);
  ISEND_RECV(is_send, ocl_util_vars->multipole_radius_sq, ocl_util_vars->n_atoms, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->multipole_radius_sq, double);
  ISEND_RECV(is_send, ocl_util_vars->l_hartree_max_far_distance, ocl_util_vars->n_atoms, MPI_INT, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->l_hartree_max_far_distance, int);
  ISEND_RECV(is_send, ocl_util_vars->outer_potential_radius, (ocl_util_vars->l_pot_max + 1) * ocl_util_vars->n_atoms, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->outer_potential_radius, double);
  ISEND_RECV(is_send, ocl_util_vars->multipole_c, ocl_util_vars->n_cc_lm_ijk_l_pot_max * ocl_util_vars->n_atoms, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->multipole_c, double);

  // ISEND_RECV(is_send, MV(hartree_potential_storage, rho_multipole_index), n_atoms, MPI_INT, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->rho_multipole_index);
  // ISEND_RECV(is_send, MV(hartree_potential_storage, rho_multipole), (l_pot_max + 1) * (l_pot_max + 1) * (n_max_radial + 2) * n_rho_multipole_atoms, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->rho_multipole);

  // ISEND_RECV(is_send, MV(hartree_f_p_functions, fp_function_spline), (lmax_Fp + 1) * n_max_spline * Fp_max_grid, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->fp_function_spline);
  // ISEND_RECV(is_send, MV(hartree_f_p_functions, fpc_function_spline), (lmax_Fp + 1) * n_max_spline * Fp_max_grid, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->fpc_function_spline);
  // ISEND_RECV(is_send, MV(hartree_f_p_functions, ewald_radius_to), 11, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->ewald_radius_to);
  // ISEND_RECV(is_send, MV(hartree_f_p_functions, inv_ewald_radius_to), 2, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->inv_ewald_radius_to);
  // ISEND_RECV(is_send, MV(hartree_f_p_functions, p_erfc_4), 6, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->p_erfc_4);
  // ISEND_RECV(is_send, MV(hartree_f_p_functions, p_erfc_5), 7, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->p_erfc_5);
  // ISEND_RECV(is_send, MV(hartree_f_p_functions, p_erfc_6), 8, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->p_erfc_6);

  // MPI_Waitall(request_count, request, status);

  free(request);
  // free(status);
}

void opencl_util_mpi_arrays_results_(OCL_UTIL_VARS* ocl_util_vars, int is_send, int send_or_recv_rank) {
  MPI_Request *request = (MPI_Request *)malloc(sizeof(MPI_Request) * 1024);
  // MPI_Status *status = (MPI_Status *)malloc(sizeof(MPI_Status) * 1024);

  int request_count = 0;
  int mpi_tag_num = 6544;

  ISEND_RECV2(is_send, ocl_util_vars->delta_v_hartree, ocl_util_vars->n_full_points_work_sumup, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->delta_v_hartree, double);
  ISEND_RECV2(is_send, ocl_util_vars->rho_multipole, ocl_util_vars->n_full_points_work_sumup, MPI_DOUBLE, send_or_recv_rank, mpi_tag_num++, MPI_COMM_WORLD, &request[request_count++], ocl_util_vars->rho_multipole, double);

  free(request);
  // free(status);
}



// #include <hipblas.h>
// #include <hip/hip_profile.h>
// #include <hip/hip_runtime.h>

// extern int MV(opencl_util, mpi_platform_relative_id);
// #define mpi_platform_relative_id MV(opencl_util, mpi_platform_relative_id)
// void hip_and_hipblas_init_(){
//   hipSetDevice(mpi_platform_relative_id / 8);
//   hipDeviceReset();
// }