#ifndef CUDA_INTEGRATEHAM_H
#define CUDA_INTEGRATEHAM_H

#include "gpuMacro.h"
#include <stdio.h>
#include <stdlib.h>

extern "C"
void FORTRAN(hamiltonian_create_gpu)(
      int* n_max_compute_ham,
      int* n_max_batch_size,
      int* ld_hamiltonian,
      int* n_spin,
      int* use_meta_gga,
      int* use_ZORA,
      int* index_on_gpu);

extern "C"
void FORTRAN(hamiltonian_destroy_gpu)();

extern "C"
void FORTRAN(get_hamiltonian_shell_gpu)(
      double *hamiltonian_shell, // matrix (nCompute, nCompute)
      int *n_compute); // # non zero elements

extern "C"
void FORTRAN(set_hamiltonian_gpu)(
      double *hamiltonian,
      int* dim_hamiltonian);

extern "C"
void FORTRAN(get_hamiltonian_gpu)(
      double *hamiltonian,
      int* dim_hamiltonian);

extern "C"
void FORTRAN(evaluate_hamiltonian_shell_gpu)(
      int* n_points, // # integration points
      const double* partition, // array (nPoints)
      int* n_compute, // # non zero elements
      const double* h_times_psi, // matrix (nBasisList,nPoints)
      int* n_basis_list, // # basis functions
      double* wave, // matrix (nBasisList,nPoints)
      double* hamiltonian_shell); // matrix (nCompute, nCompute)

extern "C"
void FORTRAN(mgga_contribution_gpu)(
      int* n_compute_1,
      int* n_compute_2,
      int* n_points,
      double* left_side_of_mgga_dot_product, // matrix (n_compute_1,3*n_points)
      double* gradient_basis_wave_store); // matrix (n_compute_2,3*n_points)

extern "C"
void FORTRAN(add_zora_matrix_gpu)(
      double *zora_vec_1, // 1. Zora tensor (nBasisList, 3, nRelPoints)
      double *zora_vec_2, // 2. Zora tensor (nBasisList, 3, nRelPoints)
      int *n_basis_list, // # basis functions
      int *n_rel_points, // # relativity points
      int *n_compute); // # non zero elements

extern "C"
void FORTRAN(update_full_matrix_via_map_gpu)(
      double* hamiltonian, // full storage in buffer 3 (1 D Object)
      int* dim_hamiltonian,
      double* hamiltonian_shell, // precalculated in buffer 2 (2 D Object)
      int* dim1_hamiltonian_shell, // number of relevant basis functions
      int* dim2_hamiltonian_shell, // number of relevant basis functions
      int* map, // array relevant basis functions
      int* dim_map);

extern "C"
void FORTRAN(update_batch_matrix_gpu)(
      int* i_spin, // spin value
      int* ld_matrix, // leading dimension of matrix
      int* ins_index, // index array for position in matrix
      int* n_compute_c); // number of nonzero basis functions

__global__ void insertInHamiltonianViaMap(
      double* hamiltonian,
      int dimHamiltonian,
      double* hamiltonianShell,
      int dim1HamiltonianShell,
      int dim2HamiltonianShell,
      int* map,
      int dimMap);

__global__ void insertShellInMatrix(
      double* matrix,
      double* matrixShell,
      int* insIndex,
      int ldMatrix,
      int iSpin,
      int ldShell);

#endif /*CUDA_INTEGRATEHAM_H*/
