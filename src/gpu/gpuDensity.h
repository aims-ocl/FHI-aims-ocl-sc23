/****************************************************************************
 ***  KS_density_densmat                                                  ***
 ****************************************************************************/

#ifndef CUDA_DENSITY_H
#define CUDA_DENSITY_H

#include "gpuMacro.h"

extern "C"
void FORTRAN(density_create_gpu)(
      int* n_max_compute_dens,
      int* n_max_batch_size,
      int* n_full_points,
      int* n_spin,
      int* use_density_gradient);

extern "C"
void FORTRAN(density_destroy_gpu)();

extern "C"
void FORTRAN(set_delta_rho_gpu)(
      double* delta_rho,
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(set_rho_gpu)(
      double* rho,
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(set_rho_change_gpu)(
      double* rho_change,
      int* n_spin);

extern "C"
void FORTRAN(set_partition_tab_gpu)(
      double* partition_tab,
      int* n_full_points);

extern "C"
void FORTRAN(set_hartree_partition_tab_gpu)(
      double* hartree_partition_tab,
      int* n_full_points);

extern "C"
void FORTRAN(set_delta_rho_gradient_gpu)(
      double* delta_rho_gradient,
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(set_rho_gradient_gpu)(
      double* rho_gradient,
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(set_density_matrix_gpu)(
      double* density_matrix, //matrix (n_compute,n_compute)
      int* n_compute); //number of relevant basis functions

extern "C"
void FORTRAN(get_delta_rho_gpu) (
      double* delta_rho,
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(get_delta_rho_gradient_gpu) (
      double* delta_rho_gradient,
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(get_rho_change_gpu) (
      double* rho_change,
      int* n_spin);

extern "C"
void FORTRAN(evaluate_ks_density_densmat_gpu)(
      int* n_points, // number of integration points
      double* waves, // matrix (n_basis_compute, n_points)
      int* n_compute, // number of relevant basis functions
      int* n_basis_compute); //maximum number of relevant basis functions

extern "C"
void FORTRAN(update_delta_rho_ks_gpu)(
      int* offset_full_points, // starting point in density difference
      int* n_batch_points, // number of Points in current Batch
      int* n_full_points, // number of all Points
      int* n_points_compute, // non zero points in rho
      int* i_spin, // actual Spin Channel
      int* n_spin); // number of Spin Channels

extern "C"
void FORTRAN(calculate_rho_change_gpu) (
      int* n_full_points,
      int* n_spin);

extern "C"
void FORTRAN(eval_density_grad_densmat_gpu)(
      int* n_points, // number of grid points
      double* grad_waves, // gradients of basis functions
      int* n_compute, // number of relevant basis functions
      int* n_basis_compute); // total number of basis functions

extern "C"
void FORTRAN(update_grad_delta_rho_ks_gpu)(
      int* n_full_points, // number of all Points
      int* n_points_compute, // non zero points in rho
      int* i_spin, // actual Spin Channel
      int* n_spin); // number of Spin Channels

template<typename T>
__global__ void gpuInitVectorDensity(
      T* vector,
      const T val,
      const ssize_t dim);

__device__ double myAtomicAddDensity(
      double* address,
      double val);

template<typename T>
__global__ void gpuDiag(
      T* result,
      const T* matrix,
      const int dim,
      const int elements);

__global__ void generateBatchPoint2FullPoint(
      int nBatchPoints,
      int offset,
      double* partitionTab,
      double* hartreePartitionTab,
      int* batch2Full);

__global__ void distributeDensityUpdate(
      int nPoints,
      int nFullPoints,
      int nSpin,
      int iSpin,
      int* batch2Full,
      double* tempRho,
      double* rho,
      double* deltaRho);

__global__ void updateRhoChange(
      int nPoints,
      int nSpin,
      double* rhoChange,
      double* deltaRho,
      double* partitionTab);

__global__ void distributeGradDensityUpdate(
      int nPoints,
      int nFullPoints,
      int nSpin,
      int iSpin,
      int* batch2Full,
      double* tempGradRho,
      double* gradRho,
      double* deltaGradRho);

#endif /*CUDA_DENSITY_H*/
