#ifndef CUDA_FORCES_H
#define CUDA_FORCES_H

#include "gpuMacro.h"
#include <stdio.h>
#include <stdlib.h>

extern "C"
void FORTRAN(forces_create_gpu)(
      int* n_max_compute_ham,
      int* n_max_batch_size,
      int* ld_dens_mat,
      int* n_spin,
      int* n_atoms,
      int* n_basis,
      int* n_basis_local,
      int* as_components,
      int* gga_forces_on,
      int* meta_gga_forces_on,
      int* use_analytical_stress,
      int* rel_atomic_ZORA,
      int* use_as_jac_in_pulay,
      int* load_balanced_matrix);

extern "C"
void FORTRAN(forces_destroy_gpu)();

extern "C"
void FORTRAN(set_h_times_psi_gpu)(
      double* h_times_psi);

extern "C"
void FORTRAN(set_d_h_times_psi_gpu)(
      double* d_h_times_psi);

extern "C"
void FORTRAN(set_wave_gpu)(
      double* wave);

extern "C"
void FORTRAN(set_d_wave_gpu)(
      double* d_wave);

extern "C"
void FORTRAN(set_permute_compute_by_atom_gpu)(
      int* permute_compute_by_atom,
      int* n_compute_c);

extern "C"
void FORTRAN(set_ins_idx_gpu)(
      int* ins_idx,
      int* n_compute_c);

extern "C"
void FORTRAN(set_xc_gradient_deriv_gpu)(
      double* xc_gradient_deriv);

extern "C"
void FORTRAN(set_xc_tau_deriv_gpu)(
      double* xc_tau_deriv);

extern "C"
void FORTRAN(set_gradient_basis_wave_gpu)(
      double* gradient_basis_wave);

extern "C"
void FORTRAN(set_hessian_basis_wave_gpu)(
      double* hessian_basis_wave);

extern "C"
void FORTRAN(set_dens_mat_gpu)(
      double* dens_mat);

extern "C"
void FORTRAN(set_partition_gpu)(
      double* partition);

extern "C"
void FORTRAN(set_sum_forces_gpu)(
      double* sum_forces,
      int* n_atoms);

extern "C"
void FORTRAN(set_as_pulay_stress_local_gpu)(
      double* as_pulay_stress_local);

extern "C"
void FORTRAN(set_as_strain_deriv_wave_gpu)(
      double* as_strain_deriv_wave);

extern "C"
void FORTRAN(set_as_jac_pot_kin_times_psi_gpu)(
      double* as_jac_pot_kin_times_psi);

extern "C"
void FORTRAN(set_as_hessian_times_xc_deriv_gga_gpu)(
      double* as_hessian_times_xc_deriv_gga);

extern "C"
void FORTRAN(set_as_hessian_times_xc_deriv_mgga_gpu)(
      double* as_hessian_times_xc_deriv_mgga);

extern "C"
void FORTRAN(set_as_strain_deriv_kinetic_wave_gpu)(
      double* as_strain_deriv_kinetic_wave);

extern "C"
void FORTRAN(set_as_strain_deriv_wave_shell_gpu)(
      double* as_strain_deriv_wave_shell);

extern "C"
void FORTRAN(get_as_strain_deriv_wave_shell_gpu)(
      double* as_strain_deriv_wave_shell);

extern "C"
void FORTRAN(get_forces_shell_gpu)(
      double* forces_shell);

extern "C"
void FORTRAN(get_sum_forces_gpu)(
      double* sum_forces,
      int* n_atoms);

extern "C"
void FORTRAN(get_as_pulay_stress_local_gpu)(
      double* as_pulay_stress_local);

extern "C"
void FORTRAN(eval_forces_shell_dpsi_h_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin);

extern "C"
void FORTRAN(eval_as_shell_dpsi_h_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin);

extern "C"
void FORTRAN(eval_as_shell_add_psi_kin_psi_shell_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin);

extern "C"
void FORTRAN(eval_forces_shell_psi_dh_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin);

extern "C"
void FORTRAN(as_evaluate_gga_stress_gpu)(
      int* n_compute_c,
      int* n_points,
      int* i_coord,
      int* i_spin);

extern "C"
void FORTRAN(evaluate_forces_shell_add_mgga_gpu)(
     int* n_points,
     int* n_compute_c,
     int* n_compute_a,
     int* i_spin);

extern "C"
void FORTRAN(eval_as_shell_psi_dkin_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin);

extern "C"
void FORTRAN(transpose_as_shell_gpu)(
      int* n_compute_c);

extern "C"
void FORTRAN(eval_gga_forces_dens_mat_gpu)(
      int* n_compute,
      int* n_points,
      int* i_dim,
      int* i_spin);

extern "C"
void FORTRAN(update_sum_forces_gpu)(
      int* n_compute,
      int* i_calculate_dimension,
      int* i_spin,
      int* n_computeForAtom);

extern "C"
void FORTRAN(as_update_sum_forces_and_stress_gpu)(
      int* n_compute,
      int* i_calculate_dimension,
      int* i_spin,
      int* n_computeForAtom);

__global__ void sumForcesAdd(
      double* address,
      double* value);

__global__ void updateSumForcesLoadBalanced(
      int nCompute,
      int iCoord,
      int* permuteComputeByAtom,
      int* indIdx,
      double* matrix,
      double* matrixShell,
      double* forceValues);

__global__ void asUpdateStressLoadBalanced(
      int nCompute,
      int iCoord,
      int* insIdx,
      double* densMat,
      double* asStrainDerivWaveShell,
      double* asStressValues);

#endif /*CUDA_FORCES_H*/

//DO NOT FORGET TO ADD ROUTINES IN CUDA_STUB.f90
