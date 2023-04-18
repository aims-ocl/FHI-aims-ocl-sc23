/* FHI-aims/gpu/gpuForces.cu
   NAME
     gpuForces
   SYNOPSIS

   PURPOSE
     A module for accelerating various contributions to the forces.  This
     module assumes that the electron density update is density-matrix-based;
     the orbital-based electron density update is used for calculations
     involving small molecules which likely won't benefit greatly from GPU
     acceleration.

     As of this writing (2018 January 9), the following force contributions are
     supported:
       o Pulay forces
       o GGA forces
       o meta-GGA forces
       o Force correction due to atomic ZORA
     The following force contributions are included as part of the batch
     integration routine which will be executed by the associated CPU code, but
     have not been GPU accelerated yet:
       o Gnonmf forces
     The following force contributions are not part of the batch integration
     routine, and thus are not accelerated here:
       o Hellman-Feynman forces
       o Electrostatic multipole corrections
       o van der Waal forces
       o EXX forces
       o and more...
   AUTHOR
     William Huhn and Bjoern Lange, Duke University
   SEE ALSO
     Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
     Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
     "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
     Computer Physics Communications (2008), submitted.
   COPYRIGHT
      Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
      e.V. Please note that any use of the "FHI-aims-Software" is subject to
      the terms and conditions of the respective license agreement."
   HISTORY
      2018 January - Updated.
*/

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "gpuError.h"
#include "gpuInterface.h"
#include "gpuMacro.h"
// The next include should be the Header file of the routine
#include "gpuForces.h"

// Global Variables
namespace gpuForces {
      // Arrays used by Pulay forces (and others)
      double* dev_partition = NULL;        // partition weights (n_max_batch_size)
      double* dev_waveComputeA = NULL;     // temporary matrix to store partition*|psi> (n_max_compute_ham, n_max_batch_size)
      double* dev_forcesShell = NULL;      // H_ij shell (n_max_compute_ham, n_max_compute_ham)
      // Arrays used exclusively for Pulay forces
      double* dev_hTimesPsi = NULL;        // H|Psi> (n_max_compute_ham, n_max_batch_size, n_spin_local)
                                           // NOTE:  If we are calculating the term involving the energy-weighted density matrix,
                                           // HTimesPsi will contain Psi instead.
      double* dev_dWave = NULL;            // |dPsi> (n_max_compute_ham, n_max_batch_size, 3)
      // Array used for pretty much everything except relativistic LDA forces-only calculations
      double* dev_wave = NULL;             // |Psi> (n_max_compute_ham, n_max_batch_size)
      // Arrays used exclusively for GGA (and meta-GGA) forces and analytical stress tensor
      double* dev_xcGradientDeriv = NULL;  // derivative of the xc Functional (3, n_spin, n_max_batch_size)
      double* dev_hessianBasisWave = NULL; // 2nd derivative of wave (3, n_spin, n_max_batch_size)
      double* dev_gradientBasisWave = NULL; // |Psi> (n_max_compute_ham * 3,  n_max_batch_size)
      double* dev_matrixTmp1 = NULL;      // temporary matrix for GGA forces (n_max_compute_ham, n_points)
      double* dev_matrixTmp2 = NULL;      // temporary matrix for GGA forces (n_points, n_max_compute_ham)
      double* dev_matrixTmp3 = NULL;      // temporary matrix for GGA analytical stress tensor (3, n_points)
      // Arrays used exclusively for meta-GGA forces
      double* dev_xcTauDeriv = NULL;       // ??? (n_spin, n_max_batch_size)
      double* dev_matrixTmp1MGGA = NULL; // temporary matrix for meta-GGA forces (n_max_compute_ham, n_points)
      double* dev_matrixTmp2MGGA = NULL; // temporary matrix for meta-GGA forces (n_points, n_max_compute_ham)
      // Arrays used exclusively for atomic ZORA corrections
      double* dev_dHTimesPsi = NULL;            // d/dr H|Psi> (n_max_compute_ham, n_max_batch_size, 3, n_spin_local)
      // Arrays used exclusively for analytical stress calculations
      double* dev_asStrainDerivWave = NULL;
      double* dev_asStrainDerivWaveShell = NULL;
      double* dev_asStrainDerivWaveShellTrans = NULL;
      double* dev_asJacPotKinTimesPsi = NULL;
      double* dev_asStrainDerivKineticWave = NULL;
      double* dev_asHessianTimesXCDerivGGA = NULL;
      double* dev_asHessianTimesXCDerivMGGA = NULL;
      // Arrays needed when indexing on the GPU
      double* dev_densMat = NULL;          // the relevant density matrix, possibly energy-weighted (n_ham_matrix_size, n_spin)
      int* dev_permuteComputeByAtom = NULL;
      int* dev_insIdx = NULL;
      double* dev_forceValues = NULL;      // temporary storage of the force contributions (n_max_compute_ham * n_max_compute_ham)
      double* dev_theNumberOneHaHaHa = NULL;  // Technically, this is a variant of The Count defined on \mathbb{R}, not \mathbb{N}
                                              // Technically technically, defined on the IEEE 754 approximation of \mathbb{R}
                                              // Who cares, GPU make code go fast
      double* dev_forceComponent = NULL;
      double* dev_sumForces = NULL;        // the resulting (partial) sum of forces (3,nAtoms)
      double* dev_asPulayStressLocal = NULL;
      // Array dimension sizes and miscellaneous settings
      int nMaxComputeHam = -1;
      int ldDensMat = -1;
      int nMaxBatchSize = -1;
      int nSpin = -1;
      int nAtoms = -1;
      int nBasis = -1;
      int nBasisLocal = -1;
      int asComponents = -1;
      bool ggaForcesOn = false;         // Are we calculating GGA *or* meta-GGA forces?
      bool metaGGAForcesOn = false;     // Are we calculating meta-GGA forces (if true, ggaForcesOn must also be true!)
      bool useAnalyticalStress = false; // Are we calculating the analytical stress tensor?
      bool relAtomicZORA = false;       // Do we need to add a correction term for atomic ZORA?
      bool useASJacInPulay = false;     // Technical parameter for analytical stress, almost always true
      bool loadBalancedMatrix = false;  // Are we using a load-balanced matrix and thus perform indexing on the GPU?
}

/*******************************************************************************
**                               CUDA CPU Code                                **
*******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                            Initialize/Finalize                             //
////////////////////////////////////////////////////////////////////////////////

// Allocate memory on GPU for variables
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
      int* load_balanced_matrix)
{
   using namespace gpuForces;

   CHECK(*n_max_compute_ham > 0);
   CHECK(*n_max_batch_size > 0);
   CHECK(*ld_dens_mat > 0);
   CHECK(*n_spin > 0);
   CHECK(*n_atoms > 0);
   CHECK(*n_basis > 0);
   // n_basis_local is only valid when load balancing enabled
   CHECK(*as_components == 6 || *as_components == 9);

   nMaxComputeHam = *n_max_compute_ham;
   nMaxBatchSize = *n_max_batch_size;
   ldDensMat = *ld_dens_mat;
   nSpin = *n_spin;
   nAtoms = *n_atoms;
   nBasis = *n_basis;
   asComponents = *as_components;

   // Set up information about runtime mode
   // We pass these flags as integers in the Fortran code, because there is no
   // standard for defining logicals in Fortran
   if (*gga_forces_on == 0) {
      ggaForcesOn = false;
   }
   else {
      ggaForcesOn = true;
   }
   if (*meta_gga_forces_on == 0) {
      metaGGAForcesOn = false;
   }
   else {
      metaGGAForcesOn = true;
   }
   if (*use_analytical_stress == 0) {
      useAnalyticalStress = false;
   }
   else {
      useAnalyticalStress = true;
   }
   if (*rel_atomic_ZORA == 0) {
      relAtomicZORA = false;
   }
   else {
      relAtomicZORA = true;
   }
   if (*use_as_jac_in_pulay == 0) {
      useASJacInPulay = false;
   }
   else {
      useASJacInPulay = true;
   }
   if (*load_balanced_matrix == 0) {
      loadBalancedMatrix = false;
   }
   else {
      loadBalancedMatrix = true;
      nBasisLocal = *n_basis_local;
   }

   // Allocate arrays which are used by Pulay force calculations, and thus
   // always applicable for all possible calculations
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_partition,
         nMaxBatchSize * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_waveComputeA,
         nMaxComputeHam * nMaxBatchSize * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_forcesShell,
         nMaxComputeHam * nMaxComputeHam * sizeof(double)));

   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_hTimesPsi,
         nMaxComputeHam * nMaxBatchSize * nSpin * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_dWave,
         nMaxComputeHam * nMaxBatchSize * 3 * sizeof(double)));

   // Allocate arrays which are needed only for certain force components
   if (ggaForcesOn || relAtomicZORA || useASJacInPulay) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_wave,
            nMaxComputeHam * nMaxBatchSize * sizeof(double)));
   }

   if (ggaForcesOn) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_xcGradientDeriv,
            3 * nSpin * nMaxBatchSize * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_hessianBasisWave,
            6 * nMaxComputeHam * nMaxBatchSize * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_gradientBasisWave,
            nMaxComputeHam * 3 * nMaxBatchSize * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_matrixTmp1,
            nMaxComputeHam * nMaxBatchSize * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_matrixTmp2,
            nMaxComputeHam * nMaxBatchSize * sizeof(double)));
   }

   if (metaGGAForcesOn) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_xcTauDeriv,
            nSpin * nMaxBatchSize * sizeof(double)));
      if (useAnalyticalStress) {
         // When calculating stress, make these work matrices three times larger
         // so they can be reused for forces and stress calculations
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_matrixTmp1MGGA,
               nMaxComputeHam * 3 * nMaxBatchSize * sizeof(double)));
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_matrixTmp2MGGA,
               nMaxComputeHam * 3 * nMaxBatchSize * sizeof(double)));
      }
      else {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_matrixTmp1MGGA,
               nMaxComputeHam * nMaxBatchSize * sizeof(double)));
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_matrixTmp2MGGA,
               nMaxComputeHam * nMaxBatchSize * sizeof(double)));
      }
   }

   if (useAnalyticalStress) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_asStrainDerivWave,
            nMaxComputeHam * nMaxBatchSize * asComponents * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_asStrainDerivWaveShell,
            nMaxComputeHam * nMaxComputeHam * sizeof(double)));
      if (ggaForcesOn) {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_matrixTmp3,
               3 * nMaxBatchSize * sizeof(double)));
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_asHessianTimesXCDerivGGA,
               nMaxComputeHam * nMaxBatchSize * nSpin * asComponents *
               sizeof(double)));
      }
      if (metaGGAForcesOn) {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_asHessianTimesXCDerivMGGA,
               nMaxComputeHam * nMaxBatchSize * nSpin * asComponents *
               sizeof(double)));
      }
      if (useASJacInPulay) {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_asJacPotKinTimesPsi,
               nMaxComputeHam * nMaxBatchSize * nSpin * asComponents *
               sizeof(double)));
      }
      if (relAtomicZORA) {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_asStrainDerivKineticWave,
               nMaxComputeHam * nMaxBatchSize * asComponents *
               sizeof(double)));
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_asStrainDerivWaveShellTrans,
               nMaxComputeHam * nMaxComputeHam * sizeof(double)));
      }
   }

   if (relAtomicZORA) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_dHTimesPsi,
            nMaxComputeHam * nMaxBatchSize * nSpin * 3 * sizeof(double)));
   }

   // Allocate arrays which are needed when indexing matrices on the GPU
   if (loadBalancedMatrix) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_insIdx,
            nBasisLocal * sizeof(int)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_densMat,
            ldDensMat * nSpin * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_permuteComputeByAtom,
            nMaxComputeHam * sizeof(int)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_forceValues,
            nMaxComputeHam * nMaxComputeHam * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_theNumberOneHaHaHa,
            sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_forceComponent,
            sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_sumForces,
            3 * nAtoms * sizeof(double)));

      // Initialize the dummy dev_theNumberOneHaHaHa variable to 1
      double theNumberOneHaHaHa = 1.0;
      HANDLE_CUBLAS(cublasSetVector(
            1, sizeof(double),
            &theNumberOneHaHaHa, 1,
            dev_theNumberOneHaHaHa, 1));

      // Initialize SumForces to Zero
      double alpha = 0.0;
      double beta = 0.0;
      HANDLE_CUBLAS(cublasDgeam(
            cublasHandle,
            CUBLAS_OP_N, CUBLAS_OP_N,
            3, nAtoms,
            &alpha,
            dev_sumForces, 3,
            &beta,
            dev_sumForces, 3,
            dev_sumForces, 3));

      if (useAnalyticalStress) {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_asPulayStressLocal,
               asComponents * sizeof(double)));
         // Initialize stress to zero
         HANDLE_CUBLAS(cublasDgeam(
               cublasHandle,
               CUBLAS_OP_N, CUBLAS_OP_N,
               1, asComponents,
               &alpha,
               dev_asPulayStressLocal, 1,
               &beta,
               dev_asPulayStressLocal, 1,
               dev_asPulayStressLocal, 1));
      }
   }
}

// Free memory on GPU for variables which was allocated previously by
// allocate_variables_forces_gpu
void FORTRAN(forces_destroy_gpu)()
{
   using namespace gpuForces;

   HANDLE_CUDA(cudaFree(dev_partition));
   HANDLE_CUDA(cudaFree(dev_waveComputeA));
   HANDLE_CUDA(cudaFree(dev_forcesShell));

   HANDLE_CUDA(cudaFree(dev_hTimesPsi));
   HANDLE_CUDA(cudaFree(dev_dWave));

   if (ggaForcesOn || relAtomicZORA || useASJacInPulay) {
      HANDLE_CUDA(cudaFree(dev_wave));
   }

   if (ggaForcesOn) {
      HANDLE_CUDA(cudaFree(dev_xcGradientDeriv));
      HANDLE_CUDA(cudaFree(dev_hessianBasisWave));
      HANDLE_CUDA(cudaFree(dev_gradientBasisWave));
      HANDLE_CUDA(cudaFree(dev_matrixTmp1));
      HANDLE_CUDA(cudaFree(dev_matrixTmp2));
   }

   if (metaGGAForcesOn) {
      HANDLE_CUDA(cudaFree(dev_xcTauDeriv));
      HANDLE_CUDA(cudaFree(dev_matrixTmp1MGGA));
      HANDLE_CUDA(cudaFree(dev_matrixTmp2MGGA));
   }

   if (useAnalyticalStress) {
      HANDLE_CUDA(cudaFree(dev_asStrainDerivWave));
      HANDLE_CUDA(cudaFree(dev_asStrainDerivWaveShell));
      if (useASJacInPulay) {
         HANDLE_CUDA(cudaFree(dev_asJacPotKinTimesPsi));
      }
      if (ggaForcesOn) {
         HANDLE_CUDA(cudaFree(dev_matrixTmp3));
         HANDLE_CUDA(cudaFree(dev_asHessianTimesXCDerivGGA));
      }
      if (metaGGAForcesOn) {
         HANDLE_CUDA(cudaFree(dev_asHessianTimesXCDerivMGGA));
      }
      if (relAtomicZORA) {
         HANDLE_CUDA(cudaFree(dev_asStrainDerivKineticWave));
         HANDLE_CUDA(cudaFree(dev_asStrainDerivWaveShellTrans));
      }
   }

   if (relAtomicZORA) {
      HANDLE_CUDA(cudaFree(dev_dHTimesPsi));
   }

   if (loadBalancedMatrix) {
      HANDLE_CUDA(cudaFree(dev_insIdx));
      HANDLE_CUDA(cudaFree(dev_densMat));
      HANDLE_CUDA(cudaFree(dev_permuteComputeByAtom));
      HANDLE_CUDA(cudaFree(dev_forceValues));
      HANDLE_CUDA(cudaFree(dev_theNumberOneHaHaHa));
      HANDLE_CUDA(cudaFree(dev_forceComponent));
      HANDLE_CUDA(cudaFree(dev_sumForces));
      if (useAnalyticalStress) {
         HANDLE_CUDA(cudaFree(dev_asPulayStressLocal));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
//                               Data Movement                                //
////////////////////////////////////////////////////////////////////////////////

// Various setter subroutines (transferring data from the CPU to the GPU).
// Pretty boilerplate stuff, won't bother commenting each one.

void FORTRAN(set_h_times_psi_gpu)(
      double* h_times_psi)
{
   using namespace gpuForces;

   CHECK(h_times_psi != NULL);
   CHECK(dev_hTimesPsi != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nMaxComputeHam * nMaxBatchSize * nSpin, sizeof(double),
         h_times_psi, 1,
         dev_hTimesPsi, 1));
}

void FORTRAN(set_d_h_times_psi_gpu)(
      double* d_h_times_psi)
{
   using namespace gpuForces;

   CHECK(d_h_times_psi != NULL);
   CHECK(dev_dHTimesPsi != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nMaxComputeHam * nMaxBatchSize * 3 * nSpin, sizeof(double),
         d_h_times_psi, 1,
         dev_dHTimesPsi, 1));
}

void FORTRAN(set_wave_gpu)(
      double* wave)
{
   using namespace gpuForces;

   CHECK(wave != NULL);
   CHECK(dev_wave != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nMaxComputeHam * nMaxBatchSize, sizeof(double),
         wave, 1,
         dev_wave, 1));
}

void FORTRAN(set_d_wave_gpu)(
      double* d_wave)
{
   using namespace gpuForces;

   CHECK(d_wave != NULL);
   CHECK(dev_dWave != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nMaxComputeHam * nMaxBatchSize * 3, sizeof(double),
         d_wave, 1,
         dev_dWave, 1));
}

void FORTRAN(set_permute_compute_by_atom_gpu)(
      int* permute_compute_by_atom,
      int* n_compute_c)
{
   using namespace gpuForces;

   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_c <= nMaxComputeHam);
   CHECK(permute_compute_by_atom != NULL);
   CHECK(dev_permuteComputeByAtom != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         *n_compute_c, sizeof(int),
         permute_compute_by_atom, 1,
         dev_permuteComputeByAtom, 1));
}

void FORTRAN(set_ins_idx_gpu)(
      int* ins_idx,
      int* n_compute_c)
{
   using namespace gpuForces;

   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_c <= nBasisLocal);
   CHECK(ins_idx != NULL);
   CHECK(dev_insIdx != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         *n_compute_c, sizeof(int),
         ins_idx, 1,
         dev_insIdx, 1));
}

void FORTRAN(set_xc_gradient_deriv_gpu)(
      double* xc_gradient_deriv)
{
   using namespace gpuForces;

   CHECK(xc_gradient_deriv != NULL);
   CHECK(dev_xcGradientDeriv != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         3 * nSpin * nMaxBatchSize, sizeof(double),
         xc_gradient_deriv, 1,
         dev_xcGradientDeriv, 1));
}

void FORTRAN(set_xc_tau_deriv_gpu)(
      double* xc_tau_deriv)
{
   using namespace gpuForces;

   CHECK(xc_tau_deriv != NULL);
   CHECK(dev_xcTauDeriv != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nSpin * nMaxBatchSize, sizeof(double),
         xc_tau_deriv, 1,
         dev_xcTauDeriv, 1));
}

void FORTRAN(set_gradient_basis_wave_gpu)(
       double* gradient_basis_wave)
{
   using namespace gpuForces;

   CHECK(gradient_basis_wave != NULL);
   CHECK(dev_gradientBasisWave != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nMaxComputeHam * 3 * nMaxBatchSize, sizeof(double),
         gradient_basis_wave, 1,
         dev_gradientBasisWave, 1));
}

void FORTRAN(set_hessian_basis_wave_gpu)(
      double* hessian_basis_wave)
{
   using namespace gpuForces;
   CHECK(hessian_basis_wave != NULL);
   CHECK(dev_hessianBasisWave != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         6 * nMaxComputeHam * nMaxBatchSize, sizeof(double),
         hessian_basis_wave, 1,
         dev_hessianBasisWave, 1));
}

void FORTRAN(set_dens_mat_gpu)(
      double* dens_mat)
{
   using namespace gpuForces;

   CHECK(dens_mat != NULL);
   CHECK(dev_densMat != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         ldDensMat * nSpin, sizeof(double),
         dens_mat, 1,
         dev_densMat, 1));
}

void FORTRAN(set_partition_gpu)(
      double* partition)
{
   using namespace gpuForces;

   CHECK(partition != NULL);
   CHECK(dev_partition != NULL);

   HANDLE_CUBLAS(cublasSetVector(
         nMaxBatchSize, sizeof(double),
         partition, 1,
         dev_partition, 1));
}

void FORTRAN(set_sum_forces_gpu)(
      double* sum_forces,
      int*    n_atoms)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         3 * *n_atoms, sizeof(double),
         sum_forces, 1,
         dev_sumForces, 1,
         "Setting dev_sumForces in set_sum_forces_gpu"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_pulay_stress_local_gpu)(
      double* as_pulay_stress_local)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         asComponents, sizeof(double),
         as_pulay_stress_local, 1,
         dev_asPulayStressLocal, 1,
         "Setting as_pulay_stress_local in set_as_pulay_stress_local_gpu"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_strain_deriv_wave_gpu)(
      double* as_strain_deriv_wave)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         nMaxComputeHam * nMaxBatchSize * asComponents, sizeof(double),
         as_strain_deriv_wave, 1,
         dev_asStrainDerivWave, 1,
         "Setting as_strain_deriv_wave in set_as_strain_deriv_wave"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_jac_pot_kin_times_psi_gpu)(
      double* as_jac_pot_kin_times_psi)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         nMaxComputeHam * nMaxBatchSize * nSpin * asComponents, sizeof(double),
         as_jac_pot_kin_times_psi, 1,
         dev_asJacPotKinTimesPsi, 1,
         "Setting as_jac_pot_kin_times_psi in set_as_jac_pot_kin_times_psi"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_hessian_times_xc_deriv_gga_gpu)(
      double* as_hessian_times_xc_deriv_gga)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         nMaxComputeHam * nMaxBatchSize * nSpin * asComponents, sizeof(double),
         as_hessian_times_xc_deriv_gga, 1,
         dev_asHessianTimesXCDerivGGA, 1,
         "Setting as_hessian_times_xc_deriv_gga in set_as_hessian_times_xc_deriv_gga"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_hessian_times_xc_deriv_mgga_gpu)(
      double* as_hessian_times_xc_deriv_mgga)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         nMaxComputeHam * nMaxBatchSize * nSpin * asComponents, sizeof(double),
         as_hessian_times_xc_deriv_mgga, 1,
         dev_asHessianTimesXCDerivMGGA, 1,
         "Setting as_hessian_times_xc_deriv_mgga in set_as_hessian_times_xc_deriv_mgga"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_strain_deriv_kinetic_wave_gpu)(
      double* as_strain_deriv_kinetic_wave)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         nMaxComputeHam * nMaxBatchSize * asComponents, sizeof(double),
         as_strain_deriv_kinetic_wave, 1,
         dev_asStrainDerivKineticWave, 1,
         "Setting as_strain_deriv_kinetic_wave in set_as_strain_deriv_kinetic_wave"));
   CHECK_FOR_ERROR();
}

void FORTRAN(set_as_strain_deriv_wave_shell_gpu)(
      double* as_strain_deriv_wave_shell)
{
   // This subroutine is intended for debugging specific subroutine calls
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsSetVector(
         nMaxComputeHam * nMaxComputeHam, sizeof(double),
         as_strain_deriv_wave_shell, 1,
         dev_asStrainDerivWaveShell, 1,
         "Setting as_strain_deriv_wave_shell in set_as_strain_deriv_wave_shell"));
   CHECK_FOR_ERROR();
}

// Getter subroutines (transferring data from the GPU to the CPU)

void FORTRAN(get_as_strain_deriv_wave_shell_gpu)(
      double* as_strain_deriv_wave_shell)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsGetVector(
         nMaxComputeHam * nMaxComputeHam, sizeof(double),
         dev_asStrainDerivWaveShell, 1,
         as_strain_deriv_wave_shell, 1,
         "Getting as_strain_deriv_wave_shell in get_as_strain_deriv_wave_shell"));
   CHECK_FOR_ERROR();
}

void FORTRAN(get_forces_shell_gpu)(
      double* forces_shell)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsGetVector(
         nMaxComputeHam * nMaxComputeHam, sizeof(double),
         dev_forcesShell, 1,
         forces_shell, 1,
         "Getting forces_shell in get_sum_forces_gpu"));
   CHECK_FOR_ERROR();
}

void FORTRAN(get_sum_forces_gpu)(
      double* sum_forces,
      int*    n_atoms)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsGetVector(
         3 * *n_atoms, sizeof(double),
         dev_sumForces, 1,
         sum_forces, 1,
         "Getting sum_forces in get_sum_forces_gpu"));
   CHECK_FOR_ERROR();
}

void FORTRAN(get_as_pulay_stress_local_gpu)(
      double* as_pulay_stress_local)
{
   using namespace gpuForces;

   HANDLE_CUBLAS(aimsGetVector(
         asComponents, sizeof(double),
         dev_asPulayStressLocal, 1,
         as_pulay_stress_local, 1,
         "Getting as_pulay_stress_local in get_as_pulay_stress_local_gpu"));
   CHECK_FOR_ERROR();
}

////////////////////////////////////////////////////////////////////////////////
//                            Computation, Forces                             //
////////////////////////////////////////////////////////////////////////////////

// Evaluate the Pulay forces for the current batch
void FORTRAN(eval_forces_shell_dpsi_h_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_points > 0);
   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_a > 0);
   CHECK(*i_coord > 0);
   CHECK(*i_coord < 4);
   CHECK(*i_spin > 0);

   int dWaveOffset = (*i_coord-1) * nMaxComputeHam * nMaxBatchSize;
   int hTimesPsiOffset = (*i_spin-1) * nMaxComputeHam * nMaxBatchSize;

   // Multiply wave with partition grid weights
   // wave x DIAG(partition)
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute_a, *n_points,
         dev_dWave + dWaveOffset, nMaxComputeHam,
         dev_partition, 1,
         dev_waveComputeA, *n_compute_a));

   double alpha = 1.0;
   double beta = 0.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_c, *n_compute_a, *n_points,
         &alpha,
         dev_hTimesPsi + hTimesPsiOffset, nMaxComputeHam,
         dev_waveComputeA, *n_compute_a,
         &beta,
         dev_forcesShell, *n_compute_c));
}

// Evaluate the correction to the forces due to atomic ZORA for the current
// batch
void FORTRAN(eval_forces_shell_psi_dh_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_points > 0);
   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_a > 0);
   CHECK(*i_coord > 0);
   CHECK(*i_coord < 4);
   CHECK(*i_spin > 0);

   // Multiply wave with partition grid weights
   // wave x DIAG(partition)
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute_a, *n_points,
         dev_wave, nMaxComputeHam,
         dev_partition, 1,
         dev_waveComputeA, *n_compute_a));

   int dHTimesPsiOffset = (*i_coord-1) * nMaxComputeHam * nMaxBatchSize;

   double alpha = 1.0;
   double beta = 0.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_a, *n_compute_c, *n_points,
         &alpha,
         dev_waveComputeA, *n_compute_a,
         dev_dHTimesPsi + dHTimesPsiOffset, nMaxComputeHam,
         &beta,
         dev_forcesShell, *n_compute_a));
}

// Evaluate the GGA and meta-GGA contribution to forces for the current batch
void FORTRAN(eval_gga_forces_dens_mat_gpu)(
      int* n_compute,
      int* n_points,
      int* i_dim,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_compute > 0);
   CHECK(*n_points > 0);
   CHECK(*i_dim > 0);
   CHECK(*i_spin > 0);
   CHECK(*i_spin <= nSpin);

   int iIndex[3][3];
   int iCounter = -1;
   int offset = -1;

   // iIndex is the indexing array from the nine overall components of the
   // Hessian to the six independent components (due to symmetry), that is,
   // 1 2 3
   // 2 4 5
   // 3 5 6
   for (int iCoord1 = 0; iCoord1 < 3; iCoord1++)  {
      for (int iCoord2 = iCoord1; iCoord2 < 3; iCoord2++)  {
         iCounter++;
         iIndex[iCoord1][iCoord2] = iCounter;
         iIndex[iCoord2][iCoord1] = iCounter;
      }
   }

   offset = (*i_dim-1) * *n_compute;

   // We omit the factor of 2.0 in dev_matrixTmp2 here; instead we move it into
   // the subsequent cublasDgemm
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute, *n_points,
         dev_gradientBasisWave + offset, 3 * nMaxComputeHam,
         dev_partition, 1,
         dev_matrixTmp2, *n_compute));

   for (int iCoord = 0; iCoord < 3; iCoord++)   {
      offset = iCoord * *n_compute;
      int xcOffset = (*i_spin - 1) * 3 + iCoord;

      HANDLE_CUBLAS(cublasDdgmm(
            cublasHandle,
            CUBLAS_SIDE_RIGHT,
            *n_compute, *n_points,
            dev_gradientBasisWave + offset, 3 * nMaxComputeHam,
            dev_xcGradientDeriv + xcOffset, 3 * nSpin,
            dev_matrixTmp1, *n_compute));

      // beta = 1.0 here, as it is assumed that this will be called after
      // calculating Pulay forces, and thus dev_forcesShells is valid data
      double alpha = 2.0;
      double beta = 1.0;
      HANDLE_CUBLAS(cublasDgemm(
            cublasHandle,
            CUBLAS_OP_N, CUBLAS_OP_T,
            *n_compute, *n_compute, *n_points,
            &alpha,
            dev_matrixTmp1, *n_compute,
            dev_matrixTmp2, *n_compute,
            &beta,
            dev_forcesShell, *n_compute));
   }

   // We omit the factor of 2.0 in dev_matrixTmp1 here; instead we move it into
   // the subsequent cublasDgemm
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute, *n_points,
         dev_wave, nMaxComputeHam,
         dev_partition, 1,
         dev_matrixTmp1, *n_compute));

   for (int iCoord = 0; iCoord < 3; iCoord++)   {
      offset = iIndex[*i_dim-1][iCoord] * nMaxComputeHam;
      int xcOffset = (*i_spin - 1) * 3 + iCoord;
      HANDLE_CUBLAS(cublasDdgmm(
            cublasHandle,
            CUBLAS_SIDE_RIGHT,
            *n_compute, *n_points,
            dev_hessianBasisWave + offset, 6 * nMaxComputeHam,
            dev_xcGradientDeriv + xcOffset, 3 * nSpin,
            dev_matrixTmp2, *n_compute));

      if (metaGGAForcesOn) {
         offset = iCoord * *n_compute;
         HANDLE_CUBLAS(cublasDdgmm(
               cublasHandle,
               CUBLAS_SIDE_RIGHT,
               *n_compute, *n_points,
               dev_gradientBasisWave + offset, 3 * nMaxComputeHam,
               dev_partition, 1,
               dev_matrixTmp1MGGA, *n_compute));

         offset = iIndex[*i_dim-1][iCoord] * nMaxComputeHam;
         xcOffset = (*i_spin - 1);
         HANDLE_CUBLAS(cublasDdgmm(
               cublasHandle,
               CUBLAS_SIDE_RIGHT,
               *n_compute, *n_points,
               dev_hessianBasisWave + offset, 6 * nMaxComputeHam,
               dev_xcTauDeriv + xcOffset, nSpin,
               dev_matrixTmp2MGGA, *n_compute));
      }

      double alpha = 2.0;
      double beta = 1.0;
      HANDLE_CUBLAS(cublasDgemm(
            cublasHandle,
            CUBLAS_OP_N, CUBLAS_OP_T,
            *n_compute, *n_compute, *n_points,
            &alpha,
            dev_matrixTmp1, *n_compute,
            dev_matrixTmp2, *n_compute,
            &beta,
            dev_forcesShell, *n_compute));

      if (metaGGAForcesOn) {
         alpha = 1.0; // The meta-GGA matrices have no missing factor of 2.0
                      // to account for here, hence a 1.0
         beta = 1.0;
         HANDLE_CUBLAS(cublasDgemm(
               cublasHandle,
               CUBLAS_OP_N, CUBLAS_OP_T,
               *n_compute, *n_compute, *n_points,
               &alpha,
               dev_matrixTmp1MGGA, *n_compute,
               dev_matrixTmp2MGGA, *n_compute,
               &beta,
               dev_forcesShell, *n_compute));
      }
   }
}

// Calculates the forces from the force shell.  The version including the
// analytical stress tensor is as_update_sum_forces_and_stress_gpu.
void FORTRAN(update_sum_forces_gpu)(
      int* n_compute,
      int* i_calculate_dimension,
      int* i_spin,
      int* n_compute_for_atom)
{
   using namespace gpuForces;

   CHECK(*n_compute > 0);
   CHECK(*n_compute <= nMaxComputeHam);
   CHECK(*i_calculate_dimension > 0);
   CHECK(*i_calculate_dimension <= 3);
   CHECK(*i_spin > 0);
   CHECK(*i_spin <= nSpin);
   CHECK(dev_forceValues != NULL);
   CHECK(dev_densMat != NULL);
   CHECK(dev_theNumberOneHaHaHa != NULL);
   CHECK(dev_forcesShell != NULL);
   if (loadBalancedMatrix) {
      CHECK(dev_insIdx != NULL);
      CHECK(n_compute_for_atom != NULL);
   }
   CHECK(nBasis >= 0);

   int densMatOffset = (*i_spin - 1) * ldDensMat;

   dim3 threads = dim3(maxThreads/2);
   dim3 grids = dim3((*n_compute * *n_compute + threads.x - 1)/threads.x);

   // Multiply the forces shell by the density matrix elementwise.
   if (loadBalancedMatrix) {
      updateSumForcesLoadBalanced<<<grids,threads>>>(
            *n_compute,
            *i_calculate_dimension,
            dev_permuteComputeByAtom,
            dev_insIdx,
            dev_densMat + densMatOffset,
            dev_forcesShell,
            dev_forceValues);
      CHECK_FOR_ERROR();
   } else {
      printf ("Unsupported Matrix Mode for forces indexing on GPU.  Exiting.");
      AIMS_EXIT();
   }

   cublasPointerMode_t pMode;
   HANDLE_CUBLAS(cublasGetPointerMode(cublasHandle, &pMode));
   HANDLE_CUBLAS(cublasSetPointerMode(cublasHandle,CUBLAS_POINTER_MODE_DEVICE));

   // Now update the appropriate coordinate for the forces
   unsigned int offsetForceValues = 0;
   unsigned int nComputeCurrent = 0;
   for (unsigned int iAtom = 1; iAtom <= nAtoms; iAtom++)  {
      // Fortunately, as a finite number of atoms will touch a batch, this
      // algorithm will asympotically be independent of the total number of
      // atoms in the system
      if (n_compute_for_atom[iAtom-1] == 0) continue;
      nComputeCurrent = n_compute_for_atom[iAtom-1];

      int target = (*i_calculate_dimension - 1) + 3 * (iAtom - 1);
      // Sum over all matrix elements contributing to the current force
      // component.  We here exploit that the fact that we previously ordered
      // matrix elements by atom, so that the summation reduces to reduction
      // over an interrupted block of GPU memory.
      // See the AS tensor indexing for an explanation on usage of cublasDdot.
      HANDLE_CUBLAS(cublasDdot(
            cublasHandle,
            nComputeCurrent * *n_compute,
            dev_forceValues + offsetForceValues, 1,
            dev_theNumberOneHaHaHa, 0,
            dev_forceComponent));
      // Finally, add the results to the current value for the force component
      sumForcesAdd<<<1,1>>>(dev_sumForces + target, dev_forceComponent);
      CHECK_FOR_ERROR();
      offsetForceValues += nComputeCurrent * *n_compute;
   }
   HANDLE_CUBLAS(cublasSetPointerMode(cublasHandle,pMode));
}


////////////////////////////////////////////////////////////////////////////////
//                      Computation, Analytical Stress                        //
////////////////////////////////////////////////////////////////////////////////

// Virtually identical to eval_forces_shell_dpsi_h_psi_gpu, but for
// Pulay contributions to analytical stress tensor
void FORTRAN(eval_as_shell_dpsi_h_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_points > 0);
   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_a > 0);
   CHECK(*i_coord > 0);
   CHECK(*i_coord < 9);
   CHECK(*i_spin > 0);

   int dWaveOffset = (*i_coord-1) * nMaxComputeHam * nMaxBatchSize;
   int hTimesPsiOffset = (*i_spin-1) * nMaxComputeHam * nMaxBatchSize;

   // Multiply wave with partition grid weights
   // wave x DIAG(partition)
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute_a, *n_points,
         dev_asStrainDerivWave + dWaveOffset, nMaxComputeHam,
         dev_partition, 1,
         dev_waveComputeA, *n_compute_a));

   double alpha = 1.0;
   double beta = 0.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_c, *n_compute_a, *n_points,
         &alpha,
         dev_hTimesPsi + hTimesPsiOffset, nMaxComputeHam,
         dev_waveComputeA, *n_compute_a,
         &beta,
         dev_asStrainDerivWaveShell, *n_compute_c));
}

void FORTRAN(eval_as_shell_add_psi_kin_psi_shell_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_points > 0);
   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_a > 0);
   CHECK(*i_coord > 0);
   CHECK(*i_coord < 9);
   CHECK(*i_spin > 0);

   // Multiply wave with partition grid weights
   // wave x DIAG(partition)
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute_a, *n_points,
         dev_wave, nMaxComputeHam,
         dev_partition, 1,
         dev_waveComputeA, *n_compute_a));

   int asJacPotOffset = ((*i_coord-1) * nSpin + (*i_spin-1)) *
         nMaxComputeHam * nMaxBatchSize;

   double alpha = 0.5;
   double beta = 1.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_c, *n_compute_a, *n_points,
         &alpha,
         dev_asJacPotKinTimesPsi + asJacPotOffset, nMaxComputeHam,
         dev_waveComputeA, *n_compute_a,
         &beta,
         dev_asStrainDerivWaveShell, *n_compute_c));
}

// Evaluate the GGA (and part of the meta-GGA) contribution to the analytical
// stress tensor
void FORTRAN(as_evaluate_gga_stress_gpu)(
      int* n_compute_c,
      int* n_points,
      int* i_coord,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_c <= nMaxComputeHam);
   CHECK(*n_points > 0);
   CHECK(*i_coord > 0);
   CHECK(*i_coord <= 9);
   CHECK(*i_spin > 0);
   CHECK(*i_spin <= nSpin);

   // The prefactor of 2 that shows up when constructing the intermediate matrix
   // in the CPU version of this subroutine is carried over to the DGEMM here.

   // WPH:  There must be a better way to do this.  In the CPU code, we loop
   //       over the three coordinates, and matrix_term is a sum over the three
   //       coordinates.  One matrix multiply is then needed.  Here, since
   //       cublasDdgmm has no accumulate functionality, I have to do a matrix
   //       multiply after every construction of dev_waveComputeA (the GPU
   //       analogue of matrix_term.)  I'm thinking a custom kernel is better
   //       suited to this task rather than contorting cuBLAS.
   double alpha = 2.0;
   double beta = 1.0;
   for (unsigned int iCoord2 = 0; iCoord2 < 3; iCoord2++) {
      alpha = 2.0;
      unsigned int xcGradDerivOffset = 3*(*i_spin-1) + iCoord2;
      HANDLE_CUBLAS(cublasDdgmm(
            cublasHandle,
            CUBLAS_SIDE_RIGHT,
            1, *n_points,
            dev_xcGradientDeriv + xcGradDerivOffset, 3*nSpin,
            dev_partition, 1,
            dev_matrixTmp3, 1));
      unsigned int gradBasisWaveOffset = iCoord2*(*n_compute_c);
      HANDLE_CUBLAS(cublasDdgmm(
            cublasHandle,
            CUBLAS_SIDE_RIGHT,
            *n_compute_c, *n_points,
            dev_gradientBasisWave + gradBasisWaveOffset, nMaxComputeHam*3,
            dev_matrixTmp3, 1,
            dev_waveComputeA, *n_compute_c));
      int dWaveOffset = (*i_coord-1) * nMaxComputeHam * nMaxBatchSize;
      HANDLE_CUBLAS(cublasDgemm(
            cublasHandle,
            CUBLAS_OP_N, CUBLAS_OP_T,
            *n_compute_c, *n_compute_c, *n_points,
            &alpha,
            dev_asStrainDerivWave + dWaveOffset, nMaxComputeHam,
            dev_waveComputeA, *n_compute_c,
            &beta,
            dev_asStrainDerivWaveShell, *n_compute_c));

      if (metaGGAForcesOn) {
         // This is a bit redundant, since we could restructure the
         // previous GGA part to calculate this piece first, but my experience
         // is that cublasDdgmm calls of these sizes are essentially free
         alpha = 1.0;
         HANDLE_CUBLAS(cublasDdgmm(
               cublasHandle,
               CUBLAS_SIDE_RIGHT,
               *n_compute_c, *n_points,
               dev_gradientBasisWave + gradBasisWaveOffset, nMaxComputeHam*3,
               dev_partition, 1,
               dev_matrixTmp1, *n_compute_c));
         int asHessXCOffset =
              (*i_coord-1) * nMaxComputeHam * nMaxBatchSize * nSpin
              + (*i_spin-1) * nMaxComputeHam * nMaxBatchSize;
         HANDLE_CUBLAS(cublasDgemm(
               cublasHandle,
               CUBLAS_OP_N, CUBLAS_OP_T,
               *n_compute_c, *n_compute_c, *n_points,
               &alpha,
               dev_asHessianTimesXCDerivMGGA + asHessXCOffset, nMaxComputeHam,
               dev_matrixTmp1, *n_compute_c,
               &beta,
               dev_asStrainDerivWaveShell, *n_compute_c));
      }
   }

   // Multiply waves with partition grid weights
   // wave x DIAG(partition)
   // The prefactor of 2 in the CPU code is carried over to the next operation.
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute_c, *n_points,
         dev_wave, nMaxComputeHam,
         dev_partition, 1,
         dev_waveComputeA, *n_compute_c));

   alpha = 2.0;
   int asHessXCOffset = (*i_coord-1) * nMaxComputeHam * nMaxBatchSize * nSpin +
                        (*i_spin-1) * nMaxComputeHam * nMaxBatchSize;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_c, *n_compute_c, *n_points,
         &alpha,
         dev_asHessianTimesXCDerivGGA + asHessXCOffset, nMaxComputeHam,
         dev_waveComputeA, *n_compute_c,
         &beta,
         dev_asStrainDerivWaveShell, *n_compute_c));
}

void FORTRAN(evaluate_forces_shell_add_mgga_gpu)(
     int* n_points,
     int* n_compute_c,
     int* n_compute_a,
     int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_points > 0);
   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_c <= nMaxComputeHam);
   CHECK(*n_compute_a > 0);
   CHECK(*n_compute_a <= nMaxComputeHam);
   CHECK(*i_spin > 0);
   CHECK(*i_spin <= nSpin);

   // dev_matrixTmp1MGGA here serves the role of left_side_of_mgga_dot_product
   // from the CPU code
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         (*n_compute_c)*3, *n_points,
         dev_gradientBasisWave, nMaxComputeHam*3,
         dev_partition, 1,
         dev_matrixTmp2MGGA, (*n_compute_c)*3));
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         (*n_compute_c)*3, *n_points,
         dev_matrixTmp2MGGA, (*n_compute_c)*3,
         dev_xcTauDeriv + (*i_spin-1), nSpin,
         dev_matrixTmp1MGGA, (*n_compute_c)*3));

   // dev_matrixTmp2MGGA here serves the role of gradient_basis_wave_compute_a
   // from the CPU code
   // Multiplying by the identity matrix may seem silly, but it serves a point;
   // gradient_basis_wave is allocated as an
   // (n_max_ham_compute, 3, n_max_batch_size) array, but it is computed
   // point-by-point, so the first two arguments are essentially fused together
   // to form a (n_max_ham_compute * 3, n_max_batch_size) in memory layout.
   // But here we want the opposite layout; we want a (n_compute_c, 3*n_points)
   // array so that we can use 3*n_points as the inner dimension in the matrix
   // multiplication.  So we reshape the matrix here by multiplying
   // gradient_basis_wave by the identity matrix, structuring the matrix
   // dimensions to get the desired contiguous layout in memory.
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         (*n_compute_a)*3, (*n_points),
         dev_gradientBasisWave, nMaxComputeHam*3,
         dev_theNumberOneHaHaHa, 0,
         dev_matrixTmp2MGGA, (*n_compute_a)*3));

   double alpha = -0.5;
   double beta = 1.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_c, *n_compute_a, 3*(*n_points),
         &alpha,
         dev_matrixTmp1MGGA, *n_compute_c,
         dev_matrixTmp2MGGA, *n_compute_a,
         &beta,
         dev_asStrainDerivWaveShell, *n_compute_c));
}

// Evaluate the correction to the analytical stress due to atomic ZORA for the
// current batch
void FORTRAN(eval_as_shell_psi_dkin_psi_gpu)(
      int* n_points,
      int* n_compute_c,
      int* n_compute_a,
      int* i_coord,
      int* i_spin)
{
   using namespace gpuForces;

   CHECK(*n_points > 0);
   CHECK(*n_compute_c > 0);
   CHECK(*n_compute_a > 0);
   CHECK(*i_coord > 0);
   CHECK(*i_coord < 9);
   CHECK(*i_spin > 0);

   // Multiply wave with partition grid weights
   // wave x DIAG(partition)
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute_a, *n_points,
         dev_wave, nMaxComputeHam,
         dev_partition, 1,
         dev_waveComputeA, *n_compute_a));

   int dKinPsiOffset = (*i_coord-1) * nMaxComputeHam * nMaxBatchSize;

   double alpha = 1.0;
   double beta = 0.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_c, *n_compute_a, *n_points,
         &alpha,
         dev_waveComputeA, *n_compute_a,
         dev_asStrainDerivKineticWave + dKinPsiOffset, nMaxComputeHam,
         &beta,
         dev_asStrainDerivWaveShell, *n_compute_a));
}

// Transpose the analytical stress shell matrix.
void FORTRAN(transpose_as_shell_gpu)(
      int* n_compute_c)
{
   using namespace gpuForces;

   CHECK(*n_compute_c > 0);

   double alpha = 1.0;
   double beta = 0.0;
   HANDLE_CUBLAS(cublasDgeam(
         cublasHandle,
         CUBLAS_OP_T, CUBLAS_OP_N,
         *n_compute_c, *n_compute_c,
         &alpha,
         dev_asStrainDerivWaveShell, *n_compute_c,
         &beta,
         NULL, *n_compute_c,
         dev_asStrainDerivWaveShellTrans, *n_compute_c));

   HANDLE_CUBLAS(cublasDcopy(cublasHandle,
         (*n_compute_c)*(*n_compute_c),
         dev_asStrainDerivWaveShellTrans, 1,
         dev_asStrainDerivWaveShell, 1));
}

// Calculates the forces and stress tensor components from the force and stress
// tensor shells, respectively.  The forces-only version is
// update_sum_forces_gpu
void FORTRAN(as_update_sum_forces_and_stress_gpu)(
      int* n_compute,
      int* i_calculate_dimension,
      int* i_spin,
      int* n_compute_for_atom)
{
   using namespace gpuForces;

   CHECK(*n_compute > 0);
   CHECK(*n_compute <= nMaxComputeHam);
   CHECK(*i_calculate_dimension > 0);
   CHECK(*i_calculate_dimension <= 9);
   CHECK(*i_spin > 0);
   CHECK(*i_spin <= nSpin);
   CHECK(dev_forceValues != NULL);
   CHECK(dev_densMat != NULL);
   CHECK(dev_theNumberOneHaHaHa != NULL);
   CHECK(dev_forcesShell != NULL);
   if (loadBalancedMatrix) {
      CHECK(dev_insIdx != NULL);
   }
   CHECK(nBasis >= 0);

   int densMatOffset = (*i_spin - 1) * ldDensMat;

   dim3 threads = dim3(maxThreads/2);
   dim3 grids = dim3((*n_compute * *n_compute + threads.x - 1)/threads.x);
   cublasPointerMode_t pMode;

   // Multiply the analytical stress shell by the density matrix elementwise.
   // I've omitted the non-packed version here, because density/force/AS update
   // via the density matrix required a packed matrix (or load balancing) be
   // used.
   if (loadBalancedMatrix) {
      asUpdateStressLoadBalanced<<<grids,threads>>>(
            *n_compute,
            *i_calculate_dimension,
            dev_insIdx,
            dev_densMat + densMatOffset,
            dev_asStrainDerivWaveShell,
            dev_forceValues);
      CHECK_FOR_ERROR();
   } else {
      printf ("Unsupported Matrix Mode for stress tensor indexing on GPU.  Exiting.");
      AIMS_EXIT();
   }

   HANDLE_CUBLAS(cublasGetPointerMode(cublasHandle, &pMode));
   HANDLE_CUBLAS(cublasSetPointerMode(cublasHandle,CUBLAS_POINTER_MODE_DEVICE));
   // Now sum up all elements calculated in the previous step to generate the
   // desired analytical stress tensor element.
   // We here use a trick to compute the reduction that a dot product, when one
   // of the vectors is set to one for all elements, is equivalent to a reduction.
   // This allows us to use cublasDdot to perform a reduction.  I've (WPH) done
   // some light testing of this trick versus the well-optimized reduction example
   // in the NVIDIA SDK, and the timings are virtually identical;  it's very likely
   // that NVIDIA is using the same algorithm for both implementations, since
   // reductions and dot products are so similar.
   // We furthermore exploit the fact we're allowed to define vectors with zero
   // stride to cut down on memory usage; our second "vector" is a single floating
   // point number which is repeatedly accessed.
   int target = *i_calculate_dimension - 1;
   HANDLE_CUBLAS(cublasDdot(
         cublasHandle,
         *n_compute * *n_compute,
         dev_forceValues, 1,
         dev_theNumberOneHaHaHa, 0,
         dev_forceComponent));
   // Finally, add the results to the current value for the force component
   sumForcesAdd<<<1,1>>>(dev_asPulayStressLocal + target, dev_forceComponent);
   CHECK_FOR_ERROR();
   HANDLE_CUBLAS(cublasSetPointerMode(cublasHandle,pMode));

   // Force components only range from 1..3, so exit if the component is higher
   if (*i_calculate_dimension > 3) return;

   FORTRAN(update_sum_forces_gpu)(
         n_compute,
         i_calculate_dimension,
         i_spin,
         n_compute_for_atom);
}

/*******************************************************************************
**                               CUDA Kernels                                 **
*******************************************************************************/

// Add the contribution to the forces from the current batch to the force
// component
// Written as a kernel to keep it on the GPU, but only executed by one thread
// (hence the absense of the usual thread/block id's)
__global__ void sumForcesAdd(
      double* address,
      double* value)
{
   *address += *value;
}

// Perform the indexing of forces from the batch matrix to the relevant
// "full" matrix when load balancing is being used
__global__ void updateSumForcesLoadBalanced(
      int nCompute,
      int iCoord,
      int* permuteComputeByAtom,
      int* insIdx,
      double* densMat,
      double* forcesShell,
      double* forceValues)
{

   int iOffset    = -1;
   int iIndexReal = -1;
   int newElement = -1;

   int element = threadIdx.x + blockIdx.x * blockDim.x;
   while (element < nCompute * nCompute)  {
      int iCompute2 = element % nCompute;
      int iCompute1 = (element - iCompute2) / nCompute;
      const int iPermBasis1 = __ldg(permuteComputeByAtom + iCompute1);
      const int iPermBasis2 = __ldg(permuteComputeByAtom + iCompute2);
      const int iBasisLoc1 = __ldg(insIdx + iCompute1);
      const int iBasisLoc2 = __ldg(insIdx + iCompute2);
      if (iBasisLoc1 <= iBasisLoc2) {
         iOffset = iBasisLoc2 * (iBasisLoc2-1)/2;
         iIndexReal = iOffset + iBasisLoc1;
      } else {
         iOffset = iBasisLoc1 * (iBasisLoc1-1)/2;
         iIndexReal = iOffset + iBasisLoc2;
      }
      const double densMatElem = __ldg(densMat + iIndexReal - 1);
      const double shellElem = __ldg(forcesShell + element);

      newElement = (iPermBasis2 - 1) + (iPermBasis1 - 1) * nCompute;
      forceValues[newElement] = densMatElem * shellElem;

      element += gridDim.x * blockDim.x;
   }
}

// Multiply the stress tensor shell by the associated element of the
// density matrix.  These quantities will be summed up later to generate
// the stress tensor component.
__global__ void asUpdateStressLoadBalanced(
      int nCompute,
      int iCoord,
      int* insIdx,
      double* densMat,
      double* asStrainDerivWaveShell,
      double* asStressValues)
{
   int iOffset    = -1;
   int iIndexReal = -1;

   int element = threadIdx.x + blockIdx.x * blockDim.x;
   while (element < nCompute * nCompute)  {
      asStressValues[element] = 0.0;
      int iCompute2 = element % nCompute;
      int iCompute1 = (element - iCompute2) / nCompute;
      const int iBasisLoc1 = __ldg(insIdx + iCompute1);
      const int iBasisLoc2 = __ldg(insIdx + iCompute2);
      if (iBasisLoc1 <= iBasisLoc2) {
         iOffset = iBasisLoc2 * (iBasisLoc2-1)/2;
         iIndexReal = iOffset + iBasisLoc1;
      } else {
         iOffset = iBasisLoc1 * (iBasisLoc1-1)/2;
         iIndexReal = iOffset + iBasisLoc2;
      }
      const double densMatElem = __ldg(densMat + iIndexReal - 1);
      const double shellElem = __ldg(asStrainDerivWaveShell + element);

      asStressValues[element] = densMatElem * shellElem;

      element += gridDim.x * blockDim.x;
   }
}
