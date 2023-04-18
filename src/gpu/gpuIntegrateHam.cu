#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "gpuError.h"
#include "gpuInterface.h"
#include "gpuMacro.h"

#include "gpuIntegrateHam.h"

// Global Variables
namespace gpuIntegrateHam {
      double* dev_hamiltonianShell;
      double* dev_wave;
      double* dev_partition;
      double* dev_hTimesPsi;
      // For meta-GGA calculations
      double* dev_leftSideOfMGGADotProduct;
      double* dev_gradientBasisWaveStore;
      // For applying ZORA
      double* dev_zoraVector1;
      double* dev_zoraVector2;
      // For indexing the Hamiltonian on the GPU
      double* dev_hamiltonian;
      int* dev_map;
};

/*******************************************************************************
**                               CUDA CPU Code                                **
*******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                            Initialize/Finalize                             //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(hamiltonian_create_gpu)(
      int* n_max_compute_ham,
      int* n_max_batch_size,
      int* ld_hamiltonian,
      int* n_spin,
      int* use_meta_gga, // zero-if-false
      int* use_ZORA, // zero-if-false
      int* index_on_gpu) // zero-if-false
{
   using namespace gpuIntegrateHam;

   // Arrays used for evaluating the integrand of the Hamiltonian integration
   // on the current batch
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_hamiltonianShell,
         (*n_max_compute_ham) * (*n_max_compute_ham) * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_wave,
         (*n_max_compute_ham) * (*n_max_batch_size) * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_partition,
         *n_max_batch_size * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_hTimesPsi,
         (*n_max_compute_ham) * (*n_max_batch_size) * sizeof(double)));

   // Arrays used in meta-GGA calculations
   if (*use_meta_gga != 0) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_leftSideOfMGGADotProduct,
            (*n_max_compute_ham) * 3 * (*n_max_batch_size) * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_gradientBasisWaveStore,
            (*n_max_compute_ham) * 3 * (*n_max_batch_size) * sizeof(double)));
   }

   // Arrays used for applying ZORA
   if (*use_ZORA != 0) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_zoraVector1,
            (*n_max_compute_ham) * 3 * (*n_max_batch_size) * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_zoraVector2,
            (*n_max_compute_ham) * 3 * (*n_max_batch_size) * sizeof(double)));
   }

   // Arrays used when indexing of the batch matrices back into the Hamiltonian
   // is done on the GPU
   if (*index_on_gpu != 0) {
      // dev_map should have the same size as dev_hamiltonianShell
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_map,
            (*n_max_compute_ham) * (*n_max_compute_ham) * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_hamiltonian,
            (*ld_hamiltonian) * (*n_spin) * sizeof(double)));
   }
}

void FORTRAN(hamiltonian_destroy_gpu)()
{
   using namespace gpuIntegrateHam;

   HANDLE_CUDA(cudaFree(dev_hamiltonianShell));
   HANDLE_CUDA(cudaFree(dev_wave));
   HANDLE_CUDA(cudaFree(dev_partition));
   HANDLE_CUDA(cudaFree(dev_hTimesPsi));

   // Meta-GGA-related quantities
   HANDLE_CUDA(cudaFree(dev_leftSideOfMGGADotProduct));
   HANDLE_CUDA(cudaFree(dev_gradientBasisWaveStore));

   HANDLE_CUDA(cudaFree(dev_zoraVector1));
   HANDLE_CUDA(cudaFree(dev_zoraVector2));

   // Hamiltonian indexing quantities
   HANDLE_CUDA(cudaFree(dev_hamiltonian));
   HANDLE_CUDA(cudaFree(dev_map));
}

////////////////////////////////////////////////////////////////////////////////
//                               Data Movement                                //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(set_hamiltonian_gpu)(
      double *hamiltonian,
      int* dim_hamiltonian)
{
   using namespace gpuIntegrateHam;

   HANDLE_CUBLAS(aimsSetVector(
         *dim_hamiltonian, sizeof(double),
         hamiltonian, 1,
         dev_hamiltonian, 1,
         "Setting dev_hamiltonian in set_hamiltonian_gpu"));
}

void FORTRAN(get_hamiltonian_shell_gpu)(
      double *hamiltonian_shell, // matrix (n_compute, n_compute)
      int *n_compute) // # non zero elements
{
   using namespace gpuIntegrateHam;

   // copy zMatrix to CPU
   HANDLE_CUBLAS(aimsGetMatrix(
         *n_compute, *n_compute, sizeof(double),
         dev_hamiltonianShell, *n_compute,
         hamiltonian_shell,    *n_compute,
         "Getting hamiltonianShell in get_hamiltonian_shell_gpu"));
}

void FORTRAN(get_hamiltonian_gpu)(
      double *hamiltonian,
      int* dim_hamiltonian)
{
   using namespace gpuIntegrateHam;

   HANDLE_CUBLAS(aimsGetVector(
         *dim_hamiltonian, sizeof(double),
         dev_hamiltonian, 1,
         hamiltonian, 1,
         "Getting hamiltonian in get_hamiltonian_gpu"));
}

////////////////////////////////////////////////////////////////////////////////
//                  Computation, Hamiltonian Integration                      //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(evaluate_hamiltonian_shell_gpu) (
      int* n_points, // # integration points
      const double* partition, // array (nPoints)
      int* n_compute, // # non zero elements
      const double* h_times_psi, // matrix (nBasisList,nPoints)
      int* n_basis_list, // # basis functions
      double* wave, // matrix (nBasisList,nPoints)
      double* hamiltonian_shell) // matrix (n_compute, n_compute)
//  INPUTS
//   o  n_points -- number of grid points in this grid batch
//   o  n_compute -- number of non-zero basis functions in this grid batch
//   o  n_basis_list -- the total number  of basis functions
//   o  partition  -- values of partition function in this grid batch
//   o  wave -- values of basis functions in this grid batch
//   o  h_times_psi -- hamiltonian times basis functions in this grid batch
//  OUTPUT
//   o  hamiltonian_shell -- ( wave * hamiltonian * wave)
{
   using namespace gpuIntegrateHam;
   // Copy Partition Grid factors to GPU
   HANDLE_CUBLAS(aimsSetVector(
         *n_points, sizeof(double),
         partition, 1,
         dev_partition, 1,
         "Setting dev_partition in evaluate_hamiltonian_shell_gpu"));

   // Copy non-zero wave to GPU
   HANDLE_CUBLAS(aimsSetMatrix(
         *n_compute, *n_points, sizeof(double),
         wave, *n_basis_list,
         dev_wave, *n_compute,
         "Setting dev_wave in evaluate_hamiltonian_shell_gpu"));

   // Copy H|Psi> into memory of dev_partition
   HANDLE_CUBLAS(aimsSetMatrix(
         *n_basis_list, *n_points, sizeof(double),
         h_times_psi, *n_basis_list,
         dev_hTimesPsi, *n_basis_list,
         "Setting dev_hTimesPsi in evaluate_hamiltonian_shell_gpu"));

   // Multiply wave with partition grid weights
   // wave x DIAG(partition)
   HANDLE_CUBLAS(cublasDdgmm(
         cublasHandle,
         CUBLAS_SIDE_RIGHT,
         *n_compute, *n_points,
         dev_wave, *n_compute,
         dev_partition, 1,
         dev_wave, *n_compute));

   // Evaluate hamiltonian_shell = 0.5 * (H_times_psi * wave_compute**T
   //                            + wave_compute*H_times_psi**T)
   double alpha = 0.5;
   double beta = 0.0;
   HANDLE_CUBLAS(cublasDsyr2k(
         cublasHandle,
         CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N,
         *n_compute, *n_points,
         &alpha,
         dev_hTimesPsi, *n_basis_list,
         dev_wave, *n_compute,
         &beta,
         dev_hamiltonianShell, *n_compute));
}

// Name was shortened from
// evaluate_mgga_contribution_and_add_to_hamiltonian_shell_gpu_ due to an XL
// bug
void FORTRAN(mgga_contribution_gpu)(
      int* n_compute_1,
      int* n_compute_2,
      int* n_points,
      double* left_side_of_mgga_dot_product, // matrix (n_compute_1,3*n_points)
      double* gradient_basis_wave_store) // matrix (n_compute_2,3*n_points)
{
   using namespace gpuIntegrateHam;

   double alpha = 1.0;
   double beta = 1.0;

   HANDLE_CUBLAS(aimsSetMatrix(
         *n_compute_1, 3 * (*n_points), sizeof(double),
         left_side_of_mgga_dot_product, *n_compute_1,
         dev_leftSideOfMGGADotProduct, *n_compute_1,
         "Setting dev_leftSideOfMGGADotProduct in mgga_contribution_gpu"));

   HANDLE_CUBLAS(aimsSetMatrix(
         *n_compute_2, 3 * (*n_points), sizeof(double),
         gradient_basis_wave_store, *n_compute_2,
         dev_gradientBasisWaveStore, *n_compute_2,
         "Setting dev_gradientBasisWaveStore in mgga_contribution_gpu"));

   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute_1, *n_compute_2, 3 * *n_points,
         &alpha,
         dev_leftSideOfMGGADotProduct, *n_compute_1,
         dev_gradientBasisWaveStore, *n_compute_2,
         &beta,
         dev_hamiltonianShell, *n_compute_1));
}

void FORTRAN(add_zora_matrix_gpu)(
      double *zora_vector1, // 1. Zora tensor (nBasisList, 3, nRelPoints)
      double *zora_vector2, // 2. Zora tensor (nBasisList, 3, nRelPoints)
      int *n_basis_list, // # basis functions
      int *n_rel_points, // # relativity points
      int *n_compute) // # non zero elements
{
   using namespace gpuIntegrateHam;

   int elemSize = sizeof(double);

   // Copy zora to GPU
   HANDLE_CUBLAS(aimsSetMatrix(
         *n_basis_list, 3 * *n_rel_points, elemSize,
         zora_vector1, *n_basis_list,
         dev_zoraVector1, *n_basis_list,
         "Setting dev_zoraVector1 in add_zora_matrix_gpu"));
   HANDLE_CUBLAS(aimsSetMatrix(
         *n_basis_list, 3 * *n_rel_points, elemSize,
         zora_vector2, *n_basis_list,
         dev_zoraVector2, *n_basis_list,
         "Setting dev_zoraVector2 in add_zora_matrix_gpu"));

   double alpha = 1.0;
   double beta = 1.0;
   HANDLE_CUBLAS(cublasDgemm(
         cublasHandle,
         CUBLAS_OP_N, CUBLAS_OP_T,
         *n_compute, *n_compute, 3 * *n_rel_points,
         &alpha,
         dev_zoraVector1, *n_basis_list,
         dev_zoraVector2, *n_basis_list,
         &beta,
         dev_hamiltonianShell, *n_compute));
}

// This is the indexing used for non-packed matrices
void FORTRAN(update_full_matrix_via_map_gpu)(
      double* hamiltonian, // full storage in buffer 3 (1 D Object)
      int* dim_hamiltonian,
      double* hamiltonian_shell, // precalculated in buffer 2 (2 D Object)
      int* dim1_hamiltonian_shell, // number of relevant basis functions
      int* dim2_hamiltonian_shell, // number of relevant basis functions
      int* map, // array relevant basis functions
      int* dim_map) // dimension of map
{
   using namespace gpuIntegrateHam;

   // Copy to GPU
   HANDLE_CUBLAS(aimsSetVector(
         *dim_map, sizeof(int),
         map, 1,
         dev_map, 1,
         "Setting dev_map in update_full_matrix_via_map_gpu"));

   dim3 threads = dim3(maxThreadsX,maxThreadsY);
   dim3 grids = dim3((*dim1_hamiltonian_shell + threads.x - 1)/threads.x,
                     (*dim2_hamiltonian_shell + threads.y - 1)/threads.y);
   //dim3 grids = dim3(32,32);

   insertInHamiltonianViaMap<<<grids,threads,0,NULL>>>(
         dev_hamiltonian,
         *dim_hamiltonian,
         dev_hamiltonianShell,
         *dim1_hamiltonian_shell,
         *dim2_hamiltonian_shell,
         dev_map,
         *dim_map);
   CHECK_FOR_ERROR();
}

// This is the indexing used when load balancing is enabled
void FORTRAN(update_batch_matrix_gpu)(
      int* i_spin, // spin value
      int* ld_matrix, // leading dimension of matrix
      int* map, // index array for position in matrix
      int* n_compute_c) // number of nonzero basis functions
{
   using namespace gpuIntegrateHam;
   HANDLE_CUBLAS(aimsSetVector(*n_compute_c, sizeof(int),
         map, 1,
         dev_map, 1,
         "Setting dev_map in update_batch_matrix_gpu"));

   nThreads = dim3(maxThreadsX,maxThreadsY);
   nBlocks  = dim3((*n_compute_c + nThreads.x - 1)/nThreads.x,
                   (*n_compute_c + nThreads.y - 1)/nThreads.y);

   insertShellInMatrix<<<nBlocks,nThreads>>>(
         dev_hamiltonian,
         dev_hamiltonianShell,
         dev_map,
         *ld_matrix,
         *i_spin,
         *n_compute_c);
   CHECK_FOR_ERROR();
}

/*******************************************************************************
**                               CUDA Kernels                                 **
*******************************************************************************/

__global__ void insertInHamiltonianViaMap(
      double* hamiltonian,
      int dimHamiltonian,
      double* hamiltonianShell,
      int dim1HamiltonianShell,
      int dim2HamiltonianShell,
      int* map,
      int dimMap)
{
   int basis2 = threadIdx.y + blockIdx.y * blockDim.y;
   while (basis2 < dim2HamiltonianShell)  {
      int basis1 = threadIdx.x + blockIdx.x * blockDim.x;
      while (basis1 <= basis2)  {
         int element = basis1 + basis2 * dim1HamiltonianShell;
         int target = __ldg(map + element) - 1;
         hamiltonian[target] += hamiltonianShell[element];
         basis1 += blockDim.x * gridDim.x;
      }
      basis2 += blockDim.y * gridDim.y;
   }
}

__global__ void insertShellInMatrix(
      double* matrix,
      double* matrixShell,
      int* insIndex,
      int ldMatrix,
      int iSpin,
      int ldShell)
{
   int iCompute = threadIdx.x + blockIdx.x * blockDim.x;

   while (iCompute < ldShell)  {
      int insertIdx = insIndex[iCompute];
      int offset = (iSpin - 1) * ldMatrix + (insertIdx*(insertIdx-1))/2;

      int jCompute = threadIdx.y + blockIdx.y * blockDim.y;
      while (jCompute <= iCompute) {
         int matIdx   = insIndex[jCompute] + offset - 1;
         int shellIdx = jCompute + iCompute*ldShell;
         matrix[matIdx] += matrixShell[shellIdx];

         jCompute += blockDim.y * gridDim.y;
      }

      iCompute += blockDim.x * gridDim.x;
   }
}
