#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "gpuError.h"
#include "gpuInterface.h"
#include "gpuMacro.h"

#include "gpuDensity.h"

namespace gpuDensity {
       double* dev_rho;
       double* dev_deltaRho;
       double* dev_tempRho;
       double* dev_rhoChange; // This should not be a device pointer
       double* dev_partitionTab;
       double* dev_hartreePartitionTab;
       double* dev_wave;
       double* dev_densityMatrix;
       int* dev_batchPoint2FullPoint;
       double* dev_work;
       double* dev_resultMat;
       // For density gradients
       double* dev_rhoGradient;
       double* dev_deltaRhoGradient;
       double* dev_tempRhoGradient;
       double* dev_gradientBasisWave;
};

/*******************************************************************************
**                               CUDA CPU Code                                **
*******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                            Initialize/Finalize                             //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(density_create_gpu)(
      int* n_max_compute_dens,
      int* n_max_batch_size,
      int* n_full_points,
      int* n_spin,
      int* use_density_gradient)
{
   using namespace gpuDensity;

   // Input density arrays
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_partitionTab,
         *n_full_points * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_hartreePartitionTab,
         *n_full_points * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_rho,
         *n_full_points * *n_spin * sizeof(double)));

   // Output (and input) density arrays
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_deltaRho,
         *n_full_points * *n_spin * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_rhoChange,
         *n_spin * sizeof(double)));

   // Temporary density arrays
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_tempRho,
         *n_max_batch_size * sizeof(double)));
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_batchPoint2FullPoint,
         *n_max_batch_size * sizeof(int)));

   // Matrices used for calculation of density gradients (needed for GGA and
   // beyond)
   if (*use_density_gradient != 0) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_tempRhoGradient,
            3 * *n_max_batch_size * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_deltaRhoGradient,
            3* *n_full_points * *n_spin * sizeof(double)));
      // Input/output density gradient array
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_rhoGradient,
            3* *n_full_points * *n_spin * sizeof(double)));
   }

   // Arrays that are updated in every batch iteration
   // Note: For system with large vacuum regions, n_max_compute_den can be zero
   //       on certain tasks!  This means that the GPU (and CPU) arrays for
   //       wave, density_matrix, and gradient_basis_wave may be unallocated,
   //       which must be taken into account.
   if (*n_max_compute_dens > 0) {
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_wave,
            *n_max_compute_dens * *n_max_batch_size * sizeof(double)));
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_densityMatrix,
            *n_max_compute_dens * *n_max_compute_dens * sizeof(double)));
      if (*use_density_gradient != 0) {
         HANDLE_CUDA(cudaMalloc(
               (void**) &dev_gradientBasisWave,
               3 * *n_max_compute_dens * *n_max_batch_size * sizeof(double)));
      }
   }

   // Intermediate work matrices for CUBLAS calls
   // These work matrices replace the "buffer" system used in the original code
   // The factors of 3 are needed to use the same work matrices for both density
   // and density gradients calculations; there is wasted space here for
   // calculations not involving density gradients, but these matrices aren't
   // the memory bottleneck (dev_densityMatrix is) and most calculations use
   // density gradients via GGA (and beyond) functionals
   HANDLE_CUDA(cudaMalloc(
         (void**) &dev_resultMat,
         (3 * *n_max_batch_size) * (*n_max_batch_size) * sizeof(double)));
   if (*n_max_compute_dens > 0) {
      // Not needed if memory ever becomes are issue; can reuse
      // dev_gradientBasisWave
      HANDLE_CUDA(cudaMalloc(
            (void**) &dev_work,
            *n_max_compute_dens * (3 * *n_max_batch_size) * sizeof(double)));
   }

   // initialize vectors
   size_t dim;
   if (*n_max_compute_dens > 0) {
     dim = *n_max_compute_dens * *n_max_batch_size;
     nThreads = dim3(maxThreads);
     nBlocks = dim3((dim + nThreads.x - 1)/nThreads.x);
     gpuInitVectorDensity<<<nBlocks,nThreads>>>(dev_wave, 0.0, dim);
     CHECK_FOR_ERROR();

     dim = *n_max_compute_dens * *n_max_compute_dens;
     nThreads = dim3(maxThreads);
     nBlocks = dim3((dim + nThreads.x - 1)/nThreads.x);
     gpuInitVectorDensity<<<nBlocks,nThreads>>>(dev_densityMatrix, 0.0, dim);
     CHECK_FOR_ERROR();

     if (*use_density_gradient != 0) {
        dim = 3 * *n_max_compute_dens * *n_max_batch_size;
        nThreads = dim3(maxThreads);
        nBlocks = dim3((dim + nThreads.x - 1)/nThreads.x);
        gpuInitVectorDensity<<<nBlocks,nThreads>>>(
              dev_gradientBasisWave, 0.0, dim);
        CHECK_FOR_ERROR();
     }
   }
}

void FORTRAN(density_destroy_gpu)()
{
   using namespace gpuDensity;

   HANDLE_CUDA(cudaFree(dev_rho));
   HANDLE_CUDA(cudaFree(dev_deltaRho));
   HANDLE_CUDA(cudaFree(dev_tempRho));
   HANDLE_CUDA(cudaFree(dev_rhoChange));
   HANDLE_CUDA(cudaFree(dev_partitionTab));
   HANDLE_CUDA(cudaFree(dev_hartreePartitionTab));
   HANDLE_CUDA(cudaFree(dev_wave));
   HANDLE_CUDA(cudaFree(dev_densityMatrix));
   HANDLE_CUDA(cudaFree(dev_batchPoint2FullPoint));
   HANDLE_CUDA(cudaFree(dev_work));
   HANDLE_CUDA(cudaFree(dev_resultMat));

   HANDLE_CUDA(cudaFree(dev_gradientBasisWave));
   HANDLE_CUDA(cudaFree(dev_tempRhoGradient));
   HANDLE_CUDA(cudaFree(dev_deltaRhoGradient));
   HANDLE_CUDA(cudaFree(dev_rhoGradient));
}

////////////////////////////////////////////////////////////////////////////////
//                               Data Movement                                //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(set_delta_rho_gpu)(
      double* delta_rho,
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         *n_full_points * *n_spin, sizeof(double),
         delta_rho, 1,
         dev_deltaRho, 1,
         "Setting dev_deltaRho in set_delta_rho_gpu"));
}

void FORTRAN(set_rho_gpu)(
      double* rho,
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         *n_full_points * *n_spin, sizeof(double),
         rho, 1,
         dev_rho, 1,
         "Setting dev_rho in set_rho_gpu"));
}

void FORTRAN(set_rho_change_gpu)(
      double* rho_change,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         *n_spin, sizeof(double),
         rho_change, 1,
         dev_rhoChange, 1,
         "Setting dev_rhoChange in set_rho_change_gpu"));
}

void FORTRAN(set_partition_tab_gpu)(
      double* partition_tab,
      int* n_full_points)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         *n_full_points, sizeof(double),
         partition_tab, 1,
         dev_partitionTab, 1,
         "Setting dev_partitionTab in set_partition_tab_gpu"));
}

void FORTRAN(set_hartree_partition_tab_gpu)(
      double* hartree_partition_tab,
      int* n_full_points)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         *n_full_points, sizeof(double),
         hartree_partition_tab, 1,
         dev_hartreePartitionTab, 1,
         "Setting dev_hartreePartitionTab in set_hartree_partition_tab_gpu"));
}

void FORTRAN(set_delta_rho_gradient_gpu)(
      double* delta_rho_gradient,
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         3 * *n_full_points * *n_spin, sizeof(double),
         delta_rho_gradient, 1,
         dev_deltaRhoGradient, 1,
        "Setting dev_deltaRhoGradient in set_delta_rho_gradient_gpu"));
}

void FORTRAN(set_rho_gradient_gpu)(
      double* rho_gradient,
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsSetVector(
         3 * *n_full_points * *n_spin, sizeof(double),
         rho_gradient, 1,
         dev_rhoGradient, 1,
         "Setting dev_rhoGradient set_rho_gradient_gpu"));
}

void FORTRAN(set_density_matrix_gpu)(
      double* density_matrix, //matrix (n_compute,n_compute)
      int* n_compute)         //number of relevant basis functions
{
   using namespace gpuDensity;

   if (*n_compute > 0) {
      HANDLE_CUBLAS(aimsSetMatrix(
            *n_compute, *n_compute, sizeof(double),
            density_matrix, *n_compute,
            dev_densityMatrix, *n_compute,
            "Setting dev_densityMatrix in push_density_matrix"));
   }
}

void FORTRAN(get_delta_rho_gpu) (
      double* delta_rho,
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsGetVector(
         *n_full_points * *n_spin, sizeof(double),
         dev_deltaRho, 1,
         delta_rho, 1,
         "Getting delta_rho in get_delta_rho_gpu"));
}

void FORTRAN(get_delta_rho_gradient_gpu) (
      double* delta_rho_gradient,
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsGetVector(
         3 * *n_full_points * *n_spin, sizeof(double),
         dev_deltaRhoGradient, 1,
         delta_rho_gradient,1,
         "Getting delta_rho_gradient in get_delta_rho_gradient_gpu"));
}

void FORTRAN(get_rho_change_gpu) (
      double* rho_change,
      int* n_spin)
{
   using namespace gpuDensity;

   HANDLE_CUBLAS(aimsGetVector(
         *n_spin, sizeof(double),
         dev_rhoChange, 1,
         rho_change,1,
         "Getting rho_change in get_rho_change_gpu"));
}

////////////////////////////////////////////////////////////////////////////////
//                           Computation, Density                             //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(evaluate_ks_density_densmat_gpu) (
      int* n_points,    // number of integration points
      double* wave,    // matrix (n_basis_compute, n_points)
      int* n_compute,   // number of relevant basis functions
      int* n_basis_compute) //maximum number of relevant basis functions
{
   using namespace gpuDensity;
   if (*n_points > 0) {
      //Safeguard if nothing to do
      if (*n_compute == 0) {
         nThreads = dim3(maxThreads);
         nBlocks = dim3((*n_points + nThreads.x - 1)/nThreads.x);
         gpuInitVectorDensity<<<nBlocks,nThreads>>>(
               dev_tempRho, 0.0, *n_points);
         CHECK_FOR_ERROR();
         return;
      }
   }

   if (*n_points > 0 && *n_basis_compute > 0) {
      HANDLE_CUBLAS(aimsSetVector(
            *n_basis_compute * *n_points, sizeof(double),
            wave, 1,
            dev_wave, 1,
            "Setting dev_wave in evaluate_ks_density_densmat_gpu"));

      double alpha = 1.0;
      double beta = 0.0;
      HANDLE_CUBLAS(cublasDsymm(
            cublasHandle,
            CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER,
            *n_compute, *n_points,
            &alpha,
            dev_densityMatrix, *n_compute,
            dev_wave, *n_basis_compute,
            &beta,
            dev_work, *n_compute));

      alpha = 1.0;
      beta  = 0.0;
      HANDLE_CUBLAS(cublasDgemm(
            cublasHandle,
            CUBLAS_OP_T, CUBLAS_OP_N,
            *n_points, *n_points, *n_compute,
            &alpha,
            dev_work, *n_compute,
            dev_wave, *n_basis_compute,
            &beta,
            dev_resultMat, *n_points));

      nThreads = dim3(maxThreadsX, maxThreadsY);
      nBlocks  = dim3((*n_points + nThreads.x - 1) / nThreads.x,
                      (*n_points + nThreads.y - 1) / nThreads.y);
      gpuDiag<<<nBlocks,nThreads>>>(
            dev_tempRho, dev_resultMat, *n_points,*n_points);
      CHECK_FOR_ERROR();
   }
}

void FORTRAN(update_delta_rho_ks_gpu) (
      int* offset_full_points,  // starting point in density difference
      int* n_batch_points,      // number of Points in current Batch
      int* n_full_points,       // number of all Points
      int* n_points_compute,    // non zero points in rho
      int* i_spin,             // actual Spin Channel
      int* n_spin)             // number of Spin Channels
{
   using namespace gpuDensity;

   CHECK(dev_batchPoint2FullPoint != NULL);
   CHECK(dev_partitionTab != NULL);
   CHECK(dev_hartreePartitionTab != NULL);
   CHECK(dev_tempRho != NULL);
   CHECK(dev_rho != NULL);
   CHECK(dev_deltaRho != NULL);

   // Generate Mapping between iPoint index of local_rho in batch
   // and fullPoint index of delta rho
   generateBatchPoint2FullPoint<<<1,1>>>(
         *n_batch_points, *offset_full_points,
         dev_partitionTab, dev_hartreePartitionTab,
         dev_batchPoint2FullPoint);
   CHECK_FOR_ERROR();

   nThreads = dim3(maxThreads);
   if (*n_points_compute > 0) {
      nBlocks = dim3((*n_points_compute + nThreads.x - 1)/nThreads.x);

      distributeDensityUpdate<<<nBlocks,nThreads>>> (
            *n_points_compute,
            *n_full_points, *n_spin, *i_spin, dev_batchPoint2FullPoint,
            dev_tempRho, dev_rho, dev_deltaRho);
      CHECK_FOR_ERROR();
   }
}

void FORTRAN(calculate_rho_change_gpu) (
      int* n_full_points,
      int* n_spin)
{
   using namespace gpuDensity;

   nThreads = dim3(maxThreads);
   if (*n_full_points > 0) {
      nBlocks = dim3((*n_full_points + nThreads.x - 1)/nThreads.x);
      updateRhoChange<<<nBlocks,nThreads>>> (
            *n_full_points, *n_spin,
            dev_rhoChange, dev_deltaRho,
            dev_partitionTab);
      CHECK_FOR_ERROR();
   }
}


////////////////////////////////////////////////////////////////////////////////
//                     Computation, Density Gradients                         //
////////////////////////////////////////////////////////////////////////////////

// Note: The following functions should not be called if density
// gradients are not used, that is, use_density_gradient = .false.
// (i.e. "standard" LDA calculations)

// Name was shortened from evaluate_density_gradient_densmat_gpu_ due to an XL
// bug
void FORTRAN(eval_density_grad_densmat_gpu) (
      int* n_points, // number of grid points
      double* gradient_basis_wave, // gradients of basis functions
      int* n_compute, // number of relevant basis functions
      int* n_basis_compute) // total number of basis functions
{
   using namespace gpuDensity;

   if (*n_points > 0) {
      //Safeguard if nothing to do
      if (*n_compute == 0) {
         nThreads = dim3(maxThreads);
         nBlocks = dim3((*n_points + nThreads.x - 1)/nThreads.x);
         gpuInitVectorDensity<<<nBlocks,nThreads>>>(
               dev_tempRhoGradient, 0.0, *n_points * 3);
         CHECK_FOR_ERROR();
         return;
      }
   }

   if (*n_points > 0 && *n_basis_compute > 0) {
      HANDLE_CUBLAS(aimsSetVector(
            *n_basis_compute * 3 * *n_points, sizeof(double),
            gradient_basis_wave, 1,
            dev_gradientBasisWave, 1,
            "Setting dev_gradientBasisWave in eval_density_grad_densmat_gpu"));

      double alpha = 2.0;
      double beta = 0.0;
      HANDLE_CUBLAS(cublasDsymm(
            cublasHandle,
            CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER,
            *n_compute, *n_points * 3,
            &alpha,
            dev_densityMatrix, *n_compute,
            dev_gradientBasisWave, *n_basis_compute,
            &beta,
            dev_work, *n_compute));

      alpha = 1.0;
      beta  = 0.0;
      HANDLE_CUBLAS(cublasDgemm(
            cublasHandle,
            CUBLAS_OP_T, CUBLAS_OP_N,
            3 * *n_points, *n_points, *n_compute,
            &alpha,
            dev_work, *n_compute,
            dev_wave, *n_basis_compute,
            &beta,
            dev_resultMat, 3 * *n_points));

      nThreads = dim3(maxThreadsX, maxThreadsY);
      nBlocks  = dim3(3 * *n_points /nThreads.x+1, *n_points / nThreads.y +1);
      gpuDiag<<<nBlocks,nThreads>>>(
            dev_tempRhoGradient, dev_resultMat,
            3 * *n_points, *n_points);
      CHECK_FOR_ERROR();
   }
}

void FORTRAN(update_grad_delta_rho_ks_gpu) (
      int* n_full_points,       // number of all Points
      int* n_points_compute,    // non zero points in rho
      int* i_spin,             // actual Spin Channel
      int* n_spin)             // number of Spin Channels
{
   using namespace gpuDensity;

   nThreads = dim3(maxThreads);
   if (*n_points_compute > 0) {
      nBlocks = dim3((*n_points_compute + nThreads.x - 1)/nThreads.x);

      distributeGradDensityUpdate<<<nBlocks,nThreads>>>(
            *n_points_compute,
            *n_full_points, *n_spin, *i_spin, dev_batchPoint2FullPoint,
            dev_tempRhoGradient, dev_rhoGradient, dev_deltaRhoGradient);
      CHECK_FOR_ERROR();
   }
}

/*******************************************************************************
**                               CUDA Kernels                                 **
*******************************************************************************/

template<typename T>
__global__ void gpuInitVectorDensity(
      T* vector,
      const T val,
      const ssize_t dim)
{
   ssize_t elem = threadIdx.x + blockDim.x * blockIdx.x;
   ssize_t stride = blockDim.x * gridDim.x;

   while (elem < dim)  {
      vector[elem] = val;
      elem += stride;
   }
}

__device__ double myAtomicAddDensity(
      double* address,
      double val)
{
   unsigned long long int* addressAsULL =
                              (unsigned long long int*)address;
   unsigned long long int old = *addressAsULL, assumed;
   do {
      assumed = old;
      double result_d = val + __longlong_as_double(assumed);
      unsigned long long int resultULL = __double_as_longlong(result_d);
      old = atomicCAS(addressAsULL, assumed, resultULL);
   } while (assumed != old);
   return __longlong_as_double(old);
}

template<typename T>
__global__ void gpuDiag(
      T* result,
      const T* matrix,
      const int nRows,
      const int nCols)
{
   // This routine gets the diagonal elements from a matrix
   // With nElements the number of elements taken from the matrix can be
   // specified. This comes in handy if the the matrix is not square.
   int nElements = nRows / nCols;
   int col = threadIdx.y + blockDim.y * blockIdx.y;
   while (col < nCols)  {
      int row = threadIdx.x + blockDim.x * blockIdx.x;
      while (row < nRows)  {
         if ((row >= col * nElements) &&
               ((row - col * nElements) < nElements)) {
            int elem = row - col * nElements;
            result[col * nElements + elem] =
                  matrix[row + col * nElements * nCols];
         }
         row += blockDim.x * gridDim.x;
      }
      col += blockDim.y * gridDim.y;
   }
}

__global__ void generateBatchPoint2FullPoint(
      int nBatchPoints,
      int offset,
      double* partitionTab,
      double* hartreePartitionTab,
      int* batch2Full)
{
   int computePointIdx = 0;
   for (int i = 0; i < nBatchPoints; i++)  {
      batch2Full[i] = -1;
      int elem = offset + i;
      if (max(partitionTab[elem],hartreePartitionTab[elem]) > 0.0) {
         batch2Full[computePointIdx] = elem;
         computePointIdx++;
      }
   }
}

__global__ void distributeDensityUpdate(
      int nPoints,
      int nFullPoints,
      int nSpin,
      int iSpin,
      int* batch2Full,
      double* tempRho,
      double* rho,
      double* deltaRho)
{
   int elem = threadIdx.x + blockIdx.x * blockDim.x;

   while(elem < nPoints)  {
      const int iFullPoint = __ldg(batch2Full + elem);
      const double myRho = __ldg(tempRho + elem);
      if (iFullPoint >= 0 && iFullPoint < nFullPoints)  {
         if (nSpin == 1) {
            int targetDeltaRho = iFullPoint + (iSpin-1) * nFullPoints;
            int targetRho = (iSpin-1) + iFullPoint * nSpin;
            const double rhoVal = __ldg(rho + targetRho);
            myAtomicAddDensity(deltaRho+targetDeltaRho, myRho - rhoVal);
         }
         if (nSpin == 2) {
            int targetDeltaRho1 = iFullPoint;
            int targetDeltaRho2 = iFullPoint + nFullPoints;
            int targetRho1 = iFullPoint * nSpin;
            int targetRho2 = 1 + iFullPoint * nSpin;
            const double rhoVal1 = __ldg(rho + targetRho1);
            const double rhoVal2 = __ldg(rho + targetRho2);
            if (iSpin == 1) {
               myAtomicAddDensity(deltaRho+targetDeltaRho1, myRho - rhoVal1);
               myAtomicAddDensity(deltaRho+targetDeltaRho2, myRho - rhoVal1);
            } else { // iSpin = 2
               myAtomicAddDensity(deltaRho+targetDeltaRho1, myRho - rhoVal2);
               myAtomicAddDensity(deltaRho+targetDeltaRho2, rhoVal2 - myRho);
            }
         }
      }
      elem += blockDim.x * gridDim.x;
   }
}

__global__ void updateRhoChange(
      int nPoints,
      int nSpin,
      double* rhoChange,
      double* deltaRho,
      double* partitionTab)
{
   __shared__ double cache[1024] ;

   for (int iSpin = 0; iSpin < nSpin; iSpin++)  {
      double value = 0.0;

      int iPoint = threadIdx.x + blockIdx.x * blockDim.x;
      while (iPoint < nPoints)  {
         int rhoElem = iPoint + iSpin * nPoints;
         const double rho = __ldg(deltaRho + rhoElem);
         const double pTab = __ldg(partitionTab + iPoint);
         value += pTab * rho * rho;
         iPoint += blockDim.x * gridDim.x;
      }
      cache[threadIdx.x] = value;

      __syncthreads();

      int half = blockDim.x / 2;

      while (half != 0) {
         if (threadIdx.x < half)
               cache[threadIdx.x] += cache[threadIdx.x + half];
         half /= 2;
         __syncthreads();
      }

      if (threadIdx.x == 0) myAtomicAddDensity(rhoChange+iSpin,cache[0]);
   }
}

__global__ void distributeGradDensityUpdate(
      int nPoints,
      int nFullPoints,
      int nSpin,
      int iSpin,
      int* batch2Full,
      double* tempGradRho,
      double* gradRho,
      double* deltaGradRho)
{
   int elem = threadIdx.x + blockIdx.x * blockDim.x;

   while(elem < nPoints)  {
      int iFullPoint = batch2Full[elem];
      if (nSpin == 1) {
         for (int iCoord = 0; iCoord < 3; iCoord++)  {
            const double myGradRho = __ldg(tempGradRho + iCoord + 3 * elem);
            int targetDeltaRho = iCoord + 3 * iFullPoint;
            int targetRho = iCoord + 3 * iFullPoint * nSpin;
            const double gradRhoVal = __ldg(gradRho + targetRho);
            myAtomicAddDensity(
                  deltaGradRho+targetDeltaRho,
                  myGradRho - gradRhoVal);
         }
      }
      if (nSpin == 2) {
         for (int iCoord = 0; iCoord < 3; iCoord++)  {
            const double myGradRho = __ldg(tempGradRho + iCoord + 3 * elem);
            int targetDeltaRho1 = iCoord + 3 * iFullPoint;
            int targetDeltaRho2 = iCoord + 3 * (iFullPoint + nFullPoints);
            int targetRho1 = iCoord + 3 * iFullPoint * nSpin;
            int targetRho2 = iCoord + 3 * (1 + iFullPoint * nSpin);
            const double gradRhoVal1 = __ldg(gradRho + targetRho1);
            const double gradRhoVal2 = __ldg(gradRho + targetRho2);
            if (iSpin == 1) {
               myAtomicAddDensity(
                     deltaGradRho+targetDeltaRho1,
                     myGradRho - gradRhoVal1);
               myAtomicAddDensity(
                     deltaGradRho+targetDeltaRho2,
                     myGradRho - gradRhoVal1);
            } else { // iSpin = 2
               myAtomicAddDensity(
                     deltaGradRho+targetDeltaRho1,
                     myGradRho - gradRhoVal2);
               myAtomicAddDensity(
                     deltaGradRho+targetDeltaRho2,
                     gradRhoVal2 - myGradRho);
            }
         }
      }
      elem += blockDim.x * gridDim.x;
   }
}
