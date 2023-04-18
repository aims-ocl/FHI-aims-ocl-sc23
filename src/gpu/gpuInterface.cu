#include <ctype.h>
#include <cublas_v2.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include "gpuError.h"
#include "gpuMacro.h"
#include <time.h>

#include "gpuInterface.h"

/*******************************************************************************
**                              Global Variables                              **
*******************************************************************************/

int nDevices = -1;

ssize_t maxThreads = -1;
ssize_t maxThreadsX = -1;
ssize_t maxThreadsY = -1;

// device specifications to perform parallel calculations
dim3 nBlocks;
dim3 nThreads;

// cublas handle
cublasHandle_t cublasHandle;

/*******************************************************************************
**                                CUDA CPU Code                               **
*******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                              Memory Management                             //
////////////////////////////////////////////////////////////////////////////////

cublasStatus_t aimsSetVector(
      const ssize_t dim,
      const int elemSize,
      const void* CPU,
      const int ld1,
      void* GPU,
      const int ld2,
      const char* info)
{
   cublasStatus_t error;
   CHECK_INFO(dim > 0, info);
   CHECK_INFO(elemSize > 0, info);
   CHECK_INFO(CPU != NULL, info);
   CHECK_INFO(GPU != NULL, info);
   CHECK_INFO(ld1 > 0, info);
   CHECK_INFO(ld2 > 0, info);

   error = cublasSetVector(dim, elemSize, CPU, ld1, GPU, ld2);
   return error;
}

cublasStatus_t aimsSetMatrix(
      const ssize_t rows,
      const ssize_t cols,
      const int elemSize,
      const void* CPU,
      const int ld1,
      void* GPU,
      const int ld2,
      const char* info)
{
   cublasStatus_t error;
   ssize_t dim = rows * cols;
   CHECK_INFO(dim > 0, info);
   CHECK_INFO(elemSize > 0, info);
   CHECK_INFO(CPU != NULL, info);
   CHECK_INFO(GPU != NULL, info);
   CHECK_INFO(ld1 > 0, info);
   CHECK_INFO(ld2 > 0, info);

   error = cublasSetMatrixAsync(rows, cols, elemSize, CPU, ld1, GPU, ld2, 0);
   return error;
}

cublasStatus_t aimsGetVector(
      const ssize_t dim,
      const int elemSize,
      const void* GPU,
      const int ld1,
      void* CPU,
      const int ld2,
      const char* info)
{
   cublasStatus_t error;
   CHECK_INFO(dim > 0, info);
   CHECK_INFO(elemSize > 0, info);
   CHECK_INFO(GPU != NULL, info);
   CHECK_INFO(CPU != NULL, info);
   CHECK_INFO(ld1 > 0, info);
   CHECK_INFO(ld2 > 0, info);

   error = cublasGetVector(dim, elemSize, GPU, ld1, CPU, ld2);
   return error;
}

cublasStatus_t aimsGetMatrix(
      const ssize_t rows,
      const ssize_t cols,
      const int elemSize,
      const void* GPU,
      const int ld1,
      void* CPU,
      const int ld2,
      const char* info)
{
   cublasStatus_t error;
   ssize_t dim = rows * cols;
   CHECK_INFO(dim > 0,info);
   CHECK_INFO(elemSize > 0, info);
   CHECK_INFO(GPU != NULL, info);
   CHECK_INFO(CPU != NULL, info);
   CHECK_INFO(ld1 > 0, info);
   CHECK_INFO(ld2 > 0, info);

   error = cublasGetMatrixAsync(rows, cols, elemSize, GPU, ld1, CPU, ld2, 0);
   return error;
}

////////////////////////////////////////////////////////////////////////////////
//                             Device Management                              //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(initialize_cuda_and_cublas)(
      int *dev_id)
{
   CHECK_VAR(*dev_id > -1, *dev_id);
   HANDLE_CUDA(cudaSetDevice(*dev_id));

   // explicitly clean up GPU
   HANDLE_CUDA(cudaDeviceReset());

   HANDLE_CUBLAS(cublasCreate(&cublasHandle));

   cudaDeviceProp *propPtr =
         (cudaDeviceProp*) malloc (sizeof(cudaDeviceProp));

   // get Properties
   HANDLE_CUDA(cudaGetDeviceProperties(propPtr,*dev_id));

   /*
   maxThreads = propPtr->maxThreadsDim[0];
   maxThreadsX = sqrt(propPtr->maxThreadsDim[0]);
   maxThreadsY = sqrt(propPtr->maxThreadsDim[0]);
   */
   maxThreads = 32;
   maxThreadsX = 32;
   maxThreadsY = 32;

   free (propPtr);
}

void FORTRAN(set_gpu)(
      int *dev_id)
{
   HANDLE_CUDA(cudaSetDevice (*dev_id));
}

void FORTRAN(finalize_cuda_and_cublas)()
{
   HANDLE_CUBLAS(cublasDestroy(cublasHandle));

   HANDLE_CUDA(cudaDeviceReset());
}

////////////////////////////////////////////////////////////////////////////////
//                          Device Information                                //
////////////////////////////////////////////////////////////////////////////////

void FORTRAN(get_num_gpus)(
      int *result)
{
   HANDLE_CUDA(cudaGetDeviceCount(result));
   nDevices = *result;
}

void FORTRAN(get_gpu_specs)(
      int *dev_id,
      int *major,
      int *minor,
      int *total_global_mem,
      int *total_constant_mem,
      int *shared_mem_per_block,
      int *clock_rate,
      int *multiprocessor_count,
      int *max_threads_per_multiprocessor,
      char* name)
{
   CHECK (dev_id != NULL);

   cudaDeviceProp *propPtr
         = (cudaDeviceProp*) malloc (sizeof(cudaDeviceProp));

   // get Properties
   HANDLE_CUDA(cudaGetDeviceProperties(propPtr,*dev_id));

   // map CUDA device properties to own struct
   *major = propPtr->major;
   *minor = propPtr->minor;
   *total_global_mem = static_cast<int>(propPtr->totalGlobalMem);
   *total_constant_mem = static_cast<int>(propPtr->totalConstMem);
   *shared_mem_per_block = static_cast<int>(propPtr->sharedMemPerBlock);
   *shared_mem_per_block = static_cast<int>(propPtr->sharedMemPerBlock);
   *multiprocessor_count = propPtr->multiProcessorCount;
   *max_threads_per_multiprocessor = propPtr->maxThreadsPerMultiProcessor;
   *clock_rate = static_cast<int>(propPtr->clockRate);
   for (int i = 0; i < 256; i++) name[i] = propPtr->name[i];

   free(propPtr);
}
