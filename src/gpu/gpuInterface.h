#ifndef CUDA_INTERFACE_H
#define CUDA_INTERFACE_H

#include "gpuMacro.h"
#include <stdio.h>
#include <stdlib.h>

// global CUDA address array

extern int nDevices;

extern ssize_t maxMem;

extern ssize_t maxThreads;
extern ssize_t maxThreadsX;
extern ssize_t maxThreadsY;

// device specifications to perform parallel calculations
extern dim3 nBlocks;
extern dim3 nThreads;

// cublas handle
extern cublasHandle_t cublasHandle;

// cuda timing, times are in milliseconds
extern cudaEvent_t start, stop;

extern cublasStatus_t aimsSetVector(
      const ssize_t dim,
      const int elemSize,
      const void* from,
      const int ld1,
      void* to,
      const int ld2,
      const char* info);

extern cublasStatus_t aimsSetMatrix(
      const ssize_t rows,
      const ssize_t cols,
      const int elemSize,
      const void* from,
      const int ld1,
      void* to,
      const int ld2,
      const char* info);

extern cublasStatus_t aimsGetVector(
      const ssize_t dim,
      const int elemSize,
      const void* from,
      const int ld1,
      void* to,
      const int ld2,
      const char* info);

extern cublasStatus_t aimsGetMatrix(
      const ssize_t rows,
      const ssize_t cols,
      const int elemSize,
      const void* from,
      const int ld1,
      void* to,
      const int ld2,
      const char* info);

extern "C"
void FORTRAN(initialize_cuda_and_cublas)(
      int *dev_id);

extern "C"
void FORTRAN(set_gpu)(
      int *dev_id);

extern "C"
void FORTRAN(finalize_cuda_and_cublas)();

extern "C"
void FORTRAN(get_num_gpus)(
      int *result);

extern "C"
void FORTRAN(get_gpu_specs)(
      int *dev_id,
      int *major,
      int *minor,
      int *total_global_mem,
      int *total_constant_mem,
      int *shared_mem_per_block,
      int *clock_rate,
      int* multiprocessor_count,
      int *max_threads_per_multiprocessor,
      char *name);

#endif /*CUDA_INTERFACE_H*/
