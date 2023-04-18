#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "gpuError.h"
#include "gpuInterface.h"
#include "gpuMacro.h"

#include "gpuMR.h"

namespace gpuMR {
    double *d_matrix_batch;
    double *d_matrix_batch_GIAO;
    double *d_r_mn_basis;
    double *d_wave;
    double *d_matrix_wave;
    double *d_matrix_packed;
    int *d_basis_glb_to_loc;
    int *d_i_basis;
    double **d_matrix_list[3]; // work array for cublasDgemmBatched
    int matrix_wave_size_factor;
};

void mr_initialize_gpu(const int &wave_size,
                       const int &matrix_wave_size,
                       const int &matrix_batch_size,
                       const int &matrix_batch_GIAO_size,
                       const int &r_mn_basis_size,
                       const int &matrix_packed_size,
                       const int &i_basis_size,
                       const int &basis_glb_to_loc_size,
                       const int *basis_glb_to_loc,
                       const int &matrix_wave_size_f)
{
    using namespace gpuMR;
    // Allocations
    cudaMalloc((void**) &d_wave, wave_size*sizeof(double));
    cudaMalloc((void**) &d_matrix_wave, matrix_wave_size*sizeof(double));
    cudaMalloc((void**) &d_matrix_batch, matrix_batch_size*sizeof(double));
    cudaMalloc((void**) &d_matrix_batch_GIAO,
               matrix_batch_GIAO_size*sizeof(double));
    cudaMalloc((void**) &d_r_mn_basis, r_mn_basis_size*sizeof(double));
    cudaMalloc((void**) &d_matrix_packed, matrix_packed_size*sizeof(double));
    cudaMalloc((void**) &d_i_basis, i_basis_size*sizeof(int));
    cudaMalloc((void**) &d_basis_glb_to_loc, basis_glb_to_loc_size*sizeof(int));
    // Initialize the packed matrix
    nThreads = dim3(maxThreadsX*maxThreadsY);
    nBlocks  = dim3(matrix_packed_size/nThreads.x);
    if (matrix_packed_size%nThreads.x) ++nBlocks.x;
    d_zero_matrix<<<nBlocks,nThreads>>>(d_matrix_packed, matrix_packed_size);
    // Copy basis_glb_to_loc to gpu
    cublasSetVector(basis_glb_to_loc_size, sizeof(int),
                    basis_glb_to_loc, 1, d_basis_glb_to_loc, 1);
    for (int i = 0; i < 3; ++i)
        cudaMalloc((void**) &d_matrix_list[i], 9*sizeof(double*));
    matrix_wave_size_factor = matrix_wave_size_f;
}

__global__
void d_zero_matrix(double *d_matrix, const int matrix_size)
{
    const int i = threadIdx.x + blockIdx.x*blockDim.x;
    if (i < matrix_size)
        d_matrix[i] = 0.0;
}

void evaluate_mr_batch(const int &n_points,
                       const int &n_compute,
                       const double *wave,
                       const double *matrix_wave,
                       const int &i_symmetrization,
                       const int &i_dir,
                       const int &max_dims)
{
    using namespace gpuMR;
    if (i_dir == 1) {
        cublasSetVector(n_compute*n_points, sizeof(double),
                        wave, 1, d_wave, 1);
        cublasSetVector(max_dims*n_compute*n_points, sizeof(double),
                        matrix_wave, 1, d_matrix_wave, 1);
    }
    const double alpha = 1.0;
    const double beta = 0.0;
    cublasDgemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                n_compute, n_compute, n_points,
                &alpha, d_wave, n_compute,
                d_matrix_wave+(i_dir-1)*n_compute*n_points, n_compute,
                &beta, d_matrix_batch, n_compute);
    nThreads = dim3(maxThreadsX, maxThreadsY);
    nBlocks  = dim3(n_compute/nThreads.x, n_compute/nThreads.y);
    if (n_compute%nThreads.x) ++nBlocks.x;
    if (n_compute%nThreads.y) ++nBlocks.y;
    if (i_symmetrization == 1) {
        symmetrize<<<nBlocks,nThreads>>>(d_matrix_batch, n_compute);
    } else if (i_symmetrization == 2) {
        antisymmetrize<<<nBlocks,nThreads>>>(d_matrix_batch, n_compute);
    }
}

void evaluate_mr_batch_no_symm(const int &n_points,
                               const int &n_compute,
                               const double *wave,
                               const double *matrix_wave,
                               const int &i_symmetrization,
                               const int &i_dir,
                               const int &max_dims)
{
    using namespace gpuMR;
    if (i_dir == 1) {
        cublasSetVector(n_compute*n_points, sizeof(double),
                        wave, 1, d_wave, 1);
        cublasSetVector(max_dims*n_compute*n_points, sizeof(double),
                        matrix_wave, 1, d_matrix_wave, 1);
    }
    const double alpha = 1.0;
    const double beta = 0.0;
    cublasDgemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                n_compute, n_compute, n_points,
                &alpha, d_wave, n_compute,
                d_matrix_wave+(i_dir-1)*n_compute*n_points, n_compute,
                &beta, d_matrix_batch, n_compute);
}

__global__
void symmetrize(double *d_matrix_batch, const int n_compute)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb > nb && mb < n_compute)
        d_matrix_batch[mb+nb*n_compute] = (d_matrix_batch[mb+nb*n_compute] +
                                           d_matrix_batch[nb+mb*n_compute])/2;
}

__global__
void antisymmetrize(double *d_matrix_batch, const int n_compute)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute) {
        if (mb > nb) {
            d_matrix_batch[mb+nb*n_compute] =
                (d_matrix_batch[mb+nb*n_compute] -
                 d_matrix_batch[nb+mb*n_compute])/2;
        } else if (mb == nb) {
            d_matrix_batch[mb+nb*n_compute] = 0.0;
        }
    }
}

void update_mr_batch_gpu(const int &n_compute,
                         const int &starting_point,
                         const int &i_basis)
{
    using namespace gpuMR;
    cublasSetVector(n_compute, sizeof(int), &i_basis, 1, d_i_basis, 1);
    nThreads = dim3(maxThreadsX, maxThreadsY);
    nBlocks  = dim3((n_compute + nThreads.x - 1)/nThreads.x,
                    (n_compute + nThreads.y - 1)/nThreads.y);
    d_batch_to_packed<<<nBlocks,nThreads>>>(d_matrix_packed,
                                            d_matrix_batch,
                                            n_compute,
                                            starting_point,
                                            d_i_basis,
                                            d_basis_glb_to_loc);
}

__global__
void d_batch_to_packed(double *d_matrix_packed,
                       const double *d_matrix_batch,
                       const int n_compute,
                       const int starting_point,
                       const int *d_i_basis,
                       const int *d_basis_glb_to_loc)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute)  {
        int i_index = d_basis_glb_to_loc[d_i_basis[mb]-1];
        i_index = (i_index*(i_index-1))/2;
        if (nb <= mb) {
            i_index += d_basis_glb_to_loc[d_i_basis[nb]-1]-1;
            d_matrix_packed[i_index+starting_point] +=
                d_matrix_batch[mb+nb*n_compute];
        }
    }
}

void update_mr_batch_gpu_full(const int &n_compute,
                              const int &starting_point,
                              const int *i_basis,
                              const int &matrix_packed_ld)
{
    using namespace gpuMR;
    cublasSetVector(n_compute, sizeof(int), i_basis, 1, d_i_basis, 1);
    nThreads = dim3(maxThreadsX,maxThreadsY);
    nBlocks  = dim3((n_compute + nThreads.x - 1)/nThreads.x,
                    (n_compute + nThreads.y - 1)/nThreads.y);
    d_batch_to_packed_full<<<nBlocks,nThreads>>>(d_matrix_packed,
                                             d_matrix_batch,
                                             n_compute,
                                             starting_point,
                                             d_i_basis,
                                             d_basis_glb_to_loc,
                                             matrix_packed_ld);
}

__global__
void d_batch_to_packed_full(double *d_matrix_packed,
                            const double *d_matrix_batch,
                            const int n_compute,
                            const int starting_point,
                            const int *d_i_basis,
                            const int *d_basis_glb_to_loc,
                            const int matrix_packed_ld)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute)  {
        int i_index = d_basis_glb_to_loc[d_i_basis[mb]-1];
        i_index = (i_index*(i_index-1))/2;
        if (nb <= mb) {
            i_index += d_basis_glb_to_loc[d_i_basis[nb]-1]-1;
            d_matrix_packed[i_index+starting_point] +=
                d_matrix_batch[mb+nb*n_compute];
            d_matrix_packed[i_index+matrix_packed_ld+starting_point] +=
                d_matrix_batch[nb+mb*n_compute];
        }
    }
}

void symm_antisymm_gpu(const int &i_symmetrization, const int &n_compute)
{
    using namespace gpuMR;
    nThreads = dim3(maxThreadsX,maxThreadsY);
    nBlocks  = dim3(n_compute/nThreads.x, n_compute/nThreads.y);
    if (n_compute%nThreads.x) ++nBlocks.x;
    if (n_compute%nThreads.y) ++nBlocks.y;
    if (i_symmetrization == 1) {
        symmetrize<<<nBlocks,nThreads>>>(d_matrix_batch, n_compute);
    } else if (i_symmetrization == 2) {
        antisymmetrize<<<nBlocks,nThreads>>>(d_matrix_batch, n_compute);
    }
}

void giao_mult_psi_gpu(const bool &compact_indexing,
                       const bool &do_transpose,
                       const bool &mult_R_mn_both,
                       const int &i_dim,
                       const int &n_dims,
                       const int &n_compute,
                       const int &n_points,
                       const double *wave,
                       const double *matrix_wave,
                       const double *r_mn_basis,
                       double *matrix_batch_GIAO,
                       double *matrix_batch,
                       const int &n_spins)
{
    using namespace gpuMR;
    const int n_batches = 3;
    // Size of the current batch in terms of basis functions (not grid
    // points)
    const int batch_size = n_compute*n_compute;
    const double alpha = 1.0;
    const double beta = 0.0;
    int strideB;
    int strideC;
    if (i_dim == 1) {
        // Copy data to gpu
        cublasSetVector(n_compute*n_points, sizeof(double), wave, 1, d_wave, 1);
        cublasSetVector(matrix_wave_size_factor*n_spins*n_compute*n_points,
                        sizeof(double), matrix_wave, 1, d_matrix_wave, 1);
        cublasSetVector(3*batch_size, sizeof(double), r_mn_basis, 1,
                        d_r_mn_basis, 1);
        // Set up nThreads and nBlocks
        nThreads = dim3(maxThreadsX,maxThreadsY);
        nBlocks  = dim3(n_compute/nThreads.x,n_compute/nThreads.y);
        if (n_compute%nThreads.x) ++nBlocks.x;
        if (n_compute%nThreads.y) ++nBlocks.y;
    }
    // These integers denote shifts of indices in certains
    // expressions. For example, have a look at d_operate4: l_1 is the
    // shift of the matrix on the left hand side on the first line;
    // r_GIAO_11 and r_GIAO_12 are shifts of the first and second GIAO
    // work arrays on the right-hand side on the first line, and so
    // on.
    int l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12;
    int l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22;
    int l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32;
    int r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42;
    if (mult_R_mn_both) {
        if (compact_indexing) {
            if (i_dim == 1) {
                double *matrix_list[3][6];
                for (int i_dir = 0; i_dir < 6; ++i_dir) {
                    matrix_list[0][i_dir] = d_wave;
                    matrix_list[1][i_dir] =
                        d_matrix_wave + i_dir*n_compute*n_points;
                    matrix_list[2][i_dir] =
                        d_matrix_batch_GIAO + i_dir*batch_size;
                }
                for (int i = 0; i < 3; ++i)
                    cublasSetVector(6, sizeof(double*),
                                    matrix_list[i], 1, d_matrix_list[i], 1);
                cublasDgemmBatched
                    (cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                     n_compute, n_compute, n_points,
                     &alpha, (const double**)d_matrix_list[0], n_compute,
                     (const double**)d_matrix_list[1],
                     n_compute, &beta, d_matrix_list[2], n_compute, 6);
                if (do_transpose) {
                    d_do_transpose<<<nBlocks,nThreads>>>
                        (6, n_compute, d_matrix_batch_GIAO);
                }
            }
            if (n_dims == 9) {
                switch(i_dim) {
                case 1:
                    // chi_11 = R_2 A - R_3 B
                    l_1       = 6*batch_size;
                    r_GIAO_11 = 5*batch_size;
                    r_mn_11   = 2*batch_size;
                    r_GIAO_12 = 2*batch_size;
                    r_mn_12   = batch_size;
                    l_2       = 7*batch_size;
                    r_GIAO_21 = batch_size;
                    r_mn_21   = 2*batch_size;
                    r_GIAO_22 = 5*batch_size;
                    r_mn_22   = batch_size;
                    l_3       = 8*batch_size;
                    r_GIAO_31 = 3*batch_size;
                    r_mn_31   = 2*batch_size;
                    r_GIAO_32 = 4*batch_size;
                    r_mn_32   = batch_size;
                    r_mn_41   = batch_size;
                    r_GIAO_41 = 6*batch_size;
                    r_mn_42   = 2*batch_size;
                    r_GIAO_42 = 7*batch_size;
                    d_operate4<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 2:
                    // chi_21 = R_3 C - R_1 A
                    r_mn_11   = 2*batch_size;
                    r_GIAO_11 = 8*batch_size;
                    r_mn_12   = 0;
                    r_GIAO_12 = 6*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 3:
                    // chi_31 = R_1 B - R_2 C
                    r_mn_11   = 0;
                    r_GIAO_11 = 7*batch_size;
                    r_mn_12   = batch_size;
                    r_GIAO_12 = 8*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 4:
                    // chi_12 = R_2 A - R_3 B
                    l_1       = 6*batch_size;
                    r_GIAO_11 = 2*batch_size;
                    r_mn_11   = 0;
                    r_GIAO_12 = 4*batch_size;
                    r_mn_12   = 2*batch_size;
                    l_2       = 7*batch_size;
                    r_GIAO_21 = 5*batch_size;
                    r_mn_21   = 0;
                    r_GIAO_22 = 3*batch_size;
                    r_mn_22   = 2*batch_size;
                    l_3       = 8*batch_size;
                    r_GIAO_31 = 4*batch_size;
                    r_mn_31   = 0;
                    r_GIAO_32 = 0;
                    r_mn_32   = 2*batch_size;
                    r_mn_41   = batch_size;
                    r_GIAO_41 = 6*batch_size;
                    r_mn_42   = 2*batch_size;
                    r_GIAO_42 = 7*batch_size;
                    d_operate4<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 5:
                    // chi_22 = R_3 C - R_1 A
                    r_mn_11   = 2*batch_size;
                    r_GIAO_11 = 8*batch_size;
                    r_mn_12   = 0;
                    r_GIAO_12 = 6*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 6:
                    // chi_32 = R_1 B - R_2 C
                    r_mn_11   = 0;
                    r_GIAO_11 = 7*batch_size;
                    r_mn_12   = batch_size;
                    r_GIAO_12 = 8*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 7:
                    // chi_13 = R_2 A - R_3 B
                    l_1       = 6*batch_size;
                    r_GIAO_11 = 4*batch_size;
                    r_mn_11   = batch_size;
                    r_GIAO_12 = 5*batch_size;
                    r_mn_12   = 0;
                    l_2       = 7*batch_size;
                    r_GIAO_21 = 3*batch_size;
                    r_mn_21   = batch_size;
                    r_GIAO_22 = batch_size;
                    r_mn_22   = 0;
                    l_3       = 8*batch_size;
                    r_GIAO_31 = 0;
                    r_mn_31   = batch_size;
                    r_GIAO_32 = 3*batch_size;
                    r_mn_32   = 0;
                    r_mn_41   = batch_size;
                    r_GIAO_41 = 6*batch_size;
                    r_mn_42   = 2*batch_size;
                    r_GIAO_42 = 7*batch_size;
                    d_operate4<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 8:
                    // chi_23 = R_3 C - R_1 A
                    r_mn_11   = 2*batch_size;
                    r_GIAO_11 = 8*batch_size;
                    r_mn_12   = 0;
                    r_GIAO_12 = 6*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 9:
                    // chi_33 = R_1 B - R_2 C
                    r_mn_11   = 0;
                    r_GIAO_11 = 7*batch_size;
                    r_mn_12   = batch_size;
                    r_GIAO_12 = 8*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                }
            } else {
                switch(i_dim) {
                case 1:
                    // chi_11 = R_2 A - R_3 B
                    l_1       = 6*batch_size;
                    r_GIAO_11 = 5*batch_size;
                    r_mn_11   = 2*batch_size;
                    r_GIAO_12 = 2*batch_size;
                    r_mn_12   = batch_size;
                    l_2       = 7*batch_size;
                    r_GIAO_21 = batch_size;
                    r_mn_21   = 2*batch_size;
                    r_GIAO_22 = 5*batch_size;
                    r_mn_22   = batch_size;
                    r_mn_31   = batch_size;
                    r_GIAO_31 = 6*batch_size;
                    r_mn_32   = 2*batch_size;
                    r_GIAO_32 = 7*batch_size;
                    d_operate3<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 2:
                    // chi_33 = R_1 B - R_2 C
                    l_1       = 6*batch_size;
                    r_GIAO_11 = 2*batch_size;
                    r_mn_11   = 0;
                    r_GIAO_12 = 4*batch_size;
                    r_mn_12   = 2*batch_size;
                    l_2       = 7*batch_size;
                    r_GIAO_21 = 4*batch_size;
                    r_mn_21   = 0;
                    r_GIAO_22 = 0;
                    r_mn_22   = 2*batch_size;
                    r_mn_31   = 2*batch_size;
                    r_GIAO_31 = 7*batch_size;
                    r_mn_32   = 0;
                    r_GIAO_32 = 6*batch_size;
                    d_operate3<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 3:
                    l_1       = 6*batch_size;
                    r_GIAO_11 = 3*batch_size;
                    r_mn_11   = batch_size;
                    r_GIAO_12 = batch_size;
                    r_mn_12   = 0;
                    l_2       = 7*batch_size;
                    r_GIAO_21 = 0;
                    r_mn_21   = batch_size;
                    r_GIAO_22 = 3*batch_size;
                    r_mn_22   = 0;
                    r_mn_31   = 0;
                    r_GIAO_31 = 6*batch_size;
                    r_mn_32   = batch_size;
                    r_GIAO_32 = 7*batch_size;
                    d_operate3<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                }
            }
        } else { // not compact indexing
            if (i_dim == 1) {
                double *matrix_list[3][9];
                for (int i_dir = 0; i_dir < 9; ++i_dir) {
                    matrix_list[0][i_dir] = d_wave;
                    matrix_list[1][i_dir] =
                        d_matrix_wave + i_dir*n_compute*n_points;
                    matrix_list[2][i_dir] =
                        d_matrix_batch_GIAO + i_dir*batch_size;
                }
                for (int i = 0; i < 3; ++i)
                    cublasSetVector(9, sizeof(double*),
                                    matrix_list[i], 1, d_matrix_list[i], 1);
                cublasDgemmBatched
                    (cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                     n_compute, n_compute, n_points,
                     &alpha, (const double**)d_matrix_list[0], n_compute,
                     (const double**)d_matrix_list[1],
                     n_compute, &beta, d_matrix_list[2], n_compute, 9);
                if (do_transpose) {
                    d_do_transpose<<<nBlocks,nThreads>>>
                        (9, n_compute, d_matrix_batch_GIAO);
                }
            }
            if (n_dims == 9) {
                switch(i_dim) {
                case 1:
                    // chi_11 = R_2 A - R_3 B
                    l_1       = 9*batch_size;
                    r_GIAO_11 = 5*batch_size;
                    r_mn_11   = 2*batch_size;
                    r_GIAO_12 = 8*batch_size;
                    r_mn_12   = batch_size;
                    l_2       = 10*batch_size;
                    r_GIAO_21 = 4*batch_size;
                    r_mn_21   = 2*batch_size;
                    r_GIAO_22 = 7*batch_size;
                    r_mn_22   = batch_size;
                    l_3       = 11*batch_size;
                    r_GIAO_31 = 3*batch_size;
                    r_mn_31   = 2*batch_size;
                    r_GIAO_32 = 6*batch_size;
                    r_mn_32   = batch_size;
                    r_mn_41   = batch_size;
                    r_GIAO_41 = 9*batch_size;
                    r_mn_42   = 2*batch_size;
                    r_GIAO_42 = 10*batch_size;
                    d_operate4<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 2:
                    // chi_21 = R_3 C - R_1 A
                    r_mn_11   = 2*batch_size;
                    r_GIAO_11 = 11*batch_size;
                    r_mn_12   = 0;
                    r_GIAO_12 = 9*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 3:
                    // chi_31 = R_1 B - R_2 C
                    r_mn_11   = 0;
                    r_GIAO_11 = 10*batch_size;
                    r_mn_12   = batch_size;
                    r_GIAO_12 = 11*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 4:
                    // chi_12 = R_2 A - R_3 B
                    l_1       = 9*batch_size;
                    r_GIAO_11 = 8*batch_size;
                    r_mn_11   = 0;
                    r_GIAO_12 = 2*batch_size;
                    r_mn_12   = 2*batch_size;
                    l_2       = 10*batch_size;
                    r_GIAO_21 = 7*batch_size;
                    r_mn_21   = 0;
                    r_GIAO_22 = batch_size;
                    r_mn_22   = 2*batch_size;
                    l_3       = 11*batch_size;
                    r_GIAO_31 = 6*batch_size;
                    r_mn_31   = 0;
                    r_GIAO_32 = 0;
                    r_mn_32   = 2*batch_size;
                    r_mn_41   = batch_size;
                    r_GIAO_41 = 9*batch_size;
                    r_mn_42   = 2*batch_size;
                    r_GIAO_42 = 10*batch_size;
                    d_operate4<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 5:
                    // chi_22 = R_3 C - R_1 A
                    r_mn_11   = 2*batch_size;
                    r_GIAO_11 = 11*batch_size;
                    r_mn_12   = 0;
                    r_GIAO_12 = 9*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 6:
                    // chi_32 = R_1 B - R_2 C
                    r_mn_11   = 0;
                    r_GIAO_11 = 10*batch_size;
                    r_mn_12   = batch_size;
                    r_GIAO_12 = 11*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 7:
                    // chi_13 = R_2 A - R_3 B
                    l_1       = 9*batch_size;
                    r_GIAO_11 = 2*batch_size;
                    r_mn_11   = batch_size;
                    r_GIAO_12 = 5*batch_size;
                    r_mn_12   = 0;
                    l_2       = 10*batch_size;
                    r_GIAO_21 = batch_size;
                    r_mn_21   = batch_size;
                    r_GIAO_22 = 4*batch_size;
                    r_mn_22   = 0;
                    l_3       = 11*batch_size;
                    r_GIAO_31 = 0;
                    r_mn_31   = batch_size;
                    r_GIAO_32 = 3*batch_size;
                    r_mn_32   = 0;
                    r_mn_41   = batch_size;
                    r_GIAO_41 = 9*batch_size;
                    r_mn_42   = 2*batch_size;
                    r_GIAO_42 = 10*batch_size;
                    d_operate4<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         l_3, r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         r_GIAO_41, r_mn_41, r_GIAO_42, r_mn_42,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 8:
                    // chi_23 = R_3 C - R_1 A
                    r_mn_11   = 2*batch_size;
                    r_GIAO_11 = 11*batch_size;
                    r_mn_12   = 0;
                    r_GIAO_12 = 9*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 9:
                    // chi_33 = R_1 B - R_2 C
                    r_mn_11   = 0;
                    r_GIAO_11 = 10*batch_size;
                    r_mn_12   = batch_size;
                    r_GIAO_12 = 11*batch_size;
                    d_operate1<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                }
            } else {
                switch(i_dim) {
                case 1:
                    // chi_11 = R_2 A - R_3 B
                    l_1       = 9*batch_size;
                    r_GIAO_11 = 5*batch_size;
                    r_mn_11   = 2*batch_size;
                    r_GIAO_12 = 8*batch_size;
                    r_mn_12   = batch_size;
                    l_2       = 10*batch_size;
                    r_GIAO_21 = 4*batch_size;
                    r_mn_21   = 2*batch_size;
                    r_GIAO_22 = 7*batch_size;
                    r_mn_22   = batch_size;
                    r_mn_31   = batch_size;
                    r_GIAO_31 = 9*batch_size;
                    r_mn_32   = 2*batch_size;
                    r_GIAO_32 = 10*batch_size;
                    d_operate3<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 2:
                    // chi_22 = R_3 C - R_1 A
                    l_1       = 9*batch_size;
                    r_GIAO_11 = 8*batch_size;
                    r_mn_11   = 0;
                    r_GIAO_12 = 2*batch_size;
                    r_mn_12   = 2*batch_size;
                    l_2       = 10*batch_size;
                    r_GIAO_21 = 6*batch_size;
                    r_mn_21   = 0;
                    r_GIAO_22 = 0;
                    r_mn_22   = 2*batch_size;
                    r_mn_31   = 2*batch_size;
                    r_GIAO_31 = 10*batch_size;
                    r_mn_32   = 0;
                    r_GIAO_32 = 9*batch_size;
                    d_operate3<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                case 3:
                    // chi_33 = R_1 B - R_2 C
                    l_1       = 9*batch_size;
                    r_GIAO_11 = batch_size;
                    r_mn_11   = batch_size;
                    r_GIAO_12 = 4*batch_size;
                    r_mn_12   = 0;
                    l_2       = 10*batch_size;
                    r_GIAO_21 = 0;
                    r_mn_21   = batch_size;
                    r_GIAO_22 = 3*batch_size;
                    r_mn_22   = 0;
                    r_mn_31   = 0;
                    r_GIAO_31 = 9*batch_size;
                    r_mn_32   = batch_size;
                    r_GIAO_32 = 10*batch_size;
                    d_operate3<<<nBlocks,nThreads>>>
                        (n_compute, d_r_mn_basis,
                         l_1, r_GIAO_11, r_mn_11, r_GIAO_12, r_mn_12,
                         l_2, r_GIAO_21, r_mn_21, r_GIAO_22, r_mn_22,
                         r_GIAO_31, r_mn_31, r_GIAO_32, r_mn_32,
                         d_matrix_batch, d_matrix_batch_GIAO);
                    break;
                }
            }
        }
    } else {
        // Copy data to gpu
        if (compact_indexing) {
            double *matrix_list[3][2];
            switch(i_dim) {
            case 1:
                // Do the dgemm
                matrix_list[0][0] = d_wave;
                matrix_list[0][1] = d_wave;
                matrix_list[1][0] = d_matrix_wave;
                matrix_list[1][1] = d_matrix_wave + n_compute*n_points;
                matrix_list[2][0] = d_matrix_batch;
                matrix_list[2][1] = d_matrix_batch_GIAO;
                for (int i = 0; i < 3; ++i)
                    cublasSetVector(2, sizeof(double*),
                                    matrix_list[i], 1, d_matrix_list[i], 1);
                cublasDgemmBatched
                    (cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                     n_compute, n_compute, n_points,
                     &alpha, (const double**)d_matrix_list[0], n_compute,
                     (const double**)d_matrix_list[1],
                     n_compute, &beta, d_matrix_list[2], n_compute, 2);
                // Cross products
                d_cross_product_modified<<<nBlocks,nThreads>>>
                    (n_compute, d_r_mn_basis, d_matrix_batch_GIAO,
                     batch_size, 2*batch_size, false, d_matrix_batch);
                break;
            case 2:
                // Do the dgemm
                matrix_list[1][0] = d_matrix_wave+2*n_compute*n_points;
                matrix_list[1][1] = d_matrix_wave+3*n_compute*n_points;
                cublasSetVector(2, sizeof(double*),
                                matrix_list[1], 1, d_matrix_list[1], 1);
                cublasDgemmBatched
                    (cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                     n_compute, n_compute, n_points,
                     &alpha, (const double**)d_matrix_list[0], n_compute,
                     (const double**)d_matrix_list[1],
                     n_compute, &beta, d_matrix_list[2], n_compute, 2);
                // Cross products
                d_cross_product_modified<<<nBlocks,nThreads>>>
                    (n_compute, d_r_mn_basis, d_matrix_batch_GIAO,
                     2*batch_size, 0, true, d_matrix_batch);
                break;
            case 3:
                // Do the dgemm
                matrix_list[1][0] = d_matrix_wave+4*n_compute*n_points;
                matrix_list[1][1] = d_matrix_wave+5*n_compute*n_points;
                cublasSetVector(2, sizeof(double*),
                                matrix_list[1], 1, d_matrix_list[1], 1);
                cublasDgemmBatched
                    (cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                     n_compute, n_compute, n_points,
                     &alpha, (const double**)d_matrix_list[0], n_compute,
                     (const double**)d_matrix_list[1],
                     n_compute, &beta, d_matrix_list[2], n_compute, 2);
                // Cross products
                d_cross_product_modified<<<nBlocks,nThreads>>>
                    (n_compute, d_r_mn_basis, d_matrix_batch_GIAO,
                     0, batch_size, false, d_matrix_batch);
                break;
            }
            cublasGetVector(batch_size, sizeof(double),
                            d_matrix_batch, 1, matrix_batch, 1);
            cublasGetVector(batch_size, sizeof(double),
                            d_matrix_batch_GIAO, 1, matrix_batch_GIAO, 1);
        } else {
            switch(i_dim) {
            case 1:
            case 4:
            case 7:
                // Do the dgemm
                strideB = n_compute*n_points;
                strideC = batch_size;
                cublasDgemmStridedBatched
                    (cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                     n_compute, n_compute, n_points,
                     &alpha, d_wave, n_compute, 0,
                     d_matrix_wave+(i_dim-1)*n_compute*n_points,
                     n_compute, strideB, &beta, d_matrix_batch_GIAO, n_compute,
                     strideC, n_batches);
                // do_tranpose if necessary
                if (do_transpose) {
                    d_do_transpose<<<nBlocks,nThreads>>>
                        (3, n_compute, d_matrix_batch_GIAO);
                }
                // Cross products
                d_cross_product<<<nBlocks,nThreads>>>
                    (n_compute, d_r_mn_basis, d_matrix_batch_GIAO,
                     batch_size, 2*batch_size, d_matrix_batch);
                break;
            case 2:
            case 5:
            case 8:
                d_cross_product<<<nBlocks,nThreads>>>
                    (n_compute, d_r_mn_basis, d_matrix_batch_GIAO,
                     2*batch_size, 0, d_matrix_batch);
                break;
            case 3:
            case 6:
            case 9:
                d_cross_product<<<nBlocks,nThreads>>>
                    (n_compute, d_r_mn_basis, d_matrix_batch_GIAO,
                     0, batch_size, d_matrix_batch);
                break;
            }
        }
    }
}

__global__
void d_operate1(const int n_compute,
                const double *d_r_mn_basis,
                const int r_GIAO_1,
                const int r_mn_1,
                const int r_GIAO_2,
                const int r_mn_2,
                double *matrix_batch,
                double *matrix_batch_GIAO)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute && nb < n_compute) {
        const int mn = mb + nb*n_compute;
        matrix_batch[mn] =
            d_r_mn_basis[mn+r_mn_1]*matrix_batch_GIAO[mn+r_GIAO_1] -
            d_r_mn_basis[mn+r_mn_2]*matrix_batch_GIAO[mn+r_GIAO_2];
    }

}

__global__
void d_operate3(const int n_compute,
                const double *d_r_mn_basis,
                const int l_1,
                const int r_GIAO_11,
                const int r_mn_11,
                const int r_GIAO_12,
                const int r_mn_12,
                const int l_2,
                const int r_GIAO_21,
                const int r_mn_21,
                const int r_GIAO_22,
                const int r_mn_22,
                const int r_GIAO_31,
                const int r_mn_31,
                const int r_GIAO_32,
                const int r_mn_32,
                double *matrix_batch,
                double *matrix_batch_GIAO)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute && nb < n_compute) {
        const int mn = mb + nb*n_compute;
        matrix_batch_GIAO[mn+l_1] =
            matrix_batch_GIAO[mn+r_GIAO_11]*d_r_mn_basis[mn+r_mn_11] -
            matrix_batch_GIAO[mn+r_GIAO_12]*d_r_mn_basis[mn+r_mn_12];
        matrix_batch_GIAO[mn+l_2] =
            matrix_batch_GIAO[mn+r_GIAO_21]*d_r_mn_basis[mn+r_mn_21] -
            matrix_batch_GIAO[mn+r_GIAO_22]*d_r_mn_basis[mn+r_mn_22];
        matrix_batch[mn] =
            d_r_mn_basis[mn+r_mn_31]*matrix_batch_GIAO[mn+r_GIAO_31] -
            d_r_mn_basis[mn+r_mn_32]*matrix_batch_GIAO[mn+r_GIAO_32];
    }

}

__global__
void d_operate4(const int n_compute,
                const double *d_r_mn_basis,
                const int l_1,
                const int r_GIAO_11,
                const int r_mn_11,
                const int r_GIAO_12,
                const int r_mn_12,
                const int l_2,
                const int r_GIAO_21,
                const int r_mn_21,
                const int r_GIAO_22,
                const int r_mn_22,
                const int l_3,
                const int r_GIAO_31,
                const int r_mn_31,
                const int r_GIAO_32,
                const int r_mn_32,
                const int r_GIAO_41,
                const int r_mn_41,
                const int r_GIAO_42,
                const int r_mn_42,
                double *matrix_batch,
                double *matrix_batch_GIAO)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute && nb < n_compute) {
        const int mn = mb + nb*n_compute;
        matrix_batch_GIAO[mn+l_1] =
            matrix_batch_GIAO[mn+r_GIAO_11]*d_r_mn_basis[mn+r_mn_11] -
            matrix_batch_GIAO[mn+r_GIAO_12]*d_r_mn_basis[mn+r_mn_12];
        matrix_batch_GIAO[mn+l_2] =
            matrix_batch_GIAO[mn+r_GIAO_21]*d_r_mn_basis[mn+r_mn_21] -
            matrix_batch_GIAO[mn+r_GIAO_22]*d_r_mn_basis[mn+r_mn_22];
        matrix_batch_GIAO[mn+l_3] =
            matrix_batch_GIAO[mn+r_GIAO_31]*d_r_mn_basis[mn+r_mn_31] -
            matrix_batch_GIAO[mn+r_GIAO_32]*d_r_mn_basis[mn+r_mn_32];
        matrix_batch[mn] =
            d_r_mn_basis[mn+r_mn_41]*matrix_batch_GIAO[mn+r_GIAO_41] -
            d_r_mn_basis[mn+r_mn_42]*matrix_batch_GIAO[mn+r_GIAO_42];
    }

}

__global__
void d_do_transpose(const int n_batches, const int n_compute, double *matrix)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb > nb && mb < n_compute) {
        for (int i_batch = 0; i_batch < n_batches; ++i_batch) {
            const int mn = mb + nb*n_compute + i_batch*n_compute*n_compute;
            const int nm = nb + mb*n_compute + i_batch*n_compute*n_compute;
            matrix[mn] += matrix[nm];
            matrix[nm] = matrix[mn];
        }
    }
}

__global__
void d_cross_product(const int n_compute,
                     const double *d_r_mn_basis,
                     const double *d_matrix_batch_aux,
                     const int shift1,
                     const int shift2,
                     double *d_matrix_batch)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute && nb < n_compute) {
        const int mn = mb + nb*n_compute;
        d_matrix_batch[mn] =
            d_r_mn_basis[mn+shift1]*d_matrix_batch_aux[mn+shift2] -
            d_r_mn_basis[mn+shift2]*d_matrix_batch_aux[mn+shift1];
    }
}

__global__
void d_cross_product_modified(const int n_compute,
                              const double *d_r_mn_basis,
                              const double *d_matrix_batch_aux,
                              const int shift1,
                              const int shift2,
                              const bool negative,
                              double *d_matrix_batch)
{
    const int mb = threadIdx.x + blockIdx.x*blockDim.x;
    const int nb = threadIdx.y + blockIdx.y*blockDim.y;
    if (mb < n_compute && nb < n_compute) {
        const int mn = mb + nb*n_compute;
        if (negative) {
            d_matrix_batch[mn] =
                d_r_mn_basis[mn+shift1]*d_matrix_batch[mn] -
                d_r_mn_basis[mn+shift2]*d_matrix_batch_aux[mn];
        } else {
            d_matrix_batch[mn] =
                d_r_mn_basis[mn+shift1]*d_matrix_batch_aux[mn] -
                d_r_mn_basis[mn+shift2]*d_matrix_batch[mn];
        }
    }
}

void get_matrix_packed_gpu(double *matrix_packed, const int &matrix_packed_size)
{
    using namespace gpuMR;
    cublasGetVector(matrix_packed_size, sizeof(double), d_matrix_packed, 1,
                    matrix_packed, 1);
}

void mr_destroy_gpu()
{
    using namespace gpuMR;
    cudaFree(d_matrix_batch);
    cudaFree(d_matrix_batch_GIAO);
    cudaFree(d_r_mn_basis);
    cudaFree(d_wave);
    cudaFree(d_matrix_wave);
    cudaFree(d_matrix_packed);
    cudaFree(d_basis_glb_to_loc);
    cudaFree(d_i_basis);
    for (int i = 0; i < 3; ++i)
        cudaFree(d_matrix_list[i]);
}
