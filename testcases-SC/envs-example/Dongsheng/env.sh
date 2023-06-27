module unload mpi/hpcx/2.11.0/gcc-7.3.1
module unload compiler/rocm/dtk/22.10.1
module load mpi/intelmpi/2017.4.239
module load compiler/rocm/dtk/21.04
module load compiler/cmake/3.23.1
source /public/software/compiler/intel/oneapi/mkl/2021.3.0/env/vars.sh

export MY_ROCM_PATH="/public/software/compiler/dtk/21.04"

export PATH="$MY_ROCM_PATH/opencl/bin:$PATH"
export LD_LIBRARY_PATH="$MY_ROCM_PATH/lib64:$MY_ROCM_PATH/lib:$MY_ROCM_PATH/llvm/lib:$MY_ROCM_PATH/opencl/lib:$LD_LIBRARY_PATH"
