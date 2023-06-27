#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --gres=dcu:1

export FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE=512
export FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM=OFF
export FHI_AIMS_OCL=ON
export FHI_AIMS_OCL_DFPT_POLAR_RHO=ON
export FHI_AIMS_OCL_DFPT_DIEL_SUMUP=ON
export FHI_AIMS_OCL_DFPT_POLAR_H=ON
export MPI_PER_NODE=1
export GPU_PER_NODE=1
export FHI_AIMS_OCL_DEBUG_IO=ON


echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo $SLURM_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

module unload mpi/hpcx/2.11.0/gcc-7.3.1
module unload compiler/rocm/dtk/22.10.1
module load compiler/rocm/dtk/21.04
module load mpi/intelmpi/2017.4.239

source /public/software/compiler/dtk/21.04/env.sh
source /public/software/compiler/intel/oneapi/mkl/2021.3.0/env/vars.sh
export MY_ROCM_PATH="/public/software/compiler/dtk/21.04"

export PATH="$MY_ROCM_PATH/opencl/bin:$PATH"
export LD_LIBRARY_PATH="$MY_ROCM_PATH/lib64:$MY_ROCM_PATH/lib:$MY_ROCM_PATH/llvm/lib:$MY_ROCM_PATH/opencl/lib:$LD_LIBRARY_PATH"

export OMP_NUM_THREADS=1

filename="../../build/aims.191127.scalapack.mpi.x"

echo ${filename}
ldd ${filename}
whereis mpirun

mpirun ${filename}
