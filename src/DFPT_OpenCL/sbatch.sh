#!/bin/bash
#SBATCH --partition=dongsheng
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --gres=dcu:1

##SBATCH --ntasks-per-node=16
##SBATCH --ntasks=4
##SBATCH --cpus-per-task=1
##SBATCH --cpus-per-gpu=8
##SBATCH --ntasks=32
##SBATCH --mem-per-cpu=3200M

# source /public/software/compiler/intel/oneapi/mkl/2021.3.0/env/vars.sh
# source /public/software/compiler/intel/oneapi/mpi/2021.3.0/env/vars.sh
# mpirun -np 32 /public/home/autopar/FHI-aims-test/fhi-aims_MPE_O3_local_index_final/bin/aims3.191127.scalapack.mpi.x

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

echo $PATH
echo $LD_LIBRARY_PATH

# echo "SLURM_PROCID=$SLURM_PROCID, SLURM_ARRAY_TASK_MAX=$SLURM_ARRAY_TASK_MAX, SLURM_ARRAY_TASK_MIN=$SLURM_ARRAY_TASK_MIN"

# echo "SLURM_PROCID=$SLURM_PROCID, SLURM_NODEID=$SLURM_NODEID, SLURM_NTASKS=$SLURM_NTASKS"

# export GPU_DEVICE_ORDINAL=0,1

# /opt/hpc/software/mpi/hpcx/v2.7.4/gcc-7.3.1/bin/mpirun -np 16 ./a.out
# srun --mpi=pmi2 ./a.out
mpirun -np 1 ./a.out
# srun ./a.out
# ./aims.191127.scalapack.mpi.x
# srun ../aims.191127.scalapack.mpi.x


