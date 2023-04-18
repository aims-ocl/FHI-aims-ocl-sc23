#!/bin/bash
#SBATCH --partition=dongsheng
#SBATCH --nodes=8
#SBATCH --gres=dcu:4
#SBATCH --tasks-per-node=32
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1

#//SBATCH --time=01:00:00

export OMP_NUM_THREADS=1

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo $SLURM_NODELIST
echo $SLURM_CLUSTER_NAME
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

filename="../../../bin/aims.191127.scalapack.mpi.x"

echo "filename=$filename"

ldd ${filename}

mpirun ${filename}
