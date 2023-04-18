#!/bin/bash -eu

module purge
module load intel-compilers-14.0
module load intel-mkl-11.1.1
module load mpich-3.1.4
module load cmake-3.1.3
module load python3.3

export OMP_NUM_THREADS=1
export MPICH_FC=ifort
ulimit -s unlimited

echo "==================="
echo "= Default Modules ="
echo "==================="
module list
echo

echo "===================="
echo "= GNU Make Version ="
echo "===================="
make --version
echo

echo "========================="
echo "= Environment Variables ="
echo "========================="
echo "MPIEXE = $MPIEXE"
echo "N_TASKS_BUILD = $N_TASKS_BUILD"
echo "N_TASKS_TEST = $N_TASKS_TEST"
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "MPICH_FC = $MPICH_FC"
echo "ulimit -s = $(ulimit -s)"
