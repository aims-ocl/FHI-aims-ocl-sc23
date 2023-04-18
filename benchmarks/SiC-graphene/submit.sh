#!/bin/bash                                                                                                                          
#PBS -N benchmark
#PBS -l select=10:ncpus=24:mpiprocs=24                                                                                                               
#PBS -l walltime=4:00:00                                                                                                            
#PBS -q normal

#PBS -m ea                                                                                                                           
#PBS -M <email>@<organization>
#PBS -j oe                                                                                                                           
#PBS -P <program>

# Example submit script for the PBS batch queueing system as installed
# at the South African Centre for High-Performance Computing (CHPC). This
# submit script does give the principle of settings needed to run FHI-aims
# on typical HPC machines but will have to be adjusted in ways specific 
# to any other supercomputer / queueing system. The #PBS flags above 
# specify a run on 10 nodes with 24 cores and 24 MPI tasks per node 
# (i.e., a total of 240 cores and MPI tasks).

# Making sure that the right compiler/runtime, MPI and linear algebra 
# libraries are in the $LD_LIBRARY_PATH so that the code finds them
# is essential. The following line is just an example that works to
# set the environment variables for the Intel 2016 suite on the
# high-performance computer Lengau at the South African Centre for CHPC. 
# The line must be altered appropriately for any other supercomputer.
module load chpc/parallel_studio_xe/16.0.1/2016.1.150

# PBS-specific way of determining the number of processors used
# in this PBS submit script
nproc=`cat $PBS_NODEFILE | wc -l`

# FHI-aims version variable to call the right binary (see below).
# Adjust <version-stamp> for your own setting.
aims_version="<version-stamp>.scalapack.mpi"

# Adjust working directory for your own needs (e.g., some temporary
# directory prescribed by your computing facility)
workdir=$HOME/work/Ac-Lys-Ala19-H.$aims_version.$nproc

# VERY important setting for ANY FHI-aims run. Please ensure that this
# setting is correctly propagated through to ALL nodes used by the mpirun
# command or its equivalent. If this variable is not altered from its
# standard settings, FHI-aims may be stopped by the operating system
# (segmentation fault) seemingly at random. A larger stack size limit than 
# the default (ideally, unlimited) will avoid this behavior. If FHI-aims 
# is compiled with a CC variable (C compiler) specified in make.sys, then 
# the code will write out the stack size found on each processor / MPI 
# task near the beginning of its standard output.
ulimit -s unlimited

# Very important Intel MKL-specific settings. Since FHI-aims does not
# use shared memory parallelization (OpenMP) internally - MPI produces
# similar of better scaling in our tests on different platforms - 
# we must make sure that Intel MPI does not incorrectly adds another 
# parallelization layer. (For instance, when MPI-parallelizing 
# over 24 cores, Intel MKL might try to use 24 threads per MPI task by
# default. This would lead to flooding the computer with 576 threads.
# The following variables prevent this behavior.)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

# Directory from which PBS' qsub command was called.
cd $PBS_O_WORKDIR

mkdir -p $workdir
cp geometry.in $workdir
cp control.in $workdir

cd $workdir

  # Check that we are indeed using an mpirun or equivalent command
  # from the right location - in the case of this example, Intel's mpirun.
  which mpirun

  # Execution.
  # Please substitute the correct path to your FHI-aims binary here.
  mpirun -np $nproc <location of FHI-aims binary>/aims.$aims_version.x < /dev/null > aims.out.$aims_version 2>&1

cd $PBS_O_WORKDIR
mv $workdir .
