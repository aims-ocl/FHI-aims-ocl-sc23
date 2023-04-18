#######################################################
#example SGE batch script for a plain MPI program using
#512 cores (16 nodes x 32 cores per node) with 32 MPI 
#tasks per node
#executed under the intel runtime 
#######################################################

## run in /bin/bash
#$ -S /bin/bash
## do not join stdout and stderr
#$ -j n
## name of the job
#$ -N esp
## execute job from the current working directory
#$ -cwd
## do not send mail
#$ -m n
## request 16 nodes (x 32 cores), must be a multiple of 32  
#$ -pe impi_hydra 192
## run for xx minutes
#$ -l h_rt=24:00:00

module load intel impi mkl
ulimit -s unlimited
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_HOME/lib/intel64/ 

## to save memory at high core counts, the "connectionless user datagram protocol" 
## can be enabled (might come at the expense of speed)
## see https://software.intel.com/en-us/articles/dapl-ud-support-in-intel-mpi-library
# export I_MPI_DAPL_UD=1

##gather MPI statistics to be analyzed with itac's mps tool
# export I_MPI_STATS=all

##gather MPI debug information (high verbosity)
# export I_MPI_DEBUG=5

filename="potential_esp_1.cube"
filename2="ACF.dat"

ind=$(seq 1.1 0.1 3 )

esp_max=4.0
R_c=20.0
km=7
rm=7
beta=1.0

for i in $ind; do 

esp_min=$i




output="esp_charges_"$esp_min"_"$esp_max".out" 

mpiexec -n $NSLOTS /scratch/biebj/esp_charges/Fit_esp_charges $filename $esp_min $esp_max $R_c $rm $km> $output

#mpiexec -n $NSLOTS /scratch/biebj/esp_charges/Fit_esp_charges $filename $esp_min $esp_max $filename2 $beta> $output

done