#!/bin/bash -l
#$ -j y
#$ -cwd
#$ -m n
#$ -N test
#$ -pe impi 96
#$ -l h_rt=86400.0
#$ -l h_vmem=45G
module purge
module load intel mkl impi
module list
export TMPDIR=/tmp

filename="potential_esp_32.cube"
filename2="ACF.dat"
ind=$(seq 1.1 0.1 3 )

esp_max=4.0
R_c=20.0
km=7
rm=7
beta=1.0

for i in $ind; do 

esp_min=$i




output="esp_charges_bader_"$esp_min"_"$esp_max".out"  


mpiexec -ppn 12 -envlist LD_LIBRARY_PATH -n $NSLOTS /mnt/lxfs2/scratch/bieniek/esp_charges/Fit_esp_charges $filename $esp_min $esp_max $R_c $rm $km> $output;
done
#mpiexec -ppn 12 -envlist LD_LIBRARY_PATH -n $NSLOTS /mnt/lxfs2/scratch/bieniek/esp_charges/Fit_esp_charges $filename $esp_min $esp_max $filename2 $beta> $output

