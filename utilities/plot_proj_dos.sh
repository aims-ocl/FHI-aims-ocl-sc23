#!/bin/bash

if [[ "$1" == "" ]]
then
	echo "A simple script for plotting projected densities of state (either atom or species) using xmgrace."
	echo "To run, specify one of the projected DOS output files (atom_projected_dos_*_.dat or *_l_proj_dos.dat) as a command line parameter."
	echo "This script has been hacked together to serve my (WPH's ) needs, but is provided to the community all the same."
	echo "Caveat computor."
	echo "If a KS_DOS_total.dat file is in the working directory, will plot it alongside the projected DOS, but is not required"
	echo "(note this assumes that you're plotting the non-raw projected DOS.)"
	exit
fi

echo "This script has been hacked together to serve my (WPH's ) needs, but is provided to the community all the same."
echo "Caveat computor."
echo 

NCOLUMNS=$(sed -n '4p' $1 | wc | awk '{print $2}')

if [[ -a total_proj_dos ]]
then
	echo "File \"total_proj_dos\" already exists, will not overwrite, exiting."
	exit
fi

echo "Creating intermediate files."
for i in $(seq 2 1 $NCOLUMNS)
do
	if [[ -a l=$((i-3)) ]] 
	then
		echo "File \"l=$((i-3))\" already exists, will not overwrite, exiting."
		exit
	fi 
done


for i in $(seq 2 1 $NCOLUMNS)
do
	if [[ $i -eq 2 ]] 
	then
		tail -n +4 $1  | awk "{print \$1, \$$i}" > total_proj_dos 
	else 
		tail -n +4 $1  | awk "{print \$1, \$$i}" > l=$((i-3))		
	fi
done

if [[ -a KS_DOS_total.dat ]]
then
	echo "Plotting the *shifted* total DOS alongside the projected DOS."
	xmgrace -legend load -nxy KS_DOS_total.dat total_proj_dos l=*
else
	echo "Could not find a \"KS_DOS_total.dat\" file, will only plot the projected DOS."
	xmgrace -legend load -nxy species_total_dos l=*
fi

echo "Deleting intermediate files."
rm -f total_proj_dos
rm -f l=*
