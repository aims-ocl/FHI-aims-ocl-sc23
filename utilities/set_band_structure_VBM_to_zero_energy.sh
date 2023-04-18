#!/bin/bash

# This script shifts the band structure output by aims (band****.out) so that the highest occupied state coincides with
# the zero of energy.  This is necessary for comparison between different calculations, as the default zero of energy for
# band structures in aims is taken to be the Fermi level, which will differ from calculation to calculation as parameters
# are tweaked.  Even for metals, this will cause noticeable misalignments.
#
# This script runs in the folder containing the band structure, with no parameters needed.  If you provide a number as a
# parameter, it will shift the band structure by that number instead of zeroing out the highest occupied state.
#
# Note this script will not shift the DOS, because I (WPH) have never needed that functionality.  This should be trivial 
# to add in, if you are so inclined.
#
# While this script is provided without warranty, as always in the utility aims folders, this is a slight modification of
# a script that I have used thousands (yes, thousands) of times for a number of different types of band structures.  
# Spin-polarized, non-spin-polarized, spin-orbit-coupled, non-spin-orbit-coupled, insulating, metallic. I'm pretty sure it 
# works.

if [[ ! -a band1001.out ]] && [[ ! -a band2001.out ]]
then
  echo "Could not find either band1001.out or band2001.out, presuming there are no band structure files, exiting."
  exit -1
fi

if [[ -a band1001.out.original ]]
then
  echo "band1001.out.original already exists in this folder, suggesting that you've already shifted this band structure"
  echo "before.  Cowardly refusing to overwrite the backup band structure."
  exit -1
fi

if [ $# -ne 0 ]
then
  FERMI=$1
else
  FERMI="-999999999999"
  if [[ ! -a band1001.out ]] && [[ ! -a band2001.out ]]
  then
    echo "Could not find either band1001.out or band2001.out, presuming there are no band structure files, exiting."
    exit -1
  fi

  for i in band*.out 
  do
    TEMP_FERMI=$(python -c "
import sys

# This script finds the highest occupied state for a given segment of the band structure.

theFile = open(sys.argv[1],'r').readlines()
currentFermi=-9999999999
fullOcc = float(theFile[0].split()[4]) # The 1s state of the first k-point will always be fully occupied, use it to establish full occupancy.

for line in theFile:
  line = line.split()
  for iter,val in enumerate(line):
    if iter > 3 and (iter - 4) % 2 == 0:
      if float(val) >= fullOcc/2. and float(line[iter+1]) > currentFermi:
        currentFermi = float(line[iter+1])
print \"%20.10f\" % (currentFermi)" $i )

    # Over all segments of the band structure, find the highest value for the highest occupied state
    if [ $(echo "($TEMP_FERMI)>($FERMI)" | bc -l) -eq 1 ]; then
      FERMI=$TEMP_FERMI
    fi
  done
fi

for i in band*.out 
do
  python -c '
import sys

# This script takes in a band*.out and outputs, in an identical fixed format, all bands shifted 
# by some given constant

theFile1 = open(sys.argv[1], "r").readlines()
en1 = float(sys.argv[2])

for i in xrange(len(theFile1)):
  line1 = [float(k) for k in theFile1[i].split()]
  print "%4d%17.7f%15.7f%15.7f" % (line1[0],line1[1],line1[2],line1[3]),
  for j in xrange(len(line1)):
    if j >= 4:
      if (j - 4) % 2 == 0:
        print "%11.5f" % line1[j],
      else:
        print "%14.5f" % (float(line1[j])-en1),
  print' $i $FERMI > temp

  mv $i $i.original
  mv temp $i
done
