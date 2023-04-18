#!/bin/bash 
#===============================================================================
#
#          FILE:  aims2crest.sh [AIMS dir(s)]
# 
#         USAGE:  aims2crest.sh
# 
#   DESCRIPTION:  Takes output of aims calculation and creates a directory from
#                 which aims can get information
# 
#       CREATED:  03/02/2014 13:45:00 PM IDT
#        AUTHOR:  Ofer Sinai
#===============================================================================

set -o nounset                              # Treat unset variables as an error

USAGE="Usage: aims2crest.sh [AIMS dir(s)]"

# Definitions
sevcleanup=/home/ofersin/scripts/sevcleanup.py
target_dir=atom2crest
source_U=plane_average_lrp.cube
target_U=U.dat
source_DOS=Sorted_Eigenvalues.dat
target_DOS=DOS.dat
source_aims_output=aims.out
target_Ef0=Ef0.dat


# If no arguments use current dir
if [[ $# -eq 0 ]]; then
  calcdirs="."
else
  calcdirs="$@"
fi

for d in $calcdirs; do
  # Avoid uneccesary work
  if [[ ! -d $d ]]; then
    echo "$(basename $0): $d not found, skipping"
    continue
  fi
  if [[ ! -f $d/$source_U ]] || [[ ! -f $d/$source_DOS ]] || [[ ! -f $d/$source_aims_output ]]; then
    echo "$(basename $0): Some source files not found in $d, skipping"
    continue
  fi

  # Create dir
  if [[ ! -d $d/$target_dir ]]; then
    mkdir $d/$target_dir
  fi

  # Comment first line of U and copy to dir
  if [[ -f $d/$source_U ]]; then
    sed '1 s/^/# /g' $d/$source_U > $d/$target_dir/$target_U
  fi

  if [[ -f $d/$source_aims_output ]]; then
    # Extract Fermi level
    Ef0=$(grep '(Fermi level)' $d/$source_aims_output | tail -1 | cut -d: -f2 | sed 's/ //g')
    if [[ ! $Ef0 = "" ]]; then
      # Save Ef0 to file
      python -c "print('{0:.12f}'.format($Ef0))" > $d/$target_dir/$target_Ef0

      # Convert Sorted_Eigenvalues to usable DOS
      if [[ -f $d/$source_DOS ]]; then
        $sevcleanup --ef $Ef0 $d/$source_DOS > $d/$target_dir/$target_DOS
      fi
    fi
  fi
done
