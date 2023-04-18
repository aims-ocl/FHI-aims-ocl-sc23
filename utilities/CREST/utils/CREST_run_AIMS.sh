#!/bin/bash
#===============================================================================
#
#          FILE:  CREST_run_AIMS.sh
# 
#         USAGE:  ./CREST_run_AIMS.sh Qsheet(e/cm^2) [divisor_for_Qsheet] [use_existing_iterdirs_upto]
#       OPTIONS:  - Divisor for the given charge (e.g. if the charge is spread
#                   over several point charges). Default 1.
#                 - Instruction to use existing iteration dirs, i.e., take
#                   results from previous iterations up to and including given
#                   number (default 0)
# 
#   DESCRIPTION:  The command CREST executes in order to run AIMS and translate
#                 the result to CREST's format.
# 
#       CREATED:  02/20/2014 12:01:02 PM IST
#       REVISED:  Wed Sep 10 11:40:12 IDT 2014
#        AUTHOR:  Ofer Sinai
#===============================================================================

set -o nounset                              # Treat unset variables as an error

# Definitions
aims2crest=/home/ofersin/scripts/aims2crest.sh
atom2crest_dir=atom2crest
crest_in=crest.in
AIMS_control_in=control.in
AIMS_geometry_in=geometry.in
AIMS_out=aims.out

# Make sure all files exist
for f in $aims2crest $crest_in $AIMS_control_in; do
  if [[ ! -f $f ]]; then
    echo "$(basename $0): Path   $f   not found (see definitions at the beginning of this file)" 1>&2
    exit 1
  fi
done

if [[ $# -lt 1 ]]; then
  echo "$(basename $0): Error: no sheet charge given" 1>&2
  exit 1
fi

# First argument must be sheet charge to insert into AIMS control.in. Given in e/cm^2
given_Qsheet=$(echo $1 | sed 's/+//g')
# Second (optional) argument is a divisor for the given charge
divisor=${2:-1}
# Third (optional) argument is an instruction to use existing iteration dirs
useOldCalcsUpto=${3:-0}

# Create file to keep track of iterations if necessary
if [[ $useOldCalcsUpto -gt 0 ]] && [[ ! -f "__LASTITER" ]]; then echo "0" > __LASTITER; fi

# Resolve current iteration
if [[ -f "__LASTITER" ]]; then
  current_iter_num=$[$(cat __LASTITER) + 1]
  # Update file
  echo "$current_iter_num" > __LASTITER
else
  # Find last numbered iteration dir
  last_iter_num=$(ls -d * | grep -oE "^iteration[0-9][0-9]([^0-9]|$)" | tail -1 | sed -r 's:[^0-9]|n0::g')
  if [[ $last_iter_num = "" ]]; then
    current_iter_num=1
  else
    current_iter_num=$[$last_iter_num + 1]
  fi
fi

# Two-digit iteration numbers
if [[ $current_iter_num -lt 10 ]]; then
  iter_num_str="0$current_iter_num"
else
  iter_num_str="$current_iter_num"
fi

# Generate iteration dir / load data as necessary
if [[ $current_iter_num -gt $useOldCalcsUpto ]]; then
  # Name and create iteration directory
  iter_dir="iteration$iter_num_str""_$(date +'%Y-%m-%d_%H-%M-%S')"
  mkdir $iter_dir
  
  # Extract area of surface cell from crest input file
  surface_cell_area_AA=$(grep "surface_cell_area" $crest_in | sed -r -e 's/ //g' -e 's/#.*$//g' -e 's/^.*=//g')     # \AA^2
  
  # Convert given Qsheet from surface charge to electronic charge
  Qsheet_e=$(python -c "print('{0:0.12f}'.format($given_Qsheet / $divisor * $surface_cell_area_AA * 1e-16))")     # e/cm^2 * AA^2 * cm^2/AA^2 = e
  
  # Substitute Qsheet_e into control.in
  sed -r "s:\S+\s+#xxQsheetxx:$Qsheet_e          #xxQsheetxx:g" $AIMS_control_in > $iter_dir/$AIMS_control_in
  # Copy geometry file to dir
  cp -af $AIMS_geometry_in $iter_dir/
  # If Qsheet is 0 remove "sheet" point-charge atoms from geometry entirely. Assumes these atoms are denoted "CP"
  if [[ $Qsheet_e = "0.000000000000" ]]; then
    sed -i -r "/ +CP *$/d" $iter_dir/$AIMS_geometry_in
  fi
  
  # Run AIMS
  echo "$(basename $0): Running AIMS with the given Qsheet = $given_Qsheet eV/cm^2"
  cd $iter_dir
  mpirun aims > $AIMS_out
  cd ..
  
  # Ensure AIMS run completed
  if [[ $(grep -c "nice day" $iter_dir/$AIMS_out) -eq 0 ]]; then
    echo "$(basename $0): AIMS run seems not to have completed. Stopping"
    exit 1
  fi

  # Save output from calculation in CREST-readable format
  $aims2crest $iter_dir
  
else
  iter_dir=$(ls -d * | grep -E "^iteration$iter_num_str([^0-9]|$)")

  if [[ $iter_dir = "" ]]; then
    echo "$(basename $0): Error: Dir for iteration no. $current_iter_num not found" 1>&2
    exit 1
  fi
fi

# Copy the dir to main dir
cp -r $iter_dir/$atom2crest_dir ./
