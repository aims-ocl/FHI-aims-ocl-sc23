#!/bin/sh
#===============================================================================
#
#          FILE:  crest_run.sh
# 
#         USAGE:
USAGE="Usage: crest_run.sh <crest_input_file> <crest_output_file>"
# 
#   DESCRIPTION:  Sets up MCR environment variables before running the CREST
#                 Matlab standalone executable, giving it the provided
#                 input file and output file.
# 
#       CREATED:  Wed Aug  6 16:12:08 IDT 2014
#        AUTHOR:  Ofer Sinai
#===============================================================================

#======================================================================
  ### CHANGE THIS VARIABLE TO POINT TO THE MCR ROOT PATH ###
  MCRROOT="/usr/local/matlab/R2012a"

  ### CHANGE THIS TO POINT TO THE CREST MATLAB STANDALONE EXECUTABLE
  CREST_BIN="/home/ofersin/apps/CREST_main"
#======================================================================

if [[ $# -lt 2 ]]; then
  echo $USAGE 1>&2
  exit 1
fi

crest_input_file=$(readlink -f $1)
crest_output_file=$(readlink -f $2)
# Make sure given input file exists
if [[ ! -f $crest_input_file ]]; then
  echo "$(basename $0): File $f not found"
  exit 1
fi

# Set up environment variables
echo "Setting up environment variables"
echo "---"
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
      MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export LD_LIBRARY_PATH;
export XAPPLRESDIR;
echo "LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};"

# Run CREST
${CREST_BIN} $crest_input_file $crest_output_file

exit 0
