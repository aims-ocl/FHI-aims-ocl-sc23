#!/usr/bin/python
#===============================================================================
#
#          FILE:  sevcleanup.py
# 
#         USAGE:  sevcleanup.py [--ef fermi_level] aims_dir(s) > STDOUT
# 
#   DESCRIPTION:  Takes "Sorted_Eigenvalues.dat" file and sums it up such that each
#                 energy appears only once and the weights are summed up
# 
#       CREATED:  02/02/2014 18:23:00 PM IDT
#        AUTHOR:  Ofer Sinai
#===============================================================================

from __future__ import print_function
import sys
import re
import os

# Main action
def main(args=sys.argv[1:]):
  # Basic command line argument parsing code.
  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  USAGE = 'Usage: sevcleanup.py [--ef fermi_level] aims_dir(s) > STDOUT'

  # Parse options
  Efermi = 0.0
  if '--ef' in args:
    opt_ind = args.index('--ef')
    opt_args = args[opt_ind:]
    if len(opt_args) < 2:
      print(USAGE, file=sys.stderr)
      sys.exit(1)
    try:
      Efermi = float(opt_args[1])
    except ValueError:
      print(USAGE, file=sys.stderr)
      sys.exit(1)
    del args[opt_ind:opt_ind+2]

  # At least 1 argument must be given
  if len(args) < 1:
    print(USAGE, file=sys.stderr)
    sys.exit(1)


  # Format output
  floatformat = '  {0:10.3f}           {1:2.9f}'

  for path in args:
    # Check path and see if it pertains to a file or a directory
    if os.path.isfile(path):
      filepath = path
    elif os.path.isdir(path):
      filepath = os.path.join(path, 'Sorted_Eigenvalues.dat')
    else:
      print('Path ' + path + ' not found', file=sys.stderr)
      continue

    # Check existence of file
    if not os.path.isfile(filepath):
      print('File ' + filepath + ' not found', file=sys.stderr)
      continue

    # Get data from file
    f = open(filepath, 'rU')
    lines = f.readlines()
    f.close()
    
    # Go over lines by order and sum weights. Use dictionary
    eigenvals_weights_hash = {}
    for line in lines:
      # Print comment lines as is
      match = re.search(r'^\s*#', line)
      if match:
        print(line.split('\n')[0])
        continue
    
      # Not a comment. Get data
      words = line.split()
      try:
        numbers = [   float(n) for n in words    ]
      except ValueError:
        # Comment out, print as is and continue
        print('# ' + line.split('\n')[0])
        continue

      if numbers[0] not in eigenvals_weights_hash.keys():
            eigenvals_weights_hash[numbers[0]] = 0
    
      eigenvals_weights_hash[numbers[0]] += numbers[1]
    
    # Print out all keys and values, shifted back by Ef
    for eigenval in sorted(eigenvals_weights_hash.keys()):
      weight = eigenvals_weights_hash[eigenval]
      print(floatformat.format(eigenval + Efermi, weight))



# Standard boilerplate  
if __name__ == "__main__":
  main()

