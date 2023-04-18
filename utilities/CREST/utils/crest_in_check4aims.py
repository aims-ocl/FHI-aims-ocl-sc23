#!/usr/bin/python
#===============================================================================
#
#          FILE:  crest_in_check4aims.py
# 
#         USAGE:  crest_in_check4aims.py [crest_input]
# 
#   DESCRIPTION:  Goes over crest input file and checks what info it can vs. the
#                 AIMS input files in the directory
# 
#       CREATED:  Wed Apr  2 08:10:20 IDT 2014
#        AUTHOR:  Ofer Sinai
#===============================================================================

from __future__ import print_function
import sys
import os
import re

# Main action
def main(args=sys.argv[1:]):

  USAGE = 'Usage: crest_in_check4aims.py [crest_input]'

  # Find crest input file path
  if len(args) > 0:
    crest_input_file = os.path.abspath(args[0])
  else:
    crest_input_file = os.path.abspath('crest.in')
  # Find AIMS input files
  run_dir = os.path.dirname(crest_input_file)
  aims_control_in = os.path.join(run_dir, 'control.in')
  aims_geometry_in = os.path.join(run_dir, 'geometry.in')
  # Check existence of files
  for infile in [aims_control_in, aims_geometry_in]:
    if not os.path.isfile(infile):
      print('File ' + infile + ' not found', file=sys.stderr)
      sys.exit(1)
  if not os.path.isfile(crest_input_file):
    print('CREST input file not found. Calculating from aims only.')
    crest_input_file = None

  # Generic messages
  warning_msg = 'Warning: {0!s} may be wrong:'
  ok_msg = 'Value for flag {0!s} seems OK:'
  not_found_msg = 'Flag {0!s} not specified in CREST input'
  report_msg = '... in CREST input: {0:0.14f}\n'\
      + '... calculated from AIMS input: {1:0.14f}'
  
  # Check flags
  for flag in ['surface_cell_area', 'intrinsic_num_electrons']:
    if crest_input_file is not None:
      crest_value = find_flag_value(flag, crest_input_file)
    else:
      crest_value = 0
    if crest_value is not None:
      if flag == 'surface_cell_area':
        aims_value = calculate_surface_cell_area(aims_geometry_in)
      elif flag == 'intrinsic_num_electrons':
        aims_value = calculate_intrinsic_nelect(aims_control_in, aims_geometry_in)
      if abs(crest_value - aims_value) > 1e-8:
        print(warning_msg.format(flag))
      else:
        print(ok_msg.format(flag))
      print(report_msg.format(crest_value, aims_value))
    else:
      print(not_found_msg.format(flag))



def find_flag_value(flagname, file):
  '''Find value of flag in given file. If file does not exist return None'''
  fid = open(file, 'r')
  filetext = fid.read()
  fid.close()

  match = re.search('\n\s*' + flagname + '\s*=\s*(\S+)', filetext)
  if not match:
    return None
  else:
    return float(match.group(1))



def calculate_surface_cell_area(aims_geometry_in):
  '''Calculate surface area from first 2 lattice vectors'''
  fid = open(aims_geometry_in, 'r')
  filetext = fid.read()
  fid.close()

  lvsstr = re.findall('\n\s*lattice_vector\s+(\S+)\s+(\S+)\s+(\S+)', filetext)
  lvs = [  [  float(n) for n in v  ] for v in lvsstr  ]
  return abs(lvs[0][0] * lvs[1][1] - lvs[0][1] * lvs[1][0])



def calculate_intrinsic_nelect(aims_control_in, aims_geometry_in):
  '''Calculate intrinsic number of electrons from AIMS input'''
  # Get numbers of ion types from geometry input
  fid = open(aims_geometry_in, 'r')
  filetext = fid.read()
  fid.close()
  atom_line_matches = re.findall('\n\s*atom(_frac)?(\s+\S+){3}\s+(\S+)', filetext)
  num_atoms = {}
  for match_groups in atom_line_matches:
    if match_groups[2] not in num_atoms.keys():
      num_atoms[match_groups[2]] = 0
    num_atoms[match_groups[2]] += 1
  
  # Get nuclear charges associated with each type and sum
  fid = open(aims_control_in, 'r')
  filetext = fid.read()
  fid.close()
  lines = filetext.split('\n')
  atom_charges = {}
  line_num = 0
  cur_type=''
  while line_num < len(lines):
    words = lines[line_num].split()
    if len(words) > 0:
      if words[0] == 'species':
        cur_type = words[1]
      elif words[0] == 'nucleus' and cur_type != '':
        if cur_type == 'CP':
          cur_charge = 0
        else:
          cur_charge = float(words[1])
        atom_charges[cur_type] = cur_charge
        cur_type=''
    line_num += 1

  total_charge = 0
  for atom_name in num_atoms:
    if atom_name not in atom_charges.keys():
      aims_input_error_msg = 'Could not find charge of atom {0!s} in control.in'
      print(aims_input_error_msg.format(atom_name))
    else:
      total_charge += num_atoms[atom_name] * atom_charges[atom_name]
  return total_charge



# Standard boilerplate  
if __name__ == "__main__":
  main()

