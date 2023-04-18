#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script to convert a FHI-aims geometry.in file of a periodic structure
# from direct coordinates to fractional coordinates
# before printing all atoms are moved into the unit cell.
#
# Usage: ./real2frac.py geometry.in
# The output written to frac_geometry.in

import sys
import string
from numpy import *


def wrap_vec(v):
   relvec = dot(super_invmat,v)
   wrapvec = (relvec+1e-5) % 1.0 - 1e-5
   return dot(super_mat,wrapvec)
 
##########################################  
# Lese Datei ein 
##########################################  

if len(sys.argv) <= 1: 
  print 'Usage: ./real2frac.py geometry_file'  
  sys.exit()  
  
fnames = sys.argv[1]

print "Read geometry file %s "%(fnames)   

f=open(fnames)   
fdata = [] 
element = []  
lattice = []  
for line in f:   
   t = string.split(line) 
   if len(t) == 0:
    continue
   if t[0] == '#':
    continue
   if t[0] =='constrain_relaxation':
    continue 
   if t[0] == 'lattice_vector': 
    lattice += [( float(t[1]),float(t[2]), float(t[3]) )] 
   elif t[0] =='atom':
    fdata += [( float(t[1]),float(t[2]), float(t[3]) )]   
    element += [( str(t[4]) )] 
   else:
    continue

fdata = array(fdata)   
lattice = array(lattice)  
element = array(element)

##########################################  
# define transformation matrix 
##########################################
lattice_inv_mat = linalg.inv(lattice)
super_mat = lattice.transpose()
super_invmat = linalg.inv(super_mat)

##########################################  
# Gebe Datei aus 
##########################################
print "Write geometry file: frac_geometry.in"

outputfile=open("frac_geometry.in","w")
outputfile.write("""#
lattice_vector """+ ((" %.8f"*3)%tuple(lattice[0,:])) + """
lattice_vector """+ ((" %.8f"*3)%tuple(lattice[1,:])) + """
lattice_vector """+ ((" %.8f"*3)%tuple(lattice[2,:])) + """
#
""")
for i in range(0,len(fdata)):
  trans_fdata = wrap_vec(fdata[i])
  trans_fdata = dot(super_invmat,trans_fdata)
  outputfile.write("atom_frac" + ((" %.8f"*3)%tuple(trans_fdata)) + " " + element[i] + "\n")

outputfile.close()
