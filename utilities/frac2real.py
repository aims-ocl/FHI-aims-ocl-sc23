#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Script to convert a FHI-aims geometry.in file of a periodic structure in fractional coordinates to direct coordinates
#  before printing all atoms are moved into the unit cell.
#  Usage:
#  frac2real.py geometry.in
#  output written to new_geometry.in

import sys
import string
from numpy import *


def wrap_vec(v):
   wrapvec = (v+1e-5) % 1.0 - 1e-5
   return dot(super_mat,wrapvec)

if len(sys.argv) <= 1: 
  print 'Usage: ./frac2real.py geometry_file'  
  sys.exit()  
  
fnames = sys.argv[1]

print "Read geometry file: %s"%(fnames)   

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
   elif t[0] =='atom_frac':   
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
print "Write geometry file: new_geometry.in"

outputfile=open("new_geometry.in","w")
outputfile.write("""#
lattice_vector """+ ((" %.8f"*3)%tuple(lattice[0,:])) + """
lattice_vector """+ ((" %.8f"*3)%tuple(lattice[1,:])) + """
lattice_vector """+ ((" %.8f"*3)%tuple(lattice[2,:])) + """
#
""")
for i in range(0,len(fdata)):
  trans_fdata = wrap_vec(fdata[i])
  outputfile.write("atom" + ((" %.8f"*3)%tuple(trans_fdata)) + " " + element[i] + "\n")

outputfile.close()
