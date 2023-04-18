# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:52:55 2014

@author: Bjoern Lange

Visual structure comparison
"""

from optparse import OptionParser
from aims_geometry import aims_geometry
import numpy as np
import os
import shutil

def voigtToMat(index):
  if (index == 0):
    return 0,0
  elif (index == 1):
    return 1,1
  elif (index == 2):
    return 2,2
  elif (index == 3):
    return 1,2
  elif (index == 4):
    return 0,2
  elif (index == 5):
    return 0,1
  else:
    print "Error in voigtToMat"
    exit(1)      

# main
if __name__ =='__main__':
    parser = OptionParser()
    parser.add_option("-i", help="structure", action="store", 
          dest="structure", default="geometry.in", type=str)
    parser.add_option("-n", help="points in each direction", action="store", 
          dest="maxStruct", default=5, type=int)
    parser.add_option("--maxGamma", help="", action="store", 
          dest="maxGamma", default=0.02, type=float)
    parser.add_option("--e1", help="strain xx", action="store_true", 
          dest="e1", default=False)
    parser.add_option("--e2", help="strain yy", action="store_true", 
          dest="e2", default=False)
    parser.add_option("--e3", help="strain zz", action="store_true", 
          dest="e3", default=False)
    parser.add_option("--e4", help="strain yz", action="store_true", 
          dest="e4", default=False)
    parser.add_option("--e5", help="strain xz", action="store_true", 
          dest="e5", default=False)
    parser.add_option("--e6", help="strain xy", action="store_true", 
          dest="e6", default=False)
 
    (options, args) = parser.parse_args()
    structure = aims_geometry()
    structure.read(options.structure)
    structure.toFrac()
      
    # set up strain tensor  
    strain = np.zeros([3,3])
    
    if(options.e1): 
      strain[0,0] = 1
    if(options.e2): 
      strain[1,1] = 1
    if(options.e3): 
      strain[2,2] = 1
    if(options.e4): 
      strain[1,2] = 0.5
      strain[2,1] = 0.5
    if(options.e5): 
      strain[0,2] = 0.5
      strain[2,0] = 0.5
    if(options.e6): 
      strain[0,1] = 0.5
      strain[1,0] = 0.5
      
    
    print strain
    
    #Give equation d2E/dgamma^2 = ...
    coeff = np.zeros([6,6])
    for i in range(6):
      ir,ic = voigtToMat(i)
      for j in range(i,6):
        jr,jc = voigtToMat(j)
        if (strain[ir,ic] and strain[jr,jc]):
          coeff[i,j] = 1
          coeff[j,i] = 1
    
    formula = "1/V d^2E/d\gamma^2 = "
    first = True      
    for i in range(6):
      for j in range(i,6):
        if (coeff[i,j]):
          if (not first): 
            formula += " + "
          first = False 
          if (i==j):
            formula += str(int(coeff[i,j])) + "c_" + str(i+1) + str(j+1)
          else:
            formula += str(int(coeff[i,j]+coeff[j,i])) + "c_" + str(i+1) + str(j+1)

    print formula
            
    # Set up structures
    cellRef = 1.0 * structure.cell
    delta = options.maxGamma / options.maxStruct

    for i in range(-options.maxStruct,options.maxStruct+1):
      gamma = i * delta
      defTensor = np.identity(3) + gamma * strain
      directory = "gamma="+str(i*delta)
      if not os.path.exists(directory):
        os.makedirs(directory)
      filename = directory + "/geometry.in"
      structure.cell = defTensor.dot(cellRef)
      structure.write(filename)
      shutil.copyfile("control.in",directory + "/control.in")
      shutil.copyfile("submit.sh", directory + "/submit.sh")
