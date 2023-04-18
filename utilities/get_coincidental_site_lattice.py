#!/usr/bin/env python

#  Script that takes two different geometry files and find possible 
#  coincidental site lattices with respect to a small epsilon
#
#  > get_coincidental_site_lattice.py <file_1> <file_2> <epsilon>
#
#  Output is repitition and rotation angle of the lattices 

__author__ = 'blange'

import numpy as np
from optparse import OptionParser

import sys
def printf(format, *args):
    sys.stdout.write(format % args)

class structure:

   def __init__(self,nDim):
      self.cell = []
      self.nDim = nDim

   def read(self, fileName):
      print "Reading lattice vectors from ", fileName
      for line in file(fileName):
         line = line.split("#")[0]
         words = line.split()
         if len(words) == 0:
            continue
         if words[0] == "lattice_vector":
            if len(words) != 4:
               raise Exception(fileName+": Syntax error in line '"+line+"'")
            self.cell += [ np.array(map(float,words[1:4])) ]

      if len(self.cell) != 3:
         raise Exception(fileName+": Must contain exactly 3 lattice vectors")

      self.cell = np.asarray(self.cell)
      
   def printLattice(self):
      print "Lattice vectors:"
      for i in range(3):
         print self.cell[i,:]
      print
   
   def getDist(self,vec):
       x = vec[0]
       y = vec[1]
       z = vec[2]
       result = np.sqrt(x*x + y*y + z*z)
       return result
       
   def getVec(self,N1,N2,N3):
       x = N1 * self.cell[0,0] + N2 * self.cell[1,0] + N3 * self.cell[2,0]
       y = N1 * self.cell[0,1] + N2 * self.cell[1,1] + N3 * self.cell[2,1]
       z = N1 * self.cell[0,2] + N2 * self.cell[1,2] + N3 * self.cell[2,2]
       result = [x,y,z]
       return result
      
   def getEuler(self, vec1, vec2):
              
       #angle in xy plane
       gamma = 0.
       dotXY = (vec1[0]*vec2[0]+vec1[1]*vec2[1])
       normXY = np.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]) \
              * np.sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1])
       if (normXY > 1e-6):
          cosa = dotXY/normXY
          if (np.fabs(cosa - 1) < 1e-6): cosa = 1.0
          gamma = np.arccos(cosa) * 180 / np.pi
          
       #angle in xz plane
       beta = 0.
       dotXZ = vec1[0]*vec2[0]+vec1[2]*vec2[2]
       normXZ = np.sqrt(vec1[0]*vec1[0]+vec1[2]*vec1[2]) \
              * np.sqrt(vec2[0]*vec2[0]+vec2[2]*vec2[2])
       if (normXZ > 1e-6):
          cosa = dotXZ/normXZ
          if (np.fabs(cosa - 1) < 1e-6): cosa = 1.0
          beta = np.arccos(cosa) * 180 / np.pi
          
       #angle in yz plane
       alpha = 0.
       dotYZ = vec1[1]*vec2[1]+vec1[2]*vec2[2]
       normYZ = np.sqrt(vec1[1]*vec1[1]+vec1[2]*vec1[2]) \
              * np.sqrt(vec2[1]*vec2[1]+vec2[2]*vec2[2])
       if (normYZ > 1e-6):
          cosa = dotYZ/normYZ
          if (np.fabs(cosa - 1) < 1e-6): cosa = 1.0
          alpha = np.arccos(cosa) * 180 / np.pi
                
       return [alpha,beta,gamma]
          
          
   def getCoincidental(self,vec,eps):
       
       result = []
       
       dist = self.getDist(vec)
       N1Max = int(np.ceil(dist/self.getDist(self.getVec(1,0,0))))
       N2Max = int(np.ceil(dist/self.getDist(self.getVec(0,1,0))))
       N3Max = int(np.ceil(dist/self.getDist(self.getVec(0,0,1))))
       
       if (self.nDim < 3): N3Max = 0
       if (self.nDim < 2): N2Max = 0  
       
       for mult1 in range(0,N1Max+1):
          for mult2 in range(0,N2Max+1):
             for mult3 in range(0,N3Max+1):
                 myVec = self.getVec(mult1,mult2,mult3)
                 myDist = self.getDist(myVec)
                 
                 if (np.fabs(myDist-dist) < eps):
                     euler = self.getEuler(vec,myVec)
                     result += [np.array([mult1,mult2,mult3,myDist,
                                          euler[0],euler[1],euler[2]])]
                     
       return np.asarray(result)



# main
if __name__ =='__main__':

    parser = OptionParser("get_coincidental_site_lattice [options]")

    parser.add_option("--file1", dest="file1", help="First lattice file")
    parser.add_option("--file2", dest="file2", help="Second lattice file")
    parser.add_option("-e","--epsilon", dest="epsilon", default=0.05, 
                      help="delta for identical point dist") 
    parser.add_option("--NMax1", dest="NMax1", default=1,
                      help="maximal elongation along first lattice")
    parser.add_option("--NMax2", dest="NMax2", default=1,
                      help="maximal elongation along second lattice")
    parser.add_option("--NMax3", dest="NMax3", default=1,
                      help="maximal elongation along third lattice")
    parser.add_option("--nDim", dest="nDim", default=3,
                      help="dimension of the crystal: 1=x 2=xy 3=xyz")
            
    (options, args) = parser.parse_args()        

    NMax1 = int(options.NMax1)
    NMax2 = int(options.NMax2)
    NMax3 = int(options.NMax3)    
    
    nDim = int(options.nDim)
    if (nDim < 3): NMax3 = 0
    if (nDim < 2): NMax2 = 0
    struct1 = structure(nDim)
    struct2 = structure(nDim)
    
    struct1.read(options.file1)
    struct2.read(options.file2)
    
    struct1.printLattice()
    print
    struct2.printLattice()
    print
    
    printf('Coincidental lattices:\n')
    for mult1 in range(0,NMax1+1):
        for mult2 in range(0,NMax2+1):
            for mult3 in range(0,NMax3+1):
                vec = struct1.getVec(mult1,mult2,mult3)
                dist = struct1.getDist(vec)
                coincidental = struct2.getCoincidental(vec,float(options.epsilon))
                if (len(coincidental) > 0):
                   printf ('%ix%ix%i:\n',mult1,mult2,mult3)
                   for data in coincidental: 
                       printf ('%ix%ix%i with rotation (deg) ',data[0],data[1],data[2])
                       printf ('x=%6.2f y=%6.2f z=%6.2f\n',data[4],data[5],data[6])
                   print
                
                

