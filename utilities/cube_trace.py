#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 13:58:08 2014

@author: lange
"""
import numpy as np
import numpy.linalg as lg
from optparse import OptionParser

class cube_data:

  def __init__(self,filename):
    self.data = []
    self.atomZ = []
    self.coords = []
    self.grid = np.zeros(3)
    self.spacing = np.zeros((3,3))
    self.volOrigin = np.zeros(3)

    self.read(filename)

  def read(self,filename):
    lineCounter = 0

    #read Header
    for line in file(filename):
      lineCounter += 1
      words = line.split()
      if (lineCounter == 3):
        # Numer of atoms + center of volumetric data
        self.nAtoms = int(words[0])
        self.volOrigin = np.array(map(float,words[1:4]))
      if (lineCounter == 4):
        self.grid[0] = int(words[0])
        self.spacing[0] = map(float,words[1:4])
      if (lineCounter == 5):
        self.grid[1] = int(words[0])
        self.spacing[1] = map(float,words[1:4])
      if (lineCounter == 6):
        self.grid[2] = int(words[0])
        self.spacing[2] = map(float,words[1:4])

    # read coordinates and data
    lineCounter = 0
    for line in file(filename):
      lineCounter += 1
      words = line.split()
      if ((lineCounter > 6) and (lineCounter <= 6+self.nAtoms)):
        self.atomZ.append(int(words[0]))
        self.coords.append(map(float,words[2:5]))
      if (lineCounter > 6+self.nAtoms):
        for i in range(len(words)):
          self.data.append(float(words[i]))

    self.data = np.asarray(self.data)
    self.data = self.data.reshape((self.grid[0],self.grid[1],self.grid[2]))
    self.coords = np.asarray(self.coords)

  def print_data(self,centerAtom):
    center = self.coords[centerAtom]
    spacingInv = lg.inv(self.spacing.T)
    centerRel = np.round(np.dot(spacingInv,center-self.volOrigin))

    f=open("traceX.dat","w")

    for i in range(int(self.grid[0])):
      coord = np.array([i,centerRel[1],centerRel[2]])
      coordx = (np.dot(self.spacing.T,coord) + self.volOrigin)[0]
      f.write("%10.6f %10.6f\n" % (coordx,self.data[i,centerRel[1],centerRel[2]]))

    f.close()

    f=open("traceY.dat","w")

    for i in range(int(self.grid[1])):
      coord = np.array([centerRel[0],i,centerRel[2]])
      coordy = (np.dot(self.spacing.T,coord) + self.volOrigin)[1]
      f.write("%10.6f %10.6f\n" % (coordy,self.data[centerRel[0],i,centerRel[2]]))

    f.close()


    f=open("traceZ.dat","w")

    for i in range(int(self.grid[2])):
      coord = np.array([centerRel[0],centerRel[1],i])
      coordz = (np.dot(self.spacing.T,coord) + self.volOrigin)[2]
      f.write("%10.6f %10.6f\n" % (coordz,self.data[centerRel[0],centerRel[1],i]))

    f.close()


#main
if __name__ =='__main__':

  parser = OptionParser()
  parser.description = "This tool collects cube data and plots a trace along the three defininge cube vectors."
  parser.add_option("--plot", action="store_true", dest="plot", default=False)
  parser.add_option("-i", action="store", dest="fileIn", type=str)
  parser.add_option("--centerAtom", action="store", dest="centerAtom", type=int, default=0)
  (options, args) = parser.parse_args()

  data = cube_data(options.fileIn)

  data.print_data(options.centerAtom)