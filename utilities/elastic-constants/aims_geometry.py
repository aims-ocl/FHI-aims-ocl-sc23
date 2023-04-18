#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:26:39 2014

@author: lange
"""

import numpy as np
import numpy.linalg as lg
import ase.data as pse


class aims_geometry:
  #Variables
  def __init__(self):
    self.cell = []
    self.cellConstraints = np.zeros([3,3])
    self.coords = []
    self.constraints = []
    self.options = []
    self.labels = []
    self.speciesLabels = []
    self.atomsPerSpecies = []
    self.isFrac = False
    self.atomtype={ 'H':0, 'He':1, 'Li':2, 'Be':3, 'B':4, 'C':5, 'N':6, \
       'O':7, 'F':8, 'Ne':9, 'Na':10, 'Mg':11, 'Al':12, 'Si':13, 'P':14, \
          'S':15, 'Cl':16, 'Ar':17, 'K':18, 'Ca':19, 'Sc':20, 'Ti':21, \
          'V':22, 'Cr':23, 'Mn':24, 'Fe':25, 'Co':26, 'Ni':27, \
          'Cu':28, 'Zn':29, 'Ga':30, 'Ge':31, 'As':32, 'Se':33, 'Br':34, 'Kr':35, 'Rb':36, \
          'Sr':37, 'Y':38, 'Zr':39, 'Nb':40, 'Mo':41, 'Tc':42, 'Ru':43, 'Rh':44, 'Pd':45, \
          'Ag':46, 'Cd':47, 'In':48, 'Sn':49, 'Sb':50, 'Te':51, 'I':52, 'Xe':53, 'Cs':54, \
          'Ba':55, 'La':56, 'Ce':57, 'Pr':58, 'Nd':59, 'Pm':60, 'Sm':61, 'Eu':62, 'Gd':63, \
          'Tb':64, 'Dy':65, 'Ho':66, 'Er':67, 'Tm':68, 'Yb':69, 'Lu':70, 'Hf':71, 'Ta':72, \
          'W':73, 'Re':74, 'Os':75, 'Ir':76, 'Pt':77, 'Au':78, 'Hg':79, 'Tl':80, 'Pb':81, \
          'Bi':82, 'Po':83, 'At':84, 'Rn':85, 'Fr':86, 'Ra':87, 'Ac':88, 'Th':89, 'Pa':90, \
          'U ':91, 'Np':92, 'Pu':93, 'Am':94, 'Cm':95, 'Bk':96, 'Cf':97, 'Es':98, 'Fm':99, \
          'Md':100, 'No':101, 'Lr':102}
          

    #"Covalent radii" taken from the 1st column in ABSTREE.H from the RasMol source code
    self.covradius=[0.320,1.600,0.680,0.352,0.820,0.680,0.752,0.728,0.720,1.120,0.972, \
              1.100,1.352,1.200,1.048,1.020,0.992,1.568,1.328,0.992,1.440,1.472, \
              1.328,1.352,1.352,1.340,1.328,1.500,1.150,1.448,1.220,1.168,1.208, \
              1.220,1.208,1.600,1.472,1.120,1.780,1.560,1.480,1.472,1.352,1.400, \
              1.448,1.500,1.592,1.688,1.628,1.460,0.620,1.700,1.400,1.700,1.672, \
              1.340,1.020,1.032,0.900,0.992,0.980,0.960,1.092,0.940,0.920,0.912, \
              0.888,0.880,0.872,0.860,0.848,0.780,0.680,0.700,0.720,0.880,0.680, \
              1.300,1.340,1.100,0.952,1.600,0.960,0.672,0.620,1.900,1.800,1.432, \
              1.180,1.020,0.888,0.968,0.952,0.928,0.920,0.912,0.900,0.888,0.880, \
              0.872,0.860,0.848,0.840];

    #"Atom radii" taken from Mendeleyev.c (2nd number in the line with the RGB color definitions)
    self.atomradius=[0.435,1.400,1.520,1.143,0.975,0.655,0.750,0.730,0.720,1.600,1.858, \
               1.605,1.432,1.176,1.060,1.020,0.990,1.900,2.262,1.976,1.655,1.476, \
               1.309,1.249,1.350,1.241,1.254,1.246,1.278,1.333,1.350,1.225,1.200, \
               1.160,1.140,2.000,2.470,2.151,1.824,1.616,1.432,1.363,1.367,1.353, \
               1.345,1.375,1.445,1.489,1.666,1.538,1.400,1.960,1.330,2.200,2.632, \
               2.171,1.873,1.824,1.836,1.829,1.809,1.804,1.984,1.818,1.800,1.795, \
               1.789,1.779,1.769,1.940,1.752,1.597,1.428,1.371,1.380,1.368,1.357, \
               1.387,1.442,1.502,1.728,1.750,1.460,1.460,1.450,1.430,2.500,2.140, \
               1.877,1.798,1.609,1.568,1.000,1.000,1.000,1.000,1.000,1.000,1.000, \
               1.000,1.000,1.000,1.000]

    self.vdwradius=[1.100,2.200,1.220,0.628,1.548,1.548,1.400,1.348,1.300,2.020,2.200, \
              1.500,1.500,2.200,1.880,1.808,1.748,2.768,2.388,1.948,1.320,1.948, \
              1.060,1.128,1.188,1.948,1.128,1.240,1.148,1.148,1.548,3.996,0.828, \
              0.900,1.748,1.900,2.648,2.020,1.608,1.420,1.328,1.748,1.800,1.200, \
              1.220,1.440,1.548,1.748,1.448,1.668,1.120,1.260,1.748,2.100,3.008, \
              2.408,1.828,1.860,1.620,1.788,1.760,1.740,1.960,1.688,1.660,1.628, \
              1.608,1.588,1.568,1.540,1.528,1.400,1.220,1.260,1.300,1.580,1.220, \
              1.548,1.448,1.980,1.708,2.160,1.728,1.208,1.120,2.300,3.240,2.568, \
              2.120,1.840,1.600,1.748,1.708,1.668,1.660,1.648,1.640,1.628,1.620, \
              1.608,1.600,1.588,1.580]

    self.colr=[0.800,0.643,0.700,0.643,0.900,0.350,0.200,0.800,0.700,0.643,0.600,0.600, \
         0.643,0.690,0.100,0.950,0.150,0.643,0.500,0.800,0.643,0.643,0.643,0.000, \
         0.643,0.518,0.643,0.707,0.950,0.643,0.900,0.643,1.000,0.643,0.500,0.643, \
         0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643, \
         0.643,0.643,0.643,0.643,0.500,0.643,0.643,0.643,0.643,0.800,0.643,0.643, \
         0.643,0.643,0.643,1.000,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643, \
         0.643,0.643,0.643,0.643,0.643,0.643,0.900,0.643,0.643,0.643,0.643,0.643, \
         0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643,0.643, \
         0.643,0.643,0.643,0.643,0.643,0.643,0.643];

    self.colg=[0.800,0.667,0.700,0.667,0.400,0.350,0.200,0.200,0.850,0.667,0.600,0.600, \
         0.667,0.769,0.700,0.900,0.500,0.667,0.500,0.800,0.667,0.667,0.667,0.800, \
         0.667,0.576,0.667,0.733,0.790,0.667,0.000,0.667,1.000,0.667,0.080,0.667, \
         0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667, \
         0.667,0.667,0.667,0.667,0.100,0.667,0.667,0.667,0.667,0.800,0.667,0.667, \
         0.667,0.667,0.667,0.843,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667, \
         0.667,0.667,0.667,0.667,0.667,0.667,0.800,0.667,0.667,0.667,0.667,0.667, \
         0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667,0.667, \
         0.667,0.667,0.667,0.667,0.667,0.667,0.667]

    self.colb=[0.800,0.678,0.700,0.678,0.000,0.350,0.800,0.200,0.450,0.678,0.600,0.700, \
         0.678,0.871,0.300,0.200,0.100,0.678,0.500,0.700,0.678,0.678,0.678,0.000, \
         0.678,0.653,0.678,0.746,0.014,0.678,1.000,0.678,0.300,0.678,0.120,0.678, \
         0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678, \
         0.678,0.678,0.678,0.678,0.500,0.678,0.678,0.678,0.678,0.000,0.678,0.678, \
         0.678,0.678,0.678,0.000,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678, \
         0.678,0.678,0.678,0.678,0.678,0.678,0.000,0.678,0.678,0.678,0.678,0.678, \
         0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678,0.678, \
         0.678,0.678,0.678,0.678,0.678,0.678,0.678]


  def read(self,fileName):
    cell = []
    coords = []
    labels = []
    # read lattice vectors
    c_lattice = 0
    for line in file(fileName):
      line = line.split("#")[0]
      words = line.split()
      if len(words) == 0:
        continue
      if words[0] == "lattice_vector":
        if len(words) != 4:
          raise Exception(fileName+": Syntax error in line '"+line+"'")
        cell.append(map(float,words[1:4]))
        c_lattice += 1
      if (words[0] == "constrain_relaxation" and len(words) == 2 and c_lattice > 0):
        if (words[1] == "x" and c_lattice <= 3):
          self.cellConstraints[(c_lattice-1),0] = 1
        if (words[1] == "y" and c_lattice <= 3):
          self.cellConstraints[(c_lattice-1),1] = 1
        if (words[1] == "z" and c_lattice <= 3):
          self.cellConstraints[(c_lattice-1),2] = 1
        if (words[1] == ".true." and c_lattice <= 3):
          self.cellConstraints[(c_lattice-1),:] = 1
      if (words[0] == "atom" or words[0] == "atom_frac"):
        c_lattice = 10


    # check that there are exact 3 lattice vectors
    if ((len(cell) != 3) and (len(cell) != 0)):
      raise Exception(fileName+": Must contain exactly 3 lattice vectors or none!")

    # cell hast to be transposed for Cell^frac_coord = abs_coord
    if (len(cell) == 3):
      self.cell = np.asarray(cell).T

    # read atom positions
    for line in file(fileName):
      line = line.split("#")[0]
      words = line.split()
      if len(words) == 0:
        continue
      # absolute coordinates
      if words[0] == "atom":
        if len(words) != 5:
          raise Exception(fileName+": Syntax error in line '"+line+"'")
        coords.append(map(float,words[1:4]))
        labels.append(words[4])
      # fractional coordinates
      if words[0] == "atom_frac":
        if (len(cell) != 3):
          raise Exception(fileName+": Must contain 3 lattice vectors to support fractional coordfinates!")
        if len(words) != 5:
          raise Exception(fileName+": Syntax error in line '"+line+"'")

        coords.append(self.cell.dot(np.array(map(float,words[1:4]))))
        labels.append(words[4])

      self.coords = np.asarray(coords)
      self.labels = labels
      self.calcSpeciesInformation()
    
    c_atom = 0
    nAtoms = len(self.labels)
    self.constraints = np.zeros([nAtoms,3])
    for line in file(fileName):
      words = line.split()
      if (len(words) < 2):
         continue
      if (words[0] == "atom" and len(words) == 5 ):
        c_atom += 1
      if (words[0] == "atom_frac" and len(words) == 5 ):
        c_atom += 1
      if (words[0] == "constrain_relaxation" and len(words) == 2 and c_atom != 0):
        if (words[1] == "x"):
          self.constraints[(c_atom-1),0] = 1
        if (words[1] == "y"):
          self.constraints[(c_atom-1),1] = 1
        if (words[1] == "z"):
          self.constraints[(c_atom-1),2] = 1
        if (words[1] == ".true."):
          self.constraints[(c_atom-1),:] = 1
      if (words[0] == "homogeneous_field" and len(words) == 4):
        self.options.append(line)
      if (words[0] == "set_vacuum_level" and len(words) == 2):
        self.options.append(line)


  def write(self, filename):
    f = open(filename, 'w')
    # Write lattice information if available
    if (len(self.cell) == 3):
      f.write("# Lattice vectors\n")
      for i in range(3):
        f.write("%-16s %13.6f %13.6f %13.6f\n" % ("lattice_vector", self.cell[0][i], self.cell[1][i], self.cell[2][i]))
        if (self.cellConstraints[i,0] and self.cellConstraints[i,1] and self.cellConstraints[i,2]):
          f.write("constrain_relaxation .true.\n")
        else:
          if (self.cellConstraints[i,0]):
            f.write("constrain_relaxation x\n")
          if (self.cellConstraints[i,1]):
            f.write("constrain_relaxation y\n")
          if (self.cellConstraints[i,2]):
            f.write("constrain_relaxation z\n")
      
    f.write("# Atomic coordinates\n")
    for i in range(self.getNTotalAtoms()):
      if (self.isFrac):
        f.write("%-16s %13.6f %13.6f %13.6f %-2s\n" % ("atom_frac", 
          self.coords[i][0], self.coords[i][1], self.coords[i][2], self.labels[i])) 
      else:
        f.write("%-16s %13.6f %13.6f %13.6f %-2s\n" % ("atom", 
          self.coords[i][0], self.coords[i][1], self.coords[i][2], self.labels[i]))
      if (self.constraints[i,0] and self.constraints[i,1] and self.constraints[i,2]):
        f.write("constrain_relaxation .true.\n")
      else:
        if (self.constraints[i,0]):
          f.write("constrain_relaxation x\n")
        if (self.constraints[i,1]):
          f.write("constrain_relaxation y\n")
        if (self.constraints[i,2]):
          f.write("constrain_relaxation z\n")
    for i in self.options:
      f.write(i)    
        
    f.close()
  
  def appendStructure(self,structure2):
    if (self.isFrac != structure2.isFrac):
      if self.isFrac:
        structure2.toFrac()
      else:
        structure2.toAbsolute()
    self.coords = np.concatenate((self.coords,structure2.coords))
    self.labels.extend(structure2.labels)
    self.constraints = np.concatenate((self.constraints,structure2.constraints))
    self.calcSpeciesInformation()
    
  def toFrac(self):
    if(not self.isFrac):
      cellInv = np.linalg.inv(self.cell)
      for i in range(self.getNTotalAtoms()):
        self.coords[i] = cellInv.dot(self.coords[i])
      self.isFrac = True      
  
  def toAbsolute(self):
    if(self.isFrac):
      for i in range(self.getNTotalAtoms()):
        self.coords[i] = self.cell.dot(self.coords[i])
      
  def getCellVolume(self):
    result = 0.0    
    if (len(self.cell) == 3):
      a = self.cell[0]
      b = self.cell[1]
      c = self.cell[2]
      result = np.dot(np.cross(a,b),c)
    
    return result

  def getCellAreas(self):
    ab = 0.0
    bc = 0.0
    ca = 0.0    
    if (len(self.cell) == 3):
      a = self.cell[0]
      b = self.cell[1]
      c = self.cell[2]
      ab = np.linalg.norm(np.cross(a,b))
      bc = np.linalg.norm(np.cross(b,c))
      ca = np.linalg.norm(np.cross(c,a))
    
    return ab, bc, ca

  def calcSpeciesInformation(self):
    nAtoms = len(self.labels)
    self.speciesLabels = []
    self.atomsPerSpecies = []
    nSpecies = 0
    for iAtom in range(nAtoms):
      if (nSpecies == 0):
        #First Atom defines new Species
        nSpecies += 1
        self.speciesLabels.append(self.labels[iAtom])
        self.atomsPerSpecies.append(1)
      else:
        found = False
        for iSpecies in range(nSpecies):
          if (self.labels[iAtom] == self.speciesLabels[iSpecies]):
            found = True
            self.atomsPerSpecies[iSpecies] += 1
            break
        if not found:
          nSpecies += 1
          self.speciesLabels.append(self.labels[iAtom])
          self.atomsPerSpecies.append(1)

  def getNTotalAtoms(self):
    return len(self.coords)

  def getNAtoms(self,iSpecies):
    return self.atomsPerSpecies[iSpecies]

  def getNAtomsLabel(self,label):
    result = 0
    for iSpecies in range(self.getNSpecies()):
      if self.speciesLabels[iSpecies][0] == label:
        result = self.getNAtoms(iSpecies)
        break

    return result

  def getNSpecies(self):
    return len(self.speciesLabels)
    
  def getSpeciesNumber(self,species):
      atomnumber = self.atomtype[species] + 1
      return atomnumber
    
  def getSpeciesLabel(self,iAtom):
    return self.labels[iAtom] 
    
  def getTotalMass(self):
    mTot   = 0.0
    nAtoms = self.getNTotalAtoms()
      
    for iAtom in range(nAtoms):
      species    = self.getSpeciesLabel(iAtom)
      atomNumber = self.getSpeciesNumber(species)
      mTot      += pse.atomic_masses[atomNumber]
          
    return mTot
      
  def getCenterOfMass(self):
    centerOfMass = np.zeros(3)
    nAtoms       = self.getNTotalAtoms()
    mTot         = self.getTotalMass()
      
    for iAtom in range(nAtoms):
      species       = self.getSpeciesLabel(iAtom)
      atomNumber    = self.getSpeciesNumber(species)
      m             = pse.atomic_masses[atomNumber]
      centerOfMass += m * self.coords[iAtom]
      
    return centerOfMass / mTot
         
  def toCenterOfMass(self):
    centerOfMass = self.getCenterOfMass()
    nAtoms       = self.getNTotalAtoms()
    
    for iAtom in range(nAtoms):
      self.coords[iAtom] -= centerOfMass
      
  def toCenterOfMassX(self):
    centerOfMass = self.getCenterOfMass()
    nAtoms       = self.getNTotalAtoms()

    centerOfMass[1] = 0.0
    centerOfMass[2] = 0.0    
    
    for iAtom in range(nAtoms):
      self.coords[iAtom] -= centerOfMass
      
  def toCenterOfMassY(self):
    centerOfMass = self.getCenterOfMass()
    nAtoms       = self.getNTotalAtoms()
    
    centerOfMass[0] = 0.0
    centerOfMass[2] = 0.0 
    
    for iAtom in range(nAtoms):
      self.coords[iAtom] -= centerOfMass

  def toCenterOfMassZ(self):
    centerOfMass = self.getCenterOfMass()
    nAtoms       = self.getNTotalAtoms()
    
    centerOfMass[0] = 0.0
    centerOfMass[1] = 0.0 
    
    for iAtom in range(nAtoms):
      self.coords[iAtom] -= centerOfMass  
  
  def shiftStructure(self,shift):
    nAtoms       = self.getNTotalAtoms()
    
    for iAtom in range(nAtoms):
      self.coords[iAtom] += shift
    
  def getInertiaTensor(self):
    centerOfMass = self.getCenterOfMass()
    self.shiftStructure(-centerOfMass)
    result       = np.zeros([3,3])
    nAtoms       = self.getNTotalAtoms()
    
    for iAtom in range(nAtoms):
      species       = self.getSpeciesLabel(iAtom)
      atomNumber    = self.getSpeciesNumber(species)
      m             = pse.atomic_masses[atomNumber]
      for i in range(3):
        for j in range(3):
          if (i == j):
            idx1 = (i+1) % 3
            idx2 = (i+2) % 3
            result[i,i] += m * (self.coords[iAtom,idx1]**2 
                              + self.coords[iAtom,idx2]**2)
          else:
            result[i,j] -= 0.5 * m * self.coords[iAtom,i] * self.coords[iAtom,j]
            result[j,i] -= 0.5 * m * self.coords[iAtom,j] * self.coords[iAtom,i]
          
      self.shiftStructure(centerOfMass)

    return result
    
  def getDiagonalInertiaTensor(self):
    tensor = self.getInertiaTensor()
    
    vals,vecs = lg.eig(tensor) 
    
    idx = vals.argsort()   
    vals = vals[idx]
    vecs = vecs[:,idx]
    
    return vals
    
    
  def transformCoords(self,mat):
    nAtoms       = self.getNTotalAtoms()
    
    for iAtom in range(nAtoms):
      self.coords[iAtom] = mat.dot(self.coords[iAtom]) 
      
  def mainAxesTransform(self):
      self.toCenterOfMass()
      IMoment = self.getInertiaTensor()
      
      print "Inertia Tensor for Transform"
      print IMoment
      print self.getCenterOfMass()
    
      vals,vecs = lg.eig(IMoment)
      
      idx = vals.argsort()   
      vals = vals[idx]
      vecs = vecs[:,idx]
      
      print "Trafotest"
      print vals
      print vecs
      print vecs.transpose()
      print lg.inv(vecs)
      print vecs.transpose().dot(IMoment.dot(vecs))    
      
      for i in range(3):
        print lg.norm(vecs[:,i])
      
      self.transformCoords(vecs.transpose())
      
      print "Inertia Tensor after Transform"
      IMoment = self.getInertiaTensor()
      print IMoment
      print self.getCenterOfMass()
      
         

  def writeDX(self,filename):
    f = open(filename,'w')
    nAtoms = self.getNTotalAtoms()
    f.write("object 1 class array type float rank 1 shape 3 items %i data follows\n" % nAtoms)
    for iAtom in range(nAtoms):
      f.write(" % 10.4f  % 10.4f  % 10.4f\n" % (self.coords[iAtom,0],
                 self.coords[iAtom,1],self.coords[iAtom,2]))
     
    f.write("object 2 class array type float rank 0 items %i data follows\n" % nAtoms)
    for iAtom in range(nAtoms):
      species = self.labels[iAtom]
      iSpecies = self.atomtype[species]
      f.write(" % 5.3f\n" % (self.atomradius[iSpecies]))
    f.write(" attribute \"dep\" string \"positions\"\n")

    f.write("object 3 class array type float rank 1 shape 3 item %i data follows\n" % nAtoms)
    for iAtom in range(nAtoms):
      species = self.labels[iAtom]
      iSpecies = self.atomtype[species]
      f.write(" %.3f %.3f %.3f\n" % (self.colr[iSpecies], 
                 self.colg[iSpecies], self.colb[iSpecies]))

    # Evaluate Bondings
    nBonds = 0
    bondFrom = []
    bondTo = []
    for iAtom in range(nAtoms-1):
      for jAtom in range(iAtom+1, nAtoms):
        d = np.linalg.norm(self.coords[jAtom] -
                            self.coords[iAtom])
        dCov = self.covradius[jAtom] + self.covradius[iAtom]
        if (d < dCov + 0.4):
          bondFrom.append(iAtom)
          bondTo.append(jAtom)
          nBonds += 1
             
    f.write("object 4 class array type int rank 1 shape 2 item %i data follows\n" % nBonds)
    for iBond in range(nBonds):
      f.write(" %6i %6i\n" % (bondFrom[iBond], bondTo[iBond]))
    f.write(" attribute \"element type\" string \"lines\"\n")
    f.write(" attribute \"ref\" string \"positions\"\n")
     
    f.write("object \"structure\" class field\n")   
    f.write("  component \"positions\" value 1\n")
    f.write("  component \"data\" value 2\n")
    f.write("  component \"colors\" value 3\n")
    f.write("  component \"connections\" value 4\n")
    f.write("end\n")

    f.close()





# main
if __name__ =='__main__':
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("-i", help="filename", action="store", dest="filename",
                    default="geometry.in", type=str)

  (options, args) = parser.parse_args()
  structure = aims_geometry()
  structure.read(options.filename)

  print "cell\n", structure.cell, "\n"
  print
  print "coordinates\n", structure.coords, "\n"
  print
  print "labels\n", structure.labels, "\n"
  
  print "nAtoms = ", structure.getNTotalAtoms(), "\n"

  print
  print ("Cell volume: %.6f A^3\n" % (structure.getCellVolume()))

  print ("Cell areas: \n xy: %.6f A^2\n yz: %.6f A^2\n xz: %.6f A^2"
    % (structure.getCellAreas()))
  
  structure.mainAxesTransform()
  print "coordinates\n", structure.coords, "\n"
  structure.write("result.in")