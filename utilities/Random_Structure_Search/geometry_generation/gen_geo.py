#!python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:26:39 2014

@author: lange & jank
"""

import os
import time
import random
#import subprocess as sp
import numpy as np

class aims_geometry:
  def __init__(self):
    self.cell = np.empty((3,3))
    self.cell.fill(np.nan)
    self.coords = []
    self.labels = []
    self.speciesLabels = []
    self.atomsPerSpecies = []
    self.ad_structs = []

  def get_inputs(self,inputs):
    f=open(inputs)
    content = f.readlines()
    f.close()
    for line in content:
      line = line.split('%')[0]
      word = line.split()
      if(word[0] == 'max_wiggle_len'):
       max_wiggle_len = np.float64(word[1])
      if(word[0] == 'thetax'):
        thetax = np.float64(word[1])
      if(word[0] == 'thetay'):
        thetay = np.float64(word[1])
      if(word[0] == 'thetaz'):
        thetaz = np.float64(word[1])
      if(options.isperiodic):
        if(word[0] == 'min_z_dist'):
          min_z_dist = np.float64(word[1])
        if(word[0] == 'max_z_height'):
          max_z_height = np.float64(word[1])
      else:
        min_z_dist = 2
        max_z_height = 2
      if(word[0] == 'min_atom_dist'):
        min_atom_dist = np.float64(word[1])
      if(word[0] == 'maxtries'):
        maxtries = np.float64(word[1])
    return max_wiggle_len, thetax, thetay, thetaz, min_z_dist, max_z_height, min_atom_dist, maxtries

  def read(self,fileName,min_z_dist):
    cell = np.empty((3,3))
    cell.fill(np.nan)
    ad_structs = []
    coords = []
    labels = []
    f=open(fileName,'r')
    content = f.readlines()
    f.close()
    k=0
    for line in content:
      line = line.split("#")[0]
      words = line.split()
      if len(words) == 0:
        continue
      if(options.isperiodic):
        if(words[0] == 'box_dims'):
          print("   --- bulk_geometry.in may NOT contain box dimensions in the periodic case.")
          quit()
        if(words[0] == "lattice_vector"):
          if(len(words)!=4):
            print(fileName+": Syntax error in line '"+line+"'")
            quit()
          cell[0][k] = words[1]
          cell[1][k] = words[2]
          cell[2][k] = words[3]
          k+=1
      else:
        if(words[0] == 'lattice_vector'):
          print("   --- bulk_geometry.in may NOT contain lattice vectors in the non-periodic case.")
          quit()
        if(words[0] == 'box_dims'):
          if(len(words)!=4):
            print(fileName+": Syntax error in line '"+line+"'")
            quit()
          cell.fill(0)
          cell[0][0] = words[1]
          cell[1][1] = words[2]
          cell[2][2] = words[3]

    if(options.isperiodic):
      if ((len(cell) != 3) and (len(cell) != 0)):
        raise Exception(fileName+": Must contain exactly 3 lattice vectors or none!")

    if (len(cell) == 3):
      self.cell = np.array(cell).T

    for line in content:
      line = line.split("#")[0]
      words = line.split()
      if len(words) == 0:
        continue
      if words[0] == "atom":
        if len(words) != 5:
          raise Exception(fileName+": Syntax error in line '"+line+"'")
        coords.append(list(map(float,words[1:4])))
        labels.append(words[4])
      if words[0] == "atom_frac":
        if (len(cell) != 3):
          raise Exception(fileName+": Must contain 3 lattice vectors to support fractional coordinates!")
        if len(words) != 5:
          raise Exception(fileName+": Syntax error in line '"+line+"'")
        coords.append(list(map(float,np.dot(self.cell,np.array(words[1:4])))))
        labels.append(words[4])

    pwd = os.getcwd()
    files = os.listdir(pwd)

    adfiles = []
    k=0
    for line in files:
      if(line.find('ad_geo')!=-1):
        adfiles.append(line)

    num_adfiles = len(adfiles)

    max_atoms = 0
    for line in files:
      if(line.find('ad_geo')!=-1):
        f=open(line,'r')
        content = f.readlines()
        f.close()
        for subline in content:
          if(line.find('atom')):
            max_atoms+=1

    ad_structs = np.empty((num_adfiles,max_atoms,3))
    ad_structs.fill(np.nan)
    num_adds = np.empty((num_adfiles))
    num_adds.fill(np.int(1))
    ad_layers = np.empty((num_adfiles+1,10,2))
    ad_layers.fill(0)
    ad_layers_present = [False for i in range(num_adfiles)]

    ad_name = [['Null' for j in range(max_atoms)] for i in range(num_adfiles)]
    family_name = ['NoName' for i in range(num_adfiles)]
    k=0
    m=0
    for line in files:
      if(line.find('ad_geo')!=-1):
        f=open(line,'r')
        content = f.readlines()
        f.close()
        l=0
        for subline in content:
          subline = subline.split("%")[0]
          words = subline.split()
          if(len(words)==0):
            continue
          if(words[0]=='#'):
            continue
          if(words[0]=='lattice_vector'):
            raise Exception(line+": May NOT contain any lattice vectors!")
          if(words[0]=='atom_frac'):
            raise Exception(line+": May NOT contain any fractional coordinates!")
          if(words[0]=='number_of_adds'):
            num_adds[m] = np.int(words[1])
            continue
          if(words[0]=='atom'):
            ad_structs[m][l][0] = words[1]
            ad_structs[m][l][1] = words[2]
            ad_structs[m][l][2] = words[3]
            ad_name[m][l] = words[4]
            l+=1
          if(words[0]=='family'):
            family_name[m]=words[1]
          if(words[0]=='layer_1'):
            if(not options.isperiodic):
              print("       ----   Layered structure generation NOT supported without periodic boundaries. --- ")
              quit()
            ad_layers_present[m]=True
            ad_layers[m,0,0] = words[1]
            ad_layers[m,0,1] = words[2]
          if(words[0]=='layer_2'):
            ad_layers[m,1,0] = words[1]
            ad_layers[m,1,1] = words[2]
          if(words[0]=='layer_3'):
            ad_layers[m,2,0] = words[1]
            ad_layers[m,2,1] = words[2]
          if(words[0]=='layer_4'):
            ad_layers[m,3,0] = words[1]
            ad_layers[m,3,1] = words[2]
          if(words[0]=='layer_5'):
            ad_layers[m,4,0] = words[1]
            ad_layers[m,4,1] = words[2]
          if(words[0]=='layer_6'):
            ad_layers[m,5,0] = words[1]
            ad_layers[m,5,1] = words[2]
        m+=1

    for i in range(num_adfiles):
      if(ad_layers_present[i]):
        if(np.sum(ad_layers[i,:,0])!=num_adds[i]):
           print('\n\n *****  ERROR: The number of atoms to be added and')
           print('      the number of atoms in the layers do not match.')
           quit()

    for i in range(1,10):
      for j in range(num_adfiles):
        ad_layers[j,i,:] = ad_layers[j,i,:] + ad_layers[j,i-1,:]

    ad_layers[:,-1,1] = 2 #min_z_dist

    if(options.isperiodic):
      self.coords = np.asarray(coords)
    else:
      self.coords = np.asarray(cell)

    self.labels = labels
    self.calcSpeciesInformation()
    self.ad_structs = ad_structs
    self.ad_name = ad_name
    self.max_atoms = max_atoms
    self.num_adds = num_adds
    self.family_name = family_name
    self.ad_layers = ad_layers
    self.num_adfiles = num_adfiles
    self.ad_layers_present = ad_layers_present

  def getCellVolume(self):
    a = self.cell[0]
    b = self.cell[1]
    c = self.cell[2]
    result = np.dot(np.cross(a,b),c)
    return result

  def getCellAreas(self):
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

class Manipulator:
  def wiggleCoords(ad_coords,ad_name,max_wiggle_len):
    new_ad_coords = np.empty((len(ad_coords[:,0]),len(ad_coords[0,:])))
    new_ad_coords.fill(np.nan)
    for i in range(len(ad_coords[:,0])):
      for j in range(len(ad_coords[0,:])):
        if(ad_name[i]!='Null'):
          myseed = time.time()
          random.seed(myseed)
          new_ad_coords[i,j] = ad_coords[i,j] - max_wiggle_len/2 * random.random()
    return new_ad_coords

  def rotateCoords(ad_coords,ad_name,thetax,thetay,thetaz):
    theta = np.empty((3))
    theta.fill(np.nan)

    rotx = np.empty((3,3))
    rotx.fill(np.nan)
    roty = np.empty((3,3))
    roty.fill(np.nan)
    rotz = np.empty((3,3))
    rotz.fill(np.nan)

    for i in range(3):
      myseed = time.time()
      random.seed(myseed)
      theta[i] = random.random()

    theta[0] = theta[0] * thetax
    theta[1] = theta[1] * thetay
    theta[2] = theta[2] * thetaz

    rotx = [[1, 0, 0],[0, np.cos(theta[0]), -np.sin(theta[0])], [0, np.sin(theta[0]), np.cos(theta[0])]]
    roty = [[np.cos(theta[1]), 0, np.sin(theta[1])],[0, 1, 0], [-np.sin(theta[1]), 0, np.cos(theta[1])]]
    rotz = [[np.cos(theta[2]), -np.sin(theta[2]), 0], [np.sin(theta[2]), np.cos(theta[2]), 0], [0, 0, 1]]
    
    totrot = np.dot(rotx,np.dot(roty,rotz))

    j=0
    natoms=0
    while(j<len(ad_name[:])):
      if(ad_name[j]!='Null'):
        natoms+=1
      j+=1

    xctr = np.sum(ad_coords[0:natoms,0])/natoms
    yctr = np.sum(ad_coords[0:natoms,1])/natoms
    zctr = np.sum(ad_coords[0:natoms,2])/natoms

    ctr = np.array([xctr,yctr,zctr])

    if(ad_name[i]!='Null'):
      ad_coords[i,:] = np.dot(totrot,ad_coords[i,:]-ctr)
    return ad_coords

  def centerCoords(str_coords,ad_coords,ad_name,min_z_dist,max_z_height,numspec,count):
    bmin = np.empty((3))
    bmin.fill(np.nan)
    bmax = np.empty((3))
    bmax.fill(np.nan)
    bulk_dims = np.empty((3))
    bulk_dims.fill(np.nan)
    new_ctr = np.empty((3))
    new_ctr.fill(np.nan)

    bmin[0] = np.amin(str_coords[:,0])
    bmin[1] = np.amin(str_coords[:,1])
    bmin[2] = np.amin(str_coords[:,2])
    bmax[0] = np.amax(str_coords[:,0])
    bmax[1] = np.amax(str_coords[:,1])
    bmax[2] = np.amax(str_coords[:,2])

    bulk_dims[0] = bmax[0] - bmin[0]
    bulk_dims[1] = bmax[1] - bmin[1]
    bulk_dims[2] = bmax[2] - bmin[2]

    if(options.isperiodic):
      for k in range(2):
        myseed = time.time()
        random.seed(myseed)
        new_ctr[k] = bmin[k] + bulk_dims[k] * random.random()
      myseed = time.time()
      random.seed(myseed)
      new_ctr[2] = bmin[2] - min_z_dist - max_z_height * random.random()
      if(structure.ad_layers_present[numspec]):
        for i in range(10):
          if(count<=structure.ad_layers[numspec,i,0]):
            max_layer_thick = structure.ad_layers[numspec,i,1]-structure.ad_layers[numspec,i-1,1]
            myseed = time.time()
            random.seed(myseed)
            new_ctr[2] = bmin[2] - structure.ad_layers[numspec,i-1,1] - max_layer_thick * random.random()
            break
    else:
      for k in range(3):
        myseed = time.time()
        random.seed(myseed)
        new_ctr[k] = bmin[k] + bulk_dims[k] * random.random()
      myseed = time.time()
      random.seed(myseed)
      if(structure.ad_layers_present[numspec]):
        for i in range(10):
          if(count<=structure.ad_layers[numspec,i,0]):
            myseed = time.time()
            random.seed(myseed)
            break

    for i in range(len(ad_coords[:,0])):
      ad_coords[i,:] = ad_coords[i,:] + new_ctr[:]
    return ad_coords, bmin[2]

  def checkDistances(all_new,tmp_new,min_atom_dist,min_z_dist,bmin):
    valid = 1
    for o in range(len(all_new[:,0,0])):
      for i in range(len(all_new[0,:,0])):
          cart_dist = np.linalg.norm(all_new[o,i,:]-tmp_new[i,:])
          if(cart_dist<min_atom_dist):
            valid = 0
          if(options.isperiodic):
            cart_dist = np.linalg.norm(all_new[o,i,:]-tmp_new[i,:]+structure.cell[0])
            if(cart_dist<min_atom_dist):
              valid = 0
            cart_dist = np.linalg.norm(all_new[o,i,:]-tmp_new[i,:]+structure.cell[1])
            if(cart_dist<min_atom_dist):
              valid = 0
            z_layer_dist = bmin-tmp_new[i,2]
            if(z_layer_dist<min_z_dist):
              valid = 0
    return valid

  def putAdsOnStructure(str_coords,ad_coords,ad_name,num_adds,family_name,thetax,thetay,thetaz,max_wiggle_len,min_z_dist,max_z_height,min_atom_dist,maxtries):
    maxdim = 0
    for i in range(len(ad_coords[:,0,0])):
      maxdim += num_adds[i]

    all_new = np.empty((maxdim,len(ad_coords[0,:,0]),3))
    all_new.fill(np.nan)
    all_name = [['Null' for j in range(len(structure.ad_name[0]))] for i in range(np.int(maxdim))]

    o=0; runs = 0
    for k in range(len(ad_coords[:,0,0])):
      for m in range(np.int(num_adds[k])):
        print("Adding %14s  Family:  Item %3.1i of %3.1i." % (family_name[k], m+1, np.int(num_adds[k])) )
        atom_set = 0
        while(atom_set==0):
          tmp_new = Manipulator.wiggleCoords(structure.ad_structs[k,:,:],structure.ad_name[k],max_wiggle_len)
          tmp_new = Manipulator.rotateCoords(tmp_new,structure.ad_name[k],thetax,thetay,thetaz)
          tmp_new, bmin = Manipulator.centerCoords(str_coords,tmp_new,ad_name[k],min_z_dist,max_z_height,k,m+1)
          atom_set = Manipulator.checkDistances(all_new,tmp_new,min_atom_dist,min_z_dist,bmin)
          runs+=1
          if(runs>maxtries):
            Manipulator.MaxtiresError(runs)
        all_new[o,:,:] = tmp_new[:,:]
        all_name[o] = structure.ad_name[k]
        o+=1
    return all_new, all_name

  def writeNewGeo(ad_new,ad_name,file):
    f=open(options.input_filename)
    content = f.readlines()
    f.close()

    f=open(file,'w')
    if(options.isperiodic):
      for line in content:
        f.write(line)
    for k in range(len(ad_new[:,0,0])):
      for i in range(len(ad_new[0,:,0])):
        if(ad_name[k][i]!='Null'):
          f.write("atom   %20.10f   %20.10f   %20.10f   %s\n" % (ad_new[k,i,0], ad_new[k,i,1], ad_new[k,i,2], ad_name[k][i]) )
    f.close()

  def MaxtiresError(tries):
    print(" ***** Error  -  Could not find any suitable space within %7.1i tries, aborting now." % tries)
    print(" *****        Check number of adatoms, layer thicknesses and molecule dimensions and retry.")
    quit()

if __name__ =='__main__':
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("-i", help="The Input Filename", action="store", dest="input_filename", default="bulk_geometry.in", type=str)
  parser.add_option("-o", help="The Output Filename", action="store", dest="output_filename", default="geometry.in", type=str)
  parser.add_option("--non-periodic", help="Option to trigger non periodic cluster output", action="store_false", dest="isperiodic", default=True)
  parser.add_option('--pars', help='The file containing parameters', action='store', dest='para_file', default='parameters')

  (options, args) = parser.parse_args()

  structure = aims_geometry()

  (max_wiggle_len, thetax, thetay, thetaz, min_z_dist, max_z_height, min_atom_dist, maxtries) = structure.get_inputs(options.para_file)

  structure.read(options.input_filename,min_z_dist)

  new_str, new_name = Manipulator.putAdsOnStructure( \
    structure.coords, \
    structure.ad_structs, \
    structure.ad_name, \
    structure.num_adds, \
    structure.family_name, \
    thetax, thetay, thetaz, \
    max_wiggle_len, \
    min_z_dist, \
    max_z_height, \
    min_atom_dist, \
    maxtries )
  Manipulator.writeNewGeo(new_str,new_name,options.output_filename)
