#!/usr/bin/python

import numpy as np
import random
import math
import urllib.request
import time

x_width = 3
y_width = 3

offset = 50

def makegeo():

  def getnewcoords(tmpatom,ii):
    myseed = time.time()
    random.seed(myseed)
    xr = random.random()
    yr = random.random()
    zr = random.random()
    xnew = xmin + xr * xlen
    ynew = ymin + yr * ylen
    if(ii<=si_atoms):
      zlen=1
      z_min_layerdist = 1.5
    elif(ii<=si_atoms+o_atoms):
      zlen=1
      z_min_layerdist = 2.5
    else:
      zlen=1
      z_min_layerdist = 6
    znew = zmin - zr * zlen - z_min_layerdist
    tmpatom[0] = xnew
    tmpatom[1] = ynew
    tmpatom[2] = znew
    print("Generating new %1s with %1s %1s %1s" % (ii,xnew,ynew,znew))
    return tmpatom

#  f=open("frac_geo","r")
#  self.endfiles = f.readlines()
#  f.close()

#  k=0
#  for line in self.endfiles:
#    name = line.split("    ")[0]
#    f = open(name,"r")
#    content = f.readlines()
#    f.close()
#    m=0
#    lj = 0
#    gj = 0
#    for cont in content:
#      cont = cont.split("%")[0]
#      entrs = cont.split()
#      if(len(entrs)==0):
#        continue
#      if(len(entrs)==4):
#        lat_data[lj][0] = entrs[1]
#        lat_data[lj][1] = entrs[2]
#        lat_data[lj][2] = entrs[3]
#        lj += 1
#      if(len(entrs)==5):
#        geo_data[gj][0] = entrs[1]
#        geo_data[gj][1] = entrs[2]
#        geo_data[gj][2] = entrs[3]
#        gj += 1


  # 4-layer H-terminated SiC:
  A = np.array([[1.55246281,-2.68894446,0.00000000],[1.55247732,2.68896959,0.00000000],[0.00000000,0.00000000,100.0000000]])

  lx = np.array([1.55246281,-2.68894446,0.00000000])
  ly = np.array([1.55247732,2.68896959,0.00000000])
  lz = np.array([0.00000000,0.00000000,100.0000000])

  atom_name = [ '' for i in range(9)]

  atom = np.empty((9,3))
  atom.fill(np.nan)

  atom[0][0] = 0.66661982
  atom[0][1] = 0.33157201
  atom[0][2] = 0.02414528
  atom_name[0] = 'Si' #np.array([0.66661982,0.33157201,0.02414528]) #Si

  atom[1][0] = -0.00004667
  atom[1][1] = -0.00175763
  atom[1][2] = 0.02070601
  atom_name[1] = 'C' #np.array([-0.00004667,-0.00175763,0.02070601]) #C

  atom[2][0] = 0.66669400
  atom[2][1] = 0.33149693
  atom[2][2] = 0.04418312
  atom_name[2] = 'C' #np.array([0.66669400,0.33149693,0.04418312]) #C

  atom[3][0] = 0.00003432
  atom[3][1] = -0.00183981
  atom[3][2] = 0.04970283
  atom_name[3] = 'Si' #np.array([0.00003432,-0.00183981,0.04970283]) #Si

  atom[4][0] = 0.66672093
  atom[4][1] = 0.33146954
  atom[4][2] = 0.07499307
  atom_name[4] = 'Si' #np.array([0.66672093,0.33146954,0.07499307]) #Si

  atom[5][0] = 0.00005678
  atom[5][1] = -0.00186249
  atom[5][2] = 0.06886178
  atom_name[5] = 'C' #np.array([0.00005678,-0.00186249,0.06886178]) #C

  atom[6][0] = 0.66675572
  atom[6][1] = 0.33143285
  atom[6][2] = 0.09403052
  atom_name[6] = 'C' #np.array([0.66675572,0.33143285,0.09403052]) #C

  atom[7][0] = 0.00009715
  atom[7][1] = -0.00190501
  atom[7][2] = 0.10011432
  atom_name[7] = 'Si' #np.array([0.00009715,-0.00190501,0.10011432]) #Si

  atom[8][0] = 0.00010097
  atom[8][1] = -0.00190658
  atom[8][2] = 0.11511007
  atom_name[8] = 'H' #np.array([0.00010097,-0.00190658,0.11511007]) #H

  f = open("tmp","w")

  f.write("lattice_vector  %1s  %1s  %1s\n" % (x_width*lx[0], y_width*lx[1], lx[2]))
  f.write("lattice_vector  %1s  %1s  %1s\n" % (x_width*ly[0], y_width*ly[1], ly[2]))
  f.write("lattice_vector  %1s  %1s  %1s\n" % (x_width*lz[0], y_width*lz[1], lz[2]))
  f.write("constrain_relaxation .true.\n\n")

  for n in range(0,x_width):
    for m in range(0,y_width):
      for k in range(9):
        f.write("atom  %1s  %1s  %1s  %1s\n" % (n*lx[0]+m*ly[0]+lz[0]+A.T.dot(atom[k])[0], n*lx[1]+m*ly[1]+lz[1]+A.T.dot(atom[k])[1], n*lx[2]+m*ly[2]+A.T.dot(atom[k])[2]+offset, atom_name[k]))
        if(k>3):
          f.write("constrain_relaxation .true.\n")
  f.close()

  f = open("tmp","r")
  content = f.readlines()
  f.close()

  xmin = 1000.0
  xmax = 0.0
  ymin = 1000.0
  ymax = 0.0
  zmin = 1000.0
  for line in content:
    a = line.split("  ")[0]
    if(a=="atom"):
      for n in range(x_width):
        for m in range(y_width):
          for k in range(9):
            x = n*lx[0]+m*ly[0]+lz[0]+A.T.dot(atom[k])[0]
            y = n*lx[1]+m*ly[1]+lz[1]+A.T.dot(atom[k])[1]
            z = n*lx[2]+m*ly[2]+A.T.dot(atom[k])[2]+offset
            if(x>xmax):
              xmax = x
            if(x<xmin):
              xmin = x
            if(y>ymax):
              ymax = y
            if(y<ymin):
              ymin = y
            if(z<zmin):
              zmin = z

  xlen = xmax - xmin
  ylen = ymax - ymin

  si_atoms = 5 # Attention: take one less because the first is generated extra
  c_atoms = 20
  o_atoms = 9
  adatoms = si_atoms + o_atoms + c_atoms

  min_cart_dist = 1.45

  si = [[0 for i in range(3)] for j in range(adatoms)]
  tmpatom = [0 for i in range(3)]

  getnewcoords(tmpatom,0)
  si[0][0]=tmpatom[0]
  si[0][1]=tmpatom[1]
  si[0][2]=tmpatom[2]

  for i in range(1,adatoms):
    si_set = 0
    k=0
    while(si_set==0):
      getnewcoords(tmpatom,i)
      k+=1
      if(i<=si_atoms):
        min_cart_dist = 1.8
      elif(i<=si_atoms+o_atoms):
        min_cart_dist = 1.7
      else:
        min_cart_dist = 1.6
      if(k>10000):
        print("Count not find a space for the next atom in 10000 steps, aborting")
#        quit()
      if(tmpatom[0]>xlen):
        print("x >> xlen")
        quit()
      if(tmpatom[1]>ylen):
        print("y >> ylen")
        quit()
      si[i][0]=tmpatom[0]
      si[i][1]=tmpatom[1]
      si[i][2]=tmpatom[2]
      si_set = 1
      for j in range(0,i):
        cart_dist = math.sqrt(math.pow((si[j][0]-tmpatom[0]),2) + math.pow((si[j][1]-tmpatom[1]),2) + math.pow((si[j][2]-tmpatom[2]),2) )
        if( cart_dist <= min_cart_dist):
          si_set = 0
        cart_dist = math.sqrt(math.pow((si[j][0]-tmpatom[0]+lx[0]),2) + math.pow((si[j][1]-tmpatom[1]+lx[1]),2) + math.pow((si[j][2]-tmpatom[2]),2) )
        if( cart_dist <= min_cart_dist):
          si_set = 0
        cart_dist = math.sqrt(math.pow((si[j][0]-tmpatom[0]+ly[0]),2) + math.pow((si[j][1]-tmpatom[1]+ly[1]),2) + math.pow((si[j][2]-tmpatom[2]),2) )
        if( cart_dist <= min_cart_dist):
          si_set = 0

  f = open("geometry.in","w")

  for line in content:
    f.write(line)

  for i in range(0,adatoms):
    if(i<=si_atoms):
      line = "atom  " + repr(si[i][0]) + "  " + repr(si[i][1]) + "  " + repr(si[i][2]) + "  Si\n"
      f.write(line)
    elif(i<=si_atoms+o_atoms):
      line = "atom  " + repr(si[i][0]) + "  " + repr(si[i][1]) + "  " + repr(si[i][2]) + "  O\n"
      f.write(line)
    else:
      line = "atom  " + repr(si[i][0]) + "  " + repr(si[i][1]) + "  " + repr(si[i][2]) + "  C\n"
      f.write(line)
  f.close()

makegeo()
