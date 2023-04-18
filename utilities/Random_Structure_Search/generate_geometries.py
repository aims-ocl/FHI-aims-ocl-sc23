#!/usr/bin/python
__author__ = 'jkloppenburg and blange'
import numpy as np
import numpy.linalg as linalg
import random as random
import time as time
from optparse import OptionParser
from aims_geometry import aims_geometry

# Specify the number of atoms to add and their chemical symbol here as a list,
# e.g.
nAdd = [4,2]
nAddLabels = ["Mo","C"]

def die(error_message):
	print ("Error:")
	print(error_message)
	exit()


class structure:

  def __init__(self,zmin,zmax,mindist,latRange):
    # the basis cell
    self.cell = []
    # the basis coords
    self.coords = []
    # the basis chemical symbols
    self.labels = []
    # coords of the additional atoms
    self.addAtoms = []
    # chemical symbols of the additional atoms
    self.addAtomLabels = []
    # minimal distance allowed between additional atoms
    self.mindist = mindist
    # minimal z value to distribute additional atoms
    self.zmin = zmin
    # maximal z value to distribute additional atoms
    self.zmax = zmax
    # lateral restriction in fractional coordinate
    self.latRange = latRange
    # box to distribute additional atoms
    self.positionCell = np.array([[0.0,0.0,0.0],
                                  [0.0,0.0,0.0],
                                  [0.0,0.0,zmax-zmin]])

    # Initialize random number generator
    self.initializeRandom()

    # create expanded list of additional atoms chemical symbols
    for i in range(len(nAdd)):
      for j in range(nAdd[i]):
        self.addAtomLabels.append(nAddLabels[i])

  # read basis geometry
  def read(self, fileName):
    structure = aims_geometry()
    structure.read(fileName)

    self.cell = structure.cell
    self.coords = structure.coords
    self.labels = structure.labels

    self.positionCell[0] = self.cell[0]
    self.positionCell[1] = self.cell[1]

    area = np.linalg.norm(np.cross(self.positionCell[0],self.positionCell[1]))
    vol = np.dot(np.cross(self.positionCell[0],self.positionCell[1]),self.positionCell[2])
    vol = vol * self.latRange**2
    nMaxCubes = vol/self.mindist**3
    print ("maximal Number of atoms for Volume: ", np.floor(nMaxCubes))
    print ("Number of atoms to generate: ", len(self.addAtomLabels))
    zopt = self.mindist**3 / area / self.latRange**2 * 5.0
    print ("Choose optimum z: ", zopt)
    self.positionCell[2][2] = zopt

  def write(self, fileName):

     f = open(fileName,"w")
     for i in range(len(self.cell)):
       #cell is transposed
       f.write("lattice_vector %10.6f %10.6f %10.6f\n" % \
                       (self.cell[0][i],self.cell[1][i],self.cell[2][i]))

     for i in range(len(self.coords)):
       f.write("atom %10.6f %10.6f %10.6f %2s\n" % \
                       (self.coords[i][0], self.coords[i][1],
                        self.coords[i][2], self.labels[i]))

     for i in range(len(self.addAtoms)):
       coords = self.positionCell.dot(self.addAtoms[i])
       coords[2] += self.zmin
       f.write("atom %10.6f %10.6f %10.6f %2s\n" % \
                       (coords[0], coords[1], coords[2], self.addAtomLabels[i]))

     f.close()

  def initializeRandom(self):
    seed = time.time()
    random.seed(seed)


  def getRandomPos(self):
    relPos = np.array([0.0,0.0,0.0])
    # random gives a number [0,1), so relative coordinates are natural
    relPos[0] = random.random() * self.latRange
    relPos[1] = random.random() * self.latRange
    relPos[2] = random.random()
    return relPos

  def checkDistance(self,coords):
    # False is no vioalation
    # True is violation of the minimal distance requirements
    result = False
    for i in range(len(self.addAtoms)):
      d = self.addAtoms[i] - coords
      # in periodic boundary conditions the next nearest atom
      # might be in the next cell so we add or substract one "cell"
      # if the distance is large than 0.5 in one direction

      for iCoord in range(3):
        if d[iCoord] > 0.5:
          d[iCoord] -= 1.0
        if d[iCoord] < -0.5:
          d[iCoord] += 1.0

      dist = linalg.norm(self.positionCell.dot(d))
      if (dist < self.mindist):
        result = True

    return result

  def mapCoordsIntoCell(self,coords):
    for i in range(len(coords)):
      if (coords[i] < 0):
        coords[i] += 1.0
      if (coords[i] >= 1.0):
        coords[i] -= 1.0

    return coords

  def constructRandom(self):
    # get a random position
    self.addAtoms = []
    self.addAtoms.append(self.getRandomPos())

    for i in range(1,len(self.addAtomLabels)):
      coords = self.getRandomPos()
      # if we violate the minimal distance requirements
      # displace the atom randomly until the requirement are fulfilled
      nTrials = 0
      while (self.checkDistance (coords) == True and nTrials < 50):
        coords = self.getRandomPos()
        nTrials += 1

      if (nTrials >= 50):
        print ("WARNING: minimum distance requirement could not be fulfilled.")
      self.addAtoms.append(coords)

    self.addAtoms = np.asarray(self.addAtoms)
    print self.addAtomLabels
    print self.addAtoms

# main
if __name__ =='__main__':
  parser = OptionParser()
  parser.description = "This tool constructs random structures " + \
    "for a random structure search."
  parser.add_option("-i", dest="file", help="basis geometry file",
    default = "geometry.in")
  parser.add_option("-o", dest="fileOut", help="resulting geometry file",
    default = "geometry.out")
  parser.add_option("--zmin", dest="zmin", help="minimal z value", type=float)
  parser.add_option("--zmax", dest="zmax", help="maximal z value", type=float)
  parser.add_option("--dmin", dest="dmin", help="minimal atomic distance",\
    type=float)
  parser.add_option("--latRange", dest="latRange", \
    help="lateral restiction 0 < latRange < 1.0",\
    default=1.0,
    type=float)

  (options, args) = parser.parse_args()

  if (options.zmin == None):
    die("zmin has to be specified")

  if (options.zmax == None):
    die("zmax has to be specified")

  if (options.dmin == None):
    die("dmin has to be specified")

  if ((options.latRange > 1.0) or (options.latRange < 0.0)):
    die("latRange must be in the interval [0,1]")

  struct = structure(options.zmin,options.zmax,options.dmin,options.latRange)
  struct.read(options.file)
  struct.constructRandom()
  struct.write(options.fileOut)