#!/usr/bin/python
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys

print "This script has been hacked together to serve my (WPH's) needs, but is provided to the community all the same."
print "I'm still actively working on it, so it's probably pretty buggy right now."
print "Caveat computor."
print
print "To see a list of options (ex. changing the length or thickness of k-vectors), use the -h or --help flags."
print

# Have an option to select atom instead of offset, selectable via flags
skip = 0
OK = np.array([0.5,0.5,0.5]) # Origin for k-vectors (aka shift) (direct coordinates)
atom_scale = 1
kvec_scale = 1
atom_centered = -1
userDefinedDimMin = False
userDefinedDimMax = False
plotAtoms = True
linescale = 1
for i in xrange(1,len(sys.argv)):
	if skip != 0:
		skip -= 1
	else:
		if sys.argv[i] == "-k_shift":
			skip = 3
			OK = np.array([float(i) for i in sys.argv[i+1:i+4]]) # Origin for k-vectors (aka shift) (direct coordinates)
		elif sys.argv[i] == "-atom_centered":
			skip = 1
			atom_centered = int(sys.argv[i+1])
		elif sys.argv[i] == "-atom_scale":
			skip = 1
			atom_scale = float(sys.argv[i+1])
		elif sys.argv[i] == "-kvec_scale":
			skip = 1
			kvec_scale = float(sys.argv[i+1])
		elif sys.argv[i] == "-line_scale":
			skip = 1
			linescale = float(sys.argv[i+1])
		elif sys.argv[i] == "-min":
			skip = 1
			userDefinedDimMin = True
			minDimUser = float(sys.argv[i+1])
		elif sys.argv[i] == "-max":
			skip = 1
			userDefinedDimMax = True
			maxDimUser = float(sys.argv[i+1])
		elif sys.argv[i] == "-noatoms":
			skip = 0
			plotAtoms = False
		elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
			print "This tool is intended to aid in conceptual analysis of band structures.  It draws the starting and ending k-vectors of bands (i.e. the k-points in an output band command in control.in) on top of the geometry."
			print "Options:"
			print "-k_shift:          Shift the origin for the k-vectors (in direct coordinates) [default 0.5 0.5 0.5]"
			print "-atom_centered:    Specify an atom (as an integer starting with 1) to set as the origin of the k-vectors, instead of a point in space [default none]"
			print "-atom_scale:       Change the radii scale of the atoms [default 1]"
			print "-kvec_scale:       Change the length scale of the k-vectors [default 1]"
			print "                   By default, the longest k-vector is half the length of the longest direct lattice vector."
			print "                   For highly asymmetric cells, this may need to be changed."
			print "-line_scale:       Change the thickness of the k-vectors [default 1]"
			print "-min:              Specify the minimum coordinate of the bounding box [default set by matplotlib]"
			print "-max:              Specify the maximum coordinate of the bounding box [default set by matplotlib]"
			print "-noatom:           Do not plot the atoms, only the direct lattice and k-vectors"
			print "-h or --help:      This menu"
			sys.exit(0)
		else:
			print "Invalid flag provided, exiting."
			sys.exit(1)

#if len(sys.argv) >= 4:
#		OK = np.array([float(i) for i in sys.argv[1:5]]) # Origin for k-vectors (aka shift) (direct coordinates)
#else:
#		OK = np.array([0.5,0.5,0.5]) # Origin for k-vectors (aka shift) (direct coordinates)

O = np.array([0,0,0]) # Origin for vertices vectors (Cartesian coordinates)

# TODO
# Lookup table for atomic radii
# Label atoms

# Get lattice vectors and atoms from geometry.in
geometrydotin = open("geometry.in","r")
lat = []
atoms = []
species = []
alreadyCart = [] # Tracks whether the atom location is specified in Cartesian coordinates or direct coordinates
		 # Needed for plotting the positions of the atoms
for line in geometrydotin.readlines():
	line = line.split()
	if len(line) > 0:
		if line[0] == "lattice_vector":
			lat.append(np.array([float(i) for i in line[1:4]]))
		elif line[0] == "atom_frac":
			if line[4] not in species:
				species.append(line[4])
			index = species.index(line[4]) # This will determine the color
			atoms.append([np.array([float(i) for i in line[1:4]]),index])
			alreadyCart.append(False)
			if len(atoms) == atom_centered: # We've found the atom that the user wants to set as the origin of the k-vectors
				OK = atoms[-1][0]
		elif line[0] == "atom":
			if line[4] not in species:
				species.append(line[4])
			index = species.index(line[4]) # This will determine the color
			atoms.append([np.array([float(i) for i in line[1:4]]),index])
			alreadyCart.append(True)
			if len(atoms) == atom_centered: # We've found the atom that the user wants to set as the origin of the k-vectors
				OK = atoms[-1][0]

if len(lat) != 3:
	print "There were not three (and only three) lattice vectors in geometry.in, exiting."
	sys.exit()

# From lattice, get reciprocal lattice
rLat = []
tempVec = np.cross(lat[1],lat[2])
temp = np.dot(tempVec,lat[0])
rLat.append(2*np.pi*tempVec/temp)
tempVec = np.cross(lat[2],lat[0])
temp = np.dot(tempVec,lat[1])
rLat.append(2*np.pi*tempVec/temp)
tempVec = np.cross(lat[0],lat[1])
temp = np.dot(tempVec,lat[2])
rLat.append(2*np.pi*tempVec/temp)

# Get the maximum length of the (real) lattice vectors
maxLengthLat = 0
for i in lat:
	temp = np.dot(i,i)
	if temp > maxLengthLat:
		maxLengthLat = np.sqrt(temp)

# Get k-points from control.in
controldotin = open("control.in","r")
kpoints = []
names = [] # To make sure k-points aren't added more than once
for line in controldotin.readlines():
	line = line.split()
	# If line starts with "output band", it will contain the desired k-points
	if len(line) > 0 and line[0] == "output" and line[1] == "band":
		kpoint1 = [np.array([float(i) for i in line[2:5]]),line[9]]
		kpoint2 = [np.array([float(i) for i in line[5:8]]),line[10]]
		if kpoints == []: # if no k-points exist, add the first one you come across
			kpoints.append(kpoint1)
			names.append(kpoint1[1])
		else:
			if kpoint1[1] not in names:
				kpoints.append(kpoint1)
				names.append(kpoint1[1])
		if kpoint2[1] not in names:
			kpoints.append(kpoint2)
			names.append(kpoint2[1])

# Get the maximum length of the k-points
maxLengthKPoint = 0 
print "K-Points:"
for i in kpoints:
	print i[1], ": ", i[0]
	temp = np.dot(i[0],i[0])
	if temp > maxLengthKPoint:
		maxLengthKPoint = np.sqrt(temp)

# Calculate Cartesian coordinates of atoms
for index, item in enumerate(atoms):
	if alreadyCart[index] == False: # We need to convert the specified direct coordinates to Cartesian coordinates
		temp = O + item[0][0]*lat[0] + item[0][1]*lat[1] + item[0][2]*lat[2]
	else: # Already in Cartesian coordinates, can use the coordinates directly (heh)
		temp = O + item[0]
	atoms[index] = [item[0],temp,item[1]]

# Having read in k-points and calculated reciprocal lattice, now do conversion
OK = OK[0]*lat[0] + OK[1]*lat[1] + OK[2]*lat[2]
for index, item in enumerate(kpoints):
	temp = item[0][0]*rLat[0] + item[0][1]*rLat[1] + item[0][2]*rLat[2]
	kpoints[index] = [item[0],temp,item[1]]

# Begin the plotting stuff
mpl.rcParams['legend.fontsize'] = 10

# Lattice vectors for computational (real-space) cell
A = lat[0]
B = lat[1]
C = lat[2]

vertices = [O,O+A,O+B,O+C,O+A+B,O+A+C,O+B+C,O+A+B+C] # Lattice points for computational cell
# Connectivity of vertices points (i.e. which are neighboring)
connect = [[0,1],[1,4],[4,2],[2,0], # base
		[0,3],[1,5],[2,6],[4,7], # sides
		[3,5],[3,6],[5,7],[6,7]] # top

fig = plt.figure(figsize=plt.figaspect(1.0)*1.5) # Initially make the aspect ratio of the figure 1:1, so that the unit cell it isn't distorted
ax = Axes3D(fig)

# Form the computational cell
# Loop over all connected vertices points, draw a line between them
for pair in connect:
	ax.plot([vertices[pair[0]][0],vertices[pair[1]][0]],[vertices[pair[0]][1],vertices[pair[1]][1]],[vertices[pair[0]][2],vertices[pair[1]][2]],'k-')

# Now draw in all the atoms
if plotAtoms:
	colors = ["k","r","g","b","y","c","m"]
	if (len(colors) < len(species)):
		print "Not enough colors specified to plot all species, please add more!"
		sys.exit()
	for atom in atoms:
		ax.plot([atom[1][0]],[atom[1][1]],[atom[1][2]],'o', ms=atom_scale*20,color=colors[atom[2]])
	# Getting labels for 3D scatters is a pain in matplotlib, this hack does it
	for i in xrange(len(species)):
		ax.plot([], [], 'o', c=colors[i], label=species[i])

# We modify the kvec_scale so that the longest k-point gets plotted as a vector half the length
# of the longest reciprocal lattice vector
kvec_scale = kvec_scale * 0.5 * (maxLengthLat/maxLengthKPoint)

# Finally, plot k-vectors 
if kpoints != []:
	for item in kpoints[:]:
		kpoint_recip, kpoint, label = item
		ax.plot([OK[0],OK[0]+kvec_scale*kpoint[0]],[OK[1],OK[1]+kvec_scale*kpoint[1]],[OK[2],OK[2]+kvec_scale*kpoint[2]],linewidth=2*linescale)
		ax.text(OK[0]+kvec_scale*kpoint[0], OK[1]+kvec_scale*kpoint[1], OK[2]+kvec_scale*kpoint[2], label)

ax.legend()


# We now have to force matplotlib to plot things nicely, since it always wants to fill the entire window
# We've already set it to initially have a 1:1 aspect ratio initially
# Now manually reset the default bounding box to make it cubic
minDim = 9999999999999999999999999999999
minDim = min(ax.get_xlim3d()[0],minDim)
minDim = min(ax.get_ylim3d()[0],minDim)
minDim = min(ax.get_zlim3d()[0],minDim)
maxDim = -9999999999999999999999999999999
maxDim = max(ax.get_xlim3d()[1],maxDim)
maxDim = max(ax.get_ylim3d()[1],maxDim)
maxDim = max(ax.get_zlim3d()[1],maxDim)
if userDefinedDimMin == True:
	minDim = minDimUser
if userDefinedDimMax == True:
	maxDim = maxDimUser
print
print "The initial minimum and maximum coordinates are ", minDim, " and ", maxDim, " (matplotlib may play around with this.)"
ax.auto_scale_xyz([minDim,maxDim],[minDim,maxDim],[minDim,maxDim])

print
print "Right now, changing the window size will destroy the aspect ratio."
print "Fun fact:  To zoom in and out in matplotlib, hold down the right mouse button and move the mouse up and down."
print "Fun fact:  I don't know how Macbook mice work, so you're on your own there."

# CHANGE THIS TO CHANGE RESOLUTION.
mpl.rcParams['savefig.dpi'] =  250
print
print "The resolution for saving figures is set to ", mpl.rcParams['savefig.dpi'], " dpi." 

plt.show()
