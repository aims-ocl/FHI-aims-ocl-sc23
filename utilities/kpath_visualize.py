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

# Future update:   Add arrows via https://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector

# Get lattice vectors from geometry.in
geometrydotin = open("geometry.in","r")
lat = []
for line in geometrydotin.readlines():
	line = line.split()
	if len(line) > 0:
		if line[0] == "lattice_vector":
			lat.append(np.array([float(i) for i in line[1:4]]))
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

# Get k-point paths from control.in
controldotin = open("control.in","r")
segments = []
kpoint = []
kpoint_name = []
for line in controldotin.readlines():
	line = line.split()
	# If line starts with "output band", it will contain the desired segments
	if len(line) > 0 and line[0] == "output" and line[1] == "band":
		kpoint1 = [np.array([float(i) for i in line[2:5]]),line[9]]
		kpoint1[0] = kpoint1[0][0]*rLat[0] + kpoint1[0][1]*rLat[1] + kpoint1[0][2]*rLat[2]
		kpoint2 = [np.array([float(i) for i in line[5:8]]),line[10]]
		kpoint2[0] = kpoint2[0][0]*rLat[0] + kpoint2[0][1]*rLat[1] + kpoint2[0][2]*rLat[2]
		segments.append([kpoint1[0],kpoint2[0],kpoint1[1] + "-" + kpoint2[1]])
		if kpoint1[1] not in kpoint_name:
			kpoint.append(kpoint1[0])
			kpoint_name.append(kpoint1[1])
		if kpoint2[1] not in kpoint_name:
			kpoint.append(kpoint2[0])
			kpoint_name.append(kpoint2[1])
		

# Begin the plotting stuff
mpl.rcParams['legend.fontsize'] = 10

# Lattice vectors for reciprocal cell (not Brillouin!)
O = np.array([0,0,0])
A = rLat[0]
B = rLat[1]
C = rLat[2]

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

# Plot k-vectors 
for item in segments[:]:
	ax.plot([item[0][0],item[1][0]],[item[0][1],item[1][1]],[item[0][2],item[1][2]],label=item[2],linewidth=2)

# Plot the k-point labels
for index, item in enumerate(kpoint):
	ax.text(item[0], item[1], item[2], kpoint_name[index])	

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
#if userDefinedDimMin == True:
#	minDim = minDimUser
#if userDefinedDimMax == True:
#	maxDim = maxDimUser
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
