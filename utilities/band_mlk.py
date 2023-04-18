#!/usr/bin/env python
# To run this script, provide energy shift
# ./*.py shift png_name species substates(0, 1, 2...)

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
from os import listdir
from os.path import isfile, join
import numpy as np
from numpy import pi
import operator
import matplotlib.patches as mpatches

########################################
# Figure settings
########################################
#matplotlib.rc('font', family='serif') 
#matplotlib.rc('font', serif='Times New Roman')
#matplotlib.rc('text', usetex='false') 
matplotlib.rcParams.update({'font.size': 34})

#matplotlib.rcParams.update({'font.size': 34})
matplotlib.rcParams['axes.linewidth'] = 2
#matplotlib.rcParams['font.serif'] = ['Times New Roman']
#matplotlib.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman']
#plt.rcParams['font.family'] = 'serif'

fig = plt.gcf()
fig.set_size_inches(10, 6)
linewidth = 1
black_bands = 1
special_line = 3
#energyshift = float(sys.argv[1])
ymax = 5
ymin = -2
pythoncolors = 'bgrcmy'
markersizeunit = 5

########################################
# Data
########################################
atoms = []
species_id = {} # map species to its index in geometry.in: Pb --> 1
energyshift = float(sys.argv[1])
filename = sys.argv[2]
spec_to_plot = sys.argv[3]
if len(sys.argv) == 5:
    substate = int(sys.argv[4]) # 0, s; 1, p; 2, d; ...
line_energy_id = 1
line_mlk_start_id = 4
line_atom_id = 3

############################################
# geometry.in
############################################
latvec = []
rlatvec = []
for line in open("geometry.in"):
    words = line.strip().split()
    if len(words) == 0:
        continue
    if words[0] == "lattice_vector":
        latvec.append([float(i) for i in words[1:4]])
    if line.strip().startswith("atom"):
    	atoms.append(words[-1])
print "Lattice vectors:"
for j in range(3):
    print latvec[j]
print "atoms:"
print atoms
species = list(set(atoms))
#species = ['C', 'Br', 'H', 'N', 'Pb', 'S'] # orders for AE4TPbBr4
print "species:"
print species
for i in range(len(atoms)):
	if atoms[i] not in species_id:
		species_id[atoms[i]] = []
        species_id[atoms[i]].append(i)
print "species_id:"
print species_id

species_color = {}
for s in range(len(species)):
	species_color[species[s]] = pythoncolors[s]
print "species_color:", species_color

#Calculate reciprocal lattice vectors
volume = (np.dot(latvec[0],np.cross(latvec[1],latvec[2])))
rlatvec.append(2*pi*np.cross(latvec[1],latvec[2])/volume)
rlatvec.append(2*pi*np.cross(latvec[2],latvec[0])/volume)
rlatvec.append(2*pi*np.cross(latvec[0],latvec[1])/volume)
print "Reciprocal lattice vectors:"
for j in range(3):
    print rlatvec[j]
print

############################################
# control.in
############################################
kpoint = []
for line in open('control.in'):
	if line.strip().startswith('output band'):
		words = line.strip().split()
		kpoint.append([float(i) for i in words[2:8]])
		n_sample = int(words[-3])
print 'K path'
print kpoint
band_len = []
xvals = []
for i in kpoint:
	kvec = []
	for j in range(3):
		kvec.append(i[j+3] - i[j])
	temp = math.sqrt(sum([k * k for k in list(np.dot(kvec, rlatvec))]))
	#print temp
	step = temp / (n_sample - 1)
	xval = []
	for i in range(n_sample):
		xval.append(i * step)
	xvals.append(xval)
	band_len.append(temp)
band_len_tot = []
for i in range(len(band_len)):
	if i == 0:
		band_len_tot.append(0)
	else:
		band_len_tot.append(sum(band_len[:i]))
#print band_len_tot
for i in range(len(xvals)):
	xvals[i] = [j + band_len_tot[i] for j in xvals[i]]
print 'x values'
# print xvals, [file 1][1 : 21 points], [file 2][1 : 21 points], etc....

############################################
# band10**.out
############################################
onlyfiles = [f for f in listdir('.') if isfile(join('.', f))]
bandfiles = [f for f in onlyfiles if f.startswith('band') and '.out' in f and len(f) == 12]
bandfiles.sort()
print 'band output files'
print bandfiles
n_files = len(bandfiles)
bandmlkfiles = [f for f in onlyfiles if f.startswith('bandmlk') and '.out' in f]
bandmlkfiles.sort()
print 'band mulliken output files'
print bandmlkfiles

for file_id in range(len(bandfiles)):
	file = bandfiles[file_id]
	print file
	energys = []
	f = open(file)
	for line in f:
		words = line.strip().split()
		energy = []
		occ_ene = words[4:]
		#print occ_ene
		for i in range(0, len(occ_ene), 2):
			energy.append(float(occ_ene[i+1]) - energyshift)
		energys.append(energy)
        # len(energy) = 171, number of states; len(energys) = 21, number of k points
		#band_points.append(len(k_index))
	bands = []
	for i in range(len(energy)):
		band = []
		for j in range(len(energys)):
			band.append(energys[j][i])
		bands.append(band)
	#print bands
	#for xval in xvals:
	for band_id in range(len(bands)):
		band = bands[band_id]
                plt.plot(xvals[file_id], band, color='k', lw=black_bands)
############################################
# mlk
############################################
	mlkfile = bandmlkfiles[file_id]
	f = open(mlkfile)
	print "open " + mlkfile
	currentMlk = []
	currentAtom = []
	for line in f:
		if line.strip() == '':
			continue
		if line.strip().startswith("k point number:"):
			currentK = int(line.strip().split()[3][:-1])
		elif line.strip().startswith("State") and len(line.strip().split()) == 2:
			# print "previous mlk", currentMlk
			words = line.strip().split()
			currentState = int(words[1])
			# currentMlkPerSpec = {}
			currentMlk = []
			currentAtom = []
		elif line.strip()[0].isdigit():
			words = line.strip().split()
			currentAtom.append(words[line_atom_id])
                        if len(sys.argv) == 5:
                            currentMlk.append(float(words[line_mlk_start_id + 1 + substate]))
                        else:
			    currentMlk.append(float(words[line_mlk_start_id]))
		if len(currentAtom) == len(atoms): # 12 atoms in total, current state is done, output or post process
			# print "current mlk", currentMlk
			spec_mlk = {}
			for i_spec in species:
				if i_spec not in spec_mlk:
					spec_mlk[i_spec] = 0
				i_spec_id = species_id[i_spec]
				for i in i_spec_id:
					spec_mlk[i_spec] += currentMlk[i]
			# stats = {'a':1000, 'b':3000, 'c': 100}
			# print currentState, sum(spec_mlk.values())
			# print spec_mlk
			# print currentMlk
			# maxSpec = max(spec_mlk.iteritems(), key=operator.itemgetter(1))[0]
                        # if currentK == 1:
                        #         print currentState, maxSpec
			plt.plot(xvals[file_id][currentK - 1], energys[currentK - 1][currentState - 1], species_color[spec_to_plot]+'o', markersize=markersizeunit*spec_mlk[spec_to_plot], markeredgecolor=species_color[spec_to_plot])

for kpoint_x in band_len_tot[1:]:
	plt.axvline(kpoint_x, color='k', linestyle = '--', lw=linewidth).set_dashes([5,5])
plt.axhline(0, color='k', linestyle = '--', lw=linewidth).set_dashes([5,5])

plt.yticks(range(ymin, ymax + 1), [])
plt.xticks([])
#plt.ylabel('Energy(ev)')
plt.axis([0, max([abs(i) for i in xvals[-1]]), ymin, ymax])
plt.savefig(filename, dpi = 300, bbox_inches='tight')
