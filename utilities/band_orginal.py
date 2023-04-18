#!/usr/bin/env python
# To run this script, ./band_original.py energy_shift filename

import sys
import matplotlib.pyplot as plt
import matplotlib
import math
from os import listdir
from os.path import isfile, join
import numpy as np
from numpy import pi

matplotlib.rcParams.update({'font.size': 40})
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman']
fig = plt.gcf()
fig.set_size_inches(10, 6)
linewidth = 1
black_bands = 1
special_line = 3
#energyshift = float(sys.argv[1])
ymax = 5
ymin = -2

inorghomo = []
orghomo = []
# orghomo.append(int(sys.argv[4]))
orglumo = []
inorglumo = []
energyshift = float(sys.argv[1])
filename = sys.argv[2]
# AEQTPbBr4

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
print "Lattice vectors:"
for j in range(3):
    print latvec[j]
print
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
#n_sample = 5
k_label = []
kpoint = []
for line in open('control.in'):
	if line.strip().startswith('output band'):
		words = line.strip().split()
		kpoint.append([float(i) for i in words[2:8]])
		n_sample = int(words[-3])
		k_label.append(words[-2:])
print 'K path'
print kpoint
band_len = []
print k_label
k_label_reduce = []
for k_pair in k_label:
	for i in range(len(k_pair)):
		k = k_pair[i]
		if len(k_label_reduce) == 0:
			k_label_reduce.append(k)
		elif i == 0 and k != k_label_reduce[-1]:
			k_label_reduce[-1] = k_label_reduce[-1] + "|" + k
		elif i == 0 and k == k_label_reduce[-1]:
			continue
		else:
			k_label_reduce.append(k)
print k_label_reduce
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
# print xvals
band_len_tot.append(xvals[-1][-1])

############################################
# band10**.out
############################################
onlyfiles = [f for f in listdir('.') if isfile(join('.', f))]
bandfiles = [f for f in onlyfiles if f.startswith('band') and '.out' in f and len(f) == 12]
print 'band output files'
print bandfiles
n_files = len(bandfiles)
bands_all_files = []

#energys = []
#bandfiles = ['band1001.out', 'band1002.out', 'band1003.out']
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
		if band_id in orghomo:
			#continue
			plt.plot(xvals[file_id], band, color='g', lw=special_line)
		elif band_id in inorghomo:
			#continue
			plt.plot(xvals[file_id], band, color='r', lw=special_line)
		elif band_id in orglumo:
			if file_id in solid:
				plt.plot(xvals[file_id], band, 'b', lw=special_line)
			else:
				plt.plot(xvals[file_id], band, 'b--', lw=special_line)
		elif band_id in inorglumo:
			if file_id in solid:
				plt.plot(xvals[file_id], band, 'r', lw=special_line)
			else:
				plt.plot(xvals[file_id], band, 'r:', lw=special_line)
		else:
			#continue
			plt.plot(xvals[file_id], band, color='k', lw=black_bands)
for kpoint_x in band_len_tot[1:]:
	plt.axvline(kpoint_x, color='k', linestyle = '--', lw=linewidth).set_dashes([5,5])
plt.axhline(0, color='k', linestyle = '--', lw=linewidth).set_dashes([5,5])
# plot mulliken point
 # blue for org., r for inorg
markers = ['o', 'o', 's', 's']
markersize = 10
markersizesquare = 10
marksersSize = [markersize, markersize, markersizesquare, markersizesquare]
markercolor = ['b', 'r', 'b', 'r']
plt.xticks([])
plt.xticks(band_len_tot, k_label_reduce, size=18)
plt.ylabel('Energy(ev)')
plt.axis([0, max([abs(i) for i in xvals[-1]]), ymin, ymax])
plt.savefig(filename, dpi = 300, bbox_inches='tight')
