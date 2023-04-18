#!/usr/bin/env python
#
#  Script to plot band structure and DOS calculated with FHI-aims. Requires the control.in/geometry.in as well as the output
#  of the calculation to be in the same directory ...
#
#  To achieve labelling of the special points along the band structure plot, add two arguments to the "output band"
#  command in the control.in, using the following syntax:
#
#  output band <start> <end> <npoints> <starting_point_name> <ending_point_name>
#
#  Example: to plot a band with 100 points from Gamma to half way along one of the reciprocal lattice vectors, write (in control.in)
#
#  output band 0.0 0.0 0.0 0.5 0.0 0.0 100 Gamma <End_point_name>
#

from pylab import *
from numpy.linalg import *
from numpy import dot,cross,pi
from scipy.interpolate import spline

import os,sys

###########
# OPTIONS #
###########

# The following settings are the "standard" settings:
print_resolution = 250		# The DPI used for printing out images
default_line_width = 1		# Change the line width of plotted bands and k-vectors, 1 is default
font_size = 12                  # Change the font size.  12 is the default.
should_spline = False		# Turn on spline interpolation for band structures
spline_factor = 10		# If spline interpolation turned on, the sampling factor (1 is the original grid)
output_x_axis = True		# Whether to output the x-axis (e.g. the e=0 line) or not
color_array = 'rgbycmk'         # The colors used when outputting multiple band structures, specified as matplotlib.colors builtins
# The following are settings that better suited for publications/presentations:
#print_resolution =  500	 # The DPI used for printing out images
#default_line_width = 3		 # Change the line width of plotted bands and k-vectors, 1 is default
#font_size = 36                  # Change the font size.  12 is the default.
#should_spline = True		 # Turn on spline interpolation for band structures
#spline_factor = 10		 # If spline interpolation turned on, the sampling factor (1 is the original grid)
#output_x_axis = True		 # Whether to output the x-axis (e.g. the e=0 line) or not
#color_array = 'rbycgmk'         # The colors used when outputting multiple band structures, specified as matplotlib.colors builtins

########################

nPlots = int(sys.argv[1])
# This is cumbersome in practice, change to input file?
if len(sys.argv) <=2 or nPlots <= 0:
    print "The format for using aimsplot_compare.py is:"
    print "aimsplot_compare.py N_PLOTS DIRECTORY TITLE ENERGY_OFFSET ... yMin yMax"
    print "where N_PLOTS is the number of band structures to compare (max 7)"
    print "after which each band structure has a triplet DIRECTORY TITLE ENERGY_OFFSET, containing, respectively,"
    print "DIRECTORY - the directory for the aims run of the band structure,"
    print "TITLE - the legend title for the band structure, and"
    print "ENERGY_OFFSET -an energy offset for offseting vertical alignment."
    print "The last two parameters, yMin and yMax, are optional and set the minimum and maximum energy"
    print "values (i.e. y-axis)"
    print "NOTE:  This script assumes all geometry.in and k-path's are the same.  It will still run otherwise, but it will"
    print "likely look ugly."
    sys.exit()
if nPlots > 7:
    print "At max 7 band structures are supported.  No special reason; only seven colors were specified"
    print "by the Google search I did for matplotlib colors." 
    sys.exit()

########################

print "Plotting bands for FHI-aims!"
print "============================"
print

matplotlib.rcParams['lines.linewidth'] = default_line_width

directory = []
plotName = []
energy_offset = []
for i in xrange(nPlots):
    directory.append(sys.argv[3*i+2])
    plotName.append(sys.argv[3*i+3])
    energy_offset.append(float(sys.argv[3*i+4]))

CUSTOM_YLIM = False
if len(sys.argv) >= 4+3*nPlots:
    CUSTOM_YLIM = True
    ylim_lower = float(sys.argv[3*nPlots+2])
    ylim_upper = float(sys.argv[3*nPlots+3])

# Information from the files
PLOT_BANDS = [False for i in xrange(nPlots)]
PLOT_SOC = [False for i in xrange(nPlots)]  # This is needed because there will only be one "spin" channel output,
                           # but collinear spin may (or may not) be turned on, so the "spin 
		           # collinear" setting needs to be overridden
latvec = [[] for i in xrange(nPlots)]
species = [[] for i in xrange(nPlots)]
max_spin_channel = [1 for i in xrange(nPlots)]
band_segments = [[] for i in xrange(nPlots)]
band_totlength = [0.0 for i in xrange(nPlots)] # total length of all band segments
rlatvec = [[] for i in xrange(nPlots)]
thePlot = []

ax_bands = subplot(1,1,1)
for i in xrange(nPlots):
    print
    print "Reading lattice vectors from " + directory[i] + "/geometry.in ..."

    for line in file(directory[i] + "/geometry.in"):
        line = line.split("#")[0]
        words = line.split()
        if len(words) == 0:
            continue
        if words[0] == "lattice_vector":
            if len(words) != 4:
                raise Exception("geometry.in: Syntax error in line '"+line+"'")
            latvec[i] += [ array(map(float,words[1:4])) ]

    if len(latvec[i]) != 3:
        raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

    latvec[i] = asarray(latvec[i])

    print "Lattice vectors:"
    for j in range(3):
        print latvec[i][j,:]
    print

    #Calculate reciprocal lattice vectors                                                                                                
    volume = (np.dot(latvec[i][0,:],np.cross(latvec[i][1,:],latvec[i][2,:])))
    rlatvec[i].append(array(2*pi*cross(latvec[i][1,:],latvec[i][2,:])/volume))
    rlatvec[i].append(array(2*pi*cross(latvec[i][2,:],latvec[i][0,:])/volume))
    rlatvec[i].append(array(2*pi*cross(latvec[i][0,:],latvec[i][1,:])/volume))
    rlatvec[i] = asarray(rlatvec[i])

    #rlatvec = inv(latvec) Old way to calculate lattice vectors
    print "Reciprocal lattice vectors:"
    for j in range(3):
        print rlatvec[i][j,:]
    print

    ########################

    print "Reading information from " + directory[i] + "/control.in ..."

    for line in file(directory[i] + "/control.in"):
        words = line.split("#")[0].split()
        nline = " ".join(words)

	if nline.startswith("spin collinear") and not PLOT_SOC[i]:
            max_spin_channel[i] = 2

        if nline.startswith("calculate_perturbative_soc") or nline.startswith("include_spin_orbit") or nline.startswith("include_spin_orbit_sc"):
    	    PLOT_SOC[i] = True
            max_spin_channel[i] = 1

        if nline.startswith("output band "):
            if len(words) < 9 or len(words) > 11:
                raise Exception("control.in: Syntax error in line '"+line+"'")
            start = array(map(float,words[2:5]))
            end = array(map(float,words[5:8]))
            length = norm(dot(rlatvec[i],end) - dot(rlatvec[i],start))
            band_totlength[i] += length
            npoint = int(words[8])
            startname = ""
            endname = ""
            if len(words)>9:
                startname = words[9]
            if len(words)>10:
                endname = words[10]
            band_segments[i] += [ (start,end,length,npoint,startname,endname) ]

    #######################

    if PLOT_SOC[i]:
        max_spin_channel[i] = 1    

    #######################

    print "Plotting %i band segments..."%len(band_segments[i])

    if output_x_axis:
        ax_bands.axhline(0,color=(0.,0.,0.),linestyle=":")

    prev_end = band_segments[i][0][0]
    distance = band_totlength[i]/30.0 # distance between line segments that do not coincide

    iband = 0
    xpos = 0.0
    labels = [ (0.0,band_segments[i][0][4]) ]

    for start,end,length,npoint,startname,endname in band_segments[i]:
	iband += 1

        if any(start != prev_end):
            xpos += distance
            labels += [ (xpos,startname) ]

        xvals = xpos+linspace(0,length,npoint)
        xpos = xvals[-1]

        labels += [ (xpos,endname) ]

        prev_end = end
        prev_endname = endname

        for spin in range(1,max_spin_channel[i]+1):
            fname = directory[i] + "/band%i%03i.out"%(spin,iband)
            idx = []
            kvec = []
            band_energies = []
            band_occupations = []
            for line in file(fname):
                words = line.split()
                idx += [ int(words[0]) ]
                kvec += [ map(float,words[1:4]) ]
                band_occupations += [ map(float,words[4::2]) ]
                band_energies += [ map(float,words[5::2]) ]
       		# Apply energy offset (if specified) to all band energies
   		band_energies[-1] = [x - energy_offset[i] for x in band_energies[-1]]
	    maxLength = 0
	    for x in band_energies:
                if maxLength < len(x):
	            maxLength = len(x)
            for x in band_energies:
		for y in xrange(maxLength-len(x)):
		    x.append(x[-1])
            assert(npoint) == len(idx)
            band_energies = asarray(band_energies)
            # Now perform spline interpolation on band structure, if requested
	    if should_spline == True:
                xvals_smooth = np.linspace(xvals.min(),xvals.max(),spline_factor*len(xvals) ) # Interpolated x axis for spline smoothing
		new_band_energies = []
	        for b in range(band_energies.shape[1]): # Spline every band, one by one
		    new_band_energies.append(spline(xvals, band_energies[:,b], xvals_smooth))
                band_energies = asarray(new_band_energies).transpose() # recombine the bands back into the original data format
                xvals_plot = xvals_smooth # and use the interpolated x axis for plotting
            else:
                xvals_plot = xvals # Use the "raw" x axis for plotting

            for b in range(band_energies.shape[1]):
		# Make sure legend entry only appears once
		if b == 0 and iband == 1:
			if max_spin_channel[i] > 1 and spin == 1:
                		ax_bands.plot(xvals_plot,band_energies[:,b],color=color_array[i],label=plotName[i] + " up")
			elif max_spin_channel[i] > 1 and spin == 2:
                		ax_bands.plot(xvals_plot,band_energies[:,b],color=color_array[i],linestyle="--",label=plotName[i] + " dn")
			else:
                		ax_bands.plot(xvals_plot,band_energies[:,b],color=color_array[i],label=plotName[i])
		else:
			if max_spin_channel[i] > 1 and spin == 2:
                		ax_bands.plot(xvals_plot,band_energies[:,b],color=color_array[i],linestyle="--")
			else:
                		ax_bands.plot(xvals_plot,band_energies[:,b],color=color_array[i])

    tickx = []
    tickl = []
    for xpos,l in labels:
        ax_bands.axvline(xpos,color='k',linestyle=':')
        tickx += [ xpos ]
        if len(l)>1:
            l = "$\\"+l+"$"
        tickl += [ l ]
    for x, l in zip(tickx, tickl):
        print "| %8.3f %s" % (x, repr(l))

    ax_bands.set_xlim(labels[0][0],labels[-1][0])
    ax_bands.set_xticks(tickx)
    ax_bands.set_xticklabels(tickl)

#######################

if CUSTOM_YLIM:
    ax_bands.set_ylim(ylim_lower,ylim_upper)
else:
    ax_bands.set_ylim(-20,20) # just some random default -- definitely better than the full range including core bands
legend = ax_bands.legend()
if legend:
    frame = legend.get_frame()
    frame.set_linewidth(default_line_width)

#######################

matplotlib.rcParams['savefig.dpi'] =  print_resolution
matplotlib.rcParams['font.size'] = font_size

print
print "The resolution for saving figures is set to ", matplotlib.rcParams['savefig.dpi'], " dpi." 

if should_spline:
	print
	print "Spine interpolation has been used on the band structure, with an interpolation factor of ", spline_factor
	print "You should check this band structure against the un-interpolated version, as splining may cause some small artifacts not present in the original band structure." 

def on_q_exit(event):
    if event.key == "q": sys.exit(0)
connect('key_press_event', on_q_exit)

show()
