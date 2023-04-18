#!/usr/bin/python
# Reads an output file of FHI-aims, gets the dipole moments and calculates
# various time correlation functions
# MR - 2013 (latest version)

import numpy
from sys import argv, exit, stderr
from math import pi, cos, sqrt, exp
from os import system
from optparse import OptionParser

def error(msg):
   """ Writes error message and quits
   """
   stderr.write(msg + "\n")
   exit(3)



def gaussian(length, broad):
 """ Calculates a Gaussian centered at length/2 and with broad=2*sigma
 """
 x0 = int(length/2)
 return [1/sqrt(pi*broad**2)*exp(-((i-x0)/(broad))**2) for i in numpy.arange(0, length)]


def _main():

   usage = """ %prog control.autocorr.in
        Reads an output file of FHI-aims and calculates either an IR spectrum or
        VDOS, based on the Fourier transform of the dipole time derivative or 
        the velocity, respectively. 

        Please provide a file control.autocorr.in with the following specifications (without the comment lines!!!) 
        (the numbers are just examples, choose them according to your particular case)
        choice dipole # this field takes either 'dipole' or 'velocity'
        sampling_interval 1 # interval between MD time steps that you consider
        broadening 2. # sets the 'broad' parameter (2*sigma), in cm-1, for the convolution of the spectra with a Gaussian
        cut_traj 20000 # cuts your trajectory after a certain number of MD steps. -1 corresponds to not cutting anything
        n_beads 16 # number of beads, if you are analyzing a PIMD run. 
        prefix pimd_h5o2 # name of the output files, without numbered extension. Both for PI runs or for averages all outputs should have the same starting name, followed by "_00", "_01", etc. 
        time_step 0.0005 # time step used in the simulaiton, in units of ps
        """
   parser = OptionParser(usage=usage)
#   parser.add_option("-e","--exp",
#                     default="exp.dat",dest="exp",metavar="filename",
#                     help="input file experimental spectrum [default: %default]")

   options, args = parser.parse_args()

   if (len(args) != 1):
     error("Please provide control.autocorr.in and execute: python auto-correlate-PI.py control.autocorr.in")
   input_file=open(args[0])

   for line in input_file:
      if "choice" in line:
          choice=line.split()[-1]
      if "sampling_interval" in line:
          delta_t = int(line.split()[-1])
      if "broadening" in line:
          broadening = float(line.split()[-1])
      if "cut_traj" in line:
          cut = int(line.split()[-1])
      if "n_beads" in line:
          nbeads=int(line.split()[-1])
      if "prefix" in line:
          pref=str(line.split()[-1])
      if "time_step" in line:
          print "Attention, you set the time step explicitly in the control.autocorr.in file"
          MD_time_step=float(line.split()[-1])    
# Check if all input flags exist:
   try:
     choice
     delta_t
     broadening
     cut 
     nbeads
     pref 
   except NameError:
    print 'There is a flag missing in your input file. Please correct and try again'
    exit()

#build list of PI files
   list_of_files = [(pref+'_%02d' % i) for i in range(nbeads)]


   print " "
   print "Now wait for the files autocorr.dat raw_fourier_transform.dat and convoluted_fourier_transform.dat to be generated"

# find all dipole moments and keep them in a list


   if choice == 'dipole':
           dipoles = [[] for i in range(len(list_of_files))]
           for i,data in enumerate(list_of_files):
               file=open(data)
               for line in file:
                  if "Molecular dynamics time step" in line:
                     MD_time_step=float(line.split()[-2]) # always in pico-seconds -> 1.10^{-12} s
                  if "Total dipole moment" in line:
                     dipoles[i].append(map(float, line.split()[-3:]))
               file.close()
	# if no dipoles were found, exit the program
               if not len(dipoles[i]):
                   print "No dipoles found in file!!!!"
                   exit()
           final_data=dipoles
   elif choice == 'velocity':
        velocities = [[] for i in range(len(list_of_files))]
        for i,data in enumerate(list_of_files):
            file=open(data)
            for line in file:
		if "Number of atoms" in line:
			n_atoms = int(line.split()[-1])
		if "Molecular dynamics time step" in line:
			MD_time_step=float(line.split()[-2]) # always in pico-seconds -> 1.10^{-12} s
		if "velocity" in line:
			velocities[i].append(map(float, line.split()[-3:]))
            file.close()
	# if no velocities were found, exit the program
       	    if not len(velocities):
         	print "No velocities found in file!!!!"
	        exit()

	#rearrange velocities array so that all velocities from the same time step are in a single line
	velocities_arranged=[[] for i in range(len(list_of_files))]
        for inst, data in enumerate(velocities):
	  for i in range(len(data)/n_atoms):
		temp=[0]
		for j in range(i*n_atoms, i*n_atoms+n_atoms, 1):
			temp=temp+data[j]
		velocities_arranged[inst].append(temp[1:])
        final_data=velocities_arranged
   else:
	print "We don't compute this quantity. Bye bye!"
	exit()




# now define t0's and calculate the correlation for all of them, storing in another list.
# for now I will use all of them, can change this later...

   dt = 0 #time steps of the MD run
   file = open("autocorr.dat", "w")
# calculate the estimator:
   quantity_array=numpy.average(final_data, axis=0)[:cut]
   n = len(quantity_array)

# take derivatives for dipoles, just normal quantities for velocities
   if choice == 'dipole':
      quantity_array_deriv=[(quantity_array[i+1] - quantity_array[i]) for i in range(n-1)]
      avg = numpy.average(quantity_array_deriv, axis=0)
      n=n-1
      quantity_array=numpy.array(quantity_array_deriv)
   if choice=='velocity':
      avg = numpy.average(quantity_array, axis=0) 

# now build the autocorrelation
   time_correlations=numpy.array([numpy.sum(numpy.mean(((quantity_array[:(n)-i:delta_t]- avg )*(quantity_array[i:(n):delta_t]- avg)), 0)) for i in range(int(n/2.))])
# Now we window the autocorrelation funciton and pad with zeroes
# window the function using the Hann window as w0 (maximum at 0)
#window_time_correlations = [time_corr_av[i]*0.5*(1+cos(2.*pi*i/(2*cutoff-1.))) for i in range(len(time_corr_av))] 
# or a triangular window
   window_time_correlations = [time_correlations[i]*(1.-float(i)/len(time_correlations)) for i in range(len(time_correlations))] 
# pad with 10000 zeroes
   window_time_correlations=numpy.append(window_time_correlations, numpy.zeros(10000))
   time_correlations=numpy.append(time_correlations, numpy.zeros(10000))
   n_total=len(window_time_correlations)
#write autocorrelation function
   temp = [file.write("%.10f %.12f %.12f\n" % (j*MD_time_step, i, z)) for i,j,z in zip(time_correlations, range(len(time_correlations)), window_time_correlations)]
   file.close()

#perform the fft
#call fortran program that does it

   if (system('./home_made_ft.x') is not 0):
      print "Please put the binary home_made_ft.x into this folder"
      exit

#normalize by total simulation time and broaden
   file=open('raw_fourier_transform.dat')
   ft = [map(float, i.split()[:]) for i in file.readlines()[:]]
# IS THIS RIGHT? CHECK!
   ft_norm=[i[1]*len(ft)/(n_total) for i in ft]
   file.close()

   f2=open('norm_raw_fourier_transf.dat', "w")
   temp = [f2.write("%.3f %.10f\n" % (i[0], j)) for i,j in zip(ft,ft_norm)]
   f2.close()

#get only second column for broadening
   ft1=[i[1] for i in ft]
   my_gaussian= gaussian(100, broadening)
# convolute with gaussian
   convoluted=numpy.convolve(ft1, my_gaussian, 'same')
   file = open("convoluted_fourier_transform.dat", "w")
   temp = [file.write("%.6f %.12f\n" % (i[0], j)) for i,j in zip(ft,convoluted)]

   file.close()

   print "It is done, happy analysis!"

if __name__ == "__main__":
      _main()

