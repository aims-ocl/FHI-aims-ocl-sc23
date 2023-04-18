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

   usage = """ %prog -l list control.autocorr.in
        Reads an output file of FHI-aims and calculates either an IR spectrum or
        VDOS, based on the Fourier transform of the dipole time derivative or 
        the velocity, respectively. 

        Please provide a file control.autocorr.in with the following specifications (without the comment lines!!!) 
        (the numbers are just examples, choose them according to your particular case):
        choice dipole # this field takes either 'dipole' or 'velocity'
        sampling_interval 1 # interval between MD time steps that you consider
        broadening 2. # sets the 'broad' parameter (2*sigma), in cm-1, for the convolution of the spectra with a Gaussian
        cut_traj 20000 # cuts your trajectory after a certain number of MD steps. -1 corresponds to not cutting anything
        time_step 0.0005 # (optional, if it is already written in the FHI-aims output file) time step used in the simulaiton, in units of ps

        Please provide also a file containing the path of each file of the trajectories to be averaged
        using the -l option. Specifying more than one file means that the trajectories will be averaged! Example:
        file1.dat
        file2.dat
        """
   parser = OptionParser(usage=usage)
   parser.add_option("-l","--list",
                     default="list.dat",dest="list",metavar="filename",
                     help="File containing the path of each file of the trajectories to be averaged [default: %default]")

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
      if "time_step" in line:
          print "Attention, you set the time step explicitly in the control.autocorr.in file"
          MD_time_step=float(line.split()[-1])    
# Check if all input flags exist:
   try:
     choice
     delta_t
     broadening
     cut 
   except NameError:
     print 'There is a flag missing in your input file!! Please correct and try again'
     exit()

#build list of runs 
   list_file=open(options.list)
   list_of_files = [str(i)[:-1] for i in list_file]
   list_file.close()


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
   time_corr_array=[]

   for run in final_data:
      quantity_array = numpy.array(run[:cut])
      n=len(quantity_array)
      if choice == 'dipole':
         quantity_array_deriv=[(quantity_array[i+1] - quantity_array[i]) for i in range(n-1)]
         avg = numpy.average(quantity_array_deriv, axis=0)
         n=n-1
         quantity_array=numpy.array(quantity_array_deriv)
      if choice=='velocity':
         avg = numpy.average(quantity_array, axis=0)
      time_correlations=numpy.array([numpy.sum(numpy.mean((quantity_array[:n-i:delta_t]- avg )*(quantity_array[i:n:delta_t]- avg), 0)) for i in range(int(n/2.))])
      #normalize
      #time_correlations = time_correlations/time_correlations[0]
      time_corr_array.append(time_correlations)

# Now we cut for the shortest run and average
   dimensions = numpy.array([len(i) for i in time_corr_array])
   time_corr_cut = numpy.array([j[:min(dimensions)] for j in time_corr_array])
   time_corr_av=numpy.average(time_corr_cut, axis=0)


# window the function using the Hann window as w0 (maximum at 0)
#window_time_correlations = [time_corr_av[i]*0.5*(1+cos(2.*pi*i/(2*cutoff-1.))) for i in range(len(time_corr_av))] 
# or a triangular window
   window_time_correlations = [time_corr_av[i]*(1.-float(i)/len(time_corr_av)) for i in range(len(time_corr_av))]
# Pad with zeroes
   window_time_correlations=numpy.append(window_time_correlations, numpy.zeros(10000))
   time_corr_av=numpy.append(time_corr_av, numpy.zeros(10000))
   n_total=len(time_corr_av)
#write autocorrelation function
   temp = [file.write("%.10f %.12f %.12f\n" % (j*MD_time_step, i, z)) for i,j,z in zip(time_corr_av, range(len(time_corr_av)), window_time_correlations)]
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
   temp = [f2.write("%.4f %.12f\n" % (i[0], j)) for i,j in zip(ft,ft_norm)]
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

