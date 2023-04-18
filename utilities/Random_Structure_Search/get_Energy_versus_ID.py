#!/usr/bin/python
__author__ = 'jkloppenburg and blange'

from optparse import OptionParser
import os


class randomSearchData:

	def __init__(self):
		self.jobDirs = []
		self.converged =[]
		self.energies = []
		
  
	def getJobDirs(self):
		path = os.getcwd()
		
		for jobDir in os.listdir(path):
			if (os.path.isdir(jobDir)):
				self.jobDirs += [jobDir]
    
	def getConverged(self):
		
		for dirs in self.jobDirs:
				converged = False
				os.chdir(dirs)
				for line in file("result.dat"):
					string = "Present geometry is converged."
					if (line.find(string) > 0):
						converged = True
						
				if (not converged):
					print ("%s is not converged" % (str(dirs)))
				self.converged += [converged]
				os.chdir("..")
				
	def getEnergies(self):
		for iJob in range(len(self.jobDirs)):
			if (self.converged[iJob]):
				os.chdir(self.jobDirs[iJob])
				for line in file("result.dat"):
					string = "Total energy of the DFT / Hartree-Fock s.c.f. calculation "
					if (line.find(string) > 0):
						words = line.split()
						energy = float(words[11])
						self.energies += [energy]
				os.chdir("..")
			else:
				self.energies += [0.0]
				
	def printEnergies(self):
		f = open("EnergyVsJob.dat","w")
		f.write ("# Jobnumber total Energy [eV]\n")
		lowestEnergy = min(self.energies)
		for iJob in range(len(self.jobDirs)):
			f.write ("%5s %15.8f\n" % (self.jobDirs[iJob],self.energies[iJob]-lowestEnergy))
			
		f.close()
			
							
				
		
		
# main
if __name__ =='__main__':

	parser = OptionParser()
	parser.description = "This tool collects data from a random search run."    
    
	#parser.add_option("-i", dest="file", help="input file", default = "geometry.in")
    
    
            
	(options, args) = parser.parse_args()        
 
	data = randomSearchData()
	data.getJobDirs ()
	data.getConverged()
	data.getEnergies()
	data.printEnergies()
	

	


