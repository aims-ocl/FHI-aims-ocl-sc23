#!/usr/bin/python 
import sys

if len(sys.argv) < 3:
	print "WIEN2k_Band_to_aims.py takes in WIEN2k output file from either x lapw1 (-up/dn) -band or x lapwso"
        print "-band (respectively, <case>.output1, <case>.output1up/dn, or <case>.outputso) and converts them into"
        print "an FHI-aims file format suitable for plotting with aimsplot_compare.py.  The syntax is "
        print "                          WIEN2k_Band_to_aims.py <file> <spin> <Fermi>"
	print "where:"
	print "<file> is the associated WIEN2k file,"
	print "<spin> denotes the spin associated with the file (0 for non-spin-polarized OR SOC regardless of"
	print "initial print spin polarization, 1 for up, and 2 for dn), and" 
	print "<Fermi> is the Fermi level used for occupations."
	print "Note that the Fermi level is assumed to be in units of Rydbergs, as is outputted by WIEN2k"
	print "in the <case>.output2 file."
	print "Also note that no smearing is applied;  levels below will be fully occupied and levels above fully"
	print "unoccupied."
	print "The output will be stored into various band* files and the skeleton of a control.in needed for "
	print "plotting with aimsplot_compare.py will be output to stdout.  To use aimsplot_compare.py, you will"
	print "need to supply your own geometry.in skeleton containing, at bare minimum, lattice vectors."
	print "NOTE:  THIS IS VERY EXPERIMENTAL CODE.  MAKE SURE YOU CHECK THE OUTPUT."
	print "ALSO, BECAUSE WIEN2k OUTPUTS AN ENERGY WINDOW OF THE BAND STRUCTURE, AS"
	print "OPPOSED TO A SET NUMBER OF BANDS AS IS ASSUMED BY aimsplot.py, YOU "
	print "WILL NEED TO USE aimsplot_compare.py TO VIEW THE RESULTING BAND STRUCTURE."
	exit()

theFile = open(sys.argv[1],'r').readlines()
spin = int(sys.argv[2])
if len(sys.argv) < 4:
	fermi = "NA" # Will also set all occupations to 0
else:
	fermi = float(sys.argv[3])

counter = 0
i_point = []
RY_TO_eV = 27.2113845/2.
isStartAndStop = False
isSpinOrbitCoupled = False
nBand = 0
k_points = []
energies = []
bandStartName = []
bandEndName = []
# Look for the gamma point
while True:
	i = theFile[counter].split()
	pointName = ""
	if len(i) > 2 and i[1] == "spin-orbit" and i[2] == "calculation":
		isSpinOrbitCoupled = True
	if len(i) > 3 and i[0] == 'K=':
		# Possible k-point, save coordinate information just in case
		# This check is necessary because beginning of file contains a lot of "K=" lines that are not relevant
		k_x = i[1]
		k_y = i[2]
		k_z = i[3]
		# If this k-point is labeled, it is presumed either the beginning or ending of the current band
		if len(i) > 4:
			pointName = i[4]
		else:
			pointName = ""

		# Jump ahead two lines;  if this is a k-point, then this line will be "EIGENVALUES ARE:", so we
		# can begin to parse it.  Otherwise, ignore this fake k-point and move on
		counter += 2
		i = theFile[counter].split()
		if i[0] == "EIGENVALUES":
			# If k-point was labeled, begin a new band and end the current band
			# The only exceptions are when this is the very first k-point, and 
			# when the next k-point is also labeled, which implies the k-path is disconnected
			if pointName != "":
				if nBand == 0:
					nBand += 1
					i_point.append(0)
					k_points.append([])
					energies.append([])
					bandStartName.append(pointName)
				elif i_point[-1] == 1:
					i_point[-1] = 0
					k_points[-1] = []
					energies[-1] = []
					bandStartName.pop()
					bandStartName.append(pointName)
				else:	
					isStartAndStop = True
					nBand += 1
					i_point.append(0)
					k_points.append([])
					energies.append([])
					bandEndName.append(pointName)
					bandStartName.append(pointName)
			i_point[-1] += 1
			k_points[-1].append([k_x, k_y, k_z])
			energies[-1].append([])
			counter += 1
			if isStartAndStop:
				i_point[-2] += 1
				k_points[-2].append([k_x, k_y, k_z])
				energies[-2].append([])
			# Now parse all eigenvalues for this k-point
			while True:
				i = theFile[counter].split()
				# If an empty line is encountered, ignore
				if len(i) == 0:
					counter += 1
					continue
				# Finished parsing eigenvalues
				if len(i) > 1 and i[1] == "EIGENVALUES":
					counter += 1
					break
				# We set the given Fermi level to be the zero of energy here
				for j in i:
					if fermi != "NA":
						energies[-1][-1].append((float(j)-fermi)*RY_TO_eV)
						if isStartAndStop:
							energies[-2][-1].append((float(j)-fermi)*RY_TO_eV)
					else:
						energies[-1][-1].append((float(j))*RY_TO_eV)
						if isStartAndStop:
							energies[-2][-1].append((float(j))*RY_TO_eV)
				counter += 1
			isStartAndStop = False				
	counter += 1
	if counter == len(theFile):
		# Throw away the last band if it's empty (i.e. contains only the last labeled point)
		if i_point[-1] == 1:
			nBand -= 1

		if spin == 0:
			print "spin none"
		else:
			print "spin collinear"
		if isSpinOrbitCoupled:
			print "calculate_perturbative_soc"
		for i in xrange(nBand):
			if spin == 0:
				fileName = "band%i%03i.out"%(1,i+1)
			else:
				fileName = "band%i%03i.out"%(spin,i+1)
            		bandFile = open(fileName,"w")
			for j in xrange(i_point[i]):
				bandFile.write(str(j+1) + " ")
				bandFile.write(k_points[i][j][0] + " " + k_points[i][j][1] + " " + k_points[i][j][2] + " ")
				for k in energies[i][j]:
					if fermi != "NA" and float(k) <= 0.:
						bandFile.write(" 1.0000 " + str(k) )
					else:
						bandFile.write(" 0.0000 " + str(k) )
				bandFile.write("\n")
			bandFile.close()
			print "output band ", k_points[i][0][0], k_points[i][0][1], k_points[i][0][2], k_points[i][-1][0], k_points[i][-1][1], k_points[i][-1][2], i_point[i], bandStartName[i][0], bandEndName[i][0], "#", bandStartName[i], bandEndName[i]
		exit()
