from sys import argv
from numpy import array, sqrt, histogram

name=argv[1]
inp_file = open(name)

def distance(a,b):
    """ Defines euclidean distances between a and b 
    """
    return sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)

oo_out=open('rdf_oo_'+name+'.dat', 'w')
oh_out=open('rdf_oh_'+name+'.dat', 'w')
hh_out=open('rdf_hh_'+name+'.dat', 'w')

all_oo=[]
all_oh=[]
all_hh=[]

while True:
    line = inp_file.readline()
    if not line:
       break
    if "Number of atoms" in line:
        natoms=int(line.split()[-1]) # this should always happen before next structure, so this should be fine...
    if "Time step" in line:
        tstep=float(line.split()[-1])
    if "Atomic structure that was used" in line:
        counter=0
        allp=[]
        alltype=[]
        while counter<natoms:
           line=inp_file.readline()
           if not line:
              break
           if "atom" in line:
              type=str(line.split()[-1])
              pos=map(float, line.split()[1:4])
              counter=counter+1
              allp.append(pos)
              alltype.append(type)
        for species, position in zip(alltype, allp):
             for species2, position2 in zip(alltype, allp):
                dist=distance(position, position2)
                if dist>0:
                   if species=='O' and species2=='O':
                      all_oo.append(dist)
                   if species=='H' and species2=='H':
                      all_hh.append(dist)
                   if species=='O' and species2=='H':
                      all_oh.append(dist)


histog_oo, bin_oo=histogram(array(all_oo), bins=500, range=(0, 5))
histog_hh, bin_hh=histogram(array(all_hh), bins=500, range=(0, 5))
histog_oh, bin_oh=histogram(array(all_oh), bins=500, range=(0, 5))


for i,j in zip(histog_oo, bin_oo[1:]):
    oo_out.write('%f\t%f\n' % (j, float(i)/len(all_oo)) )
for i,j in zip(histog_oh, bin_oh[1:]):
    oh_out.write('%f\t%f\n' % (j, float(i)/len(all_oh)) )
for i,j in zip(histog_hh, bin_hh[1:]):
    hh_out.write('%f\t%f\n' % (j, float(i)/len(all_hh)) )         

oo_out.close()
oh_out.close()
hh_out.close()
