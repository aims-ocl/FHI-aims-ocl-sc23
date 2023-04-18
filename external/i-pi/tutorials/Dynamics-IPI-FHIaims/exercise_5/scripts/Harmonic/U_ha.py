import numpy as np
import sys
import os

"""
From a list of frequencies, We compute the  quantum harmonic total energy (Kinetic + Potential) "
Syntax:

   python U_ha.py <Input-Freq_file> <Max-Temperature>  <number of zero modes>


"""

def coth(x):
    """
    coth(x)
    Uses numpy's tanh(x).
    """
    return 1./np.tanh(x)


#np.set_printoptions(precision=3, suppress=True)
##Constants
h    = 6.62607e-34 #Jxs
Kb   = 1.38065e-23 #J/K
c    = 3.0e10      #cm/s
c    = 29979245800.0 #cm/s
R    = 8.314       #J/(kmol)
j2eV = 6.2415091e+18 #

#Open and read
name_input = sys.argv[1]
temp       = float(sys.argv[2])
zero       = 6 

freq = np.zeros(1500)
n = 0
if (os.path.exists(name_input)):
    filee=open(name_input,'r')
    print '# We are reading %s' %name_input
    print '# We consider that the freq are in cm^-1'
    print '# In a two column format file'
    line=filee.readline().split() #comment line

    while True:
        line=filee.readline().split()
        if not line:
             break
        freq[n]=float(line[1])
        n+=1
    filee.close()
else:
   print '# We can not open %s' %name_input
   sys.exit()

#Compute

print '# We have %i frequencies' %n
print '# We are discarding %i of them' %zero


ZPE=0.0
for i in range(zero,n):
    ff   = freq[i]*c*h #freq in J
    ZPE  += (  ff/2.0 )*j2eV
print '# ZPE   %f ev' %ZPE

#Compute U=K+V  for each temperature:

print '# T(K)  U_c(eV)  U_q(eV)' 


for t in range(50,int(temp+1),25):
    beta = 1.0/(float(t)*Kb) #J^-1 
    K    = 0.0
    V    = 0.0
    U_cl = 0.0
    for i in range(zero,n):
      ff    = freq[i]*c*h #freq in J

      f     = beta*ff
      K    += (ff/4.0)*coth(f/2.0)
      V    += (ff/4.0)*coth(f/2.0)
      U_cl += Kb*float(t)
    print '%i %f %f' %(t,U_cl*j2eV, (K+V)*j2eV)

