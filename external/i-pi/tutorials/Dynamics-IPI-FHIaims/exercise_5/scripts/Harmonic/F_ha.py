import numpy as np
import sys
import os

np.set_printoptions(precision=3, suppress=True)
"""
From a list of frequencies, We compute the quantum harmonic free  energy  "
Syntax:

   python F_ha.py <Input-Freq_file> <Temperature>  

"""
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




#Compute F_c,F_q for each temperature:

#print '# T(K)  F_c(eV)  F_q(eV)' 

print '# T(K)   F_q(eV)' 


#Way one
#for t in range(30,int(temp+1)):
#   beta = 1.0/(float(t)*Kb) #J^-1 
#   Q_c  = 1.0
#   Q_q  = 1.0

#   for i in range(zero,n):
#      ff    = freq[i]*c*h #freq in J
#      f     = beta*ff
#      Q_c   = Q_c /f
#      Q_q   = Q_q* np.e**( (-f/2.0) ) / (1 - np.e**(-f) )
#   F_cc = -Kb*int(t)*np.log(Q_c)*j2eV
#   F_qq = -Kb*int(t)*np.log(Q_q)*j2eV #+extra
#
#   print '%i %f %f ' %(t,F_cc,F_qq )



################################################################
#Way two
for t in range(1,int(temp+1)):
   beta = 1.0/(float(t)*Kb) #J^-1 
   F_c = 0.0
   F_q = 0.0
   for i in range(zero,n):
      ff   = freq[i]*c*h #freq in J
      f    = beta*ff
      F_c += -Kb*int(t)*np.log(1.0/f)
      F_q += (  ff/2.0    + Kb*int(t)*np.log( 1 - np.e**(-f)  ) ) 
   F_c = F_c *j2eV
   F_q = F_q *j2eV

  # print '%i %f %f    ' %(t,F_c,F_q)
   print '%i  %f    ' %(t,F_q)

