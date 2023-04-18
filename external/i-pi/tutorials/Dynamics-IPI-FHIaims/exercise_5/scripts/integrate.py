import numpy as np
import sys
import os
from scipy.integrate import simps as simpson
from scipy.integrate import trapz as trap

#np.set_printoptions(precision=3, suppress=True)
"""
Syntax: python integrate.py <name_input> <number of points> 
"""

#Open and read
name_input  = sys.argv[1]
m           = int(sys.argv[2])

#Hard input
##Constants
#t    = 1.0/315774.66 #convert K to hartree
#h    = 6.62607e-34 #Jxs
#Kb   = 1.38065e-23 #J/K
#c    = 3.0e10      #cm/s
#c    = 29979245800.0 #cm/s
#R    = 8.314       #J/(kmol)
#j2eV = 6.2415091e+18 #
#

T          = np.zeros(m)
U          = np.zeros(m)   #This is K+V
DF         = np.zeros(m)

#Read
if (os.path.exists(name_input)):
    filee=open(name_input,'r')
    print '#We are reading %s' %name_input
    i=0
    icontrol=0
    while i <  m:
        icontrol+= 1
        line=filee.readline()
        if not line or '#' in line:
          if icontrol >25:
            print 'We are not able to find %i points' %m
            sys.exit()
          continue
        line=line.split()
        T[i]      = line[0]
        U[i]      = line[1]
        i+=1

    filee.close()
else:
    print '#We can not find  %s' %name_input
    sys.exit()
#Compute
x=T[0]
y=U[0]/(T[0]**2)

print "#  "
print "#Integration starting at %f K" %T[0]
print "#First column temperature"
print "#Second column anharmonic contribution to the Free Energy"
print "#  "

print T[0],0.0 
for i in range(1,m):
    x=np.append(x,T[i])
    y=np.append(y,U[i]/(T[i])**2)

    DF[i]=simpson(y,x)
    print T[i],-DF[i]*T[i]

