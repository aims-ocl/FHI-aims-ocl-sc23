#V1.0
import numpy as np
import sys
import os
import math

#Syntax
# python proc.py <name_input>  <mode>  


#Modes
# T = temperature
# V = Potential energy
# K = Kinetic energy
# K+V = K+V = U = Internal energy

#We consider this output format###
# column   1     --> step : The current simulation time step.
# column   2     --> time{picosecond} : The elapsed simulation time.
# column   3     --> temperature{kelvin} : The current temperature, as obtained from the MD kinetic energy.
# column   4     --> potential{electronvolt} : The physical system potential energy.
# column   5     --> kinetic_md{electronvolt} : The kinetic energy of the (extended) classical system.
###  

#Hard input
ndim =100000
np.set_printoptions(precision=5, suppress=True)
name_input  = sys.argv[1]
mode = sys.argv[2]
var=np.array([0])
     
#Open files
if (os.path.exists(name_input)):
      print 'We are reading %s' %name_input 
      f=open(name_input,'r')
else:
      print 'We can not open the input file'
      sys.exit()  

#Skip comments line
i=0
while True:
     line = f.readline()
     
     if '#' not in line:
         break
     i+=1

f.seek(0)
for ii in range(i):
    f.readline()

#Read
ncount = 0
while True:
   line = f.readline().split() 
   if not line:
       break
   ncount+=1
   if mode == 'T':
      var= np.append(var,float(line[2]))
   elif mode == 'V':
      var= np.append(var,float(line[3]))
   elif mode == 'K':
      var= np.append(var,float(line[4]))
   elif mode == 'K+V':
      var= np.append(var, ( float(line[3]) +float(line[4]) ) )
   else:
      print 'We can not recognize the mode'
      sys.exit()

#Compute average and Standar desviation:
mean = np.mean(var)
std  = np.std(var)
print '   '
print '   '
print 'We are done for now '
print 'We have %d data points' %ncount
print 'The selected mode is %s' %mode
print 'The mean value is        %f' %mean
print 'The standar deviation is %f' %std
print '   '
print '   '
print ' %f  %f' %(mean,  std)

