#!/usr/bin/python
from numpy import *
from matplotlib.pyplot import *

#--------------bechmark compared with origin--------------------
#X & (1/2,0,0)  &  &  &\\
#Y & (0,1/2,0)  &  &  &\\
#Z & (0,0,1/2)  &  &  &\\
# m$^{*}$ /m$_{e}$             & $\Gamma$-X&    $\Gamma$-Y   &  $\Gamma$-Z  \\ 
#(3 point near $\Gamma$)       &       0.33     &      4.0   &  --- &\\
#(by Origin, 2nd diff method)  &       0.35     &      3.8   &     ---&  \\
#(by Origin, fitting 10 points method)& 0.35    &      3.8   &  --- &   \\
#([aims-effect-mass.py],fitting 10 points) & 0.34 &    3.8   &  ---& ---\\
#--------------end benchmark compared with origin---------------


print "==================================================================================="
print "Finding effect mass for FHI-aims! " 
print "need files: output  geometry.in  band100*.out "
print "                   ---shanghui,2014.09@Berlin "
print 

print "Which band do you want to use? ---> edit in this script, this is the easiest way at present..."
print "in geneal we see three direction:" 
print " band1001.out === Gamma-X "
print " band1002.out === Gamma-Y "
print " band1003.out === Gamma-Z "
print ""
print "We suppose you have 100 points along each path, and we pick up 10 points near gamma point"
n=10       # here we choose 10 points 
n_bands=3 # here we chosse 3 bands
print "*********************************begin warning******************************************"
print "For     mass~3,   10  points in [0,0.05]   is enough,the out_plot_band.f90 can be F15.5"
print "but for mass~0.2, 10  points in [0,0.0005] is needed,the out_plot_band.f90 should be F18.8"
print "you need to check by seeing the fitting plot"
print "*********************************end warning******************************************"
print 
print "An example contron.in for effect mass is: "
print "********************************begin contron.in*****************************************"
print "output band  0.0  0.0  0.0   0.005000000000  0.0000000000   0.000000000   100    Gamma  M "
print "output band  0.0  0.0  0.0   0.003333333333  0.003333333333   0.0000000000  100 Gamma  K"
print "output band  0.0  0.0  0.0   0.0000000000  0.0000000000   0.500000000     100  Gamma  A"
print "********************************end contron.in********************************************"
print 
print "===================================================================================="
print

#=========================(1) read HOMO=============================================
print "------(1) Reading n_HOMO and n_LUMO from 'output' ------"

for line in open("output"):
 if "      0.00000       " in line:
   words = line.split()
   n_energy=int(words[0])
   break

n_HOMO = 4+2*(n_energy-1)-1
n_LUMO = 4+2*n_energy-1
print 'n_HOMO=',n_HOMO
print 'n_LUMO=',n_LUMO
print 

#=======================(2) read lattic===============================================
print "------(2) Reading real lattic vector from geomtry.in ------"
#this code comes from aimsplot.py:
latvec = []

for line in file("geometry.in"):
    line = line.split("#")[0]
    words = line.split()
    if len(words) == 0:
        continue
    if words[0] == "lattice_vector":
        if len(words) != 4:
            raise Exception("geometry.in: Syntax error in line '"+line+"'")
        latvec += [ array(map(float,words[1:4])) ]

if len(latvec) != 3:
    raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

latvec = asarray(latvec)

print "Lattice vectors: a (Ang)"
for i in range(3):
    print latvec[i,:]
print


a1=latvec[0,:]
a2=latvec[1,:]
a3=latvec[2,:]

b1=2*pi*cross(a2,a3)/dot(a1,cross(a2,a3))
b2=2*pi*cross(a3,a1)/dot(a2,cross(a3,a1))
b3=2*pi*cross(a1,a2)/dot(a3,cross(a1,a2))
print "Reciprocal lattice vectors: b=2pi/a (Ang -1)"
print b1
print b2
print b3
print 

#=======================(3) read band===============================================
print "------(3) Reading band data from band and print data------"

#--------inintal all the arrays we need------------------
mass_hole_array = zeros(n_bands)
mass_electron_array = zeros(n_bands)
x_array = zeros((n,n_bands))
y_array_hole = zeros((n,n_bands))
y_array_hole_fit = zeros((n,n_bands))
y_array_electron = zeros((n,n_bands))
y_array_electron_fit = zeros((n,n_bands))
#--------end inintal all the arrays we need------------------

def read_and_plot(file_name):
    data = loadtxt(file_name)
    x=zeros(n)
    #x=data[0:n,1]
    x_vec=data[0:n,1:4]
    for i in range(n):
        x[i]=linalg.norm(x_vec[i,0]*b1+x_vec[i,1]*b2+x_vec[i,2]*b3)*0.529177 # b2
    y_hole=data[0:n,n_HOMO]
    y_hole=y_hole/27.2113845  # change energy from eV to a.u.
    fit_results=polyfit(x, y_hole,2,full=True)
    coefficients=fit_results[0]
    log=fit_results[1:]
    print 'effect hole mass = ', 1.0/(2.0*coefficients[0]), 'with SSR=',log[0]
    mass_hole=fabs(1.0/(2.0*coefficients[0]))
    polynomial = poly1d(coefficients)
    y_hole_fit = polynomial(x)
    
    y_electron=data[0:n,n_LUMO]
    y_electron=y_electron/27.2113845  # change energy from eV to a.u.
    fit_results=polyfit(x, y_electron,2,full=True)
    coefficients=fit_results[0]
    log=fit_results[1:]
    print 'effect electron mass = ', 1.0/(2.0*coefficients[0]), 'with SSR=',log[0]
    mass_electron=fabs(1.0/(2.0*coefficients[0]))
    polynomial = poly1d(coefficients)
    y_electron_fit = polynomial(x)

    return mass_hole, mass_electron, x[:], y_hole[:], y_hole_fit[:] ,y_electron[:],y_electron_fit[:]
    


#---------main call for read_and_plot subroutine--------------------
for i_band in range(n_bands):
    file_name = 'band100'  + str(i_band+1) + '.out'
    mass_hole_array[i_band],mass_electron_array[i_band],x_array[:,i_band], y_array_hole[:,i_band], y_array_hole_fit[:,i_band], y_array_electron[:,i_band], y_array_electron_fit[:,i_band] = read_and_plot(file_name)
    
f, ((ax1_hole, ax2_hole, ax3_hole), (ax1_electron, ax2_electron, ax3_electron)) = subplots(2, 3, sharex='col')

#ax1_hole.set_title('Hole effect mass')
ax1_hole.plot(x_array[:,0], y_array_hole[:,0],label='hole-X')
ax2_hole.plot(x_array[:,1], y_array_hole[:,1],label='hole-Y')
ax3_hole.plot(x_array[:,2], y_array_hole[:,2],label='hole-Z')
ax1_hole.plot(x_array[:,0], y_array_hole_fit[:,0],label='fit-X')
ax2_hole.plot(x_array[:,1], y_array_hole_fit[:,1],label='fit-Y')
ax3_hole.plot(x_array[:,2], y_array_hole_fit[:,2],label='fit-Z')
ax1_hole.legend(loc=0)
ax2_hole.legend(loc=0)
ax3_hole.legend(loc=0)

#ax1_electron.set_title('Eletron effect mass')
ax1_electron.plot(x_array[:,0], y_array_electron[:,0],label='electron-X')
ax2_electron.plot(x_array[:,1], y_array_electron[:,1],label='electon-Y')
ax3_electron.plot(x_array[:,2], y_array_electron[:,2],label='electron-Z')
ax1_electron.plot(x_array[:,0], y_array_electron_fit[:,0],label='fit-X')
ax2_electron.plot(x_array[:,1], y_array_electron_fit[:,1],label='fit-Y')
ax3_electron.plot(x_array[:,2], y_array_electron_fit[:,2],label='fit-Z')
ax1_electron.legend(loc=0)
ax2_electron.legend(loc=0)
ax3_electron.legend(loc=0)

show()



#----------print summary-----------------
print 
print "------(4) Print summary data----------------------------------"
print "==============================================================="
print " direction           G-X                G-Y         G-Z  " 
print "---------------------------------------------------------------"
print "  hole            %8.2f       %8.2f      %8.2f" % (mass_hole_array[0],mass_hole_array[1],mass_hole_array[2]  )
print "  electron        %8.2f       %8.2f      %8.2f" % (mass_electron_array[0],mass_electron_array[1],mass_electron_array[2]  )
print "==============================================================="
