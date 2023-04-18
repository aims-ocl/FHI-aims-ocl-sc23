from sys import argv
import numpy as np

angstrom_m = 1e-10
ps_s=1.e-12
bohr_to_m = 5.2917721e-11
atomic_time = 2.4188843e-17
mass_el = 9.1093819e-31
kb=1.3806488e-23 #J/K
amu_kg=1.66053886e-27
mass_dict={'H':1.0079, 'He':4.0026, 'Li':6.941, 'Be':9.0122, 'B':10.811, 'C':12.0107, 'N':14.0067, 'O':15.9994, 'F':18.9984, 'Ne':20.1797, 'Na':22.9897, 'Mg':24.305, 'Al':26.9815, 'Si':28.0855, 'P':30.9738, 'S':32.065, 'Cl':35.453, 'K':39.0983, 'Ar':39.948, 'Ca':40.078, 'Sc':44.9559, 'Ti':47.867, 'V':50.9415, 'Cr':51.9961, 'Mn':54.938, 'Fe':55.845, 'Ni':58.6934, 'Co':58.9332, 'Cu':63.546, 'Zn':65.39, 'Ga':69.723, 'Ge':72.64, 'As':74.9216, 'Se':78.96, 'Br':79.904, 'Kr':83.8, 'Rb':85.4678, 'Sr':87.62, 'Y':88.9059, 'Zr':91.224, 'Nb':92.9064, 'Mo':95.94, 'Tc':98, 'Ru':101.07, 'Rh':102.905, 'Pd':106.42, 'Ag':107.868, 'Cd':112.411, 'In':114.818, 'Sn':118.71, 'Sb':121.76, 'I':126.904, 'Te':127.6, 'Xe':131.293, 'Cs':132.905, 'Ba':137.327, 'La':138.905, 'Ce':140.116, 'Pr':140.908, 'Nd':144.24, 'Pm':145, 'Sm':150.36, 'Eu':151.964, 'Gd':157.25, 'Tb':158.925, 'Dy':162.5, 'Ho':164.930, 'Er':167.259, 'Tm':168.934, 'Yb':173.04, 'Lu':174.967, 'Hf':178.49, 'Ta':180.948, 'W':183.84, 'Re':186.207, 'Os':190.23, 'Ir':192.217, 'Pt':195.078, 'Au':196.966, 'Hg':200.59, 'Tl':204.383, 'Pb':207.2, 'Bi':208.980, 'Po':209, 'At':210, 'Rn':222, 'Fr':223, 'Ra':226, 'Ac':227, 'Pa':231.036, 'Th':232.038, 'Np':237, 'U':238.029, 'Am':243, 'Pu':244, 'Cm':247, 'Bk':247, 'Cf':251, 'Es':252, 'Fm':257, 'Md':258, 'No':259, 'Rf':261, 'Lr':262, 'Db':262, 'Bh':264, 'Sg':266, 'Mt':268, 'Rg':272, 'Hs':277}

name=argv[1]
inp_file = open(name)
atomtype = [str(argv[2]), str(argv[3])]

def inst_temperature(type, p):
    #T=mass_dict[type]*amu_kg*(v[0]**2 + v[1]**2 + v[2]**2)*(angstrom_m/ps_s)**2/(3.*kb) # from velocities
    T=(p[0]**2 + p[1]**2 + p[2]**2)*(mass_el*bohr_to_m/atomic_time)**2/(3.*kb*mass_dict[type]*amu_kg) # from momenta
    return T

output=open('temperatures_'+name+'.dat', 'w')

output.write('# Time_step\tT(%s)\tT(%s)\n' % (atomtype[0], atomtype[1]))

for i in range(13):
  inp_file.readline()
while True:
    line = inp_file.readline()
    if not line:
       break
    tstep=line.split()[0]
    t1=[]
    t2=[]
    for i in xrange(6,11,3): # Oxygen atoms
       momenta=np.array(map(float, line.split()[i:i+3]))
       #velocities=momenta/float(mass_dict[atomtype[0]])
       #t1.append(inst_temperature(atomtype[0], velocities))
       t1.append(inst_temperature(atomtype[0], momenta))
    for i in xrange(12,26,3): # Hydrogen atoms
       momenta=np.array(map(float, line.split()[i:i+3]))
       #velocities=momenta/float(mass_dict[atomtype[1]])
       #t2.append(inst_temperature(atomtype[1], velocities))
       t2.append(inst_temperature(atomtype[1], momenta))
    T1=sum(np.array(t1))/len(t1)
    T2=sum(np.array(t2))/len(t2)
    output.write('%f\t%f\t%f\n' % (float(tstep), T1, T2))

output.close()
