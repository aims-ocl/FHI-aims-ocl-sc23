'''
Created on Sep 16, 2013

@author: supady
'''
import glob
import numpy
from random import randint, choice
import shutil

from core import user_input
from structures.structure import Structure
from utilities import check_geometry, conversions
from utilities.angle import *
from utilities.energy import *
from utilities.mmff94 import *
from utilities.new_molecule_parametrize import parametrize


def main(stoic, seed):
    config = user_input.get_config()    
    smile = config.get('random_gen' , 'smile')
    custom = config.get('random_gen' , 'custom')
    (atoms, bonds, torsion, cistrans, matrix_of_connections, custom_torsion) = parametrize(smile, custom)

    s = Structure()
    s.build_geo_whole_atom_format(aims_geo)
    return s

def create_population(popsize, cistrans, torsion, molc, atoms, bonds, directory, sourcedir, matrix_of_connections, energy_function, black_dir, machine):
    population = numpy.zeros((popsize, len(cistrans) + len(torsion)))
    random_tor_all = []
    random_cistrans_all = []
    relax_tor_all = []
    relax_cistrans_all = []
    global energy_all
    energy_all = []
    n = 0  # initializes the number of chromosomes
    cnt = 1  # intialized the number of trials 
    while n < popsize:
        print " "
        print "    New trial " + str(cnt)
        
        random_tor_n = []  # array for random torsions
        random_cistrans_n = []  # array for random 0/180 value of the cis/trans bond
        for x in range(len(cistrans)):
            # random_cistrans_n.append(choice([0,0]))
            random_cistrans_n.append(choice([180, 180]))
            # random_cistrans_n.append(choice([0,60,120,180,240,300]))

        for y in range(len(torsion)):
            # random_tor_n.append(randint(0, 179)) # for symmetric molecules
            random_tor_n.append(randint(0, 359) - 179)  # standard
            # random_tor_n.append(random.choice([0,60,120,180,240,300]))
        
        angle_set(molc, torsion, cistrans, random_tor_n, random_cistrans_n)
        conversions.sdf_to_xyz.WriteFile(molc, "molc_start_" + str(n) + ".xyz")
        check = check_geometry("molc_start_" + str(n) + ".xyz", atoms, bonds, matrix_of_connections)
        if check == True:
            
            if energy_function == "aims":
                cnt_black = len(glob.glob(black_dir + "/*.xyz"))
                shutil.copy("molc_start_" + str(n) + ".xyz", black_dir)
                old_name = os.path.join(black_dir, "molc_start_" + str(n) + ".xyz")
                new_name = os.path.join(black_dir, "black_" + str(cnt_black) + ".xyz")
                os.rename(old_name, new_name)
                
                energyc = energy_aims("molc_start_" + str(n) + ".xyz", directory, sourcedir, n, atoms, black_dir, machine)
            if energy_function == "mmff94":
                energyc = energy_mmff94(molc, "molc_start_" + str(n) + ".xyz", n)
            print "Energy of the candidate after the relaxation: " + str(energyc)
            print "Accepted. Chromosome " + str(n)
            
            angle_measure("molca_" + str(n) + ".sdf", torsion, cistrans)
            relax_tor_n = relax_tor()
            relax_cistrans_n = relax_cistrans()
                                
            energy_all.append(energyc)
            relax_tor_all.append(relax_tor_n)
            relax_cistrans_all.append(relax_cistrans_n)
                
            if len(cistrans) != 0  :
                print "Setting cis/trans isomerisation for chromosome " + str(n) + " (cis=0, trans=180)   " + str(random_cistrans_n)
                print "Cis/trans isomerisation for chromosome " + str(n) + " after the relaxation  " + str(relax_cistrans_n)
            else: 
                print "No cis/trans bond"
            if len(torsion) != 0  :
                print "Setting random torsion values of rotatable bonds for chromosome " + str(n) + str(random_tor_n)
                print "Torsion values of rotatable bonds for chromosome " + str(n) + " after the relaxation " + str(relax_tor_n)
            else: 
                print "No rotatable bond"
                
            n += 1
            cnt += 1
            
        else:
            print "Rejected. Sorry, the geometry cannot be accepted." 
            os.remove("molc_start_" + str(n) + ".xyz")
            cnt += 1
        
    ratio = float(float(n - 1) / float(cnt - 1))
    print "Success ratio: " + str(ratio)

    for j in range(popsize):
        population[j, 0:len(cistrans)] = relax_cistrans_all[j]
        population[j, len(cistrans):(len(cistrans) + len(torsion))] = relax_tor_all[j]

    return population
def energy_population():
    return energy_all
