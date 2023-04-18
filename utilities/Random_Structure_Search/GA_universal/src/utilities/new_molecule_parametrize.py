'''
Created on Sep 16, 2013

@author: supady
'''
import numpy
from operator import itemgetter

import openbabel
from pybel import *
import pybel


def parametrize(smile, custom):
    with open('input.dat', 'r') as inf:
        dict_from_file = eval(inf.read())
    
    mol = readstring("smi", smile)
    smart_torsion = Smarts(dict_from_file['smart_torsion'])  # Matches the 4 atom long chains with a rotatable bond inside = all possible torsions 
    smart_cistrans = Smarts(dict_from_file['smart_cistrans'])  # Customized smart for a candidate cis/trans (to be checked!)
    torsion = sorted(list(smart_torsion.findall(mol)), key=itemgetter(1))
    cistrans = sorted(list(smart_cistrans.findall(mol)), key=itemgetter(1))
    to_del_tor = []
    to_del_cistrans = []
    custom_torsion = []
    
    if custom == True:
        smart_custom = Smarts(dict_from_file['smart_custom'])  # for the CF3 
        custom = sorted(list(smart_custom.findall(mol)), key=itemgetter(1))
        to_del_custom = []

        def compare(list1, list2):
            ln = []
            for i in list1:
                if i not in list2:
                    ln.append(i)
            return ln
        
        custom_torsion = compare(torsion, custom)
        for x in reversed(range(len(custom_torsion))):  # maps repeated   to delete
            for y in reversed(range(x)):
                if (itemgetter(1)(custom_torsion[x]) == itemgetter(1)(custom_torsion[y]) and itemgetter(2)(custom_torsion[x]) == itemgetter(2)(custom_torsion[y])) or (itemgetter(1)(custom_torsion[x]) == itemgetter(2)(custom_torsion[y]) and itemgetter(2)(custom_torsion[x]) == itemgetter(1)(custom_torsion[y])):
                    to_del_custom.append(y)
        to_del_custom.reverse()
        for d in reversed(list(set(to_del_custom))):  # deletes the redundand cis/trans
            del custom_torsion[d]
    
    for x in reversed(range(len(cistrans))):  # maps repeated  cis/trans to delete
        for y in reversed(range(x)):
            if (itemgetter(1)(cistrans[x]) == itemgetter(1)(cistrans[y]) and itemgetter(2)(cistrans[x]) == itemgetter(2)(cistrans[y])) or (itemgetter(1)(cistrans[x]) == itemgetter(2)(cistrans[y]) and itemgetter(2)(cistrans[x]) == itemgetter(1)(cistrans[y])):
                to_del_cistrans.append(y)
    to_del_cistrans.reverse()
    for d in reversed(list(set(to_del_cistrans))):  # deletes the redundand cis/trans
        del cistrans[d]
        
    for x in reversed(range(len(torsion))):  # maps repeated torsions to delete
        for y in reversed(range(x)):
             if (itemgetter(1)(torsion[x]) == itemgetter(1)(torsion[y]) and itemgetter(2)(torsion[x]) == itemgetter(2)(torsion[y])) or (itemgetter(1)(torsion[x]) == itemgetter(2)(torsion[y]) and itemgetter(2)(torsion[x]) == itemgetter(1)(torsion[y])):
                to_del_tor.append(y)
    to_del_tor.reverse()
    for d in reversed(list(set(to_del_tor))):  # deletes the redundand torsions
        del torsion[d]
    
    
    mol.make3D()  # creates a 3D representation applying force field!
    atoms = mol.OBMol.NumAtoms()
    bonds = mol.OBMol.NumBonds()
    mol.write("sdf", "mol.sdf", overwrite=True)  # creates the sdf of the molecule of interest
    mol.write("ct", "mol.ct", overwrite=True)
    
    w = open('mol.ct')
    lines = w.readlines()
    con = numpy.array([[i for i in line.strip().split()] for line in lines[(atoms + 2):]])  # reads only the lines with indicies of connected atoms 
    w.close()
    connectivity = con[0:, 0:2]  # deletes the information about the kind of bond
    con1 = con[0:, 0:1] 
    con2 = con[0:, 1:2]
    matrix_of_connections = [] 
    for i in xrange(atoms):  # builds a atoms*atoms matrix with zeros
        matrix_of_connections.append([])
        for j in xrange(atoms):
            matrix_of_connections[i].append(0)
    x = 0
    while x < bonds:  # the atoms that are connected get 1 instead of 0
        k = int(con1[x][0]) - 1
        v = int(con2[x][0]) - 1
        matrix_of_connections[k][v] = 1
        matrix_of_connections[v][k] = 1
        x = x + 1
    
    
    
    return (atoms, bonds, torsion, cistrans, matrix_of_connections, custom_torsion)
    

# mol = readstring("smi", "N(=Nc1ccccc1)c2ccccc2") # molecule
# mol = readstring("smi", "C1CCCCC1") # cyklohexane
# mol = readstring("smi", "CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=O)C)C") # retinal
# mol = readstring("smi", "c1ccc(c(c1)N=Nc1ccccc1)NC(=S)Nc1cc(cc(c1)C(F)(F)F)C(F)(F)F") # molecule 8
# mol = readstring("smi", "CCC(=O)N[C@]([C@H](OC(C)(C)C)C)(CO)C(=O)OC") #dawmoe
# mol=readstring("smi", "CCCCCC") # hexane
# mol=readstring("smi", "CCCCCCCCCCCCCCCCCC") # C18H38
# mol = readstring("smi", "c1cccc(c1)N=Nc1ccc(cc1)NC(=S)Nc1cc(cc(c1)C(F)(F)F)C(F)(F)F") # molecule 10
# mol = readstring("smi", "c1ccc(cc1CN)N=Nc1c(cc(c(c1C)NC(=S)Nc1cc(cc(c1)C(F)(F)F)C(F)(F)F)C)C") # catalyst
# mol = readstring ("smi", "CCCC") # butan
# mol = readstring ("smi", "C") # metan
# mol= readstring("smi", "N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)O)C)C)C") # Ac-Ala-Ala-Ala-NMe
# mol = readstring ("smi", "c1ccccc1COC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NC") 
# mol=readstring("smi", "N[C@@H](Cn1oc(=O)[nH]c1=O)C(O)=O")
# mol=readstring("smi", "O=C(O)C(N)C") #alanine
# mol=readstring("smi", "c1cc(ccc1)NC(=S)Nc1ccccc1") #thiocore
