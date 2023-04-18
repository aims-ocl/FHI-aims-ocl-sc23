'''
Created on Sep 16, 2013

@author: supady
'''
from operator import itemgetter
import openbabel
import conversions
import os, sys

def angle_measure(filename, torsion, cistrans):
    molc = openbabel.OBMol()
    conversions.sdf_to_sdf.ReadFile(molc, filename)
    global relax_tor_n
    relax_tor_n=[]
    global relax_cistrans_n
    relax_cistrans_n=[]
    if ( len(cistrans)!=0 ) : 
        for x in range(len(cistrans)):
            relaxcistrans=float(molc.GetTorsion(itemgetter(0)(cistrans[x]),itemgetter(1)(cistrans[x]),itemgetter(2)(cistrans[x]),itemgetter(3)(cistrans[x])))
            relax_cistrans_n.append('{0:.2f}'.format(relaxcistrans))
    if ( len(torsion)!=0) :
        for y in range(len(torsion)):
            relaxtorsions_pre=molc.GetTorsion(itemgetter(0)(torsion[y]),itemgetter(1)(torsion[y]),itemgetter(2)(torsion[y]),itemgetter(3)(torsion[y]))
            relaxtorsions=float(relaxtorsions_pre)
            relax_tor_n.append('{0:.2f}'.format(relaxtorsions))
    if ( len(cistrans)==0 and len(torsion)==0) :
        sys.exit("There are no rotatable bonds and no cis/trans bonds")
def relax_tor():
    return relax_tor_n
def relax_cistrans():
    return relax_cistrans_n    
    
    
def angle_set(molc, torsion, cistrans, values_tor_n, values_cistrans_n):
    
    if ( len(cistrans)!=0 and len(torsion)!=0) : # creates the candidate molecule with previously chosen values 
        for x in range(len(cistrans)):
            for y in range(len(torsion)):
                molc.SetTorsion(itemgetter(0)(torsion[y]),itemgetter(1)(torsion[y]),itemgetter(2)(torsion[y]),itemgetter(3)(torsion[y]),itemgetter(y)(values_tor_n)*0.0174)
                molc.SetTorsion(itemgetter(0)(cistrans[x]),itemgetter(1)(cistrans[x]),itemgetter(2)(cistrans[x]),itemgetter(3)(cistrans[x]),itemgetter(x)(values_cistrans_n)*0.0174)
    elif ( len(cistrans)==0 and len(torsion)!=0) : # only torsions
        for y in range(len(torsion)):
            molc.SetTorsion(itemgetter(0)(torsion[y]),itemgetter(1)(torsion[y]),itemgetter(2)(torsion[y]),itemgetter(3)(torsion[y]),itemgetter(y)(values_tor_n)*0.0174)
    elif ( len(cistrans)!=0 and len(torsion)==0) : # only cis/trans
        for x in range(len(cistrans)):    
            molc.SetTorsion(itemgetter(0)(cistrans[x]),itemgetter(1)(cistrans[x]),itemgetter(2)(cistrans[x]),itemgetter(3)(cistrans[x]),itemgetter(x)(values_cistrans_n)*0.0174)
    else:                 # program closes if there are no degrees of freedom
        sys.exit("There are no rotatable bonds and no cis/trans bonds")
    
