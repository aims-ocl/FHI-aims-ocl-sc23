'''
Created on Sep 16, 2013

@author: supady
'''
"""energy.py sample script 
"""

import openbabel, pybel
from pybel import *
import sys

def energy(molecule):
    """Calculates the force field energy (MMFF 94) of the molecule
    Returns a number.    """
    
    ff = pybel._forcefields["mmff94"]
    
    try:  
        success = ff.Setup(molecule.OBMol)    # assigns the force field
    except:
        success = ff.Setup(molecule)
        
    if not success:
        sys.exit("Cannot set up forcefield")
    return ff.Energy()
if __name__ == "__main__": 
    mol=  readstring("smi", "N(=Nc1ccccc1)c2ccccc2") # molecule
    #mol.addh()
    mol.make3D()
    print "If the function is implemented correctly, you can read the energy of azobenzene: " + str(energy(mol))
    mol1=  readstring("smi", "c1cccc(c1)N=Nc1ccc(cc1)NC(=S)Nc1cc(cc(c1)C(F)(F)F)C(F)(F)F") # molecule
    mol1.make3D()
    print "If the function is implemented correctly, you can read the energy of catalyst: " + str(energy(mol1))
