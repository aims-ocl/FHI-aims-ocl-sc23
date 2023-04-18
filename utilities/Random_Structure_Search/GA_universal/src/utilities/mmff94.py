'''
Created on Sep 16, 2013

@author: supady
'''
import openbabel, pybel
from pybel import *
import conversions
from energy import *

def energy_mmff94(molc,filename,n):
    """Calculates the force field energy (MMFF 94) of the molecule
    Returns a number.    """
    def isint(x):
        try:
            a = int(x)
        except ValueError:
            return False
        else:
            return True
    ind=filename.split(".")[-2].split("_")[-1] 
    #molc = openbabel.OBMol()
    #conversions.xyz_to_fhiaims.ReadFile(molc, filename)
    ff = openbabel.OBForceField.FindForceField("MMFF94")
    ff.Setup(molc)
    #ff.SteepestDescent(2500, 1.0e-6)
    ff.ConjugateGradients(1000, 1.0e-10) 
    ff.UpdateCoordinates(molc)     
    
    if isint(ind)==True:
        conversions.sdf_to_sdf.WriteFile(molc, "molca_"+str(n)+".sdf")
    elif str(ind)=='child1':
        conversions.sdf_to_sdf.WriteFile(molc, "molca_child1.sdf")                    
    elif str(ind)=='child2':
        conversions.sdf_to_sdf.WriteFile(molc, "molca_child2.sdf")    
        
    energyc=float('{0:.4f}'.format(energy(molc)))   
    return energyc
    
if __name__ == "__main__": 
    mol=  readstring("smi", "N(=Nc1ccccc1)c2ccccc2") # molecule
    #mol.addh()
    mol.make3D()
    print "If the function is implemented correctly, you can read the energy of azobenzene: " + str(energy(mol))
    mol1=  readstring("smi", "c1cccc(c1)N=Nc1ccc(cc1)NC(=S)Nc1cc(cc(c1)C(F)(F)F)C(F)(F)F") # molecule
    mol1.make3D()
    print "If the function is implemented correctly, you can read the energy of catalyst: " + str(energy(mol1))


