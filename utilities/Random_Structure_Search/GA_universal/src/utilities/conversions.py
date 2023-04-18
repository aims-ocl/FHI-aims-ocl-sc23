'''
Created on Sep 16, 2013

@author: supady
'''
import openbabel
#
# Definitions of different type conversions from openbabel 
#

sdf_to_xyz = openbabel.OBConversion()    
sdf_to_xyz.SetInAndOutFormats("sdf", "xyz")

sdf_to_fhiaims = openbabel.OBConversion()
sdf_to_fhiaims.SetInAndOutFormats("sdf", "fhiaims")

sdf_to_mol2 = openbabel.OBConversion()   # not in use
sdf_to_mol2.SetInAndOutFormats("sdf", "mol2")

sdf_to_sdf = openbabel.OBConversion()
sdf_to_sdf.SetInAndOutFormats("sdf", "sdf")


fhiaims_to_xyz = openbabel.OBConversion()
fhiaims_to_xyz.SetInAndOutFormats("fhiaims", "xyz")

sdf_to_ct = openbabel.OBConversion()
sdf_to_ct.SetInAndOutFormats("sdf", "ct")

xyz_to_sdf = openbabel.OBConversion() # not in use
xyz_to_sdf.SetInAndOutFormats("xyz", "sdf")

xyz_to_fhiaims = openbabel.OBConversion() # not in use
xyz_to_fhiaims.SetInAndOutFormats("xyz", "fhiaims")
