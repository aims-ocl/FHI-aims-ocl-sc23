# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:52:55 2014

@author: Bjoern Lange

Visual structure comparison
"""

from optparse import OptionParser
from aims_geometry import aims_geometry

# main
if __name__ =='__main__':
    parser = OptionParser()
    parser.add_option("-r", help="reference sructure", action="store",
          dest="filename1", default="geometry.in", type=str)
    parser.add_option("-i", help="structure", action="store", 
          dest="filename2", default="geometry.in.next_step", type=str)
    parser.add_option("-o", help="output file", action="store",
          dest="filenameOut", default="displacement.dat", type=str)

    (options, args) = parser.parse_args()
    structRef = aims_geometry()
    structRef.read(options.filename1)
    structure = aims_geometry()
    structure.read(options.filename2)  

    diff = structure.coords - structRef.coords   
    
    x = structRef.coords[:,0]
    y = structRef.coords[:,1]
    z = structRef.coords[:,2]
    
    dx = diff[:,0]
    dy = diff[:,1]
    dz = diff[:,2]
    
    #Write result

    structRef.writeDX("geometry.dx")
    f = open(options.filenameOut, 'w')    
    
    nAtoms = structRef.getNTotalAtoms()
    for iAtom in range(nAtoms):
        f.write(" % 10.4f % 10.4f % 10.4f % 10.4f % 10.4f % 10.4f\n" % \
              (x[iAtom],y[iAtom],z[iAtom],dx[iAtom],dy[iAtom],dz[iAtom]))
    

    f.close ()
    
    f = open("displacement.general","w")
    f.write("file = ./%s\n" % options.filenameOut)   
    f.write("points = %i\n" % nAtoms)
    f.write("format = ascii\n")
    f.write("interleaving = field\n")
    f.write("field = locations, displacements\n")
    f.write("structure = 3-vector, 3-vector\n")
    f.write("type = float, float\n")
    f.write("\n")
    f.write("end")

    f.close()
