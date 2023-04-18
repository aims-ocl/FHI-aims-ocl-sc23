'''
Created on Sep 16, 2013

@author: supady
'''
import numpy
#import scipy 
#import sympy



## Returns 'True' if the geometry is acceptable.

def check_geometry(filename, atoms, bonds, matrix_of_connections):
    
    array=numpy.loadtxt(filename, comments='#', delimiter=None, converters=None, skiprows=2, usecols=xrange(1,4), unpack=False )
        
    #dist=scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(array))
    
    dist=numpy.zeros( (len(range(atoms)),len(range(atoms))))
    def distance(x,y):
        """ Defines euclidean distances between a and b 
        """
        return numpy.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)
    
    for x in range(atoms):
        for y in range(atoms):
            a=numpy.array(array[float(x)])
            b=numpy.array(array[float(y)])
            #dist1[x][y]=distance(a,b)
            dist[x][y]=distance(a,b)
    
    matrix=matrix_of_connections
    tocheck=dist+matrix 
       
    
    b=0
    c=0
    condition_clashes=True
    for b in range(atoms):
        for c in range(atoms): 
            dis=tocheck[b][c]
        
            if (dis < 1.2 and dis >0): # if a distance between to atoms lower 1.2 is found the loop breaks and the condtition is set to False 
                condition_clashes=False
                break 
            else: 
                continue
            break
    if (condition_clashes==True):
        #A=sympy.Matrix(matrix)
        #B=sympy.Matrix(dist)
        #tocheck2=A.multiply_elementwise(B)
      
        tocheck2=numpy.multiply(matrix, dist)
       
        
      
        e=0
        f=0
        condition_bonded=True
        for e in range(atoms):
            for f in range(atoms):
                dis2=tocheck2[e,f]
                if (dis2>2.0): # the distance between two bonded atoms shouldn't be loonger than 2.0
                    condition_bonded=False
                    break
                else:
                    continue
                break
        if ( condition_bonded==True) :
            check=True
                    
            
        else:
            check=False
    else:
        check=False
    return check
            







