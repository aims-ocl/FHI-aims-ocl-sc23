'''
Created on Aug 5, 2013

@author: newhouse
'''
from __future__ import division

from copy import deepcopy
from math import cos, sin
import math
import numpy
import operator
from random import choice
import time

from core import user_input, output
from structures.structure import StoicDict, Structure
from utilities import periodic_check


def main(list_of_structures, targit_stoic, replica):
    '''
    Every crossover module must have a method named "main" that will take
    as arguments a list of 2 or more Structures and a target stoichiometry (in the form of
    a StoicDict). It must return a single Structure with defined geometry
    or return False if the crossover fails and new structures are needed.
    
    This example uses a class to handle the crossover process.
    
    See the Structure class documentation for instructions on easily
    building a Structure's geometry.
    '''
    cross_object = Crossover(list_of_structures[0], list_of_structures[1], targit_stoic, replica)
    return_struct = cross_object.cross()
    if cross_object.verbose and not return_struct is False: 
        cross_object.output('geometries to cross: \n' + list_of_structures[0].get_geometry_atom_format() + \
                       '\n' + list_of_structures[1].get_geometry_atom_format())
        cross_object.output('crossed_geometry: \n' + return_struct.get_geometry_atom_format())
    return return_struct
    
class Crossover(object):
    '''
    Crossover takes two structures, randomly rotates, and splices 
    two halves together to create a new un-relaxed geometry
    '''
    

    def __init__(self, struct_a, struct_b, stoic, replica):
        '''
        __init__ will always run when a class is initialized. 
        This is also how the class can take arguments.
        '''
        self.ui = user_input.get_config()
        self.lattice_vectors = struct_a.get_lattice_vectors()
        # deepcopy makes copy of data instead of altering the obect being referenced
        self.replica = replica
        self.set_geometries(struct_a, struct_b)
        self.target_stoic = deepcopy(stoic)
        self.verbose = self.ui.get_eval('run_settings', 'verbose')

    def output(self, message): output.local_message(message, self.replica)

    def set_geometries(self, struct_a, struct_b):
        self.substrate = Structure()
        self.struct_a = Structure()
        self.struct_b = Structure()
        for atom in struct_a.get_geometry():
            if atom['fixed'] == True:
                self.substrate.build_geo_by_whole_atom(atom)
            else:
                self.struct_a.build_geo_by_whole_atom(atom)
        for atom in struct_b.get_geometry():
            self.struct_b.build_geo_by_whole_atom(atom)
    
    def cross(self):
        '''
        Rotates and combines the structures so that the resulting structure
        matches the target stoichiometry. 
        '''
        counter = 0
        numpy.random.seed(int(time.time() * 10000))
        while True:
            if counter == 1000: return False  # too many attempts
            temp_a = deepcopy(self.struct_a.get_geometry())
            temp_b = deepcopy(self.struct_b.get_geometry())
            
            # generate random cutting plane
            cut_plane_vector = (float(numpy.random.rand(1) * choice([-1, 1])),
                                float(numpy.random.rand(1) * choice([-1, 1])),
                                float(numpy.random.rand(1) * choice([-1, 1])))
            cut_plane_point = list(choice(list(temp_a) + list(temp_a)))[:3]
            # take one hemisphere from each geometry
            temp_a_cut = Structure()
            for atom in temp_a:
                if numpy.dot(cut_plane_vector, numpy.subtract(list(atom)[:3], cut_plane_point)) > 0:
                    temp_a_cut.build_geo_by_whole_atom(atom)
            temp_b_cut = Structure()
            for atom in temp_b:
                if numpy.dot(cut_plane_vector, numpy.subtract(list(atom)[:3], cut_plane_point)) <= 0:
                    temp_b_cut.build_geo_by_whole_atom(atom)
            # must take at least one atom from each structure
            if len(temp_a_cut.get_geometry()) == 0 or len(temp_b_cut.get_geometry()) == 0: counter +=1; continue
            return_struct = Structure()
            return_struct.set_lattice_vectors(self.lattice_vectors)
            return_struct.build_geo_whole(numpy.concatenate((temp_a_cut.get_geometry(), temp_b_cut.get_geometry(), self.substrate.get_geometry())))
            crossed_stoic = return_struct.get_stoic()
            if crossed_stoic == self.target_stoic:
                if periodic_check.is_acceptable(return_struct):
                    return return_struct
            counter += 1
            

if __name__ == '__main__':
    struct_a = Structure()
    struct_a.build_geo_whole_atom_format('''
##########geometry.in############################
 lattice_vector   6.010000      0.0000000E+00  0.0000000E+00
 lattice_vector  0.0000000E+00   12.02000      0.0000000E+00
 lattice_vector  0.0000000E+00  0.0000000E+00   20.00000

 atom   5.581003     -0.3513795       6.050800     H
 atom   5.632699     -0.2467149       7.580442     H
 atom   5.560969      0.2955519       6.776486     O
 atom   2.927149       4.545666       5.755867     H
 atom   3.141190       4.209540       7.237297     H
 atom   2.942419       4.950531       6.639831     O
 atom   4.763853       2.783684       8.411676     H
 atom   4.804386       1.641465       7.388386     H
 atom   4.753829       1.813700       8.344066     O
 atom   1.766136       4.027981       7.288442     H
 atom   1.934766       3.467659       8.706551     H
 atom   1.681143       3.208160       7.804399     O
 atom   1.592302       10.35316       8.659322     H
 atom  0.5542359       9.916296       7.617691     H
 atom   1.509411       9.906470       7.799592     O
 atom   5.165367       12.00099       6.536237     H
 atom   5.344592       10.90052       7.589934     H
 atom   4.663804       11.36354       7.072549     O
 atom  0.6463394       11.84118       8.402984     H
 atom   1.057370       12.59310       7.130547     H
 atom   1.401787       11.98685       7.808326     O
    atom        0.00000        0.00000        5.10000        Mg
    atom        3.00500        0.00000        5.10000        Mg
    atom        0.00000        3.00500        5.10000        Mg
    atom        3.00500        3.00500        5.10000        Mg
    atom        0.00000        6.01000        5.10000        Mg
    atom        3.00500        6.01000        5.10000        Mg
    atom        0.00000        9.01500        5.10000        Mg
    atom        3.00500        9.01500        5.10000        Mg
    atom        1.50250        1.50250        5.10000        O
    atom        4.50750        1.50250        5.10000        O
    atom        1.50250        4.50750        5.10000        O
    atom        4.50750        4.50750        5.10000        O
    atom        1.50250        7.51250        5.10000        O
    atom        4.50750        7.51250        5.10000        O
    atom        1.50250        10.51750    5.10000        O
    atom        4.50750        10.51750    5.10000        O
    atom        0.00000        0.00000        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        0.00000        3.00500        O
        constrain_relaxation    .true.
    atom        0.00000        3.00500        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        3.00500        3.00500        O
        constrain_relaxation    .true.
    atom        0.00000        6.01000        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        6.01000        3.00500        O
        constrain_relaxation    .true.
    atom        0.00000        9.01500        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        9.01500        3.00500        O
        constrain_relaxation    .true.
    atom        1.50250        1.50250        3.00500        Mg
        constrain_relaxation    .true.
    atom        4.50750        1.50250        3.00500        Mg
        constrain_relaxation    .true.
    atom        1.50250        4.50750        3.00500        Mg
        constrain_relaxation    .true
    atom        4.50750        4.50750        3.00500        Mg
        constrain_relaxation    .true.
    atom        1.50250        7.51250        3.00500        Mg
        constrain_relaxation    .true.
    atom        4.50750        7.51250        3.00500        Mg
        constrain_relaxation    .true.
    atom        1.50250        10.51750    3.00500        Mg
        constrain_relaxation    .true.
    atom        4.50750        10.51750    3.00500        Mg
        constrain_relaxation    .true. 
    ''')

    struct_b = Structure()
    struct_b.build_geo_whole_atom_format('''
##########geometry.in############################
 lattice_vector   6.010000      0.0000000E+00  0.0000000E+00
 lattice_vector  0.0000000E+00   12.02000      0.0000000E+00
 lattice_vector  0.0000000E+00  0.0000000E+00   20.00000

 atom   5.581003     -0.3513795       6.050800     H
 atom   5.632699     -0.2467149       7.580442     H
 atom   5.560969      0.2955519       6.776486     O
 atom   2.927149       4.545666       5.755867     H
 atom   3.141190       4.209540       7.237297     H
 atom   2.942419       4.950531       6.639831     O
 atom   4.763853       2.783684       8.411676     H
 atom   4.804386       1.641465       7.388386     H
 atom   4.753829       1.813700       8.344066     O
 atom   1.766136       4.027981       7.288442     H
 atom   1.934766       3.467659       8.706551     H
 atom   1.681143       3.208160       7.804399     O
 atom   1.592302       10.35316       8.659322     H
 atom  0.5542359       9.916296       7.617691     H
 atom   1.509411       9.906470       7.799592     O
 atom   5.165367       12.00099       6.536237     H
 atom   5.344592       10.90052       7.589934     H
 atom   4.663804       11.36354       7.072549     O
 atom  0.6463394       11.84118       8.402984     H
 atom   1.057370       12.59310       7.130547     H
 atom   1.401787       11.98685       7.808326     O
    atom        0.00000        0.00000        5.10000        Mg
    atom        3.00500        0.00000        5.10000        Mg
    atom        0.00000        3.00500        5.10000        Mg
    atom        3.00500        3.00500        5.10000        Mg
    atom        0.00000        6.01000        5.10000        Mg
    atom        3.00500        6.01000        5.10000        Mg
    atom        0.00000        9.01500        5.10000        Mg
    atom        3.00500        9.01500        5.10000        Mg
    atom        1.50250        1.50250        5.10000        O
    atom        4.50750        1.50250        5.10000        O
    atom        1.50250        4.50750        5.10000        O
    atom        4.50750        4.50750        5.10000        O
    atom        1.50250        7.51250        5.10000        O
    atom        4.50750        7.51250        5.10000        O
    atom        1.50250        10.51750    5.10000        O
    atom        4.50750        10.51750    5.10000        O
    atom        0.00000        0.00000        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        0.00000        3.00500        O
        constrain_relaxation    .true.
    atom        0.00000        3.00500        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        3.00500        3.00500        O
        constrain_relaxation    .true.
    atom        0.00000        6.01000        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        6.01000        3.00500        O
        constrain_relaxation    .true.
    atom        0.00000        9.01500        3.00500        O
        constrain_relaxation    .true.
    atom        3.00500        9.01500        3.00500        O
        constrain_relaxation    .true.
    atom        1.50250        1.50250        3.00500        Mg
        constrain_relaxation    .true.
    atom        4.50750        1.50250        3.00500        Mg
        constrain_relaxation    .true.
    atom        1.50250        4.50750        3.00500        Mg
        constrain_relaxation    .true
    atom        4.50750        4.50750        3.00500        Mg
        constrain_relaxation    .true.
    atom        1.50250        7.51250        3.00500        Mg
        constrain_relaxation    .true.
    atom        4.50750        7.51250        3.00500        Mg
        constrain_relaxation    .true.
    atom        1.50250        10.51750    3.00500        Mg
        constrain_relaxation    .true.
    atom        4.50750        10.51750    3.00500        Mg
        constrain_relaxation    .true. 
    ''')
    crossed = main([struct_a, struct_b], struct_a.get_stoic(), 1)
    print crossed.get_geometry_atom_format()