'''
Created on Dec 16, 2013

@author: newhouse
'''
from __future__ import division

from copy import deepcopy
import copy
import math
import numpy
import os
import time

from core import user_input
from core.file_handler import read_data, cwd
from structures import structure
from structures.structure import Structure, StoicDict
from utilities.element_radii import radii
from utilities import periodic_check


def main(stoic, seed):
    # decides what type of mutation to execute with a weighted random
    rand_struct = RandomStructure(stoic, seed)
    return rand_struct.create_random_structure()
    

def scalmultv(vector, scalar): return [vector[0]*scalar, vector[1]*scalar, vector[2]*scalar]
def addv(*argv):
    rv = [0, 0, 0]
    for vector in argv:
        for i in range(3):
            rv[i] += vector[i]
    return tuple(rv)

class RandomStructure(object):
    '''
    classdocs
    TODO Refactor me!
    '''

    def __init__(self, stoic, seed):
        self.ui = user_input.get_config()
        self.closeness_ratio = self.ui.get_eval('random_gen', 'min_distance_ratio')
        self.stoic = stoic
        self.seed = abs(int(seed.__hash__()))
        
    def get_model(self):
        model_structure_file = self.ui.get_eval('random_gen', 'model_structure')
        model_struct = Structure()
        if model_structure_file is not None:
            model_struct.build_geo_from_atom_file(os.path.join(cwd, model_structure_file))
            new_geo = numpy.array(filter(lambda x: x['fixed'] == False, model_struct.get_geometry()))
            model_struct.geometry = new_geo
        else:
            model_struct.add_lattice_vector(self.ui.get_eval('periodic', 'lattice_vector_a'))
            model_struct.add_lattice_vector(self.ui.get_eval('periodic', 'lattice_vector_b'))
            model_struct.add_lattice_vector(self.ui.get_eval('periodic', 'lattice_vector_c'))
        return model_struct

    def get_atom_list(self, stoic, subtract_stoic):
        atom_list = []
        for element in stoic:
            to_add = [element] * (stoic[element] - subtract_stoic['element']) 
            atom_list.extend(to_add)
        atom_list.sort
        return atom_list  # form: list:[Mg, Mg, O, O, O, O, O]

    def create_random_structure(self):
        substrate = self.get_model()
        self.atom_list = self.get_atom_list(self.stoic, subtract_stoic=substrate.get_stoic())
        self.set_random_array()
        counter = 0
        while True:
            # I'm unsure why the structures sometimes are not filled properly
            # But this loop should account for any anomalies 
            if counter > 10000: return False
            struct = self.generate_structure(substrate)
            if struct == False: counter+=1; continue
            for atom in struct.geometry: 
                if atom['element'] == '': continue
            if periodic_check.is_acceptable(struct) is False: continue
            # if all checks pass
            return struct   
            
    def generate_structure(self, substrate):
        struct = deepcopy(substrate)
        lv_a, lv_b, lv_c = struct.get_lattice_vectors()
        counter = 0
        while True:
            for element in self.atom_list:
                # x, y, z, element
                a_coord = scalmultv(lv_a, self.rand_coordinate())
                b_coord = scalmultv(lv_b, self.rand_coordinate())
                c_coord = scalmultv(lv_c, self.rand_coordinate())
                atom = addv(a_coord, b_coord, c_coord)

                if self.acceptable_distance(atom, element, struct):
                    struct.build_geo_by_atom(atom[0], atom[1], atom[2], element)
                    counter = 0
                else:
                    counter += 1
                if counter > 1000 * len(self.atom_list):
                    # nothing will work. return failure and start again
                    return False
            return struct
    
    def set_random_array(self):
        self.array_index = 0
        # multiply by some number to avoid clashes
        numpy.random.seed(self.seed + int(time.time()*10000))
        self.random_array = numpy.random.random(size=len(self.atom_list) * 3 * 3)
        # times 3 to account for unacceptable coordinates
        
    def rand_coordinate(self):
        try: coordinate = self.random_array[self.array_index]
        except: self.set_random_array(); coordinate = self.random_array[self.array_index]
        self.array_index += 1
        return coordinate
        
    def stoic_to_list(self, stoic):
        '''
        retrieves user input structure in the form Mg:2 O:5 
        and returns elements in usable list form
        '''
        # form: StoicDict(<type 'int'>, {'Mg': 2, 'O': 5})
        atom_list = []
        for element in stoic:
            to_add = [element] * int(stoic.get(element))
            atom_list.extend(to_add)
        atom_list.sort
        return atom_list  # form: list:[Mg, Mg, O, O, O, O, O]
        
    def acceptable_distance(self, atom, element, struct):
        rad_a = radii.get(element)
        for atom_b in struct.get_geometry():
            rad_b = radii.get(atom_b['element'])
            dist = math.sqrt((atom[0] - atom_b['x']) ** 2 + 
                             (atom[1] - atom_b['y']) ** 2 + 
                             (atom[2] - atom_b['z']) ** 2)
            threshold = self.closeness_ratio * (rad_a + rad_b)
            if dist < threshold: return False
        return True
       
if __name__ == '__main__':
    struct = Structure()
    struct.build_geo_whole_atom_format('''
    lattice_vector 3.5 0.0 0.0 
    lattice_vector 0.0 4.0 0.0 
    lattice_vector 0.0 0.0 44.0 
    atom -0.811486 -0.615322 -1.51538 Mg
    atom -0.697116 0.0386657 1.39756 Mg
    atom 0.883206 1.59085 -0.477114 Mg
    atom 1.54927 -1.27103 0.154104 Mg
    atom -0.972383 1.45032 -0.453429 O
    atom -0.732895 -1.63033 0.578495 O
    atom 1.09637 0.519226 1.59259 O
    atom 1.5156 -0.421332 -1.49254 O
    ''')

    struct = main(struct.get_stoic(), 100)
    print struct.get_geometry_atom_format()