'''
Created on Aug 7, 2013

@author: newhouse
'''
from __future__ import division

from copy import deepcopy
import numpy
import random
import time

from core import user_input
from structures.structure import Structure


def main(struct, r_stoic, replica):
    '''
    Every mutation module must have a method named "main" that will take
    as arguments one Structure and a target stoichiometry (in the form of
    a StoicDict). It must return a single Structure with defined geometry
    or return False if the mutation fails and a new Structure is needed.
    
    This example uses classes to handle the mutation process.
    A decision can be made in the main() function about which mutation to execute.
    
    See the Structure class documentation for instructions on easily
    building a Structure's geometry.
    '''
    # deepcopy makes copy of data instead of altering the obect being referenced
    input_structure = deepcopy(struct)
    replica_stoic = deepcopy(r_stoic)
    
    # decides what type of mutation to execute with a weighted random
    mutate_object = select_mutator(input_structure, replica_stoic, replica)
    return input_structure
    return mutate_object.mutate()

def select_mutator(input_structure, replica_stoic, replica):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method is responsible for reading the preferences set
    by the user and selecting which mutation to empoly, or no mutation at all.
    Expects: Structure, target_stoic
    Returns: Class
    '''
    # TODO: implement
    # SCAFFOLDING
    # read_UI_preferences
    # check if structure stoic is different than target stoic
        # if different, assign to adder or subtracter
    # decide based on random weighted with preferences
    return RandomTranslationMutation(input_structure, replica_stoic, replica)

def get_target_stoic(struct):
    '''
    This method is optional and intended for the case where mutation will add or remove
    one or more atoms from the input structure. The decided upon target stoichiometry
    will be returned to the genetic algorithm module for use in the crossover module.
    
    For example, if one intends to add 2 oxygen atoms to the resulting crossover,
    this method must return the stoichiometry of the input structure MINUS 2 oxygens
    
    The intent is to maintain that each replica of the GA is responsible for one stoichiometry
    
    If this method is not implemented, the GA will simply use its defined replica stoic in crossover
    '''
    # TODO: implement
    # SCAFFOLDING
    # read_UI_preferences
    # decide based on random weighted with preferences
    return struct.get_stoic()
    

class RandomTranslationMutation(object):

    def __init__(self, input_structure, target_stoic, replica):
        '''
        __init__ will always run when a class is initialized. 
        This is also how the class can take arguments.
        '''
        
        self.geometry = deepcopy(input_structure.get_geometry())
        self.ui = user_input.get_config()
        self.input_structure = input_structure
        # based on observation of perl code
        # TODO: check if can be mutated (forbidden_eles) etc.
        
    def mutate(self):
        return self.random_translation()  # TODO: make this a test based on mutation weights
#         return self.switch_species()
        pass

    def random_translation(self):
        ''' 
        make a random translation within reasonable bounds and returns Structure
        '''
        temp_geo = center_geometry(self.geometry)  # center the geometry
#         jmol(temp_geo, 'a')
        # find most external atom
        max_dist_x = 0
        max_dist_y = 0
        max_dist_z = 0
        for atom in temp_geo:
            if abs(atom['x']) > max_dist_x: max_dist_x = abs(atom['x'])
            if abs(atom['y']) > max_dist_y: max_dist_y = abs(atom['y'])
            if abs(atom['z']) > max_dist_z: max_dist_z = abs(atom['z'])
        max_dist_x = max_dist_x * 1.2  # increase the maximum sphere slightly. subject to change
        max_dist_y = max_dist_y * 1.2  # increase the maximum sphere slightly. subject to change
        max_dist_z = max_dist_z * 1.2  # increase the maximum sphere slightly. subject to change

        counter = 0
        while True:
            if counter > 10000: 
                # mutation failed
                return False  
            if counter > 1000:
                # try bigger sphere
                max_dist_x = max_dist_x * 1.2  # increase the maximum sphere slightly. subject to change
                max_dist_y = max_dist_y * 1.2  # increase the maximum sphere slightly. subject to change
                max_dist_z = max_dist_z * 1.2  # increase the maximum sphere slightly. subject to change
            counter += 1
            # select a random index to translate
            translate_index = random.randint(0, len(temp_geo) - 1)
            # skip forbidden elements
            if temp_geo[translate_index]['element'] in self.ui.get_list('mutation', 'forbidden_elements'): continue
            
            # this is strange, but with simple float i get an error
            seed = (abs(self.geometry[0][0]*23) + time.time()/1041).__hash__()            
            # pick a random coordinate within +/- max distance 
            r_x, r_y, r_z = [rand_coordinate(max_dist_x, seed), 
                             rand_coordinate(max_dist_y, seed+33), 
                             rand_coordinate(max_dist_z, seed+159)]  
            
            # this avoids comparing the new atom to its original position in closeness calculations
            compare_geo = numpy.delete(temp_geo, translate_index)
            
            # checks if too close
            if not is_too_close(r_x, r_y, r_z, compare_geo, float(self.ui.get('mutation', 'minimum_bond_length'))):
                temp_geo[translate_index]['x'] = r_x
                temp_geo[translate_index]['y'] = r_y
                temp_geo[translate_index]['z'] = r_z
                self.geometry = temp_geo
#                 jmol(self.geometry,'b')
                new_struct = Structure()
                new_struct.build_geo_whole(temp_geo)
                new_struct.set_lattice_vectors(self.input_structure.get_lattice_vectors())
                return new_struct

    def switch_species(self):
        forbidden_elements = self.ui.get_list('mutation', 'forbidden_elements')
        temp_geo = self.geometry
        indices = list(range(len(temp_geo)))
        random.shuffle(indices)
        
        to_switch_1 = indices.pop(0)
        while True:
            if temp_geo[to_switch_1]['element'] not in forbidden_elements:
                break
            to_switch_1 = indices.pop(0)
        
        to_switch_2 = -1
        for item in indices:
            if temp_geo[item]['element'] not in forbidden_elements and \
            temp_geo[item]['element'] not in temp_geo[to_switch_1]['element']:  # not switching the same element
                to_switch_2 = item
                break
        if to_switch_2 == -1: return False

        temp_geo[to_switch_1]['element'], temp_geo[to_switch_2]['element'] = \
        temp_geo[to_switch_2]['element'], temp_geo[to_switch_1]['element']
        
        self.geometry = temp_geo
        return temp_geo
        pass

    def rotate_selection(self):
        pass

def rand_coordinate(max_dist, seed):
    numpy.random.seed(seed)
    coordinate = numpy.random.random(size=None) * max_dist
    sign = 1 * numpy.random.randint(0, 2)
    if sign < 1: sign = -1
    return coordinate * sign    
        
def center_geometry(geometry):
    '''
    the centers the origin in relation to the max and min of each axis
    '''
    for i in range(2):  # x, y, and z
        coordinate_sum = 0
        counter = 0
        for atom in geometry:
            coordinate_sum += atom[i]
            counter += 1
        # average axis value
        average = coordinate_sum / counter 
        for atom in geometry:
            # shift all towards center
            atom[i] = atom[i] - average 
    return geometry

def is_too_close(ax, ay, az, geometry, min_dist):
    is_too_close = False
    for item in geometry:
        if item['element'] == '': continue
        bx = item['x']
        by = item['y']
        bz = item['z']
        dist = numpy.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
        if dist < min_dist: 
            is_too_close = True
    return is_too_close
