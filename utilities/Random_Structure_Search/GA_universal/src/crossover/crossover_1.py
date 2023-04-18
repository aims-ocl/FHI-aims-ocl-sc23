'''
Created on Aug 5, 2013

@author: newhouse
'''
from __future__ import division

from copy import deepcopy
from math import cos, sin
import math
import numpy

from core import user_input, output
from structures.structure import StoicDict, Structure


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
#     if cross_object.verbose and not return_struct is False: 
#         cross_object.output('geometries to cross: \n' + list_of_structures[0].get_geometry_atom_format() +\
#                        '\n' + list_of_structures[1].get_geometry_atom_format())
#         cross_object.output('crossed_geometry: \n' + return_struct.get_geometry_atom_format())
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
        self.replica = replica
        # deepcopy makes copy of data instead of altering the obect being referenced
        self.geometry_a = deepcopy(struct_a.get_geometry())
        self.geometry_b = deepcopy(struct_b.get_geometry())
        self.target_stoic = deepcopy(stoic)
        self.verbose = self.ui.get_eval('run_settings', 'verbose')

    def output(self, message): output.local_message(message, self.replica)

    def cross(self):
        '''
        Rotates and combines the structures so that the resulting structure
        matches the target stoichiometry. 
        '''
        counter = 0
        while True:
            # generate random rotation angles
            temp_a = self.randomly_rotate_structure(deepcopy(self.geometry_a))
            temp_b = self.randomly_rotate_structure(deepcopy(self.geometry_b))

            # center the structures
            self.center_geometry(temp_a)
            self.center_geometry(temp_b)
            
            # order atoms by z coordinate
            temp_a.sort(order='z')
            temp_b.sort(order='z')

            # alternates contraction and expansion of atoms.
            to_contract = 'b'
            to_expand = 'a'
            
            # take one hemisphere from each geometry
            temp_a_cut = temp_a[temp_a['z'] > 0]
            temp_b_cut = temp_b[temp_b['z'] <= 0]
            
            while True:
                if counter == 1000: return False  # too many attempts


                # calculate resulting stoichiometries
                stoic_a = self.calc_stoic(temp_a_cut)
                stoic_b = self.calc_stoic(temp_b_cut)
                
                # decides whether to expand or contract the structures to match target stoic
                expand = False
                contract = False
                # to avoid bais towards a structure
                rand_decision = 1 * numpy.random.randint(0, 2)
                for element in self.target_stoic:
                # checks stoichiomentry for each element
                    combined = stoic_a.get(element, 0) + stoic_b.get(element, 0)
                    if self.target_stoic.get(element) < combined: contract = True
                    if self.target_stoic.get(element) > combined: expand = True
                
                # Acceptable stoichometry
                if expand == contract == False:
                    # checks if the atoms at interface are too close
                    if not self.is_too_close(temp_a_cut, temp_b_cut):
                        # crossover successful. return the joined geometries as a structure ############ SUCCESS ############
                        new_struct = Structure()
                        new_struct.build_geo_whole(numpy.concatenate((temp_a_cut, temp_b_cut)))
                        return new_struct
                    else:
#                         # if too close attempt to increase the distance slightly
                        counter += 1
                        # separate slightly
                        if rand_decision < 1: temp_a = self.separate_slightly(temp_a)
                        else: temp_b = self.separate_slightly(temp_b)
                        # TODO: is this really the best way? this includes
                        # new atoms if they are too close to the interface. 
                        # Maybe I should not do it this way
                        continue
                    
                # stoichometries conflict. at least one has too few while one has too many.         
                if expand == contract == True: 
                    counter += 1
                    break  # restart crossover with new rotation  
                
                # Needs more atoms
                if expand == True and contract == False: 
                    counter += 1                    
                    if to_expand == 'a':  # alternate which structure to modify
                        to_expand = 'b'
                        temp_a_cut = self.increment_z(temp_a_cut)
                    else: 
                        to_expand = 'b'
                        temp_b_cut = self.decrement_z(temp_b_cut)
                    if temp_a_cut == False or temp_b_cut == False: break 
                    continue
                
                # Needs fewer atoms
                if contract == True and expand == False:
                    counter += 1
                    if to_contract == 'a':  # alternate which structure to modify
                        to_contract = 'b'
                        temp_a_cut = self.decrement_z(temp_a_cut)
                    else: 
                        to_contract = 'a'
                        temp_b_cut = self.increment_z(temp_b_cut)
                    if temp_a_cut == False or temp_b_cut == False: break 
                    continue
        
    
    def randomly_rotate_structure(self, geometry):
        '''
        uses a randomly generated rotation matrix to rotate each atom
        of a geometry in 3 dimensions
        '''
        rot_matrix = self.get_rot_matrix()
        for i in range(len(geometry)):
            coordinate = numpy.matrix([    [geometry[i][0]],
                                        [geometry[i][1]],
                                        [geometry[i][2]]
                                   ])
            rot_coordinate = rot_matrix * coordinate
            geometry[i][0] = float(rot_coordinate[0])
            geometry[i][1] = float(rot_coordinate[1])
            geometry[i][2] = float(rot_coordinate[2])
        return geometry

    def get_rot_matrix(self):
        '''
        returns a rotation matrix used in crossover calculations
        '''
        phi = numpy.random.rand(1) * math.pi * 2
        psi = numpy.random.rand(1) * math.pi * 2
        theta = numpy.random.rand(1) * math.pi * 2
        # TODO: should theta really be 2pi?

        rotation_matrix = numpy.matrix([ ((cos(theta) * cos(psi)),
                             (-cos(phi) * sin(psi)) + (sin(phi) * sin(theta) * cos(psi)),
                             (sin(phi) * sin(psi)) + (cos(phi) * sin(theta) * cos(psi))),
                            
                             ((cos(theta) * sin(psi)),
                             (cos(phi) * cos(psi)) + (sin(phi) * sin(theta) * sin(psi)),
                             (-sin(phi) * cos(psi)) + (cos(phi) * sin(theta) * sin(psi))),
                            
                             ((-sin(theta)),
                             (sin(phi) * cos(theta)),
                             (cos(phi) * cos(theta)))])
        return rotation_matrix
        
    def center_geometry(self, geometry):
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
        
    def calc_stoic(self, geo):
        '''
        returns a dictionary representing the stoichiometries
        '''
        stoic = StoicDict(int)
        for item in geo:
            stoic[item[3]] += 1
        return stoic    

    def is_too_close(self, geometry_a, geometry_b):
        '''
        only compares the distances between structures and not within each half.
        relaxed geometries may have distances smaller than the minimum distance.
        '''
        is_too_close = False
        min_dist = float(self.ui.get('crossover', 'crossover_minimum_interface_distance'))
        for atom_a in geometry_a:
            ax = atom_a[0]
            ay = atom_a[1]
            az = atom_a[2]
            for atom_b in geometry_b:
                bx = atom_b[0]
                by = atom_b[1]
                bz = atom_b[2]
                dist = numpy.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
                if dist < min_dist: 
                    is_too_close = True
        return is_too_close
    
    def separate_slightly(self, geometry):
        '''
        designed to move away from z axis whether positive or negative
        '''
        fraction_of_min_dist = float(self.ui.get('crossover', 'crossover_minimum_interface_distance')) * 0.1
        z_sum = 0
        for atom in geometry:
            z_sum += atom['z']
        if sum < 0: sign = -1
        else: sign = 1
        for atom in geometry:
            atom['z'] += fraction_of_min_dist * sign
        return geometry

    def increment_z(self, geometry):
        # find first atom below 0
        search_index = 0
        for atom in geometry:
            if atom['z'] < 0: break
            search_index += 1
        if search_index >= len(geometry): return False
        z_increase = -geometry[search_index]['z']
        for atom in geometry:
            atom['z'] += z_increase
        return geometry
    
    def decrement_z(self, geometry):
        # find first atom above 0
        search_index = 0
        geometry = geometry[::-1]  # reverse
        for atom in geometry:
            if atom['z'] > 0: break
            search_index += 1
        if search_index >= len(geometry): return False
        z_increase = -geometry[search_index]['z']
        for atom in geometry:
            atom['z'] += z_increase
        return geometry[::-1]
