'''
Created on Oct 16, 2013

@author: newhouse
'''
from __future__ import division

from collections import defaultdict
from compiler.ast import Print
import math
import numpy
from random import choice
import random
import time

from core import user_input, output
from core.file_handler import INITIAL_POOL_REFID
from structures import structure_collection
from structures.structure import StoicDict


def main(structure_supercoll, replica_stoic, replica):
    '''
    For a structure selection module, the method main() must be implemented.
    
    It must take as arguments a 'Structure supercollection' (a StoicDict of structure_collection)
    and a 'replica_stoic' which defines the main stoichiometry of the replica calling this module.
    
    It must return a list of 2 or more Structure objects intended for crossing.
    
    This module has the responsibility of selecting stoichiometries in cross-stoic crossings,
    calculating fitness for those stoichiometries, and selecting a structure from each.
    
    The following is an example of a possible implementation
    '''
    struct_sel = StructureSelection(replica)
    structure_collection.update_supercollection(structure_supercoll)
    while True:
        try:
            struct_sel.select_collections(structure_supercoll, replica_stoic)
            # for the special case where there are more simultaneous replicas than initial pool structures
            if len(struct_sel.structure_coll_a) == 0 or len(struct_sel.structure_coll_b) == 0: 
                output.local_message('No structures in collection, selecting random initial pool structures', struct_sel.replica)
                struct_a = choice(structure_supercoll.get((struct_sel.structure_coll_a.get_stoic(), INITIAL_POOL_REFID)).get_structures().values())
                struct_b = choice(structure_supercoll.get((struct_sel.structure_coll_b.get_stoic(), INITIAL_POOL_REFID)).get_structures().values())
                break
            # else proceed as normally
            [struct_a, struct_b] = struct_sel.get_structures()
            break
        except: pass  # may look in a collection which has nothing in it, select collections again
    return [struct_a, struct_b]


class StructureSelection():
    
    def __init__(self, replica):
        self.ui = user_input.get_config()
        self.replica = replica

    def select_collections(self, structure_supercoll, replica_stoic):
        '''filters for a particular input reference so structures selected will be of the same group'''
        control_list = self.ui.get_list('control', 'control_in_filelist')
        self.MAX_CASCADE = len(control_list) - 1
        self.key_list = filter(lambda x : x[1] == self.MAX_CASCADE, structure_supercoll.iterkeys()) 
        self.replica_stoic = replica_stoic       
        self.structure_coll_a = structure_supercoll[self.select_stoichiometry()]
        self.structure_coll_b = structure_supercoll[self.select_stoichiometry()]        
        
    def select_stoichiometry(self):
        '''
        this method will return a key to a structure supercollection. The decision will be based on a probability
        distribution and the availibility of other stoichiometries. In the single-stoic case the structure matching
        the replica stoichiometry will be selected.
        '''
        if len(self.key_list) == 1: return self.key_list[0]
#         return choice(self.key_list)
        # create a collection of information of possible elements and thier possibe quantities
        possible_elements = defaultdict(set)
        for stoic, input_ref in self.key_list:
            for element, quantity in stoic.items():
                quantity_set = possible_elements[element]
                quantity_set.update([quantity])
                possible_elements[element] = quantity_set
        del element, quantity
        
        while True:
        # do until a suitable stoichiometry is chosen
            stoic_to_choose = StoicDict()
            # select quantity of an element
            stdev = self.ui.get_eval('selection', 'stoic_stdev')
            for element, quantity_set in possible_elements.items():
                # use a normal distribution to specify weights
                weights = {}
                if element in self.replica_stoic: mean = self.replica_stoic[element]
                else: mean = 0
                for x in quantity_set:
                    y = math.exp(-((x - mean) ** 2) / (2 * stdev ** 2)) / (math.sqrt(2 * math.pi) * stdev)
                    weights[x] = y
                quantity = self.weighted_choice(weights)
                if not quantity == 0: stoic_to_choose[element] = quantity
            for key in self.key_list:
                if stoic_to_choose == key[0]: return key
            
    def weighted_choice(self, choices):
        total = sum(w for w in choices.values())
        r = random.uniform(0, total)
        upto = 0
        for c, w in choices.items():
            if upto + w > r:
                return c
            upto += w
        assert False, "Shouldn't get here"
            
    def get_structures(self):
        # select which stoichiometries will be considered. will be same if only one stoic present.
        if self.structure_coll_a == self.structure_coll_b:
            # if structure collection has one structure, return it twice
            if len(self.structure_coll_a) == 1: return [self.structure_coll_a.get_struct(self.structure_coll_a.structures.keys()[0])] * 2
            # if two structures, return them both
            if len(self.structure_coll_a) == 2: return [self.structure_coll_a.get_struct(self.structure_coll_a.structures.keys()[0]),
                                                        self.structure_coll_a.get_struct(self.structure_coll_a.structures.keys()[1])]
        # select structures from each stoichiometry
        struct_a, fit_a = self.select_structure(self.structure_coll_a)
        while True: 
            struct_b, fit_b = self.select_structure(self.structure_coll_b)
            if not struct_a == struct_b: break  # remove this to allow double-selection
        return [struct_a, struct_b]

    
    def select_structure(self, structure_coll):
        '''
        Will calculate fitness for each structure in a collection and select one structure based on a distribution
        
        An optional feature can select from a lower fitness structure. Yet to be implemented.
        '''
        # take care of single-structure case. return only structure
        if len(structure_coll.structures) == 1: 
            return (structure_coll.get_struct(0), 0)
        else: 
            fitness_dict = self.get_energy_fitness(structure_coll)  # this can be altered to select for different fitness
            return self.select_best_from_fitness(fitness_dict)  # TODO: change the selection method
            
    def get_energy_fitness(self, structure_coll):
        '''
        Takes a structure collection of 2 or more structures and returns a dictionary of fitness based on energy
        '''
        # decides whether to reverse fitness values
        numpy.random.seed(int(time.time() * 83 + len(structure_coll.structures) * 19))
        reverse = numpy.random.random() < self.ui.get_eval('selection', 'fitness_reversal_probability')
        # find min and max energy
        # random energy in collection just a base for comparison
        min_e = max_e = structure_coll.itervalues().next().get_property('energy')
        for index, struct in structure_coll:
            try: e = float(struct.get_property('energy'))
            except: 
                message = 'no energy: structure ' + str(struct.get_struct_id()) + '\n'
                message += str(struct.properties)
                output.error(message, self.replica)
            if max_e < e: max_e = e
            if min_e > e: min_e = e
    
        fitness = {}
        for index, struct in structure_coll:
            try: energy = float(struct.get_property('energy'))
            except: pass
            rho = (max_e - energy) / (max_e - min_e)
            if reverse: rho = 1 - rho
            if self.ui.get('selection', 'fitness_function') == 'standard':
                fitness[struct] = rho
            if self.ui.get('selection', 'fitness_function') == 'exponential':
                fitness[struct] = math.exp(-self.ui.get('selection', 'alpha') * rho)  # TODO: implement fitness
        return fitness
        
    def sorted_fitness(self, fitness_dict):
        '''
        returns fitness as a sorted list of tuples.
        '''
        sorted_fitness = sorted(fitness_dict.iteritems(), key=lambda x:x[1])
        # sorted_fitness = sorted_fitness[::-1] # sort fitness 1 to 0
        return sorted_fitness
    
    def normalized_fitness(self, sorted_fit):
        '''
        takes a sorted list of fitness tuples and returns the list of tuples normalized
        '''
        # sum all fitness
        f_sum = 0
        # reduce all fitnesses
        for index, fitness in sorted_fit: f_sum += fitness
        tmp = []
        for index, fitness in sorted_fit: tmp.append((index, fitness / f_sum))
        
        total = 0
        normalized_fit = []
        for index, t_fitness in tmp: 
            normalized_fit.append((index, t_fitness + total))
            total = t_fitness + total
        return normalized_fit
    
    def select_best_from_fitness(self, fitness_dict):
        fitness = self.sorted_fitness(fitness_dict)
        fitness = self.normalized_fitness(fitness)
        
        numpy.random.seed(int(time.time() * 87 + len(fitness) * 11))
        random_num = numpy.random.random()
        # selects first element to be greater than random number
        return list(filter((lambda x : x[1] > random_num), fitness))[0]
        
                      
    
    
    
    
