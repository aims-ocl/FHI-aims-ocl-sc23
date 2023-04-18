'''
Created on Oct 16, 2013

@author: newhouse
'''
from __future__ import division

import datetime
import math
import numpy
import time

from core import user_input, output


def main(struct, structure_coll, replica):
    '''
    For a structure comparison module, the method main() must be implemented.
    
    The module will take a new structure and compare it to the given structure collection.
    It will return True if the structure passes a test for uniqueness and for energy
    '''
#     a = datetime.datetime.now()
    comp = Comparison(struct, structure_coll, replica)
    result = comp.is_acceptable()
#     b = datetime.datetime.now()
#     if comp.verbose: self.output("time taken to compare: " + str(struct.struct_id) + ' -- ' + str(b - a))
    return result

class Comparison:
    def __init__(self, struct, structure_coll, replica):
        if struct is False or struct is None: raise Exception
        self.replica = replica
        self.struct_a = struct
        self.structure_coll = structure_coll
        self.ui = user_input.get_config()
        self.verbose = self.ui.get_eval('run_settings', 'verbose')

    def output(self, message): output.local_message(message, self.replica)
        
    def is_acceptable(self):
        # structure has not yet been relaxed. energy is unknown
        # unwise to spend time comparing to existing structures, likely random crossover is not contained here
        # comment out to remove this behaviour
        if self.struct_a.get_property('energy') == None: return True 

        if self.acceptable_energy() is False: 
            return False
        if self.is_unique() is False: return False
        return True
    
    def is_unique(self):
        '''
        determines whether a structure is unique to the current pool
        '''
        list_to_compare_1 = self.structure_coll.structures.values()
        # filter by energy criteria
        list_to_compare_2 = self.energy_window(self.struct_a, list_to_compare_1)
        # filter by histogram
        list_to_compare_3 = self.histogram_comparison(self.struct_a, list_to_compare_2)
        # compare at fine-level the distance array
        if self.verbose: 
            try: self.output('Comparing ' + str(int(len(list_to_compare_3) / len(list_to_compare_1) * 100)) + '% of structures at fine level')
            except: pass
        result = self.compare_distance_array(self.struct_a, list_to_compare_3)
        return result
    
    def compare_distance_array(self, struct_a, list_to_compare):
        '''
        sums the difference in the sorted distances and compares to the user-definied minimum
        '''
        dist_array_tolerance = self.ui.get_eval('comparison', 'dist_array_tolerance')
        da_a = self.calc_distance_array(struct_a)
        for struct_b in list_to_compare:
            if struct_a == struct_b: raise Exception # debugging feature added the wrong structure to comparison list
            da_b = self.calc_distance_array(struct_b)
            dist_sum = 0.0
            if len(da_a) != len(da_b): 
                continue  # different distance array lengths. obviously not the same structure
            for i in range(len(da_a)):
                dist_sum += abs(da_a[i] - da_b[i])
#             if self.verbose: self.output('distance sum = ' + str(dist_sum) + ' dist_array_tolerance = ' + str(dist_array_tolerance))
            if float(dist_sum) < float(dist_array_tolerance):
                if self.verbose: self.output('not unique')
                return False
        if self.verbose: self.output('structure unique')
        return True
        
    def energy_window(self, struct_a, list_to_compare):
        '''
        reduces the list of structures that warrant comparison by filtering those not within
        a certain window of energy defined by the user
        '''
        e_window = self.ui.get_eval('comparison', 'energy_window')
        if e_window == None: return list_to_compare
        return_list = []
        e_a = float(struct_a.get_property('energy'))
        for struct_b in list_to_compare:
            e_b = float(struct_b.get_property('energy'))
            if e_b <= e_a + e_window and e_b >= e_a + e_window: return_list.append(struct_b) 
        return return_list
    
    def histogram_comparison(self, struct_a, list_to_compare):
        '''
        returns list of structures with histograms that are equal to input strucure.
        unpacks the histogram of each structure and compares. can be used for allowing tolerance.
        also allows to make decision based on bin size. string comparison may me quicker. 
        '''
        return_list = []
        histogram_tolerance = self.ui.get_eval('comparison', 'histogram_tolerance')
        h_a = self.get_histogram(struct_a)
        for struct_b in list_to_compare:
            h_b = self.get_histogram(struct_b)
            h_sum = 0
            for i in range(len(h_a)):
                h_sum += abs(h_a[i] - h_b[i])
            if h_sum < histogram_tolerance:
                return_list.append(struct_b)
        return return_list
    
    def get_histogram(self, struct):
        if not 'dist_histogram' in struct.properties: 
            struct.set_property('dist_histogram', self.fuzzy_histogram(self.calc_distance_array(struct)))
        return struct.get_property('dist_histogram')
    
    def acceptable_energy(self):
        '''
        For each structure in the collection, check if the new energy is greater than the maximum existing energy
        '''
        e_a = self.struct_a.get_property('energy')
        if e_a == None: return True  # structure has not yet been relaxed. energy is unknown
        energy_list = []
        for struct in self.structure_coll.itervalues():
            energy_list.append(struct.get_property(('energy')))
        all_energies = int(len(energy_list))
        if all_energies == 0: return True
        # filter for energies lower than e_a
        energy_list_2 = [e_b for e_b in energy_list if e_a > e_b]
        lower_energies = int(len(energy_list_2))
        ratio = lower_energies / all_energies 
#         if ratio < 1.0: # accept energy as long as it's not greater than max
#         if ratio < 0.5: # accept energy as long as it's not greater than half of the existing energies
#         if ratio <= 0.0: # accept energy only if it's lower than all the energies
        if ratio <= self.ui.get_eval('comparison', 'energy_comparison'):  
            return True
        else: 
            if self.verbose: self.output('Energy too high')
            return False
        
    def get_sorted_distance_array(self, struct):
        distance_array = self.get_dist_array(struct)
        distance_array.sort()
        return distance_array
    
#     def get_dist_array(self, struct):
#         '''
#         Returns the distance array from the structures properties.
#         If not found, calculates and adds to properties
#         '''
#         struct.set_property('dist_array', self.calc_distance_array(struct))
#         return struct.get_property('dist_array')
            
    def calc_distance_array(self, struct):  # TODO: this a particular instance of uniqueness checking. move it to a module.
        '''
        creates an UNSORTED list of distances between atoms in the structure
        used for existence checks
        '''
        geo = struct.get_geometry() 
        distance_array = []
        
        for i in range(len(geo)):
            for j in range(len(geo)):
                if i > j:
                    distance = (geo[i]['x'] - geo[j]['x']) ** 2 
                    + (geo[i]['y'] - geo[j]['y']) ** 2 
                    + (geo[i]['z'] - geo[j]['z']) ** 2
                    distance = math.sqrt(distance)
                    # distance_array[i][j] = distance
                    distance_array.append(distance)
        # important to sort distance list
        return distance_array
    
    def fuzzy_histogram(self, distance_array):  # TODO: should the histogram be normalized?
        BIN_SIZE = self.ui.get_eval('comparison', 'bin_size')
        N_BINS = self.ui.get_eval('comparison', 'n_bins')
        clipped_distance_array = numpy.clip(distance_array, 0, BIN_SIZE * N_BINS)  # brings outliers to boundaries
        bins = numpy.arange(0, N_BINS * BIN_SIZE, BIN_SIZE)  # create bins based on user input
        hist = numpy.zeros(N_BINS + 2)  # initialize histogram (first bin is < 0) # bins[i] <=> hist[i+1]
        def fit_in_bin(element, bins):  # find the proper bin in which each element fits
            if element >= bins[-1]: return len(bins) - 1  # gt max
            for i in range(len(bins)):
                if element < bins[i]: return i - 1
        for element in clipped_distance_array:
            closest_bin = fit_in_bin(element, bins)
            center_of_bin = bins[closest_bin] + BIN_SIZE / 2
            # gives element and bin information to a function that determines which bins to increase
            array_to_add = self.hist_function(element, center_of_bin, BIN_SIZE)
            hist[closest_bin] += array_to_add[0]
            hist[closest_bin + 1] += array_to_add[1]
            hist[closest_bin + 2] += array_to_add[2]
        return list(hist)    
    
    def hist_function(self, element, center_of_bin, BIN_SIZE):
        # y = m*x+b # piecewise linear function 
        x = abs(center_of_bin - element)  # symetrical piecewise function
        b = 1
        m = -1 / BIN_SIZE
        k = m * x + b
        l = 1 - k
        return_array = []  # [bin_to_left, center_bin, bin_to_right]
        if element < center_of_bin: return_array = [l, k, 0.0]
        if element > center_of_bin: return_array = [0.0, k, l]
        if element == center_of_bin: return_array = [0.0, 1.0, 0.0]
        return return_array
