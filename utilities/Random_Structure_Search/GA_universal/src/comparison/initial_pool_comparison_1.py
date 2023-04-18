'''
Created on Oct 16, 2013

@author: newhouse
'''
from comparison.structure_comparison_1 import Comparison


def main(struct, structure_coll, replica):
    '''
    For a structure comparison module, the method main() must be implemented.
    
    The module will take a new structure and compare it to the given structure collection.
    It will return True if the structure passes a test for uniqueness and for energy
    '''
    comp = IPComparison(struct, structure_coll, replica)
    return comp.is_acceptable()

class IPComparison(Comparison):
    '''
    This class only checks the uniqueness and not the energy 
    by extending Comparison and overriding is_acceptable()
    '''
        
    def is_acceptable(self):
        if self.ui.get_eval('comparison', 'always_pass_initial_pool') is True: return True
        return self.is_unique() # and self.acceptable_energy()
    
    def energy_window(self, struct_a, list_to_compare):
        return list_to_compare