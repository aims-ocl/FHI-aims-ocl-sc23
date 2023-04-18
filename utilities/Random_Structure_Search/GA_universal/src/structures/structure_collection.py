'''
Created on Aug 2, 2013

@author: newhouse
'''
from ast import literal_eval
import os
import sqlite3
import time

from core import user_input
from core.file_handler import structure_dir, mkdir_p, write_data, read_data, \
    get_random_index
from structures.structure import Structure, StoicDict


class StructureCollection(object):
    '''
    This class handles the functionality for groups of structures. It is an aggregation of Structure
    ''' 
    # TODO: may need to restructure supercollection a little and have it add a level of dict. so: 
    # supercollection[STOIC][INPUTREF][STRUCTURE] or supercollection[INPUTREF][STOIC][STRUCTURE]    
    def __init__(self, stoichiometry, input_ref):
        '''
        Constructor
        '''
        # structures stored as a dictionary
        self.structures = {}  # {int, Structure}
        # a structure collection is unique to its stoichiometry and the input file common to each
        # structue. The input files are stored in the database.
        self.stoichiometry = stoichiometry
        self.input_ref = int(input_ref)
        self.path = os.path.join(structure_dir, stoichiometry.get_string(), str(input_ref))
        if not os.path.exists(self.path): mkdir_p(self.path)
        # read allupdate_local_from_db data into local memory
        self.update_local()
        
    # partially extend dictionary functionality. Returns iterators for structure dictionary
    def __iter__(self): return self.structures.iteritems()
    def iteritems(self): return self.structures.iteritems()
    def iterkeys(self): return self.structures.iterkeys()
    def itervalues(self): return self.structures.itervalues()
    def __len__(self): return self.structures.__len__()
        
    def get_struct(self, struct_index): return self.structures.get(struct_index)
    def get_structures(self): return self.structures
    def get_stoic(self): return self.stoichiometry
    def get_stoic_str(self): return self.stoichiometry.get_string()
    def get_key(self): return (self.get_stoic, int(self.input_ref))
    def get_input_ref(self): return self.input_ref
    
    def update_local(self):
        '''
        Reads structures from database and updates local object list with new structures
        This is vital for parallelism when other replicas are populating the same list
        '''
        dirs = [ name for name in os.listdir(self.path) \
                if os.path.isdir(os.path.join(self.path, name)) ]
        # for new indicies, add information to local memory
        for sdir in dirs: 
            struct_id = sdir
            if struct_id not in self.structures.keys():
                counter = 0
                while True:  # other structure may not have finished writing all of it's data. wiat for a second
                    try: 
                        self.__add_structure_to_local(struct_id)
                        break
                    except: 
                        counter += 1
                        if counter>1000: raise Exception
                        time.sleep(0.1)

    def __add_structure_to_local(self, struct_id):
        """
        Uses json loads to read data from file and create a structure with those properties
        """
        # TODO: may need to surround in try/except to catch when file is not fully written
        struct = Structure()
        struct.loads(read_data(os.path.join(self.path, str(struct_id)), 'data.json'))
        self.structures[struct_id] = struct

    def get_unrelaxed_structure(self):
        for struct in self.structures.itervalues():
            # check local memory
            if struct.get_property('to_be_relaxed'):
                # update and check again
                struct.update_local_data()
                if struct.get_property('to_be_relaxed'):
                    # update local and shared data and return the structure
                    struct.set_property('to_be_relaxed', False)
                    struct.update_shared_data()
                    return struct
            else: continue
        return False  # none found
    

def add_structure(struct, stoichiometry, input_ref):
    '''
    Inserts a structure and its properties into the shared filesystem
    '''
    # just to ensture backwards references are kept
    struct.input_ref = input_ref
    struct.stoichiometry = stoichiometry
    
    path = os.path.join(structure_dir, stoichiometry.get_string(), str(input_ref))
    # gets a new index while file-locked to avoid overwriting
    index = get_random_index(os.path.join(path, 'index.dat'))
    struct.struct_id = index
    # writes data for sharing and geometry for human readability
    mkdir_p(os.path.join(path, str(index)))
    write_data(os.path.join(path, str(index)), 'data.json', struct.dumps())
    write_data(os.path.join(path, str(index)), 'geometry.in', struct.get_geometry_atom_format())
    return index
        
def update_supercollection(structure_supercoll):
    '''
    updates and initializes entire supercollection, reading all database data in to local memory
    if filesystem collection is not in supercollection, it is added
    '''
    key_list = []
    # read all stoic in filesystem
    list_of_stoics = [name for name in os.listdir(structure_dir) \
                      if os.path.isdir(os.path.join(structure_dir, name))]
    for stoic_string in list_of_stoics:
        # read all cascade levels in filesystem
        list_of_inputs = [name for name in os.listdir(os.path.join(structure_dir, stoic_string)) \
                          if os.path.isdir(os.path.join(structure_dir, stoic_string, name))]
        for input_ref in list_of_inputs:
            # make a list of structure collection keys to update
            key_list.append((string_to_stoic(stoic_string), int(input_ref)))
        
    # update each collection or create new structure_collection in supercollection and update that  
    for key in key_list:
        if not key in structure_supercoll:
            structure_supercoll[key] = StructureCollection(key[0], key[1])
        structure_supercoll[key].update_local()
        
def string_to_stoic(stoic_string):
    '''
    Takes a string in the form Mg:2_O:5 and returns the StoicDict representation
    '''
    stoic_dict = StoicDict(int)
    for item in stoic_string.split('_'):
        [key, val] = item.split(':')
        stoic_dict[key] = int(val)
    return stoic_dict

if __name__ == '__main__':
    pass
