'''
Created on Aug 2, 2013

@author: newhouse
'''
from ast import literal_eval
import os
import sqlite3
import time

from core.file_handler import db_file
from structures.structure import Structure, StoicDict


class StructureCollection(object):
    '''
    This class handles the functionality for groups of structures. It is an aggregation of Structure
    '''    
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
    def get_stoic_str(self): return str(self.stoichiometry)
    def get_key(self): return (self.get_stoic_str(), int(self.input_ref))
    def get_conn(self):
        ''' makes a connection to the database'''
        conn = sqlite3.connect(db_file, detect_types=sqlite3.PARSE_DECLTYPES)
        conn.text_factory = str
        return conn    
#     def get_key(self):
#         return (self.stoic, self.input_ref)

    def update_local(self):
        '''
        Reads structures from database and updates local object list with new structures
        This is vital for parallelism when other replicas are populating the same list
        '''
        conn = self.get_conn()
        cursor = conn.cursor()
        # query structures of this stoic and control input
        while True:
            try: # in case database is locked
                cursor.execute('SELECT struct_id FROM structure \
                                WHERE stoic = ? AND input_id = ?',
                                (self.stoichiometry.get_string(), self.input_ref))
                database_indicies = cursor.fetchall()
                if database_indicies == None: raise Exception
            except Exception,e: print str(e); continue
            break
        cursor.close()
        conn.close()
        # for new indicies, add information to local memory
        for s_index in database_indicies: 
            if s_index[0] not in self.structures.keys():
                self.add_structure_to_local(s_index[0])

    def add_structure_to_local(self, s_index):
        conn = self.get_conn()
        cursor = conn.cursor()
        # get single structure info
        while True:
            try: 
                cursor.execute('SELECT * FROM structure \
                                WHERE struct_id = ?', (s_index,))
            except Exception,e: print str(e); continue
            break
        [struct_id, input_id, stoic_string, geo, prev_struct_id] = cursor.fetchall()[0]  # @UnusedVariable
        struct = Structure()
        # build structure
        struct.build_geo_whole(geo)
        struct.input_id = input_id
        struct.prev_struct_id = prev_struct_id
        struct.struct_id = struct_id
        struct.get_stoic() #calculates
        # get and set properties
        while True:
            try: 
                cursor.execute('SELECT * FROM property \
                                WHERE struct_id = ?', (s_index,))
            except Exception,e: print str(e); continue
            break
        properties = cursor.fetchall()
        
        cursor.close()
        conn.close()
        
        for item in properties:
            try: value = literal_eval(item[2]) 
            except: value = str(item[2])
            struct.set_property(item[1], value) 
        # add to structures dictionary
        self.structures[struct_id] = struct

    def add_structure(self, struct):
        '''
        iserts a structure and its properties into a sqlite database and returns the struct_id
        '''
        conn = self.get_conn()
        cursor = conn.cursor()
        
        # insert structure
        new_id = -1  # initialize
        while True:
            # select current max id number and add one
            new_id = cursor.execute('SELECT max(struct_id) FROM structure').fetchone()[0] + 1
            # attemts to insert. If structure with id already exists, error is returned and loop restarts
            try:
                cursor.execute('INSERT INTO structure (struct_id, input_id, stoic, geo) \
                                VALUES (?, ?, ?, ?)',
                                (new_id, self.input_ref, self.stoichiometry.get_string(), struct.get_geometry()))
                
                prop_list = []
                for prop in struct.properties.iteritems():
                    prop_list.append((new_id, prop[0], str(prop[1]),))
                cursor.executemany('INSERT INTO property (struct_id, key, value) \
                                VALUES (?, ?, ?)', prop_list)
                conn.commit()
            except Exception,e: print str(e); continue # structure id clash, re-evaluate new id
            break # if successful
        
        # insert attributes
        conn.commit()
    
        cursor.close()
        conn.close()
        
        if new_id > 0:  # no errors 
            struct.index = new_id
            print 'structure added to DB with ID: ' + str(new_id)
            return new_id
        else: raise Exception  # insertion problerm
        
class InitialPoolStructureCollection(StructureCollection):
    '''
    This class extends the StructureCollection class with additional functionality for 
    initial-pool specific operations
    
    This acts as a self-implemented type casting for StructureCollection
    '''
    def __init__(self, structure_coll):
        self.structure_coll = structure_coll

    def __getattr__(self, attr):
        return getattr(self.structure_coll, attr)
    
    def get_list_of_ui_filenames(self):
        conn = self.get_conn()
        cursor = conn.cursor()
        while True:
            try: 
                cursor.execute('SELECT value FROM property \
                                WHERE key = "u_def_filename" \
                                GROUP BY value')
            except Exception,e: print str(e); continue
            break
        data = cursor.fetchall()
        cursor.close()
        conn.close()
        filenames = [element[0] for element in data]
        return filenames
    
    def get_unrelaxed_structure(self):
        for struct in self.structures.itervalues():
            if struct.get_property('to_be_relaxed') == True: 
                self.mark_as_relaxed(struct)
                return struct
            else: continue
        return False  # none found
    
    def mark_as_relaxed(self, struct):
        conn = self.get_conn()
        cursor = conn.cursor()
        while True:
            try: 
                cursor.execute('UPDATE property\
                                SET value="False"\
                                WHERE struct_id = ? AND key = "to_be_relaxed"',
                                (struct.struct_id,))
            except Exception,e: print str(e); continue
            break
        conn.commit()
        cursor.close()
        conn.close()
        struct.set_property('to_be_relaxed', False)
        
    
def update_supercollection(structure_supercoll):
    '''
    updates and initializes entire supercollection, reading all database data in to local memory
    '''

    conn = sqlite3.connect(db_file, detect_types=sqlite3.PARSE_DECLTYPES)
    conn.text_factory = str  # use utf8 instead of unicode
    cursor = conn.cursor()
    # read stoichiometry/input_ref tuples from database
    while True:
            try: 
                cursor.execute('SELECT stoic, input_id FROM structure \
                                WHERE stoic IS NOT NULL \
                                GROUP BY stoic, input_id')
            except Exception,e: print str(e); continue
            break
    result = cursor.fetchall()
    cursor.close()
    conn.close()
    # update each collection or create new structure_collection in supercollection and update that  
    for stoic_string, input_ref in result:
        stoic = string_to_stoic(stoic_string)
        key = (stoic, input_ref)
        if not key in structure_supercoll:
            structure_supercoll[key] = StructureCollection(stoic, input_ref)
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
