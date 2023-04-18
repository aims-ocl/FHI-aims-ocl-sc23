#!/usr/bin/python

'''
Created on Jul 19, 2013

@author: newhouse
'''
from __future__ import division

import datetime
import os
import sys
import time
# add source directory to python path
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(src_dir)

from core import user_input, data_tools, output
from file_handler import *
from kill import *
from structures import structure_collection
from structures.structure import Structure
from structures.structure_collection import StructureCollection, string_to_stoic
from utilities.stoic_model import determine_stoic



# from external_libs import asizeof


# from profilestats import profile  # for time profiling
# @profile
def main(stoic):
    # setup
    mkdir_p(tmp_dir)
    mkdir_p(structure_dir)
    set_unkill()
    # check if going to be multiprocess
    ui = user_input.get_config()
    n_processors = ui.get_eval('run_settings', 'parallel_on_core')
    if not isinstance(n_processors, (int, long)) : begin(stoic)
    else: begin_multiprocess(stoic, int(n_processors))

def begin_multiprocess(stoic_string, n_processors):
    from multiprocessing import Process
    processes_to_run = []
    for i in range(n_processors):  # @UnusedVariable
        time.sleep(0.1)  # slight separation between replicas
        p = Process(target=begin, args=(stoic,))
        processes_to_run.append(p)
        p.start()

def begin(stoic):
    replica = get_random_index()
    # run genetic algorithm
    ga = RunGA(replica, stoic)
    # catch crashes
    if ga.ui.get_eval('run_settings', 'recover_from_crashes') is not True: 
        ga.start()  # no exception catching
    else:
        while True:
            try:
                ga.start() 
                break
            except Exception, e: 
                print output.error(e, replica)
    output.move_to_shared_output(replica)
    
class RunGA():
    '''
    DOCUMENT
    '''
    def __init__(self, replica, stoic):
        '''
        Document
        '''

        # Define class fields
        self.replica = replica
        self.ui = user_input.get_config()
        self.replica_stoic = stoic
        self.structure_supercoll = {}
        self.number_of_structures = int(self.ui.get('run_settings', 'number_of_structures'))
        self.working_dir = os.path.join(tmp_dir, str(self.replica))
        self.verbose = self.ui.get_eval('run_settings', 'verbose')        
        # read in all other structures to local memory from storage
        self.control_list = self.ui.get_list('control', 'control_in_filelist')
        self.max_cascade = len(self.control_list) - 1
        for i in range(INITIAL_POOL_REFID, len(self.control_list)):  # initializes the structure collections necessary for cascade
            self.structure_supercoll[(self.replica_stoic, i)] = StructureCollection(self.replica_stoic, i)
        structure_collection.update_supercollection(self.structure_supercoll)
        self.ip_coll = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID))
        self.structure_coll = self.structure_supercoll.get((self.replica_stoic, self.max_cascade))
    
    def start(self):
        '''
        Performs main genetic algorithm operations
        Loads necessary modules based on UI at runtime
        '''  # TODO: improve documentation
        # Dynamically load modules outside of loop
        initial_pool_module = my_import(self.ui.get('modules', 'initial_pool_module'), package='initial_pool')
        selection_module = my_import(self.ui.get('modules', 'selection_module'), package='selection')
        crossover_module = my_import(self.ui.get('modules', 'crossover_module'), package='crossover')
        mutation_module = my_import(self.ui.get('modules', 'mutation_module'), package='mutation')
        relaxation_module = my_import(self.ui.get('modules', 'relaxation_module'), package='relaxation')
        comparison_module = my_import(self.ui.get('modules', 'comparison_module'), package='comparison')

        while True:
            ########## Beginning of Iteration Tasks ##########
            output.move_to_shared_output(self.replica)
            if self.verbose: self.output('Beginning iteration')
            self.ui = user_input.get_config()
            begin_time = datetime.datetime.now()
            
            ########## Check if finished ##########
            try: 
                if len(self.structure_coll.structures) - 1 >= self.number_of_structures: return
            except: pass
            if get_kill() == "kill": 
                print 'replica ' + str(self.replica) + ' killed'
                self.output('replica ' + str(self.replica) + ' killed')
                return
            
            ########## Fill initial pool ##########
            if self.verbose: self.output('Scanning initial pool')
            initial_pool_module.main(self.replica, self.replica_stoic)
            # update all
            structure_collection.update_supercollection(self.structure_supercoll)
            
            ########## Check for Unrelaxed Initial Pool Structures ##########
            unrelaxed_ip_struct = self.ip_coll.get_unrelaxed_structure()
            if unrelaxed_ip_struct is not False:
                if self.verbose: self.output('Found initial pool structure to be added to main pool')
                # create a new structure with the same geometry, this strips the unrelaxed structure of properties
                new_struct = Structure()
                new_struct.set_lattice_vectors(unrelaxed_ip_struct.get_lattice_vectors())
                new_struct.build_geo_whole(unrelaxed_ip_struct.get_geometry())
                # for keeping a record of parents
                structures_to_cross = [unrelaxed_ip_struct] 
            else:
            ########## Perform Genetic Algorithm Tasks ##########
            
                ########## Structure Selection ##########
                # Expects: dictionary_of_on_or_more<Stoic, StructureCollection> #Returns: list_of_2_or_more<Structure>
                if self.verbose: self.output('Beginning structure selection')
                structures_to_cross = selection_module.main(self.structure_supercoll, self.replica_stoic, self.replica)
                if structures_to_cross is False: self.output('Selection failure'); continue
                if self.verbose: 
                    stc_message = ''
                    for stc in structures_to_cross: stc_message += stc.get_path() + ' + '
                    self.output('Crossing structures: ' + stc_message[:-3])
                    
                ########## Mutation Decision ##########
                # Expects: None #Returns: stoichiometry or None
                # target_stoic in StoicDict format
                if not hasattr(mutation_module, 'get_targ_stoic'): targ_stoic = self.replica_stoic
                else: targ_stoic = mutation_module.get_targ_stoic()
                
                ########## Crossover ##########
                # Expects: list_of_2_or_more<Structure>, target stochiometry #Returns Structure or False
                if self.verbose: self.output('Beginning crossover')
                new_struct = crossover_module.main(structures_to_cross, targ_stoic, self.replica)
                if new_struct is False: 
                    if self.verbose: self.output('Crossover failure')
                    continue  # crossover failed, start with new selection
        
                ########## Mutation Execution ##########
                # Expects: Structure, target_stoichiometry [decision] #Returns: Structure
                new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
                if new_struct is False: 
                    if self.verbose: self.output('Mutation failure')
                    continue  # mutation failed, start with new selection
    
            ########## Begin Cascade ##########
            cascade_counter = 0
            if self.max_cascade > 0: self.output('Beginning cascade')
            ########## Comparison ##########
            # Expects: Structure, Structure:collection #Returns: Boolean
            # Checks if structure is unique or does not meet constraints (energy, etc.)
            structure_collection.update_supercollection(self.structure_supercoll)
            is_acceptable = comparison_module.main(new_struct, self.structure_supercoll.get((self.replica_stoic, cascade_counter)), self.replica)
            if is_acceptable is False: 
                if self.verbose: self.output('Structure to be relaxed is not acceptable')
                continue  # structure not acceptable start with new selection
    
            structures_to_add = {}
            while True:
                input_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                         self.control_list[cascade_counter])
                ########## Relaxation ##########
                # Expects: Structure, working_dir, input_path #Returns: Structure
                self.output('Beginning relaxation')
                struct = relaxation_module.main(new_struct, self.working_dir, input_string, self.replica)
                if struct is False: 
                    if self.verbose: self.output('Relaxation failure')
                    break  # optimization failed, start with new selection
                    
                ########## Comparison ##########
                structure_collection.update_supercollection(self.structure_supercoll)
                is_acceptable = comparison_module.main(struct, self.structure_supercoll.get((struct.get_stoic(), cascade_counter)), self.replica)
                if is_acceptable is False:
                    if self.verbose: self.output('Newly relaxed structure is not acceptable') 
                    break  # structure not acceptable start with new selection

                # Add structure to a list of structures to be put in to collection
                struct.input_ref = cascade_counter
                self.set_parents(structures_to_cross, struct)
                structures_to_add[(struct.get_stoic(), cascade_counter)] = struct

                # End of cascade tasks
                cascade_counter = cascade_counter + 1  # move to next input
                new_struct = struct  # set new structure and return to beginning of cascade
                if cascade_counter > self.max_cascade: break
                
            if is_acceptable is False or struct is False: continue  # start with new selection
            else: self.end_of_iteration_tasks(structures_to_add)
            end_time = datetime.datetime.now()
            if self.verbose: self.output("Iteration time: -- " + str(end_time - begin_time))

    def end_of_iteration_tasks(self, structures_to_add):
        prev_struct_index = None
        for key, struct in structures_to_add.items():
            # adds to shared pool if structure is acceptable
            struct.set_property('prev_struct_id', prev_struct_index)  # tracks cascade sequence
            index = structure_collection.add_structure(struct, key[0], key[1])
            prev_struct_index = str(key) + str(index)
            message = 'Success: \n  stoichiometry-- ' + key[0].get_string() + \
                      '\n  cascade-- ' + str(key[1]) + \
                      '\n  structure index-- ' + str(index) + \
                      '\n  replica-- ' + str(self.replica)
            print message
            data_tools.write_energy_hierarchy(self.structure_coll)
            self.output(message)
#             print str(asizeof.asizeof(self.structure_supercoll) / 1048576) + ' MB'
            
    def set_parents(self, structures_to_cross, struct):
        # set the the parent structures as a property for family tree tracing
        for i in range(len(structures_to_cross)): 
            par_st = structures_to_cross[i]
            struct.set_property('parent_' + str(i), par_st.get_stoic_str() + '/' \
                                + str(par_st.get_input_ref()) + '/' + str(par_st.get_struct_id()))

    def output(self, message): output.local_message(message, self.replica)


if __name__ == '__main__':
    '''
    This command is important. If the module is run directly instead of imported,
    it will execute the main() method. This allows for a single replica to join 
    a current genetic algorithm search.
    '''
    try: 
        stoic_filename = sys.argv[1]
        stoic = determine_stoic(os.path.join(cwd, stoic_filename))
    except: 
        stoic = determine_stoic()
    if stoic == None: raise Exception
        
    main(stoic)

