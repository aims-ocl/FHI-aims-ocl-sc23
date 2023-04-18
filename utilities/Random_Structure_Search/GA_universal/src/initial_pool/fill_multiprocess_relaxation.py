'''
Created on Aug 13, 2013

@author: newhouse
'''
from multiprocessing import Process, Queue
import os
import time

from core import user_input
from core.file_handler import cwd, set_progress, my_import, tmp_dir, read_data
from external_libs.filelock import FileLock
from structures import structure_collection
from structures.structure import get_geo_from_file, Structure
from structures.structure_collection import StructureCollection


def main(replica, stoic):
    # create a structure collection for initial pool generation. -1 specifies initial pool
    fip = FillInitialPool(replica, stoic)
    
    if len(fip.structure_coll.structures) < fip.num_initial_pool_structures: 
        fip.add_random_structures()
    
    # check is done internally
    fip.copy_user_structures()

    set_progress('initial_pool_filled')

class FillInitialPool():
    
    def __init__(self, replica, stoic):
        self.replica = replica
        self.ui = user_input.get_config()
        self.INITIAL_POOL_REFID = -1
        self.structure_coll = StructureCollection(stoic, self.INITIAL_POOL_REFID)
        ########## Random Structure Generation ##########
        self.random_gen_module = my_import(self.ui.get('modules', 'random_gen_module'), package='random_gen')
        ########## Relaxation ##########
        self.initial_pool_relaxation_module = my_import(self.ui.get('modules', 'initial_pool_relaxation_module'), package='relaxation')
        ########## Comparison ##########
        self.comparison_module = my_import(self.ui.get('modules', 'initial_pool_comparison_module'), package='comparison')
        
        self.stoic = stoic
        self.num_initial_pool_structures = int(self.ui.get('initial_pool', 'num_initial_pool_structures'))
        self.control_in = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                    self.ui.get('control', 'initial_pool'))
        self.verbose = self.ui.get_eval('run_settings', 'verbose')

    def add_random_structures(self):
    
        # read number of processoors available
        self.batch_size = int(self.ui.get('initial_pool', 'number_of_processors'))
        batch_counter = 0
                
        while True:
            # Records the structure_coll. Shared across each process.
            q = Queue() 
            # keeps track of the running processes
            processes_to_run = []
            # make structure_coll on as many threads as there are available processors
            for i in range(self.batch_size):  # @UnusedVariable
                p = Process(target=self.make_initial_structure, args=(batch_counter, q))
                processes_to_run.append(p)
                p.start()
                batch_counter += 1
            # check if all threads are complete
            while True:
                all_done = True
                for p in processes_to_run:
                    p.join()
                    if p.is_alive():
                        all_done = False 
                if all_done: break
                else:
                    time.sleep(0.01) 
            # add current batch to structure_coll
            initial_pool_full = False
            while not q.empty():  # check if required random structure_coll are fulfilled
                self.structure_coll.update_local()
                if len(self.structure_coll.structures) >= self.num_initial_pool_structures:
                    initial_pool_full = True 
                    break
                # add new structure from queue (if unique)
                else:
                    # get next structure from the queue and check 
                    struct = q.get()
                    if not struct.get_stoic() == self.stoic: continue
                    if self.comparison_module.main(struct, self.structure_coll, self.replica):
                        structure_collection.add_structure(struct, struct.get_stoic(), self.INITIAL_POOL_REFID)                        
            if initial_pool_full: 
                break
        print 'initial pool filled'

    def make_initial_structure(self, counter, q):
        '''
        Called by a multiprocessing thread. Creates a random structure and optimizes it.
        '''
        seed = time.time() + 17 * counter + int(self.replica[:-4], 16)
        rand_struct = self.random_gen_module.main(self.stoic, seed)
    
        working_dir = os.path.join(tmp_dir, 'initial_pool', str(self.replica), str(counter % self.batch_size)) 
        struct = self.initial_pool_relaxation_module.main(rand_struct, working_dir, self.control_in, self.replica)

        # add structure to shared queue
        if not struct == False: 
            struct.set_property('to_be_relaxed', True)
            struct.input_ref = self.INITIAL_POOL_REFID
            q.put(struct)
        return
            
    def copy_user_structures(self):
        '''
        Adds structures defined by the users to the shared storage.
        Only adds if filename is not already present in shared storage.
        Can be called at every iteration to check for new files
        '''
        # setup
        user_structures_dir = os.path.join(cwd, self.ui.get('initial_pool', 'user_structures_dir'))
        added_user_structures = os.path.join(tmp_dir, 'added_user_structures.dat')
        open(added_user_structures, 'a').close()  # creates if non-existant
        # read files
        try:files = os.listdir(user_structures_dir)
        except Exception: return  # no directory
        if len(files) == 0:
            return
        # copy files and add name to list of copied files to avoid duplicates
        files_to_add = []
        with FileLock(added_user_structures):
            ffile = open(added_user_structures, 'r')
            file_list = ffile.read()
            ffile.close()
            ffile = open(added_user_structures, 'a')
            for item in files:
                filename = os.path.join(user_structures_dir, item)
                if not filename in file_list:
                    ffile.write(filename)
                    files_to_add.append(filename)
            ffile.close()
            
        for file in files_to_add:
            # Tests passed. add structure to collection.
            try:
                struct = Structure()
                struct.build_geo_from_atom_file(file)
            except: continue  # could not build structure from geometry
            # keeps track of which user defined structures have been added
            struct.set_property('to_be_relaxed', True)
            struct.set_property('u_def_filename', filename)
            structure_collection.add_structure(struct, struct.get_stoic(), self.INITIAL_POOL_REFID)
            print 'user initial pool structure added'

if __name__ == '__main__':
    main()
