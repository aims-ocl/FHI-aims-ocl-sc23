'''
Created on Jul 29, 2013

@author: newhouse
'''
import os
import subprocess

from core import user_input, output
from core.file_handler import mkdir_p_clean
from structures.structure import Structure
from copy import deepcopy


def main(input_structure, working_dir, control_in_string, replica):
    '''
    For a relaxation module, the method main() must be implemented.
    
    It must take as arguments a 'Structure', a string which defines 
    in which working directory the relaxation will take place, and a string
    representation of the full control.in file
    
    It must return a Structure object. (Please see structure documentation
    for information on how to easily setup a Structure object) 
    '''
    r = FHIAimsRelaxation(input_structure, working_dir, control_in_string, replica)
    r.setup()
    r.execute()
    if r.is_successful(): return r.extract()
    else: return False

class FHIAimsRelaxation():
    '''
    Creates an optimizer object using FHI-Aims to execute a geometry.
    
    This is an example of using a class to define a state, in this case the
    state of the relaxation process throughout setup, execution, and extraction
    '''
    def __init__(self , input_structure, working_dir, control_in_string, replica):
        '''
        __init__ runs when class is initialized. 
        '''
        # this gives easy access to retrieve preferences from the user_input file
        self.ui = user_input.get_config()
        self.replica = replica
        # defining 'self' fields gives access to these values throughout the class 
        self.input_structure = input_structure
        self.working_dir = working_dir
        self.control_in_string = control_in_string
    
    def output(self, message): output.local_message(message, self.replica)    

    def setup(self):
        '''
        Creates the environment for the optimizer and constructs the input files.
        In this case, FHI-aims will take a geometry.in file and an comtrol.in file.
        '''
        mkdir_p_clean(self.working_dir)
        self.create_geometry_in()
        self.create_control_in()
        
        
    def create_geometry_in(self):
        '''
        Writes the geometry.in file called by FHI-aims. 
        Calls from FHI-aims control settings found in res directory.
        Returns: None
        '''
        geometry_file = open(os.path.join(self.working_dir, 'geometry.in'), 'w')
        geometry_file.write(self.input_structure.get_geometry_atom_format())
        geometry_file.close()
        
    def create_control_in(self):
        '''
        Writes the control.in input file required by FHI-aims.
        Returns: None
        '''
        control_file = open(os.path.join(self.working_dir, 'control.in'), 'w')
        control_file.write(self.control_in_string)
        control_file.close()
        
    def execute(self):
        '''
        Creates an executable shell script with proper FHI-aims syntax 
        then executes in a subprocess
        '''
        # creates shell script with proper syntax and runs it as a subprocess
        exe_string = ''
        exe_string += '#!/bin/bash\n'
        # load required libraries
        exe_string += 'source /opt/intel/composer_xe_2013.1.117/bin/ifortvars.sh intel64\n'
        exe_string += 'source /opt/intel/composer_xe_2013.1.117/mkl/bin/mklvars.sh intel64\n'
        exe_string += 'ulimit -s unlimited\n'
        exe_string += 'ulimit -c unlimited\n'
        exe_string += 'export PATH=/usr/local/bin:$PATH\n'
        exe_string += 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH\n'
        # the mpi command
        exe_string += 'mpirun -np ' + self.ui.get('FHI-aims', 'number_of_processors') + ' ' + self.ui.get('FHI-aims', 'path_to_aims_executable') + ' > ' + 'aims.out'
        
        exe_script = open(os.path.join(self.working_dir, 'exe.sh'), 'w')
        exe_script.write(exe_string)
        exe_script.close()
        
        p = subprocess.Popen(['/bin/bash', 'exe.sh'], cwd=self.working_dir)
        p.wait()    

    def extract(self):
        '''
        Reads the output files from relaxation and specifies properties of the structure (e.g. energy)
        
        Returns: Structure or False
        
        Note: there is a quicker way to parse these files, but is left as is for readabilities sake
        '''
        # creates empty structure
        self.result_struct = Structure()
        
        try: self.result_struct.set_lattice_vectors(self.input_structure.get_lattice_vectors())
        except: pass
        
        # reads output file to extract energy
        energy = self.extract_energy()
        if energy == False: return False  # did not converge
        self.result_struct.set_property('energy', energy)
        
        # sets the result structure's geometry based on output file
        extract_success = self.extract_geometry()
        # structure not relaxed, no final atomic structure, keep original geometry
        if not extract_success: self.result_struct.build_geo_whole(deepcopy(self.input_structure.get_geometry()))
        if not self.input_structure.get_stoic() == self.result_struct.get_stoic(): return False
        
        # extract the spin moment
        if self.ui.get('FHI-aims','initial_moment') == 'custom': self.extract_spin_moment()
        
        return self.result_struct

    def extract_energy(self):
        '''
        reads the resulting energy from the FHI-aims log
        Returns: float if successful, False if unsuccessful
        '''
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        while True:
            line = aims_out.readline()
            if not line: return False  # energy not converged
            if '  | Total energy corrected        :' in line:
                tokens = line.split()
                energy = float(tokens[5])  # converts from SI string to float
                return energy
        

    def extract_geometry(self):
        '''
        Reads the FHI-aims output file and builds the structure's geometry
        Returns: True if successful
        '''
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        while True:
            line = aims_out.readline()
            if not line: return False  # not converged
            if 'Final atomic structure' in line: break
        aims_out.readline()  # skip one line
        
        atom_string = ''
        while True:
            line = aims_out.readline()
            if '----' in line: break
            else: atom_string += line
            
        self.result_struct.build_geo_whole_atom_format(atom_string)
        return True
    
    def extract_spin_moment(self):
        try:
            os.system("grep '|   Hirshfeld spin moment' aims.out > spin_moment.out")
            spin_file = open('spin_moment.out', 'r')
            total_spin_list = spin_file.read()
            spin_file.close()
        except: return False
        
        reduced_spin_list = total_spin_list[-len(self.result_struct.get_geometry()):]
        for i in range(len(reduced_spin_list)):
            self.result_struct.geometry[i]['spin'] = float(reduced_spin_list[i].split()[-1])
        return True
        
    def extract_charge(self):
        try:
            os.system("grep '|   Hirshfeld charge' aims.out > charge.out")
            charge_file = open('charge.out', 'r')
            total_charge_list = charge_file.read()
            charge_file.close()
        except: return False
        
        reduced_charge_list = total_charge_list[-len(self.result_struct.get_geometry()):]
        for i in range(len(reduced_charge_list)):
            self.result_struct.geometry[i]['charge'] = float(reduced_charge_list[i].split()[-1])
        return True
    
    def is_successful(self):
        '''
        checks if relaxation/optimization was successful
        '''
        # TODO: need proper check
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        while True:
            line = aims_out.readline()
            if not line: return False  # energy not converged
            if 'Have a nice day' in line:
                return True
        return False
    
    
        