'''
Created on Jul 29, 2013

@author: newhouse
'''
import os
import re
from shutil import copyfile
import shutil
import subprocess
import time

from core import user_input, output
from core.file_handler import mkdir_p_clean, res_dir
from structures.structure import Structure
from utilities.element_masses import masses


def main(input_structure, working_dir, control_in_string, replica):
    '''
    For a relaxation module, the method main() must be implemented.
    
    It must take as arguments a 'Structure' and an 'int' which defines 
    in which working directory the relaxation will take place.
    
    It must return a Structure object. (Please see structure documentation
    for information on how to easily setup a Structure object) 
    '''
    try: 
        r = ReaxRelaxation(input_structure, working_dir, control_in_string, replica)
        r.setup()
        r.execute()
        return r.extract()
    except Exception, e: print str(e); return False # if there is an error in lammps this will prevent crashing
    pass

class ReaxRelaxation():
    '''
    Creates an relaxation object using lammps to execute a geometry.
    
    This is an example of using a class to define a state, in this case the
    state of the relaxation process throughout setup, execution, and extraction
    '''
    def __init__(self, input_structure, working_dir, control_in_string, replica):
        '''
        __init__ will always run when a class is initialized. 
        This is also how the class can take arguments.
        '''
        # defining 'self' fields gives access to these values throughout the class 
        self.replica = replica
        self.working_dir = working_dir
        self.input_structure = input_structure
        self.control_in_string = control_in_string
        # this gives easy access to retrieve preferences from the user_input file
        self.ui = user_input.get_config()

    def output(self, message): output.local_message(message, self.replica)

    def setup(self):
        '''
        Creates the environment for the optimizer and constructs the input files.
        In this case, LAMMPS will take a geometry.dat file and an in.dat file.
        '''
        mkdir_p_clean(self.working_dir)
        self.set_element_list()
        self.create_geometry_in()
        self.create_input()
        shutil.copy(os.path.join(res_dir, 'ffield', 'reax.ffield'), os.path.join(self.working_dir, 'ffield'))
        
    def execute(self):
        '''
        Creates an executable shell script with proper LAMMPS syntax 
        then executes in a subprocess
        
        Subprocess is very powerful and useful in this context
        The only reason a shell script is used in the intermediate is because of
        LAMMPS' "read in" style with "<"
        '''
        # creates shell script with proper syntax and runs it as a subprocess
        exe_file = open(os.path.join(self.working_dir, 'exe.sh'), 'w')
        exe_file.write(self.ui.get('lammps', 'path_to_lammps_executable') + ' < in.dat' + '  >' + os.path.join(self.working_dir, 'log'))  # write log        
        exe_file.close()

        p = subprocess.Popen(['sh', 'exe.sh'], cwd=self.working_dir)
        p.wait()    
        
    def extract(self):
        '''
        Reads the output files from relaxation and specifies properties of the structure (e.g. energy)
        
        Returns: Structure or False
        '''
        # creates empty structure
        self.result_struct = Structure()
        
        # reads output file to extract energy
        energy = self.extract_energy()
        if energy == False: return False  # did not converge
        self.result_struct.set_property('energy', energy)
        
        # sets the result structure's geometry based on output file
        extract_success = self.extract_geometry()
        if not extract_success: return False
        if not self.input_structure.get_stoic() == self.result_struct.get_stoic(): return False
        
        return self.result_struct
        
    '''
    The following are helper methods that do the actual work of setup and extraction.
    Helper methods are not required by the interface but may be implemented by the developer
    to make the module more organized and readable.
    '''
            
    def set_element_list(self):
        '''
        Keeps track of the elements in the structure.
        This is because LAMMPS uses indices instead of strings to specify elements.
        '''
        element_list = []
        for i in range(self.input_structure.get_n_atoms()):
            if self.input_structure.geometry[i]['element'] not in element_list:
                element_list.append(self.input_structure.geometry[i]['element'])
        element_list.sort()
        self.element_list = element_list
        
    def create_geometry_in(self):
        '''
        Writes the data (geometry) file called by LAMMPS. Builds with no template.
        Format is very particular.
        Returns: None
        '''
        BOUNDARY = '-100.000 100.000'  # decimals mater #TODO: make this vary with cluster size
        geo = self.input_structure.geometry
        
        data_string = ''
        # set up specifications
        data_string += '# '
        for i in range(len(self.element_list)):
            data_string += self.element_list[i] + ' ';
        data_string += '\n\n'
        data_string += str(geo.size) + ' atoms\n'
        data_string += str(len(self.element_list)) + ' atom types\n\n'
        data_string += BOUNDARY + ' xlo xhi\n'
        data_string += BOUNDARY + ' ylo yhi\n'
        data_string += BOUNDARY + ' zlo zhi\n'
        # search for element masses
        data_string += '\nMasses\n\n'
        for i in range(len(self.element_list)):
            if self.element_list[i] in masses and masses[self.element_list[i]]:
                data_string += str(i + 1) + ' ' + str(masses[self.element_list[i]]) + '\n';
        # write element coordinates
        # line example:  1 1 0.0 -0.290828 -0.177916 0.0508826
        # definitions:   index element_index charge x y z
        data_string += '\nAtoms\n\n'
        counter = 1
        for i in range(len(self.element_list)):
            for j in range(geo.size):
                if geo[j]['element'] == self.element_list[i]:
                    data_string += '    ' + str(counter) + ' ' + str(i + 1) + ' 0.0 ' 
                    data_string += str(geo[j]['x'])\
                     + ' ' + str(geo[j]['y']) + ' ' + str(geo[j]['z'])
                    data_string += '\n'
                    counter += 1
        # write file
        data_file = open(os.path.join(self.working_dir, 'geo.dat'), 'w')
        data_file.write(data_string)
        data_file.close()
        
    def create_input(self):
        '''
        Writes the input file required by LAMMPS.
        Returns: None
        '''
        input_file = open(os.path.join(self.working_dir, 'in.dat'), 'w')
        input_file.write(self.control_in_string)
        input_file.close()

    def extract_geometry(self):
        '''
        Reads the dump file in xyz format and builds the structure's geometry
        Returns: True if successful
        '''
        # reads file
        dump = open(os.path.join(self.working_dir, 'dump.xyz'), 'r')
        dump_list = []
        for line in dump:
            dump_list.append(line)
        # reverses file to read from bottom up
        dump_list.reverse()
        for i in range(len(self.input_structure.geometry)):
            # example line:  1 -0.290828 -0.177916 0.0508826
            # definitions:   element_index x_coord y_coord z_coord
            tokens = dump_list[i].split()
            self.result_struct.build_geo_by_atom(float(tokens[1]),
                                                 float(tokens[2]),
                                                 float(tokens[3]),
                                                 self.element_list[int(tokens[0]) - 1])
             
        # sort for the sake of cosistancy. not very necessary.
        self.result_struct.geometry.sort(order='element')
        return True
    
    def extract_energy(self):
        '''
        reads the resulting energy from the LAMMPS log
        Returns: float if successful, False if unsuccessful
        '''
        loglammps = open(os.path.join(self.working_dir, 'log.lammps'), 'r')

        while True:
            line = loglammps.readline()
            if 'Energy initial, next-to-last, final =' in str(line):
                energies = loglammps.readline().split()
                return float(energies[2])  # final energy
            if not line: return False

