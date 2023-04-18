'''
Created on Jul 29, 2013

@author: newhouse
'''
import os

from utilities import element_masses
from relaxation import reax_relaxation

def main(input_structure, working_dir, control_in_string, replica):
    '''
    For a relaxation module, the method main() must be implemented.
    
    It must take as arguments a 'Structure' and an 'int' which defines 
    in which working directory the relaxation will take place.
    
    It must return a Structure object. (Please see structure documentation
    for information on how to easily setup a Structure object) 
    '''
    try:
        r = LJRelaxation(input_structure, working_dir, control_in_string, replica)
        r.setup()
        r.execute()
        return r.extract()
    except Exception, e: print str(e); return False # if there is an error in lammps this will prevent crashing
    pass

class LJRelaxation(reax_relaxation.ReaxRelaxation):
    '''
    This example takes advantage of the similarities between reax and lj relaxations.
    The LJRelaxation class extends the ReaxRelaxation class 
    and now has access to all of its methods
    Those that do not apply to LJ can simply be overridden in this class
    
    Strictly speaking, this is an abuse of subclassing. 
    One can not reasonably say "LJRelaxation is-a ReaxRelaxation"
    However it is beneficial for this case of module separability.
    '''
        
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
            if self.element_list[i] in element_masses.masses and element_masses.masses[self.element_list[i]]:
                data_string += str(i + 1) + ' ' + str(element_masses.masses[self.element_list[i]]) + '\n';
        # write element coordinates
        # line example:  1 1 1.0 2.0 3.0
        # definitions:   index element_index x y z
        data_string += '\nAtoms\n\n'
        counter = 1
        for i in range(len(self.element_list)):
            for j in range(geo.size):
                if geo[j]['element'] == self.element_list[i]:
                    data_string += str(counter) + ' ' + str(i + 1) + ' ' # DIFFERENCE 
                    data_string += str(geo[j]['x'])\
                     + ' ' + str(geo[j]['y']) + ' ' + str(geo[j]['z'])
                    data_string += '\n'
                    counter += 1
        # write file
        data_file = open(os.path.join(self.working_dir, 'geo.dat'), 'w')
        data_file.write(data_string)
        data_file.close()
