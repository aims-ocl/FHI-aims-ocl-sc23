'''
Created on Dec 16, 2013

@author: newhouse
'''
import math

from core import user_input
from element_radii import radii
from structures.structure import Structure


def is_acceptable(struct):
    # TODO: This MUST be made more efficent. a huge point of lag here
    ui = user_input.get_config()
    closeness_ratio = ui.get_eval('periodic', 'min_distance_ratio')
    superstructure = create_superstructure(struct)
    result = check_closeness(superstructure, closeness_ratio)
    return result

def addv(*argv):
    rv = [0, 0, 0]
    for vector in argv:
        for i in range(3):
            rv[i] += vector[i]
    return tuple(rv)

def create_superstructure(original_struct):
    lv = original_struct.get_lattice_vectors()
    if lv == False: raise Exception
    superstructure = Structure()
    # initial image
    augment_structure(superstructure, original_struct, (0, 0, 0))
    augment_structure(superstructure, original_struct, lv[0])
    augment_structure(superstructure, original_struct, lv[1])
    augment_structure(superstructure, original_struct, lv[2])
    augment_structure(superstructure, original_struct, addv(lv[0], lv[1]))
    augment_structure(superstructure, original_struct, addv(lv[1], lv[2]))
    augment_structure(superstructure, original_struct, addv(lv[2], lv[0]))
    augment_structure(superstructure, original_struct, addv(lv[0], lv[1], lv[2]))
    return superstructure

def augment_structure(structure_to_augment, structure_to_add, vector):
    for atom in structure_to_add.get_geometry():
        structure_to_augment.build_geo_by_atom(atom['x'] + vector[0], atom['y'] + vector[1], atom['z'] + vector[2], atom['element'])
    
def check_closeness(structure, closeness_ratio):
    c_a = 0
    for atom_a in structure.get_geometry():
        c_b = 0
        rad_a = radii.get(atom_a['element'])
        for atom_b in structure.get_geometry():
            if c_a == c_b: c_b += 1; continue
            c_b += 1
            rad_b = radii.get(atom_b['element'])
            dist = math.sqrt((atom_a['x'] - atom_b['x']) ** 2 + 
                             (atom_a['y'] - atom_b['y']) ** 2 + 
                             (atom_a['z'] - atom_b['z']) ** 2)
            threshold = closeness_ratio * (rad_a + rad_b)
            if dist < threshold: return False
        c_a += 1
    return True
