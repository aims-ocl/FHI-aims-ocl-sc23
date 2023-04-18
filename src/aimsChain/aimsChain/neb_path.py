"""
This module defines the NEB Path
"""

import numpy as np
from aimsChain.utility import vmag, vunit, vproj
from aimsChain.path import Path


class NebPath(Path):
    
    def __init__(self, 
             nodes = [],
             control = None):
        super(NebPath,self).__init__(nodes,control)


    def move_nodes(self):
        """ Move the nodes according to NEB rules."""

        import os

        positions = []
        forces = []
        new_pos = []

        for node in self.nodes:
            positions.append(node.positions)
            forces.append(node.spring_forces)      

        forces = np.array(forces)
        positions = np.array(positions)        

        if self.control.global_opt:
            save = os.path.join(self.nodes[0].dir_pre, "path.opt")
            new_pos = self.g_opt(
                self.control.optimizer,
                positions,
                forces,
                save)[0]
        else:
            new_pos = self.nong_opt(
                self.control.optimizer,
                self.nodes,
                positions,
                forces,
                ".opt")[0]

        for i,pos in enumerate(new_pos):
            self.nodes[i].positions = pos            

        forces = np.reshape(forces, (-1,3))
        forces = np.sum(forces**2,1)**0.5

        return np.nanmax(forces)
   
 
    def move_climb(self):
        """ Move the climbing nodes."""

        import os

        moving_nodes = []
        forces = []
        positions = []
        climb_mode = self.control.climb_mode

        if climb_mode == 3:
            self.find_climb(True)

        for node in self.nodes:
            if not node.fixed:
                moving_nodes.append(node)

        # Get all the forces and positions.
        for i,node in enumerate(moving_nodes):
            forces.append(node.climb_forces)
            positions.append(node.positions)

        forces = np.array(forces)
        positions = np.array(positions)

        # Move nodes, either by global or non-global optimizer.
        if self.control.climb_global_opt:
            save = os.path.join(moving_nodes[0].dir_pre, "climbing.opt")
            new_pos = self.g_opt(
                self.control.climb_optimizer,
                positions,
                forces,
                save)[0]
        else:
            new_pos = self.nong_opt(
                self.control.climb_optimizer,
                moving_nodes,
                positions,
                forces,
                ".climb.opt")[0]

        for i, position in enumerate(new_pos):
            moving_nodes[i].positions = position        

        forces = np.reshape(forces, (-1,3))
        forces = np.sum(forces**2,1)**0.5        

        return np.nanmax(forces)
