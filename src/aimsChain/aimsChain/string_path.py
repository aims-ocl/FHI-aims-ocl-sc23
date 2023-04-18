"""
This module defines the String Path
"""

import numpy as np
from aimsChain.utility import vmag, vunit, vproj
from aimsChain.path import Path
from aimsChain.optimizer.optimize import FDOptimize

class StringPath(Path):

    def __init__(self, 
             nodes = [],
             control = None):
        super(StringPath,self).__init__(nodes,control)
            

    def move_nodes(self):
        """
        Move the node according to normal force
        the original string method.
        """

        from aimsChain.interpolate import spline_pos
        import os

        positions = []
        forces = []
        new_t = []
        new_pos = []

        for node in self.nodes:
            # List all the new params:
            new_t.append(node.param)
            positions.append(node.positions)
            forces.append(node.normal_forces)

        forces = np.array(forces)
        positions = np.array(positions)

        if self.control.global_opt:
            save = os.path.join(self.nodes[0].dir_pre, "path.opt")
            new_pos,opt = self.g_opt(
                self.control.optimizer,
                positions,
                forces,
                save)
        else:
            new_pos,opt = self.nong_opt(
                self.control.optimizer,
                self.nodes,
                positions,
                forces,
                ".opt")

        if not (isinstance(opt, FDOptimize) and  opt.finite_diff):
            new_pos = spline_pos(new_pos, new_t)

        for i,position in enumerate(new_pos):
            self.nodes[i].positions = position

        forces = np.reshape(forces, (-1,3))
        forces = np.sum(forces**2,1)**0.5

        return np.nanmax(forces)
            
