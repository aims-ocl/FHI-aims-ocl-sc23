""" The growing string path.
"""

import numpy as np
from aimsChain.utility import vmag, vunit, vproj
from aimsChain.path import Path
from aimsChain.node import Node
from aimsChain.optimizer.optimize import FDOptimize
import cPickle as cp

class GrowingStringPath(Path):

    def __init__(self, 
             nodes = [],
             control = None):

        super(GrowingStringPath,self).__init__(nodes,control)

        self.lower_end = 0

    @property
    def upper_end(self):
        return self.lower_end + 1
    @property
    def linsep(self):
#        print "using " + str(self.control.nimage) + " images\n"
        return 1.0/(self.control.gs_nimage + 1)

    @property
    def params(self):
        t = []
        for node in self.nodes:
            t.append(node.param)
        return t
    @property
    def lower(self):
        return self.nodes[self.lower_end]
    @property
    def upper(self):
        return self.nodes[self.upper_end]


    def add_lower(self):

        from aimsChain.interpolate import spline_pos
        import copy
        import os

        if ((self.upper.param - self.lower.param)
            <= self.linsep * 1.5):
            return False

        self.insert_node(self.lower.param + self.linsep)
        self.lower_end += 1

#        if self.control.gs_global_optimizer:
#            save = os.path.join(self.nodes[0].dir_pre, "pgs.opt")
#            opt = self.get_optimizer("dampedBFGS", save)
#            opt.initialize()
 #           opt.load()
 #           opt.insert_node(self.lower_end, #loc to insert
 #                           len(self.nodes[0].geometry.atoms)) #num of atoms

        if ((self.upper.param - self.lower.param)
            <= self.linsep * 1.5):
            return False

        return True

        
    def add_upper(self):    

        import copy
        import os
        from aimsChain.interpolate import spline_pos

        if ((self.upper.param - self.lower.param)
            <= self.linsep * 1.5):
            return False

        self.insert_node(self.upper.param - self.linsep)        

        if ((self.upper.param - self.lower.param)
            <= self.linsep * 1.5):
            return False

 #       if self.control.gs_global_optimizer:
 #           save = os.path.join(self.nodes[0].dir_pre, "pgs.opt")
 #           opt = self.get_optimizer("dampedBFGS", save)
 #           opt.initialize()
 #           opt.load()
 #           opt.insert_node(self.lower_end, #loc to insert
 #                           len(self.nodes[0].geometry.atoms)) #num of atoms 
        return True


    def lower_tangent(self):

        from aimsChain.interpolate import spline_pos

        if (self.upper.param - self.lower.param) <= self.linsep*1.5:
            return self.lower.get_tangent()

        positions = []
        param = self.params
        loc = self.lower_end

        for i in self.nodes:
            positions.append(i.positions)          

        derv = spline_pos(positions, param, param, 3, 1)[loc]
        tangent = derv        

        return vunit(tangent)


    def upper_tangent(self):

        from aimsChain.interpolate import spline_pos

        if (self.upper.param - self.lower.param) <= self.linsep*1.5:
            return self.upper.get_tangent()

        positions = []
        param = self.params
        loc = self.upper_end

        for i in self.nodes:
            positions.append(i.positions)          

        derv = spline_pos(positions, param, param, 3, 1)[loc]
        tangent = derv        

        return vunit(tangent)


    def lower_force(self):

        forces = self.lower.forces
        tangent = self.lower_tangent()
        forces -= vproj(forces,tangent)

        return forces


    def upper_force(self):

        forces = self.upper.forces
        tangent = self.upper_tangent()
        forces -= vproj(forces,tangent)

        return forces


    def move_nodes(self):

        forces = None

        if (self.upper.param - self.lower.param) <= self.linsep*1.5:
            forces = self.move_nodes_combined()
        else:
            forces = self.move_nodes_separated()

        return forces


    def move_nodes_separated(self):

        from aimsChain.interpolate import spline_pos, get_total_length, get_t
        import os
        
        positions = []
        forces = []
        new_pos = []
        new_pos2= []
        total_length = 0
        low_force = self.lower_force()
        high_force = self.upper_force()

        for node in self.nodes:
            positions.append(node.positions)
            if node is self.upper:
                forces.append(self.upper_force())
            elif node is self.lower:
                forces.append(self.lower_force())
            else:
                forces.append(node.normal_forces)

        forces = np.array(forces)
        positions = np.array(positions)

        if self.control.gs_global_optimizer:
            save = os.path.join(self.nodes[0].dir_pre, "gs.opt")
            new_pos,opt = self.g_opt(
                self.control.gs_optimizer,
                positions,
                forces,
                save)
        else:
            new_pos,opt = self.nong_opt(
                self.control.gs_optimizer,
                self.nodes,
            positions,
                forces,
                ".gs.opt")

        lower_length = get_total_length(new_pos[:self.lower_end+1])
        upper_length = get_total_length(new_pos[self.upper_end:])
        total_length = get_total_length(spline_pos(np.array(new_pos), np.linspace(0,1,100)))

        # Decide whether to reparm locally or globally.
        lower_ideal = self.lower.param*total_length
        upper_ideal = total_length - self.upper.param*total_length

        if (abs(lower_length-lower_ideal) > 0.1*lower_ideal or 
            abs(upper_length - upper_ideal) > 0.1 * upper_ideal):
            if not (isinstance(opt, FDOptimize) and opt.finite_diff):
                new_pos2 = spline_pos(np.array(new_pos), self.params)
            else:
                new_pos2 = new_pos
        else:
            temp_pos = new_pos[:self.lower_end+1]
            old_t = get_t(temp_pos)
            new_t = np.linspace(0.0,1.0,len(temp_pos))
            
            if not (isinstance(opt, FDOptimize) and opt.finite_diff):
                temp_pos = spline_pos(temp_pos, new_t, old_t)
            else:
                pass
            new_pos2.extend(temp_pos)
            temp_pos = new_pos[self.upper_end:]
            old_t = get_t(temp_pos)
            new_t = np.linspace(0.0,1.0,len(temp_pos))
            if not (isinstance(opt, FDOptimize) and opt.finite_diff):
                temp_pos = spline_pos(temp_pos, new_t, old_t)
            else:
                pass
            new_pos2.extend(temp_pos)

        for i,position in enumerate(new_pos2):
            self.nodes[i].positions = position

        high_force = np.reshape(high_force, (-1,3))
        high_force = np.sum(high_force**2,1)**0.5
        low_force = np.reshape(low_force, (-1,3))
        low_force = np.sum(low_force**2,1)**0.5
        forces = np.reshape(forces, (-1,3))
        forces = np.sum(forces**2,1)**0.5

        return np.nanmax(forces),np.nanmax(low_force),np.nanmax(high_force)


    def move_nodes_combined(self):

        from aimsChain.interpolate import spline_pos, get_total_length
        import os
        
        positions = []
        forces = []
        new_t = self.params
        new_pos = []
        low_force = self.lower_force()
        high_force = self.upper_force()

        for i,node in enumerate(self.nodes):
            positions.append(node.positions)
            forces.append(node.normal_forces)

        forces = np.array(forces)
        positions = np.array(positions)

        if self.control.gs_global_optimizer:
            save = os.path.join(self.nodes[0].dir_pre, "pgs.opt")
            new_pos,opt = self.g_opt(
                self.control.gs_optimizer,
                positions,
                forces,
                save)
        else:
            new_pos,opt = self.nong_opt(
                self.control.gs_optimizer,
                self.nodes,
                positions,
                forces,
                ".gs.opt")

        if not (isinstance(opt, FDOptimize) and opt.finite_diff):
            new_pos = spline_pos(new_pos, new_t)

        for i,position in enumerate(new_pos):
            self.nodes[i].positions = position

        high_force = np.reshape(high_force, (-1,3))
        high_force = np.sum(high_force**2,1)**0.5
        low_force = np.reshape(low_force, (-1,3))
        low_force = np.sum(low_force**2,1)**0.5
        forces = np.reshape(forces, (-1,3))
        forces = np.sum(forces**2,1)**0.5

        return np.nanmax(forces),np.nanmax(low_force),np.nanmax(high_force)


    def write_path(self, filename="path.dat"):
        """ Write the current path into a file

        File contains:
         - current runs
         - each node's param, dir, and whether they are fixed or not
        
        This was preferred over pickling for the sake of easy debugging
        May change to pickling in the production level
        """

        data = {}
        data["runs"] = self.runs
        params = []
        dirs = []
        fix = []
        climb = []

        for node in self.nodes:
            params.append(node.param)
            dirs.append(node.dir)
            fix.append(node.fixed)
            climb.append(node.climb)

        data["param"] = params
        data["dirs"] = dirs
        data["fix"] = fix
        data["climb"] = climb
        data["lower"] = self.lower_end

        with open(filename,'w') as save:
            cp.dump(data,save)


    def read_path(self, filename="path.dat"):
        """ Read the file into current path."""

        nodes = []

        with open(filename, 'r') as save:
            print 
            data = cp.load(save)

        self.runs = data["runs"]
        params = data["param"]
        dirs = data["dirs"]
        fix = data["fix"]
        climb = data["climb"]
        self.lower_end = data["lower"]

        for i,param in enumerate(params):
            tmp_node = Node(param, path=self)
            tmp_node.dir = dirs[i]
            tmp_node.fixed = fix[i]
            tmp_node.climb = climb[i]            
            nodes.append(tmp_node)    

        self.nodes = nodes
