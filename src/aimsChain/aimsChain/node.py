""" This module defines the Node.

Node is intended to be the fundamental unit in string method.
Unlike geometry, it maintains parametric values, which remains 
unchanged through the entire calculation.
"""

import numpy as np
from aimsChain.utility import vmag, vunit, vproj
from aimsChain.config import Control

class Node(object):
    """ Class for a single Node."""

    def __init__(self, 
                 param = 0,
                 geometry = None,
                 path = None,
                 dir = None,
                 fixed = False):
        from aimsChain.atom import Atoms
        self.__param = float(param)
        self.__geometry = geometry
        if self.__geometry == None:
            self.__geometry = Atoms()
        self.__path = path
        self.__dir = dir
        self.__fixed = fixed
        self.__climb = False
        self.__dir_pre = "iterations"
        self.__previous_dir = None

    @property
    def param(self):
        return self.__param
    @param.setter
    def param(self, t):
        self.__param = float(t)
    @property
    def geometry(self):
        return self.__geometry
    @geometry.setter
    def geometry(self, atoms):
        self.__geometry = atoms
    @property
    def climb(self):
        return self.__climb
    @climb.setter
    def climb(self, climb):
        self.__climb = bool(climb)
    
    @property
    def control(self):
        if self.path:
            return self.path.control
        else:
            return Control()

    @property
    def path(self):
        return self.__path
    @path.setter
    def path(self, path):
        self.__path = path

    @property
    def dir(self):
        return self.__dir
    @dir.setter
    def dir(self, dir):
        self.__dir = dir
    
    @property
    def dir_pre(self):
        return self.__dir_pre
    @dir_pre.setter
    def dir_pre(self,dir_pre):
        self.__dir_pre = dir_pre
    @property
    def ener(self):
        return float(self.__geometry.ener)
    @ener.setter
    def ener(self, ener):
        self.__geometry.ener = float(ener)
    
    @property
    def positions(self):
        return self.__geometry.positions
    @positions.setter
    def positions(self, positions):
        if not self.__fixed:
            self.__geometry.positions = positions
    
    @property
    def forces(self):
        if self.fixed:
            return np.zeros(np.shape(self.geometry.forces))
        else:
            return self.geometry.forces
    @forces.setter
    def forces(self, forces):
        self.__geometry.forces = forces

    @property
    def climb_forces(self):
        if self.fixed:
            return np.zeros(np.shape(self.forces))
        elif self.climb:
            forces = self.forces
            tangent = self.get_tangent(for_climb=True, test_tangent = False)
            return forces - 2*vproj(forces,tangent)
        else:
            if self.control.method == "neb":
                return self.spring_forces
            else:
                return self.normal_forces

    @property
    def normal_forces(self):
        if self.fixed:
            return np.zeros(np.shape(self.geometry.forces))
        elif self.prev == None or self.next == None:
            return self.forces
        forces = self.forces
        tangent = self.get_tangent()
#        print ('forces', forces)
#        print ('tangent', tangent)
#        print ('proj', vproj(forces,tangent))
        forces -= vproj(forces, tangent)
        return forces

    @property
    def spring_forces(self):

       	k = self.control.spring_k

        if self.fixed:
            return np.zeros(np.shape(self.positions))
        elif self.prev==None or self.next==None:
            return self.forces

        forces = self.normal_forces
        tangent = self.get_tangent()
        tan1 = self.next.positions - self.positions
        tan2 = self.positions - self.prev.positions
        mag = (vmag(tan1) - vmag(tan2))
        forces += k * mag * tangent

        return forces
    
   
    @property
    def ener(self):
        return self.__geometry.ener
    @ener.setter
    def ener(self, energy):
        self.__geometry.ener = energy

    @property
    def fixed(self):
        return self.__fixed
    @fixed.setter
    def fixed(self, value):
        self.__fixed = bool(value)

    @property
    def prev(self):
        """ Get the previous element in path."""
        if not self.path:
            return None
        i = self.path.nodes.index(self)
        if i == 0:
            return None
        else:
            return self.path.nodes[i-1]
    @property
    def next(self):
        """ Get the next element in path."""
        if not self.path:
            return None
        i = self.path.nodes.index(self)
        if i+1 == len(self.path.nodes):
            return None
        else:
            return self.path.nodes[i+1]

    @property
    def previous_dir(self):
        """ Store the previous dir for this node."""
        return self.__previous_dir
    @previous_dir.setter
    def previous_dir(self, value):
        self.__previous_dir = value


    def get_tangent(self, for_climb=False, unit=True, test_tangent = False):

        from aimsChain.utility import vunit
        from aimsChain.interpolate import get_t, spline_pos

        prev = self.prev
        next = self.next
        tangent = None

        if not prev:
            tangent = (next.positions - self.positions)
        elif not next:
            tangent = (self.positions - prev.positions)
        else:
            tan1 = next.positions - self.positions
            tan2 = self.positions - prev.positions
            if next.ener >= self.ener and self.ener >= prev.ener and not for_climb:
                tangent = tan1
            elif next.ener < self.ener and self.ener < prev.ener and not for_climb:
                tangent = tan2
            elif not test_tangent:
                max = np.nanmax((abs(next.ener-self.ener), abs(self.ener-prev.ener)))
                min = np.nanmin((abs(next.ener-self.ener), abs(self.ener-prev.ener)))
                if next.ener >= prev.ener:
                    tangent = (max * tan1 + min * tan2)
                else:
                    tangent = (min * tan1 + max * tan2)
            else:
                positions = [prev.positions, self.positions, next.positions]
                loc = 1
                if prev.prev:
                    positions.insert(0,prev.prev.positions)
                    loc += 1
                if next.next:
                    positions.append(next.next.positions)
                param = get_t(positions)
                derv = spline_pos(positions, param, param, 3, 1)[loc]
                tangent = derv
        if unit:
            return vunit(tangent)
        return tangent


    def write_node(self, write_fixed = False, control_file = "control.in"):
        """ Create a directory for each node including geometry and control."""

        import glob, os, shutil
        from aimsChain.aimsio import write_aims

        dir = os.path.join(self.dir_pre, self.dir)

        if (not self.fixed) or write_fixed:
            if not os.path.isdir(dir):
                os.makedirs(dir)
            write_aims(os.path.join(dir, "geometry.in"), self.geometry)
            shutil.copy(control_file, os.path.join(dir, "control.in"))
            if (self.path and self.control.aims_restart
                and self.previous_dir):
                if self.previous_dir != self.dir:
                    restart_file = self.control.aims_restart + "*"
                    files = glob.glob(os.path.join(self.dir_pre, self.previous_dir, restart_file))
                    for file in files:
                        if os.path.isfile(file):
                            shutil.copy(file, dir)
            return dir
        else:
            return None

    
    def update_dir(self, write_fixed = False):
        """ Update the dirs for all nodes.

        If write_fixed = True, even fixed nodes will get updates.
        """

        if (not self.fixed) or write_fixed:
            self.previous_dir = self.dir
            self.dir = "iteration%04d/" % self.path.runs + "aims-chain-node-%.5f" % (self.param)
                

