#!/usr/bin/env python

import subprocess
import os
import distutils.dir_util as dir_util
import shutil
import numpy as np
import cPickle as cp

from aimsChain.path import Path
from aimsChain.node import Node
from aimsChain.aimsio import read_aims
from aimsChain.config import Control
from aimsChain.interpolate import get_t
from aimsChain.aimsio import write_mapped_aims, write_xyz, write_aims


def initial_interpolation():
    """ Performs the initial interpolation for aimsChain
    
    Reads the initial and final geometry files and interpolates based off the
    given scheme. Ugly but works.
    """

    global control
    global path

    # Set the initial and final image.
    ininode = Node(param = 0.0, 
                   geometry = read_aims(control.ini))
    finnode = Node(param = 1.0,
                   geometry = read_aims(control.fin))

    #does periodic transformation between initial and final node
    #minimize distance between initial and final image atom-wise
    if control.periodic_interp and finnode.geometry.lattice is not None:
        lattice = finnode.geometry.lattice
        initial_pos = ininode.positions
        curr_pos = finnode.positions
        fin_pos = []
        for i,atom_pos in enumerate(curr_pos):
            dis = 1000000
            pos = None
            for a in [-1,0,1]:
                for b in [-1,0,1]:
                    for c in [-1,0,1]:
                        pos_tmp = (atom_pos +
                                   a*lattice[0] + b*lattice[1] + c*lattice[2])
                        dis_tmp = np.sum(np.array(pos_tmp - initial_pos[i])**2)**0.5
                        if dis_tmp <= dis:
                            pos = pos_tmp
                            dis = dis_tmp
            fin_pos.append(pos)
        finnode.positions = fin_pos                    

    # Fix the initial and final node.
    ininode.fixed = True
    finnode.fixed = True

    # Parse the externl geometry.
    try:
        nodes = [ininode]
        if control.ext_geo and os.path.isfile(str(control.ext_geo)):
            with open(control.ext_geo) as geo:
                lines = geo.readlines()
            for line in lines:
                inp = line.split()
                if inp == []:
                    continue
                elif inp[0][0] == '#':
                    continue
                else:
                    if os.path.isfile(inp[0]):
                        nodes.append(Node(param = 0.5, 
                                          geometry = read_aims(inp[0]), 
                                          fixed = False))
            nodes.append(finnode)
            geos = []
            for node in nodes:
                geos.append(node.positions)
            t = get_t(geos)
            for i in range(len(t)):
                nodes[i].param = t[i]
            path.nodes = nodes
            # If resample is turned on and n_image does not equal to 
            # number of images from external source,
            # then resample the path.
            if control.resample and control.nimage != (len(nodes)-2):
                path.interpolate(control.nimage)

    except:
        print '!Error interprating the external geometries\n'
        print '!Using standard interpolation method for initial geometries\n'
        raise

    # If there were no external geometry, linear interpolate the image.
    if len(nodes) <= 2:
        path.nodes = [ininode, finnode]
        #interpolate the images
        path.interpolate(control.nimage)

def save_restart(path):
    """ Creates the restat file."""

    global restart_stage
    global force

    with open("iterations/restart_file", 'w') as restart:
        cp.dump((path, restart_stage, force), restart)

def read_restart():
    """ Loads the restart file."""

    global path_to_run
    global force
    global restart_stage
    global path

    file_exist = os.path.isfile("iterations/restart_file")

    if file_exist:
        with open("iterations/restart_file", 'r') as restart:
            path_to_run, restart_stage, force = cp.load(restart)
            path.read_path("iterations/path.dat")
            path.load_nodes()

    return file_exist

control = Control()
path = Path(control=control)

if os.path.isdir("interpolation"):
    shutil.rmtree("interpolation")

os.mkdir("interpolation")
initial_interpolation()
write_xyz(os.path.join("interpolation", "path.xyz"), path, control.xyz_lattice)

for i,node in enumerate(path.nodes):
    i += 1
    file_name = os.path.join("interpolation", "image%03d.in" % i)
    write_aims(file_name, node.geometry)

