#!/usr/bin/env python
"""
Script that runs aimsChain.

Module-wide global variables are used throughout the script to avoid passing
args around. Care should be taken when modifying these variables.

"""

import subprocess
import os
import sys
import distutils.dir_util as dir_util
import shutil
import numpy as np
import cPickle as cp

from aimsChain.string_path import StringPath
from aimsChain.neb_path import NebPath
from aimsChain.gs_path import GrowingStringPath
from aimsChain.node import Node
from aimsChain.aimsio import read_aims
from aimsChain.config import Control
from aimsChain.interpolate import get_t
from aimsChain.aimsio import write_mapped_aims, write_xyz, write_aims


def run_aims(paths):
    """ Does a single run of FHI-aims inside the given directory. """

    global control

    while len(paths) > 0:
        path = paths[0]
        if control.restart:
            save_restart(paths)
        # Remove the ending slash if it exist.
        if path[-1] == "/":
            path = path[:-1]
        # Generate the name for output.
        filename = path[path.rfind('/')+1:]+'.out'      
        print 'Starting FHI-aims in directory ' + path
        command = 'cd ' + path + ';' + control.run_aims + ' > ' + filename
        # Directly call shell and aborts the python script if aims crashes.
        subprocess.check_call(command, shell=True)
        paths.remove(path)
        if control.restart:
            save_restart(paths)        

    paths = []


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

    # Does periodic transformation between initial and final node
    # minimize distance between initial and final image atom-wise
    if control.periodic_interp and finnode.geometry.lattice != None:
        lattice = finnode.geometry.lattice
        initial_pos = ininode.positions
        curr_pos = finnode.positions
        fin_pos = []
        for i,atom_pos in enumerate(curr_pos):
            #large separation to start with
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

    if control.use_gs:
        path.nodes = [ininode, finnode]
        path.add_lower()
        path.add_upper()
        return

    # Parse the external geometry
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
            #if resample is turned on and n_image does not equal to 
            #number of images from external source,
            #then resample the path
            if control.resample and control.nimage != (len(nodes)-2):
                path.interpolate(control.nimage)
    # Revert to linear interpolation if anything failed. This could be:
    # -file doesn't exist
    # -wrong format
    # -interpolation error
    except:
        print '!Error interpolating the external geometries\n'
        print '!Using standard interpolation method for initial geometries\n'
        nodes = [ininode, finnode]

    # If there were no external geometry, linear interpolate the image.
    if len(nodes) <= 2:
        path.nodes = [ininode, finnode]
        #interpolate the images
        path.interpolate(control.nimage)


def save_restart(path):
    """ Creates the restart file."""

    global restart_stage
    global force

    print 'Writing restart file to ./iterations/restart_file'
    with open("iterations/restart_file", 'w') as restart:
        cp.dump((path, restart_stage, force), restart)


def read_restart():
    """ Loads the restart file. """

    global path_to_run
    global force
    global restart_stage
    global path

    file_exist = os.path.isfile("iterations/restart_file")

    if file_exist:
        print 'Loading restart file in ./iterations/restart_file'
        with open("iterations/restart_file", 'r') as restart:
            path_to_run, restart_stage, force = cp.load(restart)
            path.read_path("iterations/path.dat")
            path.load_nodes()

    return file_exist


def write_current(final = False):
    """ Writes the current path."""

    global path

    if not final:
        dir_name = "iteration%04d" % path.runs
        dir_name = os.path.join("paths", dir_name)
    else:
        dir_name = "optimized"        

    try:
        os.mkdir(dir_name)
        print 'Making ' + os.path.join('.', str(dir_name))
    except OSError:
        print 'OSError: Failed to make ' + os.path.join('.', str(dir_name)) 

    ener = []
    climb = []
    fixed = []
    print 'Writing path.xyz to ' + os.path.join('.', str(dir_name))  
    write_xyz(os.path.join(dir_name, "path.xyz"), path, control.xyz_lattice)

    for i,node in enumerate(path.nodes):
        ener.append(node.ener)
        climb.append(node.climb)
        fixed.append(node.fixed)
        i += 1
        file_name = os.path.join(dir_name, "image%03d.in" % i)
        if not final:
            if control.map_unit_cell:
                write_mapped_aims(file_name, node.geometry)
            else:
                write_aims(file_name, node.geometry)

    ener = np.array(ener) - ener[0]

    with open(os.path.join(dir_name, "ener.lst"), 'w') as enerfile:
        enerfile.write("#image num.\tenergy(eV)\t status\n")
        for i,energy in enumerate(ener):
            enerfile.write("image%03d\t%.10f" % (i+1, energy))
            if climb[i]:
                enerfile.write("\t CLIMB")
            elif fixed[i]:
                enerfile.write("\t FIXED")
            else:
                enerfile.write("\t NORMAL")
            enerfile.write("\n")   

    with open(os.path.join(dir_name, "path.lst") ,'w') as pathfile:
        pathfile.write("#image num.\tpath\t\t\t\t\t\t\t status\n")
        for i,item in enumerate(path.get_paths()):
            pathfile.write("image%03d\t%s" % (i+1,item))
            if climb[i]:
                pathfile.write("\t CLIMB")
            elif fixed[i]:
                pathfile.write("\t FIXED")
            else:
                pathfile.write("\t NORMAL")
            pathfile.write("\n")

    if not final:
        print '------------------------------------------------------------' 
        print 'End of iteration%04d' %path.runs
        print '------------------------------------------------------------\n'


def print_control():
    """ Echoes the chain.in file to the output."""

    print '------------------------------------------------------------'
    print 'Parsing chain.in' 
    print 'The contents of chain.in will be repeated verbatim below.'
    print '------------------------------------------------------------' 

    with open("chain.in", 'r') as control: 
        chain_input = control.read()

    print chain_input
    print '------------------------------------------------------------' 
    print 'End of chain.in input file.'
    print '------------------------------------------------------------\n\n\n'



print '------------------------------------------------------------'
print '       Invoking aimsChain ...'
print '       For any questions about aimsChain please consult the '
print '       FHI-aims usual manual.'
print '------------------------------------------------------------\n\n\n'
print_control()
    
force = 10.0
control = Control()
print '------------------------------------------------------------' 
print 'Begining aimsChain calculation...'
print '------------------------------------------------------------' 

if control.method == "neb":
    path = NebPath(control=control)
else:
    path = StringPath(control=control)

restart_stage = "mep"

is_restart = control.restart and read_restart() 

if control.use_gs:
    if not is_restart or restart_stage == "growing":
        path = GrowingStringPath(control=control)
        read_restart()
        restart_stage = "growing"
        growing = True
        gsforcelog = open("growing_forces.log", 'a')

# Check if the system is a restart.
if not is_restart:
    for directory in ["paths", "iterations", "optimized"]:
        if os.path.isdir(directory):
            shutil.rmtree(directory)
    os.mkdir("paths")
    initial_interpolation()
    # Write directory for images.
    path_to_run = path.write_all_node()
    path.write_path("iterations/path.dat")
    if control.use_gs:
        gsforcelog.write("Iteration\tResidual force\t\tLower end force\t\tUpper end force \n")
        gsforcelog.flush()

if restart_stage == "growing":
    while growing:
        run_aims(path_to_run)
        path.load_nodes()
        force,low_force,high_force = path.move_nodes()
        write_current()
        curr_runs = path.runs
        if low_force < control.gs_thres:
            growing = path.add_lower() and growing
        if high_force < control.gs_thres:
            growing = path.add_upper() and growing
        path.add_runs()
        path_to_run = path.write_node()
        gsforcelog.write('iteration%04d\t%16.8f\t%16.8f\t%16.8f \n' % 
                       (curr_runs,force,low_force, high_force))
        gsforcelog.flush()
        path.write_path("iterations/path.dat")
    run_aims(path_to_run)
    write_current()
    path.write_path("iterations/path.dat")
    force = 10.0
    print "Path is grown, see growing_forces.log."
    gsforcelog.write("Path is grown.\n")   
    gsforcelog.close()
    path.write_path("iterations/path.dat")
    try:
        os.mkdir('grownstring')
    except OSError:
        print 'OSError: failed to make grownstring directory.'
    for i,dir in enumerate(path.get_paths()):
        i+=1
        dir_util.copy_tree(dir, os.path.join('grownstring', "image%03d" % i))
    write_xyz("grownstring/path.xyz", path, control.xyz_lattice)
    restart_stage = "grown"

if restart_stage == "grown":
    if control.method == "neb":
        path = NebPath(control=control)
    else:
        path = StringPath(control=control)
    path.read_path("iterations/path.dat")
    path.load_nodes()
    # Resample the path.
    path.interpolate(control.nimage)
    path.add_runs()
    path_to_run = path.write_all_node()    
    restart_stage = "mep"

forcelog = open("forces.log", 'a')

if not is_restart:
    forcelog.write("#Residual Forces in the system:\n")
    forcelog.flush()

if restart_stage == "mep":
    while force > control.thres:
        run_aims(path_to_run)
        path.load_nodes()
        force = path.move_nodes()
        write_current()
        curr_runs = path.runs
        if force > control.thres:
            path.add_runs()
            path_to_run = path.write_node()
        forcelog.write('iteration%04d\t%16.10f\n' % (curr_runs,force))
        forcelog.flush()
        path.write_path("iterations/path.dat")
    force = 10.0
    forcelog.write("System has converged.\n")
    print "\nSystem has converged, see forces.log for more information.\n"

forcelog.close()

if control.use_climb:
    forcelog = open("climbing_forces.log", 'a')
    if restart_stage != "CI":
        path.find_climb()
        path.add_runs()
        if control.climb_control != "control.in":
            path_to_run = path.write_all_node(control.climb_control)
        else:
            path_to_run = path.write_node()
        forcelog.write("#Residual Forces in the Climbing image:\n")
        forcelog.flush()
        restart_stage = "CI"
    while force > control.climb_thres:
        run_aims(path_to_run)
        path.load_nodes()
        force = path.move_climb()
        write_current()
        curr_runs = path.runs
        if force > control.climb_thres:
            path.add_runs()
            path_to_run = path.write_node(control.climb_control)
        forcelog.write('iteration%04d\t%16.16f \n' % (curr_runs, force))
        forcelog.flush()
        path.write_path("iterations/path.dat")
    forcelog.write('Climbing image has converged.\n')
    print "\nClimbing image has converged, see climbing_forces.log for more information.\n"
    forcelog.close()

for i,dir in enumerate(path.get_paths()):
    i+=1
    dir_util.copy_tree(dir, os.path.join('optimized', "image%03d" % i))

write_current(True)

print '\n------------------------------------------------------------'
print '             aimsChain is finished running.'
print '             Have a nice day.'
print '------------------------------------------------------------'
