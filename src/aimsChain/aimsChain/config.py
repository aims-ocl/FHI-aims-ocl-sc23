"""
This module provides the parser for config file.

Config should be named "chain.in"

"""

import os 
import string

class Control(object):
    """ Class for storing and parsing config.
    
    A very barebones class for storing and parsing config that lacks robustness.
    Should be sufficient for our use case.
    """

    def __init__(self):
        #initial geometry
        self.ini = "ini.in"
        #final geometry
        self.fin = "fin.in"
        #number of images using
        self.nimage = 5
        #number of gs images using
        self.gs_nimage = None
        #periodic interpolation
        self.periodic_interp = False
        #threshold for convergence
        self.thres = 0.2
        #threshold for convergence
        self.gs_thres = None
        #threshold for climbing image convergence
        self.climb_thres = None
        #climbing image or not
        self.use_climb = False
        #use growing string or not
        self.use_gs = False
        #climbing by interpolation
        self.climb_interp = True
        #climbing mode
        #1 for single node climb
        #2 for 3 node climb
        #3 for all node climb
        self.climb_mode = 2
        #command to run aims
        self.run_aims = "mpiexec -ppn 8 -n $NSLOTS ~/bin/aims.081912.scalapack.mpi.x"
        #using external starting geometry
        self.ext_geo = None
        #method, string or neb
        self.method = "string"
        #global optimizer
        self.global_opt = True
        #global optimizer for climbing image
        self.climb_global_opt = True
        #spring constant
        self.spring_k = 20.0
        #restart file
        self.aims_restart = None
        #resample external geometry
        self.resample = False
        #restart or not
        self.restart = True
        #optimizer for evolving path
        self.optimizer = "trm"
        #optimizer for evolving path
        self.gs_optimizer = "trm"
        #global optimizer for gs
        self.gs_global_optimizer = False
        #optimizer for climbing image
        self.climb_optimizer = "fire"
        #control file for climbing image
        self.climb_control = "control.in"
        #lbfgs parameters
        self.lbfgs_alpha = 120.0
        self.lbfgs_memory = 30
        self.lbfgs_maxstep = 0.04
        #bfgs parameters
        self.bfgs_alpha = 120.0
        self.bfgs_maxstep = 0.04
        #fire parameters
        self.fire_dt = 0.02
        self.fire_maxstep = 0.04
        self.fire_dtmax = 1.0
        self.fire_nmin = 5
        self.fire_finc = 1.1
        self.fire_fdec = 0.5
        self.fire_astart = 0.1
        self.fire_fa = 0.99
        self.fire_a = 0.1
        #map back to the central cell?
        self.map_unit_cell = False
        #lattice view for 
        self.xyz_lattice = [2,2,1]    
        self.read()

    
    def read(self, filename = "chain.in"):
        if os.path.isfile(filename):
            with open(filename, 'r') as control:
                lines = control.readlines()
            for line in lines:
                inp = line.split()
                if inp == []:
                    continue
                elif inp[0][0] == '#':
                    continue
                elif inp[0] == "initial_file":
                    self.ini = inp[1]
                elif inp[0] == "final_file":
                    self.fin = inp[1]
                elif inp[0] == "n_images" or inp[0] == "n_image":
                    self.nimage = int(inp[1])
                elif inp[0] == "gs_n_images" or inp[0] == "gs_n_image":
                    self.gs_nimage = int(inp[1])
                elif inp[0] == "force_thres":
                    self.thres = float(inp[1])
                elif inp[0] == "climb_thres":
                    self.climb_thres = float(inp[1])
                elif inp[0] == "gs_thres":
                    self.gs_thres = float(inp[1])
                elif inp[0] == "use_climb":
                    self.use_climb = parse_bool(inp[1])
                elif inp[0] == "use_gs_method":
                    self.use_gs = parse_bool(inp[1])
                elif inp[0] == "climb_mode":
                    self.climb_mode = int(inp[1])
                elif inp[0] == "run_aims":
                    self.run_aims = string.join(inp[1:])
                elif inp[0] == "external_geometry":
                    self.ext_geo = inp[1]
                elif inp[0] == "method":
                    if inp[1].lower() == "neb":
                        self.method = "neb"
                    else:
                        self.method = 'string'
                elif inp[0] == "global_optimizer":
                    self.global_opt = parse_bool(inp[1])
                elif inp[0] == "interpolated_climb":
                    self.climb_interp = parse_bool(inp[1])
                elif inp[0] == "neb_spring_constant":
                    self.spring_k = float(inp[1])
                elif inp[0] == "periodic_interpolation":
                    self.periodic_interp = parse_bool(inp[1])
                elif inp[0] == "aims_restart":
                    self.aims_restart = str(inp[1])
                elif inp[0] == "resample":
                    self.resample = parse_bool(inp[1])
                elif inp[0] == "restart":
                    self.restart = parse_bool(inp[1])
                elif inp[0] == "climb_global_optimizer":
                    self.climb_global_opt = parse_bool(inp[1])
                elif inp[0] == "optimizer":
                    self.optimizer = str(inp[1])
                elif inp[0] == "gs_optimizer":
                    self.gs_optimizer = str(inp[1])
                elif inp[0] == "gs_global_optimizer":
                    self.gs_global_optimizer = parse_bool(inp[1])
                elif inp[0] == "climb_optimizer":
                    self.climb_optimizer = str(inp[1])
                elif inp[0] == "climb_control":
                    self.climb_control = str(inp[1])
                elif inp[0] == "lbfgs_alpha":
                    self.lbfgs_alpha = float(inp[1])
                elif inp[0] == "lbfgs_memory":
                    self.lbfgs_memory = int(inp[1])
                elif inp[0] == "lbfgs_maxstep":
                    self.lbfgs_maxstep = float(inp[1])
                elif inp[0] == "bfgs_alpha":
                    self.bfgs_alpha = float(inp[1])
                elif inp[0] == "bfgs_maxstep":
                    self.bfgs_maxstep = float(inp[1])
                elif inp[0] == "fire_dt":
                    self.fire_dt = float(inp[1])
                elif inp[0] == "fire_maxstep":
                    self.fire_maxstep = float(inp[1])
                elif inp[0] == "fire_dtmax":
                    self.fire_dtmax = float(inp[1])
                elif inp[0] == "fire_nmin":
                    self.fire_nmin = int(inp[1])
                elif inp[0] == "fire_finc":
                    self.fire_finc = float(inp[1])
                elif inp[0] == "fire_fdec":
                    self.fire_fdec = float(inp[1])
                elif inp[0] == "fire_astart":
                    self.fire_astart = float(inp[1])
                elif inp[0] == "fire_fa":
                    self.fire_fa = float(inp[1])
                elif inp[0] == "fire_a":
                    self.fire_a = float(inp[1])               
                elif inp[0] == "map_unit_cell":
                    self.map_unit_cell = parse_bool(inp[1])
                elif inp[0] == "xyz_lattice":
                    self.xyz_lattice = [int(inp[1]), int(inp[2]), int(inp[3])]
                    if len(self.xyz_lattice) != 3:
                        self.xyz_lattice = [2,2,1]        
        # Assign climbing thres if it is not set
        if self.climb_thres == None:
            self.climb_thres = self.thres
        if self.gs_thres == None:
            self.gs_thres = self.thres*1.5
        if self.gs_nimage == None:
            self.gs_nimage = self.nimage
        if self.ext_geo != None:
            self.use_gs = False


def parse_bool(input_string):
    """ Parse a string into boolean.
    
    Returns true for some common variations, returns false otherwise.
    """

    if input_string.lower() in ['true','.true.', 'y', 'yes', '1', 't','on']:
        return True
    else:
        return False
