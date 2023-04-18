import numpy as np
from numpy.linalg import eigh, solve

import cPickle as cp

class choleskyBFGS(object):
    def __init__(self, restart="hess", maxstep=0.04, 
                 alpha = 70):
        """BFGS optimizer.while force > 0.01:

        Parameters:

        restart: string
            Pickle file used to store hessian matrix. 
            load() will attempt to load existing Hessian from it.
            dump() will dump the current Hessian to it. 
        Maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.04 A).
        """
        
        if maxstep is not None:
            if maxstep > 1.0:
                raise ValueError('You are using a much too large value for ' +
                                 'the maximum step size: %.1f A' % maxstep)
            self.maxstep = maxstep
        self.restart = restart
        self.alpha = alpha
        self.save_hess = True
        self.writelog = True


    def log(self, str_in=""):
        if not self.writelog:
            return
        log = open('bfgs.log','a')
        log.write("CholeskyIteration[%04d]: " % self.iteration)
        log.write(str_in + '\n')
        log.close()


    def initialize(self):
        self.H = None
        self.r0 = None
        self.f0 = None
        self.inserted = False
        self.inserted_counter = 0
        self.iteration = 0

        
    def insert_node(self, nth, n_atoms):
        self.log("Inserting new node")

        if self.H == None:
            return
        loc = int(nth*n_atoms*3)
        H = self.H

        for _ in range(int(n_atoms*3)):
            H = np.insert(H,loc,0,axis=0)

        for _ in range(int(n_atoms*3)):
            H = np.insert(H,loc,0,axis=1)

        new_n = len(H)
        I = np.eye(new_n) * (self.alpha)

        for i in range(new_n):
            for j in range(new_n):
                if H[i][j] == 0:
                    H[i][j] = I[i][j]
        self.inserted = True
        self.inserted_counter = 0
        self.H = H
        self.dump()
            

    def load(self):
        """ Load saved hessian, positions, and forces."""

        import os.path as path

        if path.isfile(self.restart):
            with open(self.restart, 'r') as hess:
                try:
                    (self.H, self.r0, self.f0, self.inserted, 
                     self.inserted_counter, self.iteration) = cp.load(hess)
                except:
                    self.initialize()


    def dump(self):
        """ Dump necessary values for future reference."""

        with open(self.restart, 'w') as hess:
            cp.dump((self.H, self.r0, self.f0, self.inserted, 
                     self.inserted_counter, self.iteration),
                    hess)


    def step(self, r, f):
        """ Take a single step."""

        import time

        self.log("steping the nodes")
        start_t = time.time()
        r = np.array(r)
        f = np.array(f).flatten()        

        if (self.f0 == None) or len(self.f0) == len(f):
            self.update(r, f)
        chole = np.linalg.cholesky(self.H)
        tempy = np.linalg.solve(chole, f)
        dr = np.linalg.solve(chole.T.conj(),tempy)            
        dr = dr.reshape(-1,3)
        dr = self.determine_step(dr)
        dr = np.reshape(dr, np.shape(r))
        self.r0 = r.flatten()
        self.f0 = f.flatten()
        finish_t = time.time() - start_t
        self.log("Used " + str(finish_t) + " seconds in steping")

        if self.save_hess:
            filename = self.restart + ".hess"
            f_handle = open(filename, 'a')
            np.savetxt(f_handle, self.H)
            f_handle.write("\n\n")
            f_handle.close()        
        return r+dr


    def determine_step(self, dr):
        """Determine step to take according to maxstep
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        
        steplengths = (dr**2).sum(1)**0.5
        maxsteplength = np.max(steplengths)
        if maxsteplength >= self.maxstep:
            dr *= self.maxstep / maxsteplength        
        return dr


    def update(self, r, f):
        self.iteration += 1
        self.inserted_counter += 1
        r = r.flatten()
        f = f.flatten()

        if self.H is None:
            self.H = np.eye(len(r)) * self.alpha
            return

        r0 = self.r0
        f0 = self.f0
        a1 = abs(np.dot(f,f0))
        a2 = np.dot(f0,f0)

        if ((a1 >= a2) or (a2 == 0.0)) and (self.inserted_counter > 2):
            self.log("resetting the hessian")
            self.H = np.eye(len(r)) * (self.alpha) * 0.5
            return

        sk = r.reshape(-1) - r0.reshape(-1)
        yk = f0.reshape(-1) - f.reshape(-1) #force is negative
        rhok=np.dot(yk,sk)
        if np.abs(sk).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        #local curvature is convex
        if rhok <= 0:
            return

        #BFGS update formula
        try:
            rhok = 1.0/rhok
        except ZeroDivisionError:
            rhok = 1e5

        if np.isinf(rhok) or float(rhok) >= 1e5:
            self.log("damping rhok")
            rhok = 1e5

        I = np.eye(len(r))
        tk = np.dot(self.H,sk)        
        self.H = (self.H + np.outer(yk,yk) * rhok 
                  - np.outer(tk,tk) / np.dot(sk,tk))
        
