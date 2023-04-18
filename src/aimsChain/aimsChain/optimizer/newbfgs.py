import numpy as np
from numpy.linalg import eigh, solve
import cPickle as cp

class BFGS(object):
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


    def initialize(self):
        self.H = None
        self.r0 = None
        self.f0 = None
        

    def load(self):
        """ Load saved hessian, positions, and forces."""
        
        import os.path as path

        if path.isfile(self.restart):
            with open(self.restart, 'r') as hess:
                try:
                    self.H, self.r0, self.f0  = cp.load(hess)
                except:
                    self.initialize()


    def dump(self):
        """ Dump necessary values for future reference."""

        with open(self.restart, 'w') as hess:
            cp.dump((self.H, self.r0, self.f0),
                    hess)


    def step(self, r, f):
        """ Take a single step."""

        r = np.array(r)
        f = np.array(f).flatten()
        self.update(r, f)
        dr = self.H.dot(f)
        dr = dr.reshape(-1,3)
        dr = self.determine_step(dr)
        dr = np.reshape(dr, np.shape(r))
        self.r0 = r.flatten()
        self.f0 = f.flatten()

        return r+dr



    def determine_step(self, dr):
        """ Determine step to take according to maxstep.
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """

        steplengths = (dr**2).sum(1)**0.5
        maxsteplength = np.max(steplengths)
        if maxsteplength >= self.maxstep:
            dr *= self.maxstep / maxsteplength
        
        return dr



    def update(self, r, f):

        r = r.flatten()
        f = f.flatten()
        if self.H is None or np.shape(r) != np.shape(self.r0):
            self.H = np.eye(len(r)) * (1.0/self.alpha)
            self.r0 = r
            self.f0 = f
            return

        r0 = self.r0
        f0 = self.f0
        sk = r.reshape(-1) - r0.reshape(-1)

        if np.abs(sk).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        yk = f0.reshape(-1) - f.reshape(-1)
        try:
            rhok = 1.0/(np.dot(yk,sk))
        except ZeroDivisionError:
            rhok = 1000.0
        if np.isinf(rhok):
            rhok = 1000.0

#        if rhok >= 0:
        I = np.eye(len(r))
        A1 = I - sk[:,np.newaxis] * yk[np.newaxis,:] * rhok
        A2 = I - yk[:,np.newaxis] * sk[np.newaxis,:] * rhok        
        self.H = (np.dot(A1,np.dot(self.H,A2)) 
                  + rhok * sk[:,np.newaxis] * sk[np.newaxis,:])


