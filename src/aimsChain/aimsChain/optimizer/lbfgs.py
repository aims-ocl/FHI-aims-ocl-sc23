import sys
import numpy as np
import cPickle as cp


class LBFGS(object):
    """Limited memory BFGS optimizer.
    
    A limited memory version of the bfgs algorithm. Unlike the bfgs algorithm
    used in bfgs.py, the inverse of Hessian matrix is updated.  The inverse
    Hessian is represented only as a diagonal matrix to save memory

    Modified from ASE's LBFGS optimizer

    """
    def __init__(self, restart="hess",
                 maxstep=None, memory=25, 
                 damping = 1.0, alpha = 70.0):
        """
        Parameters:

        restart: string
            Pickle file used to store vectors for updating the inverse of Hessian
            matrix. If set, file with such a name will be searched and information
            stored will be used, if the file exists.

        maxstep: float
            How far is a single atom allowed to move. This is useful for DFT
            calculations where wavefunctions can be reused if steps are small.
            Default is 0.04 Angstrom.

        memory: int
            Number of steps to be stored. Default value is 100. Three numpy
            arrays of this length containing floats are stored.

        damping: float
            The calculated step is multiplied with this number before added to
            the positions. 

        alpha: float
            Initial guess for the Hessian (curvature of energy surface). A
            conservative value of 70.0 is the default, but number of needed
            steps to converge might be less if a lower value is used. However,
            a lower value also means risk of instability.
            
        """

        if maxstep is not None:
            if maxstep > 1.0:
                raise ValueError('You are using a much too large value for ' +
                                 'the maximum step size: %.1f Angstrom' % maxstep)
            self.maxstep = maxstep
        else:
            self.maxstep = 0.04

        self.memory = memory
        self.H0 = 1. / alpha  # Initial approximation of inverse Hessian
                            # 1./70. is to emulate the behaviour of BFGS
                            # Note that this is never changed!
        self.damping = damping
        self.p = None
        self.restart = restart


    def initialize(self):
        """ Initalize everything so no checks have to be done in step."""

        self.iteration = 0
        self.s = []
        self.y = []
        self.rho = [] # Store also rho, to avoid calculationg the dot product
                      # again and again
        self.r0 = None
        self.f0 = None


    def load(self):
        """ Load saved arrays to reconstruct the Hessian."""

        import os.path as path

        if path.isfile(self.restart):
            with open(self.restart, 'r') as hess:
                self.iteration, self.s, self.y, self.rho, \
                    self.r0, self.f0 = cp.load(hess)


    def dump(self):
        """ Dump necessary values for future reference."""

        with open(self.restart, 'w') as hess:
            cp.dump((self.iteration, self.s, self.y, 
                       self.rho, self.r0, self.f0), hess)


    def step(self, r, f):
        """Take a single step
        
        Use the given forces, update the history and calculate the next step --
        then take it.
        """

        f = np.array(f)
        r = np.array(r)    
        self.update(r, f, self.r0, self.f0)        
        s = self.s
        y = self.y
        rho = self.rho
        H0 = self.H0
        loopmax = len(self.s)
        a = np.empty((loopmax,), dtype=np.float64)

        ### The algorithm itself:
        q = - f.reshape(-1) 
        for i in range(loopmax - 1, -1, -1):
            a[i] = rho[i] * np.dot(s[i], q)
            q -= a[i] * y[i]
        z = H0 * q
        
        for i in range(loopmax):
            b = rho[i] * np.dot(y[i], z)
            z += s[i] * (a[i] - b)

        self.p = - z.reshape((-1, 3))
        ###
        dr = self.determine_step(self.p * self.damping) 
        dr = np.reshape(dr, np.shape(r))
        self.r0 = r.copy()
        self.f0 = f.copy()

        return r+dr


    def approxdr(self, f):
        """ Approx. the change in r according to f hessian is not updated."""

        f = np.array(f)    
        s = self.s
        y = self.y
        rho = self.rho
        H0 = self.H0
        loopmax = len(self.s)
        a = np.empty((loopmax,), dtype=np.float64)
        ### The algorithm itself:
        q = - f.reshape(-1) 
        for i in range(loopmax - 1, -1, -1):
            a[i] = rho[i] * np.dot(s[i], q)
            q -= a[i] * y[i]
        z = H0 * q
        
        for i in range(loopmax):
            b = rho[i] * np.dot(y[i], z)
            z += s[i] * (a[i] - b)

        self.p = - z.reshape((-1, 3))
        ###
#        dr = self.p * self.damping
        dr = self.determine_step(self.p) * self.damping
        dr = np.reshape(dr, np.shape(f))
        return dr


    def determine_step(self, dr):
        """Determine step to take according to maxstep.
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """

        steplengths = (dr**2).sum(1)**0.5
        longest_step = np.max(steplengths)
        if longest_step >= self.maxstep:
            dr *= self.maxstep / longest_step        
        return dr


    def update(self, r, f, r0=None, f0=None):
        """Update everything that is kept in memory.

        This function is mostly here to allow for replay_trajectory.
        """

        r = np.array(r)
        f = np.array(f)
        if r0 == None:
            if self.r0 == None:
                self.r0 = r.copy()
                self.f0 = f.copy()
            else:
                r0 = self.r0
                f0 = self.f0

        #reset direction        
        if np.shape(f) != np.shape(f0):
            self.initialize()
        if self.iteration > 0:
            a1 = abs(np.dot(f.reshape(-1),f0.reshape(-1)))
            a2 = np.dot(f0.reshape(-1),f0.reshape(-1))
            if (a1 >= 0.5*a2) or (a2 == 0.0):
                self.initialize()        
        if self.iteration > 0:
            s0 = r.reshape(-1) - r0.reshape(-1)
            if np.abs(s0).max() < 1e-7:
                # Same configuration again
                return
            # We use the gradient which is minus the force!
            y0 = f0.reshape(-1) - f.reshape(-1)
            rho0 = 1.0 / np.dot(y0, s0)            
            self.y.append(y0)
            self.s.append(s0)
            self.rho.append(rho0)
        if len(self.y) > self.memory:
            self.s.pop(0)
            self.y.pop(0)
            self.rho.pop(0)
        self.iteration += 1
