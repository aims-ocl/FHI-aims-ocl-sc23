import numpy as np
import cPickle as cp
from aimsChain.utility import vmag, vunit


class FIRE(object):
    def __init__(self, restart="fire.dat",
                 dt=0.02, maxstep=0.04, dtmax=1.0, Nmin=5, finc=1.1, fdec=0.5,
                 astart=0.1, fa=0.99, a=0.1):

        self.dt = dt
        self.Nsteps = 0
        self.maxstep = maxstep
        self.dtmax = dtmax
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa
        self.a = a
        self.restart = restart


    def initialize(self):
        self.v = None


    def load(self):
        """ Load saved velocity and dt from file."""

        import os.path as path

        if path.isfile(self.restart):
            with open(self.restart,'r') as save:
                try:
                    self.v, self.dt, self.a, self.Nsteps = cp.load(save)
                except:
                    self.initialize()


    def dump(self):
        """ Dump necessary values for future reference."""

        with open(self.restart, 'w') as save:
            cp.dump((self.v, self.dt, self.a, self.Nsteps), save)
       

    def step(self,r,f):
        r = np.array(r)
        f = np.reshape(f,(-1,3))

        if self.v is None:
            self.v = np.zeros((len(f), 3))
        else:
            try:
                vf = np.vdot(self.v,f)
            except ValueError:
                self.v = np.zeros((len(f), 3))
                vf = np.vdot(self.v,f)
            if vf > 0.0:
                self.v = ((1.0 - self.a) * self.v + 
                          (self.a * vunit(f) * vmag(self.v)))
                if self.Nsteps > self.Nmin:
                    self.dt = min(self.dt * self.finc, self.dtmax)
                    self.a *= self.fa
                self.Nsteps += 1
            else:
                self.v *= 0.0
                self.a = self.astart
                self.dt *= self.fdec
                self.Nsteps = 0

        self.v += self.dt * f
        dr = self.dt * self.v
        dr = self.determine_step(dr)

        return r + dr.reshape(r.shape)

    
    def determine_step(self, dr):
        steplengths = (dr**2).sum(1)**0.5
        maxsteplength = np.max(steplengths)

        if maxsteplength >= self.maxstep:
            dr *= self.maxstep/maxsteplength

        return dr
