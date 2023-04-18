import numpy as np
import cPickle as cp
from aimsChain.utility import vunitproj,vunit

class newFIRE(object):
    def __init__(self, restart="fire.dat",
                 dt=0.0005, maxmove=0.04, dtmax=0.048, Nmin=6, finc=1.05, fdec=0.5,
                 astart=0.00, fa=1.0, a=0.00):

        self.dt = dt
        self.Nsteps = 0
        self.maxmove = maxmove
        self.dtmax = dtmax
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa
        self.a = a
        self.restart = restart
        self.r0 = None
        self.delta_r = []


    def initialize(self):
        self.v = None


    def load(self):
        """ Load saved velocity and dt from file."""

        import os.path as path

        if path.isfile(self.restart):
            with open(self.restart,'r') as save:
                self.v, self.dt, self.a, self.Nsteps, self.r0, self.delta_r = cp.load(save)


    def dump(self):
        """ Dump necessary values for future reference."""

        with open(self.restart, 'w') as save:
            cp.dump((self.v, self.dt, self.a, self.Nsteps, self.r0, self.delta_r), save)
      
 
    def step(self,r,f):
        r = np.array(r)
        f = np.reshape(f,(-1,3))

        if self.v is None:
            self.v = f
            self.r0 = r.reshape(-1,3)
        else:
            delta_r = r.reshape(-1,3) - self.r0
            self.delta_r.append(vunit(delta_r))
            delta_r = vunit(np.array(self.delta_r[-5:]).sum(0))
            vf = np.vdot(f,delta_r)
            if vf > 0.0:
                self.v = ((1.0-self.a) * f + 
                          (self.a * delta_r *
                           np.sqrt(np.vdot(f,f))))
#                self.v = ((1.0 - self.a) * self.v + 
#                          (self.a * vunit(delta_r) * 
#                           np.sqrt(np.vdot(self.v, self.v))))
                if self.Nsteps > self.Nmin:
                    self.dt = min(self.dt * self.finc, self.dtmax)
#                    self.a *= self.finc
#                    self.a *= self.fa
                self.Nsteps += 1
            else:
#                self.v = ((1.0-self.a) * f - 
#                          (self.a * vunit(delta_r) *
#                           np.sqrt(np.vdot(f,f))))
                self.v = f
                self.a = self.astart
                self.dt = max(self.dt* self.fdec, 0.000001)
                self.Nsteps = 0

#        self.v += self.dt * f
        dr = self.dt * self.v
#        dr,dt = self.determine_step(dr,self.dt)
        return r+dr.reshape(r.shape)


    def determine_step(self, dr, dt):
        steplengths = (dr**2).sum(1)**0.5
        maxsteplength = np.max(steplengths)
        if maxsteplength >= self.maxmove:
            dr *= self.maxmove / maxsteplength
            dt *= self.maxmove / maxsteplength
        return dr,dt
