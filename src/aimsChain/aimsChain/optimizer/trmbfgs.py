import numpy as np
import cPickle as cp
from numpy.linalg import eigh, solve

from aimsChain.utility import vmag,vunit

class trm(object):
    def __init__(self, restart="hess",
                 maxstep=0.04,alpha=70):
        self.restart = restart
        self.maxstep = maxstep
        self.inistep=maxstep
        self.alpha = alpha
        self.epsilon = 1e-10
        self.writelog = False


    def log(self, str_in=""):
        if not self.writelog:
            return
        with open('trm.log','a') as log:
            log.write("Iteration[%04d]: " % self.iteration)
            log.write(str_in + '\n')


    def initialize(self):
        self.H = None
        self.r0 = None
        self.f0 = None
        self.delta = self.inistep
        self.iteration = 0


    def load(self):

        import os.path as path

        if path.isfile(self.restart):
            with open(self.restart, 'r') as hess:
                try:
                    self.H, self.r0, self.f0,self.delta,self.direc,self.trmpred,self.iteration  = cp.load(hess)
                except:
                    self.initialize()


    def dump(self):
        """ Dump necessary values for future reference."""

        with open(self.restart, 'w') as hess:
            cp.dump((self.H, self.r0, self.f0, self.delta, self.direc, self.trmpred, self.iteration), hess)


    def step(self,r,f):
        """ Take a single trust radius step."""

        r = np.array(r)
        f = np.array(f).reshape(-1)
        self.update(r,f)
        D,S = eigh(self.H)
        StF = np.dot(f,S)
        
        #This is from the Fortran version, but is  useless for us.
        #The original purpose was to report details about hessian
        #As well as deal with zero modes & machine epsilon 
        #But for a Chain these makes no sense, since we are not
        #interested in the details of the Hessian
        """
        trm_examine_Hess(D, StF)
        """
        self.log("steping the nodes")
        self.log("Delta: " + str(self.delta))
        (self.trmpred,StDX) = self.trm_TRM(D, StF, self.delta)
        dr = np.dot(S,StDX)
        #store this as direction for next iteration
        self.direc = dr
        dr = dr.reshape(-1,3)
        dr = self.determine_step(dr)
        dr = np.reshape(dr, np.shape(r))        
        self.r0 = r.flatten()
        self.f0 = f.flatten()
        self.iteration += 1

        return r+dr


    def update(self,r,f):
        r = r.flatten()
        f = f.flatten()

        if self.H is None or len(r) != len(self.r0):
            self.H = np.eye(len(f)) * self.alpha
            self.r0 = r
            self.f0 = f
            return
        r0 = self.r0
        f0 = self.f0        
        sk = r.reshape(-1) - r0.reshape(-1)
        yk = f0.reshape(-1) - f.reshape(-1) #f is negative gradient
        if np.abs(sk).max() < 1e-7:
            #same config
            return

        Hsk = np.dot(self.H, sk)
        self.H = (self.H + np.outer(yk,yk)/np.dot(sk,yk) 
                  - np.outer(Hsk,Hsk)/np.dot(sk,Hsk))
        quality = np.dot(vunit(self.direc), vunit(f))
        self.log("quality: " + str(quality))
        if quality > 0.5:
            self.delta *= 1.1
        else:
            self.delta *= 0.5


    def trm_TRM(self, D, StF, Delta):
        StDX = np.zeros(np.shape(D))
        for i,d in enumerate(D):
            if d > 0:
                StDX[i] = StF[i]/d
        minD = np.nanmin([np.nanmin(D), 
                         np.nanmin(np.clip(StF,0.0,np.inf))])
        
        #If matrix is positive and definite
        #And expected step size is less than delta
        #Then simply take it
        self.log("minD: " + str(minD))
        if (minD >= 0.0) and (vmag(StDX) < Delta):
            self.log("direct trival step")
#        if True:
            return (False,StDX)

        lambda_lower = np.nanmax([0.0, -minD])
        lambda_upper = 1e30
        Lambda = np.nanmin([lambda_lower + 0.5, 0.5*(lambda_upper+lambda_lower)])
        
        for _ in range(100):
            StDX = StF / (D+Lambda)
            g = np.sum(StDX**2) - Delta**2
            dg = -2.0 * np.sum(StDX**2/(D+Lambda))
            if (np.abs(g/dg) < self.epsilon) or (np.abs(g) < 1e-13):
                return (True,StDX)

            if g < 0.0:
                lambda_upper = min(Lambda, lambda_upper)
            else:
                lambda_lower = max(Lambda, lambda_lower)

            if (dg>0.0) or (lambda_lower > lambda_upper):
                self.log("Assertion failed")

            Lambda = Lambda - g/dg
            
            if(Lambda <= lambda_lower) or (Lambda >= lambda_upper):
                Lambda = 0.5* (lambda_lower + lambda_upper)
            if _ == 100:
                self.log("Newton's method failed")
            
        return (False,StDX)


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
