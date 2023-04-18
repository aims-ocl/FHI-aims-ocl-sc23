import numpy as np

class Optimize(object):
    def __init__(self):
        self.maxstep = 0.04

    
    def determine_step(self, dr):
        """Determine step to take according to maxstep
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        
        steplengths = (dr**2).sum(1)**0.5
        maxsteplength = np.max(steplengths)
        factor = 1
        if maxsteplength >= self.maxstep:
            factor = self.maxstep / maxsteplength
            dr *= factor
        
        return dr, factor

class FDOptimize(Optimize):
    def __init__(self):
        self.maxstep = 0.04
        self.finitie_diff = False
