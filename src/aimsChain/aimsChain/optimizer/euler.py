import numpy as np

class EULER(object):
   """ A basic forward Euler method."""

   def __init__(self, tstep = 0.01):
      self.dt = tstep
   
   #placeholder functions for uniformity
   def load(self):
      pass


   def dump(self):
      pass


   def initialize(self):
      pass


   def step(self, r, f):

      import numpy as np      

      new_pos = None      
      positions = np.array(r)
      oldshape = np.shape(positions)
      positions = np.reshape(positions, (-1,3))
      forces = np.array(f)
      forces = np.reshape(forces, (-1,3))      

      if positions.shape == forces.shape:
         new_pos = positions + self.dt*forces
         new_pos = np.reshape(new_pos, oldshape)

      return new_pos
      
