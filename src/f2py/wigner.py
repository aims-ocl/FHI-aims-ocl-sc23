#!/usr/bin/python

"""Test the Wigner rotation matrices."""


import sys

import numpy as np

import ylm

import _wigner_rot

def wigner_small_d(L, Theta):
    """Get Wigner small-d rotation matrix.

    
    Example (prepare a rotation and a vector):
    >>> import numpy.random
    >>> Theta = np.pi * np.random.random()
    >>> L = np.random.random_integers(10)
    >>> Tmat = np.array([[ np.cos(Theta), 0, np.sin(Theta)],
    ...                  [ 0,             1, 0],
    ...                  [-np.sin(Theta), 0, np.cos(Theta)]])
    ...
    >>> vec = np.random.random(3)
    >>> d = wigner_small_d(L, Theta)

    The rotation matrix should of course be orthgonal:
    >>> np.allclose(np.dot(d, d.T), np.identity(2*L+1))
    True

    Application (transposed) should give the same as Tmat on the argument:
    >>> Ys_vec = ylm.Y(vec, L, np.arange(-L, L+1))
    >>> Ys_Tvec = ylm.Y(np.dot(Tmat, vec), L, np.arange(-L, L+1))
    >>> dTYs_vec = np.dot(d.T, Ys_vec)
    >>> np.allclose(Ys_Tvec, dTYs_vec)
    True
    """
    return _wigner_rot.wigner_small_d(L, Theta)



def _test():
    import doctest
    print __doc__.rstrip()
    print "... if silent, no problems have occured ..."
    print
    print "Test docstrings"
    doctest.testmod()

if __name__ == "__main__":
    _test()
