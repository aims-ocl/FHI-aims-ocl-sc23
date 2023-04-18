"""
Provides various frequently used functions used in vector manipulation.
"""

import numpy as np

def vmag(v):
    """ Returns the magnitude of v computed with a flattened vector"""

    return np.sqrt(np.vdot(v,v))


def vproj(v1,v2):
    """ Returns the projection of v1 onto v2."""

    mag2 = vmag(v2)
#    if mag2 == 0:
#        print "Can't project onto a zero vector\n"
#        return v1
    return np.vdot(v1,v2)/mag2 *v2


def vunitproj(v1,v2):
    """ Return the projection of v1 onto unit vector of v2."""

    mag2 = vmag(v2)
#    if mag2 == 0:
#        print "Can't project onto a zero vector\n"
#        return v1
    vunit2 = vunit(v2)
    return np.vdot(v1, vunit2) * vunit2


def vunit(v):
    """ Return the unit vector of v."""

    mag = vmag(v)
 #   if mag == 0:
 #       print "can't normalize a zero vector\n"
 #       return v

    return v / mag
