#!/usr/bin/python

"""Script to test some relations among spherical harmonics."""

import sys

import numpy as np
import numpy.random

import ylm

import _triple_y

def valid_3jm(l1, l2, l3, m1, m2, m3):
    """Returns True if finite 3jm, False if zero 3jm and raises for invalid."""
    def isint(x):
        return np.allclose(round(x), x)
    if not all([isint(x) for x in [2*l1, 2*l2, 2*l3, 2*m1, 2*m2, 2*m3]]):
        raise ValueError("Some value is not (half-)integer")
    if not all([isint(x) for x in [l1-m1, l2-m2, l3-m3]]):
        raise ValueError("Half-integer mismatch between l and m")
    if not all([abs(m) <= l for l, m in [(l1, m1), (l2, m2), (l3, m3)]]):
        raise ValueError("Some m out of range")
    if not isint(l1 + l2 + l3):
        return False
    if int(round(m1 + m2 + m3)) != 0:
        return False
    if l1 < abs(l2 - l3) or l1 > l2 + l3:
        return False
    if (l1 + l2 + l3) % 2 == 1:
        return False
    return True
    
def valid_3jm_real(l1, l2, l3, m1, m2, m3):
    """Same as valid_3jm, but for real-valued spherical harmonics."""
    def isint(x):
        return np.allclose(round(x), x)
    if not all([isint(x) for x in [2*l1, 2*l2, 2*l3, 2*m1, 2*m2, 2*m3]]):
        raise ValueError("Some value is not (half-)integer")
    if not all([isint(x) for x in [l1-m1, l2-m2, l3-m3]]):
        raise ValueError("Half-integer mismatch between l and m")
    if not all([abs(m) <= l for l, m in [(l1, m1), (l2, m2), (l3, m3)]]):
        raise ValueError("Some m out of range")
    if not isint(l1 + l2 + l3):
        return False
    if l1 < abs(l2 - l3) or l1 > l2 + l3:
        return False
    if (l1 + l2 + l3) % 2 == 1:
        return False
    ma, mb, mc = sorted([abs(m) for m in m1, m2, m3])
    if ma + mb != mc:
        return False   # Not combinable to 0 otherwise.
    if len([None for m in (m1, m2, m3) if m < 0]) % 2 == 1:
        return False   # would be purely imaginary => 0.
    return True


def three_jm_from_drcjj(l1, l2, l3, m1, m2, m3):
    """Return 3jm symbol calculated with drcjj.f"""
    if not valid_3jm(l1, l2, l3, m1, m2, m3):
        return 0.
    l1min, l1max, thrcof, ier = _triple_y.drc3jj(l2, l3, m2, m3)
    assert ier == 0, ier
    return thrcof[int(round(l1 - l1min))]

def three_jm_from_drcjm(l1, l2, l3, m1, m2, m3):
    """Return 3jm symbol calculated with drcjm.f"""
    if not valid_3jm(l1, l2, l3, m1, m2, m3):
        return 0.
    m2min, m2max, thrcof, ier = _triple_y.drc3jm(l1, l2, l3, m1)
    assert ier == 0, ier
    return thrcof[int(round(m2 - m2min))]

def C3_from_3jm(l1, l2, l3, m1, m2, m3, f_3jm=three_jm_from_drcjj):
    """Calculate triple Y integral \int dOmega Y_l1m1 Y_l2m2 Y_l3m3."""
    presq = (2*l1+1) * (2*l2+1) * (2*l3+1) / (4*np.pi)
    return (np.sqrt(presq) *
            f_3jm(l1, l2, l3, 0, 0, 0) *
            f_3jm(l1, l2, l3, m1, m2, m3))

def C3_from_fortran(l1, l2, l3, m1, m2, m3):
    """Calculate triple Y integral \int dOmega Y_l1m1 Y_l2m2 Y_l3m3."""
    if not valid_3jm(l1, l2, l3, m1, m2, m3):
        return 0.
    triple_Y = _triple_y.triple_y.triple_y_cmplx(l1, l2, l3)
    return triple_Y[l1 + m1, l2 + m2]

def C3r_from_3jm(l1, l2, l3, m1, m2, m3):
    """Calculate triple Y integral \int dOmega Y^r_l1m1 Y^r_l2m2 Y^r_l3m3."""
    if not valid_3jm_real(l1, l2, l3, m1, m2, m3):
        return 0.
    # sort them:
    all_lm = [(l1, m1), (l2, m2), (l3, m3)]
    def key(lm):
        l, m = lm
        return - abs(m)
    (l1, m1), (l2, m2), (l3, m3) = sorted(all_lm, key=key)
    if abs(m1) != abs(m2) + abs(m3):   # Only because |m1| is largest.
        return 0.
    if m2 < 0 and m3 < 0:
        sign = -1
    else:
        sign = 1
    for m in (m1, m2, m3):
        if m < 0: sign *= (-1)**m
    if m3 == 0:
        w = 1.
    else:
        w = 1. / np.sqrt(2.)
    C3 = C3_from_3jm(l1, l2, l3, -abs(m1), abs(m2), abs(m3))
    return sign * (-1)**m1 * w * C3

def C3r_from_fortran(l1, l2, l3, m1, m2, m3):
    """Calculate triple Y integral \int dOmega Y^r_l1m1 Y^r_l2m2 Y^r_l3m3."""
    if not valid_3jm_real(l1, l2, l3, m1, m2, m3):
        return 0.
    triple_Yr, Ms = _triple_y.triple_y.triple_y_real(l1, l2, l3)
    m3s = list(Ms[l1 + m1, l2 + m2])
    if m3 in m3s:
        M_type = m3s.index(m3)
        return triple_Yr[l1 + m1, l2 + m2, M_type]
    else:
        return 0.

def C3_from_integral(l1, l2, l3, m1, m2, m3):
    """Numerically integrate \int dOmega Y_l1m1 Y_l2m2 Y_l3m3."""
    def YYY(theta, phi):
        vec = ylm.to_kartesian(theta, phi)
        return ylm.Y(vec, l1, m1) * ylm.Y(vec, l2, m2) * ylm.Y(vec, l3, m3)
    return ylm._integrate_Omega(YYY)


def C3r_from_integral(l1, l2, l3, m1, m2, m3):
    """Numerically integrate \int dOmega Y^r_l1m1 Y^r_l2m2 Y^r_l3m3."""
    def YYY(theta, phi):
        vec = ylm.to_kartesian(theta, phi)
        return ylm.Yr(vec, l1, m1) * ylm.Yr(vec, l2, m2) * ylm.Yr(vec, l3, m3)
    return ylm._integrate_Omega(YYY)

def _test():
    import doctest
    print "Script to test relations among spherical harmonics."
    print "... if silent, no problems have occured ..."
    print
    print "Test docstrings"
    doctest.testmod()
    print "Test C3 (complex Y)"
    for l1 in xrange(5):
        l2 = np.random.random_integers(0, 5)
        l3 = np.random.random_integers(abs(l1-l2), l1+l2)
        m1 = np.random.random_integers(-l1, l1)
        m2 = np.random.random_integers(-l2, l2)
        m3 = - (m1 + m2)
        if abs(m3) > l3:
            m3 = np.random.random_integers(-l3, l3)
        C1 = C3_from_3jm(l1, l2, l3, m1, m2, m3,
                         f_3jm=three_jm_from_drcjj)
        C2 = C3_from_3jm(l1, l2, l3, m1, m2, m3,
                         f_3jm=three_jm_from_drcjm)
        C3 = C3_from_integral(l1, l2, l3, m1, m2, m3)
        C4 = C3_from_fortran(l1, l2, l3, m1, m2, m3)
        assert (np.allclose(C1, C2) and np.allclose(C1, C3) and
                np.allclose(C1, C4)), (
                (l1, m1), (l2, m2), (l3, m3), C1, C2, C3, C4)

    print "Test C3r (real Y^r)"
    l1 = np.random.random_integers(0, 5)
    l2 = np.random.random_integers(0, 10)
    l3 = np.random.random_integers(abs(l1-l2), l1+l2)
    m1 = np.random.random_integers(-l1, l1)
    m2 = np.random.random_integers(-l2, l2)
    m3 = np.random.random_integers(-l3, l3)
    C1 = C3r_from_3jm(l1, l2, l3, m1, m2, m3)
    C2 = C3r_from_integral(l1, l2, l3, m1, m2, m3)
    C3 = C3r_from_fortran(l1, l2, l3, m1, m2, m3)
    assert np.allclose(C1, C2) and np.allclose(C1, C3), (
        (l1, m1), (l2, m2), (l3, m3), C1, C2, C3)


if __name__ == "__main__":
    _test()
