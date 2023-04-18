#!/usr/bin/python

"""Script for providing spherical harmonics."""


import sys

import numpy as np
import numpy.random

import _ylm

try:
    import scipy.special
    import scipy.integrate
    use_scipy = True
except ImportError:
    use_scipy = False


def find_index(l, m):
    """Get index in ylm(:) for given l, m.

    >>> vals = [(0,0, 0), (1,-1, 1), (1,0, 2), (1,1, 3), (2,0, 6)]
    >>> all([find_index(l, m) == lm for l, m, lm in vals])
    True
    """
    l = np.asarray(l)
    m = np.asarray(m)
    if np.any(np.abs(m) > l): raise ValueError("|m| > l")
    return l**2 + (m + l)

def index_to_lm(lm):
    """Inverse of find_index().

    >>> lm = np.random.random_integers(500)
    >>> lm == find_index(*index_to_lm(lm))
    True
    """
    l = int(np.floor(np.sqrt(lm)))
    m = lm - l**2 - l
    return l, m


def array_of_scalars(scals):
    """Regularize shape of scalars to (n).

    Example:
    >>> scals, ret_shp, n_scals = array_of_scalars([[1,2,3], [4,5,6]])
    >>> scals.shape
    (6,)
    >>> ret_shp
    (2, 3)
    """
    scals = np.asarray(scals)
    shp = np.shape(scals)
    sze = np.size(scals)
    scals = np.reshape(scals, (sze,))
    return scals, shp, sze

def array_of_vectors(rvecs, vec_dim=3):
    """Regularize shape of rvecs to (n, vec_dim).

    Example:
    >>> rvecs, ret_shp, n_vecs = array_of_vectors([1,2,3])
    >>> rvecs.shape
    (1, 3)
    >>> rvecs, ret_shp, n_vecs = array_of_vectors([[[1,2,3], [4,5,6]],
    ...                                           [[7,8,9],[10,11,12]]])
    >>> rvecs.shape
    (4, 3)
    >>> ret_shp
    (2, 2)
    """
    rvecs = np.asarray(rvecs)
    shp = np.shape(rvecs)
    if vec_dim is None:
        vec_dim = shp[-1]
    elif shp[-1] != vec_dim:
        raise ValueError("Need %i-vectors for rvec" % vec_dim)
    ret_shp = shp[:-1]
    n_vec = np.product(ret_shp)
    rvecs = np.reshape(rvecs, (n_vec, vec_dim))
    return rvecs, ret_shp, n_vec

def to_polar(rvecs):
    """Convert kartesian to spherical coordinates.

    Example:
    >>> np.allclose(to_polar([1,0,0]), (1., np.pi/2., 0.))
    True
    >>> np.allclose(to_polar([[0,1,0], [1e-10,0,2]]),
    ...             ([1.,2.], [np.pi/2.,0.], [np.pi/2., 0.]))
    True
    """
    rvecs, ret_shp, _ = array_of_vectors(rvecs)
    xs, ys, zs = rvecs.T
    rho_sq = xs**2 + ys**2
    rho = np.sqrt(rho_sq)
    r = np.sqrt(rho_sq + zs**2)
    phi = np.arctan2(ys, xs)
    theta = np.arctan2(rho, zs)
    r, theta, phi = [np.reshape(xxx, ret_shp) for xxx in (r, theta, phi)]
    return r, theta, phi

def to_kartesian(theta, phi):
    """Return 3-component vector on S^2.

    >>> np.allclose(to_kartesian([0., np.pi/2], [0., np.pi]),
    ...             [[0,0,1], [-1,0,0]])
    True
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    rvec = np.array([x,y,z])                         # x,y,z by 0-index
    return np.rollaxis(rvec, 0, np.rank(rvec)) # put 0-index last


if use_scipy:
    def _integrate_Omega(func):
        """\\int_{S^2} func(theta, phi) * sin(theta) dtheta dphi.

        >>> np.allclose(_integrate_Omega(lambda theta, phi: 1.), 4*np.pi)
        True
        """
        def integrand(theta, phi):
            return np.sin(theta) * func(theta, phi)
        def low_theta(phi): return 0.
        def high_theta(phi): return np.pi
        I, Ierr = scipy.integrate.dblquad(integrand, 0., 2*np.pi,
                                          low_theta, high_theta)
        assert Ierr < 1e-5
        return I

def Y(vec, l, m):
    """Get single value of Y_lm as defined in ylm.f.

    >>> def YY(theta, phi):
    ...    vec = to_kartesian(theta, phi)
    ...    return Y(vec, l1, m1)*np.conj(Y(vec, l2, m2))
    >>> theta = np.pi * numpy.random.random()
    >>> phi = 2*np.pi * numpy.random.random()
    >>> vec = to_kartesian(theta, phi)
    >>> if use_scipy:
    ...    for l1 in xrange(3):
    ...       for m1 in xrange(-l1, l1+1):
    ...          # Consistent with scipy.special.sph_harm():
    ...          sp_res = scipy.special.sph_harm(m1, l1, phi, theta)
    ...          assert np.allclose(Y(vec, l1, m1), sp_res)
    ...          # Orthonormal:
    ...          l2 = numpy.random.random_integers(0, l1)
    ...          m2 = numpy.random.random_integers(-l2, l2)
    ...          I = _integrate_Omega(YY)
    ...          expect = 0. if (l1, m1) != (l2, m2) else 1.
    ...          assert np.allclose(I, expect), (l1, m1, l2, m2)
    """
    lm, lm_shp, n_lm = array_of_scalars(find_index(l, m))
    vecs, vec_shp, n_vecs = array_of_vectors(vec)
    ret_shp = vec_shp + lm_shp
    ylms = np.empty((n_vecs, n_lm), dtype=np.complex_)
    for i, rvec in enumerate(vecs):
        these_ylms = _ylm.ylm(rvec, np.max(l))
        ylms[i,:] = these_ylms[lm]
    return np.reshape(ylms, ret_shp)
        

def Yr_cache(vec, lmax):
    """Get all Yrs up to lmax in a find_index() like structure."""
    return _ylm.ylm_real(vec, lmax)


def Yr(vec, l, m):
    """Get single value of Y^r_lm as defined in ylm_real.f.

    >>> def YYr(theta, phi):
    ...    vec = to_kartesian(theta, phi)
    ...    return Yr(vec, l1, m1) * Yr(vec, l2, m2)
    >>> for l1 in xrange(3):
    ...    for m1 in xrange(-l1, l1+1):
    ...       l2 = numpy.random.random_integers(0, l1)
    ...       m2 = numpy.random.random_integers(-l2, l2)
    ...       I = _integrate_Omega(YYr)
    ...       expect = 0. if (l1, m1) != (l2, m2) else 1.
    ...       assert np.allclose(I, expect), (I, expect, (l1, m1), (l2, m2))
    """
    lmax = np.max(l)
    lm, lm_shp, n_lm = array_of_scalars(find_index(l, m))
    vecs, vec_shp, n_vecs = array_of_vectors(vec)
    ret_shp = vec_shp + lm_shp
    ylms = np.empty((n_vecs, n_lm), dtype=np.float_)
    for i, rvec in enumerate(vecs):
        these_ylms = _ylm.ylm_real(rvec, lmax)
        ylms[i,:] = these_ylms[lm]
    return np.reshape(ylms, ret_shp)

def Yr_py(vec, l, m, from_which=0):
    """Get single value of Y^r_lm from complex Y_lm in ylm.f.

    >>> for l in xrange(10):
    ...    for m in xrange(-l, l+1):
    ...       vec = np.random.randn(3)
    ...       for from_which in xrange(2):
    ...           from_aims = Yr(vec, l, m)
    ...           this_one = Yr_py(vec, l, m, from_which)
    ...           assert np.allclose(from_aims, this_one), (
    ...              str(vec), (l, m), from_which, from_aims, this_one)
    """
    if from_which == 0:   # from Y_{l, m} & Y_{l, -|m|}
        Yplus = Y(vec, l, m)
        Yminus = Y(vec, l, -m)
        if m > 0:
            cval = np.sqrt(0.5) * (Yplus + (-1)**m * Yminus)
        elif m == 0:
            cval = Yplus
        else:
            cval = 1j * np.sqrt(0.5) * (Yplus - (-1)**m * Yminus)
        assert np.allclose(cval.imag, 0.)
        return cval.real
    elif from_which == 1:
        Ylm = Y(vec, l, m)
        if m > 0:
            return np.sqrt(2.) * np.real(Ylm)
        elif m == 0:   # from Y_{l, m} & Y_{l, -|m|}
            return np.real(Ylm)
        else:
            return - np.sqrt(2.) * np.imag(Ylm)

def _test():
    import doctest
    print __doc__.rstrip()
    print "... if silent, no problems have occured ..."
    print
    print "Test docstrings"
    doctest.testmod()
    
if __name__ == "__main__":
    _test()
