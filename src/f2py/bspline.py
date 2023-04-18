#!/usr/bin/python

import sys

import numpy as np
from numpy import newaxis as NA

import ylm

import _bspline

periodic_bspline = _bspline.bspline.cubic_bspline_periodic
fundamental_bspline = _bspline.bspline.cubic_bspline_fundamental
notaknot_bspline = _bspline.bspline.cubic_bspline_notaknot
natural_bspline = _bspline.bspline.cubic_bspline_natural

class ScaledBspline(object):
    """B-spline with arbitrary domain.

    >>> import numpy.random
    >>> n = 10; n_func = 3; rr_shp = (7, 8)
    >>> spl_param = np.random.random((n_func, n+2))
    >>> spl = ScaledBspline(spl_param)
    >>> rr = np.random.random(rr_shp)
    >>> val = spl(rr)
    >>> val.shape ==  (n_func,) + rr_shp
    True
    """
    def __init__(self, spl_param, a=0., b=1.):
        self.spl_param, self.shp, self.n_func \
            = ylm.array_of_vectors(spl_param, vec_dim=None)
        self.n = np.size(self.spl_param, -1) - 2
        self.a = a
        self.b = b
        assert np.shape(self.spl_param) == (self.n_func, self.n + 2)

    def __call__(self, r, deriv=0, if_outside=None):
        rs, ret_shp, n_r = ylm.array_of_scalars(r)
        if if_outside is None:
            if np.any(np.logical_or(rs < self.a-1e-11, self.b+1e-11 < rs)):
                raise ValueError("Outlier")
        bval = _bspline.val_vec_bspline(self.a, self.b, rs, self.spl_param.T,
                                        deriv=deriv, if_outside=if_outside).T
        assert np.shape(bval) == (n_r, self.n_func)
        if if_outside is not None:
            bval = np.where(np.logical_and(self.a-1e-11 < rs[:,NA],
                                           rs[:,NA] < self.b+1e-11),
                            bval,
                            if_outside)
        assert np.shape(bval) == (n_r, self.n_func)
        return np.reshape(bval.T, self.shp + ret_shp)


def spline_array(ff, a=0., b=1.,
                 bspline_type='notaknot', dfa=None, dfb=None):
    """Calculate a B-spline with arbitrary domain.

    The function values ff should be given on
    tt = np.linspace(a, b, len(ff), endpoint=True).

    Example:
    >>> a, b = 0, 2*np.pi
    >>> tt = np.linspace(a, b, 201, endpoint=True)
    >>> ff = np.cos(tt)
    >>> spl = spline_array(ff, a, b, 'fundamental', dfa=0., dfb=0.)
    >>> ttt = np.linspace(a, b, 51, endpoint=True)
    >>> np.allclose(spl(ttt), np.cos(ttt))
    True
    >>> np.allclose(spl(ttt, deriv=1), -np.sin(ttt), atol=1e-6)
    True
    """
    ff, shp, n_func = ylm.array_of_vectors(ff, vec_dim=None)
    n = np.size(ff, -1)
    if bspline_type == 'notaknot':
        spl_param = notaknot_bspline(ff.T).T
    elif bspline_type == 'natural':
        spl_param = natural_bspline(ff.T).T
    elif bspline_type == 'fundamental':
        dfa, _, _ = ylm.array_of_scalars(dfa)
        dfb, _, _ = ylm.array_of_scalars(dfb)
        spl_param = fundamental_bspline(ff.T, dfa, dfb).T
    else:
        raise ValueError("Invalid spline type %s" % bspline_type)
    assert np.shape(spl_param) == (n_func, n + 2)
    param_shp = shp + (n+2,)
    spl_param = np.reshape(spl_param, param_shp)
    return ScaledBspline(spl_param, a, b)


class PeriodicBspline(object):
    """Periodic B-spline with arbitrary domain.

    >>> import numpy.random
    >>> n = 10; n_func1 = 2; n_func2 = 3; rr_shp = (7, 8)
    >>> spl_param = np.random.random((n_func1, n_func2, n+2))
    >>> spl = PeriodicBspline(spl_param)
    >>> rr = np.random.random(rr_shp)
    >>> val = spl(rr)
    >>> val.shape == (n_func1, n_func2) + rr_shp
    True
    >>> np.allclose(val, spl(rr + 1.))   # Periodicity
    True
    """
    def __init__(self, spl_param, a=0., b=1.):
        self.spl_param, self.shp, self.n_func \
            = ylm.array_of_vectors(spl_param, vec_dim=None)
        self.n = np.size(self.spl_param, -1) - 2
        self.a = a
        self.b = b
        assert np.shape(self.spl_param) == (self.n_func, self.n + 2)

    def __call__(self, r, deriv=0):
        rs, ret_shp, n_r = ylm.array_of_scalars(r)
        bval = _bspline.val_vec_bspline(self.a, self.b, rs,
                                        self.spl_param.T, deriv=deriv,
                                        is_periodic=True)
        assert np.shape(bval) == (self.n_func, n_r), bval.shape
        return np.reshape(bval, self.shp + ret_shp)


def spline_periodic_array(ff, a=0., b=1.,
                 bspline_type='notaknot', dfa=None, dfb=None):
    """Calculate a periodic B-spline with arbitrary domain.

    The function values ff should be given on
    tt = np.linspace(a, b, len(ff), endpoint=False).

    Example:
    >>> a, b = 0, 2*np.pi
    >>> tt = np.linspace(a, b, 200, endpoint=False)
    >>> ff = np.sin(tt)
    >>> spl = spline_periodic_array(ff, a, b)
    >>> ttt = np.linspace(a, b, 400, endpoint=True)
    >>> np.allclose(spl(ttt), np.sin(ttt))
    True
    >>> np.allclose(spl(ttt, deriv=1), np.cos(ttt), atol=1e-6)
    True
    """
    ff, shp, n_func = ylm.array_of_vectors(ff, vec_dim=None)
    n = np.size(ff, -1)
    spl_param = periodic_bspline(ff.T).T
    assert np.shape(spl_param) == (n_func, n + 2)
    param_shp = shp + (n+2,)
    spl_param = np.reshape(spl_param, param_shp)
    return PeriodicBspline(spl_param, a, b)


def get_subspline(scl_spl, i0, iN1, drop_fac):
    """Get a (possibly coarsend) subspline and its error.

    Example:
    >>> a, b = 0, 2*np.pi
    >>> tt = np.linspace(a, b, 201, endpoint=True)
    >>> ff = np.cos(tt)
    >>> spl = spline_array(ff, a, b, 'fundamental', dfa=0., dfb=0.)
    >>> dspl, max_err = get_subspline(spl, 30, 170, 4)

    The splines are still rather accurate:
    >>> max_err < 1e-6
    True
    >>> ttt = np.linspace(0.5*np.pi, 1.5*np.pi, 140)
    >>> np.allclose(dspl(ttt), np.cos(ttt), atol=1e-7)
    True

    The return parameter max_err is actually sensible:
    >>> np.abs(np.log10(np.max(np.abs(dspl(ttt) - np.cos(ttt)))) -
    ...        np.log10(max_err)) < 1.
    True
    """
    assert scl_spl.n_func == 1
    spl_param_in = scl_spl.spl_param.flatten()
    get_sub = _bspline.bspline.get_subspline
    n_out = (iN1 - i0) / drop_fac + 1
    iN1 = i0 + drop_fac * (n_out-1)

    spl_param, max_err = get_sub(i0+1, iN1+1, drop_fac, spl_param_in, n_out, 0.)
    spl_param = spl_param.flatten()

    a_out = i0  * (scl_spl.b - scl_spl.a) / (scl_spl.n-1)
    b_out = iN1 * (scl_spl.b - scl_spl.a) / (scl_spl.n-1)
    scl_out = ScaledBspline(spl_param, a_out, b_out)
    return scl_out, max_err


def _test():
    import doctest
    print "Script to test numerical grids."
    print "... if silent, no problems have occured ..."
    print
    print "Test docstrings"
    doctest.testmod()


if __name__ == "__main__":
    _test()
