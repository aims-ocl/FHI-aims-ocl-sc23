#!/usr/bin/python3
"""
This module provides a set of basic comparison functions for the evaluation of
the test parameters extracted from the regression tests.
To use one of these functions in your TestCase, simply define your TestParameter
like this:

    <TestParameter name="dummy" comparison="comparisons.string_equal/>

Where the first part of the argument is the module (file) name and the second
is the name of the function you want to call. For sophisticated comparisons,
there are also functions with additional parameters available, these are called
with the following syntax;

    <TestParameter name="dummy" comparison="comparisons.float_absdiff(0.01)/>

where the additional arguments are in parantheses, as the maximal absolute
difference between two floating point numbers in the example.

Feel free to add your own comparison functions if the standard set does not
suffice for your test cases.

Each function should return two values: a boolean to indicate if the equality
is satisfied or not and a string with additional information about the test.
The second return value may be the empty string.
"""

def abs_diff(tolerance, precision=3):
    """
    returns true, if the two objects can be converted to floats and are their
    absolute difference is smaller than the given tolerance. The optional
    precision argument determines the number of digits behind the decimal dot.
    """
    tolerance = float(tolerance)
    prec = int(precision)
    msg = r"abs diff = {:.%iE} (tol = {:.%iE})"%(prec, prec)
    msg_file = r"abs diff = {:.%iE} (tol = {:.%iE}) [file={}]"%(prec, prec)
    def inner_float_absdiff(file, obj1, obj2):
        """inner function of closure"""
        diff = abs(float(obj1)-float(obj2))
        if file:
            return diff <= tolerance, msg_file.format(diff, tolerance, file)
        else:
            return diff <= tolerance, msg.format(diff, tolerance)
    return inner_float_absdiff

def absvalue_diff(tolerance, precision=3):
    """
    returns true, if the two objects can be converted to floats and the
    difference of their absolute values is smaller than the given tolerance.
    The optional precision argument determines the number of digits behind the
    decimal dot.
    """
    tolerance = float(tolerance)
    prec = int(precision)
    msg = r"absvalue diff = {:.%iE} (tol = {:.%iE})"%(prec, prec)
    msg_file = r"absvalue diff = {:.%iE} (tol = {:.%iE}) [file={}]"%(prec, prec)
    def inner_float_absdiff(file, obj1, obj2):
        """inner function of closure"""
        diff = abs(abs(float(obj1))-abs(float(obj2)))
        if file:
            return diff <= tolerance, msg_file.format(diff, tolerance, file)
        else:
            return diff <= tolerance, msg.format(diff, tolerance)
    return inner_float_absdiff

def rel_diff(tolerance, precision=3, reference=0):
    """
    returns true, if the two objects can be converted to floats and are their
    relative difference is smaller than the given tolerance (in %). The optional
    precision argument determines the number of digits behind the decimal dot.
    (reference=0 divdes by ref, 1 by test)
    """
    tolerance = float(tolerance)
    prec = int(precision)
    msg = r"rel diff = {:.%iF}%% (tol = {:.%iF}%%)"%(prec, prec)
    msg_file = r"rel diff = {:.%iF}%% (tol = {:.%iF}%%) [file={}]"%(prec, prec)
    def inner_float_absdiff(file, obj1, obj2):
        """inner function of closure"""
        ref = abs(float(obj1)) if reference == 0 else abs(float(obj2))
        diff = abs(float(obj1)-float(obj2)) / ref * 100
        if file:
            return diff <= tolerance, msg_file.format(diff, tolerance, file)
        else:
            return diff <= tolerance, msg.format(diff, tolerance)
    return inner_float_absdiff
