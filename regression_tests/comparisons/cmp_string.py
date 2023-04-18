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

def stripped_equal(file, obj1, obj2):
    """returns true, if the stripped input strings are equal"""
    return obj1.strip() == obj2.strip(), "string comparison"

def equal(file, obj1, obj2):
    """returns true, if the input strings are exactly equal"""
    return obj1 == obj2, "string comparison"
