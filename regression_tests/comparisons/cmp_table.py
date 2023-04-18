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
from .abstractcomparison import AbstractComparison

class TableFloatDiff(AbstractComparison):
    """
    compares a N-column table extracted from a file by their absolute or
    relative diff
    """

    def __init__(self, columns, tolerance, precision=3, comments=None, #pylint: disable=R0913
                 headerrows=1, threshold=0.0, relative="False"):
        """
        mandatory constructor arguments:
            -columns: the number of columns this table has
            -tolerance: the threshold, the absolute differences of all table
                        cells must be smaller than ot
        optional constructor arguments:
            -precision: digits after the decimal dot in output
            -comments: an iterable containing all characters that indicate the
                       start of a comment line (without leading whitespace)
            -headerrows: the number of rows that consitute the table header,
                         this INCLUDES the line that matches the regex
            -threshold: skip value pairs that are smaller than this value
            - relative: if true, perform a relative difference comparison
        """
        super().__init__()
        self._tolerance = float(tolerance)
        self._cols = int(columns)
        self._texts = dict()
        self._comments = comments if comments else self._comments
        self._maxdiffs = dict()
        self._headerrows = int(headerrows)
        self._threshold = float(threshold)
        self._reldiff = relative.lower == "true"
        precision = int(precision)
        prefix = "{}/{} rows equal - "
        location = " in row {}, col {}"
        first = "first diff: "
        maxdiff = "max diff: "
        text = "{:.%iE} (abs tol = {:.%iE})"%(precision, precision)+location
        self._texts["abs_max"] = prefix + maxdiff + text
        self._texts["abs_first"] = prefix + first + text
        text = "{:.%iF}%% (rel tol = {:.%iF}%%)"%(precision, precision)+location
        self._texts["rel_max"] = prefix + maxdiff + text
        self._texts["rel_first"] = prefix + first + text

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _compare(self, obj1, obj2):
        if self._position["tag"] not in self._maxdiffs:
            data = {"line":-1, "column":-1, "diff":0.0, "failed":False,
                    "line_first":-1, "column_first":-1, "diff_first":0.0,
                    "rows":0, "rows_failed":0}
            self._maxdiffs[self._position["tag"]] = data
        if self._is_comment(obj1) and self._is_comment(obj2):
            return True, "comment lines"
        if self._position["line"] <= self._headerrows:
            return True, "header rows"
        diffdata = self._maxdiffs[self._position["tag"]]
        diffdata["rows"] += 1
        extract = lambda x: [float(y) for y in x.split()]
        data_ref, data_test = self._extract_data(obj1, obj2, extract)
        line_failed = False
        for i in range(self._cols):
            if i >= len(data_ref) or i >= len(data_test):
                raise ValueError
            data = data_ref[i], data_test[i]
            if all((abs(x) < self._threshold for x in data)):
                continue
            if self._reldiff:
                diff = abs(data_ref[i] - data_test[i]) / abs(data_ref[i]) * 100
            else:
                diff = abs(data_ref[i] - data_test[i])
            if diff > diffdata["diff"]:
                diffdata["diff"] = diff
                diffdata["line"] = self._position["line"]
                diffdata["column"] = i + 1
            if diff > self._tolerance and not diffdata["failed"]:
                diffdata["failed"] = True
                diffdata["diff_first"] = diff
                diffdata["line_first"] = self._position["line"]
                diffdata["column_first"] = i + 1
            if diff > self._tolerance and not line_failed:
                diffdata["rows_failed"] += 1
                line_failed = True
        return True, "interim analysis"

    def _statistics(self):
        """returns the number of successful and total row tests"""
        tag = self._position["tag"]
        total = self._maxdiffs[tag]["rows"]
        return total - self._maxdiffs[tag]["rows_failed"], total

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def largest_error(self):
        passed, total = self._statistics()
        data = self._maxdiffs[self._position["tag"]]
        if self._reldiff:
            text = self._texts["rel_max"]
        else:
            text = self._texts["abs_max"]
        state = False if data["failed"] else True
        return state, text.format(passed, total, data["diff"], self._tolerance,
                                  data["line"], data["column"])

    def first_error(self):
        passed, total = self._statistics()
        data = self._maxdiffs[self._position["tag"]]
        if data["failed"]:
            if self._reldiff:
                text = self._texts["rel_first"]
            else:
                text = self._texts["abs_first"]
            text = text.format(passed, total, data["diff_first"],
                               self._tolerance, data["line_first"],
                               data["column_first"])
            return False, text
        return self.largest_error()
