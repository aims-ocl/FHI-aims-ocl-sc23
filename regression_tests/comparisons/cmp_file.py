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

from collections import OrderedDict

class TableFileDiff(AbstractComparison):
    """
    compare a simple table file, either by absolute or relative differences,
    optionally a threshold criterion can be defined to ignore noise level values
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, columns, tolerance, precision=3, comments=None, #pylint: disable=R0913
                 threshold=0.0, relative="False"):
        """
        mandatory constructor arguments:
            -columns: number of columns this table file has
            -tolerance: RMSD of the two files must be smaller than this value
        optional constructor arguments:
            -precision: digits after decimal dot in output
            -comments: a string with all characters that mark comment lines
            -threshold: ignore value pairs below this (absolute) threshold
            -relative: if true, make a relative difference comparison
        """
        super().__init__()
        self._tolerance = float(tolerance)
        self._cols = int(columns)
        self._threshold = float(threshold)
        self._reldiff = relative.lower() == "true"
        self._comments = comments if comments != None else self._comments
        self._maxdiffs = OrderedDict()
        self._texts = dict()
        precision = int(precision)
        prefix = "{}/{} files equal - "
        location = " in '{}', line {}, col {}"
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
                    "line_first":-1, "column_first":-1, "diff_first":0.0}
            self._maxdiffs[self._position["tag"]] = data
        if self._is_comment(obj1) and self._is_comment(obj2):
            return True, "comment lines"
        diffdata = self._maxdiffs[self._position["tag"]]
        extract = lambda x: [float(y) for y in x.split()]
        data_ref, data_test = self._extract_data(obj1, obj2, extract)
        for i in range(self._cols):
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
                if not diffdata["failed"] and diff > self._tolerance:
                    diffdata["failed"] = True
                    diffdata["diff_first"] = diff
                    diffdata["line_first"] = self._position["line"]
                    diffdata["column_first"] = i + 1
        return True, "interim analysis"

    def _statistics(self):
        """returns the number of successful and total file tests"""
        failed = len([x for x in self._maxdiffs.items() if x[1]["failed"]])
        total = len(self._maxdiffs)
        return total-failed, total

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def largest_error(self):
        passed, total = self._statistics()
        keyfunc = lambda x: x[1]["diff"]
        file, data = sorted(self._maxdiffs.items(), key=keyfunc)[-1]
        if self._reldiff:
            text = self._texts["rel_max"]
        else:
            text = self._texts["abs_max"]
        state = False if data["diff"] >= self._tolerance else True
        return state, text.format(passed, total, data["diff"], self._tolerance,
                                  file, data["line"], data["column"])

    def first_error(self):
        passed, total = self._statistics()
        for file, data in self._maxdiffs.items():
            if data["failed"]:
                if self._reldiff:
                    text = self._texts["rel_first"]
                else:
                    text = self._texts["abs_first"]
                text = text.format(passed, total, data["diff_first"],
                                   self._tolerance, file, data["line_first"],
                                   data["column_first"])
                return False, text
        return self.largest_error()

class CubeFileRMSD(AbstractComparison):
    """compare a cubefile by means of an RMSD-analysis"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, tolerance, precision=3, tolerance_header=1E-5):
        """
        mandatory constructor arguments:
            -tolerance: RMSD of the two files must be smaller than this value
        optional constructor arguments:
            -precision. digits after decimal dot in output
            -tolerance_header: allowed tolerance for values in the file headers
        """
        super().__init__()
        self._tolerance = float(tolerance)
        self._tolerance_header = float(tolerance_header)
        self._rmsd = OrderedDict()
        self._texts = dict()
        precision = int(precision)
        prefix = "{}/{} files equal - "
        location = " in '{}'"
        first = "first diff: "
        maxdiff = "max diff: "
        text = "{:.%iE} (RMSD tol = {:.%iE})"%(precision, precision)+location
        self._texts["rmsd_max"] = prefix + maxdiff + text
        self._texts["rmsd_first"] = prefix + first + text
        text = "'{}' line {} col {}: abs diff = {:.%iE} (tol = {:.%iE})"%(
            precision, precision)
        self._texts["fail_header"] = text

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
        if self._position["tag"] not in self._rmsd:
            data = {"points":0, "maxpoints":1, "pointsum":0.0, "result":None,
                    "volume":0.0, "vectors":[(), (), ()], "headersize":6}
            self._rmsd[self._position["tag"]] = data
        headersize = self._rmsd[self._position["tag"]]["headersize"]
        if self._position["line"] <= headersize:
            return self._eval_header(obj1, obj2)
        elif self._position["line"] >= headersize:
            return self._eval_body(obj1, obj2)

    def _eval_header(self, obj1, obj2):
        """evaluation and comparison of the cubefile header"""
        file = self._position["tag"]
        line = self._position["line"]
        if self._position["line"] <= 2:
            return True, "comment lines"
        else:
            extract = lambda x: [int(x.split()[0])] + [float(y)
                                                       for y in x.split()[1:]]
            data_ref, data_test = self._extract_data(obj1, obj2, extract)
            if line == 3:
                self._rmsd[file]["headersize"] += data_ref[0]
            if line >= 4 and line <= 6:
                self._rmsd[file]["maxpoints"] *= data_ref[0]
                vector = tuple(x/data_ref[0] for x in data_ref[1:])
                self._rmsd[file]["vectors"][line-4] = vector
            if line == 6:
                vecs = self._rmsd[file]["vectors"]
                cross = (vecs[1][1] * vecs[2][2] - vecs[1][2] * vecs[2][1],
                         vecs[1][2] * vecs[2][0] - vecs[1][0] * vecs[2][2],
                         vecs[1][0] * vecs[2][1] - vecs[1][1] * vecs[2][0])
                volume = sum((vecs[0][i] * cross[i] for i in range(3)))
                self._rmsd[file]["volume"] = volume
            for i in range(len(data_ref)):
                diff = abs(data_ref[i] - data_test[i])
                if diff > self._tolerance_header:
                    text = self._texts["fail_header"]
                    return False, text.format(file, line, i+1, diff,
                                              self._tolerance_header)
            return True, "interim analysis"

    def _eval_body(self, obj1, obj2):
        """comparison of the cubefile body"""
        rmsd = self._rmsd[self._position["tag"]]
        extract = lambda x: [float(y) for y in x.split()]
        data_ref, data_test = self._extract_data(obj1, obj2, extract)
        for i in range(max(len(data_ref), len(data_test))):
            rmsd["points"] += 1
            rmsd["pointsum"] += (data_ref[i] - data_test[i])**2
        return True, "interim analysis"

    def _statistics(self):
        """returns the number of successful and total file tests"""
        failed = 0
        for file, data in self._rmsd.items(): #pylint: disable=W0612
            assert data["points"] == data["maxpoints"]
            result = data["pointsum"] * data["volume"]**2 / data["points"]
            data["result"] = result**0.5
            if data["result"] >= self._tolerance:
                failed += 1
        total = len(self._rmsd)
        return total-failed, total

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def largest_error(self):
        passed, total = self._statistics()
        keyfunc = lambda x: x[1]["result"]
        file, data = sorted(self._rmsd.items(), key=keyfunc)[-1]
        text = self._texts["rmsd_max"]
        state = data["result"] < self._tolerance
        return state, text.format(passed, total, data["result"],
                                  self._tolerance, file)

    def first_error(self):
        passed, total = self._statistics()
        for file, data in self._rmsd.items():
            if data["result"] >= self._tolerance:
                text = self._texts["rmsd_first"]
                return False, text.format(passed, total, data["result"],
                                          self._tolerance, file)
        return self.largest_error()

class MullikenFileDiff(AbstractComparison):
    """this class takes care of diffing mulliken charge analysis files"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, tolerance, precision=3):
        super().__init__()
        self._tolerance = float(tolerance)
        self._maxdiffs = dict()
        self._texts = dict()
        precision = int(precision)
        prefix = "{}/{} files equal - "
        location = " in '{}', line {}, col {}"
        first = "first diff: "
        maxdiff = "max diff: "
        text = "{:.%iE} (abs tol = {:.%iE})"%(precision, precision)+location
        self._texts["abs_max"] = prefix + maxdiff + text
        self._texts["abs_first"] = prefix + first + text
        text = "'{}' line {}: section headers do not match"
        self._texts["fail_header"] = text
        text = "'{}' line {}: state index mismatch - {} != {}"
        self._texts["fail_index"] = text

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
        try:
            in_head = not type(int(obj1.split()[0])) == int
        except (ValueError, IndexError):
            in_head = True
        if in_head:
            return self._eval_header(obj1, obj2)
        else:
            return self._eval_body(obj1, obj2)

    def _eval_header(self, obj1, obj2):
        """evaluate the per atom/spin header"""
        file = self._position["tag"]
        line = self._position["line"]
        if line < 2:
            return True, "comment lines"
        elif self._is_comment(obj1) and self._is_comment(obj2):
            return True, "comment lines"
        elif obj1 == obj2:
            return True, "header lines"
        else:
            return False, self._texts["fail_header"].format(file, line)

    def _eval_body(self, obj1, obj2):
        """evaluate the body of a Mulliken output file"""
        if self._position["tag"] not in self._maxdiffs:
            data = {"line":-1, "column":-1, "diff":0.0, "failed":False,
                    "line_first":-1, "column_first":-1, "diff_first":0.0}
            self._maxdiffs[self._position["tag"]] = data
        extract = lambda x: [int(x.split()[0])] + [float(y)
                                                   for y in x.split()[1:]]
        data_ref, data_test = self._extract_data(obj1, obj2, extract)
        diffdata = self._maxdiffs[self._position["tag"]]
        if data_ref[0] != data_test[0]:
            text = self._texts["fail_index"]
            file = self._position["tag"]
            line = self._position["line"]
            return False, text.format(file, line, data_ref[0], data_test[0])
        for i in range(max(len(data_ref[1:]), len(data_test[1:]))):
            diff = abs(data_ref[i+1] - data_test[i+1])
            if diff > diffdata["diff"]:
                diffdata["diff"] = diff
                diffdata["line"] = self._position["line"]
                diffdata["column"] = i + 1
                if not diffdata["failed"] and diff > self._tolerance:
                    diffdata["failed"] = True
                    diffdata["diff_first"] = diff
                    diffdata["line_first"] = self._position["line"]
                    diffdata["column_first"] = i + 1
        return True, "interim analysis"

    def _statistics(self):
        """returns the number of successful and total file tests"""
        failed = len([x for x in self._maxdiffs.items() if x[1]["failed"]])
        total = len(self._maxdiffs)
        return total-failed, total

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def largest_error(self):
        passed, total = self._statistics()
        keyfunc = lambda x: x[1]["diff"]
        file, data = sorted(self._maxdiffs.items(), key=keyfunc)[-1]
        text = self._texts["abs_max"]
        state = False if data["diff"] >= self._tolerance else True
        return state, text.format(passed, total, data["diff"], self._tolerance,
                                  file, data["line"], data["column"])

    def first_error(self):
        passed, total = self._statistics()
        for file, data in self._maxdiffs.items():
            if data["failed"]:
                text = self._texts["abs_first"]
                text = text.format(passed, total, data["diff_first"],
                                   self._tolerance, file, data["line_first"],
                                   data["column_first"])
                return False, text
        return self.largest_error()
