#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -AbstractComparison: the template for custom comparison functions

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

class AbstractComparison(object): #pylint: disable=R0921
    """the abstract base for any comparison function"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self):
        self._position = {"tag":"", "line":0, "taglist":set()}
        self._comments = list()

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _compare(self, obj_ref, obj_test):
        """
        compare the two given objects in regex mode, inputs:
            obj_ref - the extracted value from the reference file, as string
            obj_test - the extracted value from the test file, as string

        This function must return a two-tuple (Result, Message) with:
            -Result - Boolean value indicating equality (true) or difference
            -Message - informative string about the test
        """
        raise NotImplementedError

    def _eval_errors(self, errors):
        """generalized error formatting for failed value extraction"""
        if "file" in self._position:
            file = self._position["file"]
        else:
            file = "Unknown file"

        if "line" in self._position:
            line = self._position["line"]
        else:
            line = "Unknown line"

        err = ""
        if "obj1" in errors:
            err = "ref '%s': line %i invalid"%(file, line)
        if "obj2" in errors and len(err) > 0:
            err += ", "
        if "obj2" in errors:
            err += "test '%s': line %i invalid"%(file, line)
        return err

    def _extract_data(self, line_ref, line_test, func_ref, func_test=None):
        """
        given the extraction functions, this helper function will extract the
        desired data from the input lines and takes care of error handling
        if the extraction fails
        """
        func_test = func_ref if func_test == None else func_test
        errors = dict()
        try:
            data_ref = func_ref(line_ref)
        except (Exception) as exception: #pylint: disable=W0703
            errors["ref"] = exception
        try:
            data_test = func_test(line_test)
        except (Exception) as exception: #pylint: disable=W0703
            errors["test"] = exception
        if len(errors) > 0:
            raise ValueError(self._eval_errors(errors))
        return data_ref, data_test

    def _is_comment(self, line):
        """
        checks if a given line is a comment, i.e. either empty or the first
        non-space character is contained in comments
        """
        return len(line.strip()) == 0 or line.strip()[0] in self._comments

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def compare(self, tag, obj_ref, obj_test):
        """
        compare the two given objects in regex mode, inputs:
            file - the (relative) path of the current file
            obj_ref - the extracted value from the reference file, as string
            obj_test - the extracted value from the test file, as string
        """
        self._position["tag"] = tag
        self._position["taglist"].add(tag)
        try:
            return self._compare(obj_ref, obj_test)
        except ValueError as err:
            return False, str(err)


    def compare_file(self, tag, index, line_ref, line_test):
        """
        compare the two given objects in file or table mode, inputs:
            file - the (relative) path of the current file
            index - line number of the currently compared line
            line_ref - the extracted line from the reference file, as string
            line_test - the extracted line from the test file, as string
        """
        self._position["tag"] = tag
        self._position["taglist"].add(tag)
        self._position["line"] = index
        try:
            return self._compare(line_ref, line_test)
        except ValueError as err:
            return False, str(err)

    def first_error(self):
        """
        after analyzing the complete table or file, this function returns the
        info message containing details about the first encountered difference
        """
        raise NotImplementedError

    def largest_error(self):
        """
        after analyzing the complete table or file, this function returns the
        info message containing details about the most severe difference found
        """
        raise NotImplementedError
