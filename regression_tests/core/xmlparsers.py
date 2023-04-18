#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -XMLParserSuites: parsing of main input file, creates test hierarchy
    -XMLParserDefaults: parses the file with the default parameters and makes
     registers them for all tests

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import sys
from xml.parsers import expat

from . import RegressionSuite, TestSuite, TestCase, TestParameter

class XMLParserSuites(object): #pylint: disable=R0903
    """a class to parse the input file with the test definitions"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, params):
        self._regression_suite = RegressionSuite(params)
        self._current = None
        self._parser = expat.ParserCreate()
        self._parser.StartElementHandler = self._start_element
        self._parser.EndElementHandler = self._end_element

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _start_element(self, name, attr):
        """
        handles XML-element starts in the XML-file and constructs the
        RegressionSuite from them
        """
        if name == "regressionsuite" and self._current == None:
            self._current = self._regression_suite
        elif name == "testsuite":
            item = TestSuite(attr, self._current)
            self._current.add_suite(item) #pylint: disable=E1103
            self._current = item
        elif name == "testcase":
            item = TestCase(attr, self._current)
            self._current.add_test(item)
            self._current = item
        elif name == "testparameter":
            item = TestParameter(attr, self._current)
            self._current.add_test_parameter(item)
            self._current = item
        else:
            print("unknown or invalid element: ", name, " attributes: ", attr)
            sys.exit(1)

    def _end_element(self, name):
        """end of XML element handler"""
        if name != "regressionsuite":
            self._current = self._current.parent

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def xmlparse(self, inputfile):
        """parse the given XML-file and return the generated RegressionSuite"""
        self._parser.ParseFile(open(inputfile, mode="rb"))
        return self._regression_suite

class XMLParserDefaults(object): #pylint: disable=R0903
    """
    a class to parse the input file with the default TestParameter which apply
    to all testcases
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self):
        self._current = None
        self._parser = expat.ParserCreate()
        self._parser.StartElementHandler = self._start_element
        self._parser.EndElementHandler = self._end_element

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _start_element(self, name, attr):
        """
        handles XML-element starts in the XML-file and constructs the
        RegressionSuite from them
        """
        if name == "defaults" and self._current == None:
            self._current = True
        elif name == "testparameter" and self._current == True:
            TestCase.add_default_param(attr)
            self._current = attr
        else:
            print("unknown or invalid element: ", name, " attributes: ", attr)
            sys.exit(1)

    def _end_element(self, name):
        """end of XML element handler"""
        if name == "testparameter":
            self._current = True

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def xmlparse(self, inputfile):
        """
        parse the given XML-file and store the defaults in the TestCase class
        """
        self._parser.ParseFile(open(inputfile, mode="rb"))
