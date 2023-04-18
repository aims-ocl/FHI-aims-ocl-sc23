#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -TestSuite: a grouping object for similar TestCases, allows to easily
                restrict the regression test to specific subsets

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

from .abstractbase import AbstractBase
from .states import StateCase, StateSuite
from .utilities import ColoredString

class TestSuite(AbstractBase):
    """
    groups a set of similar testcases to allow group defaults and disabling
    of test groups

    This constructor may raise:
    KeyError - if no name is given
    InvalidAttributeException - if one of the given attributes cannot be
        transformed into the required data type
    UnknownAttributeException - if the passed attr dict contains unknown
        parameters, which might be typos of known parameters
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, attr, parent):
        super().__init__(attr, parent)
        self._tests = set()
        self._consistency_check(attr)

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    @property
    def command(self):
        """the console command to execute (splitted into a list of strings)"""
        return self._parent.command

    @property
    def folder_base(self):
        """the folder where the regression test input files are located"""
        return self._parent.folder_base

    @property
    def folder_test(self):
        """the working directory where the calculations should be performed"""
        return self._parent.folder_test

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def add_test(self, testcase):
        """
        adds the given TestCase to this TestSuite, duplicates will raise a
        ValueError
        """
        if testcase in self._tests:
            msg = "Test case '%s' was already defined for test suite '%s'"%(
                testcase.name, self.name)
            raise ValueError(msg)
        self._tests.add(testcase)

    def evaluate(self):
        self.output(ColoredString("TestSuite [b]{0}[/]:", (self.name,)), 0)
        self.output(ColoredString("", tuple()), 1)
        status = {StateCase.passed:0, StateCase.optfail:0, StateCase.failed:0}
        for test in sorted(self._tests):
            result = test.evaluate()
            status[result] += 1
            self._parent.interaction()
        finalstatus = StateSuite.elevate_state(status)
        text = "TestSuite [b]{0}[/]: {1} [{2} passed, {3} optfails, {4} fails]"
        self.output(ColoredString(text, (self.name, finalstatus,
            status[StateCase.passed], status[StateCase.optfail],
            status[StateCase.failed])), 0)
        self.output(ColoredString("", tuple()), 0)
        return ColoredString(None, (self.name, finalstatus,
            str(status[StateCase.passed]), str(status[StateCase.optfail]),
            str(status[StateCase.failed])))

    def execute(self):
        for test in sorted(self._tests):
            test.execute()

    def prepare(self):
        for case in sorted(self._tests):
            case.prepare()

    def prepare_remote(self, remote):
        for case in sorted(self._tests):
            case.prepare_remote(remote)
