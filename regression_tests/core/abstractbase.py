#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -AbstractBase: the template for most other classes

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

from functools import total_ordering

from .exceptions import InvalidAttributeException, UnknownAttributeException

@total_ordering #pylint: disable=R0921
class AbstractBase(object):
    """
    abstract base class defining magic functions and common properties
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, attr, parent):
        self._parent = parent
        self._name = self._parse_name(attr.pop("name", ""))

    def _consistency_check(self, attr):
        """
        after initializing all attributes, make sure that the attr dict is
        empty, i.e. no unprocessed keywords remain
        """
        if len(attr) > 0:
            raise UnknownAttributeException(self.__class__, self.name,
                                            self._parent, attr)

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.name.strip().lower() == other.name.strip().lower()

    def __hash__(self):
        return hash((self.__class__, self.name.lower().strip()))

    def __lt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.name < other.name

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    @property
    def mode(self):
        """the current operational mode"""
        return self._parent.mode

    @property
    def name(self):
        """the display name of this object"""
        return getattr(self, "_name", "------")

    @property
    def parent(self):
        """the object which contains this object"""
        return self._parent

    @property
    def reference(self):
        """the reference, either an executable or a folder with results"""
        return self._parent.reference

    @property
    def testsubject(self):
        """the test values, either an executable or a folder with results"""
        return self._parent.testsubject

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _parse_name(self, name):
        """make sure the specified name is not empty or whitespace only"""
        if name.strip() == "":
            data = {"attr":"name", "value":name, "valid":"Any non-empty string"}
            raise InvalidAttributeException(self.__class__, self.name,
                                            self.parent, data)
        else:
            return name.strip()

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def evaluate(self):
        """evaluate the results of the last run for this object"""
        raise NotImplementedError

    def execute(self):
        """execute this object on the local machine, i.e. run the aims-job"""
        raise NotImplementedError

    def prepare(self):
        """
        prepare this object for being executed (fhi-aims run) within a not
        specified runtime environment (can be either local or remote)
        """
        raise NotImplementedError

    def prepare_remote(self, remote):
        """
        additional preparations for the run that should be executed if the
        regression test is executed on a remote system
        """
        raise NotImplementedError

    def output(self, string, level):
        """formatted output, simultaneously to terminal and logfile"""
        self.parent.output(string, level)
