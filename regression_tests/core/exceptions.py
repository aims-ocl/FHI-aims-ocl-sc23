#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

class MatchingException(Exception):
    """
    an exception which is raised when the scanned file ends before the desired
    number of matches was found
    """

    def __init__(self, found, expected, name, file):
        super().__init__()
        self._found = found
        self._expected = expected
        self._name = name
        self._file = file

    def __eq__(self, other):
        if type(other) == MatchingException:
            return str(self) == str(other)
        else:
            return NotImplemented

    def __hash__(self):
        return hash((MatchingException, self._found, self._expected, self._name,
            self._file))

    def __str__(self):
        return ("Found %i/%i matches in file '%s' for parameter '%s'"%
            (self._found, self._expected, self._file, self._name))

    def __repr__(self):
        return "%i/%i matches found"%(self._found, self._expected)

class UnknownAttributeException(Exception):
    """
    an exception which is raised when a constructor recieves a attribute dict
    with unknown entries (to prevent silent ignorance of typos)
    """

    def __init__(self, cls, name, parent, attr):
        super().__init__()
        self._cls = cls.__name__
        self._attr = attr
        self._name = name
        self._parent = parent
        self.__cause__ = None

    def __str__(self):
        if self._parent:
            msg = "Constructor for %s '%s' [parent: %s '%s'] "%(self._cls,
                self._name, self._parent.__class__.__name__, self._parent.name)
        else:
            msg = "Constructor for %s '%s' "%(self._cls, self._name)
        msg += "recieved unknown arguments: %s"%(self._attr)
        return msg

class InvalidAttributeException(Exception):
    """
    an exception which is raised when a constructor recieves a attribute dict
    which contains an invalid value
    """

    def __init__(self, cls, name, parent, data):# pylint: disable=W0231
        self._cls = cls.__name__
        self._name = name
        self._parent = parent
        self._attr = data["attr"]
        self._value = data["value"]
        self._valid = data["valid"]
        self.__cause__ = None

    def __str__(self):
        if self._parent and self._parent.name != None:
            msg = "Constructor for %s '%s' [parent: %s '%s'] "%(self._cls,
                self._name, self._parent.__class__.__name__, self._parent.name)
        else:
            msg = "Constructor for %s '%s' "%(self._cls, self._name)
        msg += "recieved an invalid argument '%s' for attribute '%s'"%(
            self._value, self._attr)
        if self._valid:
            msg += " (valid arguments for this attribute are: %s)"%(self._valid)
        return msg

class CriticalError(Exception):
    """
    exception-subclass for critical, but expected errors, these will fail
    gracefully without stack trace because their origin is known
    """

    def __init__(self, msg):
        super().__init__()
        self._message = msg

    def __str__(self):
        return "\x1b[31mERROR:\x1b[0m %s"%(self._message)
