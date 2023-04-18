#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

from functools import total_ordering

class AbstractEnum(object): #pylint: disable=R0903
    """abstract enum implementation, used to avoid python3.4 dependency"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, sort=False):
        self._entries = dict()
        self._sorted = None
        for item, value in self._vars().items():
            self._entries[item] = EnumEntry(self, item, value)
            setattr(self, item, self._entries[item])
        if sort:
            self._sorted = sorted(self._entries.values(), key=lambda x: x.value)

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __call__(self, searchvalue):
        for key, entry in self._entries.items():
            if entry.value == searchvalue:
                return key
        raise ValueError("value %s does not exist in Enumeration %s"%(
            searchvalue, self.__class__.__name__))

    def __getitem__(self, key):
        return self._entries[key]

    def __iter__(self):
        if self._sorted:
            return iter(self._sorted)
        else:
            return iter(self._entries.values())

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    @classmethod
    def _vars(cls):
        """
        returns a filtered dict of the class attributes, listing only public
        attributes in a name->value map
        """
        filterfunc = lambda x: x[0] != "_" and type(getattr(cls, x)) == int
        return {x:getattr(cls, x) for x in dir(cls) if filterfunc(x)}

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def elevate_state(self, other):
        """
        given a state from the lower level, translate it into its equivalent on
        the higher level
        """
        raise NotImplementedError

@total_ordering
class EnumEntry(object):
    """a single entry in an enumeration object"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, parent, name, value):
        self._parent = parent
        self._type = parent.__class__
        self._name = name
        self._value = value

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __add__(self, other):
        return self._parent(self.value + other)

    def __eq__(self, other):
        if self.enumtype == other.enumtype:
            return self.value == other.value
        return NotImplemented

    def __hash__(self):
        return hash((self.enumtype, self.name, self.value))

    def __lt__(self, other):
        if self.enumtype == other.enumtype:
            return self.value < other.value
        return NotImplemented

    def __repr__(self):
        return "<%s.%s: %s>"%(self._type.__name__, self.name, self.value)

    def __sub__(self, other):
        return self._parent(self.value - other)

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    @property
    def enumtype(self):
        """the enumeration this entry belongs too"""
        return self._parent

    @property
    def name(self):
        """the name of this entry"""
        return self._name

    @property
    def value(self):
        """the value of this entry"""
        return self._value

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

class StateTest(AbstractEnum): #pylint: disable=R0903
    """
    comparison states for the value extraction/comparison.
    'failed' states indicate that less values than expected could be extracted
        from the source file
    'missing' states indicate that the source file was completely missing
    'wrong' indicates that the exact issue does not matter
    """
    equal = 1
    opt_notequal = 2
    opt_failed = 3
    opt_missing = 4
    info_wrong = 5
    notequal = 12
    failed = 13
    missing = 14

    def elevate_state(self, other):
        return NotImplemented

class StateCase(AbstractEnum): #pylint: disable=R0903
    """comparison states for the value extraction/comparison"""
    passed = 1
    optfail = 2
    failed = 3

    def elevate_state(self, other):
        assert type(other) == EnumEntry and other.enumtype == StateTest
        if other.value == 1:
            return self.passed
        elif other.value >= 2 and other.value <= 4:
            return self.optfail
        elif other.value == 5:
            return self.passed
        elif other.value >= 12 and other.value <= 14:
            return self.failed
        else:
            return NotImplemented

class StateSuite(AbstractEnum): #pylint: disable=R0903
    """states for TestSuites"""
    passed = 1
    optfail = 2
    failed = 3

    def elevate_state(self, other):
        assert type(other) == dict
        if other[StateCase.failed] > 0:
            return self.failed
        elif other[StateCase.optfail] > 0:
            return self.optfail
        elif other[StateCase.passed] > 0:
            return self.passed
        else:
            return NotImplemented

class StateOptional(AbstractEnum): #pylint: disable=R0903
    """
    comparison states for the optional state of a TestParameter:
    mandatory  - this test must finish successfully to pass the TestCase
    optional   - this test is completely optional and does not affect the final
                 result of the TestCase
    consistent - specifies that the test is optional, if the value extraction
                 for both files ends with the same result (fail/success).
    """
    mandatory = 3001
    consistent = 3002
    optional = 3003
    informative = 3004

    @classmethod
    def elevate_state(cls, other):
        return NotImplemented

#This hack is required to make the imported objects behave like instances
#pylint: disable=C0103
StateOptional = StateOptional(sort=True)
StateSuite = StateSuite(sort=True)
StateCase = StateCase(sort=True)
StateTest = StateTest(sort=True)
#pylint: enable=C0103
