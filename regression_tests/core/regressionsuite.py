#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -RegressionSuite: orchestrates the evaluation of all defined tests

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import shutil

from .abstractbase import AbstractBase
from .utilities import ColoredString, ColoredTable

class RegressionSuite(AbstractBase):
    """
    contains a set of all testsuites to execute (with possible exclusions) and
    runtime relevant parameters (working directories, MPI process count, etc)

    All arguments for the constructor are already checked in 'validate_params'
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, params):
        super().__init__({"name":"Regression Test"}, None)
        self._config = {"mode":params["mode"], "log":"", "color":False,
            "strict":False}
        if params["mode"] in ("full", "analysis"):
            self._config["interactive"] = params["interactive"]
        if params["mode"] == "analysis":
            self._base_folder = ""
            self._runinfo = None
        else:
            self._base_folder = params["basefolder"]
            self._runinfo = (params["test"], params["ref"], params["command"])
        if params["mode"] in ("full", "analysis"):
            self._config["log"] = params["logfile"]
            self._config["color"] = params["colors"]
            self._config["strict"] = params["strict"]
            self._config["firstfail"] = params["firstfail"]
            self._config["termwidth"] = params["termwidth"]
        self._workdir = params["workdir"]
        self._exclusions = params["exclusions"]
        self._suites = set()
        self._logfile = None
        self._config["termwidth"] = self._terminal_size()
        self._failed = False

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __enter__(self):
        self._logfile = open(self._config["log"], "w")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._logfile.close()

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    @property
    def command(self):
        """the console command to execute (splitted into a list of strings)"""
        return self._runinfo[2]

    @property
    def mode(self):
        return self._config["mode"]

    @property
    def testsubject(self):
        """the executable to be used for the test calculations"""
        return self._runinfo[0]

    @property
    def reference(self):
        """the executable to be used for the reference calculations"""
        return self._runinfo[1]

    @property
    def folder_base(self):
        """the folder where the regression test input files are located"""
        return self._base_folder

    @property
    def folder_test(self):
        """the working directory where the calculations should be performed"""
        return self._workdir

    @property
    @classmethod
    def name(cls):
        """
        this is a fake function to allow TestSuites to have a parent without
        failing to throw exceptions (due to missing parent name)
        """
        return None

    @property
    def strict(self):
        """true if informative tests should be treated as optional ones"""
        return self._config["strict"]

    @property
    def firstfail(self):
        """
        true if the first encountered difference, not the largest one is desired
        """
        return self._config["firstfail"]

    @property
    def suites_failed(self):
        """
        true if at least one non-optional tests in the regression suite has
        failed
        """
        return self._failed

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    @classmethod
    def _indent(cls, indentlevel):
        """returns an indented version of the input message"""
        if indentlevel > 0:
            msg = indentlevel * "|   "
            return msg, "\033[90m" + msg + "\033[0m"
        else:
            return "", ""

    def _terminal_size(self):
        """try to determine the size of the terminal"""
        try:
            if self._config["termwidth"] != None:
                return self._config["termwidth"]
            import fcntl, termios, struct
            data = struct.unpack('HHHH',
                fcntl.ioctl(0, termios.TIOCGWINSZ,
                struct.pack('HHHH', 0, 0, 0, 0)))
            return data[1]
        except Exception:
            return 80

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def add_suite(self, suite):
        """
        adds the given TestCase to this TestSuite, duplicates will raise a
        ValueError
        """
        if suite in self._suites:
            msg = "Test suite '%s' was already defined"%(suite.name)
            raise ValueError(msg)
        self._suites.add(suite)

    def evaluate(self):
        head = ColoredString("Final Summary of the RegressionTest", tuple())
        aligns = ["<", "^", "^", "^", "^"]
        statustable = ColoredTable(5, head, alignments=aligns)
        statustable.add_row(ColoredString(None, ("TestSuite", "Status",
            "passed", "optfails", "failed")))
        statustable.add_separator()
        for suite in sorted(self._suites):
            if suite.name.lower().strip() not in self._exclusions:
                statustable.add_row(suite.evaluate())
        self.output(ColoredString("", tuple()), 0)
        for row, indent in statustable.formatted_output(0):
            self.output(row, indent)
            # Check to see if any non-optional tests have failed.  If at least
            # one non-optional test has failed, mark the entire regression suite
            # as failed and stop checking.
            if not self._failed and len(row._data_plain) > 4:
                self._failed = row._data_plain[1] == 'FAILED' and int(row._data_plain[4]) > 0

    def execute(self):
        for suite in sorted(self._suites):
            if suite.name.lower().strip() not in self._exclusions:
                suite.execute()

    def prepare(self):
        try:
            shutil.rmtree(self.folder_test)
        except OSError:
            pass
        for suite in sorted(self._suites):
            if suite.name.lower().strip() not in self._exclusions:
                suite.prepare()

    def prepare_remote(self, remote):
        remote.process_regression(self)
        for suite in sorted(self._suites):
            if suite.name.lower().strip() not in self._exclusions:
                suite.prepare_remote(remote)

    def output(self, string, level):
        """formatted output, simultaneously to logfile and terminal"""
        assert type(string) == ColoredString
        indents = self._indent(level)
        plain = string.format_plain().rstrip()
        if self._config["color"]:
            termindent = indents[1]
            terminal = string.format_terminal().rstrip()
            length = self._config["termwidth"] + string.get_padwidth()
        else:
            termindent = indents[0]
            terminal = plain
            length = self._config["termwidth"]
        length -= len(indents[0])
        if len(terminal) > length:
            terminal = terminal[:length-3] + "..."
        print(termindent + terminal)
        if self._logfile:
            self._logfile.write(indents[0] + plain + "\n")

    def interaction(self):
        """
        if interactive mode was requested, stop after each testcase and wait for
        user input before proceeding
        """
        if self._config["interactive"]:
            input("Press enter to proceed to the next test case...")
