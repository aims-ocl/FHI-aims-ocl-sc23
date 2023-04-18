#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -TestCase: encapsulates a single aims calculation (or a logical group, like
               like a fragment calculation)

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import os
import shutil
import subprocess

from .abstractbase import AbstractBase
from .states import StateTest, StateCase
from .utilities import ColoredString, ColoredTable
from .exceptions import InvalidAttributeException
from .testparameter import TestParameter

class TestCase(AbstractBase):
    """
    this class encapsulates a single regression test with the location of the
    input files and any special variables/files to be compared
    The parameters for this constructors are:
        name - the name of this test case
        folder - the path to folder where the input files are located, relative
        to the folder_base from the RegressionSuite
        subfolders - the name of subfolders in the main folder of this testcase
        skip - any default parameters to skip in this particular test
        defaultfolder - the default output folder for the default parameters
        defaultfile - the default output file for the default parameters
    This constructor may raise:
        InvalidAttributeException - if one of the given attributes cannot be
            transformed into the required data type
        UnknownAttributeException - if the passed attr dict contains unknown
            parameters, which might be typos of known parameters
    """
    _defaultparams = list()

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    @classmethod
    def add_default_param(cls, attr):
        """
        add a attribute dict for instancing a TestParameter to the default list,
        i.e. create a TestParameter from it and append it to all CONSECUTIVELY
        created TestCases instances
        """
        cls._defaultparams.append(attr)

    def __init__(self, attr, parent):
        super().__init__(attr, parent)
        self._folder = self._parse_folder(attr.pop("folder", ""))
        self._subfolders = self._parse_subfolders(attr.pop("subfolders", ""))
        self._defaultfolder = self._parse_subfolders(attr.pop("defaultfolder", ""))[0]

        for param in self._defaultparams:
            param["file"] = os.path.join(self._defaultfolder,
                                         attr.pop("defaultfile", "aims.out"))

        self._parameters = set()
        skips = tuple(x.strip() for x in attr.pop("skip", "").split(','))
        for param in [x for x in self._defaultparams if x["name"] not in skips]:
            self.add_test_parameter(TestParameter(param.copy(), self))
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
    def exe_ref(self):
        """the executable to be used for the reference calculations"""
        return self._parent.exe_ref

    @property
    def exe_test(self):
        """the executable to be used for the test calculations"""
        return self._parent.exe_test

    @property
    def folder_base(self):
        """the folder where the regression test input files are located"""
        path = os.path.join(self._parent.folder_base, self._folder)
        return os.path.abspath(os.path.realpath(path))

    @property
    def folder_test(self):
        """
        the working directory where the calculations should be performed, given
        as a normalized absolute path (with resolved symbolic links)
        """
        return os.path.abspath(os.path.realpath(self._parent.folder_test))

    @property
    def folder_case(self):
        """the name of the folder for this TestCase"""
        return self._folder

    @property
    def subfolders(self):
        """all subfolders in this job that need to computed"""
        return self._subfolders

    @property
    def defaultfolder(self):
        """the default output folder for all checks that need to be performed"""
        return self._defaultfolder

    @property
    def defaultfile(self):
        """the default output file for all checks that need to be performed"""
        return self._defaultfile
    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    @classmethod
    def _run_calc(cls, command, folder):
        """
        run the given command in the specified subfolder, redirecting
        stdout and stderr to the specified streams
        """
        with open(os.path.join(folder, "aims.out"), "w") as out:
            with open(os.path.join(folder, "aims.err"), "w") as err:
                subprocess.call(command, stderr=err, stdout=out, cwd=folder)

    def _parse_folder(self, folder):
        """
        parse the folder name and make sure that it exists within the specified
        base folder (at the RegressionSuite level)
        """
        inputpath = os.path.join(self._parent.folder_base, folder)
        if self._parent.mode == "analysis":
            return folder
        if os.path.isabs(folder) or not os.path.exists(inputpath):
            raise InvalidAttributeException(TestCase, self.name, self._parent,
                {"attr":"folder", "value":folder, "valid":"Any existing "
                + "path, relative to the folder_base from the Regressionsuite"})
        return folder

    def _parse_subfolders(self, subfolders):
        """
        parse subfolders with additional input files, order of appearance
        equals the order of calculation. Note that defining subfolders will skip
        the main directory, unless it is listed too.
        """
        if subfolders != "":
            subs = tuple((x.strip() for x in subfolders.split(",")))
        else:
            subs = (".",)
        if self._parent.mode == "analysis":
            return subs
        exists = lambda x: not os.path.exists(os.path.join(self.folder_base, x))
        if any([os.path.isabs(x) or exists(x) for x in subs]):
            raise InvalidAttributeException(TestCase, self.name, self._parent,
                {"attr":"subfolders", "value":subfolders, "valid":"A comma-"
                + "separated list of subfolders with calculations, relative to "
                + "the main folder of this TestCase."})
        return subs

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def add_test_parameter(self, item):
        """
        generates a TestParameter from the given XML-attributes and adds
        it to this TestCase
        """
        if item in self._parameters:
            self.output(ColoredString("Overriding Parameter \"" + item.name +
                "\" for TestCase \"" + self.name + "\"", (self.name, "dummy")), 0)
            self._parameters.remove(item)
        self._parameters.add(item)

    def evaluate(self):
        head = ColoredString("TestCase [b]{0}[/]: [folder '{1}']", (self.name,
            self.folder_case))
        statustable = ColoredTable(6, head, ["<", ">", ">", "<", "<", "<"])
        head = ColoredString(None, ("Property", "Ref-Value", "Test-Value",
            "Status", "Importance", "Details"))
        statustable.add_row(head)
        statustable.add_separator()
        status = StateTest.equal
        for test in sorted(self._parameters):
            for result, row in test.evaluate():
                statustable.add_row(row)
                status = max(status, result)
        for row, indent in statustable.formatted_output(1, newline=False):
            self.output(row, indent)
        status = StateCase.elevate_state(status)
        self.output(ColoredString("TestCase [b]{0}[/]: {1}", (self.name,
            status)), 1)
        self.output(ColoredString("", tuple()), 1)
        return status

    def execute(self):
        self.output(ColoredString("Now running: {}", (self.name,)), 0)
        if not os.path.isdir(self.reference):
            path = os.path.join(self.folder_test, "reference", self._folder)
            for folder in [os.path.join(path, x) for x in self._subfolders]:
                command = self.command + [self.reference]
                self._run_calc(command, folder)
        path = os.path.join(self.folder_test, "test", self._folder)
        for folder in [os.path.join(path, x) for x in self._subfolders]:
            command = self.command + [self.testsubject]
            self._run_calc(command, folder)

    def prepare(self):
        subpath = os.path.join("reference", self._folder)
        inputpath = os.path.join(self._parent.folder_base, self._folder)
        if os.path.isdir(self.reference):
            shutil.copytree(os.path.join(self.reference, subpath),
                os.path.join(self._parent.folder_test, subpath))
        else:
            shutil.copytree(inputpath,
                os.path.join(self._parent.folder_test, subpath))
        subpath = os.path.join("test", self._folder)
        shutil.copytree(inputpath,
                os.path.join(self._parent.folder_test, subpath))

    def prepare_remote(self, remote):
        remote.process_testcase(self)

