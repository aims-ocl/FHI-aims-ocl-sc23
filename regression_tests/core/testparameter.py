#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -TestParameter: a single quantity of interest that should be compared
                    between the test and reference calculations

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import importlib
import os
import re
from glob import glob
from itertools import zip_longest

from .abstractbase import AbstractBase
from .states import StateTest, StateOptional
from .utilities import ColoredString
from .exceptions import InvalidAttributeException, MatchingException

class TestParameter(AbstractBase):
    """
    a comparison test declaration with the following components:
    -name - the name of this test component
    -file - the file to search (defaults to the main output file)
    -regex - the regular expression to find, if none is given, the whole file
        will be checked for equality
    -occurences - if regex is given, this defines the "range" of the comparison:
        given the notation "b:e:s" this means start with the b-th hit and
        compare each s-th occurence thereafter until (inclusive) reaching the
        e-th occurence (defaults to first occurence only)
    -comparison - if given, import the given function from another module and
        use it as a equality operator, use function closures to add additional
        parameters to the function call (defaults to stripped string comparison)
    -optional - if set to True, this parameter is not required to be equal to
        succeed the test, it might even be missing completely

    This constructor may raise:
        KeyError - if no name is given
        InvalidAttributeException - if one of the given attributes cannot be
            transformed into the required data type
        UnknownAttributeException - if the passed attr dict contains unknown
            parameters, which might be typos of known parameters
    """
    regex_importcheck = re.compile(r"(?P<mod>[A-z0-9_]+).(?P<func>[A-z0-9_]+)"
                        + r"[ ]*(?:\((?P<args>.*)\)|)$")

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, attr, parent):
        """
        Note that the different parameters are initialised indepentenly of each
        other, it is therefore sufficient to test the individual parsing
        functions in isolation.
        """
        super().__init__(attr, parent)
        self._file = attr.pop("file", "aims.out")
        self._comparison = self._import_comparison(attr.pop("comparison"))
        self._tablerows = self._parse_tablerows(attr.pop("tablerows", None))
        self._regex = self._parse_regex(attr.pop("regex", None))
        self._optional = self._parse_importance(attr.pop("importance", ""))
        self._occurences = self._parse_occurence(attr.pop("occurence", "1"))
        self._consistency_check(attr)

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    @property
    def optional(self):
        """Is this test optional to succeed the testcase?"""
        return self._optional

    @property
    def parent(self):
        """the parent object of this item"""
        return self._parent

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _analyse_comparison_state(self, result):
        """
        transform the boolean returned by the comparison into a StateTest
        instance, depending on the importance of this test and global strictness
        """
        if not result and self._optional == StateOptional.optional:
            state = StateTest.opt_notequal
        elif (not result and self._optional == StateOptional.informative and
            not self.parent.parent.parent.strict):
            state = StateTest.info_wrong
        elif (not result and self._optional == StateOptional.informative and
            self.parent.parent.parent.strict):
            state = StateTest.opt_notequal
        elif not result and not self._optional == StateOptional.optional:
            state = StateTest.notequal
        else:
            state = StateTest.equal
        return state

    def _diff_files(self, filepath_ref, filepath_test):
        """
        given the two files (or sets of files), diff them line-by-line* with
        the specified function

        * the use of instance variables on your comparison classes allows for
        more complex behaviour if you need it
        """
        reffiles = self._get_filelist(filepath_ref, "reference")
        testfiles = self._get_filelist(filepath_test, "test")
        for file in sorted(set(reffiles.keys()).union(set(testfiles.keys()))):
            errors = dict()
            ref = "<FILE>"
            test = "<FILE>"
            try:
                refvalues = open(reffiles[file], "r")
            except (IOError, KeyError) as exception:
                errors["ref"] = exception
                ref = "------"
            try:
                testvalues = open(testfiles[file], "r")
            except (IOError, KeyError) as exception:
                errors["test"] = exception
                test = "------"
            if len(errors) > 0:
                yield self._error_analysis(errors, self.name, ref, test)
                raise StopIteration
            for i, (ref, test) in enumerate(zip_longest(refvalues, testvalues)):
                if isinstance(self._comparison, object):
                    difffunc = self._comparison.compare_file
                else:
                    difffunc = self._comparison
                result, msg = difffunc(file, i+1, ref, test)
                state = self._analyse_comparison_state(result)
                if state != StateTest.equal:
                    yield state, ColoredString(None, (self.name, "<FILE>",
                        "<FILE>", state, self._optional, msg))
                    raise StopIteration
        if self.parent.parent.parent.firstfail:
            result, msg = self._comparison.first_error()
        else:
            result, msg = self._comparison.largest_error()
        state = self._analyse_comparison_state(result)
        yield state, ColoredString(None, (self.name, "<FILE>", "<FILE>",
            state, self._optional, msg))

    def _diff_regex(self, filepath_ref, filepath_test):
        """
        given the two files, extract data according to given regex and compare
        them with the specified function
        """
        reffiles = self._get_filelist(filepath_ref, "reference")
        testfiles = self._get_filelist(filepath_test, "test")
        fileset = set(reffiles.keys()).union(set(testfiles.keys()))
        for i, file in enumerate(sorted(fileset)):
            desc = self.name + (" #1" if len(self._occurences) > 1 else "")
            if len(fileset) > 1:
                desc += " [file #%i]"%i
            refvalues = self._extract_regex(reffiles[file])
            testvalues = self._extract_regex(testfiles[file])
            file = os.path.relpath(reffiles[file], os.path.join(
                self.parent.folder_test, "reference", self.parent.folder_case))
            if len(fileset) == 1:
                file = None
            errors = dict()
            ref = None
            test = None
            for index in enumerate(self._occurences):
                try:
                    ref = next(refvalues)
                except(IOError, MatchingException) as exception:
                    errors["ref"] = exception
                try:
                    test = next(testvalues)
                except(IOError, MatchingException) as exception:
                    errors["test"] = exception
                if len(errors) > 0:
                    yield self._error_analysis(errors, desc, ref, test)
                    break
                result, msg = self._comparison(file, ref, test)
                state = self._analyse_comparison_state(result)
                yield state, ColoredString(None, (desc, ref, test, state,
                                                  self._optional, msg))
                desc = self.name + " #%i"%(index[0]+2)
                desc += " [file #%i]"%i if len(fileset) > 1 else ""

    def _diff_table(self, filepath_ref, filepath_test):
        """
        given the two files, extract data tables according to given regex and
        compare them with the specified function

        if a tablesize was defined for this testparameter, this function will
        work in the table mode, i.e. the input to the comparison function is
        not a single value, but a list of strings, containing the header of the
        table section and the specified number of following lines
        """
        reffiles = self._get_filelist(filepath_ref, "reference")
        testfiles = self._get_filelist(filepath_test, "test")
        fileset = set(reffiles.keys()).union(set(testfiles.keys()))
        for i, file in enumerate(sorted(fileset)):
            desc = self.name + (" #1" if len(self._occurences) > 1 else "")
            if len(fileset) > 1:
                desc += " [file #%i]"%i
            refvalues = self._extract_table(reffiles[file])
            testvalues = self._extract_table(testfiles[file])
            file = os.path.relpath(reffiles[file], os.path.join(
                self.parent.folder_test, "reference", self.parent.folder_case))
            errors = dict()
            ref = "<TABLE>"
            test = "<TABLE>"
            if isinstance(self._comparison, object):
                difffunc = self._comparison.compare_file
            else:
                difffunc = self._comparison
            for index in enumerate(self._occurences):
                try:
                    ref = next(refvalues)
                except(IOError, MatchingException) as exception:
                    errors["ref"] = exception
                try:
                    test = next(testvalues)
                except(IOError, MatchingException) as exception:
                    errors["test"] = exception
                if len(errors) > 0:
                    yield self._error_analysis(errors, desc, "<TABLE>",
                        "<TABLE>")
                    break
                for i, rows in enumerate(zip(ref, test)):
                    result, msg = difffunc(index, i+1, rows[0], rows[1])
                    state = self._analyse_comparison_state(result)
                    if state != StateTest.equal:
                        yield state, ColoredString(None, (desc, "<TABLE>",
                                                          "<TABLE>", state,
                                                          self._optional, msg))
                        break
                if state != StateTest.equal:
                    continue
                if self.parent.parent.parent.firstfail:
                    result, msg = self._comparison.first_error()
                else:
                    result, msg = self._comparison.largest_error()
                state = self._analyse_comparison_state(result)
                desc = self.name
                desc += " #%i"%(index[0]+1) if len(self._occurences) > 1 else ""
                desc += " [file #%i]"%i if len(fileset) > 1 else ""
                yield state, ColoredString(None, (desc, "<TABLE>", "<TABLE>",
                                                  state, self._optional, msg))

    def _error_analysis(self, errors, desc, ref, test):
        """expects a set of (origin, error) pairs to analyze what went wrong"""
        msg = ""
        fails = {"ref":False, "test":False}
        data = {"ref":ref, "test":test}
        for comp in ("ref", "test"):
            if comp in errors:
                data[comp] = "------"
                fails[comp] = type(errors[comp])
                if msg != "":
                    msg += ", "
                if type(errors[comp]) == MatchingException:
                    msg += "%s: %s"%(comp, repr(errors[comp]))
                elif type(errors[comp]) == IOError:
                    msg += "%s: file '%s' missing"%(comp,
                        os.path.basename(errors[comp].filename))
                elif type(errors[comp]) == KeyError:
                    msg += "%s: file '%s' missing"%(comp,
                        os.path.basename(errors[comp].args[0]))
        samefail = fails["ref"] == fails["test"]
        if self._optional == StateOptional.consistent and samefail:
            state = StateTest.equal
        elif self._optional == StateOptional.optional:
            state = StateTest.opt_missing
        elif (self.optional == StateOptional.informative and not
            self.parent.parent.parent.strict):
            state = StateTest.info_wrong
        elif (self.optional == StateOptional.informative and
            self.parent.parent.parent.strict):
            state = StateTest.opt_missing
        else:
            state = StateTest.missing
        return state, ColoredString(None, (desc, data["ref"], data["test"],
                                           state, self._optional, msg))

    def _extract_regex(self, filepath):
        """
        generator which will extract the desired occurences of the regex from
        the given file, raises an Exception if the file is missing or yielded
        too few matches
        """
        hits = 0
        lastmatch = self._occurences[-1]
        if lastmatch > 0:
            with open(filepath, 'r', encoding='utf-8') as infile:
                for line in infile:
                    match = self._regex.search(line)
                    if match:
                        hits += 1
                    if match and hits in self._occurences:
                        yield match.group(1)
                    if hits == lastmatch:
                        break
            if hits != lastmatch:
                raise MatchingException(hits, lastmatch, self.name, filepath)
        else:
            matches = list()
            with open(filepath, 'r', encoding='utf-8') as infile:
                for line in infile:
                    match = self._regex.search(line)
                    if match:
                        matches.append(match.group(1))
            try:
                for occurence in self._occurences:
                    yield matches[occurence]
            except IndexError:
                raise MatchingException(len(matches), self._occurences[0]
                    - self._occurences[-1], self.name, filepath)

    def _extract_table(self, filepath):
        """
        generator which will extract the desired table occurences of the regex
        from the given file, raises an Exception if the file is missing or
        yielded too few matches
        """
        hits = 0
        lastmatch = self._occurences[-1]
        if lastmatch > 0:
            with open(filepath, 'r') as infile:
                for line in infile:
                    match = self._regex.search(line)
                    if match:
                        hits += 1
                    if match and hits in self._occurences:
                        rowcount = range(self._tablerows)
                        table = [infile.readline() for x in rowcount]
                        yield [line,] + table
                    if hits == lastmatch:
                        break
            if hits != lastmatch:
                raise MatchingException(hits, lastmatch, self.name, filepath)
        else:
            matches = list()
            with open(filepath, 'r') as infile:
                for line in infile:
                    match = self._regex.search(line)
                    if match:
                        rowcount = range(self._tablerows)
                        table = [infile.readline() for x in rowcount]
                        matches.append([line,] + table)
            try:
                for occurence in self._occurences:
                    yield matches[occurence]
            except IndexError:
                raise MatchingException(len(matches), self._occurences[0]
                    - self._occurences[-1], self.name, filepath)

    def _get_filelist(self, filepattern, middlepart):
        """given a glob-compatible expression, build the required path-dict"""
        basefolder = os.path.join(self.parent.folder_test, middlepart,
            self.parent.folder_case)
        relpath = lambda x: os.path.relpath(x, basefolder)
        files = {relpath(x):x for x in glob(filepattern)}
        if len(files) == 0:
            files = {relpath(filepattern):filepattern}
        return files

    def _import_comparison(self, importname):
        """
        import a function from a given module and construct the closure from
        the additional arguments, if any
        """
        try:
            match = self.regex_importcheck.match(importname)
            modname = "comparisons.%s"%(match.group("mod"))
            module = importlib.import_module(modname)
            if match.group("args") != None:
                paramlist = match.group("args").split(",")
                args = tuple(x for x in paramlist if "=" not in x)
                keywords = tuple(x.split("=") for x in paramlist if "=" in x)
                keywords = {x[0].strip():x[1] for x in keywords}
                return getattr(module, match.group("func"))(*args, **keywords) #pylint: disable=W0142
            else:
                return getattr(module, match.group("func"))
        except (AttributeError, ImportError):
            data = {"attr":"comparison", "value":importname,
                "valid":"any function (closure) from the comparison submodules"}
            raise InvalidAttributeException(self.__class__, self.name,
                self._parent, data)

    def _parse_importance(self, importance):
        """transform the specied importance string into an StateOptional"""
        try:
            return StateOptional[importance.lower()]
        except KeyError:
            data = {"attr":"importance", "value":importance,
                    "valid":"'mandatory', 'optional' or 'consistent'"}
            raise InvalidAttributeException(self.__class__, self.name,
                self._parent, data)

    def _parse_occurence(self, occurences):
        """
        parse the given occurence range and make sure it is a sane, nonempty set
        """
        try:
            occ = [int(x) for x in occurences.split(":")]
            if len(occ) == 1:
                occ.append(occ[0])
                occ.append(1)
            elif len(occ) == 2:
                occ.append(1)
            if occ[0] == 0:
                occ[0] = 1
            elif occ[1] == 0:
                occ[1] = -1
            occrange = range(occ[0], occ[1] + 1, occ[2])
            if len(occ) != 3 or len(occrange) == 0:
                raise ValueError
            if any(x < 0 for x in occrange) and any(x > 0 for x in occrange):
                raise ValueError
            return occrange
        except ValueError:
            data = {"attr":"occurences", "value":occurences,
                "valid":"a non-empty range given in the format 'start:stop:"
                + "step' with either positive (first matches) or negative "
                + "(last matches) elements"}
            raise InvalidAttributeException(self.__class__, self.name,
                self._parent, data)

    def _parse_regex(self, regex):
        """compile the specified regex, if any"""
        if regex != None:
            try:
                pattern = re.compile(regex)
                if (pattern.groups != 1 and self._tablerows == None) or (
                    pattern.groups != 0 and self._tablerows != None):
                    raise re.error
                return pattern
            except re.error:
                data = {"attr":"regex", "value":regex,
                "valid":"any valid regular expression where the desired value "
                + "must be in the only matching group, which are denoted by "
                + "embracing (), any normal () must be escaped by a \\, regex "
                + "patterns for table headers must not have any groups"}
                raise InvalidAttributeException(self.__class__, self.name,
                self._parent, data)
        else:
            return None

    def _parse_tablerows(self, size):
        """verify that the size of the table to extract is a positive integer"""
        if size != None:
            try:
                size = int(size)
                assert size >= 1
                return size
            except ValueError:
                data = {"attr":"tablesize", "value":size,
                "valid":"any positive integer > 0"}
                raise InvalidAttributeException(self.__class__, self.name,
                self._parent, data)
        else:
            return None

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def evaluate(self):
        filepath_ref = os.path.join(self._parent.folder_test, "reference",
            self._parent.folder_case, self._file)
        filepath_test = os.path.join(self._parent.folder_test, "test",
            self._parent.folder_case, self._file)
        if self._regex and self._tablerows:
            return self._diff_table(filepath_ref, filepath_test)
        elif self._regex:
            return self._diff_regex(filepath_ref, filepath_test)
        else:
            return self._diff_files(filepath_ref, filepath_test)

    def execute(self):
        pass

    def prepare(self):
        pass

    def prepare_remote(self, remote):
        pass
