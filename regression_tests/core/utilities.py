#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -indent function
    -ColoredTable class
    -ColoredString class

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import re
from .states import StateCase, StateOptional, StateSuite, StateTest
import os.path
import shutil

class ColoredTable(object):
    """
    this object is a container for accumulating the results from a testcase and
    formatted/aligned output
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, cols, headline, alignments=None):
        assert type(cols) == int and cols >= 1
        assert type(headline) == ColoredString or headline == None
        assert alignments == None or (all(x in "<^>" for x in alignments)
            and len(alignments) == cols)
        self._cols = cols
        self._rows = list()
        self._headline = headline
        self._aligns = alignments if alignments else ["<" for x in range(cols)]

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __len__(self):
        return len(self._rows)

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _calc_columnwidths(self):
        """compute the maximal width each column needs"""
        widths = [0 for x in range(self._cols)]
        for row in (x.get_data_widths() for x in self._rows if x and x.empty):
            for i in range(self._cols):
                widths[i] = max(widths[i], row[i])
        return widths

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def add_row(self, row):
        """add a new row in form of a colored string to this table"""
        assert type(row) == ColoredString
        self._rows.append(row)

    def add_separator(self):
        """ add a line of ---- as a separator to the table"""
        self._rows.append(None)

    def formatted_output(self, indentlevel, newline=True):
        """
        generate a sequence of colorized strings of this table for the
        commandline, optionally color coding can be disabled as a command line
        parameter
        """
        assert type(indentlevel) == int and indentlevel >= 0
        widths = self._calc_columnwidths()
        tableformat = " | ".join(["{:%s%i}" for x in widths])
        fmtdata = zip(self._aligns, widths)
        sepformat = "-|-".join(["{:-%s%i}"%(x, y) for x, y in fmtdata])
        if self._headline:
            if self._headline.empty:
                self._headline.set_template(tableformat, self._aligns, widths)
            yield (self._headline, indentlevel)
        for row in self._rows:
            if row == None:
                row = ColoredString(sepformat, ["" for x in range(self._cols)])
            elif row.empty:
                row.set_template(tableformat, self._aligns, widths)
            yield (row, indentlevel+1)
        if newline:
            yield (ColoredString("", tuple()), 0)

class ColoredString(object):
    """
    utility class that repesents a string which can use colors and other
    highlights when printed to a terminal, but is represented as a plain string
    when written to a textfile.
    """
    _fmtdict = dict()
    _replacers = dict()
    _re_plain = re.compile(r"\x1b\[[^m]+m")

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, string, data):
        """
        Input:
        string - either a format string with %s placeholders or None for a
                 simple table row
        data   - the values that should be inserted into the placeholders
        """
        assert type(string) == str or string == None
        if string != None:
            self._string_color = self._colorize((string,))[0]
            self._string_plain = self._simplify((self._string_color,))[0]
        else:
            self._string_color = None
            self._string_plain = None
        self._data_color = self._colorize(data)
        self._data_plain = self._simplify(self._data_color)

    @classmethod
    def setup_fmt_dict(cls):
        """define all known format replacements"""
        passed = "[ok]PASSED[/]"
        failed = "[alarm]FAILED[/]"
        optfail = "[warn]FAILED[/]"
        optfail_full = "[warn]FAILED (OPTIONAL ONLY)[/]"
        fmtdata = {
            StateTest.equal:passed,
            StateTest.opt_notequal:optfail,
            StateTest.opt_failed:optfail,
            StateTest.opt_missing:optfail,
            StateTest.info_wrong:optfail,
            StateTest.notequal:failed,
            StateTest.failed:failed,
            StateTest.missing:failed,
            StateCase.passed:passed,
            StateCase.failed:failed,
            StateCase.optfail:optfail_full,
            StateSuite.passed:passed,
            StateSuite.failed:failed,
            StateSuite.optfail:optfail,
            StateOptional.mandatory:"mandatory",
            StateOptional.optional:"optional",
            StateOptional.consistent:"consistent",
            StateOptional.informative:"informative"
            }
        cls._fmtdict.update(fmtdata)
        cls._replacers["[alarm]"] = "\x1b[31m"
        cls._replacers["[warn]"] = "\x1b[33m"
        cls._replacers["[ok]"] = "\x1b[32m"
        cls._replacers["[/]"] = "\x1b[0m"
        cls._replacers["[b]"] = "\x1b[1m"
        cls._replacers["[/b]"] = "\x1b[0m"

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __eq__(self, other):
        if type(other) == ColoredString:
            return self.format_terminal() == other.format_terminal()
        else:
            return NotImplemented

    def __repr__(self):
        return """ColoredString: "%s" <-- %s"""%(self._string_plain,
            str(self._data_plain))

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    @property
    def empty(self):
        """true, if this instance contains no valid string for formatting"""
        return self._string_color == None

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _colorize(self, inputdata):
        """returns a colorized version of the given data tuple"""
        data = [str(self._fmtdict.get(x, x)) for x in inputdata]
        for i, item in enumerate(data):
            for key, val in self._replacers.items():
                item = item.replace(key, val)
            data[i] = item
        return data

    def _simplify(self, inputdata):
        """returns a decolorized version of the given colored data tuple"""
        return [self._re_plain.sub("", x) for x in inputdata]

    def _get_padding(self):
        """get the extra lengths needed for the colorized string"""
        data = zip(self._data_color, self._data_plain)
        return [len(x)-len(y) for x, y in data]

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def format_plain(self):
        """plain representation of this string"""
        return self._string_plain.format(*self._data_plain)

    def format_terminal(self):
        """
        colored/formatted version of this string, if template is given it will
        be used instead of the stored string
        """
        return self._string_color.format(*self._data_color)

    def get_padwidth(self):
        """returns the number of invisible characters in this string"""
        return sum(self._get_padding())

    def get_data_widths(self):
        """compute the width in characters of each data element"""
        return [len(x) for x in self._data_plain]

    def set_template(self, template, aligns, widths):
        """overwrite the format string of this instance"""
        colorwidths = [x+y for x, y in zip(widths, self._get_padding())]
        data_color = zip(aligns, colorwidths)
        data_plain = zip(aligns, widths)
        self._string_color = template%tuple(x for y in data_color for x in y)
        self._string_plain = template%tuple(x for y in data_plain for x in y)


def find_program(cmd):
    try:
        program = shutil.which(cmd)
    except AttributeError:
        import subprocess
        try:
            program = subprocess.check_output(['which', cmd]).decode().strip()
        except subprocess.CalledProcessError:
            program = None
    if program:
        return os.path.realpath(program)
