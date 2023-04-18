#!/usr/bin/env python3
"""
This module provides a regression test infrastructure for the FHI-aims all-
electron ab initio electronic structure code. The actual test suite is loaded
from an external XML-file. See the testsuite.xml distributed with this script
as an example.
"""
import sys
#python interpreter version check
if sys.version_info[0] != 3 or sys.version_info[1] < 2:
    print("This script requires at least Python version 3.2.\n"
        + "If this python version is not available, "
        + "use the old regression script instead.")
    sys.exit(1)

from subprocess import CalledProcessError

from core.xmlparsers import XMLParserSuites, XMLParserDefaults
from core.commandparser import CommandParser
from core.utilities import ColoredString
from core.exceptions import CriticalError

def main():
    """
    main function, parses command-line parameters and executes regression
    tests accordingly
    """
    try:
        ColoredString.setup_fmt_dict()
        cmdparser = CommandParser()
        arguments = cmdparser.parse_commandline(sys.argv[1:])
        cmdparser.validate_params(arguments)
        xmlparser = XMLParserDefaults()
        xmlparser.xmlparse(cmdparser["defaults"])
        xmlparser = XMLParserSuites(cmdparser)
        suite = xmlparser.xmlparse(cmdparser["configuration"])
        if cmdparser["mode"] != "analysis":
            suite.prepare()
        if cmdparser["mode"] in ("local", "full"):
            suite.execute()
        elif cmdparser["mode"] == "remote":
            remote = cmdparser["servertype"](cmdparser, cmdparser["workdir"])
            with remote:
                suite.prepare_remote(remote)
                remote.finalize_preparation()
            remote.submit()
        if cmdparser["mode"] in ("analysis", "full"):
            with suite:
                suite.evaluate()
            # Output whether the regression tests passed or failed and exit with
            # the proper exit status
            if suite.suites_failed:
                print()
                print("\x1b[31mONE OR MORE TESTS FAILED, CHECK YOUR COMPILATION SETTINGS!\x1b[0m")
                print()
                sys.exit(2)
            else:
                print()
                print("\x1b[32mALL TESTS PASSED, YOU'RE GOOD TO GO!\x1b[0m")
                print()
                sys.exit(0)
    except CalledProcessError as err:
        print("\x1b[31mERROR:\x1b[0m %s"%(err.output))
        sys.exit(1)
    except CriticalError as err:
        print(str(err))
        sys.exit(1)

if __name__ == "__main__":
    main()
