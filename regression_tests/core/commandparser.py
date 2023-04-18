#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import os.path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys

from .remotes.abstractremote import AbstractRemote
from .exceptions import CriticalError
from .utilities import find_program

class CommandParser(object):
    """this is the central parser for command line arguments"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self):
        self._data = dict()
        self._parser = ArgumentParser(description="""Validate the functionality
            of FHI-aims with a regression test. This programm can be performed
            in several modes, you can either run it completely local or you can
            outsource the computational part to a remote system. To obtain a
            detailed call this program with the desired execution mode followed
            by -h""", fromfile_prefix_chars="@", epilog="""You can save your
            preferred set of command line argumenents to a file and then load
            it in future runs by writing @filename as argument.""",
            formatter_class=ArgumentDefaultsHelpFormatter)
        subparsers = self._parser.add_subparsers(dest="mode")
        parser_full = subparsers.add_parser("full", help="""Local Execution and
            Analysis of the Regression Tests""", description="""execute the
            regression tests on your local machine and directly evaluate them
            after all calculations have been performed""",
            formatter_class=ArgumentDefaultsHelpFormatter)
            #conflict_handler='resolve')
        parser_local = subparsers.add_parser("local", help="""Local Execution of
            the Regression Tests, no subsequent analysis""",
            description="""execute the regression tests on your local machine
            without evaluating them afterwards""",
            formatter_class=ArgumentDefaultsHelpFormatter)
        parser_remote = subparsers.add_parser("remote", help="""Remote execution
            of the Regression Tests, without subsequent analysis""",
            description="""prepare the regression tests for execution on a
            remote system, includes archive packing, data transmission,
            unpacking and registration with the queueing system""",
            formatter_class=ArgumentDefaultsHelpFormatter)
        parser_analysis = subparsers.add_parser("analysis", help="""Analysis of
            a previously executed Regression Tests""", description="""analyse
            the results of previously performed regression test. This mode is
            intended to be used as postprocessing for the remote mode""",
            formatter_class=ArgumentDefaultsHelpFormatter)
        self._add_testsuite_group(parser_local)
        self._add_testsuite_group(parser_full)
        self._add_testsuite_group(parser_remote)
        self._add_testsuite_group(parser_analysis)
        self._add_execute_group(parser_local)
        self._add_execute_group(parser_full)
        self._add_execute_group(parser_remote)
        self._add_remote_group(parser_remote)
        self._add_output_group(parser_analysis)
        self._add_output_group(parser_full)

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __getitem__(self, key):
        return self._data[key]

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    @classmethod
    def _add_execute_group(cls, parser):
        """
        add a group with execution specific parameters to the specified parser
        """
        group = parser.add_argument_group("EXECUTION PARAMETERS", """this
            group of settings controls how the regression test is executed""")
        parser.add_argument("ref", action="store", type=str, metavar="REF",
            help="""the reference, this can either be an executable or a folder
            in the local filesystem with the results from a previous regression
            test. If a folder is specified, the data in its "reference"
            subfolder will be used, otherwise the specified executable will be
            used to generate the reference data""")
        parser.add_argument("test", action="store", type=str, metavar="TEST",
            help="""the executable which should be regression tested""")
        group.add_argument("--sourcedir", action="store", type=str,
            metavar="DIR", default="./testcases", help="""path to the input
            files of the regression test""")
        group.add_argument("--mpiexe", action="store", default="mpiexec.hydra",
            type=str, metavar="BIN", help="""the executable of the MPI
            environment, if set to NULL will execute a serial run""")
        group.add_argument("--cpus", action="store", default=4, type=int,
            metavar="N", help="number of CPUs to use")
        group.add_argument("-f", "--force", action="store_true", default=False,
            dest="force", help="""allows the working directory to exist (this
            will purge the folder at the beginning of the regression test!)""")
        group.add_argument("--cml", action="store", default="", type=str,
            metavar="CML", help="""full command line to invoke an FHI-aims
            calculation (this will ignore --mpiexe and --cpus!)""")

    @classmethod
    def _add_output_group(cls, parser):
        """
        add a group with output specific parameters to the specified parser
        """
        group = parser.add_argument_group("OUTPUT PARAMETERS", """this
            group of settings allows you to customize the analysis output""")
        group.add_argument("--nocolors", action="store_true", help="""disable
            color coding in the terminal output""")
        group.add_argument("--logfile", action="store", type=str,
            metavar="FILE", default="./regressionlog.txt", help="""logfile for
            this regressiontest""")
        group.add_argument("--termsize", action="store", type=int,
            metavar="INT", default=None, help="""specify a terminal size for
            output cropping instead of using autodetection""")
        group.add_argument("--batch", action="store_true", default=False,
            help="""if specified, do not stop after each test""")
        group.add_argument("--strict", action="store_true", default=False,
            help="""if specified, failed informative tests count as optional
            fails""")
        group.add_argument("--firstfail", action="store_true", default=False,
            help="""if specified, table and file comparisons will report the
            first difference instead of the largest one""")

    @classmethod
    def _add_remote_group(cls, parser):
        """
        add a group with output specific parameters to the specified parser
        """
        group = parser.add_argument_group("REMOTE CONFIGURATION", """this
            group of settings allows you to configure on which remote server
            you want to execute the regression test""")
        group.add_argument("--ssh", action="store", type=str, required=True,
            metavar="REMOTE", help="""the address of the server (usually for
            ssh)""")
        subs = [x.typename() for x in AbstractRemote.__subclasses__()] #pylint: disable=E1101
        types = [x[0] for x in subs]
        desc = ["%s [%s]"%(x[0], x[1]) for x in subs]
        group.add_argument("--servertype", action="store", default="SGE",
            type=str, metavar="TYPE", help="""the type of server, defined by its
            queuing system, supports: """+", ".join(desc), choices=types)
        group.add_argument("--remoteworkdir", action="store", type=str,
            default="regressiontests", metavar="DIR", help="""the working
            directory on the server (will be created if necessary)""")
        group.add_argument("--remoteconfig", action="store", type=str,
            help="""path to the configuration file for the specified server
            type, its expected content may vary depending on the chosen server
            type""", default="templates/thnec.conf")
        group.add_argument("--archive", action="store", type=str,
            metavar="FILE", default="regression", help="""base name (without
            extention) for the generated archive""")

    @classmethod
    def _add_testsuite_group(cls, parser):
        """
        add a group with the basic parameter like working directory etc to the
        specified parser
        """
        group = parser.add_argument_group("TEST CONFIGURATION PARAMETERS",
            """this group of settings allows you to customize the Regression
            Test to your personal requirements""")
        group.add_argument("--configuration", action="store", type=str,
            metavar="FILE", default="./testsuite.xml", help="""path to the
            test configuration file that should be used""")
        group.add_argument("--defaulttests", action="store", type=str,
            metavar="FILE", default="./defaultparameters.xml", help="""path to
            the file specifying the default tests parameters which apply to all
            test cases""")
        group.add_argument("--workdir", action="store", type=str,
            metavar="DIR", default="./workspace", help="""working directory
            for calculations (must not exist, unless -f is specified""")
        group.add_argument("--exclude", dest="exclusions", type=str,
            default=[], action="append", metavar="TESTSUITE", help="""exclude
            the testsuite with the given name (may be specified more than once)
            """)

    def _validate_analysis(self, args):
        """validate settings for a local regression test"""
        if not os.path.exists(args.workdir) and args.mode == "analysis":
            msg = "Working directory '%s' not found, "%args.workdir
            msg += "requested an analysis of existing data!"
            raise CriticalError(msg)
        self._data["colors"] = not args.nocolors
        if not os.path.exists(os.path.dirname(args.logfile)):
            msg = "Logfile folder '%s' does not exist!"%(args.logfile)
            raise CriticalError(msg)
        self._data["logfile"] = os.path.realpath(args.logfile)
        self._data["interactive"] = not args.batch
        self._data["strict"] = args.strict
        self._data["firstfail"] = args.firstfail
        if args.termsize != None and args.termsize < 20:
            raise CriticalError("No valid terminal size given (must be >= 20)")
        self._data["termwidth"] = args.termsize

    def _validate_execute(self, args):
        """
        validate settings for a local regression test, internally we convert the
        executable pathes to absolute ones because this allows us not to use the
        shell-flag for subprocesses
        """
        if os.path.exists(args.workdir) and not args.force:
            raise CriticalError("Working directory '%s' does "%args.workdir
                + "already exist and --force was not specified!")
        if not os.path.exists(args.ref) and args.mode != "remote":
            raise CriticalError("Reference '%s' does not exist!"%args.ref)
        self._data["ref"] = os.path.realpath(args.ref)
        if not os.path.exists(args.test) and args.mode != "remote":
            executable = find_program(args.test)
            if executable:
                args.test = executable
            else:
                msg = "Test executable '%s' does not exist!" % args.test
                raise CriticalError(msg)
        self._data["test"] = os.path.realpath(args.test)
        if not os.path.exists(args.sourcedir):
            raise CriticalError("Source folder '%s' not found!"%args.sourcedir)
        self._data["basefolder"] = os.path.realpath(args.sourcedir)
        if args.cml:
            self._data["command"] = args.cml.split()
            self._data["test"] = ""
        elif args.mpiexe != "NULL":
            if args.cpus < 1:
                raise CriticalError("%i processes per node does not "%args.cpus
                + "make any sense! Please choose a positive integer.")
            self._data["command"] = [args.mpiexe, "-n", str(args.cpus)]
        else:
            self._data["command"] = []

    def _validate_remote(self, args):
        """validate settings for a local regression test"""
        self._data["ssh"] = args.ssh
        self._data["remoteworkdir"] = args.remoteworkdir
        subs = [(x, x.typename()[0]) for x in AbstractRemote.__subclasses__()] #pylint: disable=E1101
        for remote, shorthand in subs:
            if shorthand == args.servertype:
                self._data["servertype"] = remote
                break
        self._data["archivename"] = args.archive
        if not os.path.exists(args.remoteconfig):
            raise CriticalError("Non-existing remote configuration specified")
        self._data["remoteconfig"] = args.remoteconfig

    def _validate_testsuites(self, args):
        """validate settings that make up the test configuration"""
        self._data["workdir"] = os.path.realpath(args.workdir)
        if not os.path.exists(args.configuration):
            raise CriticalError("Configuration file '%s' does not exist!"%(
                args.configuration))
        self._data["configuration"] = os.path.realpath(args.configuration)
        if not os.path.exists(args.defaulttests):
            msg = "Default test config file '%s' does not exist!"%(
                args.defaulttests)
            raise CriticalError(msg)
        self._data["defaults"] = os.path.realpath(args.defaulttests)
        self._data["exclusions"] = [x.lower().strip() for x in args.exclusions]

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def parse_commandline(self, args):
        """read arguments from the command line"""
        if len(args) == 0:
            print("You have not specified any parameters. Please use the -h "
                + "option to see the built-in help.")
            sys.exit()
        return self._parser.parse_args(args)

    def validate_params(self, options):
        """validate the given command line parameters"""
        self._data["mode"] = options.mode
        self._validate_testsuites(options)
        if options.mode != "analysis":
            self._validate_execute(options)
        if options.mode in ["full", "analysis"]:
            self._validate_analysis(options)
        if options.mode == "remote":
            self._validate_remote(options)
