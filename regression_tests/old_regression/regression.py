#!/usr/bin/python
from __future__ import print_function, with_statement

import sys
import os
import re
import subprocess
import optparse
import shutil
import glob
from xml.parsers import expat

class XMLWrapper(object):

    def __init__(self, targetdict):
        self.current = None
        self.target = targetdict

    def xml_start_element(self, name, attr):
        """
        fill TESTDICT with data from new regression test, discarding all unused
        new information
        """
        if name == "testcase":
            #case distinction for the one test that had numdiff arguments
            diff = None if attr["folder"] != "ScN.plus_u" else "-a 1e-5 -r 1e-3"
            self.target[attr["name"]] = (attr["folder"], list(), diff)
            self.current = attr["name"]
        elif name == "testparameter":
            if "regex" not in attr and "file" in attr:
                self.target[self.current][1].append(attr["file"])

    def xml_end_element(self, name):
        if name == "testcase":
            self.current = None

    def parse_configuration(self, configfile):
        parser = expat.ParserCreate()
        parser.StartElementHandler = self.xml_start_element
        parser.EndElementHandler = self.xml_end_element
        with open(configfile, "r") as infile:
            parser.ParseFile(infile)

TESTDICT = dict()
XMLWrapper(TESTDICT).parse_configuration("../testsuite.xml")

USAGE = """%prog (run|diff|all) BINARY [options]

Test the FHI-aims executable BINARY.  The flags "-n NP" or "-s" can be used to
explicitely specify parallel or serial mode.  There are three modes to run
this script; either only do the test runs (run), only diff results of previous
runs (diff) or do both subsequently (all)."""

# === Helpers

def simplify_binary(binary):
    import os.path
    abspath = os.path.abspath(binary)
    # Does not work because '~' is not recognized.
    # home = os.environ.get('HOME')
    # if home:
    #     abspath = abspath.replace(home, '~')
    return abspath


def outstem_from_binary(binary):
    m = re.match(r'aims\.([^/]*)\.x$', os.path.basename(binary))
    if m:
        return m.group(1)
    else:
        import datetime
        return datetime.datetime.now().strftime('noname-%m%d%y')

def eval_np(mpi_np, binary):
    """Deduce MPI setting from explicit option and binary name."""
    binname = os.path.basename(binary)
    if mpi_np == 'auto':
        if re.search(r'\.mpi\.', binname):
            # print "Use implicit mpi (np=2, disable by '--mpi-np none')"
            np = 8
        elif re.search(r'\.serial\.', binname):
            # print "Use implicit serial (enable mpi by '--mpi-np <n>')"
            np = None
        else:
            parser.error('Cannot infer --mpi-np setting from name %s' % binname)
    elif mpi_np == 'none':
        np = None
    else:
        try:
            np = int(mpi_np)
            if np <=0:
                raise ValueError()
        except ValueError:
            parser.error("Arg to --mpi-np is either none, auto or <n>")
    return np


# === Parse argument list

formatter = optparse.IndentedHelpFormatter(max_help_position=33)
parser = optparse.OptionParser(usage=USAGE, formatter=formatter)

# Test options
parser.add_option("-l", "--list", action='store_true',
                  help="only list tests and exit")
parser.add_option("-x", "--exclude", dest="exclude",
                  help="tests to exclude (comma separated)", metavar="TESTS")
parser.add_option("-t", "--tests", dest="tests",
                  help="only these tests (comma separated)", metavar="TESTS")

parser.add_option("-n", "--mpi-np", dest="mpi_np", default="auto", metavar="NP",
                  help="""use mpiexec ['auto' (deduce from BINARY)]""")
parser.add_option("-s", "--serial", dest="mpi_np", action='store_const',
                  const='none', help="""alias for '--mpi-np none'""")

parser.add_option("--outstem", metavar="OUT", default=None,
                  help="stem for output file names [from BINARY]")
parser.add_option("-r", "--refstem", metavar="REF", default='reference',
                  help="stem for reference file names ['reference']")

parser.add_option("--differ", metavar="PROGRAM", default='../ediff',
                  help="program for diffs")
parser.add_option("-O", "--numdiff-opts", dest="nd_opts", metavar="OPTS",
                  default="", help="options for numdiff.py")


(options, args) = parser.parse_args()


# === Which tests to be done

if options.tests:
    tests = [test.strip() for test in options.tests.split(",")]
    for test in tests:
        if test not in TESTDICT:
            parser.error("Do not know about test %s" % test)
else:
    tests = [test for test in TESTDICT.keys()]

if len(set(tests)) != len(tests):
    parser.error("Duplicate tests in list")

if options.exclude:
    exclude_tests = [test.strip() for test in options.exclude.split(",")]
    for test in exclude_tests:
        try:
            tests.remove(test)
        except ValueError:
            parser.error("Test %s to be excluded not in list" % test)

if options.list:
    sys.stdout.write("Test\tDirectory\n")
    sys.stdout.write("------\t---------------\n")
    for test in tests:
        dir_, aux, ndadd = TESTDICT[test]
        sys.stdout.write("%s\t%s\n" % (test, dir_))
    sys.exit(0)


# === Which binaries to use

try:
    task = args.pop(0)
    binary = args.pop(0)
except IndexError:
    parser.error("Need both (run|diff|both) and BINARY.")

if task not in ["all", "run", "diff"]:
    parser.error("Task should be either run, diff, or all.")

if args:
    parser.error("Need only (run|diff|both) and BINARY.")

# # add ../../bin to PATH
# os_path = os.environ['PATH'].split(':')
# os_path.insert(0, os.path.abspath('../../bin'))
# os.environ['PATH'] = ':'.join(os_path)

if options.outstem:
    outstem = options.outstem
else:
    outstem = outstem_from_binary(binary)

refstem = options.refstem

# === Test subs

def do_test(test, dir_, binary, mpi_np, outfile):
    """Run given binary in given dir with np procs."""
#    print 78*"="
    np = eval_np(mpi_np, binary)
    abs_binary = simplify_binary(binary)
    if np is None:
        cmd = [abs_binary]
    else:
        cmd = ["mpirun", "-n", str(np), abs_binary]
#    print "| %8s: %-47s |" % (test, " ".join(cmd))
    out = open(os.path.join(dir_, outfile), "w")
    null = open('/dev/null', "r")
    subprocess.call(cmd, stdin=null, stdout=out, cwd=dir_)
    out.close()
    null.close()

def look_at_test(dir_, oufile, reffile):
    """Quick comparison of outfile with reffile."""
    subprocess.call(["../compare.pl", outfile, reffile], cwd=dir_)

def diff_test(dir_, differ, ndopts, outfile, reffile, ndadd=None):
    """Detailed comparison (side-by-side diff) of outfile with reffile.

    As numdiff.py also takes directories, this function is also usable
    with directories of auxiliary output.
    """
    cwd = os.getcwd()
    os.chdir(dir_)
    if ndadd is None:
        ndadd = ""
    try:
        cmd = "../numdiff.py -c ../numdiff.rc --differ \"%s\" %s %s %s %s" % (
            differ, ndopts, outfile, reffile, ndadd)
        print(cmd)
        subprocess.call(cmd, shell=True)
    finally:
        os.chdir(cwd)

# === Do work

for test in tests:
    dir_, aux, ndadd = TESTDICT[test]

    # Get binary and outfile
    outfile = "%s.%s.out" % (dir_, outstem)
    reffile = "%s.%s.out" % (dir_, refstem)
    outdir = "%s.%s.data" % (dir_, outstem)
    refdir = "%s.%s.data" % (dir_, refstem)

    # Run test if needed
    if task in ["run", "all"]:
        if not os.path.exists(dir_):
            os.mkdir(dir_)
        for path in glob.glob("../testcases/%s/*"%(dir_)):
            shutil.copy(path, dir_)
        do_test(test, dir_, binary, options.mpi_np, outfile)

        if len(aux) > 0:
            subprocess.call("rm -rf %s/%s" % (dir_, outdir), shell=True)
            subprocess.call("mkdir %s/%s" % (dir_, outdir), shell=True)
            for file_ in aux:
                subprocess.call("mv %s/%s %s/%s/" %
                                (dir_, file_, dir_, outdir), shell=True)

    #copy reference data, if no explicit binary was specified
    if options.refstem == "reference":
        if not os.path.exists(os.path.join(dir_, refdir)) and len(aux) > 0:
            os.mkdir(os.path.join(dir_, refdir))
        refpath = "../references_lastrelease/reference/%s/*"%(dir_)
        for path in glob.glob(refpath):
            if os.path.basename(path) == "aims.out":
                copypath = os.path.join(dir_, reffile)
                subprocess.call("cp %s %s" % (path, copypath), shell=True)
            else:
                copypath = os.path.join(dir_, refdir)
                if not os.path.exists(copypath):
                   os.mkdir(copypath)
                subprocess.call("cp -r %s %s" % (path, copypath + "/."), shell=True)

    # Unconditionally glance at test result
    look_at_test(dir_, outfile, reffile)

    # Diff test if needed
    if task in ["diff", "all"]:
        diff_test(dir_, options.differ, options.nd_opts, outfile, reffile,
                  ndadd)
        # Diff auxs if present.
        if len(aux) > 0:
            diff_test(dir_, options.differ, options.nd_opts, outdir, refdir)
