#!/usr/bin/env python

import sys
import re
import optparse
import fileinput
import StringIO

USAGE = """%prog [options] <f90 file names>

Write Fortran 90 module dependencies in given source files to
standard output in GNU make readable format.\
"""

def uc(s):
    return s.upper()
def lc(s):
    return s.lower()

MODULE_RE = re.compile(r'^\s*module\s+(\w+)', re.I)
USE_RE = re.compile(r'^\s*use\s+(\w+)', re.I)
PROGRAM_RE = re.compile(r'^\s*program', re.I)
HASH_RE = re.compile(r'^\#')

# Configurable:
TOCASE = lc
MODSUFFIX = ".mod"

def stem(filename, suffixes=frozenset(['f', 'f90', 'F', 'F90'])):
    """Return stem of filename ('foo.f90' -> 'foo')."""
    #stem, dot, suffix = filename.partition('.')
    i = filename.index('.')
    if i < 0 or filename[i+1:] not in suffixes:
        raise ValueError("Invalid filename %s" % filename)
    stem = filename[:i]
    return stem

def objstr(filename):
    """Return object file of filename ('foo.f90' -> '$(OBJDIR)/foo.o')."""
    return "$(OBJDIR)/%s.o" % stem(filename)

def modstr(mod):
    """Return modulse file from module name ('foo' -> '$(MODDIR)/foo.mod')."""
    return "$(MODDIR)/%s%s" % (mod, MODSUFFIX)

def get_dependson_tree(defmods_d, usemods_d):
    """Generate dependency tree of source files.

    The return value is a dictionary which maps source file names to a set
    of other source file names on which it depends.

    >>> defmods_d = {'a.f90': 'a', 'b.f90': 'b'}
    >>> usemods_d = {'b.f90': 'a'}
    >>> (get_dependson_tree(defmods_d, usemods_d) ==
    ...  {'b.f90': set(['a.f90'])})
    True
    """
    mod2filename = dict()
    for filename, defmods in defmods_d.items():
        for mod in defmods:
            if mod in mod2filename:
                raise ValueError("Module %s defined twice" % mod)
            mod2filename[mod] = filename

    dependson = dict()
    for filename, usemods in usemods_d.items():
        dependson[filename] = set()
        for mod in usemods:
            if mod in mod2filename:
                dependson[filename].add(mod2filename[mod])
    return dependson
    
def recursive_cycle_check(tree, node, trace=None):
    """Check for cycles and remove leaves from tree, start at node.

    >>> dependson = {'b': set(['a']), 'c': set(['b']), 'a': set(['c'])}
    >>> recursive_cycle_check(dependson, 'a') == ['a', 'c', 'b', 'a']
    True
    """

    # Leaf: Fast forward
    if node not in tree:
        return None

    # Update and check trace
    if trace is None:
        trace = []
    mytrace = trace + [node]
    if node in trace:
        return mytrace

    # Check subnodes:
    subnodes = tree[node].copy()   # might change
    for subnode in subnodes:
        # Recurse
        ret = recursive_cycle_check(tree, subnode, mytrace)
        if ret:
            return ret
        # Remove reference to leaf (done after(!) recursion)
        if subnode not in tree:
            tree[node].remove(subnode)
    # Leaves are signified by missing or empty bentries.  The first is better.
    if not tree[node]:   # Now empty
        tree.pop(node)
    return None


def check_dependson_for_cycles(dependson):
    """Check for cycles.

    >>> dependson = {'b.f90': set(['a.f90'])}
    >>> check_dependson_for_cycles(dependson) == None
    True
    >>> dependson['c.f90'] = set(['b.f90'])
    >>> dependson['a.f90'] = set(['c.f90'])
    >>> (set(check_dependson_for_cycles(dependson)) ==
    ...  set(['a.f90', 'b.f90', 'c.f90']))
    True
    """
    # Deep copy.
    mydeps = dict([(key, val.copy()) for key, val in dependson.items()])
    # Get files to iterate first, as mydeps might change.
    filenames = dependson.keys()
    for filename in filenames:
        ret = recursive_cycle_check(mydeps, filename)
        if ret:
            # Return trace from first to second occurence of tail.
            last = ret[-1]
            return ret[ret.index(last):]
    return None

def parse_file(fileiter):
    """Return tuple of sets defmods, usemods, and logicals is_exe, needs_fpp.

    >>> TOCASE = lc
    >>> f = StringIO.StringIO('use EXAMPLEMOD')
    >>> parse_file(f) == (set(), set(['examplemod']), False, False)
    True
    >>> f = StringIO.StringIO('''program\\nmodule examplemod''')
    >>> parse_file(f) == (set(['examplemod']), set([]), True, False)
    True
    """
    defmods = set()
    usemods = set()
    is_exe = False
    needs_fpp = False
    for line in fileiter:
        # The following is performance critical and optimized for speed.
        line = line.lstrip()
        if not line:
            continue
        if line[0] == "#":
            use_fpp = True
            continue
        word = line[0:line.find(' ')].lower()
        if word == "module":
            m = MODULE_RE.match(line)
            if m:
                defmods.add(TOCASE(m.group(1)))
        elif word == "use":
            m = USE_RE.match(line)
            if m:
                usemods.add(TOCASE(m.group(1)))
        elif word == "program":
            if PROGRAM_RE.match(line):
                is_exe = True
    # ignore "module procedure" (from interfaces)
    defmods.discard('procedure')
    # Discard usage of a module defined in this file.
    usemods.difference_update(defmods)
    return defmods, usemods, is_exe, needs_fpp

def parse_files(filenames):
    """Return dicts defmods_d, usemods_d, exes, fpp_needers"""
    defmods_d = dict()
    usemods_d = dict()
    exes = []
    fpp_needers = []
    for name in filenames:
        f = open(name, "r")
        defmods, usemods, is_exe, needs_fpp = parse_file(f)
        defmods_d[name] = defmods
        usemods_d[name] = usemods
        if is_exe:
            exes.append(name)
        if needs_fpp:
            fpp_needers.append(name)
    return defmods_d, usemods_d, exes, fpp_needers


if __name__ == "__main__":

    if "-t" in sys.argv or "--test" in sys.argv:
        import doctest
        doctest.testmod()
        sys.exit(0)

    # --- Parse command line

    formatter = optparse.IndentedHelpFormatter(max_help_position=32)
    parser = optparse.OptionParser(usage=USAGE, formatter=formatter)
    parser.add_option("-o", "--output", default="Makefile.moddeps",
                      help="File name for dependencies")
    parser.add_option("-u", "--upper",
                      dest='tocase', action='store_const', const=uc,
                      help="Use uppercase names")
    parser.add_option("-l", "--lower",
                      dest='tocase',  action='store_const', const=lc,
                      default=lc,
                      help="Use lowercase names (default)")
    parser.add_option("-m", "--suffix", metavar="MOD", default=".mod",
                      help="Suffix (default is .mod)")
    parser.add_option("-x", "--exclude", metavar="REGEX", action='append',
                      help="Do not scan files matching REGEX.")
    parser.add_option("-t", "--test", action="store_true", help="Run doctests")

    (options, args) = parser.parse_args()

    TOCASE = options.tocase
    MODSUFFIX = options.suffix

    if options.output == "-":
        out = sys.stdout
    else:
        out = open(options.output, "w")

    if options.exclude:
        for excl in options.exclude:
            args = [arg for arg in args if not re.search(excl, arg)]

    out.write("# Command line:\n")
    out.write("# " + " ".join(sys.argv) + "\n")
    out.write("\n")


    # --- Parse file

    defmods_d, usemods_d, exes, fpp_needers = parse_files(args)

    # --- Preprocessor?

    if fpp_needers:
        out.write("# List all f90 files which need a preprocessor.\n")
        out.write("# This could be used with a compile line like\n")
        out.write("#    $(FC) -c $($(@:.obj=)_FPPFLAGS) $(FFLAGS) $<\n")
        out.write("\n")
        for filename in fpp_needers:
            out.write("%s_FPPFLAGS = $(FPPFLAGS)\n" % stem(filename))
        out.write("\n\n")

    # --- List of modules

    mods = set([mod
                for filename, defmods in defmods_d.items()
                for mod in defmods
                ])

    if mods: 
        out.write("# List of all .mod files\n")
        modstrs = [modstr(mod) for mod in mods]
        out.write("MODULE_HEADERS = %s\n" % ' '.join(modstrs))
        out.write("\n")

        # --- <file>.mod: <file>.o ;

        out.write("# .mod : defining .o ; # empty rule\n")
        for filename, defmods in defmods_d.items():
            if not defmods: continue
            modfiles = [modstr(mod) for mod in defmods]
            out.write(" ".join(modfiles) + ": " + objstr(filename) + " ;\n")
        out.write("\n")

    # --- <file>.o: <used_mods>.mod

    out.write("# .o : used .mod\n")
    for filename, usemods in usemods_d.items():
        if not usemods: continue
        modfiles = [modstr(mod) for mod in usemods]
        out.write("%s: %s\n" % (objstr(filename), " ".join(modfiles)))

    # --- Check dependency tree

    dependson = get_dependson_tree(defmods_d, usemods_d)
    cycle = check_dependson_for_cycles(dependson)
    if cycle:
        cycstr = ' -> '.join(cycle)
        sys.stderr.write("*** Found dependency cycle: %s\n" % cycstr)
        sys.stderr.write(
            "*** Please remove some use statements to break this.\n")
        # This is an ugly kludge to stop make from proceeding.
        # First, remove Makefile.moddeps to allow retry.
        if options.output != "-":
            out.write("$(shell rm %s)\n" % options.output)
        # Then, bail out exactly this time.
        out.write("$(error Found dependency cycle: %s)\n" % cycstr)
