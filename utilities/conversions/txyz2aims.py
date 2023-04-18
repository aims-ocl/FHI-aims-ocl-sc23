#!/usr/bin/python

import sys
import optparse
import os
import subprocess

USAGE = """%prog [options] [tinker.xyz [outfile]]

Convert a Tinker xyz file into another format.
By default, the last geometry info is output in a format given
explicitly with the -o flag or deduced from outfile.  The default
is the FHI-aims geometry.in input format.  For formats other than
aims and xyz, babel is called and must be able to understand aims.
"""

def error(msg):
    sys.stderr.write(msg + "\n")
    sys.exit(2)

def read_tinker(in_):
    """Read Tinker-XYZ from iterator in_ and return xyz-list and comment."""
    line = in_.next()
    fields = line.split()
    n_atoms = int(fields[0])
    comment = " ".join(fields[1:])
    xyz = []
    for i_atom, line in enumerate(in_):
        fields = line.split()
        j_atom, name, x, y, z = fields[0:5]
        if int(j_atom) != i_atom+1:
            error("Wrong atom number in line %i" % (i_atom+1))
        if len(name) > 1 and name[1].islower():
            aimsname = name[0:2]
        else:
            aimsname = name[0:1]
        xyz.append((float(x), float(y), float(z), aimsname))
    if len(xyz) != n_atoms:
        error("Wrong number of atoms")
    return xyz, comment


def write_xyz(out, xyz, comment):
    out.write("%i\n" % len(xyz))
    out.write(comment + "\n")
    for x, y, z, name in xyz:
        out.write("%2s %15s %15s %15s\n" % (name, x, y, z))

def write_aims(out, xyz, comments):
    if isinstance(comments, (str, unicode)):
        out.write("# " + comments + "\n")
    else:
        out.writelines(["# " + line + "\n" for line in comments])
    for x, y, z, name in xyz:
        out.write("atom %15s %15s %15s %2s\n" % (x, y, z, name))


if __name__ == "__main__":
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-o", "--outformat", metavar="FORMAT",
                      action="store", default=None,
                      help="Output format extension (default is in for aims)")
    parser.add_option("-c", "--compat-babel", action="store_true",
                      help="Feed babel with xyz format")

    options, args = parser.parse_args()

    if len(args) >= 1:
        inname = args[0]
        in_ = open(inname)
    else:
        inname = "-"
        in_ =sys.stdin

    if len(args) >= 2:
        outname = args[1]
        out = open(outname, "w")
        _basename, _dot, outformat = outname.rpartition('.')
    else:
        outname = "-"
        out = sys.stdout
        outformat = "aims"

    if len(args) >= 3:
        parser.error("Too many arguments")

    if options.outformat is not None:
        outformat = options.outformat
    if outformat == "in":
        outformat = "aims"

    xyz, comment = read_tinker(in_)

    if inname != "-":
        comment += " from file '%s'" % inname

    if outformat == "xyz":
        write_xyz(out, xyz, comment)
    elif outformat == "aims":
        write_aims(out, xyz, comment)
    else:
        if options.compat_babel:
            cmd = ["babel", "-ixyz", "-", "-o"+outformat, outname]
        else:
            cmd = ["babel", "-iaims", "-", "-o"+outformat, outname]
        sys.stderr.write("$ " + " ".join(cmd) + "\n")
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        if options.compat_babel:
            write_xyz(p.stdin, xyz, comment)
        else:
            write_aims(p.stdin, xyz, comment)
        p.communicate()
