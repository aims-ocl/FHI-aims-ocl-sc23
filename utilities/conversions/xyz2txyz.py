#!/usr/bin/python

import sys
import optparse

USAGE = """%prog template [xyzfile [outfile]]

Convert a standard (xmakemol) xyz file to TINKER xyzfile.
The atom type and bonding information is taken from template.
Therefore, template must be a Tinker xyz file with the same
atom sequence.

Can be used in conjunction with aimsout2xyz.py:
$ aimsout2xyz.py aims.out | xyz2txyz.py template.xyz >aims_out.xyz"""

parser = optparse.OptionParser(usage=USAGE)

options, args = parser.parse_args()

if len(args) < 1 or len(args) > 3:
    parser.error("Invalid number of arguments")

templ = open(args[0])
if len(args) >= 2:
    in_ = open(args[1])
else:
    in_ = sys.stdin
if len(args) >= 3:
    out = open(args[2], "w")
else:
    out = sys.stdout

# check numbers of atoms
first_templ_line = templ.next()
n_atom = int(first_templ_line.split()[0])
first_xyz_line = in_.next()
n_xyz_atom = int(first_xyz_line)
if n_atom != n_xyz_atom:
    sys.stderr.write("Numbers of atoms differ between template and xyz.")
    sys.exit(2)
comment = in_.next()
out.write("%i " % n_atom + comment)

for tline, xline in zip(templ, in_):
    tfields = tline.split()
    i, tname, _x, _y, _z, atomtype = tfields[:6]
    xfields = xline.split()
    xname, x, y, z = xfields[:4]
    if not tname.startswith(xname):
        sys.stderr.write("Atom class mismatch in line %s\n" % i)
    out.write("%5s %4s %16s %16s %16s  %3s  %s\n" %
              (i, tname, x, y, z, atomtype, " ".join(tfields[6:])))


templ.close()
in_.close()
out.close()
