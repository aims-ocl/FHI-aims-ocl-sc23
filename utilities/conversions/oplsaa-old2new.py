#!/usr/bin/python

import sys
import optparse

USAGE = """%prog [options] tinker-xyz

Convert a TINKER OPLS-aa-(~2001) xyz to a OPLS-aa-(2008) xyz file.
"""

parser = optparse.OptionParser(usage=USAGE)
parser.add_option("-r", "--reverse", action="store_true",
                  help="Convert from newer to older format")
parser.add_option("-g", "--generate", action="store_true",
                  help="""\
Generate translation table instead.
In this case, two files are needed as parameters which contain
the same structure with the atom numbering from both the
force fields.
""")

options, args = parser.parse_args()

TRANS = {1: (77, 'CT', 'CT'),
         2: (78, 'CT', 'CT'),
         6: (82, 'HC', 'HC'),
         11: (87, 'CA', 'CA'),
         12: (88, 'HA', 'HA'),
         15: (91, 'CT', 'CT'),
         74: (163, 'CT', 'CT'),
         82: (174, 'C', 'C'),
         83: (175, 'O', 'O'),
         85: (177, 'N', 'N'),
         86: (178, 'N', 'N'),
         88: (180, 'H', 'H'),
         93: (185, 'CT', 'CT'),
         107: (206, 'C', 'C'), # _C_OOH
         108: (208, 'OH', 'OH'), # CO_O_H
         109: (207, 'O', 'O'), # C_O_OH
         110: (209, 'HO', 'HO'), # COO_H_
         111: (210, 'C', 'C'),
         112: (211, 'O2', 'O2'),
         146: (227, 'N3', 'N3'),
         151: (230, 'H3', 'H3'),
         155: (232, 'CT', 'CT'),
         181: (184, 'CT', 'CT'),
         183: (222, 'CT', 'CT')}

if options.reverse:
    newtrans = dict()
    for key, value in TRANS.items():
        newkey, atom2, atom1 = value
        newtrans[newkey] = (key, atom1, atom2)
    trans = newtrans
else:
    trans = TRANS



if options.generate:
    import pprint

    filename1, filename2 = args

    file1 = open(filename1)
    file2 = open(filename2)

    title1 = file1.next()
    title2 = file2.next()

    n1 = int(title1.split()[0])
    n2 = int(title2.split()[0])

    def error(msg):
        sys.stderr.write(msg + "\n")
        sys.exit(2)

    if n1 != n2: error("Numbers of atoms differ")

    d1to2 = dict()
    d2to1 = dict()
    for line1, line2 in zip(file1, file2):
        fields1 = line1.split()
        fields2 = line2.split()
        i1, name1, x1, y1, z1, atom1 = fields1[0:6]
        i2, name2, x2, y2, z2, atom2 = fields2[0:6]
        connects1 = fields1[6:]
        connects2 = fields2[6:]
        if i1 != i2: error("Atom number differ: %s %s(?)" % (i1, i2))
        if connects1 != connects2: error("Connectiveties differ at %s" % i1)
        atom1, atom2 = map(int, (atom1, atom2))
        if atom1 in d1to2:
            if atom2 != d1to2[atom1][0]:
                error("atom1 %s to %s and %s" % (atom1, d1to2[atom1], atom2))
        else:
            d1to2[atom1] = atom2, name1, name2
        if atom2 in d2to1:
            if atom1 != d2to1[atom2][0]:
                error("atom2 %s to %s and %s" % (atom2, d2to1[atom2], atom1))
        else:
            d2to1[atom2] = atom1, name2, name1

    pprint.pprint(d1to2)
    sys.exit(0)


if len(args) == 0:
    in_ = sys.stdin
    out = sys.stdout
if len(args) == 1:
    inname, = args
    in_ = open(inname)
    out = sys.stdout
elif len(args) == 2:
    inname, outname = args
    in_ = open(inname)
    out = open(outname, "w")
else:
    parser.error("Need one or two arguments")
    
for i, line in enumerate(in_):
    if i == 0:
        out.write(line)
    else:
        fields = line.split()
        old = int(fields[5])
        new, _name1, _name2 = trans[old]
        fields[5] = str(new)
        out.write("\t".join(fields) + "\n")

in_.close()
out.close()
