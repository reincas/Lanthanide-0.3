#!/usr/bin/env python
# -*- Python -*-

import getopt

from Lanthanide import *
from Ion import Ion
from Matrix import Matrix
from Term import Term

energyparms = [
    "H1:2",
    "H1:4",
    "H1:6",
    "e1",
    "e2",
    "e3",
    "H2",
    "H3:0",
    "H3:1",
    "H3:2",
    "H4:2",
    "H4:3",
    "H4:4",
    "H4:6",
    "H4:7",
    "H4:8",
    "Hss:0",
    "Hss:2",
    "Hss:4",
    "Hsoo:0",
    "Hsoo:2",
    "Hsoo:4",
    "H5:0",
    "H5:2",
    "H5:4",
    "H5fix",
    "H6:2",
    "H6:4",
    "H6:6",
    "H6fix",
    ]

matrixparms = [
    "Bz",
    "L2",
    "Lz",
    "S2",
    "Sz",
    "J2",
    "Jz",
    "Vkkk:2:2:2",
    "Vkkk:2:2:4",
    "Vkkk:2:4:4",
    "Vkkk:2:4:6",
    "Vkkk:4:4:4",
    "Vkkk:4:4:6",
    "Vkkk:2:6:6",
    "Vkkk:4:6:6",
    "Vkkk:6:6:6",
    "T1:1:-1",
    "T1:1:+0",
    "T1:1:+1",
    "U1:1:-1",
    "U1:1:+0",
    "U1:1:+1",
    "U1:2:-2",
    "U1:2:-1",
    "U1:2:+0",
    "U1:2:+1",
    "U1:2:+2",
    "U1:4:-4",
    "U1:4:-3",
    "U1:4:-2",
    "U1:4:-1",
    "U1:4:+0",
    "U1:4:+1",
    "U1:4:+2",
    "U1:4:+3",
    "U1:4:+4",
    "U1:6:-6",
    "U1:6:-5",
    "U1:6:-4",
    "U1:6:-3",
    "U1:6:-2",
    "U1:6:-1",
    "U1:6:+0",
    "U1:6:+1",
    "U1:6:+2",
    "U1:6:+3",
    "U1:6:+4",
    "U1:6:+5",
    "U1:6:+6",
    ]

crystalparms = [
    "cf:0:0",
    "cf:2:0",
    "cf:2:2",
    "cf:4:0",
    "cf:4:2",
    "cf:4:4",
    "cf:6:0",
    "cf:6:2",
    "cf:6:4",
    "cf:6:6",
    ]

def calcMatrix(ion, names, key="slj"):
    config = ion.config
    term = ion.term
    for name in names:
        Matrix(config, term, key, name)

def usage():
    name = os.path.basename(sys.argv[0])
    print "Usage: %s [-h] -i ion -t type" % name
    print "       -h:       help"
    print "       -i ion:   rare earth ion to prepare"
    print "       -t type:  term, energy, matrix, or crystal"
    sys.exit()

ion = ""
key = ""
try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:t:")
    for arg, val in opts:
	if arg == "-h":
            usage()
	elif arg == "-i":
	    ion = val
	elif arg == "-t":
	    key = val
	else:
            usage()
except:
    usage()
if not ion or not key:
    usage()
if not key in ("term", "energy", "matrix", "crystal"):
    usage()

re = Ion(ion)
if key == "term":
    Term(re.config)
elif key == "energy":
    calcMatrix(re, energyparms)
elif key == "matrix":
    calcMatrix(re, matrixparms)
elif key == "crystal":
    calcMatrix(re, crystalparms)
