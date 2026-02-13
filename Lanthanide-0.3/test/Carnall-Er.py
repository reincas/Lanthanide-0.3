#
# (c) 1978 by W. T. Carnall et al.
# (c) 01/2001 by Reinhard Caspary
#
# Data from:
#   W. T. Carnall, H. Crosswhite, H. M. Crosswhite
#   "Energy level structure and transition probabilities of the
#    trivalent lanthanides in LaF3"
#   ANL-78-XX-95, Argonne National Laboratory Report, 1978
#   (available as microfiche)
#
# Attention:
#   There is an error in appendix XII, table 1A. Compared to table 1,
#   term 4S3/2 has energy 18583 cm-1 instead of 8583 cm-1.
#

from Lanthanide import *
from Ion import *


# Parameters:
parms = [
    [ 250537.8,  "offset", 1 ],
    [ 100274,    "H1:2",   0 ],
    [  70555,    "H1:4",   0 ],
    [  49900,    "H1:6",   0 ],
    [   2381,    "H2",     0 ],
    [     17.88, "H3:0",   0 ],
    [   -599,    "H3:1",   0 ],
    [   1719,    "H3:2",   0 ],
    [    441,    "H4:2",   1 ],
    [     42,    "H4:3",   1 ],
    [     64,    "H4:4",   1 ],
    [   -314,    "H4:6",   1 ],
    [    387,    "H4:7",   1 ],
    [    363,    "H4:8",   1 ],
    [      4.58, "H5fix",  0 ],
    [    852,    "H6fix",  0 ],
    ]
       
       
# Measurement:


# Calculation:
energies = [
    [   0, "4I15/2",   217,  1 ],
    [   1, "4I13/2",  6712,  1 ],
    [   2, "4I11/2", 10346,  1 ],
    [   3, "4I9/2",  12597,  1 ],
    [   4, "4F9/2",  15455,  1 ],
    [   5, "4S3/2",  18583,  1 ],
    [   6, "2H11/2", 19337,  1 ],
    [   7, "4F7/2",  20715,  1 ],
    [   8, "4F5/2",  22376,  1 ],
    [   9, "4F3/2",  22712,  1 ],
    [  10, "4F9/2",  24756,  1 ],
    [  11, "4G11/2", 26631,  1 ],
    [ -12, "4G9/2",  27637,  1 ],
    [ -13, "2K15/2", 27922,  1 ],
    [ -14, "4G7/2",  28224,  1 ],
    [ -15, "2P1/2",  33319,  1 ],
    ]

Er = Ion("Er3+", parms)
Er.Term(TERM_SLJ)
#Er.Fit(energies, 30000)
Er.showLevels(energies, 40000)
