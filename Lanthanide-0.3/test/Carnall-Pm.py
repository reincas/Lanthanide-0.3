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

from Lanthanide import *
from Ion import *


# Parameters:
parms = [
    [ 52914.8,   "offset", 1 ],
    [  77000,    "H1:2",   0 ],
    [  55000,    "H1:4",   0 ],
    [  37500,    "H1:6",   0 ],
    [   1022,    "H2",     0 ],
    [     21.00, "H3:0",   0 ],
    [   -560,    "H3:1",   0 ],
    [   1400,    "H3:2",   0 ],
    [    330,    "H4:2",   1 ],
    [     41.5,  "H4:3",   1 ],
    [     62,    "H4:4",   1 ],
    [   -295,    "H4:6",   1 ],
    [    360,    "H4:7",   1 ],
    [    310,    "H4:8",   1 ],
    [      2.49, "H5fix",  0 ],
    [    440,    "H6fix",  0 ],
    ]
       
       
# Measurement:


# Calculation:
energies = [
    [  0, "5I4",   120,  1 ],
    [  1, "5I5",  1612,  1 ],
    [  2, "5I6",  3239,  1 ],
    [  3, "5I7",  4951,  1 ],
    [  4, "5I8",  6714,  1 ],
    [  5, "5F1", 12638,  1 ],
    [  6, "5F2", 13080,  1 ],
    [  7, "5F3", 13933,  1 ],
    [  8, "5S2", 14486,  1 ],
    [  9, "5F4", 14887,  1 ],
    [ 10, "5F5", 16223,  1 ],
    [ 11, "3K6", 16939,  1 ],
    [ 12, "5G2", 18053,  1 ],
    [ 13, "3H4", 18075,  1 ],
    [ 14, "3K7", 18255,  1 ],
    [ 15, "5G3", 18565,  1 ],
    [ 16, "3K8", 19862,  1 ],
    [ 17, "3H5", 20307,  1 ],
    [ 18, "5G4", 20554,  1 ],
    [ 19, "3G3", 21935,  1 ],
    [ 20, "5G5", 22475,  1 ],
    [ 21, "5G6", 22807,  1 ],
    [ 22, "3D2", 23140,  1 ],
    [ 23, "3L7", 23772,  1 ],
    [ 24, "3P1", 24216,  1 ],
    [ 25, "3H6", 24702,  1 ],
    [ 26, "3G4", 24840,  1 ],
    [ 27, "3L8", 24907,  1 ],
    [ 28, "3P0", 25811,  1 ],
    [ 29, "3D3", 25895,  1 ],
    [ 30, "3L9", 25907,  1 ],
    ]

Pm = Ion("Pm3+", parms)
Pm.Term(TERM_SLJ)
#Pm.Fit(energies, 30000)
Pm.showLevels(energies, 30000)
