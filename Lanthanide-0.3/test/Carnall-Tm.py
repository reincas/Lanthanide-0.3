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
    [ 282780.7,  "offset",   1 ],
    [ 102459,    "H1:2",     1 ],
    [  72424,    "H1:4",     1 ],
    [  51380,    "H1:6",     1 ],
    [   2640,    "H2",       1 ],
    [     17,    "H3:0",     1 ],
    [   -737,    "H3:1",     1 ],
    [   1700,    "H3:2",     1 ],
    [      4.93, "H5fix",    1 ],
    [    729.6,  "H6fix",    1 ],
    ]

# Measurement:
meas = [
    [  0, "3H6",   200,  1 ],
    [  1, "3F4",  5858,  1 ],
    [  2, "3H5",  8336,  1 ],
    [  3, "3H4", 12711,  1 ],
    [  4, "3F3", 14559,  1 ],
    [  5, "3F2", 15173,  1 ],
    [  6, "1G4", 21352,  1 ],
    [  7, "1D2", 28061,  1 ],
    [  8, "1I6", 34886,  1 ],
    [  9, "3P0", 35604,  1 ],
    [ 10, "3P1", 36559,  1 ],
    [ 11, "3P2", 38344,  1 ],
    ]

# Calculation:
energies = [
    [  0, "3H6",   175,  1 ],
    [  1, "3F4",  5818,  1 ],
    [  2, "3H5",  8391,  1 ],
    [  3, "3H4", 12721,  1 ],
    [  4, "3F3", 14597,  1 ],
    [  5, "3F2", 15181,  1 ],
    [  6, "1G4", 21314,  1 ],
    [  7, "1D2", 28001,  1 ],
    [  8, "1I6", 34975,  1 ],
    [  9, "3P0", 35579,  1 ],
    [ 10, "3P1", 36615,  1 ],
    [ 11, "3P2", 38268,  1 ],
    [ 12, "1S0", 75300,  1 ],
    ]

parms = TmParms
parms[0][0] = 282780.7

Tm = Ion("Tm3+", parms)
Tm.Term(TERM_SLJ)
#Tm.Fit(energies)
Tm.showLevels(energies)
