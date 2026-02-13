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
from CarnallPrParms import *

# Parameters:
aparms = [
    [  11393.2,  "offset",   0 ],
    [  69305,    "H1:2",     1 ],
    [  50675,    "H1:4",     1 ],
    [  32813,    "H1:6",     1 ],
    [    750.8,  "H2",       1 ],
    [     21,    "H3:0",     1 ],
    [   -842,    "H3:1",     1 ],
    [   1625,    "H3:2",     1 ],
    [      1.99, "H5fix",    1 ],
    [    200,    "H6fix",    1 ],
    ]

# Measurement:
meas = [
    [  0, "3H4",   200,  1 ],
    [  2, "3H6",  4487,  1 ],
    [  3, "3F2",  5215,  1 ],
    [  4, "3F3",  6568,  1 ],
    [  5, "3F4",  7031,  1 ],
    [  6, "1G4", 10001,  1 ],
    [  7, "1D2", 17047,  1 ],
    [  8, "3P0", 20927,  1 ],
    [  9, "3P1", 21514,  1 ],
    [ 11, "3P2", 22746,  1 ],
    [ 12, "1S0", 46986,  1 ],
    ]
 
# Calculation:
energies = [
    [  0, "3H4",   191,  1 ],
    [  1, "3H5",  2303,  1 ],
    [  2, "3H6",  4495,  1 ],
    [  3, "3F2",  5196,  1 ],
    [  4, "3F3",  6595,  1 ],
    [  5, "3F4",  7009,  1 ],
    [  6, "1G4", 10012,  1 ],
    [  7, "1D2", 17052,  1 ],
    [  8, "3P0", 20935,  1 ],
    [  9, "3P1", 21555,  1 ],
    [ 10, "1I6", 21743,  1 ],
    [ 11, "3P2", 22690,  1 ],
    [ 12, "1S0", 46986,  1 ],
    ]

parms = PrParms
parms[0][0] = 11393.2

Pr = Ion("Pr3+", parms)
Pr.Term()
#Pr.Fit(meas)
Pr.showLevels(energies)
#Pr.saveParms("CarnallPrParms.new")
