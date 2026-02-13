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
    [  31000.5,  "offset", 1 ],
    [  73036,    "H1:2",   1 ],
    [  52624,    "H1:4",   1 ],
    [  35793,    "H1:6",   1 ],
    [    884.9,  "H2",     1 ],
    [     21.28, "H3:0",   1 ],
    [   -583,    "H3:1",   1 ],
    [   1443,    "H3:2",   1 ],
    [    306,    "H4:2",   1 ],
    [     41,    "H4:3",   1 ],
    [     59,    "H4:4",   1 ],
    [   -283,    "H4:6",   1 ],
    [    326,    "H4:7",   1 ],
    [    298,    "H4:8",   1 ],
    [      2.237,"H5fix",  1 ],
    [    213,    "H6fix",  1 ],
    ]
       
# Measurement:


# Calculation:
energies = [
    [  0, "4I9/2",     235,  1 ],
    [  1, "4I11/2",   2114,  1 ],
    [  2, "4I13/2",   4098,  1 ],
    [  3, "4I15/2",   6148,  1 ],
    [  4, "4F3/2",   11621,  1 ],
    [  5, "4F5/2",   12660,  1 ],
    [  6, "2H9/2",   12768,  1 ],
    [  7, "4F7/2",   13619,  1 ],
    [  8, "4S3/2",   13691,  1 ],
    [  9, "4F9/2",   14899,  1 ],
    [ 10, "2H11/2",  16105,  1 ],
    [ 11, "4G5/2",   17428,  1 ],
    [ 12, "4G7/2",   17469,  1 ],
    [ 13, "4G7/2",   19293,  1 ],
    [ 14, "4G9/2",   19709,  1 ],
    [ 15, "2K13/2",  19785,  1 ],
    [ 17, "2D3/2",   21425,  1 ],
    [ 18, "4G11/2",  21714,  1 ],
    [ 19, "2K15/2",  21780,  1 ],
    [ 20, "2P1/2",   23458,  1 ],
    [ 21, "2D5/2",   24004,  1 ],
    [ 22, "2P3/2",   26424,  1 ],
    ]

Nd = Ion("Nd3+", parms)
Nd.Term(TERM_SLJ)
#Nd.Fit(energies, 30000)
Nd.showLevels(energies, 30000)
