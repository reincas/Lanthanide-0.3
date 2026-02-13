#
# (c) 1970 by W. T. Carnall et al.
# (c) 02/2001 by Reinhard Caspary
#
# Data from:
#   W. T. Carnall, P. R. Fields, J. Morrison, R. Sarup
#   "Absorptions Spectrum of Tm(3+):LaF3"
#   J Chem Phys 52(8), 1970, pp. 4054-4059
#

from Lanthanide import *
from Ion import *
import TermNames
import TermReduce



parms = [
    [ -3, "H6:2", 1 ],
    [ -1, "H6:4", 1 ],
    [  2, "H6:6", 1 ],
]



#-------------------------------------------------------------------------
# My own calculations:

ion = Ion("Tm3+", parms)
ion.Term()
myHci = ion.hamiltonBuild()
myHci = TriToMatrix(myHci)



#-------------------------------------------------------------------------
# Table from Carnall:

p2 = parms[0][0]/225.0
p4 = parms[1][0]/1089.0
p6 = parms[2][0]*25/184041.0

t11red = [
    0, # 1S
    0, # 1D
    0,
    0, # 1G
    0,
    0,
    0, # 1I
    0,
    0,
    0,
    -315*p2 -693*p4 -1287*p6, # 3P
    math.sqrt(15/2.0)*(137*p2 +198*p4 +143*p6),
    0,
    0,
    -150*p2 -264*p4 +858*p6,
    0, # 3F
    -math.sqrt(10)*((219/2.0)*p2 +165*p4 +(1287/2.0)*p6),
    math.sqrt(11)*(85*p2 +263*p4 +325*p6),
    0,
    0,
    -math.sqrt(14)*(95*p2 +198*p4 +143*p6),
    0, # 3H
    0,
    -math.sqrt(10)*((155/2.0)*p2 +254*p4 +(923/2.0)*p6),
    math.sqrt(13/2.0)*(105*p2 +210*p4 +423*p6),
    0,
    0,
    -math.sqrt(55)*(80*p2 +180*p4 +416*p6),
    ]
Hci = -math.sqrt(3) * Numeric.array(t11red)
Hci = TriToMatrix(Hci)
Hci = TermReduce.expandSL(ion.term, Hci, 1, MATRIX_H6)



#-------------------------------------------------------------------------
# Comparison of results:

diff = Numeric.absolute(Numeric.absolute(myHci)-Numeric.absolute(Hci))
maximum = max(Numeric.ravel(diff))
names = TermNames.short(ion.term, "slj")

#showMatrix(myHci, names)
#showMatrix(Hci, names)
#showMatrix(diff)

print "Maximum difference of Hci (f12) to Carnall: %g" % maximum
