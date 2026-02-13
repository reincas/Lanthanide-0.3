#
# (c) 1968 by B. R. Judd et al.
# (c) 02/2001 by Reinhard Caspary
#
# Data from:
#   B. R. Judd, H. M. Crosswhite, H. Crosswhite
#   "Intra-Atomic Magnetic Interactions for f Electrons"
#   Phys Rev 169(1), 1968, pp. 130-138
#

from Lanthanide import *
from Ion import *
import TermNames
import TermReduce



parms = [
    [  2, "H6:0", 1 ],
    [ -3, "H6:2", 1 ],
    [ -1, "H6:4", 1 ],
    [  1, "H6:6", 1 ],
]



#-------------------------------------------------------------------------
# My own calculations:

ion = Ion("Pr3+", parms)
ion.Term()
myHci = ion.hamiltonBuild()
myHci = TriToMatrix(myHci)



#-------------------------------------------------------------------------
# Table from Judd:

p0 = parms[0][0]
p2 = parms[1][0]/225.0
p4 = parms[2][0]/1089.0
p6 = parms[3][0]*25/184041.0

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
    -2*p0 -105*p2 -231*p4 -429*p6, # 3P
    math.sqrt(15/2.0)*(p0 +32*p2 -33*p4 -286*p6),
    0,
    0,
    -p0 -45*p2 -33*p4 +1287*p6,
    0, # 3F
    math.sqrt(10)*(-p0 -(9/2.0)*p2 +66*p4 -(429/2.0)*p6),
    math.sqrt(11)*(p0 -20*p2 +32*p4 -104*p6),
    0,
    0,
    math.sqrt(14)*(-p0 +10*p2 +33*p4 +286*p6),
    0, # 3H
    0,
    math.sqrt(10)*(-p0 +(55/2.0)*p2 -23*p4 -(65/2.0)*p6),
    math.sqrt(13/2.0)*(p0 -21*p4 -6*p6),
    0,
    0,
    math.sqrt(55)*(-p0 +25*p2 +51*p4 +13*p6),
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

print "Maximum difference of Hci (f02) to Judd: %g" % maximum
