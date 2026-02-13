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
    [ -1, "Hss:0", 1 ],
    [  3, "Hss:2", 1 ],
    [  1, "Hss:4", 1 ],
]



#-------------------------------------------------------------------------
# My own calculations:

ion = Ion("Tm3+", parms)
ion.Term()
myHss = ion.hamiltonBuild()
myHss = TriToMatrix(myHss)



#-------------------------------------------------------------------------
# Table from Judd:

m0 = parms[0][0]
m2 = parms[1][0]
m4 = parms[2][0]

T22red = [
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
    0, # 3P
    0,
    0,
    0,
    -12*m0 -24*m2 -(300/11.0)*m4,
    0, # 3F
    0,
    0,
    0,
    (8/math.sqrt(3))*(3*m0 +m2 -(100/11.0)*m4),
    (4/3.0)*math.sqrt(14)*(-m0 +8*m2 -(200/11.0)*m4),
    0, # 3H
    0,
    0,
    0,
    0,
    (8/3.0)*math.sqrt(11/2.0)*(2*m0 -(23/11.0)*m2 -(325/121.0)*m4),
    (4/3.0)*math.sqrt(143)*(m0 -(34/11.0)*m2 -(1325/1573.0)*m4),
    ]
Hss = math.sqrt(5) * Numeric.array(T22red)
Hss = TriToMatrix(Hss)
Hss = TermReduce.expandSL(ion.term, Hss, 2, MATRIX_H5)



#-------------------------------------------------------------------------
# Comparison of results:

diff = Numeric.absolute(Numeric.absolute(myHss)-Numeric.absolute(Hss))
maximum = max(Numeric.ravel(diff))
names = TermNames.short(ion.term, "slj")

#showMatrix(myHss, names)
#showMatrix(Hss, names)
#showMatrix(diff)

print "Maximum difference of Hss (f12) to Judd: %g" % maximum
