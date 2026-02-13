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



#-------------------------------------------------------------------------
# My own calculations:

parms = [
    [  2, "Hsoo:0", 1 ],
    [ -3, "Hsoo:2", 1 ],
    [  1, "Hsoo:4", 1 ],
]

ion = Ion("Pr3+", parms)
ion.Term()
myHsoo = ion.hamiltonBuild()
myHsoo = TriToMatrix(myHsoo)



#-------------------------------------------------------------------------
# Table from Judd:

m0 = parms[0][0]
m2 = parms[1][0]
m4 = parms[2][0]

T11red = [
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
    6*m0 +2*m2 +(10/11.0)*m4, # 3P
    -math.sqrt(2/15.0)*(27*m0 +14*m2 +(115/11.0)*m4),
    0,
    0,
    -36*m0 -72*m2 -(900/11.0)*m4,
    0, # 3F
    math.sqrt(2/5.0)*(23*m0 +6*m2 -(195/11.0)*m4),
    math.sqrt(11)*(-6*m0 +(64/33.0)*m2 -(1240/363.0)*m4),
    0,
    0,
    2*math.sqrt(14)*(-15*m0 -m2 +(10/11.0)*m4),
    0, # 3H
    0,
    math.sqrt(2/5.0)*(39*m0 -(728/33.0)*m2 -(3175/363.0)*m4),
    math.sqrt(26)*(-5*m0 -(30/11.0)*m2 -(375/1573.0)*m4),
    0,
    0,
    8/math.sqrt(55.0)*(-132*m0 +23*m2 +(130/11.0)*m4),
    ]
Hsoo = -math.sqrt(3) * Numeric.array(T11red)
Hsoo = TriToMatrix(Hsoo)
Hsoo = TermReduce.expandSL(ion.term, Hsoo, 1, MATRIX_H5)



#-------------------------------------------------------------------------
# Comparison of results:

diff = Numeric.absolute(Numeric.absolute(myHsoo)-Numeric.absolute(Hsoo))
maximum = max(Numeric.ravel(diff))
names = TermNames.short(ion.term, "slj")

#showMatrix(myHsoo, names)
#showMatrix(Hsoo, names)
#showMatrix(diff)

print "Maximum difference of Hsoo (f02) to Judd: %g" % maximum
