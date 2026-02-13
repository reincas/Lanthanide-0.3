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



#-------------------------------------------------------------------------
# My own calculations:

parms = [
    [  2, "Hsoo:0", 1 ],
    [ -3, "Hsoo:2", 1 ],
    [  1, "Hsoo:4", 1 ],
]

ion = Ion("Tm3+", parms)
ion.Term()
myHsoo = ion.hamiltonBuild()
myHsoo = TriToMatrix(myHsoo)



#-------------------------------------------------------------------------
# Table from Carnall:

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
    138*m0 -10*m2 -(50/11.0)*m4, # 3P
    -math.sqrt(1/30.0)*(1044*m0 -62*m2 -20*m4),
    0,
    0,
    30*m0 -78*m2 -(930/11.0)*m4,
    0, # 3F
    math.sqrt(10)*((353/5.0)*m0 -(24/5.0)*m2 -(69/11.0)*m4),
    math.sqrt(11)*(-72*m0 +(262/33.0)*m2 -(250/363.0)*m4),
    0,
    0,
    math.sqrt(14)*(36*m0 -8*m2 -(10/11.0)*m4),
    0, # 3H
    0,
    math.sqrt(1/10.0)*(738*m0 -(3436/33.0)*m2 -(16250/363.0)*m4),
    math.sqrt(26)*(-38*m0 +(3/11.0)*m2 +(1770/1573.0)*m4),
    0,
    0,
    math.sqrt(1/55.0)*(2574*m0 -146*m2 -(610/11.0)*m4),
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

print "Maximum difference of Hsoo (f12) to Carnall: %g" % maximum
