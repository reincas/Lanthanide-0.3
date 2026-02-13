# (c) 01/2001 by Reinhard Caspary

from Lanthanide import *
from Ion import *
import TermNames
import TermReduce



parms = [
    [ -3,  "H4:1", 1 ],
    [ 11,  "H4:2", 1 ],
    [ -13, "H4:3", 1 ],
    [  5,  "H4:4", 1 ],
    [ -1,  "H4:5", 1 ],
    [  2,  "H4:6", 1 ],
    [ -5,  "H4:7", 1 ],
    [ -7,  "H4:8", 1 ],
    [  7,  "H4:9", 1 ],
]



#-------------------------------------------------------------------------
# My own calculations:

ion = Ion("Nd3+", parms)
ion.Term()
myTi = ion.hamiltonBuild()
myTi = TriToMatrix(myTi)



#-------------------------------------------------------------------------
# Table from Judd:

t1 = parms[0][0]
t2 = parms[1][0]
t3 = parms[2][0]
t4 = parms[3][0]
t5 = parms[4][0]
t6 = parms[5][0]
t7 = parms[6][0]
t8 = parms[7][0]
t9 = parms[8][0]

s22 = math.sqrt(22)
s33 = math.sqrt(33)
s455 = math.sqrt(455)
s4290 = math.sqrt(4290)

elements = Numeric.array([
    [   6,       0, 288,         0,        0,         0,        0,     0, 0 ],
    [   6,    1694,   8,     -8008,        0,         0,        0,     0, 0 ],
    [   6,       0, -72,         0,        0,         0,        0,     0, 0 ],
    [   6,     616,   8,      7280,        0,         0,        0,     0, 0 ],
    [   6,   -1078,   8,     -1960,        0,         0,        0,     0, 0 ],
    [  -1,    -385, -48,         0,        0,    -30030,        0,     0, 0 ],
    [  -1,    -319,  32,     -1144,      286,     12870,    10296,     0, 0 ],
    [   0,  36*s33,   0,   468*s33,  156*s33,  -624*s33,  156*s33,     0, 0 ],
    [  -1,    -423,  -3,      3237,     -377,     -1677,    -1833,  4641, 0 ],
    [ -15,       0,   0,         0,        0,         0,        0,     0, 0 ],
    [   0, 231*s22,   0,         0,        0,         0,        0,     0, 0 ],
    [  -1,     -21,  -3,      1365,     -455,      1365,    -1365, -3315, 0 ],
    [  -1,    -116,  32,      1040,     -260,      4680,    -9360,     0, 0 ],
    [   0, 3*s4290,   0, -24*s4290, -8*s4290, -52*s4290, -8*s4290,     0, 0 ],
    [  -1,      11,  -3,     -2475,      561,      1221,     1947,  1309, 0 ],
    [  -1,     105, -48,         0,        0,      8190,        0,     0, 0 ],
    [   0,       0,   0,   84*s455, -28*s455,         0, 252*s455,     0, 0 ],
    [  -1,    -399,  -3,     -1995,      -49,     -2709,      567, -1071, 0 ],
    [  -1,     203,  32,      -280,       70,     -8190,     2520,     0, 0 ],
    [  -1,      56,  -3,      1827,      315,      -252,       21, -1071, 0 ],
    [  -1,     336,  -3,      -525,     -245,      1260,     -315,   945, 0 ],
    ], Numeric.Float)

parms = Numeric.array([
    t1 * math.sqrt(33.0/6860),
    t2 * math.sqrt(2) / 2156.0,
    t3 * 1 / math.sqrt(6720),
    t4 * 1 / (56.0 * math.sqrt(15015)),
    t5 * 3 / (49.0 * math.sqrt(17160)),
    t6 * 1 / (924.0 * math.sqrt(455)),
    t7 * 1 / (168.0 * math.sqrt(5005)),
    t8 * 1 / math.sqrt(16336320),
    t9,
    ], Numeric.Float)

elements = Numeric.matrixmultiply(elements, parms)

lsterms = [
    ("4S", "4S"),
    ("4D", "4D"),
    ("4F", "4F"),
    ("4G", "4G"),
    ("4I", "4I"),
    ("2P", "2P"),
    ("2D(1)", "2D(1)"),
    ("2D(1)", "2D(2)"),
    ("2D(2)", "2D(2)"),
    ("2F(1)", "2F(1)"),
    ("2F(1)", "2F(2)"),
    ("2F(2)", "2F(2)"),
    ("2G(1)", "2G(1)"),
    ("2G(1)", "2G(2)"),
    ("2G(2)", "2G(2)"),
    ("2H(1)", "2H(1)"),
    ("2H(1)", "2H(2)"),
    ("2H(2)", "2H(2)"),
    ("2I", "2I"),
    ("2K", "2K"),
    ("2L", "2L"),
    ]

key = "sl"
names = TermNames.short(ion.term, key)
states = ion.term[key]["states"]

Tired = Numeric.zeros((states, states), Numeric.Float)
for i in range(len(elements)):
    a = names.index(lsterms[i][0])
    b = names.index(lsterms[i][1])
    val = elements[i]
    Tired[a,b] = val
    Tired[b,a] = val
Ti = -math.sqrt(3) * Tired
Ti = TermReduce.expandSL(ion.term, Ti, 1, MATRIX_H4)



#-------------------------------------------------------------------------
# Comparison of results:

diff = Numeric.absolute(Numeric.absolute(myTi)-Numeric.absolute(Ti))
maximum = max(Numeric.ravel(diff))
names = TermNames.short(ion.term, "slj")

#span = (0,20)
#showMatrix(myTi, names, span=span)
#showMatrix(Ti, names, span=span)
#showMatrix(diff)

print "Maximum difference of Ti (f03) to Judd: %g" % maximum
