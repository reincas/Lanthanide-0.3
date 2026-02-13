# -*- coding: latin-1 -*-
# Copyright 2000,2003 Reinhard Caspary <r.caspary@tu-bs.de>
#
# This file is part of the Python package Lanthanide.
# 
# Lanthanide is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# Lanthanide is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Lanthanide; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from Lanthanide import *
from Utilities import *
import TermReduce


######################################################################

def Test0Matrix(config, term, key, opts, rw):
    l = config["l"]
    k = opts[0]
    q = opts[1]
    m1 = Matrix(config, term, key, "U1:%i:%i" % (k,+q), rw)
    m2 = Matrix(config, term, key, "U1:%i:%i" % (k,-q), rw)
    matrix = m1 + m2
    matrix *= 3
    #if q == 2:
        #matrix *= math.sqrt(7.0)/2
        #matrix *= 1000000000
    #if q == 0:
        #matrix *= math.sqrt(35.0/2)
    matrix = MatrixToSparse(matrix)
    return matrix

def Test1Matrix(config, term, key, opts, rw):
    l = config["l"]
    k = opts[0]
    matrix = SparseToTri(Matrix(config, term, key, "U1U1:%i" % k, rw))
    multadd(Matrix(config, term, key, "U1U2:%i" % k, rw), 2, matrix)
    matrix = TriToSparse(matrix)
    #matrix["data"] *= 1
    return matrix

######################################################################

def juddfactor(c):
    a,b = c
    if a < 0:
        return -math.sqrt(-1.0*a/b)
    return math.sqrt(1.0*a/b)


    
juddtable = [[
    ( 2, 2, 2, 1 ),
    (     -11,     1134),
    (     605,     5292),
    (   32761,   889056),
    (    3575,   889056),
    (  -17303,   396900),
    (   -1573,     8232),
    (  264407,   823200),
    (   21879,   274400),
    (  -46189,   231525),
    ],[
    ( 2, 2, 4, 3 ),
    (       4     , 189),
    (   -6760,    43659),
    (      33,     1372),
    (    -325,    37044),
    (     416,    33075),
    (  -15028,   305613),
    (   28717,  2778300),
    (  -37349,   926100),
    (   -8398,   694575),
    ],[
    ( 2, 4, 4, 3 ),
    (       1,      847),
    (   -1805,   391314),
    (      -4,    33957),
    (  -54925,   373527),
    (    -117,   296450),
    (    4693, 12326391),
    (-1273597, 28014525),
    (  849524,  9338175),
    ( -134368,  3112725),
    ],[
    ( 2, 4, 6, 6 ),
    (      26,     3267),
    (   -4160,   754677),
    (     -13,      264),
    (     625,    26136),
    (     256,   571725),
    (    1568,   107811),
    (     841,  1960200),
    (     -17,   653400),
    (  -15827,   245025),
    ],[
    ( 4, 4, 4, 1 ),
    (   -6877,   139755),
    (   55016,   717409),
    (   49972,   622545),
    (   92480,  1369599),
    (  178802,   978285),
    ( -297680,  5021863),
    ( -719104,  2282665),
    (  -73644,  2282665),
    (   -2584,    18865),
    ],[
    ( 4, 4, 6, 3 ),
    (     117,     1331),
    (    -195,   204974),
    (      52,     1089),
    (     529,    11979),
    (   -2025,    18634),
    (     -49,   395307),
    (   -1369,    35937),
    (      68,    11979),
    (       0,        1),
    ],[
    ( 2, 6, 6, 3 ),
    (    2275,    19602),
    (    1625,   143748),
    (     325,   199584),
    (    6889,  2195424),
    (   71*71,  198*198),
    (      -1,   223608),
    (     625,    81312),
    (    1377,    27104),
    (     323,    22869),
    ],[
    ( 4, 6, 6, 3 ),
    (   12376,   179685),
    (   88400,  1185921),
    (    -442,    12705),
    (  -10880,   251559),
    (   -1088,   179685),
    ( -174080,  8301447),
    (   -8704,  3773385),
    ( -103058,  1257795),
    (     -19,    31185),
    ],[
    ( 6, 6, 6, 1 ),
    (    4199,   539055),
    (   29393,   790614),
    (  205751,   784080),
    (  -79135,  1724976),
    (    2261,  1078110),
    (   79135,   175692),
    (   15827,   319440),
    (   -8379,   106480),
    (     -98,     1485),
    ]]


def RkMatrix(config, term, key, opts, rw):
    l = config["l"]
    states = term[key]["states"]
    dim = opts[0]

    sum = NewTri(states)
    for k in range(1, dim, 2):
        matrix = Matrix(config, term, key, "U1U1:%i" % k, rw)
        multadd(matrix, (2*k+1), sum)
        matrix = Matrix(config, term, key, "U1U2:%i" % k, rw)
        multadd(matrix, 2*(2*k+1), sum)
    sum /= dim - 2.0
    return TriToSparse(sum)

def G2Matrix(config, term, key, opts, rw):
    states = term[key]["states"]
    sum = NewTri(states)
    matrix = Matrix(config, term, key, "U1U1:1", rw)
    multadd(matrix, 3, sum)
    matrix = Matrix(config, term, key, "U1U2:1", rw)
    multadd(matrix, 2*3, sum)
    matrix = Matrix(config, term, key, "U1U1:5", rw)
    multadd(matrix, 11, sum)
    matrix = Matrix(config, term, key, "U1U2:5", rw)
    multadd(matrix, 2*11, sum)
    sum /= 4.0
    return TriToSparse(sum)

def ekMatrix(config, term, key, opts, rw):
    l = config["l"]
    states = term[key]["states"]
    if l != 3:
        raise RuntimeError, "ekMatrix: Rules available only for l=3!"
    i = opts[0]
    if not i in (1,2,3):
        raise RuntimeError, "ekMatrix: Option must be 1, 2, or 3!"
    sum = NewTri(states)
    if i == 1:
        multadd(h1Matrix(config, term, key, [2], rw),     75.0/14,  sum)
        multadd(h1Matrix(config, term, key, [4], rw),     99.0/7,   sum)
        multadd(h1Matrix(config, term, key, [6], rw),   5577.0/350, sum)
    elif i == 2:
        multadd(h1Matrix(config, term, key, [2], rw),  10725.0/14,  sum)
        multadd(h1Matrix(config, term, key, [4], rw), -12870.0/7,   sum)
        multadd(h1Matrix(config, term, key, [6], rw), 195195.0/350, sum)
    else:
        multadd(h1Matrix(config, term, key, [2], rw),    825.0/14,  sum)
        multadd(h1Matrix(config, term, key, [4], rw),    396.0/7,   sum)
        multadd(h1Matrix(config, term, key, [6], rw), -39039.0/350, sum)
    return TriToSparse(sum)
        
def DkqMatrix(config, term, key, opts, rw):
    l = config["l"]
    k = opts[0]
    q = opts[1]
    matrix = Matrix(config, term, key, "U1:%i:%i" % (k,q), rw)
    matrix *= SIGNEXP(l) * (2*l+1) * wign3j(l,k,l,0,0,0)
    return matrix

def DqMatrix(config, term, key, opts, rw):
    l = config["l"]
    if l != 2:
        raise RuntimeError, "DMatrix: Rules available only for l=2!"
    m1 = DkqMatrix(config, term, key, [4,0], rw)
    m2 = DkqMatrix(config, term, key, [4,+4], rw)
    m3 = DkqMatrix(config, term, key, [4,-4], rw)
    matrix = 21 * (m1 + math.sqrt(5.0/14) * (m2+m3))
    matrix = MatrixToSparse(matrix)
    return matrix
        
def CMatrix(config, term, key, opts, rw):
    l = config["l"]
    states = term[key]["states"]
    sum = NewTri(states)
    multadd(h1Matrix(config, term, key, [2], rw), 7.0, sum)
    multadd(h1Matrix(config, term, key, [4], rw), 63.0/5, sum)
    return TriToSparse(sum)

def BMatrix(config, term, key, opts, rw):
    l = config["l"]
    matrix = h1Matrix(config, term, key, [2], rw)
    matrix["data"] *= 49.0
    return matrix

def h1Matrix(config, term, key, opts, rw):
    l = config["l"]
    k = opts[0]
    matrix = Matrix(config, term, key, "U1U2:%i" % k, rw)
    f = (2*l+1) * wign3j(l,k,l,0,0,0)
    matrix["data"] *= f * f
    return matrix

def h2Matrix(config, term, key, opts, rw):
    l = config["l"]
    matrix = Matrix(config, term, key, "U1T1:1", rw)
    matrix["data"] *= math.sqrt(1.5 * l * (l+1) * (2*l+1))
    return matrix

def h3Matrix(config, term, key, opts, rw):
    l = config["l"]
    i = opts[0]
    if not l in (1,2,3):
        raise RuntimeError, "h3Matrix: Rules available only for l = 1,2,3!"
    if not i in range(l):
        raise RuntimeError, "h3Matrix: Option must be < l!"
    if l == 1:
        matrix = L2Matrix(config, term, key, [], rw)
    elif l == 2:
        if i == 0:
            matrix = L2Matrix(config, term, key, [], rw)
        else:
            matrix = RkMatrix(config, term, key, [5], rw)
    else:
        if i == 0:
            matrix = L2Matrix(config, term, key, [], rw)
        elif i == 1:
            matrix = G2Matrix(config, term, key, [], rw)
        else:
            matrix = RkMatrix(config, term, key, [7], rw)
    return matrix

def VkkkMatrix(config, term, key, opts, rw):
    k1 = opts[0]
    k2 = opts[1]
    k3 = opts[2]
    if (k1==k2) and (k2==k3):
        matrix = Matrix(config, term, key, "U1U2U3:%i" % k1, rw)
    elif k1==k2:
        matrix = Matrix(config, term, key, "U1U2U3:%i:%i" % (k1,k3), rw)
    elif k2==k3:
        matrix = Matrix(config, term, key, "U1U2U3:%i:%i" % (k3,k1), rw)
    else:
        matrix = Matrix(config, term, key, "U1U2U3:%i:%i:%i" % (k1,k2,k3), rw)
    matrix["data"] *= 6 * math.sqrt((2*k1+1)*(2*k2+1)*(2*k3+1))
    return matrix

def h4Matrix(config, term, key, opts, rw):
    l = config["l"]
    states = term[key]["states"]
    if l != 3:
        raise RuntimeError, "h4Matrix: Rules available only for l=3!"
    i = opts[0]
    sum = NewTri(states)
    for j in range(len(juddtable)):
        k1,k2,k3,num = juddtable[j][0]
        factor = num * juddfactor(juddtable[j][i])
        matrix = VkkkMatrix(config, term, key, [k1, k2, k3], rw)
        multadd(matrix, factor, sum)
    return TriToSparse(sum)

def hssMatrix(config, term, key, opts, rw):
    l = config["l"]
    k = opts[0]
    ck0 = -(2*l+1) * wign3j(l,k,l,0,0,0)
    ck2 = -(2*l+1) * wign3j(l,k+2,l,0,0,0)
    f = -12 * ck0*ck2 * math.sqrt((k+1)*(k+2)*(2*k+1)*(2*k+3)*(2*k+5)/5.0)
    matrix = Matrix(config, term, key, "U1U2T1T2:%i:%i:1:1:2" % (k,k+2), rw)
    matrix["data"] *= f
    return matrix

def hsooMatrix(config, term, key, opts, rw):
    l = config["l"]
    states = term[key]["states"]
    k = opts[0]
    sum = NewTri(states)
    ck0 = -(2*l+1) * wign3j(l,k,l,0,0,0)
    a = -ck0*ck0 * math.sqrt((2*l+k+2)*(2*l-k)*(k+1)*(2*k+1)*(2*k+3))
    matrix = Matrix(config, term, key, "U1U2T1T2:%i:%i:0:1:1" % (k,k+1), rw)
    multadd(matrix, a, sum)
    matrix = Matrix(config, term, key, "U1U2T1T2:%i:%i:0:1:1" % (k+1,k), rw)
    multadd(matrix, 2*a, sum)
    ck2 = -(2*l+1) * wign3j(l,k+2,l,0,0,0)
    b = -ck2*ck2 * math.sqrt((2*l+k+3)*(2*l-k-1)*(k+2)*(2*k+3)*(2*k+5))
    matrix = Matrix(config, term, key, "U1U2T1T2:%i:%i:0:1:1" % (k+2,k+1), rw)
    multadd(matrix, b, sum)
    matrix = Matrix(config, term, key, "U1U2T1T2:%i:%i:0:1:1" % (k+1,k+2), rw)
    multadd(matrix, 2*b, sum)
    matrix = TriToSparse(2 * sum)
    return matrix

def h5Matrix(config, term, key, opts, rw):
    hss  = SparseToTri(hssMatrix(config, term, key, opts, rw))
    hsoo = SparseToTri(hsooMatrix(config, term, key, opts, rw))
    return TriToSparse(hss + hsoo)

def h5fixMatrix(config, term, key, opts, rw):
    l = config["l"]
    if l != 3:
        raise RuntimeError, "h5fixMatrix: Rules available only for l=3!"
    sum = h5Matrix(config, term, key, [0], rw)
    sum = SparseToTri(sum)
    matrix = h5Matrix(config, term, key, [2], rw)
    multadd(matrix, 0.56, sum)
    matrix = h5Matrix(config, term, key, [4], rw)
    multadd(matrix, 0.38, sum)
    return TriToSparse(sum)

def h6Matrix(config, term, key, opts, rw):
    l = config["l"]
    states = term[key]["states"]
    k = opts[0]
    ck = -(2*l+1) * wign3j(l,k,l,0,0,0)
    a = ck*ck / 6
    sum = NewTri(states)
    if k > 0:
        b = math.sqrt((2*l+k+1)*(2*l-k+1)*k*(2*k-1)/(2*k+1))
        matrix = Matrix(config, term, key,
                        "U1U2T1T2:%i:%i:0:1:1" % (k,k-1), rw)
        multadd(matrix, b, sum)
    if k < 2*l:
        c = -math.sqrt((2*l+k+2)*(2*l-k)*(k+1)*(2*k+3)/(2*k+1))
        matrix = Matrix(config, term, key,
                        "U1U2T1T2:%i:%i:0:1:1" % (k,k+1), rw)
        multadd(matrix, c, sum)
    sum *= a
    matrix = TriToSparse(2 * sum)
    return matrix

def h6fixMatrix(config, term, key, opts, rw):
    l = config["l"]
    if l != 3:
        raise RuntimeError, "h6fixMatrix: Rules available only for l=3!"
    sum = h6Matrix(config, term, key, [2], rw)
    sum = SparseToTri(sum)
    matrix = h6Matrix(config, term, key, [4], rw)
    multadd(matrix, 0.75, sum)
    matrix = h6Matrix(config, term, key, [6], rw)
    multadd(matrix, 0.50, sum)
    return TriToSparse(sum)

#def cfMatrix(config, term, key, opts, rw):
#    l = config["l"]
#    k = opts[0]
#    q = opts[1]
#    matrix = Matrix(config, term, key, "U1:%i:%i" % (k,q), rw)
#    if q != 0:
#        matrix *= SIGNEXP(q) * \
#                  Matrix(config, term, key, "U1:%i:%i" % (k,-q), rw)
#    matrix = MatrixToSparse(matrix)
#    matrix["data"] *= SIGNEXP(l) * (2*l+1) * wign3j(l,k,l,0,0,0)
#    return matrix

def cfMatrix(config, term, key, opts, rw):
    l = config["l"]
    k = opts[0]
    q = abs(opts[1])
    matrix = Matrix(config, term, key, "U1:%i:%i" % (k,q), rw)
    matrix = MatrixToSparse(matrix)
    matrix["data"] *= SIGNEXP(l) * (2*l+1) * wign3j(l,k,l,0,0,0)
    return matrix

def bzMatrix(config, term, key, opts, rw):
    lz = SparseToTri(LzMatrix(config, term, key, opts, rw))
    sz = SparseToTri(SzMatrix(config, term, key, opts, rw))
    return TriToSparse(0.466945406392*(lz + 2.00231924*sz))

def L2Matrix(config, term, key, opts, rw):
    l = config["l"]
    matrix = SparseToTri(Matrix(config, term, key, "U1U1:1", rw))
    multadd(Matrix(config, term, key, "U1U2:1", rw), 2, matrix)
    matrix = TriToSparse(matrix)
    matrix["data"] *= l * (l+1.0) * (2.0*l+1.0)
    return matrix

def LzMatrix(config, term, key, opts, rw):
    l = config["l"]
    matrix = Matrix(config, term, key, "U1:1:0", rw)
    matrix = MatrixToSparse(matrix)
    matrix["data"] *= math.sqrt(l * (l+1) * (2*l+1))
    return matrix

def S2Matrix(config, term, key, opts, rw):
    matrix = SparseToTri(Matrix(config, term, key, "T1T1:1", rw))
    multadd(Matrix(config, term, key, "T1T2:1", rw), 2, matrix)
    matrix = TriToSparse(matrix)
    matrix["data"] *= 1.5
    return matrix

def SzMatrix(config, term, key, opts, rw):
    matrix = Matrix(config, term, key, "T1:1:0", rw)
    matrix = MatrixToSparse(matrix)
    matrix["data"] *= math.sqrt(1.5)
    return matrix

def J2Matrix(config, term, key, opts, rw):
    l = config["l"]
    l2 = SparseToTri(L2Matrix(config, term, key, opts, rw))
    s2 = SparseToTri(S2Matrix(config, term, key, opts, rw))
    ls = SparseToTri(Matrix(config, term, key, "U1T1:1", rw))
    multadd(Matrix(config, term, key, "U1T2:1", rw), 2, ls)
    ls *= math.sqrt(1.5 * l * (l+1) * (2*l+1))
    return TriToSparse(l2 + 2*ls + s2)

def JzMatrix(config, term, key, opts, rw):
    lz = SparseToTri(LzMatrix(config, term, key, opts, rw))
    sz = SparseToTri(SzMatrix(config, term, key, opts, rw))
    return TriToSparse(lz + sz)


K = "([0-9]+)"
Q = "([+-]?[0-9]+)"

MatrixList = [
    [ "^Test0:%s:%s$" % (K,K),     Test0Matrix  ],
    [ "^Test1:%s$" % K,            Test1Matrix  ],
    [ "^H1:%s$" % K,               h1Matrix     ],
    [ "^H2$",                      h2Matrix     ],
    [ "^H3:%s$" % K,               h3Matrix     ],
    [ "^H4:%s$" % K,               h4Matrix     ],
    [ "^H5:%s$" % K,               h5Matrix     ],
    [ "^H5fix",                    h5fixMatrix  ],
    [ "^H6:%s$" % K,               h6Matrix     ],
    [ "^H6fix",                    h6fixMatrix  ],
    [ "^cf:%s:%s$" % (K,K),        cfMatrix     ],
    [ "^e%s$" % K,                 ekMatrix     ],
    [ "^R%s$" % K,                 RkMatrix     ],
    [ "^G2$",                      G2Matrix     ],
    [ "^Bz$",                      bzMatrix     ],
    [ "^L2$",                      L2Matrix     ],
    [ "^Lz$",                      LzMatrix     ],
    [ "^S2$",                      S2Matrix     ],
    [ "^Sz$",                      SzMatrix     ],
    [ "^J2$",                      J2Matrix     ],
    [ "^Jz$",                      JzMatrix     ],
    [ "^B$",                       BMatrix      ],
    [ "^C$",                       CMatrix      ],
    [ "^Dq$",                      DqMatrix     ],
    [ "^Hss:%s$" % K,              hssMatrix    ],
    [ "^Hsoo:%s$" % K,             hsooMatrix   ],
    [ "^Vkkk:%s:%s:%s$" % (K,K,K), VkkkMatrix   ],
]

UnitList = [
    [ "^OPE$",                                   OPE,       -1 ],
    [ "^U1:%s:%s$" % (K,Q),                      MAT1_U,     1 ],
    [ "^T1:%s:%s$" % (K,Q),                      MAT1_T,     1 ],
    [ "^U1U1:%s$" % K,                           MAT1_UU,   -1 ],
    [ "^U1U1:%s:%s:%s$" % (K, K, Q),             MAT1_UUS,  -1 ],
    [ "^T1T1:%s$" % K,                           MAT1_TT,   -1 ],
    [ "^U1T1:%s$" % K,                           MAT1_UT,   -1 ],
    [ "^U1U2:%s$" % K,                           MAT2_UU,   -1 ],
    [ "^U1U2:%s:%s:%s$" % (K, K, Q),             MAT2_UUS,  -1 ],
    [ "^T1T2:%s$" % K,                           MAT2_TT,   -1 ],
    [ "^U1T2:%s$" % K,                           MAT2_UT,   -1 ],
    [ "^U1U2T1T2:%s:%s:%s:%s:%s$" % (K,K,K,K,K), MAT2_UUTT, -1 ],
    [ "^U1U2U3:%s$" % K,                         MAT3_UUUA, -1 ],
    [ "^U1U2U3:%s:%s$" % (K,K),                  MAT3_UUUB, -1 ],
    [ "^U1U2U3:%s:%s:%s$" % (K,K,K),             MAT3_UUUC, -1 ],
    ]



def matrixOpts(optlist, name):
    for item in optlist:
        pattern = item[0]
        if re.match(pattern, name):
            opts = re.split(pattern, name)[1:-1]
            for i in range(len(opts)):
                opts[i] = string.atoi(opts[i])
            break
    else:
        raise RuntimeError, "calcMatrix: Unknown matrix: %s\n" % name

    list = [ opts ]
    for i in range(1,len(item)):
        list.append(item[i])
    return list



def calcUnit(config, name):
    opts,type,qpos = matrixOpts(UnitList, name)
    quant = config["quant"]
    unit = config["unit"]
    if (qpos < 0):
        matrix = calcmatrix(quant, unit, config, type, opts)
    else:
        q = opts[qpos]
        if q >= 0:
            tri = "u"
        else:
            opts[qpos] = -q
            tri = "d"
        matrix = calcmatrix(quant, unit, config, type, opts)
        matrix = SparseToMatrix(matrix, tri)
        if q < 0:
            matrix *= SIGNEXP(q)
    return matrix



def calcMatrix(config, term, key, name, rw):
    opts,func = matrixOpts(MatrixList, name)
    return func(config, term, key, opts, rw)



def Matrix(config, term, key, name, rw=IO_READ+IO_WRITE):
    fn = dataPath(config, key, name)
    if (rw & IO_READ) and os.path.exists(fn):
        fp = gzip.open(fn, "r")
        matrix = cPickle.loads(fp.read())
        fp.close()
        if matrix["DictType"] == DICT_SMATRIX:
            matrix = SmatrixToMatrix(matrix)
    else:
        print "Matrix: calc %s" % fn
        if re.match("^[UT]1$", name[:2]) or name == "OPE":
            if key == "product":
                matrix = calcUnit(config, name)
            elif key == "sljm":
                matrix = Matrix(config, term, "product", name, rw)
                if type(matrix) == type({}):
                    matrix = transform(matrix, term["vectors"])
                else:
                    matrix = Transform(matrix, term["vectors"])
            elif key == "slo":
                matrix = Matrix(config, term, "product", name, rw)
                if type(matrix) == type({}):
                    matrix = transform(matrix, term["slovectors"])
                else:
                    matrix = Transform(matrix, term["slovectors"])
            elif key == "slj":
                matrix = Matrix(config, term, "sljm", name, rw)
                matrix = TermReduce.prepareSLJ(term, matrix)
            else:
                raise RuntimeError, "Matrix: no rules for term %s!" % key
        else:
            matrix = calcMatrix(config, term, key, name, rw)
        if rw & IO_WRITE:
            ##########################################################
            ### Dirty hack for cygwin-python:                      ###
            ### str() should be not necessary!                     ###
            if str(type(matrix)) != str(type({})):
                m = MatrixToSmatrix(matrix)
            else:
                m = matrix
            fp = gzip.open(fn, "w")
            fp.write(cPickle.dumps(m))
            fp.close()
    return matrix
