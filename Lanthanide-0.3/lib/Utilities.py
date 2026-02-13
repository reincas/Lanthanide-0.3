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

PATT_S = "([0-9]+)"
PATT_L = "([A-Z](\([0-9]+\))?)"
PATT_J = "([0-9]+(/[0-9]+)?)"
PATT_TERM = PATT_S + PATT_L + PATT_J

def splitShortTerm(term):
    vals = re.split(PATT_TERM, term)
    vals = (vals[1], vals[2], vals[4])
    return vals

# ----------------------------------------------------------------------
def parmFit(ion, element, parent, stage, method,
            meas, fmeas, cauchy, bandspan=[], modified=0):
    form     = "energy/%s/stage-%s.%s"
    parmMask = form % (element, "%s", "parm")
    fitName  = form % (element, stage, "fit")
    parmName = form % (element, stage, "parm")

    parms    = getParms(parent, parmMask)
    base     = BaseParms[element]
    parms    = setParms(parms, base, stage)

    ion.termParms(parms)
    ion.oscCauchy(cauchy)
    if bandspan:
        ion.setBandspan(bandspan)
    nparms, dk, omega, df = ion.fullFit(fitName, meas, fmeas, method)
    if bandspan:
        fmeas = ion.JO.meas
    saveParms(parmName, parms, nparms, dk, omega, df)

def saveParms(fn, oldparms, parms, dk, omega, df):
    s = []
    s.append("parms = [\n")
    for val,name,flag in parms:
        s.append("    [ %g, \"%s\", %i ],\n" % (val, name, flag))
    s.append("    ]\n")
    s.append("dk = %g\n\n" % dk)
    s.append("oldparms = [\n")
    for val,name,flag in oldparms:
        s.append("    [ %g, \"%s\", %i ],\n" % (val, name, flag))
    s.append("    ]\n\n")
    s.append("omega = [\n")
    if len(omega) == 3:
        index = [2, 4, 6]
    else:
        index = range(len(omega))
    for i in range(len(omega)):
        name = "Omega:%i" % index[i]
        s.append("    [ %g, \"%s\" ],\n" % (omega[i], name))
    s.append("    ]\n")
    s.append("df = %g\n\n" % df)
        
    fp = open(fn, "w")
    fp.writelines(s)
    fp.close()

def getParms(stage, form):
    if not stage:
        parms = []
    else:
        name = form % stage
        fp = open(name, "r")
        exec(fp.read())
        fp.close()
    return parms

def setParms(parms, base, stage):
    if not stage:
        return parms
    if not STAGES.has_key(stage):
        raise RuntimeError, "setParms: Unknown stage %s!" % stage
    return STAGES[stage](parms, base)
    
def stage0(parms, base):
    parms = []
    parms = mergeParm(parms, base, "base")
    parms = mergeParm(parms, base, "H2")
    return parms

def stage1(parms, base):
    parms = []
    parms = mergeParm(parms, base, "base")
    parms = mergeParm(parms, base, "H1:2")
    parms = mergeParm(parms, base, "H1:4")
    parms = mergeParm(parms, base, "H1:6")
    parms = mergeParm(parms, base, "H2")
    return parms

def stage2(parms, base):
    parms = mergeParm(parms, base, "H3:0")
    parms = mergeParm(parms, base, "H3:1")
    parms = mergeParm(parms, base, "H3:2")

    parms = mergeParm(parms, base, "H4:2")
    parms = mergeParm(parms, base, "H4:3")
    parms = mergeParm(parms, base, "H4:4")
    parms = mergeParm(parms, base, "H4:6")
    parms = mergeParm(parms, base, "H4:7")
    parms = mergeParm(parms, base, "H4:8")
    return parms

def stage2a(parms, base):
    parms = mergeParm(parms, base, "H3:0")
    parms = mergeParm(parms, base, "H3:1")
    parms = mergeParm(parms, base, "H3:2")
    setParm(parms, "H3:2", 2, 0)
    return parms

def stage2b(parms, base):
    parms = mergeParm(parms, base, "H3:0")
    parms = mergeParm(parms, base, "H3:1")
    parms = mergeParm(parms, base, "H3:2")
    return parms

def stage3(parms, base):
    parms = mergeParm(parms, base, "H5fix")
    parms = mergeParm(parms, base, "H6fix")
    return parms

def stage3x(parms, base):
    parms = mergeParm(parms, base, "H5fix")
    parms = mergeParm(parms, base, "H6fix")
    setParm(parms, "H6fix", 2, 0)
    return parms

def stage3a(parms, base):
    h5 = getParm(base, "H5fix", 0) / 1.94
    parms = addParm(parms, [h5,      "H5:0", 1])
    parms = addParm(parms, [h5*0.56, "H5:2", 1])
    parms = addParm(parms, [h5*0.38, "H5:4", 1])
    h6 = getParm(base, "H6fix", 0) / 2.25
    parms = addParm(parms, [h6,      "H6:2", 1])
    parms = addParm(parms, [h6*0.75, "H6:4", 1])
    parms = addParm(parms, [h6*0.50, "H6:6", 1])
    return parms

STAGES = {
    "0":   stage0,
    "1":   stage1,
    "2":   stage2,
    "2a":  stage2a,
    "2b":  stage2b,
    "3":   stage3,
    "3x":  stage3x,
    "3a":  stage3a,
    }



# ----------------------------------------------------------------------
THRES = 1e-10

def ZERO(x, thres=THRES):
    return abs(x) < thres

def LIMIT(a, thres=THRES):
    return Numeric.where(Numeric.less_equal(a,-thres) + \
                         Numeric.greater_equal(a, thres), a, 0)

def SIGN(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0

def FORSIGN(x, y):
    if y >= 0:
        return abs(x)
    return -abs(x)

def INT(x):
    if x >= 0:
        return int(x+0.1)
    return int(x-0.1)

def SIGNEXP(x):
    if (abs(INT(x)) % 2) == 0:
        return 1.0
    return -1.0

def SWAP(x, y):
    return (y, x)



# ----------------------------------------------------------------------
def Nielson(list):
    if not list:
        return 0.0
    
    primeList = [ 1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 ]
    num = 1
    for i in range(1, len(list)):
        num *= pow(1.0*primeList[i], list[i])
    return list[0] * math.sqrt(num)

def NielsonUnit(term, list):
    key   = "sl"
    terms = term[key]["states"]
    syms  = term[key]["symmetry"]
    lvals = term[key]["values"][:,syms.index("L2")]
    
    matrix = Numeric.zeros((terms,terms), Numeric.Float)
    pos = 0
    for j in range(terms):
        lb = 0.5*lvals[j]
        for i in range(j+1):
            la = 0.5*lvals[j]
            matrix[i,j] = Nielson(list[pos])
            pos += 1
            if i != j:
                matrix[j,i] = SIGNEXP(la-lb) * matrix[i,j]
    return matrix



# ----------------------------------------------------------------------
def parmIndex(parms, name):
    for i in range(len(parms)):
	if parms[i][1] == name:
	    break
    else:
	raise RuntimeError, "parmIndex: Unknown parameter %s!" % name
    return i

def delParm(parms, name):
    i = parmIndex(parms, name)
    del parms[i]
    return parms

def setParm(parms, name, index, value):
    i = parmIndex(parms, name)
    parms[i][index] = value
    return parms

def getParm(parms, name, index):
    i = parmIndex(parms, name)
    return parms[i][index]

def addParm(parms, add):
    parms.append(add)
    return parms

def mergeParm(parms, merge, name):
    j = parmIndex(merge, name)
    try:
        i = parmIndex(parms, name)
    except:
        parms.append(merge[j])
    else:
        parms[i] = merge[j]
    return parms



# ----------------------------------------------------------------------
def CalcMatrix(quant, unit, config, key, opts):
    states = config["states"]
    matrix = Numeric.zeros((states, states), Numeric.Float)
    for i in range(states):
        for j in range(states):
            pair = orderpair(config, i, j)
            matrix[i,j] = calcelement(quant, unit, pair, key, opts)
    return matrix

def Diagonalize(dummy1, dummy2, matrix):
    matrix = TriToMatrix(matrix)
    w,v = LinearAlgebra.eigenvectors(matrix)
    index = Numeric.argsort(w)
    w = Numeric.take(w, index)
    v = Numeric.take(v, index)
    return w, v, 0

def Transform(matrix, eigenvectors):
    matrix = Numeric.matrixmultiply(eigenvectors, matrix)
    matrix = Numeric.matrixmultiply(matrix, Numeric.transpose(eigenvectors))
    return matrix



# ----------------------------------------------------------------------
def dataPath(config, term, name):
    l = "spdfghi"[config["l"]]
    dir = "%s" % DATADIR
    electrons = config["electrons"]
    dir += "/%c%02i" % (l, electrons)
    if term:
        dir += "/" + term
    if not os.path.exists(dir):
        os.makedirs(dir)
    name = string.replace(name, ":", ".")
    return "%s/%s.gz" % (dir, name)



# ----------------------------------------------------------------------
def showMatrix(matrix, rowNames=[], colNames=[], span=()):
    print MatrixToString(matrix, rowNames, colNames, span)
    
def MatrixToString(matrix, rowNames=[], colNames=[], span=()):
    if type(matrix) == type({}):
        matrix = SparseToMatrix(matrix)
        
    if len(matrix.shape) == 1:
        matrix = matrix[Numeric.NewAxis,:]
    rows = matrix.shape[0]
    cols = matrix.shape[1]

    if rowNames and not colNames:
        colNames = rowNames
    if not colNames:
        for i in range(cols):
            colNames.append("%i" % i)
    if not rowNames:
        for i in range(rows):
            rowNames.append("%i" % i)

    if span:
        first,last = span

        if first >= 0:
            firstrow = first
        else:
            firstrow = rows+first
        if last > 0:
            lastrow = last
        else:
            lastrow = rows+last

        if first >= 0:
            firstcol = first
        else:
            firstcol = cols+first
        if last > 0:
            lastcol = last
        else:
            lastcol = cols+last
    else:
        firstrow = 0
        lastrow = rows
        firstcol = 0
        lastcol = cols

    s = "        |"
    for j in range(firstcol,lastcol):
        s += " %5s" % colNames[j][:5]
    s += "\n"
    s += "--------+--" + "------"*(lastcol-firstcol)
    s += "\n"
    for i in range(firstrow,lastrow):
	s += " %5s  |" % rowNames[i][:5]
	for j in range(firstcol,lastcol):
	    if abs(matrix[i,j]) >= 0.0000001:
		s += " %5.2f" % matrix[i,j]
	    else:
		s += " %5s" % "0"
        s += "\n"
    return s[:-1]


def MatrixToLaTeX(matrix, rowNames=[], colNames=[], span=(), mode=1):
    if type(matrix) == type({}):
        matrix = SparseToMatrix(matrix)
        
    if len(matrix.shape) == 1:
        matrix = matrix[Numeric.NewAxis,:]
    rows = matrix.shape[0]
    cols = matrix.shape[1]

    if rowNames and not colNames:
        colNames = rowNames
    if not colNames:
        for i in range(cols):
            colNames.append("%i" % i)
    if not rowNames:
        for i in range(rows):
            rowNames.append("%i" % i)

    if span:
        first,last = span

        if first >= 0:
            firstrow = first
        else:
            firstrow = rows+first
        if last > 0:
            lastrow = last
        else:
            lastrow = rows+last

        if first >= 0:
            firstcol = first
        else:
            firstcol = cols+first
        if last > 0:
            lastcol = last
        else:
            lastcol = cols+last
    else:
        firstrow = 0
        lastrow = rows
        firstcol = 0
        lastcol = cols

    for i in range(firstcol,lastcol):
        if mode & 1:
            S,L,J = splitShortTerm(colNames[i])
            colNames[i] = "\\term{%s}{%s}{%s}" % (S,L,J)

    for i in range(firstrow,lastrow):
        if mode & 2:
            S,L,J = splitShortTerm(rowNames[i])
            rowNames[i] = "\\term{%s}{%s}{%s}" % (S,L,J)

    s = ""
    s += "\\begin{tabular}{|r|%s|}\n" % ("c"*(lastcol-firstcol))
    s += "  \\hline\n "
    for j in range(firstcol,lastcol):
        s += " & %s" % colNames[j]
    s += " \\\\\n"
    s += "  \\hline\n"
    for i in range(firstrow,lastrow):
	s += "  %s \n  " % rowNames[i]
        count = 0
	for j in range(firstcol,lastcol):
            value = matrix[i,j]
	    if abs(value) >= 0.000001:
                if value > 0:
                    s += "&%.2f" % value
                else:
                    s += "&$\\overline{%.2f}$" % (-value)
	    else:
		s += "&$\\cdot$"
            if count >= 5:
                s += "\n  "
                count = 0
            count += 1
        s += "\\\\\n"
    s += "  \\hline\n"
    s += "\\end{tabular}\n"
    return s[:-1]



# ----------------------------------------------------------------------
def UnitSparse(states):
    sparse = {"DictType": DICT_SPARSE}
    sparse["states"] = states
    sparse["size"] = 1
    sparse["data"] = Numeric.ones(states, Numeric.Float)
    sparse["col"] = Numeric.arrayrange(states)
    sparse["row"] = Numeric.arrayrange(states) + 1
    return sparse

def NewTri(states, typecode=Numeric.Float):
    elements = states * (states+1) / 2
    return Numeric.zeros(elements, typecode)

def SparseToTri(sparse):
    states = sparse["states"]
    elements = states * (states+1) / 2
    tri = Numeric.zeros(elements, Numeric.Float)
    index = 0
    for j in range(states):
        for i in range(j+1):
            tri[index] = getelement(sparse, i, j)
            index = index + 1
    return tri
    
def SparseToMatrix(sparse, tri=""):
    states = sparse["states"]
    matrix = Numeric.zeros((states,states), Numeric.Float)
    for j in range(states):
	for i in range(j+1):
	    value = getelement(sparse, i, j)
            if (tri == "") or (tri == "u"):
                matrix[i,j] = value
            if (tri == "") or (tri == "d"):
                matrix[j,i] = value
    return matrix

def TriToSparse(tri):
    states = (math.sqrt(1+8*tri.shape[0])-1)/2
    states = int(round(states,0))
    sparse = {"DictType": DICT_SPARSE}
    sparse["states"] = states
    if tri.typecode() == Numeric.Float:
	sparse["size"] = 1
    else:
	sparse["size"] = 2
    data = []
    col = []
    row = Numeric.zeros(states, Numeric.Int)
    index = 0
    ptr = 0
    for j in range(states):
	for i in range(j+1):
	    value = tri[ptr]
	    ptr = ptr + 1
	    if not ZERO(value):
		data.append(value)
		col.append(i)
		index = index + 1
	row[j] = index
    sparse["data"] = Numeric.array(data, Numeric.Float)
    sparse["col"] = Numeric.array(col, Numeric.Int)
    sparse["row"] = row
    return sparse

def TriToMatrix(tri):
    states = (math.sqrt(1+8*tri.shape[0])-1)/2
    states = int(round(states,0))
    matrix = Numeric.zeros((states,states), Numeric.Float)
    ptr = 0
    for j in range(states):
	for i in range(j+1):
	    value = tri[ptr]
	    ptr = ptr + 1
	    matrix[i,j] = value
	    matrix[j,i] = value
    return matrix

def MatrixToSparse(matrix):
    states = matrix.shape[0]
    sparse = {"DictType": DICT_SPARSE}
    sparse["states"] = states
    if matrix.typecode() == Numeric.Float:
	sparse["size"] = 1
    else:
	sparse["size"] = 2
    data = []
    col = []
    row = Numeric.zeros(states, Numeric.Int)
    index = 0
    for j in range(states):
	for i in range(j+1):
	    value = matrix[i,j]
	    if not ZERO(value):
		data.append(value)
		col.append(i)
		index = index + 1
	row[j] = index
    sparse["data"] = Numeric.array(data, Numeric.Float)
    sparse["col"] = Numeric.array(col, Numeric.Int)
    sparse["row"] = row
    return sparse

def MatrixToSmatrix(matrix):
    states = matrix.shape[0]
    smatrix = {"DictType": DICT_SMATRIX}
    smatrix["states"] = states
    if matrix.typecode() == Numeric.Float:
	smatrix["size"] = 1
    else:
	smatrix["size"] = 2
    data = []
    col = []
    row = Numeric.zeros(states, Numeric.Int)
    index = 0
    for j in range(states):
	for i in range(states):
	    value = matrix[i,j]
	    if not ZERO(value):
		data.append(value)
		col.append(i)
		index = index + 1
	row[j] = index
    smatrix["data"] = Numeric.array(data, Numeric.Float)
    smatrix["col"] = Numeric.array(col, Numeric.Int)
    smatrix["row"] = row
    return smatrix

def SmatrixToMatrix(smatrix):
    states = smatrix["states"]
    matrix = Numeric.zeros((states,states), Numeric.Float)
    row = smatrix["row"]
    col = smatrix["col"]
    data = smatrix["data"]
    first = 0
    for j in range(states):
        last = row[j]
        for pos in range(first,last):
            matrix[col[pos],j] = data[pos]
        first = last
    return matrix



# ----------------------------------------------------------------------
def statSparse(matrix, str=""):
    states = matrix["states"]
    total = matrix["row"][states-1]
    estimation = states * MEMFACTOR
    maximum = states * (states+1) / 2
    size = total * 12 + states * 4
    length = size
    dim = ""
    if length > 1024:
        length = length / 1024
        dim = "k"
    if length > 1024:
        length = length / 1024
        dim = "M"
    if str:
        s = " (%s):" % str
    else:
        s = ":"
    print "Matrix statistics%s" % s
    print "  states:    %7i" % states
    print "  memfactor: %7i" % MEMFACTOR
    print "  elements:  %7i" % maximum
    print "  nonzero:   %7i" % total
    print "  size:      %7.1f %sByte" % (length, dim)
    if total != 0:
        print "  estimation factor:  %10.2f" % (1.0*(states*MEMFACTOR)/total)
        print "  compression factor: %10.2f" % (1.0*maximum/total)
    print ""
    sys.stdout.flush()
    return 1.0*size
