from Lanthanide import *
from Utilities import *
from Ion import Ion
import TermNames
import Matrix
import TermReduce


P2U2 = [
    [],
    [1, 2, -1],
    [1, 0, -1, 0, 1],
    [],
    [],
    [-1],
    ]



F2U2 = [
    [],
    [1, 2, 0, 0, -1],
    [-1, -1, -1, 0, -2, 2],
    [],
    [1, 2, 0, 0, -2, 1],
    [1, 0, 3, 0, -2, -1],
    [],
    [],
    [1, 1, 0, 0, -1, -1, 1],
    [1, -1, -1, 2, 0, -1, 1],
    [],
    [],
    [],
    [],
    [-1, -1, 2, 0, -1],
    [],
    [],
    [],
    [],
    [1, 1, 1, 0, -1],
    [-1, 0, -2],
    [],
    [],
    [],
    [],
    [],
    [1, 2, -2, 0, -1, 1],
    [1, -1, -2, 0, -1, 1, 1],
    ]

def reduce(term, matrix, k, q):
    key = "sljm"
    terms  = term[key]["states"]
    syms   = term[key]["symmetry"]
    values = term[key]["values"]
    jvals  = values[:,syms.index("J2")]
    mvals  = values[:,syms.index("Jz")]
    
    for i in range(terms):
        for j in range(terms):
            ja = 0.5*jvals[i]
            ma = 0.5*mvals[i]
            jb = 0.5*jvals[j]
            mb = 0.5*mvals[j]
            if not ZERO(matrix[i,j]):
		factor = SIGNEXP(ja-ma) * wign3j(ja,k,jb,-ma,q,mb)
                matrix[i,j] /= factor
    return matrix



def showUnit(term, key, unit, k):
    names  = TermNames.short(term, key)
    terms  = term[key]["states"]
    syms   = term[key]["symmetry"]
    values = term[key]["values"]
    svals  = values[:,syms.index("S2")]
    lvals  = values[:,syms.index("L2")]
    jvals  = values[:,syms.index("J2")]

    for i in range(terms):
        for j in range(terms):
            sa = 0.5*svals[i]
            la = 0.5*lvals[i]
            ja = 0.5*jvals[i]
            sb = 0.5*svals[j]
            lb = 0.5*lvals[j]
            jb = 0.5*jvals[j]
            unit[i,j] /= math.sqrt((2*ja+1)*(2*jb+1))
            unit[i,j] *= SIGNEXP(sa+lb+ja+k)
            factor = wign6j(ja,jb,k,lb,la,sa)
            if not ZERO(unit[i,j]):
                unit[i,j] /= factor
    showMatrix(unit, names)



def getUnit(config, term, key, k, q="all"):
    if q != "all":
        return Matrix.Matrix(config, term, key, "U1:%i:%+i" % (k,q))
    
    unit = []
    for q in range(-k, k+1):
        unit.append(Matrix.Matrix(config, term, key, "U1:%i:%+i" % (k,q)))
    return unit



def getSpin(config, term, key, k, q="all"):
    if q != "all":
        return Matrix.Matrix(config, term, key, "T1:%i:%+i" % (k,q))
    
    unit = []
    for q in range(-k, k+1):
        unit.append(Matrix.Matrix(config, term, key, "T1:%i:%+i" % (k,q)))
    return unit



def testUnit(config, term, k):
    span = (0, 20)
    span=()

    #key = "product"
    #names = TermNames.long(term, key)
    #print names
    #quant = config["quant"]
    #unit = config["unit"]
    #matrix = CalcMatrix(quant, unit, config, MAT1_U, [2, 1])
    #showMatrix(matrix, names, span=span)
    
    key = "sl"
    names = TermNames.short(term, key)
    matrix = NielsonUnit(term, F2U2)
    #showMatrix(matrix, names)

    key = "slj"
    names = TermNames.short(term, key)
    matrix = TermReduce.expandSL(term, matrix, k, MATRIX_UNIT)
    showMatrix(matrix, names)

    key = "sljm"
    names = TermNames.short(term, key)
    matrix = TermReduce.expandSLJ(term, matrix, k, 1)
    #showMatrix(matrix, names, span=span)

    key = "product"
    names = TermNames.short(term, key)
    matrix = Transform(matrix, Numeric.transpose(term["vectors"]))
    #showMatrix(matrix, names, span=span)



    key = "product"
    names = TermNames.short(term, key)
    unit = getUnit(config, term, key, k, 1)
    #showMatrix(unit, names, span=span)

    key = "sljm"
    names = TermNames.short(term, key)
    states = term[key]["states"]
    unit = Numeric.zeros((states, states), Numeric.Float)
    for q in range(-k, k+1):
	u = getUnit(config, term, key, k, q)
	unit += reduce(term, u, k, q)
    #showMatrix(unit, names, span=span)

    key = "slj"
    names = TermNames.short(term, key)
    unit = getUnit(config, term, key, k)
    unit = TermReduce.reducetoSLJ(term, unit, k)
    showMatrix(unit, names, span=span)
    

    #print "---------------------------------------------------------------"
    #key = "slj"
    #names = TermNames.short(term, key)
    #unit = getUnit(config, term, key, k, -2)
    #showMatrix(unit, names, span=span)
    #unit = getUnit(config, term, key, k, -1)
    #showMatrix(unit, names, span=span)
    #unit = getUnit(config, term, key, k, +0)
    #showMatrix(unit, names, span=span)
    #unit = getUnit(config, term, key, k, +1)
    #showMatrix(unit, names, span=span)
    #unit = getUnit(config, term, key, k, +2)
    #showMatrix(unit, names, span=span)
    print "---------------------------------------------------------------"

    key = "slj"
    names = TermNames.short(term, key)
    unit = TermReduce.reducetoSLJ(term, getUnit(config, term, key, k), k)
    #showMatrix(unit, names, span=span)



ion = Ion("Pr3+")
ion.Term()

testUnit(ion.config, ion.term, 2)
