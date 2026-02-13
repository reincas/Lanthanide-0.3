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
import Matrix
import TermReduce



#-------------------------------------------------------------------------
symmetry = {
    1: [ "S2", "L2", "num", "J2", "Jz" ],
    2: [ "S2", "R5", "L2", "num", "J2", "Jz" ],
    3: [ "S2", "R7", "G2", "L2", "tau", "num", "J2", "Jz" ],
    }



#-------------------------------------------------------------------------
def quantProduct(config):
    electrons = config["electrons"]
    states = config["states"]
    values = Numeric.zeros((states,3*electrons), Numeric.Int)
    ml = config["quant"]["ml"]
    ms = config["quant"]["ms"]
    for i in range(states):
        state = config["state"][i,:]
        for j in range(electrons):
            values[i,3*j] = state[j]
            values[i,3*j+1] = INT(2*ml[state[j]])
            values[i,3*j+2] = INT(2*ms[state[j]])
    syms = []
    for j in range(electrons):
        syms.append("state(%i)" % j)
        syms.append("ml(%i)" % j)
        syms.append("ms(%i)" % j)
    return syms, values



#-------------------------------------------------------------------------
def termProduct(config):
    print "termProduct: prepare"
    syms,values = quantProduct(config)

    print "termProduct: build term dictionary"
    term = {}
    term["key"]      = TERM_PRODUCT
    term["states"]   = config["states"]
    term["symmetry"] = syms
    term["values"]   = values
    return term



#-------------------------------------------------------------------------
def rescaleQuant(sym, vals, config):
    if (sym == "S2") or (sym == "L2") or (sym == "J2"):
        vals = (Numeric.sqrt(1.0+4.0*vals)-1.0)
    elif sym == "Jz":
        vals *= 2
    elif sym == "R5":
        vals *= 3
    elif sym == "R7":
        vals *= 5
    elif sym == "G2":
        vals *= 12
    elif sym == "OPE":
        #for i in range(len(vals)):
        #    vals[i] = 2 * vals[i]**2
        vals *= 10000
    elif sym[0:5] == "Test0":
        vals *= 10000
    elif sym[0:5] == "Test1":
        vals *= 10000
    #elif sym == "T0:4:0":
    #    vals += 12
    #    vals /= 5.0
    #    vals *= 10000
    elif sym == "Dq":
        if config["l"] == 2:
            num = config["electrons"]
            vals = 0.1 * (vals + 4*num)
        vals *= 10000

    for i in range(len(vals)):
        vals[i] += SIGN(vals[i]) * 0.1
    return vals.astype(Numeric.Int)



#-------------------------------------------------------------------------
def termDiag(config, term, key, sym, span, vectors):
    states = config["states"]
        
    print "termDiag: load symmetry matrix"
    matrix = Matrix.Matrix(config, term, key, sym)

    if not span:
        print "termDiag: diagonalize symmetry matrix"
        matrix = SparseToTri(matrix)
        values,vectors,info = diagonalize(DIAG_FAST, DIAG_VEC, matrix)
        span = [(0, states)]
    else:
        print "termDiag: transform symmetry matrix"
        matrix = transform(matrix, vectors)
        matrix = SparseToTri(matrix)

        print "termDiag: diagonalize all submatrices"
        values = Numeric.zeros(states, Numeric.Float)
        for i,j in span:
            if j-i == 1:
                values[i] = submatrix(matrix, i, i+1)[0]
            else:
                sub = submatrix(matrix, i, j)
                w,v,info = diagonalize(DIAG_FAST, DIAG_VEC, sub)
                vectors[i:j,:] = Numeric.matrixmultiply(v, vectors[i:j,:])
                values[i:j] = w
    values = rescaleQuant(sym, values, config)
    
    print "termDiag: search new submatrices"
    newspan = []
    for i,j in span:
        first = i
        for k in range(i+1, j):
            if values[k] != values[first]:
                newspan.append((first, k))
                first = k
        newspan.append((first, j))

    return values, vectors, newspan



#-------------------------------------------------------------------------
def termTau(config, syms, termvals, span):
    states = config["states"]
    jvals = termvals[:,syms.index("J2")]
    num = 0.0

    values = Numeric.zeros(states, Numeric.Int)
    for first,last in span:
        if last - first > 1:
            tau = 0
            for i in range(first, last):
                tau += 1
                values[i] = tau
            num += (tau-1) / (1.0+jvals[first])

    diff = abs(num - INT(num))
    if diff > 1e-6:
        raise RuntimeError, \
              "tauSplit: Noninteger number of tau terms (%g, %g)!" % \
              (num, diff)
    print "tauSplit: found %i tau terms" % num

    termvals[:,syms.index("tau")] = values
    return termvals



#-------------------------------------------------------------------------
def sortValues(syms, values, vectors, name, span=[]):
    if not span:
        span = [ (0,values.shape[0]) ]

    vals = values[:,syms.index(name)]
    index = Numeric.arrayrange(len(vals))
    for first,last in span:
        index[first:last] = Numeric.argsort(vals[first:last]) + first

    newspan = []
    for first,last in span:
        firstval = vals[index[first]]
        for i in range(first+1, last):
            if vals[index[i]] != firstval:
                newspan.append((first, i))
                first = i
                firstval = vals[index[first]]
        newspan.append((first, last))

    for first,last in newspan:
        if last-first > 1:
            index[first:last] = Numeric.sort(index[first:last])
                
    values = Numeric.take(values, index)
    vectors = Numeric.take(vectors, index)
    return values, vectors, newspan



#-------------------------------------------------------------------------
def termNumbers(config, syms, termvals, lspan, nspan):
    states = config["states"]
    fullspan = []
    i = 0
    for first,last in lspan:
        span = []
        while (i < len(nspan)) and (nspan[i][1] <= last):
            span.append(nspan[i])
            i += 1
        fullspan.append(span)
    values = Numeric.zeros(states, Numeric.Int)
    for span in fullspan:
        if len(span) > 1:
            num = 0
            for first,last in span:
                num += 1
                values[first:last] = num
    termvals[:,syms.index("num")] = values
    return termvals
                


#-------------------------------------------------------------------------
def termSLJM(config, term):
    key = "product"
    if not term.has_key(key):
        raise RuntimeError, "termSLJM: Need a product term to work on!"

    print "termSLJM: prepare"
    states = config["states"]
    syms = symmetry[config["l"]]
    numsyms = len(syms)

    print "termSLJM: test l"
    if config["l"] not in (1, 2, 3):
        raise RunTimeError, "termSLJM: Rules available only for l=1,2,3!"

    print "termSLJM: main loop"
    termvals = Numeric.zeros((states,len(syms)), Numeric.Int)
    vectors = []
    span = []
    for level in range(numsyms):
        sym = syms[level]
        if (sym != "tau") and (sym != "num"):
            values,vectors,span = termDiag(config, term, key, sym,
                                           span, vectors)
            termvals[:,level] = values
            print "termSLJM: %s completed, %i terms" % (sym, len(span))

    if syms.count("tau") != 0:
        print "termSLJM: split to different tau values"
        termvals = termTau(config, syms, termvals, span)

    print "termSLJM: reorder symmetries"
    newsyms = [ "S2", "L2" ]
    for sym in syms:
        if (sym != "S2") and (sym != "L2"):
            newsyms.append(sym)
    index = Numeric.zeros(numsyms, Numeric.Int)
    for i in range(numsyms):
        index[i] = syms.index(newsyms[i])
    syms = newsyms
    termvals = Numeric.take(termvals, index, 1)

    print "termSLJM: sort states"
    span = []
    for sym in syms:
        if sym == "J2":
            nspan = span
        termvals,vectors,span = sortValues(syms, termvals, vectors, sym, span)
        if sym == "L2":
            lspan = span
    print "termSLJM: %i SL terms" % len(nspan)

    print "termSLJM: search for unresolved states"
    sum = 0
    for first,last in span:
        sum += last - first - 1
    if sum != 0:
        raise RuntimeError, "termSLJM: %i unresolved states!" % sum

    print "termSLJM: Calculate term numbers"
    termvals = termNumbers(config, syms, termvals, lspan, nspan)

    print "termSLJM: build term dictionary"
    term = {}
    term["key"]      = TERM_SLJM
    term["states"]   = states
    term["symmetry"] = syms
    term["values"]   = termvals
    return term, vectors



#-------------------------------------------------------------------------
def termSOh(config, term, syms):
    key = "product"
    if not term.has_key(key):
        raise RuntimeError, "termSOh: Need a product term to work on!"

    print "termSOh: prepare"
    states = config["states"]
    #syms = symmetry[config["l"]]
    #syms = [ "S2", "R5", "L2", "num", "J2", "Jz" ]
    #syms = [ "S2", "C4" ]
    numsyms = len(syms)

    print "termSOh: test l"
    if config["l"] not in (1, 2, 3):
        raise RunTimeError, "termSOh: Rules available only for l=1,2,3!"

    print "termSOh: main loop"
    termvals = Numeric.zeros((states,len(syms)), Numeric.Int)
    vectors = []
    span = []
    for level in range(numsyms):
        sym = syms[level]
        if (sym != "tau") and (sym != "num"):
            values,vectors,span = termDiag(config, term, key, sym,
                                           span, vectors)
            termvals[:,level] = values
            print "termSOh: %s completed, %i terms" % (sym, len(span))

    if syms.count("tau") != 0:
        print "termSOh: split to different tau values"
        termvals = termTau(config, syms, termvals, span)

    #print "termSOh: reorder symmetries"
    #newsyms = [ "S2", "L2" ]
    #for sym in syms:
    #    if (sym != "S2") and (sym != "L2"):
    #        newsyms.append(sym)
    #index = Numeric.zeros(numsyms, Numeric.Int)
    #for i in range(numsyms):
    #    index[i] = syms.index(newsyms[i])
    #syms = newsyms
    #termvals = Numeric.take(termvals, index, 1)

    #print "termSOh: sort states"
    #span = []
    #for sym in syms:
    #    if sym == "J2":
    #        nspan = span
    #    termvals,vectors,span = sortValues(syms, termvals, vectors, sym, span)
    #    if sym == "L2":
    #        lspan = span
    #print "termSOh: %i SL terms" % len(nspan)

    #print "termSOh: search for unresolved states"
    #sum = 0
    #for first,last in span:
    #    sum += last - first - 1
    #if sum != 0:
    #    raise RuntimeError, "termSOh: %i unresolved states!" % sum

    #print "termSOh: Calculate term numbers"
    #termvals = termNumbers(config, syms, termvals, lspan, nspan)

    print "termSOh: build term dictionary"
    term = {}
    term["key"]      = TERM_SLJM
    term["states"]   = states
    term["symmetry"] = syms
    term["values"]   = termvals
    return term, vectors



#-------------------------------------------------------------------------
def termSLJ(term):
    key = "sljm"
    if not term.has_key(key):
        raise RuntimeError, "termSLJ: Need a SLJM term to work on!"

    print "termSLJ: prepare"
    syms   = copy.copy(term[key]["symmetry"])
    values = term[key]["values"]
    jvals  = values[:,syms.index("J2")]
    mvals  = values[:,syms.index("Jz")]
    mspan  = TermReduce.getSpan(term, key, "Jz")
    terms  = len(mspan)

    print "termSLJ: pick %i stretched states" % terms
    index = Numeric.zeros(terms, Numeric.Int)
    for i in range(terms):
        first,last = mspan[i]
        for j in range(first, last):
            if mvals[j] == jvals[j]:
                index[i] = j
                break
        else:
            raise RuntimeError, "termSLJ: m=j state missing!"
    termvals = Numeric.take(values, index)

    print "termSLJ: delete Jz data"
    jzpos = syms.index("Jz")
    index = range(len(syms))
    del index[jzpos]
    index = Numeric.array(index)
    termvals = Numeric.take(termvals, index, 1)
    del syms[jzpos]
    
    print "termSLJ: build term dictionary"
    term = {}
    term["key"]      = TERM_SLJ
    term["states"]   = terms
    term["symmetry"] = syms
    term["values"]   = termvals
    return term



#-------------------------------------------------------------------------
def termSL(term):
    key = "slj"
    if not term.has_key(key):
        raise RuntimeError, "termSL: Need a SLJ term to work on!"

    print "termSL: prepare"
    syms   = copy.copy(term[key]["symmetry"])
    values = term[key]["values"]
    svals  = values[:,syms.index("S2")]
    lvals  = values[:,syms.index("L2")]
    jvals  = values[:,syms.index("J2")]
    nspan  = TermReduce.getSpan(term, key, "num")
    terms  = len(nspan)

    print "termSL: pick %i stretched states" % terms
    index = Numeric.zeros(terms, Numeric.Int)
    for i in range(terms):
        first,last = nspan[i]
        for j in range(first, last):
            if jvals[j] == svals[j] + lvals[j]:
                index[i] = j
                break
        else:
            raise RuntimeError, "termSL: j=s+l state missing!"
    termvals = Numeric.take(values, index)

    print "termSL: delete J data"
    jpos = syms.index("J2")
    index = range(len(syms))
    del index[jpos]
    index = Numeric.array(index)
    termvals = Numeric.take(termvals, index, 1)
    del syms[jpos]
    
    print "termSL: build term dictionary"
    term = {}
    term["key"]      = TERM_SL
    term["states"]   = terms
    term["symmetry"] = syms
    term["values"]   = termvals
    return term



#-------------------------------------------------------------------------
def buildTerm(config):
    l = config["l"]
    term = {}
    term["product"] = termProduct(config)
    term["sljm"],term["vectors"] = termSLJM(config, term)
    for k in range(0, 2*l+1, 2):
        term["vectors"] = TermReduce.phaseSLJM(config, term, k)
    term["slj"] = termSLJ(term)
    #for k in range(0, 2*l+1, 2):
    #    term["vectors"] = TermReduce.phaseSLJ(config, term, k)
    term["sl"] = termSL(term)
    #term["slo"],v = termSOh(config, term, [ "S2", "OPE" ])
    #term["slovectors"] = v
    return term



#-------------------------------------------------------------------------
def Term(config, rw=IO_READ+IO_WRITE):
    fn = dataPath(config, "", "term")
    if (rw & IO_READ) and os.path.exists(fn):
        fp = gzip.open(fn, "r")
        term = cPickle.loads(fp.read())
        fp.close()
    else:
        print "Term: calc %s" % fn
        term = buildTerm(config)
        if rw & IO_WRITE:
            fp = gzip.open(fn, "w")
            fp.write(cPickle.dumps(term))
            fp.close()
    return term
