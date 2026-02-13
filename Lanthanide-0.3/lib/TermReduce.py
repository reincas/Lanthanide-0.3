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


#-------------------------------------------------------------------------
def getSpan(term, key, name):
    term = term[key]
    states   = term["states"]
    syms     = term["symmetry"]
    termvals = term["values"]
    
    span = [(0,states)]
    level = 0
    while syms[level] != name:
        values = termvals[:,level]
        newspan = []
        for first,last in span:
            firstval = values[first]
            for i in range(first+1, last):
                if values[i] != firstval:
                    newspan.append((first, i))
                    first = i
                    firstval = values[first]
            newspan.append((first, last))
        span = newspan
        level += 1
    return span



#-------------------------------------------------------------------------
def factorSLJ(Ja, Ma, Jb, Mb, k, q):
    factor = SIGNEXP(0.5*(Ja+Ma)) * \
             wign3j(0.5*Ja,k,0.5*Jb,-0.5*Ma,q,0.5*Mb)
    return factor



#-------------------------------------------------------------------------
def factorSL(Sa, La, Ja, Sb, Lb, Jb, k, mtype):
    if mtype == MATRIX_UNIT:
        if Sa == Sb:
            factor = SIGNEXP(0.5*(Sa+Lb+Ja)+k) * \
                     math.sqrt((Ja+1)*(Jb+1)) * \
                     wign6j(0.5*Ja,0.5*Jb,k,0.5*Lb,0.5*La,0.5*Sa)
        else:
            factor = 0.0
    elif mtype == MATRIX_H4:
        if Ja == Jb:
            factor = SIGNEXP(k)/math.sqrt(2*k+1)
        else:
            factor = 0.0
    elif (mtype == MATRIX_H5) or (mtype == MATRIX_H6):
        if Ja == Jb:
            factor = SIGNEXP(0.5*(Sb+La+Ja)+k) * \
                     wign6j(0.5*La,0.5*Lb,k,0.5*Sb,0.5*Sa,0.5*Ja) / \
                     math.sqrt(2*k+1)
        else:
            factor = 0.0
    else:
        raise RuntimeError, "factorSL: Unknown matrix type!"
    return factor



#-------------------------------------------------------------------------
# Pick all matrix elements from a SLJM matrix corresponding to
# stretched states, i.e. all states, where Ma=Ja and Mb=Jb. The
# returned SLJ matrix may be converted to a reduced matrix by
# reducePreparedSLJ().
#
def prepareSLJ(term, matrix):
    key = "sljm"
    syms   = term[key]["symmetry"]
    jvals  = term[key]["values"][:,syms.index("J2")]
    mvals  = term[key]["values"][:,syms.index("Jz")]
    mspan  = getSpan(term, key, "Jz")
    terms  = len(mspan)

    index = Numeric.zeros(terms, Numeric.Int)
    for i in range(terms):
        first,last = mspan[i]
        m = jvals[first]
        index[i] = mvals[first:last].tolist().index(m) + first

    sparse = 0
    if type(matrix) == type({}):
        sparse = 1

    if sparse:
        matrix = SparseToMatrix(matrix)
    matrix = Numeric.take(matrix, index, 1)
    matrix = Numeric.take(matrix, index, 0)
    if sparse:
        matrix = MatrixToSparse(matrix)

    return matrix



#-------------------------------------------------------------------------
# Pick all matrix elements from a SLJ matrix corresponding to
# stretched states, i.e. all states, where Ja=La+Sa and Jb=Lb+Sb. The
# returned SL matrix may be converted to a reduced matrix by
# reducePreparedSL().
#
def prepareSL(term, matrix):
    key = "slj"
    syms   = term[key]["symmetry"]
    svals  = term[key]["values"][:,syms.index("S2")]
    lvals  = term[key]["values"][:,syms.index("L2")]
    jvals  = term[key]["values"][:,syms.index("J2")]
    jspan  = getSpan(term, key, "J2")
    terms  = len(jspan)

    index = Numeric.zeros(terms, Numeric.Int)
    for i in range(terms):
        first,last = jspan[i]
        j = svals[first] + lvals[first]
        index[i] = jvals[first:last].tolist().index(j) + first

    sparse = 0
    if type(matrix) == type({}):
        sparse = 1

    if sparse:
        matrix = SparseToMatrix(matrix)
    matrix = Numeric.take(matrix, index, 1)
    matrix = Numeric.take(matrix, index, 0)
    if sparse:
        matrix = MatrixToSparse(matrix)

    return matrix



#-------------------------------------------------------------------------
# Convert a SLJ matrix containing only elements corresponding to
# stretched SLJM states, as retured by prepareSLJ() to a reduced
# SLJ matrix by application of the Wigner-Eckart theorem.
#
def reducePreparedSLJ(term, matrix, k):
    key = "slj"

    if type(matrix) == type({}):
        matrix = SparseToMatrix(matrix)

    list = 0
    if type(matrix) == type([]):
        list = 1
        
    states = term[key]["states"]
    syms   = term[key]["symmetry"]
    jvals  = term[key]["values"][:,syms.index("J2")]

    newmatrix = Numeric.zeros((states, states), Numeric.Float)
    for a in range(states):
        Ja = jvals[a]
        for b in range(states):
            Jb = jvals[b]

            q = (Ja-Jb)/2
            factor = factorSLJ(Ja, Ja, Jb, Jb, k, q)
            if factor != 0.0:
                if list:
                    newmatrix[a,b] = matrix[k+q][a,b] / factor
                else:
                    newmatrix[a,b] = matrix[a,b] / factor
    return newmatrix


    
#-------------------------------------------------------------------------
# Convert a SL matrix containing only elements corresponding to
# stretched SLJ states, as retured by prepareSL() to a reduced SL
# matrix by application of the reduce equation given by mtype.
#
def reducePreparedSL(term, matrix, k, mtype):
    key = "sl"

    if type(matrix) == type({}):
        matrix = SparseToMatrix(matrix)

    list = 0
    if type(matrix) == type([]):
        list = 1
        
    states = term[key]["states"]
    syms   = term[key]["symmetry"]
    svals  = term[key]["values"][:,syms.index("S2")]
    lvals  = term[key]["values"][:,syms.index("L2")]

    newmatrix = Numeric.zeros((states, states), Numeric.Float)
    for a in range(states):
        Sa = svals[a]
        La = lvals[a]
        Ja = Sa + La
        for b in range(states):
            Sb = svals[b]
            Lb = lvals[b]
            Jb = Sb + Lb

            factor = factorSL(Sa, La, Ja, Sb, Lb, Jb, k, mtype)
            if factor != 0.0:
                newmatrix[a,b] = matrix[a,b] / factor
    return newmatrix


    
#-------------------------------------------------------------------------
# Apply the Wigner-Eckart theorem to all elements of a SLJM matrix,
# but don't actually shrink the matrix. Used by phaseSLJM().
#
def reduceSLJM(term, matrix, k, q):
    key = "sljm"
    terms  = term[key]["states"]
    syms   = term[key]["symmetry"]
    jvals  = term[key]["values"][:,syms.index("J2")]
    mvals  = term[key]["values"][:,syms.index("Jz")]

    for i in range(terms):
        for j in range(terms):
            Ja = jvals[i]
            Ma = mvals[i]
            Jb = jvals[j]
            Mb = mvals[j]
            if not ZERO(matrix[i,j]):
		factor = factorSLJ(Ja, Ma, Jb, Mb, k, q)
                if factor == 0.0:
                    print jvals
                    print mvals
                    print Ja, Ma, Jb, Mb, k, q
                    print i, j, matrix[i,j], factor # and die ...
                matrix[i,j] /= factor
    return matrix



#-------------------------------------------------------------------------
# Apply the reduce equation given by mtype to all elements of a SLJ
# matrix, but don't actually shrink the matrix. Used by phaseSLJ().
#
def reduceSLJ(term, matrix, k, mtype):
    key = "slj"
    terms  = term[key]["states"]
    syms   = term[key]["symmetry"]
    svals  = term[key]["values"][:,syms.index("S2")]
    lvals  = term[key]["values"][:,syms.index("L2")]
    jvals  = term[key]["values"][:,syms.index("J2")]
    
    for i in range(terms):
        for j in range(terms):
            Sa = svals[i]
            La = lvals[i]
            Ja = jvals[i]
            Sb = svals[j]
            Lb = lvals[j]
            Jb = jvals[j]

            if not ZERO(matrix[i,j]):
		factor = factorSL(Sa, La, Ja, Sb, Lb, Jb, k, mtype)
                matrix[i,j] /= factor
    return matrix



#-------------------------------------------------------------------------
# Correct the phases of the product/SLJM transformation vectors for
# all SLJ terms. The phases are choosen, so that all matrix elements
# of the unit tensor operator of rank k (default k=2l) for a SLJ term
# have the same sign as the diagonal elements (if not zero), which
# automatically have the right sign.
#
def phaseSLJM(config, term, k="l"):
    if k == "l":
        k = 2*config["l"]
    
    key = "sljm"
    states = term[key]["states"]
    mspan  = getSpan(term, key, "Jz")
    terms  = len(mspan)

    unit = Numeric.zeros((states, states), Numeric.Float)
    for q in range(-k, k+1):
	u = Matrix.Matrix(config, term, key, "U1:%i:%+i" % (k,q), IO_NONE)
	unit += reduceSLJM(term, u, k, q)

    index = Numeric.ones(states, Numeric.Int)
    for i in range(terms):
        first,last = mspan[i]
        sign = SIGN(unit[first,first])
        if sign == 0:
            continue
        for col in range(first, last):
            for row in range(first, col+1):
                thissign = SIGN(unit[row,col])
                if thissign:
                    if thissign*index[row] != sign:
                        index[col] = -1
                        break

    vectors = term["vectors"]
    for i in range(states):
        if index[i] < 0:
            vectors[i,:] *= -1
    return vectors



#-------------------------------------------------------------------------
# Correct the phases of the product/SLJM transformation vectors for
# all SL terms. The phases are choosen, so that all matrix elements
# of the unit tensor operator of rank k (default k=2l) for a SL term
# have the same sign as the diagonal elements (if not zero), which
# automatically have the right sign.
#
def phaseSLJ(config, term, k="l"):
    if k == "l":
        k = 2*config["l"]
    
    key = "slj"
    states = term[key]["states"]
    jspan  = getSpan(term, key, "J2")
    terms  = len(jspan)

    unit = Numeric.zeros((states, states), Numeric.Float)
    for q in range(-k, k+1):
	unit += Matrix.Matrix(config, term, key, "U1:%i:%+i" % (k,q), IO_NONE)
    unit = reducePreparedSLJ(term, unit, k)
    unit = reduceSLJ(term, unit, k, MATRIX_UNIT)

    index = Numeric.ones(states, Numeric.Int)
    for i in range(terms):
        first,last = jspan[i]
        sign = SIGN(unit[first,first])
        if sign == 0:
            continue
        for col in range(first, last):
            for row in range(first, col+1):
                thissign = SIGN(unit[row,col])
                if thissign:
                    if thissign*index[row] != sign:
                        index[col] = -1
                        break

    key = "sljm"
    mspan  = getSpan(term, key, "Jz")
    vectors = term["vectors"]
    for i in range(states):
        if index[i] < 0:
            for j in range(mspan[i][0], mspan[i][1]):
                vectors[j,:] *= -1
    return vectors



#-------------------------------------------------------------------------
# Expand a SL matrix to a SLJ matrix by application of the inverse
# reduce equation given by mtype.
#
def expandSL(term, matrix, k, mtype):
    key = "slj"
    states = term[key]["states"]
    syms   = term[key]["symmetry"]
    svals  = term[key]["values"][:,syms.index("S2")]
    lvals  = term[key]["values"][:,syms.index("L2")]
    jvals  = term[key]["values"][:,syms.index("J2")]
    jspan  = getSpan(term, key, "J2")
    terms  = len(jspan)

    if (matrix.shape[0] != terms) or (matrix.shape[1] != terms):
        raise RuntimeError, "expandSL: Misaligned matrix!"

    newmatrix = Numeric.zeros((states, states), Numeric.Float)
    for i in range(terms):
        for a in range(jspan[i][0],jspan[i][1]):
            Sa = svals[a]
            La = lvals[a]
            Ja = jvals[a]
            for j in range(terms):
                for b in range(jspan[j][0],jspan[j][1]):
                    Sb = svals[b]
                    Lb = lvals[b]
                    Jb = jvals[b]
            
                    factor = factorSL(Sa, La, Ja, Sb, Lb, Jb, k, mtype)
                    newmatrix[a,b] = matrix[i,j] * factor
    return newmatrix



#-------------------------------------------------------------------------
# Used by expandSLJ to expand a SLJ matrix to a SLJM matrix with a
# specific q by application of the Wigner-Eckart theorem.
#
def doexpandSLJ(term, matrix, k, q):
    key = "sljm"
    states = term[key]["states"]
    syms   = term[key]["symmetry"]
    jvals  = term[key]["values"][:,syms.index("J2")]
    mvals  = term[key]["values"][:,syms.index("Jz")]
    mspan  = getSpan(term, key, "Jz")
    terms  = len(mspan)

    if (matrix.shape[0] != terms) or (matrix.shape[1] != terms):
        raise RuntimeError, "expandSLJ: Misaligned matrix!"

    newmatrix = Numeric.zeros((states, states), Numeric.Float)
    for i in range(terms):
        for a in range(mspan[i][0],mspan[i][1]):
            Ja = jvals[a]
            Ma = mvals[a]
            for j in range(terms):
                for b in range(mspan[j][0],mspan[j][1]):
                    Jb = jvals[b]
                    Mb = mvals[b]

                    factor = factorSLJ(Ja, Ma, Jb, Mb, k, q)
                    newmatrix[a,b] = matrix[i,j] * factor
    return newmatrix



#-------------------------------------------------------------------------
# Expand a SLJ matrix to a SLJM matrix by application of the
# Wigner-Eckart theorem. If no q is given, than a list of the matrices
# for all values of q is returned.
#
def expandSLJ(term, matrix, k, q="all"):
    if q != "all":
        return doexpandSLJ(term, matrix, k, q)

    list = []
    for q in range(-k, k+1):
        list.append(doexpandSLJ(term, matrix, k, q))
    return list
