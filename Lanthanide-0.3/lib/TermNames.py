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

#-------------------------------------------------------------------------
casimir = {}

def casimirRk(w):
    sum = 0
    for i in range(len(w)):
        sum += w[i] * (w[i] - 1 + 2*len(w) - 2*i)
    return sum / 2

def casimirSU(l):
    n = 0
    for i in range(len(l)):
        n += l[i]
    sum = 0
    for i in range(len(l)):
        sum += l[i] * (l[i] + 1 + len(w) - 2*i)
    sum = len(w)*sum - n*n
    return sum

def casimirG2(u):
    return u[0]*u[0] + u[1]*u[1] + u[0]*u[1] + 5*u[0] + 4*u[1]

def casimirDict():
    if casimir == {}:
        dict = {}
        for i in range(3):
            for j in range(i+1):
                w = (i,j)
                val = casimirRk(w)
                dict[val] = w
        casimir["R5"] = dict
        dict = {}
        for i in range(3):
            for j in range(i+1):
                for k in range(j+1):
                    w = (i,j,k)
                    val = casimirRk(w)
                    dict[val] = w
        casimir["R7"] = dict
        dict = {}
        for i in range(5):
            for j in range(i+1):
                u = (i,j)
                val = casimirG2(u)
                dict[val] = u
        casimir["G2"] = dict
    return



#-------------------------------------------------------------------------
def mlName(ml):
    ml /= 2
    return "%+i" % ml

def msName(ms):
    if ms < 0:
        return "m"
    return "p"

def sName(s):
    return "%i" % (s+1)

def lName(l):
    l /= 2
    spd = "SPDFGHIKLMNOPQRSTUVWXYZ"
    return spd[l]

def jName(j):
    if (j % 2) == 0:
        return "%i" % (j/2)
    return "%i/2" % j

def jzName(jz):
    if (jz % 2) == 0:
        return "[%+i]" % (jz/2)
    return "[%+i/2]" % jz

def tauName(tau):
    tau = " AB"[tau]
    if tau == " ":
        tau = ""
    return tau

def numName(num):
    if num:
        num = "(%i)" % num
    else:
        num = ""
    return num

def wuName(x):
    s = ""
    for val in x:
        s += "%i" % val
    return "(%s)" % s

def nuName(w, s):
    nu = 0
    for i in range(len(w)):
        if w[i] == 2:
            nu += 1
    return "%i" % (2*nu + s)



#-------------------------------------------------------------------------
def full(term, key):
    return long(term, key)



#-------------------------------------------------------------------------
def long(term, key):
    term = term[key]
    termvals = term["values"]
    if key == "product":
        states,electrons = termvals.shape
        electrons /= 3
        names = []
        for i in range(states):
            name = ""
            for j in range(electrons):
                ml = mlName(termvals[i,3*j+1])
                ms = msName(termvals[i,3*j+2])
                name += "," + ml + ms
            names.append(name[1:])
        return names

    casimirDict()
    syms = term["symmetry"]
    Nu   = ""
    Tau  = ""
    W    = ""
    U    = ""
    S    = ""
    L    = ""
    Num  = ""
    J    = ""
    Jz   = ""
    D    = ""
    
    names = []
    for term in range(term["states"]):
        vals = termvals[term,:]
        for i in range(len(syms)):
            sym = syms[i]
            if sym == "S2":
                s = vals[i]
                S = sName(s)
            elif sym == "L2":
                L = lName(vals[i])
            elif sym == "J2":
                J = jName(vals[i])
            elif sym == "Jz":
                Jz = jzName(vals[i])
            elif sym == "R5":
                W = wuName(casimir["R5"][vals[i]])
            elif sym == "R7":
                w = casimir["R7"][vals[i]]
                W = wuName(w)
            elif sym == "G2":
                U = wuName(casimir["G2"][vals[i]])
            elif sym == "tau":
                Tau = tauName(vals[i])
            elif sym == "num":
                Num = numName(vals[i])
            else:
                D += "[%s]" % vals[i]
                #raise RuntimeError, "termNames: Unknown symmetry name!"
        if U:
            Nu = nuName(w, s)
        name = Nu + Tau + W + U + S + L + Num + J + Jz + D
        Nu   = ""
        Tau  = ""
        W    = ""
        U    = ""
        S    = ""
        L    = ""
        Num  = ""
        J    = ""
        Jz   = ""
        D    = ""
        names.append(name)
    return names



#-------------------------------------------------------------------------
def short(term, key):
    term = term[key]
    termvals = term["values"]
    if key == "product":
        states,electrons = termvals.shape
        electrons /= 3
        names = []
        for i in range(states):
            name = ""
            for j in range(electrons):
                state = termvals[i,3*j]
                name += ",%i" % state
            names.append(name[1:])
        return names

    syms = term["symmetry"]
    J    = ""
    Jz   = ""

    names = []
    for term in range(term["states"]):
        vals = termvals[term,:]
        for i in range(len(syms)):
            sym = syms[i]
            if sym == "S2":
                S = sName(vals[i])
            elif sym == "L2":
                L = lName(vals[i])
            elif sym == "J2":
                J = jName(vals[i])
            elif sym == "Jz":
                Jz = jzName(vals[i])
            elif sym == "num":
                Num = numName(vals[i])
        name = S + L + Num + J + Jz
        names.append(name)
    return names
