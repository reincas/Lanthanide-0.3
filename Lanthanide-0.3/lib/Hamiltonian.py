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
from Matrix import Matrix
import TermNames

NONE   = 0
BASE   = 1
OFFSET = 2

class Hamiltonian:
    energies = []
    eigenvectors = []
    offsettype = NONE
    
    # ----------------------------------------------------------------------
    def __init__(self, config, term, key, parms):
        self.config = config
        self.term = term
        self.key = key
        self._splitParms_(parms)
        
    # ----------------------------------------------------------------------
    def _setNone_(self):
        self.offsettype = NONE

    # ----------------------------------------------------------------------
    def _setBase_(self, index):
        self.offsettype = BASE
        self.offsetindex = index

    # ----------------------------------------------------------------------
    def _setOffset_(self, index):
        self.offsettype = OFFSET
        self.offsetindex = index

    # ----------------------------------------------------------------------
    def getParms(self):
        return self.pvalue

    # ----------------------------------------------------------------------
    def setParms(self, parms, key=""):
        if key:
            self.key = key
        if type(parms) == type([]):
            self._splitParms_(parms)
        else:
            self.pvalue = parms
        self.energies = []
        self.eigenvectors = []

    # ----------------------------------------------------------------------
    def _splitParms_(self, parms, key=""):
        if key:
            self.key = key
        self.pvalue  = []
        self.pname   = []
        self.pmatrix = []
        for i in range(len(parms)):
            value = parms[i][0]
            name = parms[i][1]
            if name == "base":
                self._setBase_(i)
                self.pmatrix.append(name)
            elif name == "offset":
                self._setOffset_(i)
                self.pmatrix.append(name)
            else:
                self.pmatrix.append(
                    Matrix(self.config, self.term, self.key, name))
            self.pvalue.append(value)
            self.pname.append(name)
        pvalue = Numeric.array(self.pvalue) * 1.0 # Keep complex!
        self.setParms(pvalue)

    # ----------------------------------------------------------------------
    def build(self):
        states = self.term[self.key]["states"]
        elements = states * (states+1) / 2
        hamiltonian = Numeric.zeros(elements, self.pvalue.typecode())
        for i in range(len(self.pmatrix)):
            matrix = self.pmatrix[i]
            if type(matrix) != type(""):
                multadd(matrix, self.pvalue[i], hamiltonian)
        return hamiltonian

    # ----------------------------------------------------------------------
    def diagonalize(self, parms=[]):
        if parms or not self.energies or not self.eigenvectors:
            if parms:
                if len(parms) != self.numparm:
                    raise RuntimeError, \
                          "Hamiltonian.diagonalize: wrong number of parms!"
                if type(parms) == Types.ListType:
                    parms = Numeric.array(parms) * 1.0 # Keep complex!
                self.pvalue = parms

            speed = DIAG_FAST
            self.hamiltonian = self.build()
            w,v,info = diagonalize(speed, DIAG_VEC, self.hamiltonian)
            self.energies = w
            self.eigenvectors = v

        if self.offsettype == BASE:
            self.energies += self.pvalue[self.offsetindex] - self.energies[0]
        elif self.offsettype == OFFSET:
            self.energies += self.pvalue[self.offsetindex]
        return self.energies, self.eigenvectors

    # ----------------------------------------------------------------------
    def intermediate(self, mspan=[], nkey="short"):
        if not self.eigenvectors:
            self.diagonalize()
            
        states = self.term[self.key]["states"]
        syms   = self.term[self.key]["symmetry"]
        if self.key != "slo":
            mult   = self.term[self.key]["values"][:,syms.index("J2")] + 1
        if nkey == "short":
            names  = TermNames.short(self.term, self.key)
        else:
            names  = TermNames.long(self.term, self.key)
        
        index  = Numeric.absolute(self.eigenvectors)
        index  = Numeric.argsort(index)
        self.termindex = index
        index  = index[:,-1]
        fnames = []
        values = Numeric.zeros(len(index), Numeric.Int)
        for i in range(len(index)):
            fnames.append(names[index[i]])
            if self.key != "slo":
                values[i] = mult[index[i]]
        if not mspan:
            return fnames, values

        mlist = Numeric.zeros(states, Numeric.Int)
        first = 0
        level = 0
        while first < states:
            last = first + values[first]+1
            mlist[first:last] = level
            level += 1
            first = last
        if first != states:
            raise RuntimeError, "Hamiltonian.intermediate: Wrong J values!"
        return fnames, values, mlist
