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
from Hamiltonian import Hamiltonian
from Matrix import Matrix
from Simplex import Simplex
from Fit import Fit
import Show

class TermFit(Fit):
    dchi2 = ""
    oldparms = ""

    # ----------------------------------------------------------------------
    def __init__(self, config, term, key, parms, meas, monitor):
        self.config = config
        self.term = term
        self.key = key
        self._splitParms_(parms)
        self.oldparms = copy.deepcopy(self.parms)
        self.meas = meas
        self.span,self.guess,self.y,self.sigma = self.splitMeas(meas)
        self.sigma = self.invSigma(self.sigma)
        self.monitor = monitor
        eps = 1.e-16
        self.eps3 = pow(eps, 1.0/3)
        
    # ----------------------------------------------------------------------
    def _getPfit_(self):
        parms = self.H.getParms()
        return Numeric.take(parms, self.pindex)
    
    # ----------------------------------------------------------------------
    def _putPfit_(self, pfit):
        parms = self.H.getParms()
        for i in range(len(self.pindex)):
            parms[self.pindex[i]] = pfit[i]
        self.H.setParms(parms)

    # ----------------------------------------------------------------------
    def _splitParms_(self, parms):
        self.parms = copy.deepcopy(parms)
        self.H = Hamiltonian(self.config, self.term, self.key, self.parms)

        pindex = []
        for i in range(len(self.parms)):
            name = self.parms[i][1]
            if (name != "offset") and (name != "base"):
                if self.parms[i][2] == 1:
                    pindex.append(i)
        self.pindex = Numeric.array(pindex)

    # ----------------------------------------------------------------------
    def _packParms_(self):
        self.oldparms = copy.deepcopy(self.parms)
        parms = self.H.getParms()
        for i in range(len(self.parms)):
            name = self.parms[i][1]
            if (name == "offset") or (name == "base"):
                offset = self.calcChi2()[1]
                parms[i] += offset
                self.parms[i][2] = 1
            self.parms[i][0] = parms[i]
        self.H.setParms(parms)

    # ----------------------------------------------------------------------
    def Chi2(self, pfit=[]):
        if pfit:
            self._putPfit_(pfit)
        return self.calcChi2()[0]

    # ----------------------------------------------------------------------
    def calcChi2(self, parms=[]):
        if parms:
            self.H.setParms(parms)
        energy,vect = self.H.diagonalize()
        names,mult  = self.H.intermediate()
        diff = self.y - self.ySum(energy, names, mult)
        offset = Numeric.add.reduce(diff) / len(diff)
        diff -= offset
        diff *= self.sigma
        chi2  = Numeric.add.reduce(diff*diff)
        return chi2, offset
        
    # ----------------------------------------------------------------------
    def _optParms_(self, factor, loop):
        pfit = self._getPfit_()
        f = factor - factor*loop/10.0
        f = max(f, 0.2)
        S = Simplex(self.Chi2, pfit, f*pfit)
        pfit,chi2,count = S.minimize(0.001, 1000, self.monitor)
        self._putPfit_(pfit)
        return chi2

    # ----------------------------------------------------------------------
    def fit(self, factor=2):
        chi2 = self.Chi2()

        loop = 0
        lastchi2 = 2*chi2
        while (chi2 > 1e-3) and (2*(lastchi2-chi2)/(lastchi2+chi2) > 1e-3):
            lastchi2 = chi2
            if self.monitor:
                print "----------------------------------------------------"
            chi2 = self._optParms_(factor, loop)
            loop += 1

        self._packParms_()
        self.chi2 = chi2
        return (self.parms, self.chi2)

    # ----------------------------------------------------------------------
    def Show(self, maxenergy=0):
        return Show.ShowLevels(self, maxenergy)
