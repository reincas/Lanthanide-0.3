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
from Simplex import Simplex
from TermFit import TermFit
from JuddOfelt import JuddOfelt
from Show import ShowDouble, ShowDoubleLaTeX

class DoubleFit:
    # ----------------------------------------------------------------------
    def __init__(self, config, term, key, parms, energies, cauchy, fmeas,
                 bandspan=[], monitor=1):
        self.bandspan = bandspan
        self.monitor = monitor
        self.TF = TermFit(config, term, key, parms, energies, monitor)
        self.JO = JuddOfelt(config, term, parms, cauchy, fmeas, [], monitor)

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

        self.TF._packParms_()
        self.TF.chi2 = self.TF.Chi2()
        self.chi2 = chi2

        self.JO.fed = (self.JO.u22*self.JO.omega[0] + \
                       self.JO.u42*self.JO.omega[1] + \
                       self.JO.u62*self.JO.omega[2]) * self.JO.Fed
        return (self.TF.parms, self.JO.omega, chi2)

    # ----------------------------------------------------------------------
    def Chi2(self, pfit=[]):
        tfchi2 = self.TF.Chi2(pfit)
        parms = self.TF.H.getParms()
        omega,jochi2 = self.JO.fit(parms)
        return tfchi2 + jochi2

    # ----------------------------------------------------------------------
    def _optParms_(self, factor, loop):
        pfit = self.TF._getPfit_()
        f = factor + (0.2-factor)*loop/10.0
        f = max(f, 0.2)
        S = Simplex(self.Chi2, pfit, f*pfit)
        pfit,chi2,count = S.minimize(0.001, 1000, self.monitor)
        if self.bandspan:
            self.JO.bandspan = self.bandspan
            chi2 = self.Chi2()
            self.JO.bandspan = []
        self.TF._putPfit_(pfit)
        return chi2

    # ----------------------------------------------------------------------
    def Show(self, maxenergy=0):
        return ShowDouble(self.TF, self.JO, maxenergy)
    
    # ----------------------------------------------------------------------
    def LaTeX(self, maxenergy=0):
        return ShowDoubleLaTeX(self.TF, self.JO, maxenergy)
    
