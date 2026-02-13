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
import Config
import Term
from Hamiltonian import Hamiltonian
from TermFit import TermFit
from JuddOfelt import JuddOfelt
from DoubleFit import DoubleFit
from OscPlot import OscPlot
from SigmaPlot import SigmaPlot
from BaseParms import BaseParms
import Show
import TermNames


CONFIG = {
    "Cr3+": (2,  3),
    "Ce3+": (3,  1),
    "Pr3+": (3,  2),
    "Nd3+": (3,  3),
    "Pm3+": (3,  4),
    "Sm3+": (3,  5),
    "Eu3+": (3,  6),
    "Gd3+": (3,  7),
    "Tb3+": (3,  8),
    "Dy3+": (3,  9),
    "Ho3+": (3, 10),
    "Er3+": (3, 11),
    "Tm3+": (3, 12),
    "Yb3+": (3, 13),
    "d01":  (2,  1),
    "d02":  (2,  2),
    "d03":  (2,  3),
    "d04":  (2,  4),
    "d05":  (2,  5),
    "d06":  (2,  6),
    "d07":  (2,  7),
    "d08":  (2,  8),
    "d09":  (2,  9),
    "f01":  (3,  1),
    "f02":  (3,  2),
    "f03":  (3,  3),
    "f04":  (3,  4),
    "f05":  (3,  5),
    "f06":  (3,  6),
    "f07":  (3,  7),
    "f08":  (3,  8),
    "f09":  (3,  9),
    "f10":  (3, 10),
    "f11":  (3, 11),
    "f12":  (3, 12),
    "f13":  (3, 13),
    }


class Ion:
    ion = ""
    l = -1
    electrons = -1
    config = {}
    term = {}
    parms = []
    cauchy = []
    maxenergy = 0
    H  = ""
    TF = ""
    JO = ""
    DF = ""
    listMeas = []
    listCalc = []
    filemode = "w"
    bandspan = []
    
    # ----------------------------------------------------------------------
    def __init__(self, ion, monitor=1, modified=0):
        self.monitor = monitor
        self.ion = ion
        l,electrons = CONFIG[ion]
	self.l = l
	self.electrons = electrons
	self.config = Config.Config(l, electrons)
	self.term = Term.Term(self.config)
        self.modified = modified

    # ----------------------------------------------------------------------
    def confLongNames(self, key=""):
        if not key:
            key = self.key
        return TermNames.full(self.term, key)
            
    def confShortNames(self, key=""):
        if not key:
            key = self.key
        return TermNames.short(self.term, key)
            
    def setMaxenergy(self, maxenergy):
        self.maxenergy = maxenergy
        
    def setBandspan(self, bandspan):
        self.bandspan = bandspan
        
    # ----------------------------------------------------------------------
    def termParms(self, parms, key="slj"):
	self.parms = parms
        self.key = key
        if self.H:
            self.H.setParms(parms, key)
        else:
            self.H = Hamiltonian(self.config, self.term, key, parms)

    def termNew(self, meas=[]):
        self.TF = TermFit(self.config, self.term, self.key,
                          self.parms, meas, self.monitor)

    def termShow(self, fn, meas=0):
        if meas:
            self.termNew(meas)
        self.termShowData(fn, "Energy level data for %s:" % self.ion)

    def termFit(self, fn, meas=0):
        if meas:
            self.termNew(meas)
        parms,chi2 = self.TF.fit()
	self.termParms(parms)
        sigma = self.termShowData(fn, "Energy level fit of %s:" % self.ion)
        return parms, sigma

    def termShowData(self, fn, title):
	show = self.TF.Show(self.maxenergy)
        show.open(fn, self.filemode)
        show.showTitle(title)
        show.showString()
        show.showParms(1)
        show.showString()
        show.showResult()
        show.showString()
        show.showVectors()
        sigma = show.stats[2]
        show.close()
        return sigma

    # ----------------------------------------------------------------------
    def oscCauchy(self, cauchy):
        self.cauchy = cauchy
        
    def oscNew(self, fmeas):
        self.JO = JuddOfelt(self.config, self.term, self.parms,
                            self.cauchy, fmeas, self.bandspan, self.monitor,
                            modified=self.modified)

    def oscShow(self, fn, fmeas=0):
        if fmeas:
            self.oscNew(fmeas)
        self.oscShowData(fn, "Judd-Ofelt data for %s:" % self.ion)

    def oscFit(self, fn, fmeas=0):
        if fmeas:
            self.oscNew(fmeas)
        omega,chi2 = self.JO.fit(self.parms)
        self.names = self.JO.names
        self.JO.calcData()
        self.tau = self.JO.tau
        self.Aed = self.JO.Aed
        self.Amd = self.JO.Amd
        sigma = self.oscShowData(fn, "Judd-Ofelt fit of %s:" % self.ion)
        return omega, sigma

    def oscShowData(self, fn, title):
        show = self.JO.Show(self.maxenergy)
        show.open(fn, self.filemode)
        show.showTitle(title)
        show.showString()
        show.showOmega()
        show.showString()
        show.showResult()
        show.showString()
        show.showJOdata()
        sigma = show.stats[2]
        show.close()
        return sigma
    
    # ----------------------------------------------------------------------
    def doubleLaTeX(self, fn, meas=0, fmeas=0, span=0):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
	show = self.DF.LaTeX(self.maxenergy)
        show.open(fn, self.filemode)
        show.showResult(span)
        show.close()

    def vectorLaTeX(self, fn, meas=0, fmeas=0, span=()):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
	show = self.DF.LaTeX(self.maxenergy)
        show.open(fn, self.filemode)
        show.showVectorsLaTeX(span)
        show.close()

    def emissionLaTeX(self, fn, meas=0, fmeas=0, span=0):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
	show = self.DF.LaTeX(self.maxenergy)
        show.open(fn, self.filemode)
        show.showJOdata("e", span)
        show.close()

    def absorptionLaTeX(self, fn, meas=0, fmeas=0, span=0):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
	show = self.DF.LaTeX(self.maxenergy)
        show.open(fn, self.filemode)
        show.showJOdata("a", span)
        show.close()

    # ----------------------------------------------------------------------
    def doubleNew(self, meas, fmeas):
        self.DF = DoubleFit(self.config, self.term, self.key,
                            self.parms, meas, self.cauchy, fmeas,
                            self.bandspan, self.monitor)

    def doubleShow(self, fn, meas=0, fmeas=0):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
        self.doubleShowData(fn, "Double data for %s:" % self.ion)

    def doubleFit(self, fn, meas=0, fmeas=0):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
	parms,omega,chi2 = self.DF.fit()
	self.termParms(parms)
        self.TF = self.DF.TF
        self.JO = self.DF.JO
        self.JO.calcData()
        self.names = self.JO.names
        self.tau = self.JO.tau
        self.Aed = self.JO.Aed
        self.Amd = self.JO.Amd
        sigma = self.doubleShowData(fn, "Double fit of %s:" % self.ion)
        return parms, sigma[0], omega, sigma[1]

    def doubleShowData(self, fn, title):
	show = self.DF.Show(self.maxenergy)
        show.open(fn, self.filemode)
        show.showTitle(title)
        show.showString()
	show.showParms(1)
        show.showString()
        show.showOmega()
        show.showString()
        show.showResult()
        show.showString()
        show.showVectors()
        show.showString()
        show.showJOdata()
        sigma = show.stats[2]
        show.close()
        return sigma
        
    # ----------------------------------------------------------------------
    def fullFit(self, fn, meas, fmeas, double=1, method=0):
        if double:
            parms, dk, omega, df = self.doubleFit(fn, meas, fmeas)
        else:
            sigma = []
            parms, dk = self.termFit(fn, meas)
            self.filemode = "a"
            show = self.TF.Show(self.maxenergy)
            show.open(fn, self.filemode)
            show.showString()
            show.close()
            omega, df = self.oscFit(fn, fmeas)
            self.filemode = "w"
        return parms, dk, omega, df

    # ----------------------------------------------------------------------
    def plotOsc(self, fn, xmax, ymax, meas=0, fmeas=0):
        if meas or fmeas:
            self.doubleNew(meas, fmeas)
	show = self.DF.Show(self.maxenergy)
        listCalc,listMeas = show.plotLists()
        plot = OscPlot(listMeas, listCalc)
        plot.plotData(fn, xmax, ymax)

    # ----------------------------------------------------------------------
    def plotSigma(self, fn, xmax, element, base, sigma, meas=0):
        self.termNew(meas)
        energy,v = self.TF.H.diagonalize()
        names,mult = self.TF.H.intermediate()
        plot = SigmaPlot(sigma, energy, names, element, base, xmax)
        plot.plotData(fn)

    # ----------------------------------------------------------------------
    def getLevels(self, parms=[], key="slj"):
        if type(parms) == type(""):
            parms = BaseParms[parms]
            self.termParms(parms, key)            
        elif parms:
            self.termParms(parms, key)
        self.termNew()
        energy,v = self.TF.H.diagonalize()
        names,mult = self.TF.H.intermediate()
        return energy, names

