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
import TermNames
from Matrix import Matrix
import TermReduce
from Hamiltonian import Hamiltonian
from Fit import Fit
import Show


FED = 8 * math.pi*math.pi * CONST_me * CONST_c * 100 / (3 * CONST_h)
FMD = CONST_h * 100 / (6 * CONST_me * CONST_c)
AED = 64 * pow(math.pi,4) * 1e6 / (3 * CONST_h)
AMD = 4 * math.pi*math.pi * CONST_h * 1e6 / (3 * pow(CONST_me*CONST_c,2))
COULOMB = CONST_e*CONST_e / (4 * math.pi * CONST_eps0)


class JuddOfelt(Fit):
    energy = ""
    
    # ----------------------------------------------------------------------
    def __init__(self, config, term, parms, cauchy, meas,
                 bandspan=[], monitor=1, otype=1, modified=0):
        self.config = config
        self.term = term
        self.key = "slj"
        self.parms = parms
        self.H = Hamiltonian(config, term, self.key, parms)
        index = [ 2, 4, 6 ]
        if modified:
            index = [ 1, 2, 3, 4, 5, 6 ]
        self.num = len(index)
        states = term[self.key]["states"]
        self.u = Numeric.zeros((self.num, states, states), Numeric.Float)
        self.u2 = Numeric.zeros((self.num, states), Numeric.Float)
        for i in range(self.num):
            self.u[i,:,:] = self._unitReduced_(index[i])
        #self.u2 = self._unitReduced_(2)
        #self.u4 = self._unitReduced_(4)
        #self.u6 = self._unitReduced_(6)
        self.ls = self._magnReduced_()
        self.cauchy = cauchy
        self.meas = meas
        self.span,self.guess,self.fmeas,self.sigma = self.splitMeas(meas)
        self.sigma = self.invSigma(self.sigma)
        self.otype = otype
        eps = 1.e-16
        self.eps3 = pow(eps, 1.0/3)
        self.bandspan = bandspan
        self.monitor = monitor
        self.fit()

        #self.names = []
        #self.tau = []
        #self.Aed = []
        #self.Amd = []

    # ----------------------------------------------------------------------
    def _unitReduced_(self, k):
        unit = self._sumUnit_("U1:%i:%+i", k)
        unit = TermReduce.reducePreparedSLJ(self.term, unit, k)
        return unit

    # ----------------------------------------------------------------------
    def _magnReduced_(self):
        l = self.config["l"]
        lq = self._sumUnit_("U1:%i:%+i", 1) * math.sqrt(l*(l+1)*(2*l+1))
        sq = self._sumUnit_("T1:%i:%+i", 1) * math.sqrt(1.5)
        sum = lq+2*sq
        sum = TermReduce.reducePreparedSLJ(self.term, sum, 1)
        return sum

    # ----------------------------------------------------------------------
    def _sumUnit_(self, name, k):
        states = self.term[self.key]["states"]
        sum = Numeric.zeros((states, states), Numeric.Float)
        for q in range(-k, k+1):
            sum += Matrix(self.config, self.term, self.key, name % (k,q))
        return sum

    # ----------------------------------------------------------------------
    def fit(self, parms=""):
        self._optParms_(parms)
        
        if not self.bandspan:
            return self.omega, self.chi2
        
        chi2 = self.chi2
        lastchi2 = 2*chi2
        loop = 0
        while (chi2 > 1e-3) and (2*(lastchi2-chi2)/(lastchi2+chi2) > 1e-3):
            if self.monitor:
                print "jofit: iteration = %i,  chi2 = %f" % (loop, chi2)
            lastchi2 = chi2
            self.adjust()
            chi2 = self._optParms_(parms)[1]
            loop += 1
        if self.monitor:
            print "jofit: chi2 = %f" % chi2
        self.meas = self.updateMeas(self.meas, self.fmeas)
        return self.omega, self.chi2

    # ----------------------------------------------------------------------
    # obsolet (and valid only for 3-parameter JO anyway)
    #def makeParms(self):
    #    parms = []
    #    for i in range(3):
    #        name = "omega%i" % (2*i+2)
    #        parms.append([omega[i], name, 1])
    #    return parms

    # ----------------------------------------------------------------------
    def _optParms_(self, parms=""):
        self.A,self.b = self._calcMatrix_(parms)
    
        if self.otype == 1:
            for i in range(self.num):
                self.A[:,i] *= self.sigma
            self.b *= self.sigma

        #raise "Debug", str(A.shape)
        self.omega = self._solve_(self.A, self.b)
        u2 = Numeric.transpose(self.u2)
        self.fed = Numeric.matrixmultiply(u2, self.omega) * self.Fed
        #self.fed = (self.u22*self.omega[0] + \
        #            self.u42*self.omega[1] + \
        #            self.u62*self.omega[2]) * self.Fed
    
        self.chi2 = self.calcChi2(self.omega)
        return self.omega, self.chi2
    
    # ----------------------------------------------------------------------
    def _calcMatrix_(self, parms=""):
        if parms:
            self.H.setParms(parms)
        self.energy,eigenvectors = self.H.diagonalize()
        self.names,mult = self.H.intermediate()
        j = 0.5 * (mult[0]-1)

        dk = self.energy - self.energy[0]
        n = self.calcCauchy(dk)
        n[0] = 1
        Xed = n*n + 2
        Xed = Xed*Xed / (9*n)
        self.Fed = Xed * FED * dk / (2*j+1)
        self.Fed[0] = 0
        Xmd = n
        self.Fmd = Xmd * FMD * dk / (2*j+1)
        self.Fmd[0] = 0

        for i in range(self.num):
            self.u2[i,:]=self._squareMatrix_(self.u[i,:,:], eigenvectors)[0,:]
        #self.u22 = self._squareMatrix_(self.u2, eigenvectors)[0,:]
        #self.u42 = self._squareMatrix_(self.u4, eigenvectors)[0,:]
        #self.u62 = self._squareMatrix_(self.u6, eigenvectors)[0,:]
        self.ls2 = self._squareMatrix_(self.ls, eigenvectors)[0,:]
        
        A = self.u2 * self.Fed
        #A = Numeric.array([self.u22*self.Fed,
        #                   self.u42*self.Fed,
        #                   self.u62*self.Fed], Numeric.Float)
        A = Numeric.transpose(A)
        self.fmd = self.ls2*self.Fmd

        A   = self.ySum(A, self.names)
        fmd = self.ySum(self.fmd, self.names)
        
        b   = self.fmeas - fmd
        return A, b

    # ----------------------------------------------------------------------
    def adjust(self):
        meas = self.fmeas
        calc = self.ySum(self.fed+self.fmd, self.names)
        for band in self.bandspan:
            if len(band) > 1:
                msum = 0.0
                csum = 0.0
                for i in band:
                    msum += meas[i] 
                    csum += calc[i]
                factor = msum / csum
                for i in band:
                    meas[i] = factor * calc[i]
                    
    # ----------------------------------------------------------------------
    def _cauchy_(self):
        w1 = 1.0e4 / self.energy
        w2 = w1 * w1
        w4 = w2 * w2
        c = self.cauchy
        return c[0]/w4 + c[1]/w2 + c[2] + c[3]*w2 + c[4]*w4

    # ----------------------------------------------------------------------
    def calcCauchy(self, energy):
        w1 = 1.0e4 / energy
        w2 = w1 * w1
        w4 = w2 * w2
        c = self.cauchy
        return c[0]/w4 + c[1]/w2 + c[2] + c[3]*w2 + c[4]*w4

    # ----------------------------------------------------------------------
    def _squareMatrix_(self, matrix, eigenvectors):
        matrix = Transform(matrix, eigenvectors)
        return matrix*matrix

    # ----------------------------------------------------------------------
    def _solve_(self, A, b):
        D = Numeric.matrixmultiply(Numeric.transpose(A), A)
        c = Numeric.matrixmultiply(Numeric.transpose(A), b)
        x = LinearAlgebra.solve_linear_equations(D, c)
        return x

    # ----------------------------------------------------------------------
    def calcChi2(self, omega):
        chi2 = self.b - Numeric.matrixmultiply(self.A, omega)
        if self.otype == 0:
            chi2 *= self.sigma
        chi2 = Numeric.add.reduce(chi2*chi2)
        return chi2
        
    # ----------------------------------------------------------------------
    def Show(self, maxenergy=0):
        if self.otype == 0:
            return Show.ShowTradJO(self, maxenergy)
        else:
            return Show.ShowJuddOfelt(self, maxenergy)

    # ----------------------------------------------------------------------
    def calcData(self, parms=""):
        if parms:
            self.H.setParms(parms)
        self.energy,eigenvectors = self.H.diagonalize()
        self.names,mult = self.H.intermediate()
        #print "# -1-", self.energy[:4]

        states = len(self.energy)
        dk = Numeric.zeros((states,states), Numeric.Float)
        dl = Numeric.zeros((states,states), Numeric.Float)
        n  = Numeric.zeros((states,states), Numeric.Float)       
        for i in range(states):
            for f in range(states):
                value = abs(self.energy[i]-self.energy[f])
                dk[i,f] = value
                if value:
                    dl[i,f] = 1.0e7 / value
                else:
                    dl[i,f] = -1
                k = value
                if k < 1000:
                    k = 1000
                if k > 50000:
                    k = 50000
                n[i,f] = self.calcCauchy(k)

        #fed = n*n + 2
        #fed = fed*fed / (9*n)
        #fed *= FED * dk
        #fmd = n * FMD * dk
        #
        #factor = n*n * dk*dk
        #Aed = factor * fed
        #Amd = factor * fmd

        Xed = n*n + 2
        Xed = Xed*Xed / (9*n)
        Xmd = n

        fed = Xed * FED * dk
        fmd = Xmd * FMD * dk

        factor = n*n * dk*dk*dk * COULOMB
        Aed = Xed * AED * factor
        Amd = Xmd * AMD * factor

        for i in range(states):
            for f in range(states):
                if i >= f:
                    fed[i,f] = 0
                    fmd[i,f] = 0
                if i <= f:
                    Aed[i,f] = 0
                    Amd[i,f] = 0

        for i in range(states):
            fed[i,:] /= mult[i]
            fmd[i,:] /= mult[i]
            Aed[i,:] /= mult[i]
            Amd[i,:] /= mult[i]

        u2 = Numeric.zeros((states, states, self.num), Numeric.Float)
        for i in range(self.num):
            u2[:,:,i] = self._squareMatrix_(self.u[i,:,:], eigenvectors)
        #u22 = self._squareMatrix_(self.u2, eigenvectors)
        #u42 = self._squareMatrix_(self.u4, eigenvectors)
        #u62 = self._squareMatrix_(self.u6, eigenvectors)
        ls2 = self._squareMatrix_(self.ls, eigenvectors)

        u = Numeric.matrixmultiply(u2, self.omega)
        #u = self.omega[0]*u22 + self.omega[1]*u42 + self.omega[2]*u62
        fed *= u
        fmd *= ls2
        Aed *= u
        Amd *= ls2

        beta = Aed+Amd
        tau = Numeric.add.reduce(beta, 1)
        tau[1:] = 1.0 / tau[1:]
        for i in range(states):
            beta[i,:] *= tau[i]

        self.tau = tau
        self.Aed = Aed
        self.Amd = Amd
        return self.names,tau,dk,dl,u2,ls2,fed,fmd,Aed,Amd,beta
        #return self.names,tau,dk,dl,u22,u42,u62,ls2,fed,fmd,Aed,Amd,beta

    # ----------------------------------------------------------------------
    def calcDataFull(self, parms=""):
	
	self.calcData(parms)
	
        if parms:
            self.H.setParms(parms)
        self.energy,eigenvectors = self.H.diagonalize()
        self.names,mult = self.H.intermediate()
        #print "# -1-", self.energy[:4]

        states = len(self.energy)
        dk = Numeric.zeros((states,states), Numeric.Float)
        n  = Numeric.zeros((states,states), Numeric.Float)
        for i in range(states):
            for f in range(states):
                value = abs(self.energy[i]-self.energy[f])
                dk[i,f] = value
                k = value
                if k < 1000:
                    k = 1000
                if k > 50000:
                    k = 50000
                n[i,f] = self.calcCauchy(k)

        Xed = n*n + 2
        Xed = Xed*Xed / (9*n)
        Xmd = n

        fed = Xed * FED * dk
        fmd = Xmd * FMD * dk

        factor = n*n * dk*dk*dk * COULOMB
        Aed = Xed * AED * factor
        Amd = Xmd * AMD * factor

        for i in range(states):
            fed[i,:] /= mult[i]
            fmd[i,:] /= mult[i]
            Aed[i,:] /= mult[i]
            Amd[i,:] /= mult[i]

        u2 = Numeric.zeros((states, states, self.num), Numeric.Float)
        for i in range(self.num):
            u2[:,:,i] = self._squareMatrix_(self.u[i,:,:], eigenvectors)
        ls2 = self._squareMatrix_(self.ls, eigenvectors)

        u = Numeric.matrixmultiply(u2, self.omega)
        fed *= u
        fmd *= ls2
        Aed *= u
        Amd *= ls2

        return self.names,dk,n,fed,fmd,Aed,Amd

    # ----------------------------------------------------------------------
    def showEmLaTeX(self, span=0):
        names,tau,dk,dl,u2,ls2,fed,fmd,Aed,Amd,beta = self.calcData()
        #names,tau,dk,dl,u22,u42,u62,ls2,fed,fmd,Aed,Amd,beta = \
        #    self.calcData()
        if span:
            states = span
        else:
            states = len(names)
        for i in range(states):
            S,L,J = splitShortTerm(names[i])
            names[i] = "\\term{%s}{%s}{%s}" % (S,L,J)
        s = ""
        s += "\\begin{EmTab}\n"
        for i in range(states-1,0,-1):
            s += "  \\hline\n"
            s += "  & \\multicolumn{%d}{c|}{initial: %s}\n" % \
                 (3+self.num, names[i])
            #s += "  & \\multicolumn{6}{c|}{initial: %s}\n" % names[i]
            s += "  & \\multicolumn{3}{c|}" \
                 + "{$\\tau=\\einheit{%.2f}{ms}$} \\\\\n" % (tau[i]*1000)
            s += "  \\hline\n"
            for f in range(i-1,-1,-1):
                s += "  %s\n" % names[f]
                s += "  & %.0f\n" % dk[i,f]
                if dl[i,f] < 8000:
                    s += "  & %.0f\n" % dl[i,f]
                else:
                    s += "  & --\n"
                s += " "
                for x in range(self.num):
                    s += " & %.4f" % u2[i,f,x]
                s += " & %.4f\n" % ls2[i,f]
                #s += "  & %.4f & %.4f & %.4f & %.4f\n" % \
                #     (u22[i,f], u42[i,f], u62[i,f], ls2[i,f])
                s += "  & %.0f & %.0f & %.3f \\\\\n" % \
                     (Aed[i,f], Amd[i,f], beta[i,f])
        s += "\\end{EmTab}\n"
        s += "\n"
        return s[:-1]

    # ----------------------------------------------------------------------
    def showAbsLaTeX(self, span=0):
        names,tau,dk,dl,u2,ls2,fed,fmd,Aed,Amd,beta = self.calcData()
        #names,tau,dk,dl,u22,u42,u62,ls2,fed,fmd,Aed,Amd,beta = \
        #    self.calcData()
        if span:
            states = span
        else:
            states = len(names)
        for i in range(states):
            S,L,J = splitShortTerm(names[i])
            names[i] = "\\term{%s}{%s}{%s}" % (S,L,J)
        s = ""
        s += "\\begin{AbsTab}\n"
        for i in range(states-1):
            s += "  \\hline\n"
            s += "  & \\multicolumn{%d}{c|}{initial: %s}\n" % \
                 (3+self.num, names[i])
            #s += "  & \\multicolumn{6}{c|}{initial: %s}\n" % names[i]
            s += "  & \\multicolumn{2}{c|}{} \\\\\n"
            s += "  \\hline\n"
            for f in range(i+1,states):
                s += "  %s\n" % names[f]
                s += "  & %.0f\n" % dk[i,f]
                if dl[i,f] < 8000:
                    s += "  & %.0f\n" % dl[i,f]
                else:
                    s += "  & --\n"
                s += " "
                for l in range(self.num):
                    s += " & %.4f" % u2[i,f,l]
                s += " & %.4f\n" % ls2[i,f]
                #s += "  & %.4f & %.4f & %.4f & %.4f\n" % \
                #     (u22[i,f], u42[i,f], u62[i,f], ls2[i,f])
                s += "  & %.1f & %.1f \\\\\n" % \
                     (fed[i,f]*1e8, fmd[i,f]*1e8)
        s += "\\end{AbsTab}\n"
        s += "\n"
        return s[:-1]
