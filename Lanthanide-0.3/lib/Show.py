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
from Fit import Fit
import os.path
import Hamiltonian
import TermFit
import TermNames
import JuddOfelt



FREE  = 0  # not measured, not fixed
FIXED = 1  # not measured, fixed
FLOAT = 2  # measured, not fixed
MEAS  = 3  # measured, fixed
START = 4  # span: meaured and fixed (first)
SPAN  = 5  # span: meaured and fixed (cont.)


# --------------------------------------------------------------------------
class Show(Fit):
    fp = ""
    mlevel = []

    def open(self, fn, mode):
        if fn:
            path,file = os.path.split(fn)
            if path and not os.path.exists(path):
                os.makedirs(path)
            self.fp = open(fn, mode)
        else:
            self.fp = ""

    def close(self):
        if self.fp:
            self.fp.close()
            self.fp = ""

    def lastIndex(self, energy, maxenergy=0):
        if maxenergy > 0:
            for i in range(len(energy)):
                if energy[i] > maxenergy:
                    break
            self.last = i
        else:
            self.last = len(energy)

    def showTitle(self, s):
        head = " %s " % s
        self.showString(head)
        self.showString(len(head)*"=")

    def showString(self, s="", stdout=1):
        if self.stdout:
            if stdout:
                print s
        if self.fp:
            self.fp.write(s+"\n")

    def showResult(self, span=0):
        self.showLine()
	self.sortValues()
	self.showValues(span)
        self.showString()
	self.showStats()

    def showLine(self):
        s = ""
        for line in self.form:
            if len(line) == 1:
                s += "+"
            else:
                num = len(line[1] % "")
                s += num * "-"
        self.showString(s)

    def showData(self, data, head=0):
        s = ""
        pos = 0
        for line in self.form:
            if len(line) == 1:
                s += line[0]
            else:
                num = len(line[1] % "")
                s += self.showElement(line, data[pos], num)
                pos += 1
        self.showString(s)

    def showElement(self, form, data, num):
        if data == "":
            #s = num * " "
            s = form[1] % data
        elif type(data) == type(""):
            s = form[1] % data
        else:
            if type(data) == type((0,)):
                s = []
                for val in data:
                    s.append(str(val))
                s = " " + " + ".join(s)
            else:
                s = form[0] % data
        return s
    

    def sortValues(self):
        nummeas = len(self.span)
        states  = len(self.names)

        # Prepare lists status and mlevel
        status = []
        mlevel = []
        for i in range(states):
            status.append(FREE)
            mlevel.append(-1)

        # Set fixed values, not measured and measured
        nofix = []
        for i in range(nummeas):
            span = self.span[i]
            name = self.guess[i]
            if type(span) != types.TupleType:
                if span < 0:
                    nofix.append(i)
                else:
                    sigma = self.sigma[i]
                    level = self.yIndex(name, span, self.names)
                    if sigma < 0:
                        status[level]  = FIXED
                    else:
                        status[level]  = MEAS
                    mlevel[level] = i
            else:
                for j in range(len(span)):
                    n = name[j]
                    k = span[j]
                    level = self.yIndex(n, k, self.names)
                    if j == 0:
                        status[level] = START
                    else:
                        status[level] = SPAN
                    mlevel[level] = i
                    
        # Build list of free levels
        free = []
        for i in range(states):
            if status[i] == FREE:
                free.append(i)

        # Find closest value for all floating values
        while (len(free) != 0) and (len(nofix) != 0):
            ### ERROR: in double diff = calcDiff()[0] + calcDiff()[1]
            mindiff = self.calcDiff(free[0], nofix[0])[0]
            imin    = 0
            jmin    = 0
            for i in range(len(free)):
                for j in range(len(nofix)):
                    ### ERROR: in double diff = calcDiff()[0] + calcDiff()[1]
                    diff = self.calcDiff(free[i], nofix[j])[0]
                    if (diff < mindiff):
                        mindiff = diff
                        imin    = i
                        jmin    = j
            level          = free[imin]
            status[level]  = FLOAT
            mlevel[level]  = nofix[jmin]
            del nofix[jmin]
            del free[imin]

        self.status = status
        self.mlevel = mlevel


    def showValues(self, span=0):
        totvalues = 0
        fitvalues = 0
        diff      = 0

        self.spanindex = Numeric.zeros(len(self.span))
        self.showHead()
        states = self.last
        if span:
            states = span
        for i in range(states):
            j = self.mlevel[i]
            if self.status[i] == FREE:
                # not measured, not fixed
                self.showFree(i, j)
            if self.status[i] == FIXED:
                # not measured, but fixed
                self.showFixed(i, j)
            if self.status[i] == FLOAT:
                # measured, but not fixed
                diff += self.showFloat(i, j)
                totvalues += self.valperlevel
            if self.status[i] == MEAS:
                # measured and fixed
                diff += self.showMeas(i, j)
                totvalues += self.valperlevel
                fitvalues += self.valperlevel
            if self.status[i] == START:
                # span, measured and fixed (first)
                diff += self.showStart(i, j)
                totvalues += self.valperlevel
                fitvalues += self.valperlevel
                self.spanindex[j] += 1
            if self.status[i] == SPAN:
                # span, measured and fixed (cont.)
                self.showSpan(i, j)
                self.spanindex[j] += 1
        self.showLine()
        self.totvalues = totvalues
        self.fitvalues = fitvalues
        self.showDiff(diff)


    def showStats(self):
        self.showString("fit parameters:    %4i"    % self.fitparms)
        self.showString("fit values:        %4i"    % self.fitvalues)
        self.showString("total parameters:  %4i"    % self.totparms)
        self.showString("total meas values: %4i"    % self.totvalues)


    def calcStats(self, diff):
        num = len(diff) / 3
        chi2  = []
        sum   = []
        sigma = []
        for i in range(num):
            chi2.append(diff[i])
            sum.append(diff[2*num+i])
            #sigma.append(math.sqrt(diff[i] * diff[num+i])/self.fitvalues)
            sigma.append(math.sqrt(diff[i]*diff[num+i])/self.fitvalues \
                         * self.valperlevel)
        if num == 1:
            return chi2[0], sum[0], sigma[0]
        return chi2, sum, sigma



# --------------------------------------------------------------------------
class ShowJOparms:
    def showOmega(self):
        omega = self.JO.omega * 1e24

        s = ""
        form = " %-16s | %12s \n"
        s += form % ("parameter", "value")
        s += form % ("", "pm^2")
        s += 18*"-" + "+" + 14*"-" + "\n"

        index = [2, 4, 6]
        if len(omega) != 3:
            index = [1, 2, 3, 4, 5, 6]
        for i in range(len(omega)):
            name = "omega%i" % index[i]
            #name = "omega%i" % (2*i+2)
            new  = "%.3f" % omega[i]
            s += form % (name, new)
        self.showString(s)

    def showJOdata(self, LaTeX="", span=0):
        if not self.JO.energy:
            self.JO.fit()

        if LaTeX == "e":
            self.showString(self.JO.showEmLaTeX(span)[:-1])
            return
        if LaTeX == "a":
            self.showString(self.JO.showAbsLaTeX(span)[:-1])
            return
            
        s = ""
        form = " %-16s | %10s  %10s "
        vals = ["level", "kcalc", "fmd"]
        index = [2, 4, 6]
        if self.JO.num != 3:
            index = [1, 2, 3, 4, 5, 6]
        for i in range(self.JO.num):
            form += " %10s "
            vals.append("|<U%d>|^2" % index[i])
        form += "\n"
        s += form % tuple(vals)
        #form = " %-16s | %10s  %10s  %10s  %10s  %10s \n"
        #s += form % ("level", "kcalc", "fmd",
        #              "|<U2>|^2", "|<U4>|^2", "|<U6>|^2")
        vals = ["", "cm^-1", "10^-8"] + self.JO.num * [""]
        s += form % tuple(vals)
        s += 18*"-" + "+" + (2+self.JO.num)*12*"-" + "\n"
        
        for i in range(len(self.JO.names)):
            name = self.JO.names[i]
            energy = "%i" % self.JO.energy[i]
            u = []
            for x in range(self.JO.num):
                val = self.JO.u2[x,i]
                val = round(val, 4)
                if val == 0.0:
                    val = "0      "
                else:
                    val = "%.4f " % val
                u.append(val)
            #u = [ self.JO.u22[i], self.JO.u42[i], self.JO.u62[i] ]
            #for j in range(3):
            #    u[j] = round(u[j], 4)
            #    if u[j] == 0.0:
            #        u[j] = "0      "
            #    else:
            #        u[j] = "%.4f " % u[j]
            fmd = self.JO.fmd[i]*1e8
            if ZERO(fmd):
                fmd = "0  "
            else:
                fmd = "%.1f" % fmd
            vals = [name, energy, fmd]
            for x in range(self.JO.num):
                vals.append(u[x])
            s += form % tuple(vals)
            #s += form % (name, energy, fmd, u[0], u[1], u[2])
        self.showString(s[:-1])



# --------------------------------------------------------------------------
class ShowJuddOfelt(Show, ShowJOparms):
    form = [
        [ " %2i ",    " %2s "   ],
        [ " %-12s ",  " %-12s " ],
        [ "|"                   ],
        [ " %8.1f ",  " %8s "   ],
        [ "(%4.1f) ", "%6s "    ],
        [ " %8.1f ",  " %8s "   ],
        [ " %8.1f ",  " %8s "   ],
        [ " %8.1f ",  " %8s "   ],
        [ "|"                   ],
        [ " %-12s ",  " %-12s " ],
        ]
    
    def __init__(self, JO, maxenergy=0):
        self.JO     = JO
        self.JO._optParms_()
        self.names  = JO.names
        self.energy = JO.energy
        self.fed    = JO.fed * 1
        self.fmd    = JO.fmd * 1
        self.span,self.guess,self.fmeas,self.sigma = JO.splitMeas(JO.meas, 1)
        
        self.lastIndex(self.energy, maxenergy)
        
        self.totparms = len(JO.omega)
        self.fitparms = len(JO.omega)
        self.valperlevel = 1
        
        self.fed   *= 1e8
        self.fmd   *= 1e8
        self.fmeas *= 1e8
        self.sigma *= 1e8

        self.stdout = 1

    def showHead(self):
        self.showData([ "", "level",
                        "fmeas", "", "fed", "fmd", "diff",
                        "guess" ], 1)
	unit = "10^-8"
        self.showData([ "", "",
                        unit, "", unit, unit, unit,
                        "" ], 1)
        
    def showFree(self, i, j):
        self.showData([ i, self.names[i],
                        "", "", self.fed[i], self.fmd[i], "",
                        "" ])

    def showFixed(self, i, j):
        self.showData([ i, self.names[i],
                        "", "", self.fed[i], self.fmd[i], "",
                        self.guess[j] ])

    def showFloat(self, i, j):
        guess = "..%s.." % self.guess[j]
	diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.fmeas[j], self.sigma[j],
                        self.fed[i], self.fmd[i], diff[2],
                        guess ])
        return diff

    def showMeas(self, i, j):
        diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.fmeas[j], self.sigma[j],
                        self.fed[i], self.fmd[i], diff[2],
                        self.guess[j] ])
        return diff

    def showStart(self, i, j):
        guess = self.guess[j][0]
        diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.fmeas[j], self.sigma[j],
                        self.fed[i], self.fmd[i], diff[2],
                        guess ])
        return diff

    def showSpan(self, i, j):
        k = self.spanindex[j]
        guess = "%s (%s)" % (self.guess[j][k], self.guess[j][0])
        self.showData([ i, self.names[i],
                        "...", "",
                        self.fed[i], self.fmd[i], "...",
                        guess ])

    def showDiff(self, diff):
        chi2,sum,sigma = self.calcStats(diff)
        self.stats = (chi2,sum,sigma)
        self.showData([ "", "sum",
                        "", "", "", "", sum,
                        "" ])
        self.showData([ "", "sigma",
                        "", "", "", "", sigma,
                        "" ])
        self.showData([ "", "chi2",
                        "", "", "", "", chi2,
                        "" ])

    def calcDiff(self, icalc, imeas):
        span = self.span[imeas]
        if type(span) != types.TupleType:
            fcalc = self.fed[icalc] + self.fmd[icalc]
        else:
            fcalc = 0.0
            for icalc in span:
                fcalc += self.fed[icalc] + self.fmd[icalc]
        diff  = fcalc - self.fmeas[imeas]
        #diff2 = diff * diff
        diff2 = pow(self.sigma[imeas], 2)
        chi2  = pow(diff / self.sigma[imeas], 2)
        return Numeric.array([chi2, diff2, diff], Numeric.Float)


# --------------------------------------------------------------------------
class ShowTradJO(ShowJuddOfelt):
    def calcStats(self):
        raise RuntimeError, "ShowTradJO: not defined! :-("
        chi2 = self.chi2 * 1e16
        if self.fitparms < self.fitvalues:
            sigma = "%5.1f 10^-8" % sigma
        else:
            sigma = "..."

        chi2  = "%5.1f 10^-16" % chi2
        rms   = "%5.1f 10^-8" % rms
        return chi2, rms, sigma



# --------------------------------------------------------------------------
class ShowParms:
    def showParms(self, oldparms=0):
        if oldparms:
            self.showParmsLong()
            return

        s = ""
        form = " %-16s | %12s \n"
        s += form % ("parameter", "value ")
        s += form % ("", "cm^-1 ")
        s += 18*"-" + "+" + 14*"-" + "\n"

        for i in range(len(self.TF.parms)):
            value = self.TF.parms[i][0]
            name  = self.TF.parms[i][1]
            if self.TF.parms[i][2] != 0:
                new = " %.2f " % value
            else:
                new = "(%.2f)" % value
            s += form % (name, new)
        self.showString(s[:-1])

    def showParmsLong(self):
        s = ""
        form = " %-16s | %12s  %12s \n"
        s += form % ("parameter", "value ", "start ")
        s += form % ("", "cm^-1 ", "cm^-1 ")
        s += 18*"-" + "+" + 14*"-" + 14*"-" + "\n"

        for i in range(len(self.TF.parms)):
            val1 = self.TF.parms[i][0]
            val2 = self.TF.oldparms[i][0]
            name = self.TF.parms[i][1]
            if self.TF.parms[i][2] != 0:
                new = " %.2f " % val1
                old = " %.2f " % val2
            else:
                new = "(%.2f)" % val1
                old = "(%.2f)" % val2
            s += form % (name, new, old)
        self.showString(s[:-1])

    def sortVectors(self, matspan=(), LaTeX=0):
        H = self.TF.H
        if not H.eigenvectors:
            H.diagonalize()

        cnames = TermNames.short(H.term, H.key)
        syms   = H.term[H.key]["symmetry"]
        cvals  = H.term[H.key]["values"][:,syms.index("J2")]
        rnames,mult = H.intermediate()
        states = len(cnames)
        
        cindex = Numeric.argsort(cvals)
        rindex = Numeric.argsort(mult)

        matrix = H.eigenvectors * H.eigenvectors

        # calculate center of mass (cm) energy of J-multiplets
        rvals = Numeric.take(mult, rindex)
        energy = Numeric.take(H.energies, rindex)
        mean = []
        span = []
        first = 0
        sum = 0.0
        for i in range(states):
            if rvals[i] != rvals[first]:
                mean.append(sum / (i-first))
                span.append((first, i))
                first = i
                sum = 0.0
            sum += energy[i]
        mean.append(sum / (states-first))
        span.append((first, states))
        mean = Numeric.array(mean)

        ## sort column and row indices to cm energy of J-multiplets
        #rlist = []
        #clist = []
        #index = Numeric.argsort(mean)
        #for i in index:
        #    first,last = span[i]
        #    rlist += Numeric.sort(rindex[first:last]).tolist()
        #    clist += Numeric.sort(cindex[first:last]).tolist()
        #rindex = Numeric.array(rlist)
        #cindex = Numeric.array(clist)

        # sort column index to get maximum values on the diagonal
        cindex = cindex.tolist()
        cmax = Numeric.argmax(matrix)
        for i in range(states):
            j = cindex.index(cmax[rindex[i]])
            if j > i:
                cindex[i],cindex[j] = SWAP(cindex[i],cindex[j])
        cindex = Numeric.array(cindex)

        # sort sqared matrix by new row and  column index and keep signs
        matrix = H.eigenvectors * Numeric.absolute(H.eigenvectors)
        matrix = Numeric.take(matrix, rindex, 0)
        matrix = Numeric.take(matrix, cindex, 1)

        # build list of column and row names
        newcnames = []
        newrnames = []
        for i in range(states):
            newcnames.append(cnames[cindex[i]])
            newrnames.append("%i" % rindex[i])

        # build span list of sub-matrices
        rvals = Numeric.take(mult, rindex)
        span = []
        first = 0
        for i in range(states):
            if rvals[i] != rvals[first]:
                span.append((first, i))
                first = i
        span.append((first, states))

        matrix = Numeric.transpose(matrix)
        return matrix, newcnames, newrnames, span


    def showVectors(self, matspan=()):
        matrix, names, nums, span = self.sortVectors()
        s = MatrixToString(matrix, names, nums, matspan)
        self.showString(s)
        

    def showVectorsLaTeX(self, matspan=()):
        matrix, names, nums, span = self.sortVectors()
        matrix = Numeric.transpose(matrix)

        for i in range(len(names)):
            S,L,J = splitShortTerm(names[i])
            names[i] = "\\term{%s}{%s}{%s}" % (S,L,J)

        xmax = 5
        for first,last in span:
            xmax = max(xmax, last-first)
        s = ""
        s += "\\begin{VectorTab}{%i}\n" % xmax
        for first,last in span:
            num = last - first
            c = (xmax - num) * " &\phantom{0.0000}" + " \\\\\n"
            #s += "  \\hline\n"
            s += "  \\multicolumn{1}{|c}{} "
            for i in range(first, last):
                s += " & %s" % names[i]
            s += c
            s += "  \\cline{2-%i}\n" % (xmax+1)
            for i in range(first, last):
                s += "  %s\n  " % nums[i]
                for j in range(first, last):
                    value = matrix[i,j]
                    if value > 0:
                        s += "&%.4f" % value
                    else:
                        s += "&$\\overline{%.4f}$" % (-value)
                s += c
            s += "  \\hline\n"
        s += "\\end{VectorTab}\n"
        self.showString(s)
        
            
# --------------------------------------------------------------------------
class ShowLevels(Show, ShowParms):
    form = [
        [ " %2i ",    " %2s "   ],
        [ " %-12s ",  " %-12s " ],
        [ "|"                   ],
        [ " %8.0f ",  " %8s "   ],
        [ "(%4.0f) ", "%6s "    ],
        [ " %8.0f ",  " %8s "   ],
        [ " %8.1f ",  " %8s "   ],
        [ "|"                   ],
        [ " %-12s ",  " %-12s " ],
        ]
    
    def __init__(self, TF, maxenergy=0):
        self.TF       = TF
        self.energy,v = TF.H.diagonalize()
        self.names,self.mult = TF.H.intermediate()
        self.span,self.guess,self.meas,self.sigma = TF.splitMeas(TF.meas, 1)

        self.lastIndex(self.energy, maxenergy)

	self.totparms = len(TF.parms)
	self.fitparms = 0
	for i in range(self.totparms):
	    if TF.parms[i][2] > 0:
		self.fitparms += 1
        self.valperlevel = 1

        self.stdout = 1

    def showHead(self):
        self.showData([ "", "level",
                        "kmeas", "", "kcalc", "diff",
                        "guess" ], 1)
	unit = "cm^-1"
        self.showData([ "", "",
                        unit, "", unit, unit,
                        "" ], 1)
        
    def showFree(self, i, j):
	self.showData([ i, self.names[i],
                        "", "", self.energy[i], "",
                        "" ])

    def showFixed(self, i, j):
	self.showData([ i, self.names[i],
                        "", "", self.energy[i], "",
                        self.guess[j] ])

    def showFloat(self, i, j):
        guess = "..%s.." % self.guess[j]
	diff = self.calcDiff(i, j)
	self.showData([ i, self.names[i],
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[2],
                        guess ])
	return diff

    def showMeas(self, i, j):
	diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[2],
                        self.guess[j] ])
	return diff

    def showStart(self, i, j):
        guess = self.guess[j][0]
        diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[2],
                        guess ])
	return diff

    def showSpan(self, i, j):
        k = self.spanindex[j]
        guess = "%s (%s)" % (self.guess[j][k], self.guess[j][0])
        self.showData([ i, self.names[i],
                        "...", "",
                        self.energy[i], "...",
                        guess ])

    def showDiff(self, diff):
        chi2,sum,sigma = self.calcStats(diff)
        self.stats = (chi2,sum,sigma)
        self.showData([ "", "sum",
                        "", "", "", sum,
                        "" ])
        self.showData([ "", "sigma",
                        "", "", "", sigma,
                        "" ])
        self.showData([ "", "chi2",
                        "", "", "", chi2,
                        "" ])

    def calcDiff(self, icalc, imeas):
        span = self.span[imeas]
        if type(span) != types.TupleType:
            energy = self.energy[icalc]
        else:
            energy = 0.0
            mult   = 0
            for icalc in span:
                energy += self.energy[icalc] * self.mult[icalc]
                mult   += self.mult[icalc]
            energy /= mult
        diff  = energy - self.meas[imeas]
        #diff2 = diff * diff
        diff2 = pow(self.sigma[imeas], 2)
        chi2  = pow(diff / self.sigma[imeas], 2)
        return Numeric.array([chi2, diff2, diff], Numeric.Float)



# --------------------------------------------------------------------------
class ShowDouble(Show, ShowParms, ShowJOparms):
    form = [
        [ " %2i ",    " %2s "   ],
        [ " %-12s ",  " %-12s " ],
        [ "|"                   ],
        [ " %8.0f ",  " %8s "   ],
        [ "(%4.0f) ", "%6s "    ],
        [ " %8.0f ",  " %8s "   ],
        [ " %8.1f ",  " %8s "   ],
        [ "|"                   ],
        [ " %8.1f ",  " %8s "   ],
        [ "(%4.1f) ", "%6s "    ],
        [ " %8.1f ",  " %8s "   ],
        [ " %8.1f ",  " %8s "   ],
        [ " %8.1f ",  " %8s "   ],
        [ "|"                   ],
        [ " %-12s ",  " %-12s " ],
        ]
    
    def __init__(self, TF, JO, maxenergy=0):
        self.JO     = JO
        self.names  = JO.names
        self.fed    = JO.fed * 1
        self.fmd    = JO.fmd * 1
        self.span,self.guess,self.fmeas,self.fsigma = JO.splitMeas(JO.meas, 1)
        
        self.fed    *= 1e8
        self.fmd    *= 1e8
        self.fmeas  *= 1e8
        self.fsigma *= 1e8

        self.TF       = TF
        self.energy,v = TF.H.diagonalize()
        self.names,self.mult = TF.H.intermediate()
        self.span,self.guess,self.meas,self.sigma = TF.splitMeas(TF.meas, 1)

        self.lastIndex(self.energy, maxenergy)

	self.totparms = len(TF.parms)
	self.fitparms = 0
	for i in range(self.totparms):
	    if TF.parms[i][2] > 0:
		self.fitparms += 1
        self.totparms += len(JO.omega)
        self.fitparms += len(JO.omega)
        self.valperlevel = 2

        self.stdout = 1

    def showHead(self):
        self.showData([ "", "level",
                        "kmeas", "", "kcalc", "diff",
                        "fmeas", "", "fed", "fmd", "diff",
                        "guess" ], 1)
	unit1 = "cm^-1"
	unit2 = "10^-8"
        self.showData([ "", "",
                        unit1, "", unit1, unit1,
                        unit2, "", unit2, unit2, unit2,
                        "" ], 1)

    def showFree(self, i, j):
	self.showData([ i, self.names[i],
                        "", "", self.energy[i], "",
                        "", "", self.fed[i], self.fmd[i], "",
                        "" ])

    def showFixed(self, i, j):
	self.showData([ i, self.names[i],
                        "", "", self.energy[i], "",
                        "", "", self.fed[i], self.fmd[i], "",
                        self.guess[j] ])

    def showFloat(self, i, j):
        guess = "..%s.." % self.guess[j]
        diff = self.calcDiff(i, j)
	self.showData([ i, self.names[i],
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[4],
                        self.fmeas[j], self.fsigma[j],
                        self.fed[i], self.fmd[i], diff[5],
                        guess ])
	return diff

    def showMeas(self, i, j):
        diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[4],
                        self.fmeas[j], self.fsigma[j],
                        self.fed[i], self.fmd[i], diff[5],
                        self.guess[j] ])
	return diff

    def showStart(self, i, j):
        guess = self.guess[j][0]
        diff = self.calcDiff(i, j)
        self.showData([ i, self.names[i],
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[4],
                        self.fmeas[j], self.fsigma[j],
                        self.fed[i], self.fmd[i], diff[5],
                        guess ])
	return diff

    def showSpan(self, i, j):
        k = self.spanindex[j]
        guess = "%s (%s)" % (self.guess[j][k], self.guess[j][0])
        self.showData([ i, self.names[i],
                        "...", "",
                        self.energy[i], "...",
                        "...", "",
                        self.fed[i], self.fmd[i], "...",
                        guess ])

    def showDiff(self, diff):
        chi2,sum,sigma = self.calcStats(diff)
        self.stats = (chi2,sum,sigma)
        self.showData([ "", "sum",
                        "", "", "", sum[0],
                        "", "", "", "", sum[1],
                        "" ])
        self.showData([ "", "sigma",
                        "", "", "", sigma[0],
                        "", "", "", "", sigma[1],
                        "" ])
        self.showData([ "", "chi2",
                        "", "", "", chi2[0],
                        "", "", "", "", chi2[1],
                        "" ])

    def plotLists(self):
        if not self.mlevel:
            self.sortValues()
        listCalc = []
        listMeas = []
        for i in range(self.last):
            j = self.mlevel[i]
            if self.status[i] == FREE:
                # not measured, not fixed
                listCalc.append(self.listCalc(i))
            if self.status[i] == FIXED:
                # not measured, but fixed
                listCalc.append(self.listCalc(i))
            if self.status[i] == FLOAT:
                # measured, but not fixed
                listMeas.append(self.listMeas(i, j))
            if self.status[i] == MEAS:
                # measured and fixed
                listMeas.append(self.listMeas(i, j))
            if self.status[i] == START:
                # span, measured and fixed (first)
                listMeas.append(self.listMeas(i, j))
            if self.status[i] == SPAN:
                # span, measured and fixed (cont.)
                pass
        return listCalc, listMeas

    def listCalc(self, i):
        return [self.energy[i], self.fed[i]+self.fmd[i]]

    def listMeas(self, i, j):
        energy,fcalc = self.calcSum(i, j)
        return [self.meas[j], self.sigma[j], energy,
                self.fmeas[j], self.fsigma[j], fcalc]

    def calcSum(self, icalc, imeas):
        span = self.span[imeas]
        if type(span) != types.TupleType:
            energy = self.energy[icalc]
            fcalc  = self.fed[icalc] + self.fmd[icalc]
        else:
            fcalc  = 0.0
            energy = 0.0
            mult   = 0
            for icalc in span:
                fcalc  += self.fed[icalc] + self.fmd[icalc]
                energy += self.energy[icalc] * self.mult[icalc]
                mult   += self.mult[icalc]
            energy /= mult
        return energy, fcalc
    
    def calcDiff(self, icalc, imeas):
        energy,fcalc = self.calcSum(icalc, imeas)
        diff1  = energy - self.meas[imeas]
        diff2  = fcalc - self.fmeas[imeas]
        #diff12 = diff1 * diff1
        #diff22 = diff2 * diff2
        diff12 = pow(self.sigma[imeas], 2)
        diff22 = pow(self.fsigma[imeas], 2)
        chi12  = pow(diff1 / self.sigma[imeas], 2)
        chi22  = pow(diff2 / self.fsigma[imeas], 2)
        return Numeric.array([chi12, chi22, diff12, diff22, diff1, diff2],
                             Numeric.Float)

# --------------------------------------------------------------------------
class ShowDoubleLaTeX(Show, ShowParms, ShowJOparms):
    form = [
        [ "  %i &",       "  %s &"    ],
        [ " %s\n  &",     " %s\n  &"  ],
        [ " $%.0f$ &",    " %s &"     ],
        [ " $(%.0f)$ &",  " %s &"     ],
        [ " $%.0f$ &",    " %s &"     ],
        [ " $%.1f$\n  &", " %s\n  &"  ],
        [ " $%.1f$ &",    " %s &"     ],
        [ " $(%.1f)$ &",  " %s &"     ],
        [ " $%.1f$ &",    " %s &"     ],
        [ " $%.1f$ &",    " %s &"     ],
        [ " $%.1f$ \\\\", " %s \\\\"  ],
        ]
    
    def __init__(self, TF, JO, maxenergy=0):
        self.JO     = JO
        self.names  = JO.names
        self.fed    = JO.fed * 1
        self.fmd    = JO.fmd * 1
        self.span,self.guess,self.fmeas,self.fsigma = JO.splitMeas(JO.meas, 1)
        
        self.fed    *= 1e8
        self.fmd    *= 1e8
        self.fmeas  *= 1e8
        self.fsigma *= 1e8

        self.TF       = TF
        self.energy,v = TF.H.diagonalize()
        self.names,self.mult = TF.H.intermediate()
        self.span,self.guess,self.meas,self.sigma = TF.splitMeas(TF.meas, 1)

        self.lastIndex(self.energy, maxenergy)

	self.totparms = len(TF.parms)
	self.fitparms = 0
	for i in range(self.totparms):
	    if TF.parms[i][2] > 0:
		self.fitparms += 1
        self.totparms += len(JO.omega)
        self.fitparms += len(JO.omega)
        self.valperlevel = 2

        self.stdout = 0

    def showResult(self, span=0):
        self.showString("\\begin{FitTab}")
	self.sortValues()
	self.showValues(span)
        self.showString("\\end{FitTab}")

    def showHead(self):
        pass

    def showLine(self):
        self.showString("  \\hline")

    def showFree(self, i, j):
        S,L,J = splitShortTerm(self.names[i])
        name = "\\term{%s}{%s}{%s}" % (S,L,J)
	self.showData([ i, name,
                        "", "", self.energy[i], "",
                        "", "", self.fed[i], self.fmd[i], "" ])

    def showFixed(self, i, j):
        S,L,J = splitShortTerm(self.names[i])
        name = "\\term{%s}{%s}{%s}" % (S,L,J)
	self.showData([ i, name,
                        "", "", self.energy[i], "",
                        "", "", self.fed[i], self.fmd[i], "" ])

    def showFloat(self, i, j):
        S,L,J = splitShortTerm(self.names[i])
        name = "\\term{%s}{%s}{%s}" % (S,L,J)
        guess = "..%s.." % self.guess[j]
        diff = self.calcDiff(i, j)
	self.showData([ i, name,
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[4],
                        self.fmeas[j], self.fsigma[j],
                        self.fed[i], self.fmd[i], diff[5] ])
	return diff

    def showMeas(self, i, j):
        S,L,J = splitShortTerm(self.names[i])
        name = "\\term{%s}{%s}{%s}" % (S,L,J)
        diff = self.calcDiff(i, j)
        self.showData([ i, name,
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[4],
                        self.fmeas[j], self.fsigma[j],
                        self.fed[i], self.fmd[i], diff[5] ])
	return diff

    def showStart(self, i, j):
        S,L,J = splitShortTerm(self.names[i])
        name = "\\term{%s}{%s}{%s}" % (S,L,J)
        guess = self.guess[j][0]
        diff = self.calcDiff(i, j)
        self.showData([ i, name,
                        self.meas[j], self.sigma[j],
                        self.energy[i], diff[4],
                        self.fmeas[j], self.fsigma[j],
                        self.fed[i], self.fmd[i], diff[5] ])
	return diff

    def showSpan(self, i, j):
        S,L,J = splitShortTerm(self.names[i])
        name = "\\term{%s}{%s}{%s}" % (S,L,J)
        k = self.spanindex[j]
        guess = "%s (%s)" % (self.guess[j][k], self.guess[j][0])
        self.showData([ i, name,
                        "\\ldots", "",
                        self.energy[i], "\ldots",
                        "\ldots", "",
                        self.fed[i], self.fmd[i], "\ldots" ])

    def showDiff(self, diff):
        chi2,sum,sigma = self.calcStats(diff)
        self.stats = (chi2,sum,sigma)
        s = ""
        s += "  \\multicolumn{2}{|c|}{}\n"
        s += "  & \\multicolumn{1}{l}{$\dkbar$:}\n"
        s += "  & & & $%.1f$\n" % sigma[0]
        s += "  & \\multicolumn{1}{l}{$\dfbar$:}\n"
        s += "  & & & & $%.1f$ \\\\" % sigma[1]
        self.showString(s)

    def plotLists(self):
        if not self.mlevel:
            self.sortValues()
        listCalc = []
        listMeas = []
        for i in range(self.last):
            j = self.mlevel[i]
            if self.status[i] == FREE:
                # not measured, not fixed
                listCalc.append(self.listCalc(i))
            if self.status[i] == FIXED:
                # not measured, but fixed
                listCalc.append(self.listCalc(i))
            if self.status[i] == FLOAT:
                # measured, but not fixed
                listMeas.append(self.listMeas(i, j))
            if self.status[i] == MEAS:
                # measured and fixed
                listMeas.append(self.listMeas(i, j))
            if self.status[i] == START:
                # span, measured and fixed (first)
                listMeas.append(self.listMeas(i, j))
            if self.status[i] == SPAN:
                # span, measured and fixed (cont.)
                pass
        return listCalc, listMeas

    def listCalc(self, i):
        return [self.energy[i], self.fed[i]+self.fmd[i]]

    def listMeas(self, i, j):
        energy,fcalc = self.calcSum(i, j)
        return [self.meas[j], self.sigma[j], energy,
                self.fmeas[j], self.fsigma[j], fcalc]

    def calcSum(self, icalc, imeas):
        span = self.span[imeas]
        if type(span) != types.TupleType:
            energy = self.energy[icalc]
            fcalc  = self.fed[icalc] + self.fmd[icalc]
        else:
            fcalc  = 0.0
            energy = 0.0
            mult   = 0
            for icalc in span:
                fcalc  += self.fed[icalc] + self.fmd[icalc]
                energy += self.energy[icalc] * self.mult[icalc]
                mult   += self.mult[icalc]
            energy /= mult
        return energy, fcalc
    
    def calcDiff(self, icalc, imeas):
        energy,fcalc = self.calcSum(icalc, imeas)
        diff1  = energy - self.meas[imeas]
        diff2  = fcalc - self.fmeas[imeas]
        #diff12 = diff1 * diff1
        #diff22 = diff2 * diff2
        diff12 = pow(self.sigma[imeas], 2)
        diff22 = pow(self.fsigma[imeas], 2)
        chi12  = pow(diff1 / self.sigma[imeas], 2)
        chi22  = pow(diff2 / self.fsigma[imeas], 2)
        return Numeric.array([chi12, chi22, diff12, diff22, diff1, diff2],
                             Numeric.Float)
