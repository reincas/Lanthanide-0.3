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

PATT_S = "([0-9]+)"
PATT_L = "([A-Z](\([0-9]+\))?)"
PATT_J = "([0-9]+(/[0-9]+)?)"
PATT_SLJ = PATT_S + PATT_L + PATT_J
PATT_SL = PATT_S + PATT_L


class SigmaPlot:
    def __init__(self, sigma, energy, names, element, base, xmax):
        self.sigma = sigma
        self.energy = energy
        self.names = names
        self.element = element
        self.base = base
        self.xmax = xmax

    def plotHead(self):
        base = "\"ground level ^{%s}%s_{%s}\"" % self.splitTerm(self.base)

        s = []
        s.append("#set terminal postscript eps enhanced solid" \
                 + "\"Helvetica\" 14.4\n")
        s.append("set terminal postscript eps enhanced color solid " \
                 + "\"Helvetica\" 14.4\n")
        s.append("set encoding iso_8859_1\n")
        s.append("set xlabel \"wavenumber in cm^{-1}\"\n")
        s.append("set ylabel \"absorption cross section in pm^2\"\n")
        s.append("\n")
        s.append("set format x '%.0f'\n")
        s.append("set format y '%.1f'\n")
        s.append("set xrange [0:50000]\n")
        s.append("set yrange [0:2.5]\n")
        s.append("ymax = 2.5\n")
        s.append("\n")
        s.append("set nokey\n")
        s.append("#set mxtics 2\n")
        s.append("#set mytics 2\n")
        s.append("set grid mxtics\n")
        s.append("set grid xtics\n")
        s.append("set grid mytics\n")
        s.append("set grid ytics\n")
        s.append("set grid\n")
        s.append("\n")
        s.append("set label %s at graph 0.02,0.96 \\\n" % base)
        s.append("    left font \"Helvetica,12\"\n")
        s.append("\n")
        return s

    def plotData(self, fn):
        s = []
        s += self.plotHead()

        x,y,dy = self.sigma
        for i in range(1,len(self.energy)):
            val = self.energy[i]
            if val > self.xmax:
                continue
            wave = 1.0e7 / val
            if wave < 1000:
                wave = "\"%inm\"" % round(wave, 0)
            else:
                wave = "\"%gµm\"" % round(wave/1000.0, 2)
            term = "\"^{%s}%s_{%s}\"" % self.splitTerm(self.names[i])
            i = Numeric.argmin(abs(x-val))
            sigma = 0.0
            for j in range(i-1,i+2):
                sigma += y[j]
            sigma /= 3
            line = "call \"peak.inc\" %5i %5.3f \"+0  \" %-14s %-8s ymax\n" \
                 % (val, sigma, term, wave)
            s.append(line)
        s.append("\n")
        s.append("plot \"%s-sigma.data\" using 1:2 notitle with lines 3\n" \
                 % self.element)
        fp = open(fn, "w")
        fp.writelines(s)
        fp.close()

    def splitTerm(self, term):
        vals = re.split(PATT_SLJ, term)
        if len(vals) > 1:
            vals = (vals[1], vals[2], vals[4])
        else:
            vals = re.split(PATT_SL, term)
            if len(vals) > 1:
                vals = (vals[1], vals[2], "")
            else:
                vals = ("", vals[0], "")
        return vals

