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

class OscPlot:
    def __init__(self, meas, calc):
        self.listMeas = meas
        self.listCalc = calc

    def plotHead(self, x, y=0):
        s = []
        s.append("#set terminal postscript eps enhanced " \
                 + "\"Helvetica\" 14.4\n")
        s.append("set terminal postscript eps enhanced color solid " \
                 + "\"Helvetica\" 14.4\n")
        s.append("set encoding iso_8859_1\n")
        s.append("set xlabel \"wavenumber in cm^{-1}\"\n")
        s.append("set ylabel \"{/Helvetica-Oblique " \
                 + "(f_{calc} - f_{meas}) / f_{meas}}\"\n")
        s.append("set format x '%.0f'\n")
        s.append("set format y '%.1f'\n")
        s.append("set xrange [0:%f]\n" % x)
        if y:
            s.append("set yrange [%f:%f]\n" % (-y, y))
        s.append("set nokey\n")
        s.append("set mxtics 1\n")
        s.append("set mytics 2\n")
        s.append("set grid mxtics\n")
        s.append("set grid xtics\n")
        s.append("set grid mytics\n")
        s.append("set grid ytics\n")
        s.append("set grid\n")
        s.append("set style line 1 lt 3 pt 6\n")
        return s

    def plotMeas(self):
        s = []
        for data in self.listMeas:
            s.append("%f %f %f %f\n" % (data[0], 0, data[1], data[4]/data[3]))
        s.append("e\n")
        return s
    
    def plotFit(self):
        s = []
        for data in self.listMeas:
            s.append("%f %f\n" % (data[2], data[5]/data[3]-1))
        s.append("e\n")
        return s
    
    def plotCalc(self):
        s = []
        for data in self.listCalc:
            s.append("%f %f\n" % (data[0], 0))
        s.append("e\n")
        return s

    def plotData(self, fn, xmax, ymax):
        s = []
        s += self.plotHead(xmax, ymax)
        s.append("plot '-' with xyerrorbars, " +
                 "'-' with points 7, " +
                 "'-' with points ls 1, " +
                 "0\n")
        s += self.plotMeas()
        s += self.plotFit()
        s += self.plotCalc()
        fp = open(fn, "w")
        fp.writelines(s)
        fp.close()
