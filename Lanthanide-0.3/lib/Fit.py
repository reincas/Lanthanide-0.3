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

class Fit:
    # ----------------------------------------------------------------------
    def splitMeas(self, meas, all=0):
        span   = []
        guess  = []
        y      = []
        sigma  = []
        for i in range(len(meas)):
            level = meas[i][0]
            name  = meas[i][1]
            value = meas[i][2]
            delta = meas[i][3]
            if all or ((delta > 0) and \
               ((type(level) == types.TupleType) or (level >= 0))):
                span.append(level)
                guess.append(name)
                y.append(value)
                sigma.append(delta)
	y      = Numeric.array(y, Numeric.Float)
        sigma  = Numeric.array(sigma, Numeric.Float)
	return span, guess, y, sigma

    # ----------------------------------------------------------------------
    def updateMeas(self, meas, y, all=0):
        index = 0
        for i in range(len(meas)):
            level = meas[i][0]
            delta = meas[i][3]
            if all or ((delta > 0) and \
               ((type(level) == types.TupleType) or (level >= 0))):
                meas[i][2] = y[index]
                index += 1
	return meas

    # ----------------------------------------------------------------------
    def invSigma(self, sigma):
	newsigma = Numeric.array(sigma)
	for i in range(len(sigma)):
	    if sigma[i] != 0:
		newsigma[i] = 1.0 / sigma[i]
	return newsigma

    # ----------------------------------------------------------------------
    def ySum(self, calc, guess, mult=[]):
	shape = calc.shape
	if len(shape) != 2:
	    calc = calc[:, Numeric.NewAxis]
	nummeas = len(self.guess)
        sum = Numeric.zeros((nummeas, calc.shape[1]), Numeric.Float)
        for i in range(nummeas):
	    sum[i] = self.yCalc(calc, i, guess, mult)
	if len(shape) != 2:
	    sum = Numeric.reshape(sum, (sum.shape[0],))
        return sum

    # ----------------------------------------------------------------------
    def yCalc(self, calc, i, guess, mult=[]):
	name = self.guess[i]
	span = self.span[i]
	if type(span) != types.TupleType:
	    sum = calc[self.yIndex(name, span, guess),:]
	else:
	    sumy    = 0.0
	    summult = 0
	    for j in range(len(span)):
		n = name[j]
		k = span[j]
		if mult:
		    k = self.yIndex(n, k, guess)
		    sumy    += calc[k,:] * mult[k]
		    summult += mult[k]
		else:
		    k = self.yIndex(n, k, guess)
		    sumy    += calc[k,:]
	    if mult:
		sum = sumy / summult
	    else:
		sum = sumy
        return sum

    # ----------------------------------------------------------------------
    def yIndex(self, name, i, guess):
        if not name:
            return i
        if name == guess[i]:
            return i
	num = guess.count(name)
	if num == 1:
	    j = guess.index(name)
	    return j
        #print guess
        #raise RuntimeError, "gety: Unknown level %s, %i!" % (name, i)
	return i


