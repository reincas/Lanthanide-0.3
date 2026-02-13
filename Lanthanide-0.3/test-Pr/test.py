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
from Ion import Ion

from Cauchy import Cauchy
from PrMeas import Energy, OscStrength

config   = "f02"
element  = "Pr"
meas     = Energy
fmeas    = OscStrength
cauchy   = Cauchy["ZBLAN"]

ion = Ion(config, 1)

parmFit(ion, element, "",  "1", FIT_SINGLE, meas, fmeas, cauchy)
parmFit(ion, element, "1", "2a", FIT_SINGLE, meas, fmeas, cauchy)
parmFit(ion, element, "2a", "3", FIT_SINGLE, meas, fmeas, cauchy)

#xmax = 25000
#ymax = 1
#ion.plotOsc(element + "-fit.gnu", xmax, ymax, meas, fmeas)
#ion.doubleLaTeX(element + "-best.tab", meas, fmeas)
#ion.vectorLaTeX(element + "-intermediate.tab", meas, fmeas)
#ion.emissionLaTeX(element + "-emission.tab", meas, fmeas, 0)
#ion.absorptionLaTeX(element + "-absorption.tab", meas, fmeas, 0)

#print ion.confShortNames("slj")
