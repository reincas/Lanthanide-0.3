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

CONFIG = {
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
    "p2":   (1, 2),
    }

class LevelFit:
    def __init__(self, DB, meas):
	self.DB = DB
        self.l,self.electrons = CONFIG[DB["config"]]
	self.config = Config.Config(self.l, self.electrons)
	self.term = Term.Term(self.config)
	self.key = DB["key"]
	self.meas = meas
	self.numstep = len(DB["steps"])
	
	self.parmlist = []
	for parm in DB["parms"]:
	    self.parmlist.append(parm[0])
	self.H = Hamiltonian(self.config, self.term, self.key, self.parmlist)
	    
	self.fit = []
	for i in range(self.numstep):
	    step = DB["steps"][i]
 	    fit = xxx

    def getparm(self, step="a"):
	pass
#... to be continued ...
