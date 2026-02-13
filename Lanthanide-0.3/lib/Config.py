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


def buildConfig(l, electrons):
    quant  = buildquant(l)
    unit   = buildunit(quant)
    config = buildconfig(quant, electrons)
    config["quant"] = quant
    config["unit"] = unit
    return config



def Config(l, electrons, rw=IO_READ+IO_WRITE):
    config = {"l": l, "electrons": electrons}
    fn = dataPath(config, "", "config")
    if (rw & IO_READ) and os.path.exists(fn):
        fp = gzip.open(fn, "r")
        config = cPickle.loads(fp.read())
        fp.close()
    else:
        print "Config: calc %s" % fn
	config = buildConfig(l, electrons)
        if rw & IO_WRITE:
            fp = gzip.open(fn, "w")
            fp.write(cPickle.dumps(config))
            fp.close()
    return config
