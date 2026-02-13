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

import sys
import types
import copy
import os
import math
import re
import string
import types
import cPickle
import gzip
import Numeric
import LinearAlgebra

from lanthanide import *
from BaseParms import *

TERM_PRODUCT = 0
TERM_SLJM    = 1
TERM_SLJ     = 2
TERM_SL      = 3

MATRIX_UNIT  = 0
MATRIX_H4    = 1
MATRIX_H5    = 2
MATRIX_H6    = 3

DATADIR  = "matrix"
IO_NONE  = 0
IO_READ  = 1
IO_WRITE = 2

FIT_SINGLE = 0
FIT_DOUBLE = 1

CONST_L    = 6.0220e23     # 1 / mol
CONST_k    = 1.3807e-23    # J / K
CONST_eps0 = 8.8542e-12    # C / V m
CONST_me   = 9.1095e-31    # kg
CONST_e    = 1.6022e-19    # C
CONST_h    = 6.6262e-34    # J s
CONST_c    = 2.99792458e8  # m / s
