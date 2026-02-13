#!/usr/bin/env python
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

import os
from distutils.core import setup, Extension

HOME = os.environ['HOME']

setup (name = 'Lanthanide',
       version = '0.3',
       description = 'Routines for energy level calculations',
       long_description = \
  'Routines for energy level calculations: determinantal product states, \n'+\
  'LS-coupling, intermediate coupling, all common known energy operators \n'+\
  'for f-shell configurations, many operators for p-shell configurations, \n'+\
  'Judd-Ofelt calculations. \n',
       author = 'Reinhard Caspary',
       author_email = 'r.caspary@tu-bs.de',
       url = 'http://www.tu-braunschweig.de/ihf',
       license = 'GPL',

       packages = [''],
       package_dir = {'': 'lib'},
       extra_path = HOME+'/lib/python/Lanthanide',
       include_dirs = ['include'],
       ext_modules = [Extension('lanthanide',
				['src/py-module.c',
				 'src/py-types.c',
				 'src/py-calc.c',
				 'src/py-diagonalize.c',
				 'src/py-wigner.c',
				 'src/calc.c',
				 'src/matrix-one.c',
				 'src/matrix-two.c',
				 'src/matrix-three.c',
				 'src/diagonalize.c',
				 'src/wigner.c'],
				libraries=['lapack', 'blas', 'g2c', 'm'],
				#libraries=['blas', 'g2c', 'm'],
                                #extra_objects=['/usr/lib64/liblapack.so.3',
                                #               '/usr/lib64/libblas.so',
                                #               '/usr/lib64/libg2c.so',
                                #               '/usr/lib64/libm.so'],
                                extra_compile_args = ['-fPIC'],
                                extra_link_args = ['-fPIC'],
				)
		      ],
    )

