// Copyright 2000,2003 Reinhard Caspary <r.caspary@tu-bs.de>
//
// This file is part of the Python package Lanthanide.
// 
// Lanthanide is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Lanthanide is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Lanthanide; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_lanthanide
#include <Numeric/arrayobject.h>
#include "lanthanide.h"



//--------------------------------------------------------------------
// Exported Python functions from py-calc.c:
//--------------------------------------------------------------------

// def binom(n, k): value
PyObject *lanthanide_binom(PyObject *self, PyObject *args);

// def buildquant(l): quant
PyObject *lanthanide_buildquant(PyObject *self, PyObject *args);

// def buildunit(quant): unit
PyObject *lanthanide_buildunit(PyObject *self, PyObject *args);

// def buildconfig(quant, electrons): config
PyObject *lanthanide_buildconfig(PyObject *self, PyObject *args);

// def orderpair(config, i, j): pair
PyObject *lanthanide_orderpair(PyObject *self, PyObject *args);

// def calcmatrix(quant, unit, config, key, opts): matrix
PyObject *lanthanide_calcmatrix(PyObject *self, PyObject *args);

// def calcelement(quant, unit, pair, key, opts): value
PyObject *lanthanide_calcelement(PyObject *self, PyObject *args);

// def getelement(matrix, i, j): value
PyObject *lanthanide_getelement(PyObject *self, PyObject *args);

// def multadd(sparse, factor, tri)
PyObject *lanthanide_multadd(PyObject *self, PyObject *args);

// def transform(sparse, vecs): sparse
PyObject *lanthanide_transform(PyObject *self, PyObject *args);

// def transdiag(sparse, vecs): vector
PyObject *lanthanide_transdiag(PyObject *self, PyObject *args);

// def submatrix(tri, substart, subend): tri
PyObject *lanthanide_submatrix(PyObject *self, PyObject *args);



//--------------------------------------------------------------------
// Exported Python functions from py-diagonalize.c:
//--------------------------------------------------------------------

// def diagonalize(method, vecflag, matrix): (vals, vecs, info)
PyObject *lanthanide_diagonalize(PyObject *self, PyObject *args);



//--------------------------------------------------------------------
// Exported Python functions from py-wigner.c:
//--------------------------------------------------------------------

// def wign3j(a, b, c, x, y, z): value
PyObject *lanthanide_wign3j(PyObject *self, PyObject *args);

// def clebsg(a, b, c, x, y, z): value
PyObject *lanthanide_clebsg(PyObject *self, PyObject *args);

// def wign6j(a, b, c, x, y, z): value
PyObject *lanthanide_wign6j(PyObject *self, PyObject *args);

// def racahc(a, b, c, x, y, z): value
PyObject *lanthanide_racahc(PyObject *self, PyObject *args);

// def jahnuf(a, b, c, x, y, z): value
PyObject *lanthanide_jahnuf(PyObject *self, PyObject *args);

// def wign9j(a, b, c, x, y, z, g, h, p): value
PyObject *lanthanide_wign9j(PyObject *self, PyObject *args);



//--------------------------------------------------------------------
// Python dictionary representations of C data stuctures.
//--------------------------------------------------------------------

// Type numbers for the dictionary representations.
#define DICT_QUANT   0    // quantStruct
#define DICT_UNIT    1    // unitStruct
#define DICT_CONFIG  2    // configStruct
#define DICT_PAIR    3    // pairStruct
#define DICT_SPARSE  4    // sparseStruct
#define DICT_SMATRIX 5    // sparse full matrix


// functions to handle the data type quantStruct
long QuantToPy(quantStruct *quant, PyObject **dict);
long QuantFromPy(PyObject *dict, quantStruct *quant);


// functions to handle the data type unitStruct
long UnitToPy(unitStruct *unit, PyObject **dict);
long UnitFromPy(PyObject *dict, unitStruct *unit);


// functions to handle the data type configStruct
long ConfigToPy(configStruct *config, PyObject **dict);
long ConfigFromPy(PyObject *dict, configStruct *config);


// functions to handle the data type pairStruct
long PairToPy(pairStruct *pair, PyObject **dict);
long PairFromPy(PyObject *dict, pairStruct *pair);


// functions to handle the data type sparseStruct
long SparseToPy(sparseStruct *matrix, PyObject **dict);
long SparseFromPy(PyObject *dict, sparseStruct *matrix);



//--------------------------------------------------------------------
// Initalisation function for the module lanthanide.
//--------------------------------------------------------------------

void initlanthanide();
