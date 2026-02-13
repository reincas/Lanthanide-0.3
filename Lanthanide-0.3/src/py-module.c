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

#include "py-lanthanide.h"



static PyMethodDef LanthanideMethods[] = {
  {"buildquant",    lanthanide_buildquant,    METH_VARARGS},
  {"buildunit",     lanthanide_buildunit,     METH_VARARGS},
  {"buildconfig",   lanthanide_buildconfig,   METH_VARARGS},
  {"orderpair",     lanthanide_orderpair,     METH_VARARGS},
  {"calcmatrix",    lanthanide_calcmatrix,    METH_VARARGS},
  {"calcelement",   lanthanide_calcelement,   METH_VARARGS},
  {"getelement",    lanthanide_getelement,    METH_VARARGS},
  {"multadd",       lanthanide_multadd,       METH_VARARGS},
  {"transform",     lanthanide_transform,     METH_VARARGS},
  {"transdiag",     lanthanide_transdiag,     METH_VARARGS},
  {"submatrix",     lanthanide_submatrix,     METH_VARARGS},
  {"diagonalize",   lanthanide_diagonalize,   METH_VARARGS},
  {"binom",         lanthanide_binom,         METH_VARARGS},
  {"wign3j",        lanthanide_wign3j,        METH_VARARGS},
  {"clebsg",        lanthanide_clebsg,        METH_VARARGS},
  {"wign6j",        lanthanide_wign6j,        METH_VARARGS},
  {"racahc",        lanthanide_racahc,        METH_VARARGS},
  {"jahnuf",        lanthanide_jahnuf,        METH_VARARGS},
  {"wign9j",        lanthanide_wign9j,        METH_VARARGS},
  {NULL,      NULL}        /* Sentinel */
};



//--------------------------------------------------------------------
void initlanthanide()
{
  PyObject *module, *dict;

  module = Py_InitModule("lanthanide", LanthanideMethods);
  import_array();

  initwigner();

  // indices of contributions to the hamiltonian:
  dict = PyModule_GetDict(module);

// Matrix element keys:
  PyDict_SetItemString(dict, "OPE",          PyInt_FromLong(MAT1_OPE));
  PyDict_SetItemString(dict, "MAT1_U",       PyInt_FromLong(MAT1_U));
  PyDict_SetItemString(dict, "MAT1_T",       PyInt_FromLong(MAT1_T));
  PyDict_SetItemString(dict, "MAT1_UU",      PyInt_FromLong(MAT1_UU));
  PyDict_SetItemString(dict, "MAT1_UUS",     PyInt_FromLong(MAT1_UUS));
  PyDict_SetItemString(dict, "MAT1_TT",      PyInt_FromLong(MAT1_TT));
  PyDict_SetItemString(dict, "MAT1_UT",      PyInt_FromLong(MAT1_UT));
  PyDict_SetItemString(dict, "MAT2_UU",      PyInt_FromLong(MAT2_UU));
  PyDict_SetItemString(dict, "MAT2_UUS",     PyInt_FromLong(MAT2_UUS));
  PyDict_SetItemString(dict, "MAT2_TT",      PyInt_FromLong(MAT2_TT));
  PyDict_SetItemString(dict, "MAT2_UT",      PyInt_FromLong(MAT2_UT));
  PyDict_SetItemString(dict, "MAT2_UUTT",    PyInt_FromLong(MAT2_UUTT));
  PyDict_SetItemString(dict, "MAT3_UUUA",    PyInt_FromLong(MAT3_UUUA));
  PyDict_SetItemString(dict, "MAT3_UUUB",    PyInt_FromLong(MAT3_UUUB));
  PyDict_SetItemString(dict, "MAT3_UUUC",    PyInt_FromLong(MAT3_UUUC));

  // export MEMFACTOR:
  PyDict_SetItemString(dict, "MEMFACTOR",    PyInt_FromLong(MEMFACTOR));

  // export options for the diagonalisation function:
  PyDict_SetItemString(dict, "DIAG_FAST",    PyInt_FromLong(DIAG_FAST));
  PyDict_SetItemString(dict, "DIAG_SMALL",   PyInt_FromLong(DIAG_SMALL));
  PyDict_SetItemString(dict, "DIAG_NOVEC",   PyInt_FromLong(DIAG_NOVEC));
  PyDict_SetItemString(dict, "DIAG_VEC",     PyInt_FromLong(DIAG_VEC));

  // export dictionary types:
  PyDict_SetItemString(dict, "DICT_QUANT",   PyInt_FromLong(DICT_QUANT));
  PyDict_SetItemString(dict, "DICT_UNIT",    PyInt_FromLong(DICT_UNIT));
  PyDict_SetItemString(dict, "DICT_CONFIG",  PyInt_FromLong(DICT_CONFIG));
  PyDict_SetItemString(dict, "DICT_PAIR",    PyInt_FromLong(DICT_PAIR));
  PyDict_SetItemString(dict, "DICT_SPARSE",  PyInt_FromLong(DICT_SPARSE));
  PyDict_SetItemString(dict, "DICT_SMATRIX", PyInt_FromLong(DICT_SMATRIX));
}
