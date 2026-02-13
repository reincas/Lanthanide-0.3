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



//--------------------------------------------------------------------
long QuantToPy(quantStruct *quant, PyObject **dict)
{
  long i;
  int dims[1];
  double *data;
  PyArrayObject *ml, *ms;

  //printf("QuantToPy: build ml array\n");
  dims[0] = quant->max;
  if (!(ml = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_DOUBLE)))
    return -1;
  data = (double *)ml->data;
  for (i=0; i<dims[0]; i++)
    data[i] = 0.5*quant->ml[i];

  //printf("QuantToPy: build ms array\n");
  dims[0] = quant->max;
  if (!(ms = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_DOUBLE)))
    return -2;
  data = (double *)ms->data;
  for (i=0; i<dims[0]; i++)
    data[i] = 0.5*quant->ms[i];

  //printf("QuantToPy: build dictionary\n");
  if (!(*dict = Py_BuildValue("{s:i,s:i,s:i,s:O,s:O}",
			      "DictType",  DICT_QUANT,
			      "l",         quant->l / 2,
			      "max",       quant->max,
			      "ml",        ml,
			      "ms",        ms)))
    return -3;
  Py_DECREF(ml);
  Py_DECREF(ms);

  //printf("QuantToPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long QuantFromPy(PyObject *dict, quantStruct *quant)
{
  long i;
  double *data;
  PyArrayObject *array;

  //printf("QuantFromPy: set values\n");
  quant->l      = PyInt_AsLong(PyDict_GetItemString(dict,"l")) * 2;
  quant->max    = PyInt_AsLong(PyDict_GetItemString(dict,"max"));
  quant->ownbuf = 1;

  //printf("QuantFromPy: allocate memory for structures\n");
  if (!(quant->ml = malloc(quant->max * sizeof(long))))
    return -1;
  if (!(quant->ms = malloc(quant->max * sizeof(long))))
    return -2;

  //printf("QuantFromPy: fill ml\n");
  array = (PyArrayObject *)PyDict_GetItemString(dict, "ml");
  data = (double *)array->data;
  for (i=0; i<quant->max; i++)
    quant->ml[i] = DTO2L(data[i]);
  
  //printf("QuantFromPy: fill ms\n");
  array = (PyArrayObject *)PyDict_GetItemString(dict, "ms");
  data = (double *)array->data;
  for (i=0; i<quant->max; i++)
    quant->ms[i] = DTO2L(data[i]);
  
  //printf("QuantFromPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long UnitToPy(unitStruct *unit, PyObject **dict)
{
  int dims[3];
  PyArrayObject *u, *t;

  //printf("UnitToPy: build u array\n");
  dims[0] = unit->l+1;
  dims[1] = unit->max;
  dims[2] = unit->max;
  if (!(u = (PyArrayObject *) PyArray_FromDims(3, dims, PyArray_DOUBLE)))
    return -1;
  free(u->data);
  u-> data = (char *)unit->u;

  //printf("UnitToPy: build t array\n");
  dims[0] = 2;
  dims[1] = unit->max;
  dims[2] = unit->max;
  if (!(t = (PyArrayObject *)PyArray_FromDims(3, dims, PyArray_DOUBLE)))
    return -2;
  free(t->data);
  t->data = (char *)unit->t;

  //printf("UnitToPy: build dictionary\n");
  if (!(*dict = Py_BuildValue("{s:i,s:i,s:i,s:O,s:O}",
			      "DictType",  DICT_UNIT,
			      "l",         unit->l / 2,
			      "max",       unit->max,
			      "u",         u,
			      "t",         t)))
    return -3;
  Py_DECREF(u);
  Py_DECREF(t);
  unit->ownbuf = 0;

  //printf("UnitToPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long UnitFromPy(PyObject *dict, unitStruct *unit)
{
  PyArrayObject *u, *t;

  //printf("UnitFromPy: get u and t\n");
  u = (PyArrayObject *)PyDict_GetItemString(dict, "u");
  t = (PyArrayObject *)PyDict_GetItemString(dict, "t");

  //printf("UnitFromPy: set values\n");
  unit->l      = PyInt_AsLong(PyDict_GetItemString(dict,"l")) * 2;
  unit->max    = PyInt_AsLong(PyDict_GetItemString(dict,"max"));
  unit->u      = (double *)u->data;
  unit->t      = (double *)t->data;
  unit->ownbuf = 0;
  
  //printf("UnitFromPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long ConfigToPy(configStruct *config, PyObject **dict)
{
  int dims[2];
  PyArrayObject *state;

  //printf("ConfigToPy: build states array\n");
  dims[0] = config->states;
  dims[1] = config->electrons;
  if (!(state = (PyArrayObject *) PyArray_FromDims(2, dims, PyArray_LONG)))
    return -1;
  free(state->data);
  state->data = (char *)config->state;

  //printf("ConfigToPy: build dictionary\n");
  if (!(*dict = Py_BuildValue("{s:i,s:i,s:i,s:i,s:O}",
			      "DictType",  DICT_CONFIG,
			      "l",         config->l / 2,
			      "electrons", config->electrons,
			      "states",    config->states,
			      "state",     state)))
    return -2;
  Py_DECREF(state);
  config->ownbuf = 0;

  //printf("ConfigToPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long ConfigFromPy(PyObject *dict, configStruct *config)
{
  PyArrayObject *state;

  //printf("ConfigFromPy: set state\n");
  state = (PyArrayObject *)PyDict_GetItemString(dict, "state");

  //printf("ConfigFromPy: set variables\n");
  config->l         = PyInt_AsLong(PyDict_GetItemString(dict,"l")) * 2;
  config->electrons = PyInt_AsLong(PyDict_GetItemString(dict,"electrons"));
  config->states    = PyInt_AsLong(PyDict_GetItemString(dict,"states"));
  config->state     = (long *)state->data;
  config->ownbuf    = 0;

  //printf("ConfigFromPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long PairToPy(pairStruct *pair, PyObject **dict)
{
  int dims[1];
  PyArrayObject *statei, *statej;

  //printf("PairToPy: build statei\n");
  dims[0] = pair->electrons;
  if (!(statei = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_LONG)))
    return -1;
  free(statei->data);
  statei->data = (char *)pair->statei;

  //printf("PairToPy: build statej\n");
  dims[0] = pair->electrons;
  if (!(statej = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_LONG)))
    return -2;
  free(statej->data);
  statej->data = (char *)pair->statej;

  //printf("PairToPy: build dictionary\n");
  if (!(*dict = Py_BuildValue("{s:i,s:i,s:i,s:i,s:O,s:O,s:i,s:i}",
		      "DictType",  DICT_PAIR,
		      "i",         pair->i,
		      "j",         pair->j,
		      "electrons", pair->electrons,
		      "statei",    statei,
		      "statej",    statej,
		      "unpaired",  pair->unpaired,
		      "sign",      pair->sign)))
    return -3;
  Py_DECREF(statei);
  Py_DECREF(statej);
  pair->ownbuf = 0;

  //printf("PairToPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long PairFromPy(PyObject *dict, pairStruct *pair)
{
  PyArrayObject *statei, *statej;

  //printf("PairFromPy: get statei and statej\n");
  statei = (PyArrayObject *)PyDict_GetItemString(dict, "statei");
  statej = (PyArrayObject *)PyDict_GetItemString(dict, "statej");

  //printf("PairFromPy: set variables\n");
  pair->i         = PyInt_AsLong(PyDict_GetItemString(dict,"i"));
  pair->j         = PyInt_AsLong(PyDict_GetItemString(dict,"j"));
  pair->electrons = PyInt_AsLong(PyDict_GetItemString(dict,"electrons"));
  pair->statei    = (long *)statei->data;
  pair->statej    = (long *)statej->data;
  pair->unpaired  = PyInt_AsLong(PyDict_GetItemString(dict,"unpaired"));
  pair->sign      = PyInt_AsLong(PyDict_GetItemString(dict,"sign"));
  pair->ownbuf    = 0;
  
  //printf("PairFromPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long SparseToPy(sparseStruct *matrix, PyObject **dict)
{
  long total, size;
  int dims[1], type;
  PyArrayObject *data, *col, *row;

  //printf("SparseToPy: set some values\n");
  total = matrix->row[matrix->states-1];
  size = matrix->size;
  if (size == 1)
    type = PyArray_DOUBLE;
  else
    type = PyArray_CDOUBLE;

  //printf("SparseToPy: build data\n");
  dims[0] = total;
  if(!(data = (PyArrayObject *)PyArray_FromDims(1, dims, type)))
    return -1;
  free(data->data);
  data->data = (char *)matrix->data;
  
  //printf("SparseToPy: build col\n");
  dims[0] = total;
  if (!(col = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_LONG)))
    return -2;
  free(col->data);
  col->data = (char *)matrix->col;
  
  //printf("SparseToPy: build row\n");
  dims[0] = matrix->states;
  if (!(row = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_LONG)))
    return -3;
  free(row->data);
  row->data = (char *)matrix->row;

  //printf("SparseToPy: build dictionary\n");
  if (!(*dict = Py_BuildValue("{s:i,s:i,s:i,s:O,s:O,s:O}",
			      "DictType", DICT_SPARSE,
			      "states",   matrix->states,
			      "size",     matrix->size,
			      "data",     data,
			      "col",      col,
			      "row",      row)))
    return -4;
  Py_DECREF(data);
  Py_DECREF(col);
  Py_DECREF(row);
  matrix->ownbuf = 0;

  //printf("SparseToPy: return\n");
  return 0;
}



//--------------------------------------------------------------------
long SparseFromPy(PyObject *dict, sparseStruct *matrix)
{
  PyArrayObject *data, *col, *row;
  
  //printf("SparseFromPy: get data, col, and row\n");
  data = (PyArrayObject *)PyDict_GetItemString(dict, "data");
  col  = (PyArrayObject *)PyDict_GetItemString(dict, "col");
  row  = (PyArrayObject *)PyDict_GetItemString(dict, "row");

  //printf("SparseFromPy: fill variables\n");
  matrix->states = PyInt_AsLong(PyDict_GetItemString(dict,"states"));
  matrix->size   = PyInt_AsLong(PyDict_GetItemString(dict,"size"));
  matrix->data   = (double *)data->data;
  matrix->col    = (long *)col->data;
  matrix->row    = (long *)row->data;
  matrix->ownbuf = 0;

  //printf("SparseFromPy: return\n");
  return 0;
}
