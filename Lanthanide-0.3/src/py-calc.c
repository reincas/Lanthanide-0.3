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
PyObject *lanthanide_binom(PyObject *self, PyObject *args)
{
  long n, k;
  
  //printf("lanthanide_binom: parse args\n");
  if (!PyArg_ParseTuple(args, "ll", &n, &k))
    {
      PyErr_SetString(PyExc_ValueError, "binom: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_binom: calc binomial coefficient\n");
  return(PyLong_FromLong(binom(n,k)));
}



//--------------------------------------------------------------------
PyObject *lanthanide_buildquant(PyObject *self, PyObject *args)
{
  long l;
  quantStruct quant;
  PyObject *dict;

  //printf("lanthanide_buildquant: parse args\n");
  if (!PyArg_ParseTuple(args, "l", &l))
    {
      PyErr_SetString(PyExc_ValueError, "buildquant: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_buildquant: build\n");
  if (buildquant(2*l, &quant))
    {
      PyErr_SetString(PyExc_ValueError, "buildquant: build error");
      return NULL;
    }

  //printf("lanthanide_buildquant: convert quant\n");
  if (QuantToPy(&quant, &dict))
    {
      PyErr_SetString(PyExc_ValueError, "buildquant: convert quant error");
      return NULL;
    }
  
  //printf("lanthanide_buildquant: return\n");
  freeQuant(&quant);
  return dict;
}



//--------------------------------------------------------------------
PyObject *lanthanide_buildunit(PyObject *self, PyObject *args)
{
  quantStruct quant;
  unitStruct unit;
  PyObject *qdict, *udict;

  //printf("lanthanide_buildunit: parse args\n");
  if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &qdict))
    {
      PyErr_SetString(PyExc_ValueError, "buildunit: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_buildunit: convert quant\n");
  if (QuantFromPy(qdict, &quant))
    {
      PyErr_SetString(PyExc_ValueError, "buildunit: convert quant error");
      return NULL;
    }
      
  //printf("lanthanide_buildunit: build\n");
  if (buildunit(&quant, &unit))
    {
      PyErr_SetString(PyExc_ValueError, "buildunit: build error");
      return NULL;
    }

  //printf("lanthanide_buildunit: convert unit\n");
  if (UnitToPy(&unit, &udict))
    {
      PyErr_SetString(PyExc_ValueError, "buildunit: convert unit error");
      return NULL;
    }

  //printf("lanthanide_buildunit: return\n");
  freeQuant(&quant);
  freeUnit(&unit);
  return udict;
}



//--------------------------------------------------------------------
PyObject *lanthanide_buildconfig(PyObject *self, PyObject *args)
{
  long electrons;
  quantStruct quant;
  configStruct config;
  PyObject *qdict, *cdict;

  //printf("lanthanide_buildconfig: parse args\n");
  if (!PyArg_ParseTuple(args, "O!l", &PyDict_Type, &qdict, &electrons))
    {
      PyErr_SetString(PyExc_ValueError, "buildconfig: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_buildconfig: convert quant\n");
  if (QuantFromPy(qdict, &quant))
    {
      PyErr_SetString(PyExc_ValueError, "buildconfig: convert quant error");
      return NULL;
    }
      
  //printf("lanthanide_buildconfig: build\n");
  if (buildconfig(&quant, electrons, &config))
    {
      PyErr_SetString(PyExc_ValueError, "buildconfig: build error");
      return NULL;
    }

  //printf("lanthanide_buildconfig: convert config\n");
  if (ConfigToPy(&config, &cdict))
    {
      PyErr_SetString(PyExc_ValueError, "buildconfig: convert config error");
      return NULL;
    }
  
  //printf("lanthanide_buildconfig: return\n");
  freeQuant(&quant);
  freeConfig(&config);
  return cdict;
}



//--------------------------------------------------------------------
PyObject *lanthanide_orderpair(PyObject *self, PyObject *args)
{
  long i, j;
  configStruct config;
  pairStruct pair;
  PyObject *cdict, *pdict;

  //printf("lanthanide_orderpair: parse args\n");
  if (!PyArg_ParseTuple(args, "O!ll", &PyDict_Type, &cdict, &i, &j))
    {
      PyErr_SetString(PyExc_ValueError, "orderpair: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_orderpair: convert config\n");
  if (ConfigFromPy(cdict, &config))
    {
      PyErr_SetString(PyExc_ValueError, "orderpair: convert config error");
      return NULL;
    }
      
  //printf("lanthanide_orderpair: init pair\n");
  if (initpair(&config, &pair))
    {
      PyErr_SetString(PyExc_ValueError, "orderpair: init pair error");
      return NULL;
    }

  //printf("lanthanide_orderpair: order pair\n");
  if (orderpair(&config, i, j, &pair))
    {
      PyErr_SetString(PyExc_ValueError, "orderpair: order error");
      return NULL;
    }

  //printf("lanthanide_orderpair: convert pair\n");
  if (PairToPy(&pair, &pdict))
    {
      PyErr_SetString(PyExc_ValueError, "orderpair: convert pair error");
      return NULL;
    }
  
  //printf("lanthanide_orderpair: return\n");
  freeConfig(&config);
  freePair(&pair);
  return pdict;
}



//--------------------------------------------------------------------
PyObject *lanthanide_calcmatrix(PyObject *self, PyObject *args)
{
  long i;
  long key, numopts, *opts;
  quantStruct quant;
  unitStruct unit;
  configStruct config;
  sparseStruct matrix;
  PyObject *qdict, *udict, *cdict, *mdict, *list;

  //printf("lanthanide_calcmatrix: parse first part of args\n");
  if (!PyArg_ParseTuple(args, "O!O!O!lO!",
			&PyDict_Type, &qdict,
			&PyDict_Type, &udict,
			&PyDict_Type, &cdict,
			&key,
			&PyList_Type, &list))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_calcmatrix: allocate option memory\n");
  numopts = ELEMENT[key].opts;
  if(!(opts = malloc(numopts * sizeof(long))))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: error allocating memory");
      return NULL;
    }

  //printf("lanthanide_calcmatrix: parse option args\n");
  for (i=0; i<numopts; i++)
    opts[i] = 2 * PyInt_AsLong(PyList_GetItem(list, i));

  //printf("lanthanide_calcmatrix: convert quant\n");
  if (QuantFromPy(qdict, &quant))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: convert quant error");
      return NULL;
    }

  //printf("lanthanide_calcmatrix: convert unit\n");
  if (UnitFromPy(udict, &unit))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: convert unit error");
      return NULL;
    }

  //printf("lanthanide_calcmatrix: convert config\n");
  if (ConfigFromPy(cdict, &config))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: convert config error");
      return NULL;
    }

  //printf("lanthanide_calcmatrix: fill matrix\n");
  if (calcmatrix(&quant, &unit, &config, key, opts, &matrix))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: fill error");
      return NULL;
    }

  //printf("lanthanide_calcmatrix: convert matrix\n");
  if (SparseToPy(&matrix, &mdict))
    {
      PyErr_SetString(PyExc_ValueError, "calcmatrix: convert matrix error");
      return NULL;
    }
  
  //printf("lanthanide_calcmatrix: return\n");
  freeConfig(&config);
  freeUnit(&unit);
  freeSparse(&matrix);
  return mdict;
}



//--------------------------------------------------------------------
PyObject *lanthanide_calcelement(PyObject *self, PyObject *args)
{
  long i;
  long key, numopts, *opts;
  quantStruct quant;
  unitStruct unit;
  pairStruct pair;
  double val;
  PyObject *qdict, *udict, *pdict, *value, *list;
  
  //printf("lanthanide_calcelement: parse first part of args\n");
  if (!PyArg_ParseTuple(args, "O!O!O!lO!",
			&PyDict_Type, &qdict,
			&PyDict_Type, &udict,
			&PyDict_Type, &pdict,
			&key,
			&PyList_Type, &list))
    {
      PyErr_SetString(PyExc_ValueError, "calcelement: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_calcelement: allocate option memory\n");
  numopts = ELEMENT[key].opts;
  if(!(opts = malloc(numopts * sizeof(long))))
    {
      PyErr_SetString(PyExc_ValueError,
		      "calcelement: error allocating memory");
      return NULL;
    }

  //printf("lanthanide_calcelement: parse option args\n");
  for (i=0; i<numopts; i++)
    opts[i] = 2 * PyInt_AsLong(PyList_GetItem(list, i));

  //printf("lanthanide_calcelement: convert quant\n");
  if (QuantFromPy(qdict, &quant))
    {
      PyErr_SetString(PyExc_ValueError, "calcelement: convert quant error");
      return NULL;
    }

  //printf("lanthanide_calcelement: convert unit\n");
  if (UnitFromPy(udict, &unit))
    {
      PyErr_SetString(PyExc_ValueError, "calcelement: convert unit error");
      return NULL;
    }

  //printf("lanthanide_calcelement: convert pair\n");
  if (PairFromPy(pdict, &pair))
    {
      PyErr_SetString(PyExc_ValueError, "calcelement: convert pair error");
      return NULL;
    }

  //printf("lanthanide_calcelement: calc element\n");
  if (calcelement(&quant, &unit, &pair, key, opts, &val))
    {
      PyErr_SetString(PyExc_ValueError,
		      "calcelement: error calculating element");
      return NULL;
    }
  value = PyFloat_FromDouble(val);

  //printf("lanthanide_calcelement: return\n");
  freeQuant(&quant);
  freeUnit(&unit);
  freePair(&pair);
  free(opts);
  return value;
}



//--------------------------------------------------------------------
PyObject *lanthanide_getelement(PyObject *self, PyObject *args)
{
  long i, j;
  double rval;
  complex cval;
  sparseStruct matrix;
  PyObject *dict, *value;
  
  //printf("getelement: parse args\n");
  if (!PyArg_ParseTuple(args, "O!ll", &PyDict_Type, &dict, &i, &j))
    {
      PyErr_SetString(PyExc_ValueError, "getelement: wrong arguments");
      return NULL;
    }

  //printf("getelement: convert matrix\n");
  if (SparseFromPy(dict, &matrix))
    {
      PyErr_SetString(PyExc_ValueError, "getelement: convert matrix error");
      return NULL;
    }

  if (matrix.size == 1)
    {
      //printf("getelement: get double\n");
      rval = realgetelement(&matrix, i, j);
      value = PyFloat_FromDouble(rval);
    }
  else
    {
      //printf("getelement: get complex\n");
      cval = cplxgetelement(&matrix, i, j);
      value = PyComplex_FromDoubles(cval.r, cval.i);
    }

  //printf("getelement: return\n");
  freeSparse(&matrix);
  return value;
}



//--------------------------------------------------------------------
PyObject *lanthanide_multadd(PyObject *self, PyObject *args)
{
  int type;
  long size;
  double *ddata;
  complex *cdata;
  complex factor;
  sparseStruct matrix;
  PyObject *dict;
  PyArrayObject *array;
  
  //printf("lanthanide_multadd: parse args\n");
  if (!PyArg_ParseTuple(args, "O!DO!",
			&PyDict_Type, &dict,
			&factor,
			&PyArray_Type, &array))
    {
      PyErr_SetString(PyExc_ValueError, "multadd: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_multadd: set type\n");
  type = array->descr->type_num;
  if ((type != PyArray_CDOUBLE) &&
      (type != PyArray_DOUBLE))
    {
      PyErr_SetString(PyExc_ValueError, "multadd: array has wrong type");
      return NULL;
    }
  if (type == PyArray_DOUBLE)
    size = 1;
  else
    size = 2;

  //printf("lanthanide_multadd: test factor\n");
  if ((type == PyArray_DOUBLE) &&
      (factor.i != 0.0))
    {
      PyErr_SetString(PyExc_ValueError, "multadd: factor is complex");
      return NULL;
    }

  //printf("lanthanide_multadd: convert sparse matrix\n");
  if (SparseFromPy(dict, &matrix))
    {
      PyErr_SetString(PyExc_ValueError, "multadd: convert matrix error");
      return NULL;
    }

  //printf("lanthanide_multadd: test sparse matrix size\n");
  if (matrix.size > size)
    {
      PyErr_SetString(PyExc_ValueError, "multadd: matrix dimension error");
      return NULL;
    }

  //printf("lanthanide_multadd: convert dest matrix\n");
  if ((array->nd != 1) ||
      (array->dimensions[0] != matrix.states * (matrix.states+1) / 2) ||
      (array->flags & (CONTIGUOUS == 0)))
    {
      PyErr_SetString(PyExc_ValueError, "multadd: matrix type error");
      return NULL;
    }
  
  //printf("lanthanide_multadd: calc multadd\n");
  if (type == PyArray_DOUBLE)
    {
      ddata = (double *)array->data;
      realmultadd(&matrix, factor.r, ddata);
    }
  else
    {
      cdata = (complex *)array->data;
      cplxmultadd(&matrix, factor, cdata);
    }

  //printf("lanthanide_multadd: free matrix\n");
  freeSparse(&matrix);

  //printf("lanthanide_multadd: cleanup\n");
  // Return value should be Py_None, but that is defect and causes
  // spurios errors...
  return PyLong_FromLong(0);
}



//--------------------------------------------------------------------
PyObject *lanthanide_transform(PyObject *self, PyObject *args)
{
  int type;
  long size, err, terms;
  double *dtrans;
  complex *ctrans;
  sparseStruct matrix, result;
  PyObject *dict, *mdict;
  PyArrayObject *array;
  
  //printf("lanthanide_transform: parse args\n");
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyDict_Type, &dict,
			&PyArray_Type, &array))
    {
      PyErr_SetString(PyExc_ValueError, "transform: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_transform: set type\n");
  type = array->descr->type_num;
  if ((type != PyArray_CDOUBLE) &&
      (type != PyArray_DOUBLE))
    {
      PyErr_SetString(PyExc_ValueError, "transform: trans has wrong type");
      return NULL;
    }
  if (type == PyArray_DOUBLE)
    size = 1;
  else
    size = 2;

  //printf("lanthanide_transform: convert matrix\n");
  if (SparseFromPy(dict, &matrix))
    {
      PyErr_SetString(PyExc_ValueError, "transform: convert matrix error");
      return NULL;
    }

  //printf("lanthanide_transform: test matrix size\n");
  if (matrix.size > size)
    {
      PyErr_SetString(PyExc_ValueError, "transform: matrix dimension error");
      return NULL;
    }

  //printf("lanthanide_transform: convert trans matrix\n");
  if ((array->nd != 2) ||
      (array->dimensions[1] != matrix.states))
    {
      PyErr_SetString(PyExc_ValueError, "transform: trans type error");
      return NULL;
    }
  terms = array->dimensions[0];

  //printf("lanthanide_transform: make trans matrix contiguous\n");
  array = (PyArrayObject *)
    PyArray_ContiguousFromObject((PyObject *)array, 
				 array->descr->type_num, 
				 array->nd,
				 array->nd);
  if (type == PyArray_DOUBLE) 
    dtrans = (double *)array->data;
  else
    ctrans = (complex *)array->data;

  //printf("lanthanide_transform: fill matrix\n");
  result.states = terms;
  if (type == PyArray_DOUBLE) 
    err = realtransform(&matrix, dtrans, &result);
  else
    err = cplxtransform(&matrix, ctrans, &result);
  if (err != 0)
    {
      PyErr_SetString(PyExc_ValueError, "transform: calc error");
      return NULL;
    }

  //printf("lanthanide_transform: convert result\n");
  mdict = NULL;
  if (SparseToPy(&result, &mdict))
    {
      PyErr_SetString(PyExc_ValueError, "transform: convert result error");
      return NULL;
    }
  
  //printf("lanthanide_transform: free allocated memory\n");
  freeSparse(&matrix);
  freeSparse(&result);
  Py_DECREF(array);

  //printf("lanthanide_transform: return\n");
  return mdict;
}



//--------------------------------------------------------------------
PyObject *lanthanide_transdiag(PyObject *self, PyObject *args)
{
  int type, dims[1];
  long size, terms;
  double *dtrans, *dresult;
  complex *ctrans, *cresult;
  sparseStruct matrix;
  PyObject *dict;
  PyArrayObject *array, *dest;

  //printf("lanthanide_transdiag: parse args\n");
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyDict_Type, &dict,
			&PyArray_Type, &array))
    {
      PyErr_SetString(PyExc_ValueError, "transdiag: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_transdiag: set type\n");
  type = array->descr->type_num;
  if ((type != PyArray_CDOUBLE) &&
      (type != PyArray_DOUBLE))
    {
      PyErr_SetString(PyExc_ValueError, "transdiag: trans has wrong type");
      return NULL;
    }
  if (type == PyArray_DOUBLE)
    size = 1;
  else
    size = 2;

  //printf("lanthanide_transdiag: convert matrix\n");
  if (SparseFromPy(dict, &matrix))
    {
      PyErr_SetString(PyExc_ValueError, "transdiag: convert matrix error");
      return NULL;
    }

  //printf("lanthanide_transdiag: test matrix size\n");
  if (matrix.size > size)
    {
      PyErr_SetString(PyExc_ValueError, "transdiag: matrix dimension error");
      return NULL;
    }

  //printf("lanthanide_transdiag: test trans matrix\n");
  if ((array->nd != 2) ||
      (array->dimensions[1] != matrix.states))
    {
      PyErr_SetString(PyExc_ValueError, "transdiag: trans type error");
      return NULL;
    }
  terms = array->dimensions[0];

  //printf("lanthanide_transdiag: make trans matrix contiguous\n");
  array = (PyArrayObject *)
    PyArray_ContiguousFromObject((PyObject *)array, 
				 array->descr->type_num, 
				 array->nd,
				 array->nd);
  if (type == PyArray_DOUBLE)
  {
      ctrans = NULL;
      dtrans = (double *)array->data;
  }
  else
  {
      ctrans = (complex *)array->data;
      dtrans = NULL;
  }

  //printf("lanthanide_transdiag: allocate dest\n");
  dims[0] = terms;
  dest = (PyArrayObject *)PyArray_FromDims(1, dims, type);
  if (dest == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "transdiag: can't create result");
      return NULL;
    }
  if (type == PyArray_DOUBLE)
  {
      cresult = NULL;
      dresult = (double *)dest->data;
  }
  else
  {
      cresult = (complex *)dest->data;
      dresult = NULL;
  }

  //printf("lanthanide_transdiag: calc transdiag\n");
  if (type == PyArray_DOUBLE)
    realtransdiag(&matrix, dtrans, dresult, terms);
  else
    cplxtransdiag(&matrix, ctrans, cresult, terms);

  //printf("lanthanide_transdiag: free allocated memory\n");
  freeSparse(&matrix);
  Py_DECREF(array);
  return (PyObject *)dest;
}



//--------------------------------------------------------------------
PyObject *lanthanide_submatrix(PyObject *self, PyObject *args)
{
  int type, dims[1];
  long substart, subend, diff, states;
  double *tri, *result;
  PyArrayObject *pytri, *pyresult;
  
  //printf("lanthanide_submatrix: parse args\n");
  if (!PyArg_ParseTuple(args, "O!ll",
			&PyArray_Type, &pytri,
			&substart,
			&subend))
    {
      PyErr_SetString(PyExc_ValueError, "submatrix: wrong arguments");
      return NULL;
    }

  //printf("lanthanide_submatrix: test indices\n");
  if ((substart < 0) ||
      (subend < substart))
    {
      PyErr_SetString(PyExc_ValueError, "submatrix: wrong indices");
      return NULL;
    }

  //printf("lanthanide_submatrix: set type\n");
  type = pytri->descr->type_num;
  if (type != PyArray_DOUBLE)
    {
      PyErr_SetString(PyExc_ValueError, "submatrix: tri has wrong type");
      return NULL;
    }

  //printf("lanthanide_submatrix: test tri\n");
  states = (int)((sqrt(1 + 8.0*pytri->dimensions[0]) - 1.0) / 2.0);
  if ((pytri->nd != 1) ||
      (pytri->dimensions[0] != states*(states+1)/2))
    {
      PyErr_SetString(PyExc_ValueError, "submatrix: tri type error");
      return NULL;
    }

  //printf("lanthanide_submatrix: make tri contiguous\n");
  pytri = (PyArrayObject *)
    PyArray_ContiguousFromObject((PyObject *)pytri, 
				 pytri->descr->type_num, 
				 pytri->nd,
				 pytri->nd);
  tri = (double*)pytri->data;

  //printf("lanthanide_submatrix: allocate dest\n");
  diff = subend - substart;
  dims[0] = diff * (diff+1) / 2;
  pyresult = (PyArrayObject *)PyArray_FromDims(1, dims, type);
  if (pyresult == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "submatrix: can't create result");
      return NULL;
    }
  result = (double *)pyresult->data;

  //printf("lanthanide_submatrix: get submatrix\n");
  if (diff != 0)
    subtriangle(tri, states, (long)substart, (long)subend, result);

  //printf("lanthanide_submatrix: free variables\n");
  Py_DECREF(pytri);

  //printf("lanthanide_submatrix: return\n");
  return (PyObject *)pyresult;
}
