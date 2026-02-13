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
PyObject *lanthanide_diagonalize(PyObject *self, PyObject *args)
{
  long method, vecflag, n;
  int elements, type;
  PyArrayObject *matap, *matw, *matz;
  double *dap, *dz;
  complex *cap, *cz;
  double *w;
  long ldz;
  int w_dim[1], z_dim[2];
  long info;

  // parse arguments
  if (!PyArg_ParseTuple(args, "llO!",
			&method,
			&vecflag,
			&PyArray_Type, &matap))
    {
      PyErr_SetString(PyExc_ValueError, "wrong arguments");
      return NULL;
    }
  type = matap->descr->type_num;
  if ((type != PyArray_CDOUBLE) &&
      (type != PyArray_DOUBLE))
    {
      PyErr_SetString(PyExc_ValueError, "ap has wrong type");
      return NULL;
    }

  // check properties of ap
  matap = (PyArrayObject *)
    PyArray_ContiguousFromObject((PyObject *)matap, 
				 matap->descr->type_num, 
				 matap->nd,
				 matap->nd);
  elements = matap->dimensions[0];
  n = (int)((sqrt(1 + 8.0*matap->dimensions[0]) - 1.0) / 2.0);
  if ((matap == NULL) ||
      (matap->nd != 1) ||
      (n*(n+1)/2 != elements))
    {
      PyErr_SetString(PyExc_ValueError, "ap has wrong properties");
      return NULL;
    }
  dap = (double *)matap->data;
  cap = (complex *)matap->data;

  // initialize vector w
  w_dim[0] = (int)n;
  matw = (PyArrayObject *)PyArray_FromDims(1, w_dim, PyArray_DOUBLE);
  if (matw == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "can't create w");
      return NULL;
    }
  w = (double *)matw->data;

  // set minimum value of ldz
  if (vecflag == DIAG_NOVEC)
    ldz = 1;
  else
    ldz = n;

  // initialize ldz and array z(ldz,n)
  z_dim[0] = (int)ldz;
  z_dim[1] = (int)n;
  matz = (PyArrayObject *)PyArray_FromDims(2, z_dim, type);
  if (matz == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "can't create z");
      return NULL;
    }
  dz = (double *)matz->data;
  cz = (complex *)matz->data;

  // diagonalize
  if (type == PyArray_DOUBLE)
    info = diagsym(method, vecflag, n, dap, w, dz);
  else
    info = diagherm(method, vecflag, n, cap, w, cz);

  // return tuple containing w and z
  Py_DECREF(matap);
  return Py_BuildValue("(NNi)",
		       PyArray_Return(matw),
		       PyArray_Return(matz),
		       info);
}
