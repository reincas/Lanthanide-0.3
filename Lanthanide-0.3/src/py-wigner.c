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
PyObject *lanthanide_wign3j(PyObject *self, PyObject *args)
{
  double a, b, c, x, y, z;

  if (!PyArg_ParseTuple(args, "dddddd", &a, &b, &c, &x, &y, &z))
    {
      PyErr_SetString(PyExc_ValueError, "wign3j: wrong arguments");
      return NULL;
    }
  return PyFloat_FromDouble(wign3j(DTO2L(a), DTO2L(b), DTO2L(c),
				   DTO2L(x), DTO2L(y), DTO2L(z)));
}



//--------------------------------------------------------------------
PyObject *lanthanide_clebsg(PyObject *self, PyObject *args)
{
  double a, b, c, x, y, z;

  if (!PyArg_ParseTuple(args, "dddddd", &a, &b, &c, &x, &y, &z))
    {
      PyErr_SetString(PyExc_ValueError, "clebsg: wrong arguments");
      return NULL;
    }
  return PyFloat_FromDouble(clebsg(DTO2L(a), DTO2L(b), DTO2L(c),
				   DTO2L(x), DTO2L(y), DTO2L(z)));
}



//--------------------------------------------------------------------
PyObject *lanthanide_wign6j(PyObject *self, PyObject *args)
{
  double a, b, c, x, y, z;

  if (!PyArg_ParseTuple(args, "dddddd", &a, &b, &c, &x, &y, &z))
    {
      PyErr_SetString(PyExc_ValueError, "wign6j: wrong arguments");
      return NULL;
    }
  return PyFloat_FromDouble(wign6j(DTO2L(a), DTO2L(b), DTO2L(c),
				   DTO2L(x), DTO2L(y), DTO2L(z)));
}



//--------------------------------------------------------------------
PyObject *lanthanide_racahc(PyObject *self, PyObject *args)
{
  double a, b, c, x, y, z;

  if (!PyArg_ParseTuple(args, "dddddd", &a, &b, &c, &x, &y, &z))
    {
      PyErr_SetString(PyExc_ValueError, "racahc: wrong arguments");
      return NULL;
    }
  return PyFloat_FromDouble(racahc(DTO2L(a), DTO2L(b), DTO2L(c),
				   DTO2L(x), DTO2L(y), DTO2L(z)));
}



//--------------------------------------------------------------------
PyObject *lanthanide_jahnuf(PyObject *self, PyObject *args)
{
  double a, b, c, x, y, z;

  if (!PyArg_ParseTuple(args, "dddddd", &a, &b, &c, &x, &y, &z))
    {
      PyErr_SetString(PyExc_ValueError, "jahnuf: wrong arguments");
      return NULL;
    }
  return PyFloat_FromDouble(jahnuf(DTO2L(a), DTO2L(b), DTO2L(c),
				   DTO2L(x), DTO2L(y), DTO2L(z)));
}



//--------------------------------------------------------------------
PyObject *lanthanide_wign9j(PyObject *self, PyObject *args)
{
  double a, b, c, x, y, z, g, h, p;

  if (!PyArg_ParseTuple(args, "ddddddddd",
			&a, &b, &c, &x, &y, &z, &g, &h, &p))
    {
      PyErr_SetString(PyExc_ValueError, "wign9j: wrong arguments");
      return NULL;
    }
  return PyFloat_FromDouble(wign9j(DTO2L(a), DTO2L(b), DTO2L(c),
				   DTO2L(x), DTO2L(y), DTO2L(z),
				   DTO2L(g), DTO2L(h), DTO2L(p)));
}
