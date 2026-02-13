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

#include "lanthanide.h"



//--------------------------------------------------------------------
// prototypes of private functions
void _inccount(long *count, long max, long i);



//--------------------------------------------------------------------
long binom(long n, long k)
{
  long i, kk;
  double prod;

  kk = MIN(k, n-k);

  prod = 1.0;
  for(i=0; i<kk; i++)
      prod = prod * (n-i) / (kk-i);

  return(DTOL(prod));
}



//--------------------------------------------------------------------
void freeQuant(quantStruct *quant)
{
  if (quant->ownbuf)
    {
      free(quant->ml);
      free(quant->ms);
      quant->ml = NULL;
      quant->ms = NULL;
    }
}



//--------------------------------------------------------------------
long buildquant(long l, quantStruct *quant)
{
  long i, max, ml, ms;

  //printf("buildquant: set values\n");
  max = 2*l+2;
  quant->l      = l;
  quant->max    = max;
  quant->ownbuf = 1;

  //printf("buildquant: malloc structures\n");
  if (!(quant->ml = malloc(max*sizeof(long))))
    return -1;
  if (!(quant->ms = malloc(max*sizeof(long))))
    return -2;

  //printf("buildquant: fill ml and ms\n");
  i = 0;
  for (ml=l; ml>=-l; ml-=2)
    for (ms=1; ms>=-1; ms-=2)
      {
	quant->ml[i] = ml;
	quant->ms[i] = ms;
	i++;
      }

  //printf("buildquant: return\n");
  return 0;
}



//--------------------------------------------------------------------
void freeUnit(unitStruct *unit)
{
  if (unit->ownbuf)
    {
      free(unit->u);
      free(unit->t);
      unit->u = NULL;
      unit->t = NULL;
    }
}



//--------------------------------------------------------------------
long buildunit(quantStruct *quant, unitStruct *unit)
{
  long l, max, size;
  long i, j, ma, mb, k, q;

  //printf("buildunit: set values\n");
  l   = quant->l;
  max = quant->max;
  unit->l      = l;
  unit->max    = max;
  unit->ownbuf = 1;

  //printf("buildunit: malloc structures\n");
  size = (l+1) * max * max;
  if (!(unit->u = malloc(size * sizeof(double))))
    return -1;
  size = 2 * max * max;
  if (!(unit->t = malloc(size * sizeof(double))))
    return -2;

  //printf("buildunit: fill u\n");
  for (i=0; i<max; i++)
    for (j=0; j<max; j++)
      {
	ma = quant->ml[i];
	mb = quant->ml[j];
	q  = ma-mb;
	for (k=0; k<=2*l; k+=2)
	  U(unit, k, i, j) = 
	    SIGNEXP2(l-ma) * wign3j(l, k, l, -ma, q, mb);
      }

  //printf("buildunit: fill t\n");
  for (i=0; i<max; i++)
    for (j=0; j<max; j++)
      {
	ma = quant->ms[i];
	mb = quant->ms[j];
	q  = ma-mb;
	for (k=0; k<=2; k+=2)
	  T(unit, k, i, j) = 
	    SIGNEXP2(1-ma) * wign3j(1, k, 1, -ma, q, mb);
      }

  //printf("buildunit: return\n");
  return 0;
}



//--------------------------------------------------------------------
void freeConfig(configStruct *config)
{
  if (config->ownbuf)
    {
      free(config->state);
      config->state = NULL;
    }
}



//--------------------------------------------------------------------
long buildconfig(quantStruct *quant, long electrons, configStruct *config)
{
  long i, j, *count, l, max;

  //printf("buildconfig: set values\n");
  l = quant->l;
  config->l         = l;
  config->electrons = electrons;
  config->states    = binom(quant->max, electrons);
  config->ownbuf    = 1;

  //printf("buildconfig: malloc state\n");
  if (!(config->state = malloc(config->states*electrons*sizeof(long))))
    return -1;

  //printf("buildconfig: fill count\n");
  if ((count = malloc(electrons*sizeof(long))) == NULL)
    return -2;
  for (i=0; i<electrons; i++)
    count[i] = i;

  //printf("buildconfig: build states\n");
  max = quant->max;
  for (i=0; i<config->states; i++)
    {
      for (j=0; j<electrons; j++)
	STATE(config, i, j) = count[j];
      if (count[0] < max-electrons)
	_inccount(count, max, electrons-1);
    }

  //printf("buildconfig: return\n");
  free(count);
  return 0;
}



//--------------------------------------------------------------------
void _inccount(long *count, long max, long i)
{
  count[i] = count[i] + 1;
  if (count[i] >= max)
    {
      _inccount(count, max, i-1);
      count[i] = count[i-1];
      _inccount(count, max, i);
    }
  return;
}



//--------------------------------------------------------------------
void freePair(pairStruct *pair)
{
  if (pair->ownbuf)
    {
      free(pair->statei);
      free(pair->statej);
      pair->statei = NULL;
      pair->statej = NULL;
    }
}



//--------------------------------------------------------------------
long initpair(configStruct *config, pairStruct *pair)
{
  long electrons;

  //printf("initpair: set values\n");
  electrons = config->electrons;
  pair->electrons = electrons;
  pair->ownbuf    = 1;

  //printf("initpair: malloc statei\n");
  if (!(pair->statei = malloc(electrons * sizeof(long))))
    return -1;

  //printf("initpair: malloc statej\n");
  if (!(pair->statej = malloc(electrons * sizeof(long))))
    return -2;

  //printf("initpair: return\n");
  return 0;
}



//--------------------------------------------------------------------
long orderpair(configStruct *config, long i, long j, pairStruct *pair)
{
  long k, l;
  long electrons;

  electrons = config->electrons;
  
  if (electrons != pair->electrons)
    return -1;

  pair->i = i;
  pair->j = j;

  for (k=0; k<electrons; k++)
    {
      BRA(pair, k) = STATE(config, i, k);
      KET(pair, k) = STATE(config, j, k);
    }
  
  pair->sign = 1;
  for (k=0; k<electrons; k++)
    {
      for (l=0; l<electrons; l++)
	{
	  if(BRA(pair, k) == KET(pair, l))
	    {
	      if(k != l)
		{
		  SWAP(KET(pair, k), KET(pair, l));
		  pair->sign = -pair->sign;
		}
	      break;
	    }
	} 
    }
  //if ((i-j) % 2 != 0)
  //  pair->sign = -pair->sign;

  pair->unpaired = 0;
  for(k=0; k<electrons; k++)
    if(BRA(pair, k) != KET(pair, k))
      {
	if(k != pair->unpaired)
	  {
	    SWAP(BRA(pair, k), BRA(pair, pair->unpaired));
	    SWAP(KET(pair, k), KET(pair, pair->unpaired));
	  }
	pair->unpaired++;
      }

  return 0;
}



//--------------------------------------------------------------------
void freeSparse(sparseStruct *matrix)
{
  if (matrix->ownbuf)
    {
      free(matrix->data);
      free(matrix->col);
      free(matrix->row);
      matrix->data = NULL;
      matrix->col = NULL;
      matrix->row = NULL;
    }
}



//--------------------------------------------------------------------
long calcmatrix(quantStruct *quant, unitStruct *unit, configStruct *config,
		long key, long *opts, sparseStruct *matrix)
{
  long i, j, index;
  long total;
  double val;
  pairStruct pair;

  //printf("calcmatrix: some predefinitions\n");
  total   = config->states * MEMFACTOR;

  //printf("calcmatrix: allocate row memory and init structure\n");
  if (!(matrix->row = malloc(config->states * sizeof(long))))
    return -1;
  matrix->states = config->states;
  matrix->size   = 1;
  matrix->ownbuf = 1;

  //printf("calcmatrix: test for zero matrix\n");
/* Be careful: the following is wrong!
  if ((config->electrons < ELEMENT[key].electrons) ||
      (config->electrons > quant->max - ELEMENT[key].electrons))
*/
  if (config->electrons < ELEMENT[key].electrons)
    {
      matrix->data = NULL;
      matrix->col  = NULL;
      for (i=0; i<config->states; i++)
	matrix->row[i] = 0;
      return 0;
    }

  //printf("calcmatrix: allocate data and col memory\n");
  if (!(matrix->data = malloc(total * sizeof(double))))
    return -2;
  if (!(matrix->col = malloc(total * sizeof(long))))
    return -3;

  //printf("calcmatrix: init pair\n");
  if (initpair(config, &pair) != 0)
    return -4;

  //printf("calcmatrix: calculation loop\n");
  index = 0;
  for (j=0; j<matrix->states; j++)
    {
      for (i=0; i<=j; i++)
	{
	  if (orderpair(config, i, j, &pair))
	    return -5;
	  
	  switch (ELEMENT[key].electrons)
	    {
	    case 1:
	      val = elementOne(quant, unit, &pair, key, opts);
	      break;
	    case 2:
	      val = elementTwo(quant, unit, &pair, key, opts);
	      break;
	    case 3:
	      val = elementThree(quant, unit, &pair, key, opts);
	      break;
	    default:
	      return -6;
	    }

	  if (!DZERO(val))
	    {
	      matrix->data[index] = val;
	      matrix->col[index++] = i;
	    }
	}
      matrix->row[j] = index;
      
      if ((total - index) < (j+1))
	{
	  //printf("calcmatrix: realloc to larger size\n");
	  fprintf(stderr, "WARNING! calcmatrix: matrix too small!\n");
	  total = index + matrix->states * MEMFACTOR;
	  if(!(matrix->data = 
	       realloc(matrix->data, total * sizeof(double))))
	    return -7;
	  if(!(matrix->col = 
	       realloc(matrix->col, total * sizeof(long))))
	    return -8;
	}
    }

  //printf("calcmatrix: realloc to actual size\n");
  if (index > 0)
    {
      if(!(matrix->data = 
	   realloc(matrix->data, index * sizeof(double))))
	return -9;
      if(!(matrix->col = 
	   realloc(matrix->col, index * sizeof(long))))
	return -10;
    }
  else
    {
      free(matrix->data);
      free(matrix->col);
      matrix->data = NULL;
      matrix->col = NULL;
    }

  //printf("calcmatrix: free pair\n");
  freePair(&pair);

  //printf("calcmatrix: return\n");
  return 0;
}



//--------------------------------------------------------------------
long calcelement(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		 long key, long *opts, double *data)
{
  if ((pair->electrons < ELEMENT[key].electrons) ||
      (pair->electrons > quant->max - ELEMENT[key].electrons))
    {
      *data = 0.0;
      return 0;
    }

  switch (ELEMENT[key].electrons)
    {
    case 1:
      *data = elementOne(quant, unit, pair, key, opts);
      break;
    case 2:
      *data = elementTwo(quant, unit, pair, key, opts);
      break;
    case 3:
      *data = elementThree(quant, unit, pair, key, opts);
      break;
    default:
      return -1;
    }

  if (DZERO(*data))
    *data = 0.0;
  return 0;
}



//--------------------------------------------------------------------
double realgetelement(sparseStruct *matrix, long i, long j)
{
  long a, b, p;

  if (i > j)
    SWAP(i, j);
  
  a = 0;
  if (j > 0)
    a = matrix->row[j-1];
  b = matrix->row[j];

  if (a == b)
    return 0.0;
  
  while (matrix->col[a] != i)
    {
      if (b <= a+1)
	return 0.0;

      p = (a+b)/2;
      if (matrix->col[p] <= i)
	a = p;
      else
	b = p;
    }

  return matrix->data[a];
}



//--------------------------------------------------------------------
complex cplxgetelement(sparseStruct *matrix, long i, long j)
{
  long a, b, p;
  complex result;

  result.r = 0.0;
  result.i = 0.0;

  if (i > j)
    SWAP(i, j);
  
  a = 0;
  if (j > 0)
    a = matrix->row[j-1];
  b = matrix->row[j];

  if (a == b)
    return result;
  
  while (matrix->col[a] != i)
    {
      if (b <= a+1)
	return result;

      p = (a+b)/2;
      if (matrix->col[p] <= i)
	a = p;
      else
	b = p;
    }

  result.r = matrix->data[a*2];
  result.i = matrix->data[a*2+1];
  return result;
}



//--------------------------------------------------------------------
long realmultadd(sparseStruct *matrix, double factor, double *dest)
{
  long j, p, ps, k, offset;

  if (matrix->size != 1)
    return -1;

  if (factor == 0.0)
    return 0;

  ps = 0;
  for (j=0; j<matrix->states; j++)
    {
      offset = j * (j+1) / 2;
      for (p=ps; p<matrix->row[j]; p++)
	{
	  k = matrix->col[p] + offset;
	  dest[k] = dest[k] + factor * matrix->data[p];
	}
      ps = matrix->row[j];
    }

  return 0;
}



//--------------------------------------------------------------------
long cplxmultadd(sparseStruct *matrix, complex factor, complex *dest)
{
  long j, p, ps, k, offset;

  if ((matrix->size != 1) &&
      (matrix->size != 2))
    return -1;

  if ((factor.r == 0.0) && (factor.i == 0.0))
    return 0;

  if (matrix->size == 1)
    {
      ps = 0;
      for (j=0; j<matrix->states; j++)
	{
	  offset = j * (j+1) / 2;
	  for (p=ps; p<matrix->row[j]; p++)
	    {
	      k = matrix->col[p] + offset;
	      dest[k].r = dest[k].r + factor.r * matrix->data[p];
	      dest[k].i = dest[k].i + factor.i * matrix->data[p];
	    }
	  ps = matrix->row[j];
	}
    }
  else
    {
      ps = 0;
      for (j=0; j<matrix->states; j++)
	{
	  offset = j * (j+1) / 2;
	  for (p=ps; p<matrix->row[j]; p++)
	    {
	      k = matrix->col[p] + offset;
	      dest[k].r = dest[k].r +
		factor.r * matrix->data[2*p] -
		factor.i * matrix->data[2*p+1];
	      dest[k].i = dest[k].i +
		factor.r * matrix->data[2*p+1] +
		factor.i * matrix->data[2*p];
	    }
	  ps = matrix->row[j];
	}
    }

  return 0;
}



//--------------------------------------------------------------------
long realtransform(sparseStruct *matrix, double *trans,
		   sparseStruct *result)
{
  long i, j, k, l, index;
  long ps, p;
  long states, terms, total;
  double value, sum;
  double *tl, *tk;
  
  //printf("realtransform: stage 1\n");
  if (matrix->size != 1)
    return -1;
  states = matrix->states;
  terms  = result->states;
  total  = states * 5*MEMFACTOR;
  
  //printf("realtransform: stage 2\n");
  if (!(result->row = malloc(states * sizeof(long))))
    return -2;
  result->size   = 1;
  result->ownbuf = 1;

  //printf("realtransform: stage 3\n");
  if (!(result->data = malloc(total * sizeof(double))))
    return -3;
  if (!(result->col = malloc(total * sizeof(long))))
    return -4;

  //printf("realtransform: stage 4\n");
  index = 0;
  for (l=0; l<terms; l++)
    {
      tl = &trans[l*states];
      for (k=0; k<=l; k++)
	{
	  tk = &trans[k*states];
	  sum = 0.0;
	  ps = 0;
	  for (j=0; j<states; j++)
	    {
	      for (p=ps; p<matrix->row[j]; p++)
		{
		  i = matrix->col[p];
		  value = tk[i] * tl[j];
		  if (i != j)
		    value = value + tk[j] * tl[i];
		  sum += matrix->data[p] * value;
		}
	      ps = matrix->row[j];
	    }
	  if (!DZERO(sum))
	    {
	      result->data[index] = sum;
	      result->col[index++] = k;
	    }
	}
      result->row[l] = index;

      if (total-index < l+1)
	{
	  fprintf(stderr, "WARNING! realtransform: result too small!\n");
	  total = index + result->states * 5*MEMFACTOR;
	  if(!(result->data = 
	       realloc(result->data, total * sizeof(double))))
	    return -5;
	  if(!(result->col = 
	       realloc(result->col, total * sizeof(long))))
	    return -6;
	}
    }

  //printf("realtransform: stage 5\n");
  if (index > 0)
    {
      if(!(result->data = 
	   realloc(result->data, index * sizeof(double))))
	return -7;
      if(!(result->col = 
	   realloc(result->col, index * sizeof(long))))
	return -8;
    }
  else
    {
      free(result->data);
      free(result->col);
      result->data = NULL;
      result->col = NULL;
    }

  //printf("realtransform: stage 6\n");
  return 0;
}



//--------------------------------------------------------------------
long cplxtransform(sparseStruct *matrix, complex *trans,
		   sparseStruct *result)
{
  long i, j, k, l, index;
  long ps, p;
  long states, total;
  complex value, sum;
  complex *tl, *tk;

  //printf("cplxtransform: stage 1\n");
  if ((matrix->size != 1) &&
      (matrix->size) != 2)
    return -1;
  states = matrix->states;
  total  = states * MEMFACTOR;

  //printf("cplxtransform: stage 2\n");
  if (!(result->row = malloc(states * sizeof(long))))
    return -2;
  result->states = states;
  result->size   = 2;
  result->ownbuf = 1;

  //printf("cplxtransform: stage 3\n");
  if (!(result->data = malloc(total * sizeof(complex))))
    return -3;
  if (!(result->col = malloc(total * sizeof(long))))
    return -4;

  //printf("cplxtransform: stage 4\n");
  index = 0;
  j = 0;
  for (l=0; l<states; l++)
    {
      tl = &trans[l*states];
      for (k=0; k<=l; k++)
	{
	  tk = &trans[k*states];
	  sum.r = 0.0;
	  sum.i = 0.0;
	  ps = 0;
	  for (j=0; j<states; j++)
	    {
	      for (p=ps; p<matrix->row[j]; p++)
		{
		  i = matrix->col[p];
		  value.r = tk[i].r * tl[j].r - tk[i].i * tl[j].i;
		  value.i = tk[i].r * tl[j].i + tk[i].i * tl[j].r;
		  if (i != j)
		    {
		      value.r = value.r +
			tk[j].r * tl[i].r - tk[j].i * tl[i].i;
		      value.i = value.i +
			tk[j].r * tl[i].i + tk[j].i * tl[i].r;
		    }
		  if (matrix->size == 1)
		    {
		      sum.r += matrix->data[p] * value.r;
		      sum.i += matrix->data[p] * value.i;
		    }
		  else
		    {
		      sum.r +=
			matrix->data[2*p] * value.r -
			matrix->data[2*p+1] * value.i;
		      sum.i +=
			matrix->data[2*p] * value.i +
			matrix->data[2*p+1] * value.r;
		    }
		}
	      ps = matrix->row[j];
	    }
	  if (!CZERO(sum))
	    {
	      result->data[2*index] = sum.r;
	      result->data[2*index+1] = sum.i;
	      result->col[index++] = k;
	    }
	}
      result->row[l] = index;

      if (total-index < j+1)
	{
	  fprintf(stderr, "WARNING! cplxtransform: result too small!\n");
	  total = index + result->states * MEMFACTOR;
	  if(!(result->data = 
	       realloc(result->data, total * sizeof(complex))))
	    return -5;
	  if(!(result->col = 
	       realloc(result->col, total * sizeof(long))))
	    return -6;
	}
    }

  //printf("cplxtransform: stage 5\n");
  if (index > 0)
    {
      if(!(result->data = 
	   realloc(result->data, index * sizeof(complex))))
	return -7;
      if(!(result->col = 
	   realloc(result->col, index * sizeof(long))))
	return -8;
    }
  else
    {
      free(result->data);
      free(result->col);
      result->data = NULL;
      result->col = NULL;
    }

  //printf("cplxtransform: stage 6\n");
  return 0;
}



//--------------------------------------------------------------------
long realtransdiag(sparseStruct *matrix, double *trans, double *result, long terms)
{
  long k, i, j;
  long states;
  long ps, p;
  double value, sum;
  double *tk;

  if (matrix->size != 1)
    return -1;

  states = matrix->states;
  for (k=0; k<terms; k++)
    {
      tk = &trans[k*states];
      sum = 0.0;
      ps = 0;
      for (j=0; j<states; j++)
	{
	  for (p=ps; p<matrix->row[j]; p++)
	    {
	      i = matrix->col[p];
	      value = tk[i] * tk[j];
	      if (i != j)
		value = 2.0 * value;
	      sum += matrix->data[p] * value;
	    }
	  ps = matrix->row[j];
	}
      if (!DZERO(sum))
	result[k] = sum;
      else
	result[k] = 0.0;
    }

  return 0;
}



//--------------------------------------------------------------------
long cplxtransdiag(sparseStruct *matrix, complex *trans, complex *result, long terms)
{
  long k, i, j;
  long states;
  long ps, p;
  complex value, sum;
  complex *tk;

  if ((matrix->size != 1) &&
      (matrix->size != 2))
    return -1;

  states = matrix->states;
  for (k=0; k<terms; k++)
    {
      tk = &trans[k*states];
      sum.r = 0.0;
      sum.i = 0.0;
      ps = 0;
      for (j=0; j<states; j++)
	{
	  for (p=ps; p<matrix->row[j]; p++)
	    {
	      i = matrix->col[p];
	      value.r = tk[i].r * tk[j].r - tk[i].i * tk[j].i;
	      value.i = tk[i].r * tk[j].i + tk[i].i * tk[j].r;
	      if (i != j)
		{
		  value.r = 2.0 * value.r;
		  value.i = 2.0 * value.i;
		}
	      if (matrix->size == 1)
		{
		  sum.r += matrix->data[p] * value.r;
		  sum.i += matrix->data[p] * value.i;
		}
	      else
		{
		  sum.r +=
		    matrix->data[2*p] * value.r -
		    matrix->data[2*p+1] * value.i;
		  sum.i +=
		    matrix->data[2*p] * value.i +
		    matrix->data[2*p+1] * value.r;
		}
	    }
	  ps = matrix->row[j];
	}
      if (!CZERO(sum))
	{
	  result[k].r = sum.r;
	  result[k].i = sum.i;
	}
      else
	{
	  result[k].r = 0.0;
	  result[k].i = 0.0;
	}
    }

  return 0;
}



//--------------------------------------------------------------------
long subtriangle(double *tri, long states, long substart, long subend,
		 double *result)
{
  long i, j, k, d;

  k = 0;
  for (j=substart; j<subend; j++)
  {
    d = j*(j+1)/2;
    for (i=substart; i<j+1; i++)
      result[k++] = tri[d+i];
  }
  return 0;
}
