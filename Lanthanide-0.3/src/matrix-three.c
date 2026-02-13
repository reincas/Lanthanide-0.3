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
// useful shortcuts
#define L      quant->l
#define MLA    ML(quant, A)
#define MLB    ML(quant, B)
#define MLC    ML(quant, C)
#define MLD    ML(quant, D)
#define MLE    ML(quant, E)
#define MLF    ML(quant, F)
#define UAD(k) UVAL(unit, k, A, D)
#define UBE(k) UVAL(unit, k, B, E)
#define UCF(k) UVAL(unit, k, C, F)

#define S      1
#define MSA    MS(quant, A)
#define MSB    MS(quant, B)
#define MSC    MS(quant, C)
#define MSD    MS(quant, D)
#define MSE    MS(quant, E)
#define MSF    MS(quant, F)
#define TAD(k) TVAL(unit, k, A, D)
#define TBE(k) TVAL(unit, k, B, E)
#define TCF(k) TVAL(unit, k, C, F)



//--------------------------------------------------------------------
#define A  bra[0]
#define B  bra[1]
#define C  bra[2]
#define D  ket[0]
#define E  ket[1]
#define F  ket[2]
double elementThree(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		    long key, long *opts)
{
  long i, maxi, j, maxj, k, maxk;
  long a, b, swap;
  long bra[3], ket[3];
  double sum;
  elementFunc element;

  if (pair->unpaired > 3)
    return 0.0;

  element = ELEMENT[key].element;
  swap    = ELEMENT[key].swap;

  maxi = pair->electrons-2;
  if (pair->unpaired > 0)
    maxi = 1;

  maxj = pair->electrons-1;
  if (pair->unpaired > 1)
    maxj = 2;

  maxk = pair->electrons;
  if (pair->unpaired > 2)
    maxk = 3;

  sum = 0.0;
  for (i=0; i<maxi; i++)
    for (j=i+1; j<maxj; j++)
      for (k=j+1; k<maxk; k++)
	{
	  A = BRA(pair, i);
	  B = BRA(pair, j);
	  C = BRA(pair, k);
	  D = KET(pair, i);
	  E = KET(pair, j);
	  F = KET(pair, k);

	  switch (swap)
	    {
	    case 1:
	      for (b=0; b<3; b++)
		{
		  sum += element(quant, unit, bra, ket, opts);
		  SWAP(D, E);
		  sum -= element(quant, unit, bra, ket, opts);
		  SWAP(E, F);
		}
	      break;

	    case 2:
	      for (a=0; a<3; a++)
		{
		  for (b=0; b<3; b++)
		    {
		      sum += element(quant, unit, bra, ket, opts);
		      SWAP(D, E);
		      sum -= element(quant, unit, bra, ket, opts);
		      SWAP(E, F);
		    }
		  ROT(A, B, C);
		}
	      break;

	    default:
	      for (a=0; a<3; a++)
		{
		  for (b=0; b<3; b++)
		    {
		      sum += element(quant, unit, bra, ket, opts);
		      SWAP(D, E);
		      sum -= element(quant, unit, bra, ket, opts);
		      SWAP(E, F);
		    }
		  SWAP(A, B);
		  for (b=0; b<3; b++)
		    {
		      sum -= element(quant, unit, bra, ket, opts);
		      SWAP(D, E);
		      sum += element(quant, unit, bra, ket, opts);
		      SWAP(E, F);
		    }
		  SWAP(B, C);
		}
	    }
	}
  
  if (swap == 2)
    sum /= 3;
  if (swap == 3)
    sum /= 6;
  return(sum * pair->sign);
}



//--------------------------------------------------------------------
#define K123  opts[0]
double three_uuua(quantStruct *quant, unitStruct *unit,
		  long *bra, long *ket, long *opts)
{
  long newopts[3];

  newopts[0] = K123;
  newopts[1] = K123;
  newopts[2] = K123;

  return three_uuuc(quant, unit, bra, ket, newopts);
}
#undef K123



//--------------------------------------------------------------------
#define K12  opts[0]
#define K3   opts[1]
double three_uuub(quantStruct *quant, unitStruct *unit,
		  long *bra, long *ket, long *opts)
{
  long newopts[3];

  newopts[0] = K12;
  newopts[1] = K12;
  newopts[2] = K3;

  return three_uuuc(quant, unit, bra, ket, newopts);
}
#undef K12
#undef K3



//--------------------------------------------------------------------
#define K1  opts[0]
#define K2  opts[1]
#define K3  opts[2]
double three_uuuc(quantStruct *quant, unitStruct *unit,
		  long *bra, long *ket, long *opts)
{
  long q1, q2, q3;
  //double sum;
  
  //printf("three_uuuc: <%i %i %i|(%i%i%i)|%i %i %i>\n", 
  //A, B, C, K1, K2, K3, D, E, F);
  if(NEQ(MSA, MSD) ||
     NEQ(MSB, MSE) ||
     NEQ(MSC, MSF))
    return 0.0;

  q1 = MLA - MLD;
  q2 = MLB - MLE;
  q3 = MLC - MLF;
  if (q1+q2+q3 != 0)
    return 0.0;

  //printf("%+i %+i %+i: %7.4f, %7.4f, %7.4f, %7.4f\n",
  //q1,q2,q3,wign3j(K1,K2,K3,q1,q2,q3),UAD(K1),UBE(K2),UCF(K3));
  return wign3j(K1,K2,K3,q1,q2,q3) * UAD(K1) * UBE(K2) * UCF(K3);
	    
  /*
  sum = 0.0;
  for (q1=-K1; q1<=K1; q1+=2)
    for (q2=-K2; q2<=K2; q2+=2)
      {
	q3 = -q1-q2;
	if (abs(q3) <= K3)
	  {
	    printf("%+i %+i %+i: %7.4f, %7.4f, %7.4f, %7.4f\n",
		   q1,q2,q3,wign3j(K1,K2,K3,q1,q2,q3),UAD(K1),UBE(K2),UCF(K3));
	    sum += wign3j(K1,K2,K3,q1,q2,q3) * UAD(K1) * UBE(K2) * UCF(K3);
	  }
      }
  printf("%f\n", sum);
  return sum;
  */
}
#undef K1
#undef K2
#undef K3
