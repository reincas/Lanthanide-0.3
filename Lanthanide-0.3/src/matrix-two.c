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
#define UAC(k) UVAL(unit, k, A, C)
#define UBD(k) UVAL(unit, k, B, D)

#define S      1
#define MSA    MS(quant, A)
#define MSB    MS(quant, B)
#define MSC    MS(quant, C)
#define MSD    MS(quant, D)
#define TAC(k) TVAL(unit, k, A, C)
#define TBD(k) TVAL(unit, k, B, D)



//--------------------------------------------------------------------
#define A  bra[0]
#define B  bra[1]
#define C  ket[0]
#define D  ket[1]
double elementTwo(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		  long key, long *opts)
{
  long i, maxi, j, maxj;
  long swap;
  long bra[2], ket[2];
  double sum;
  elementFunc element;

  if (pair->unpaired > 2)
    return 0.0;

  element = ELEMENT[key].element;
  swap    = ELEMENT[key].swap;

  maxi = pair->electrons-1;
  if (pair->unpaired > 0)
    maxi = 1;

  maxj = pair->electrons;
  if (pair->unpaired > 1)
    maxj = 2;

  sum = 0.0;
  for (i=0; i<maxi; i++)
    for (j=i+1; j<maxj; j++)
      {
	A = BRA(pair, i);
	B = BRA(pair, j);
	C = KET(pair, i);
	D = KET(pair, j);

	switch (swap)
	  {
	  case 1:
	    sum += element(quant, unit, bra, ket, opts);
	    SWAP(C, D);
	    sum -= element(quant, unit, bra, ket, opts);
	    break;
	    
	  default:
	    sum += element(quant, unit, bra, ket, opts);
	    SWAP(C, D);
	    sum -= element(quant, unit, bra, ket, opts);
	    SWAP(A, B);
	    sum += element(quant, unit, bra, ket, opts);
	    SWAP(C, D);
	    sum -= element(quant, unit, bra, ket, opts);
	  }
      }
  if (swap == 2)
    sum /= 2;
  return(sum * pair->sign);
}



//--------------------------------------------------------------------
#define K  opts[0]
double two_uu(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MLA+MLB, MLC+MLD) ||
     NEQ(MSA, MSC) ||
     NEQ(MSB, MSD))
    return 0.0;

  return SIGNEXP2(MLA-MLC) * UAC(K) * UBD(K);
}
#undef K



//--------------------------------------------------------------------
#define K1 opts[0]
#define K2 opts[1]
#define Q  opts[2]
double two_uus(quantStruct *quant, unitStruct *unit,
	       long *bra, long *ket, long *opts)
{
  if(NEQ(MLA+MLB, MLC+MLD) ||
     NEQ(MSA, MSC) ||
     NEQ(MSB, MSD) ||
     NEQ(Q, MLA-MLC))
    return 0.0;

  return SIGNEXP2(Q) * UAC(K1) * UBD(K2);
}
#undef K1
#undef K2
#undef Q



//--------------------------------------------------------------------
#define K  opts[0]
double two_tt(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MSA+MSB, MSC+MSD) ||
     NEQ(MLA, MLC) ||
     NEQ(MLB, MLD))
    return 0.0;

  return SIGNEXP2(MSA-MSC) * TAC(K) * TBD(K);
}
#undef K



//--------------------------------------------------------------------
#define K  opts[0]
double two_ut(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MLA+MSB, MLC+MSD) ||
     NEQ(MSA, MSC) ||
     NEQ(MLB, MLD))
    return 0.0;

  return SIGNEXP2(MLA-MLC) * UAC(K) * TBD(K);
}
#undef K


//--------------------------------------------------------------------
#define K1  opts[0]
#define K2  opts[1]
#define R1  opts[2]
#define R2  opts[3]
#define K   opts[4]
double two_uutt(quantStruct *quant, unitStruct *unit,
		long *bra, long *ket, long *opts)
{
  long q, q1, q2, p1, p2;
  
  q1 = MLA - MLC;
  q2 = MLB - MLD;
  p1 = MSA - MSC;
  p2 = MSB - MSD;
  q  = q1 + q2;
  
  if(q1 + q2 + p1 + p2 != 0)
    return 0.0;

  return
    SIGNEXP2(q + K1 - K2 + R1 - R2) *
    (K + 1.0) *
    wign3j(K1, K2, K, q1, q2, -q) *
    wign3j(R1, R2, K, p1, p2, q) *
    UAC(K1) * UBD(K2) * TAC(R1) * TBD(R2);
}
#undef K1
#undef K2
#undef R1
#undef R2
#undef K
