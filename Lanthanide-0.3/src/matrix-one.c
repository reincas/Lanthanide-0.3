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
#define UAB(k) UVAL(unit, k, A, B)

#define S      1
#define MSA    MS(quant, A)
#define MSB    MS(quant, B)
#define TAB(k) TVAL(unit, k, A, B)



//--------------------------------------------------------------------
#define A  bra[0]
#define B  ket[0]
double elementOne(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		  long key, long *opts)
{
  long i, maxi;
  long bra[1], ket[1];
  double sum;
  elementFunc element;

  if (pair->unpaired > 1)
    return 0.0;

  element = ELEMENT[key].element;

  maxi = pair->electrons;
  if (pair->unpaired > 0)
    maxi = 1;

  sum = 0.0;
  for (i=0; i<maxi; i++)
    {
      A = BRA(pair, i);
      B = KET(pair, i);
      sum += element(quant, unit, bra, ket, opts);
    }
  return(sum * pair->sign);
}



//--------------------------------------------------------------------
double ope(quantStruct *quant, unitStruct *unit,
	   long *bra, long *ket, long *opts)
{
  double val;

  if NEQ(MSA, MSB)
    return 0.0;

  val = 0.0;
  if (((MLA == 4) && (MLB == 4)) ||
      ((MLA == 4) && (MLB == -4)) ||
      ((MLA == -4) && (MLB == 4)) ||
      ((MLA == -4) && (MLB == -4)))
    val = 0.5;
  if ((MLA == 0) && (MLB == 0))
      val = 1.0/sqrt(2.0);

  return val;
}



//--------------------------------------------------------------------
#define K  opts[0]
#define Q  opts[1]
double one_u(quantStruct *quant, unitStruct *unit,
	     long *bra, long *ket, long *opts)
{
  //printf("%i:%i  %i (%i, %i), %f\n", A, B, MLA-MLB, K, Q, UAB(K));
  if(NEQ(MSA, MSB) ||
     Q != MLA-MLB)
    return 0.0;
  
  //printf("done\n");
  return UAB(K);
}
#undef K
#undef Q



//--------------------------------------------------------------------
#define K  opts[0]
#define Q  opts[1]
double one_t(quantStruct *quant, unitStruct *unit,
	     long *bra, long *ket, long *opts)
{
  if(NEQ(MLA, MLB) ||
     Q != MSA-MSB)
    return 0.0;
  
  return TAB(K);
}
#undef K
#undef Q



//--------------------------------------------------------------------
#define K  opts[0]
double one_uu(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MLA, MLB) ||
     NEQ(MSA, MSB))
    return 0.0;

  return 1.0 / (L + 1.0);
  // return SIGNEXP2(L-MLA) / (L + 1.0);
}
#undef K



//--------------------------------------------------------------------
#define K1 opts[0]
#define K2 opts[1]
#define Q  opts[2]
double one_uus(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MLA, MLB) ||
     NEQ(MSA, MSB))
    return 0.0;

  return SIGNEXP2(Q) * \
    wign3j(L,K1,L,Q-MLB,-Q,MLB) * \
    wign3j(L,K2,L,Q-MLB,-Q,MLB);
}
#undef K1
#undef K2
#undef Q



//--------------------------------------------------------------------
#define K  opts[0]
double one_tt(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MSA, MSB) ||
     NEQ(MLA, MLB))
    return 0.0;

  return 0.5;
  // return SIGNEXP2(1-MSA) / 2.0;
}
#undef K



//--------------------------------------------------------------------
#define K  opts[0]
double one_ut(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts)
{
  if(NEQ(MLA+MSA, MLB+MSB))
    return 0.0;

  return SIGNEXP2(MLA-MLB) * UAB(K) * TAB(K);
}
#undef K
