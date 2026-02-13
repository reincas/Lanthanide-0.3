//
// (c) 08/1974 by Harwell Library, authors: T. Lindelof and F. James
// (c) 01/2001 by Reinhard Caspary <r.caspary@tu-bs.de>
//
// The code in this file is based on the FORTRAN program clebsg.f from
// T. Lindelof and F. James. Translation to C by Reinhard Caspary.
// 

#include <math.h>


//--------------------------------------------------------------------
// Prototypes of functions to export
//--------------------------------------------------------------------

// This function must be called before the following functions of
// wigner.c may be used. It calculates a list of factorials to speed
// up the calculation of the functions.
void initwigner();

// Calculate a Wigner 3j symbol. All arguments are integers and must
// be given as *doubled* values!
double wign3j(long k1, long k2, long k3,
	      long k4, long k5, long k6);

// Calculate a Clebsch-Gordan coefficient. All arguments are integers
// and must be given as *doubled* values!
double clebsg(long k1, long k2, long k3,
	      long k4, long k5, long k6);

// Calculate a Wigner 6j symbol. All arguments are integers and must
// be given as *doubled* values!
double wign6j(long k1, long k2, long k3,
	      long k4, long k5, long k6);

// Calculate a Racah coefficient. All arguments are integers and must
// be given as *doubled* values!
double racahc(long k1, long k2, long k5,
	      long k4, long k3, long k6);

// Calculate a U-function of Jahn. All arguments are integers and must
// be given as *doubled* values!
double jahnuf(long k1, long k2, long k5,
	      long k4, long k3, long k6);

// Calculate a Wigner 9j symbol. All arguments are integers and must
// be given as *doubled* values!
double wign9j(long k1, long k2, long k3,
	      long k4, long k5, long k6,
	      long k7, long k8, long k9);



//--------------------------------------------------------------------
// Prototypes of private functions
double _wign3js(long k1, long k2, long k3);
long _triangle(long ka, long kb, long kc, 
	      double *ay, long *iay);



//--------------------------------------------------------------------
// key values in clebsg.f
// Not used any more, but given here for historic reasons and for
// anybody who wants to compare the code to its original.
#define WIGN3J   1
#define CLEBSG   2
#define WIGN3Js  3
#define CLEBSGs  4
#define WIGN6J  11
#define RACAHC  12
#define JAHNUF  13
#define WIGN9J -10

#define NUMFACT 100
#define IPARF(x)  (long)(4*((x)/4)-(x)+1)
#define SIGN(x)   ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)



//--------------------------------------------------------------------
// Global parameters for private use in this file only. Filled by
// initwigner() and used by all others.
static double h[NUMFACT+2];
static long j[NUMFACT+2];



//--------------------------------------------------------------------
double wign3j(long k1, long k2, long k3,
	      long k4, long k5, long k6)
{
  if ((k4 == 0) && (k5 == 0) && (k6 == 0))
    return(_wign3js(k1, k2, k3));

  return(IPARF(abs(k1 - k2 - k6))
	 * clebsg(k1, k2, k3, k4, k5, -k6)
	 / sqrt(k3+1.0));
}



//--------------------------------------------------------------------
double _wign3js(long k1, long k2, long k3)
{
  long m1, m2, m3, m4;
  long ij, ijpar;
  double y, z;
  long iy, iz;
  
  if(((k1 % 2) != 0) || ((k2 % 2) != 0) || ((k3 % 2) != 0))
    return(0.0);
  
  ij = k1 + k2 + k3;
  if(IPARF(ij) <= 0)
    return(0.0);
  
  m1 = ij - k1 - k1;
  m2 = ij - k2 - k2;
  m3 = ij - k3 - k3;
  m4 = ij + 2;
  if((m1 < 0) || (m2 < 0) || (m3 < 0))
    return(0.0);
  
  m1 = m1/2 + 1;
  m2 = m2/2 + 1;
  m3 = m3/2 + 1;
  m4 = ij/2 + 2;
  
  y = sqrt(h[m1] * h[m2] * h[m3] / h[m4]);
  iy = (j[m1] + j[m2] + j[m3] - j[m4]) / 2;
  
  ij = ij/2;
  ijpar = IPARF(ij);
  ij = ij/2 + 1;
  
  m1 = m1/2 + 1;
  m2 = m2/2 + 1;
  m3 = m3/2 + 1;
  
  z = h[ij] / (h[m1] * h[m2] * h[m3]);
  iz = j[ij] - j[m1] - j[m2] - j[m3];
  iz = iz + iy;
  
  return(ijpar * y * z * pow(10.0, iz));
};



//--------------------------------------------------------------------
double clebsg(long k1, long k2, long k3,
	      long k4, long k5, long k6)
{
  long m1, m2, m3, m4, m5, m6, m7, m8, m9, m10;
  long mm1, mm2, mm3, mm4, mm5;
  long n4, n5, n5par;
  double x, y, z;
  long ix, iy;
  
  if ((k4 == 0) && (k5 == 0) && (k6 == 0))
    return(_wign3js(k1, k2, k3) * IPARF(abs(k2 - k1)) * sqrt(k3+1.0));
  
  if((k4 + k5 - k6) != 0)
    return(0.0);
  
  m1 = k1 + k2 - k3;
  m2 = k2 + k3 - k1;
  m3 = k3 + k1 - k2;
  m4 = k1 + k4;
  m5 = k1 - k4;
  m6 = k2 + k5;
  m7 = k2 - k5;
  m8 = k3 + k6;
  m9 = k3 - k6;
  m10 = k1 + k2 + k3 + 2;
  
  if((m1<0) || (m2<0) || (m3<0) || (m4<0) || 
     (m5<0) || (m6<0) || (m7<0) || (m8<0) || (m9<0))
    return(0.0);
  
  if(((m4 % 2) != 0) ||
     ((m6 % 2) != 0) ||
     ((m8 % 2) != 0) ||
     ((m10 % 2) != 0))
    return(0.0);
  
  y = k3 + 1.0;
  
  m1 = m1/2 + 1;
  m2 = m2/2 + 1;
  m3 = m3/2 + 1;
  m4 = m4/2 + 1;
  m5 = m5/2 + 1;
  m6 = m6/2 + 1;
  m7 = m7/2 + 1;
  m8 = m8/2 + 1;
  m9 = m9/2 + 1;
  m10 = m10/2 + 1;
  
  y = sqrt(y * h[m1] * h[m2] * h[m3] * h[m4] * h[m5] *
	   h[m6] * h[m7] * h[m8] * h[m9] / h[m10]);
  iy = (j[m1] + j[m2] + j[m3] + j[m4] + j[m5] +
	j[m6] + j[m7] + j[m8] + j[m9] - j[m10]) / 2;
  
  n4 = m1;
  if(n4 > m5)
    n4 = m5;
  if(n4 > m6)
    n4 = m6;
  n4 = n4 - 1;
  m2 = k2 - k3 - k4;
  m3 = k1 + k5 - k3;
  n5 = 0;
  if(n5 < m2)
    n5 = m2;
  if(n5 < m3)
    n5 = m3;
  n5par = IPARF(n5);
  n5 = n5/2;
  
  z = 0.0;
  while(SIGN(n5-n4) < 1)
    {
      mm1 = m1 - n5;
      mm2 = m5 - n5;
      mm3 = m6 - n5;
      mm4 = n5 - (m2/2) + 1;
      mm5 = n5 - (m3/2) + 1;
      
      x = 1.0 / (h[mm1] * h[mm2] * h[mm3] * h[mm4] * h[mm5] * h[n5+1]);
      ix = -j[mm1] - j[mm2] - j[mm3] - j[mm4] - j[mm5] - j[n5+1];
      
      while(SIGN(ix + iy) != 0)
	{
	  if(SIGN(ix + iy) < 0)
	    {
	      x = 0.1 * x;
	      ix = ix + 1;
	    }
	  else
	    {
	      x = 10.0 * x;
	      ix = ix - 1;
	    }
	}
      
      if(n5par < 0)
	x = -x;
      z = z + x;
      n5par = -n5par;
      n5 = n5 + 1;
    }
  
  return(z * y);
}



//--------------------------------------------------------------------
double wign6j(long k1, long k2, long k3,
	      long k4, long k5, long k6)
{
  long ijpar;
  
  ijpar = IPARF(k1 + k2 + k4 + k5);
  if(ijpar < 0)
    return(-racahc(k1, k2, k5, k4, k3, k6));
  return(racahc(k1, k2, k5, k4, k3, k6));
}



//--------------------------------------------------------------------
double racahc(long k1, long k2, long k5,
	      long k4, long k3, long k6)
{
  long m1, m2, m3, m4, m5, m6, m7;
  long n4, n5, n5par;
  double x, y, z;
  long ix, iy;
  
  y = 1.0;
  iy = 0;
  if(_triangle(k1, k2, k3, &y, &iy) == 0)
    return(0.0);
  if(_triangle(k4, k5, k3, &y, &iy) == 0)
    return(0.0);
  if(_triangle(k4, k2, k6, &y, &iy) == 0)
    return(0.0);
  if(_triangle(k1, k5, k6, &y, &iy) == 0)
    return(0.0);
  
  m1 = (k1 + k2 + k4 + k5)/2 + 2;
  m2 = (k1 + k2 - k3)/2 + 1;
  m3 = (k4 + k5 - k3)/2 + 1;
  m4 = (k1 + k5 - k6)/2 + 1;
  m5 = (k2 + k4 - k6)/2 + 1;
  m6 = k1 + k4 - k3 - k6;
  m7 = k2 + k5 - k3 - k6;
  
  n4 = m1;
  if(n4 > m2)
    n4 = m2;
  if(n4 > m3)
    n4 = m3;
  if(n4 > m4)
    n4 = m4;
  if(n4 > m5)
    n4 = m5;
  n4 = n4 - 1;
  
  n5 = 0;
  if(n5 < m6)
    n5 = m6;
  if(n5 < m7)
    n5 = m7;
  n5par = IPARF(n5);
  n5 = n5/2;
  
  m6 = m6/2 - 1;
  m7 = m7/2 - 1;
  
  z = 0.0;
  while(n5 <= n4)
    {
      x = h[m1-n5] / (h[n5+1] * h[m2-n5] * h[m3-n5] * h[m4-n5] *
		      h[m5-n5] * h[n5-m6] * h[n5-m7]);
      ix = j[m1-n5] - j[n5+1] - j[m2-n5] - j[m3-n5] -
	j[m4-n5] - j[m5-n5] - j[n5-m6] - j[n5-m7];
      
      while(SIGN(ix + iy) != 0)
	{
	  if(SIGN(ix + iy) < 0)
	    {
	      x = 0.1 * x;
	      ix = ix + 1;
	    }
	  else
	    {
	      x = 10.0 * x;
	      ix = ix - 1;
	    }
	}
      if(n5par < 0)
	x = -x;
      z = z + x;
      n5par = -n5par;
      n5 = n5 + 1;
    }
  
  return(z * y);
}



//--------------------------------------------------------------------
long _triangle(long ka, long kb, long kc, 
	       double *ay, long *iay)
{
  long ma, mb, mc, md;
  
  ma = ka + kb - kc;
  mb = ka - kb + kc;
  mc = -ka + kb + kc;
  md = ka + kb + kc + 2;
  if((ma < 0) || (mb < 0) || (mc < 0))
    return(0);
  if((md % 2) != 0)
    return(0);
  
  ma = ma/2 + 1;
  mb = mb/2 + 1;
  mc = mc/2 + 1;
  md = md/2 + 1;
  ay[0] = ay[0] * sqrt(h[ma] * h[mb] * h[mc] / h[md]);
  iay[0] = iay[0] + (j[ma] + j[mb] + j[mc] - j[md]) / 2;
  return(1);
}



//--------------------------------------------------------------------
double jahnuf(long k1, long k2, long k5,
	      long k4, long k3, long k6)
{
  return(racahc(k1, k2, k5, k4, k3, k6) * sqrt((k3+1.0) * (k6+1)));
}



//--------------------------------------------------------------------
double wign9j(long k1, long k2, long k3,
	      long k4, long k5, long k6,
	      long k7, long k8, long k9)
{
  long m1, m2;
  long k, kup;
  double ra, rb, rc;
  double result;
  
  kup = k1 + k9;
  m1 = k4 + k8;
  m2 = k2 + k6;
  if(kup > m1)
    kup = m1;
  if(kup > m2)
    kup = m2;
  
  k = k1 - k9;
  if(k < 0)
    k = -k;
  m1 = k4 - k8;
  if(m1 < 0)
    m1 = -m1;
  m2 = k2 - k6;
  if(m2 < 0)
    m2 = -m2;
  if(k < m1)
    k = m1;
  if(k < m2)
    k = m2;
  
  result = 0.0;
  while(k > kup)
    {
      ra = racahc(k1, k4, k9, k8, k7, k);
      if(ra != 0.0)
	{
	  rb = racahc(k2, k8, k6, k4, k5, k);
	  if(rb != 0.0)
	    {
	      rc = racahc(k9, k6, k1, k2, k3, k);
	      if(rc != 0.0)
		result = result + ra*rb*rc*(k+1);
	    }
	}
      k = k + 2;
    }
  return(result);
} 



//--------------------------------------------------------------------
void initwigner()
{
  long i;
  double x;
  
  h[1] = 1.0;
  j[1] = 0;
  x = 0.0;
  
  for(i=2; i<NUMFACT+2; i++)
    {
      x = x + 1.0;
      h[i] = h[i-1] * x;
      j[i] = j[i-1];
      
      while(h[i] >= 10.0)
	{
	  h[i] = 0.01 * h[i];
	  j[i] = j[i] + 2;
	}
    }
}
