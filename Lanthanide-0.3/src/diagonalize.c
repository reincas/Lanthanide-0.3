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


// external CLAPACK functions
extern int dspevd_(char *jobz, char *uplo, long *n,
		   double *ap, double *w, double *z, long *ldz,
		   double *work, long *lwork,
		   long *iwork, long *liwork,
		   long *info);

extern int dspev_(char *jobz, char *uplo, long *n,
		  double *ap, double *w, double *z, long *ldz,
		  double *work,
		  long *info);

extern int zhpevd_(char *jobz, char *uplo, long *n,
		   complex *ap, double *w, complex *z, long *ldz,
		   complex *work, long *lwork,
		   double *rwork, long *lrwork,
		   long *iwork, long *liwork,
		   long *info);

extern int zhpev_(char *jobz, char *uplo, long *n,
		  complex *ap, double *w, complex *z, long *ldz,
		  complex *work, double *rwork,
		  long *info);



//--------------------------------------------------------------------
long diagsym(long method, long vecflag, long n,
	     double *ap, double *w, double *z)
{
  double *work, worklen;
  long *iwork, iworklen;
  long ldz, lwork, liwork;
  long info;
  char jobz, uplo;
			
  //printf("diagherm: prepare jobz, uplo, and ldz\n");
  if (vecflag == DIAG_NOVEC)
    {
      jobz = 'N';
      ldz  = 1;
    }
  else
    {
      jobz = 'V';
      ldz  = n;
    }
  uplo = 'U';

  if (method == DIAG_FAST)
    {
      //printf("diagsym: ask for workspace size\n");
      lwork = -1;
      liwork = -1;
      dspevd_(&jobz, &uplo, &n, ap, w, z, &ldz,
	      &worklen, &lwork, &iworklen, &liwork, &info);
      if (info < 0)
	return info;

      //printf("diagsym: prepare work\n");
      lwork = (long)worklen;
      if(!(work = malloc(lwork * sizeof(double))))
	return -102;

      //printf("diagsym: prepare iwork\n");
      liwork = (int)iworklen;
      if(!(iwork = malloc(liwork * sizeof(long))))
	return -103;

      //printf("diagsym: fast diagonalisation\n");
      dspevd_(&jobz, &uplo, &n, ap, w, z, &ldz,
	      work, &lwork, iwork, &liwork, &info);
    }
  else
    {
      //printf("diagsym: prepare work\n");
      lwork = 3*n;
      if(!(work = malloc(lwork * sizeof(double))))
	return -104;

      //printf("diagsym: small diagonalisation\n");
      dspev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);
    }
  
  //printf("diagsym: cleanup\n");
  free(work);
  if (method == DIAG_FAST)
    free(iwork);
  return info;
}



//--------------------------------------------------------------------
long diagherm(long method, long vecflag, long n,
	      complex *ap, double *w, complex *z)
{
  complex *work, worklen;
  double *rwork, rworklen;
  long *iwork, iworklen;
  long ldz, lwork, lrwork, liwork;
  long info;
  char jobz, uplo;
			
  //printf("diagherm: prepare jobz, uplo, and ldz\n");
  if (vecflag == DIAG_NOVEC)
    {
      jobz = 'N';
      ldz  = 1;
    }
  else
    {
      jobz = 'V';
      ldz  = n;
    }
  uplo = 'U';

  if (method == DIAG_FAST)
    {
      //printf("diagherm: ask for workspace size\n");
      lwork = -1;
      lrwork = -1;
      liwork = -1;
      zhpevd_(&jobz, &uplo, &n, ap, w, z, &ldz,
	      &worklen, &lwork, &rworklen, &lrwork, &iworklen, &liwork,
	      &info);
      if (info < 0)
	return info;

      //printf("diagherm: prepare work\n");
      lwork = (long)worklen.r;
      if(!(work = malloc(lwork * sizeof(complex))))
	return -100;
      
      //printf("diagherm: prepare rwork\n");
      lrwork = (long)rworklen;
      if(!(rwork = malloc(lrwork * sizeof(double))))
	return -101;
      
      //printf("diagherm: prepare iwork\n");
      liwork = iworklen;
      if(!(iwork = malloc(liwork * sizeof(long))))
	return -102;

      //printf("diagherm: fast diagonalisation\n");
      zhpevd_(&jobz, &uplo, &n, ap, w, z, &ldz,
	      work, &lwork, rwork, &lrwork, iwork, &liwork,
	      &info);
    }
  else
    {
      //printf("diagherm: prepare work\n");
      lwork = 2*n - 1;
      if(!(work = malloc(lwork * sizeof(complex))))
	return -103;
      
      //printf("diagherm: prepare rwork\n");
      lrwork = 2*n - 2;
      if(!(rwork = malloc(lrwork * sizeof(double))))
	return -104;
      
      //printf("diagherm: small diagonalisation\n");
      zhpev_(&jobz, &uplo, &n, ap, w, z, &ldz,
	     work, rwork, &info);
    }

  free(work);
  free(rwork);
  if (method == DIAG_FAST)
    free(iwork);
  return info;
}
