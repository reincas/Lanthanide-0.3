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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



//--------------------------------------------------------------------
//  General remark:
//    As integer data type always "long" is used and as floating
//    point data type always "double". Don't use int or float!
//--------------------------------------------------------------------



//--------------------------------------------------------------------
// Some definitions.
//--------------------------------------------------------------------

// global shortcuts
#define SIGN(x)     ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
#define MIN(x,y)    ((x) < (y) ? (x) : (y))
#define MAX(x,y)    ((x) > (y) ? (x) : (y))
#define DELTA(x,y)  ((x) == (y) ? 1 : 0)
#define NEQ(x,y)    ((x) != (y))
#define DTOL(x)     (long)((x) + ((x) > 0 ? +0.1 : -0.1))
#define DTO2L(x)    (long)(2*(x) + ((x) > 0 ? +0.1 : -0.1))
#define SIGNEXP2(x) (double)((abs((x)/2) % 2) == 0 ? 1.0 : -1.0)
#define DZERO(x)    (fabs(x) < 1e-12)
#define CZERO(x)    (x.r*x.r + x.i*x.i < 1e-24)
#define SWAP(x,y)   {register long t; t = x; x = y; y = t;}
#define ROT(x,y,z)  {register long t; t = x; x = y; y = z; z = t;}


// Estimated number of nonzero matrix elements = states * MEMFACTOR.
// 32 seems to be a realistic value, but may have to be adjusted for
// better performance. calcmatrix() does *not* crash when the
// estimated value is too small, but a time consuming memory
// reallocation is done in that case.
#define MEMFACTOR 32



//--------------------------------------------------------------------
// Useful structures.
//--------------------------------------------------------------------

// Structure keeping all combinations of quantum numbers ml and ms of
// a electron with a given quantum number l.
// 
//   l: Quantum number l of the electrons.
//   max: Maximum number of different electrons (= 4l+2).
//   ml[max]: Quantum number ml for all different electrons.
//   ms[max]: Quantum number ms for all different electrons.
//   states: Number of states = binom(max, electrons).
//   ownbuf: A value of 1 indicates that the memory blocks for all
//       internal structures were allocated by the creator of the
//       structure.
typedef struct {
  long l;
  long max;
  long *ml, *ms;
  long ownbuf;
} quantStruct;
#define ML(var,i) var->ml[i]
#define MS(var,i) var->ms[i]


// Structure keeping precalculated values of all matrix elements of
// the unit tensor operators u(k) and t(k) for a given orbital quantum
// number l.
//
//   l: Quantum number l of the electrons.
//   max: Number of combinations (= 4l+2).
//   u[(2l+1)*max*max]: Precalculated matrix elements of the orbital
//       unit tensor operator. u[(k*max+i)*max+j] = <ml[i]|u(k)q|ml[j]>
//       and q = ml[i]-ml[j].
//   t[2*max*max]: Precalculated matrix elements of the spin unit
//       tensor operator. t[(k*max+i)*max+j] = <ms[i]|t(k)q|ms[j]> and
//       q = ms[i]-ms[j].
//   ownbuf: A value of 1 indicates that the memory blocks for all
//       internal structures were allocated by the creator of the
//       structure.
typedef struct {
  long l;
  long max;
  double *u, *t;
  long ownbuf;
} unitStruct;
#define U(var,k,i,j) var->u[(((k)/2)*var->max+(i))*var->max+(j)]
#define T(var,k,i,j) var->t[(((k)/2)*var->max+(i))*var->max+(j)]
#define UVAL(var,k,i,j) ((k) <= 2*var->l ? U(var,k,i,j) : 0.0)
#define TVAL(var,k,i,j) ((k) <= 2*var->l ? T(var,k,i,j) : 0.0)


// Structure keeping all combinations of different electrons in a
// given configuration, the determinantal product states.
// 
//   l: Quantum number l of the electrons.
//   electrons: Number of electrons.
//   states: Number of states (= binom(4l+2, electrons)).
//   state[states*electrons]: List of the index numbers of all
//       electrons for every state. The index number of electron k in
//       state i is states[electrons*i+k].
//   ownbuf: A value of 1 indicates that the memory blocks for all
//       internal structures were allocated by the creator of the
//       structure.
typedef struct {
  long l, electrons;
  long states;
  long *state;
  long ownbuf;
} configStruct;
#define STATE(var,i,k) var->state[electrons*(i)+(k)]


// Structure keeping a pair of ordered determinantal product states.
//
//   i, j: Indices of the two states.
//   electrons: Number of electrons.
//   statei[electrons]: List of the index numbers of all electrons of
//       state i.
//   statej[electrons]: List of the index numbers of all electrons of
//       state j.
//   unpaired: Number of unpaired electrons. These are the first ones
//       in the lists statei[] and statej[].
//   sign: A value -1/+1 indicates an odd/even permutation of
//       electrons to pair all possible electrons in statei[] and
//       statej[].
//   ownbuf: A value of 1 indicates that the memory blocks for all
//       internal structures were allocated by the creator of the
//       structure.
typedef struct {
  long i, j;
  long electrons;
  long *statei, *statej;
  long unpaired, sign;
  long ownbuf;
} pairStruct;
#define BRA(val,k) val->statei[k]
#define KET(val,k) val->statej[k]


// Structure keeping the non-zero elements of a spare symetric or
// hermitean matrix.
//
//   states: Number of quantum states. The dimension of the matrix is
//       therefore states x states.

//   size: Length of matrix element in sizeof(double). In case of a
//       symetric / hermitean matrix, the elements are real / complex,
//       therefore size = 1 / 2.
//   data[]: Values of all nonzero matrix elements in the upper right
//       triangle of the matrix. Length of this vector is size *
//       number of nonzero matrix elements.
//   col[]: Column indices i of the appropriate nonzero matrix
//       elements. Length of this vector is therefore the number of
//       nonzero matrix elements.
//   row[]: row index pointer to col: the nonzero matrix elements in
//       the row j have col indices in the range row[j-1] until
//       row[j]-1. The first nonzero matrix element in row j has
//       therefore the column index i=col[row[j-1]] and the value
//       data[row[j-1]*size]. The case row[j]==row[j-1] indicates that
//       row j contains no nonzero matrix elements in the upper right
//       triangle.
//   ownbuf: A value of 1 indicates that the memory blocks for all
//       internal structures were allocated by the creator of the
//       structure.
//
// Number of nonzero matrix elements is row[states-1].
typedef struct
{
  long states;
  long size;
  double *data;
  long *col;
  long *row;
  long ownbuf;
} sparseStruct;


// Type definition for the data type "complex".
typedef struct {
  double r, i;
} complex;



//--------------------------------------------------------------------
// Calculation of matrix elements.
//--------------------------------------------------------------------

// Function prototype for the matrix element function type. This
// prototype is necessary to keep the functions in a structure like
// elementStruct.
//
//   quant: Quantum numbers of the electrons.
//   unit: Precalculated matrix elements of the unit tensor operators.
//   bra[]: List of quantum numbers of the bra state.
//   ket[]: List of quantum numbers of the ket state.
//   opts[]: List of option values like order of tensor operators.
//
// Returned is the value of the matrix element.
typedef double (*elementFunc)(quantStruct *quant, unitStruct *unit,
			      long *bra, long *ket, long *opts);


// The matrix element functions:
double one_u(quantStruct *quant, unitStruct *unit,
	     long *bra, long *ket, long *opts);
double one_t(quantStruct *quant, unitStruct *unit,
	     long *bra, long *ket, long *opts);
double one_uu(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts);
double one_uus(quantStruct *quant, unitStruct *unit,
	       long *bra, long *ket, long *opts);
double one_tt(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts);
double one_ut(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts);
double two_uu(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts);
double two_uus(quantStruct *quant, unitStruct *unit,
	       long *bra, long *ket, long *opts);
double two_tt(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts);
double two_ut(quantStruct *quant, unitStruct *unit,
	      long *bra, long *ket, long *opts);
double two_uutt(quantStruct *quant, unitStruct *unit,
		long *bra, long *ket, long *opts);
double three_uuua(quantStruct *quant, unitStruct *unit,
		  long *bra, long *ket, long *opts);
double three_uuub(quantStruct *quant, unitStruct *unit,
		  long *bra, long *ket, long *opts);
double three_uuuc(quantStruct *quant, unitStruct *unit,
		  long *bra, long *ket, long *opts);
double ope(quantStruct *quant, unitStruct *unit,
	     long *bra, long *ket, long *opts);


// Structure keeping fixed informations on all matrix element.
//
//   electrons: Minimum number of electrons.
//   opts: Number of option values.
//   swap: Number of different operators (see matrix-xxx.c).
//   element(): Function for calculating the matrix element.
typedef struct {
  long electrons, opts, swap;
  elementFunc element;
} elementStruct;


// List of elementStruct structures to keep fixed information for all
// matrix elements at one place.
static elementStruct ELEMENT[] = {
  { 1, 2, 1, one_u      },
  { 1, 2, 1, one_t      },
  { 1, 1, 1, one_uu     },
  { 1, 3, 1, one_uus    },
  { 1, 1, 1, one_tt     },
  { 1, 1, 1, one_ut     },
  { 2, 1, 1, two_uu     },
  { 2, 3, 2, two_uus    },
  { 2, 1, 1, two_tt     },
  { 2, 1, 2, two_ut     },
  { 2, 5, 2, two_uutt   },
  { 3, 1, 1, three_uuua },
  { 3, 2, 2, three_uuub },
  { 3, 3, 3, three_uuuc },
  { 1, 2, 1, ope        },
};


// Indices for all matrix elements, as given in the list ELEMENT.
#define MAT1_U     0
#define MAT1_T     1
#define MAT1_UU    2
#define MAT1_UUS   3
#define MAT1_TT    4
#define MAT1_UT    5
#define MAT2_UU    6
#define MAT2_UUS   7
#define MAT2_TT    8
#define MAT2_UT    9
#define MAT2_UUTT 10
#define MAT3_UUUA 11
#define MAT3_UUUB 12
#define MAT3_UUUC 13
#define MAT1_OPE  14



//--------------------------------------------------------------------
// C functions for general use.
//--------------------------------------------------------------------

// Calculate a binomial coefficient.
long binom(long n, long k);


// Free all memory blocks allocated in the different structures.
void freeQuant(quantStruct *quant);
void freeUnit(unitStruct *unit);
void freeConfig(configStruct *config);
void freePair(pairStruct *pair);
void freeSparse(sparseStruct *matrix);


// Build the structure quant with all permutations of the quantum
// numbers ml and ms for electrons with the given quantum number l.
// A return value != 0 indicates an error.
long buildquant(long l, quantStruct *quant);


// Build the structure unit with all matrix elements of the unit
// tensor operators u(k) and t(k) for the given qunatum number l.
// A return value != 0 indicates an error.
long buildunit(quantStruct *quant, unitStruct *unit);


// Build the structure config with all permutations of different
// electrons with quantum numbers as given in quant. electrons is the
// number of equivalent electrons.
// A return value != 0 indicates an error.
long buildconfig(quantStruct *quant, long electrons, configStruct *config);


// Allocate memory as needed for the configuration config and set all
// fixed variables in the structure pair.
// A return value != 0 indicates an error.
long initpair(configStruct *config, pairStruct *pair);


// Build up the structure pair by ordering the electrons of states i
// and j in the configuration config. If pair already contains data
// for a pair of states with the same number of electrons as in
// config, the structure is just filled with new data. Otherwise it is
// free'd first.
// A return value != 0 indicates an error.
long orderpair(configStruct *config, long i, long j, pairStruct *pair);


// Build up the upper right triangle of a symetric matrix.
// ELEMENT[key] gives all information necessary to calculate the
// matrix elements in the configuration config and with the options
// opts. A return value != 0 indicates an error.
long calcmatrix(quantStruct *quant, unitStruct *unit, configStruct *config,
		long key, long *opts, sparseStruct *matrix);


// Calculate the value(s) of a matrix element with the function
// ELEMENT[key].element(). pair gives the ordered pair of states and
// opts the necessary options. The value(s) are saved in the buffer
// data. A return value != 0 indicates an error.
long calcelement(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		 long key, long *opts, double *data);


// Basic matrix element functions for one-, two-, and three-electron
// operators.
double elementOne(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		  long key, long *opts);
double elementTwo(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		  long key, long *opts);
double elementThree(quantStruct *quant, unitStruct *unit, pairStruct *pair,
		    long key, long *opts);



//--------------------------------------------------------------------
// Functions to calculate eigenvalues und eigenvectors.
//--------------------------------------------------------------------

#define DIAG_FAST   0
#define DIAG_SMALL  1
#define DIAG_NOVEC  0
#define DIAG_VEC    1


// Diagonalize the symetric matrix ap by using the LAPACK functions
// dspevd or dspev. The arguments have the following meaning:
//
//   method: Set to DIAG_FAST / DIAG_SMALL to use the CLAPACK
//       functions dspevd or dspev respectively.
//   vecflag: Set to DIAG_NOVEC / DIAG_VEC to calculate the
//       eigenvalues only or also the eigenvectors.
//   n: The order of the matrix. (n x n)
//   ap[]: Pointer to the matrix elements. For an upper right triangle
//       a[i,j] = ap[i + j*(j+1)/2] and for a lower left triangle
//       a[i,j] = ap[i + j*(2*n-j-1)/2]. So in both cases the elements
//       are stored in ascending columns as usual, but with all matrix
//       elements not in the appropriate triangle left out.
//   w[]: Vector of the eigenvalues in ascending order.
//   z[]: Array of the eigenvectors. The eigenvector to the eigenvalue
//       w[i] starts at z[i*n].
//
// Memory requirement for DIAG_FAST:  
//   DIAG_VEC:    20*n*n +  80*n + 20  Bytes
//   DIAG_NOVEC:   4*n*n +  28*n +  4  Bytes
//
// Memory requirement for DIAG_SMALL:
//   DIAG_VEC:    12*n*n +  36*n  Bytes
//   DIAG_NOVEC:   4*n*n +  44*n  Bytes
long diagsym(long method, long vecflag, long n,
	     double *ap, double *w, double *z);


// Diagonalize the hermitean matrix ap by using the LAPACK functions
// zhpevd or zhpev. The arguments have the following meaning:
//
//   method: Set to DIAG_FAST / DIAG_SMALL to use the CLAPACK
//       functions zhpevd or zhpev respectively.
//   vecflag: Set to DIAG_NOVEC / DIAG_VEC to calculate the
//       eigenvalues only or also the eigenvectors.
//   n: The order of the matrix. (n x n)
//   ap[]: Pointer to the matrix elements. For an upper right triangle
//       a[i,j] = ap[i + j*(j+1)/2] and for a lower left triangle
//       a[i,j] = ap[i + j*(2*n-j-1)/2]. So in both cases the elements
//       are stored in ascending columns as usual, but with all matrix
//       elements not in the appropriate triangle left out.
//   w[]: Vector of the eigenvalues in ascending order.
//   z[]: Array of the eigenvectors. The eigenvector to the eigenvalue
//       w[i] starts at z[i*n].
//
// Memory requirement for DIAG_FAST:  
//   DIAG_VEC:    40*n*n + 108*n + 20  Bytes
//   DIAG_NOVEC:   8*n*n +  56*n +  4  Bytes
//
// Memory requirement for DIAG_SMALL:
//   DIAG_VEC:    24*n*n +  72*n - 32  Bytes
//   DIAG_NOVEC:   8*n*n +  88*n - 32  Bytes
long diagherm(long method, long vecflag, long n,
	      complex *ap, double *w, complex *z);



//--------------------------------------------------------------------
// Functions to handle sparse matrices.
//--------------------------------------------------------------------

// Save the value(s) of the matrix element matrix[i][j] in the buffer
// data. A return value 0/1 indicates a nonzero/zero matrix element.
double realgetelement(sparseStruct *matrix, long i, long j);
complex cplxgetelement(sparseStruct *matrix, long i, long j);


// Add the sparse matrix multiplied by factor to the upper triangle
// matrix dest. If factor == 0.0, nothing is done. A return value != 0
// indicates an error.
long realmultadd(sparseStruct *matrix, double factor, double *dest);
long cplxmultadd(sparseStruct *matrix, complex factor, complex *dest);


// Transform the sparse matrix by the array of eigenfunctions trans.
// The result is stored as a new sparse matrix. A return value != 0
// indicates an error.
long realtransform(sparseStruct *matrix, double *trans,
		   sparseStruct *result);
long cplxtransform(sparseStruct *matrix, complex *trans,
		   sparseStruct *result);


// Transform the sparse matrix by the array of eigenfunctions trans
// and store the diagonal elements of the result. A return value != 0
// indicates an error.
long realtransdiag(sparseStruct *matrix, double *trans,
		   double *result, long terms);
long cplxtransdiag(sparseStruct *matrix, complex *trans,
		   complex *result, long terms);



//--------------------------------------------------------------------
// Functions to handle triangle matrices.
//--------------------------------------------------------------------

// Copy the states x states triangle matrix tri with rows and columns
// from substart until subend-1 to result.
long subtriangle(double *tri, long states, long substart, long subend,
		 double *result);



//--------------------------------------------------------------------
// Wigner n-j symbols and others.
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
