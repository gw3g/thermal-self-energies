#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

/*-----------------------------------------------------------------------------------------------*/
/*                                                                "global" over ALL source files */
extern size_t calls;
extern double   tol;
extern double     g;
extern    int   HTL;
/*-----------------------------------------------------------------------------------------------*/

typedef enum p_type {
  B,  // Boson
  F,  // Fermion
} p_type;

typedef enum pol {
  T,  // transverse
  L,  // longitudinal
} pol;

double 
   f(float e, p_type X),    // f(e)       for X type
  bf(float e, p_type X);    // \bar{f}(x) for X type

double *Pi_qed(double complex ,double complex ,pol);
double *Pi_qcd(double complex ,double complex ,pol);
double *Sig_qcd(double,double,pol);

double *frakJ(double,void *,int);

/* htl.c */
double   *Pi_htl(double complex,pol);   // photon
double  *Sig_htl(double complex,pol);          // quarks
double  *T_htl(double complex, double complex, pol);          // quarks

/* disp rel */
double disp(double,pol);

/* 4-vector */
struct pair {double complex o; double complex q;};
struct Qpol {double complex o; double complex q; pol X;};

