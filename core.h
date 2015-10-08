#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

/*-----------------------------------------------------------------------------------------------*/
/*                                                                "global" over ALL source files */
extern size_t calls;
extern double tol;
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
//double re_PiL(double, void *);
//double im_PiL(double, void *);
//double re_PiT(double, void *);
//double im_PiT(double, void *);

double *PI_qed(double,double,pol);
double *PI_qcd(double,double,pol);
double reL(double,double);

double *frakJ(double,void *,int);

/* htl.c */
double   *Pi_htl(double,pol);  // photon
double  *Sig_htl(double,pol);  // quarks
//double 
  //re_Pi_L( double complex, double complex),
  //im_Pi_L( double complex, double complex);
/* disp rel */
double omega_g(double,pol);

/* 4-vector */
struct pair {double complex o; double complex q;};
struct Qpol {double  o; double q; pol X;};

