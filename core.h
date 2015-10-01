#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

typedef enum p_type {
  B,  // Boson
  F,  // Fermion
} p_type;

typedef enum pol {
  T,  // transverse
  L,  // longitudinal
} pol;

double re_PiL(double, void *);
double im_PiL(double, void *);

/* htl.c */
double   *Pi_htl(double,pol);  // photon
double  *Sig_htl(double,pol);  // quarks
double 
  re_Pi_L( double complex, double complex),
  im_Pi_L( double complex, double complex);

/* 4-vector */
struct pair {double complex o; double complex q;};

