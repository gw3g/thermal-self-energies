#include "core.h"
#include <complex.h>
#include <stdio.h>

/*
 * HTL gluon self-energy
 *
 */
double *Pi_htl(double z, pol X) {
  /*  z := omega/q  */
  double z2=z*z;

  double complex 
    LD  = clog((1.+z)/(1.-z))-I*M_PI, 
    P00 = 1. - 0.5*z*LD,            // Pi_00 =: Pi_L
    Pii = z2 + (1.-z2)*0.5*z*LD;    // Pi_T

  double *P = (double *)malloc( 2*sizeof(double) );

  switch(X) {
    case T:                         // transverse = spatial
      P[0]  = creal( Pii );
      P[1]  = cimag( Pii );         break;
    case L:                         // longitudinal = temporal
      P[0]  = creal( P00 );
      P[1]  = cimag( P00 );         break;
  }

  return P;
}

/*
 * HTL quark self-energy
 *
 */
double *Sig_htl(double z, pol X) {
  /*  z := omega/q  */
  double z2=z*z;

  double complex 
    LD  = clog((1.+z)/(1.-z))-I*M_PI, 
    S0  = 0.5*LD ,                  // Sigma_0
    Si  = -1. + 0.5*z*LD;           // Sigma_i   ( dir. ~ \hat{q} )

  double *S = (double *)malloc( 2*sizeof(double) );

  switch(X) {
    case T:                         // transverse = spatial
      S[0]  = creal( Si );
      S[1]  = cimag( Si );          break;
    case L:                         // longitudinal = temporal
      S[0]  = creal( S0 );
      S[1]  = cimag( S0 );          break;
  }

  return S;                         // return array of [ Re(Sigma), Im(Sigma) ]
}

