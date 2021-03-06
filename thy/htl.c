#include "core.h"

/*
 * HTL gluon self-energy
 *
 */
double *Pi_htl(double complex z, pol X) {
  /*  z := omega/q  */
  double complex z2=z*z, g2=g*g;

  double complex 
    LD  = clog( (z+1.)/(z-1.) )-0.*I*M_PI, 
    P00 = 1. - 0.5*z*LD,            // Pi_00 =: Pi_L
    Pii = z2 + (1.-z2)*0.5*z*LD;    // Pi_T
    Pii/= 2.; P00*=-1.;             // for consistency w/ 
    Pii*= 3.; P00*= 3.;             // Eq(5) of hep-ph/9708434

  double *P = (double *)malloc( 2*sizeof(double) );

  switch(X) {
    case T:                         // transverse = spatial
      P[0]  = g2*creal( Pii );
      P[1]  = g2*cimag( Pii );         break;
    case L:                         // longitudinal = temporal
      P[0]  = g2*creal( P00 );
      P[1]  = g2*cimag( P00 );         break;
  }

  return P;
}

/*
 * HTL quark self-energy
 *
 */
double *Sig_htl(double complex z, pol X) {
  /*  z := omega/q  */

  double complex 
    LD  = clog((z+1.)/(z-1.)), 
    S0  = 0.5*LD ,                  // Sigma_0
    Si  = -1. + 0.5*z*LD;           // Sigma_i   ( dir. ~ \hat{q} )
                                    // cf Eq.(22,26) of hep-ph/9708434
                                    // T_1 ~ S_0
                                    // T_2 ~ S_0 - Si

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


double *T_htl(double complex o, double complex q, pol X) {
  double complex z = o/q;                                     //  no factor of mf^2 

  double complex                               S0 = ( Sig_htl(z,L)[0] + I*Sig_htl(z,L)[1] )/q;
  double complex                               Si = ( Sig_htl(z,T)[0] + I*Sig_htl(z,T)[1] )/q;

  double *Ti = (double *)malloc( 2*sizeof(double) );

  switch(X) {
    case L:                                     // longitudinal = T_1
      Ti[0]  = 4.*creal( S0 );
      Ti[1]  = 4.*cimag( S0 );                   break;
    case T:                                     // transverse = T_2
      Ti[0]  = 4.*creal( o*S0 - q*Si );
      Ti[1]  = 4.*cimag( o*S0 - q*Si );          break;
  }

  return Ti;
}
