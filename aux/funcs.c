#include "core.h"

double *fAUX(double complex r, double k, double complex o, double complex q, int i) {

  double complex ll = clog(k+r+o+o*I*1e-9 ),  D;
  double complex r2=r*r, r3=r*r2, r4=r*r3, 
                 q2=q*q, o2=o*o, 
                 k2=k*k, k3=k*k2, ko = k+o;

  double *res     = (double*)malloc(2*sizeof(double));    // 0-REAL, 1-IMAG

  if (i==0) {         // j D

       D   = +  k*(    + r - ko*ll                                   )*q2;           }

  else if (i==1) {    // j (k.r) D

       D   = +.50*(    -.5*r2 + r*ko - (ko*ko + k2 - q2)*ll          )*q2;           }

  else if (i==2) {    // j (q.k) (q.r) D

       D   = -.25*(    + .25*r4 - r3*ko/3. + .5*r2*(ko*ko - 2.*k2)
                            + r*(k3 - k2*o - 3.*k*o2 - o*o2 ) 
                            + ( o2*(k+ko)*(k+ko) - q2*q2 )*ll     ) ;                }

  else if (i==3) {    // D

       D   = +   ll*q2;                                                              }

  else if (i==4) {    //  j (1/k) D

       D   = +   r  - ko*ll;                                                         }

  else if (i==5) {    //  j (1/r) D

       D   = +   k*ll;                                                               }

  else if (i==6) {    //  j ( (q.r)/kr ) D

       D   = + .5*( +.5*r2 - r*ko + ( ko*ko + q2 - k2 )*ll );                        }

  else if (i==7) {    //  j ( (q.k)/kr ) D

       D   = + .5*( -.5*r2 + r*ko - ( ko*ko - q2 - k2 )*ll );                        }


  /*printf(" %.5f + i %.5f \n", creal(D), cimag(D) );*/
  /*printf(" %.5f ; %.5f \n", k, creal(r) );*/
  res[0] = creal(D); res[1] = cimag(D);                                  return res;

}

double *frakJ (double k, void *params, int i) {
  /*                          ---+---
   * k-integrand                 |
   *                          (params) = Qpol {o, q, X};
   *
   * i =    0,    1,    2,    3
   *        J00   Jii  qJij   Jij
   */

  struct Qpol * Q = (struct Qpol *)params;

  double complex o = Q->o;
  double complex q = Q->q;
  pol            X = Q->X;

  double *r_int   = (double*)malloc(2*sizeof(double));        // 0-REAL, 1-IMAG
  double *res     = (double*)malloc(2*sizeof(double));        //

  double sr, so;
  double s = 1.;

                                                              // upper/lower r-limits
  double complex rU = csqrt( (q+k)*(q+k) );
  double complex rL = csqrt( (q-k)*(q-k) );

  res[0]  = 0. ; res[1] = 0.;

  /*
   * carry out the r-integration between (rL, rU).
   */
                                                              // enumerate energy  denoms:
  for (int j=0; j<4; j++) {      sr  = (double) (2*(j%2)-1);  //  -   +   -   +   "sign of r"
                                 so  = (double) (2*(j/2)-1);  //  -   -   +   +   "sign of o"

  if ( (i==4) ) {s = so;}
  if ( (i==5) ) {s = so*sr;}
  /*if ( (i==6)&&(creal(o)>creal(q)) ) {s =  so;}*/
  if ( (i==7) ) {s = sr;}
  /*if ( (i==6)||(i==7) ) {s = so;}*/
  //                    r,    k, omega,   q,``i''
      r_int = fAUX(   sr*rU,  k,  so*o,   q,  i   );      res[0] += s*r_int[0]  ;//     \__ upper
                                                          res[1] += s*r_int[1]  ;//     /
      r_int = fAUX(   sr*rL,  k,  so*o,   q,  i   );      res[0] -= s*r_int[0]  ;//     \__ lower
                                                          res[1] -= s*r_int[1]  ;//     /

  };                                            double complex D; 
                                                               D  = res[0] + I*res[1] ; 
                                                               D *= -.5/q             ;
                                      // (jacobian)
                                                          res[0]  = creal(D)    ; 
                                                          res[1]  = cimag(D)    ;

  /*printf(" %.3f : %.5f + i %.5f \n", k, res[0], res[1] );*/
  free(r_int);
  return res;
}
