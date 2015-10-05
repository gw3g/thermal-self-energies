#include "core.h"

double *fAUX(double r, double k, double o, double q, int i) {

  double complex ll = clog(k+r+o*(1+I*1e-6)),  D;
  double r2=r*r, r3=r*r2, r4=r*r3, 
         q2=q*q, o2=o*o, 
         k2=k*k, k3=k*k2, ko = k+o;

  double *res     = (double*)malloc(2*sizeof(double));    // 0-REAL, 1-IMAG

  if (i==0) {

       D   = +  k*(    + r - ko*ll                                   );           }

  else if (i==1) {

       D   = +.50*(    -.5*r2 + r*ko - (ko*ko + k2 - q2)*ll          );           }

  else if (i==2) {

       D   = -.25*(    + .25*r4 - r3*ko/3. + .5*r2*(ko*ko - 2.*k2)
                            + r*(k3 - k2*o - 3.*k*o2 - o*o2 ) 
                            + ( o2*(k+ko)*(k+ko) - q2*q2 )*ll     )/q2 ;          }

  else if (i==3) {

       D   = +   ll;                                                              }

  res[0] = creal(D); res[1] = cimag(D);                                  return res;

}

double *frakJ (double k, void *params, int i) {
  /*
   * k-integrand
   * params = {o, q, X};
   *
   * i =    0,    1,    2,    3
   *        J00   Jii  qJij   Jij
   */

  struct Qpol * Q = (struct Qpol *)params;

  double  o = Q->o;
  double  q = Q->q;
  pol     X = Q->X;

  double *r_int   = (double*)malloc(2*sizeof(double));    // 0-REAL, 1-IMAG
  double *res     = (double*)malloc(2*sizeof(double));    //

  double sr, so;
  double rU = fabs(q+k), rL = fabs(q-k);

  res[0]  = 0. ; res[1] = 0. ;

  for (int j=0; j<4; j++) {                                       sr  = (double) (2*(j%2)-1);
                                                                  so  = (double) (2*(j/2)-1);

      r_int = fAUX(   sr*rU,  k,  so*o,   q,  i   );      res[0] +=  r_int[0]  ;//  \__ upper
                                                          res[1] +=  r_int[1]  ;//  /
      r_int = fAUX(   sr*rL,  k,  so*o,   q,  i   );      res[0] -=  r_int[0]  ;//  \__ lower
                                                          res[1] -=  r_int[1]  ;//  /
                                                                                                      };

  double fk = f(k,X);                                     res[0] *= .5*fk/q;
                                                          res[1] *= .5*fk/q;
  free(r_int);
  return res;
}
