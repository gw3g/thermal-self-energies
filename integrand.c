#include "core.h"

double fBose( float e ) {
  return 1. / (exp( e ) - 1.);
};

double fFermi( float e ) {
  return 1. / (exp( e ) + 1.);
};

double *i1(double q, double o, double k, double r) {
  double r2=r*r, q2=q*q, ko2=(2.*k+o)*(2.*k+o);

  double complex 
    D = -.5*r2 + r*(3.*k+o) - (ko2-q2)*(clog(k+r+o)) ;

  double *res     = (double *)malloc( 2*sizeof(double) );
          res[0]  = creal(D);
          res[1]  = cimag(D);

  return res;
}


double *PiL(double xi, void *params) {         // the integrand:

  /* members */
  struct pair * Q = (struct pair *)params;

  double complex o = Q->o;
  double complex q = Q->q;

  /* variable changes */
  double
    k   = xi/(1.-xi),
    ko2=(2.*k+o)*(2.*k+o),
    k2  = k*k,
    q2  = q*q,
    fk  = fFermi(k);

  double rU = q+k, rL = q-k;

  double complex res =

    (-( i1(q,o,k,rU)[0] + i1(q,o,k,-rU)[0] + i1(q,-o,k,rU)[0] + i1(q,-o,k,-rU)[0] ) +
      ( i1(q,o,k,rL)[0] + i1(q,o,k,-rL)[0] + i1(q,-o,k,rL)[0] + i1(q,-o,k,-rL)[0] ) ) + 

  I*(-( i1(q,o,k,rU)[1] + i1(q,o,k,-rU)[1] + i1(q,-o,k,rU)[1] + i1(q,-o,k,-rU)[1] ) +
      ( i1(q,o,k,rL)[1] + i1(q,o,k,-rL)[1] + i1(q,-o,k,rL)[1] + i1(q,-o,k,-rL)[1] ) )   ;

  res *=   ( 2.*fk/q )
          *( 1./(4.*M_PI*M_PI) )
          *( 1./( (1.-xi)*(1.-xi) ) ) ;

  double *Pi  = (double*)malloc(2*sizeof(double));

  Pi[0] = creal( res );
  Pi[1] = cimag( res );

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  /*printf("%.3f + i %.3f \n", creal(1./(o+k+r)), cimag(1./(o+k-r) - 1./(o-k+r)) );*/

  return Pi;
};

double re_PiL(double xi, void *params) {
  double *res = PiL(xi,params);
  double a = res[0];
  free(res);
  return a;
};
double im_PiL(double xi, void *params) {
  double *res = PiL(xi,params);
  double a = res[1];
  free(res);
  return a;
};

