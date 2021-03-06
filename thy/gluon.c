#include "core.h"
#include <gsl/gsl_integration.h>
size_t calls;
double tol;
double g;

double *Igd_PI_qcd(double xi, void *params) {         // the integrand:

  /* members */
  struct Qpol * Q = (struct Qpol *)params;

  double complex o = Q->o;
  double complex q = Q->q;
  pol    X = Q->X;

  /* variable changes */
  double
    k   = xi/(1.-xi),
    q2  = q*q,
    o2  = o*o,
    fk  = f(k,B);

  double complex res = 0.;

  double *e_int  = (double*)malloc(2*sizeof(double));

  switch(X) {
    case L:
              e_int = frakJ(k,Q,0); res += e_int[0] + I*e_int[1];
              e_int = frakJ(k,Q,1); res += e_int[0] + I*e_int[1];
              e_int = frakJ(k,Q,3); res += .5*q2*( e_int[0] + I*e_int[1] );         break;
    case T:
              e_int = frakJ(k,Q,0); res += e_int[0] + I*e_int[1];
              e_int = frakJ(k,Q,2); res -= e_int[0] + I*e_int[1];
              e_int = frakJ(k,Q,3); res += .5*(o2-q2)*( e_int[0] + I*e_int[1] );     break;
  }

  free(e_int);

  res *=  -fk                                         // thermal weights
          *6.*3.                                      // ??
          *( 1./(4.*M_PI*M_PI) )                      // angular prefactors
          *( 1./( (1.-xi)*(1.-xi) ) )                 // jacobian
          *g*g
          ; 

  double *Pi  = (double*)malloc(2*sizeof(double));

  Pi[0] = creal( res );
  Pi[1] = cimag( res );

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  /*printf("o:%.3f q:%.3f \n",o,q);*/
  /*printf("%.3f + i %.3f \n", creal(1./(o+k+r)), cimag(1./(o+k-r) - 1./(o-k+r)) );*/

  return Pi;
};


double *Pi_qcd(double complex o, double complex q, pol X) {

  double                        *Pi  = (double*)malloc(2*sizeof(double));
  gsl_integration_workspace     *WS  = gsl_integration_workspace_alloc(calls);

  struct Qpol                     Q  = {o,q,X};
  double complex                 q2  = q*q;
  double                                                                            res, err, 
                                omq  = fabs(o-q)/2., 
                                opq  = fabs(o+q)/2.;

  double                      pts[4] = {0., omq/(omq+1.), opq/(opq+1.), 1.};

  for (int i=0;i<2;i++) {

    double          PI(double xi, void *params) { return Igd_PI_qcd(xi,params)[i]; };

    gsl_function    aux         = { &PI, &Q };

    gsl_integration_qagp (&aux,  pts, 4, tol, 0, calls, WS, &res, &err);        Pi[i] = res/q2;

  }

  gsl_integration_workspace_free (WS);

  /*printf("%.3f  ::  %.3f + i %.3f \n", o, Pi[0], Pi[1]);*/
  return Pi;
}


