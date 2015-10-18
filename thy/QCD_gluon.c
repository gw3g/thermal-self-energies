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
          *6.*3.                                     // ??
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

double *Pi_qcd(double o, double q, pol X) {

  double                        *Pi  = (double*)malloc(2*sizeof(double));
  double q2 = q*q;

  /*if ( (fabs(q)<1e-3)||(fabs(o)<1e-3) ) { // generous ``safety net''*/
    /*Pi[0] = Pi_htl(o/q,X)[0]; Pi[1] = Pi_htl(o/q,X)[1]; return Pi;*/
  /*}*/
  gsl_integration_workspace     *WS1 = gsl_integration_workspace_alloc(calls);
  gsl_integration_workspace     *WS2 = gsl_integration_workspace_alloc(calls);

  double res, err; struct Qpol Q = {o,q,X};

  /*printf("%.5f  \n", Q.o);*/

  double re_PI(double xi, void *params) { return Igd_PI_qcd(xi,params)[0]; };
  double im_PI(double xi, void *params) { return Igd_PI_qcd(xi,params)[1]; };

  gsl_function  re_aux; re_aux.function=&re_PI; re_aux.params=&Q;       // nicer to package, re_aux={&...,&...}???
  gsl_function  im_aux; im_aux.function=&im_PI; im_aux.params=&Q;

  double omq = fabs(o-q)/2., opq = fabs(o+q)/2.;
  /*printf("%.3f  %.3f \n", omq, opq);*/
  /*printf("%.3f  %.3f \n", omq/(omq+1.), opq/(opq+1.));*/
  double pts[4] = {0., omq/(omq+1.), opq/(opq+1.), 1.};
  /*double                        *pts  = (double*)malloc(4*sizeof(double));*/
  /*pts[0] = 0.; pts[1] = omq/(omq+1.); pts[2] = opq/(opq+1.); pts[3] = 1.;*/
  gsl_integration_qagp (&re_aux,  pts, 4, tol, 0, calls, WS1, &res, &err);        Pi[0] = res/q2;
  gsl_integration_qagp (&im_aux,  pts, 4, tol, 0, calls, WS1, &res, &err);        Pi[1] = res/q2;
  /*WS = gsl_integration_workspace_alloc(calls);*/
  /*gsl_integration_qags (&re_aux,  0, 1, 0, tol, calls, WS1, &res, &err);        Pi[0] = res;*/
  /*gsl_integration_qags (&im_aux,  0, 1, 0, tol, calls, WS2, &res, &err);        Pi[1] = res;*/
  gsl_integration_workspace_free (WS1);
  gsl_integration_workspace_free (WS2);
  /*free(pts);*/

  /*printf("%.3f  ::  %.3f + i %.3f \n", o, Pi[0], Pi[1]);*/
  return Pi;
}
