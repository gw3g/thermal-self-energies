#include "core.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
size_t calls;
double tol;

double *Igd_PI_qcd(double xi, void *params) {         // the integrand:

  /* members */
  struct Qpol * Q = (struct Qpol *)params;

  double o = Q->o;
  double q = Q->q;
  pol    X = Q->X;

  /*printf("%d\n",X);*/

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

  res *=  fk //(-fk/q )
          *( 12./(8.*M_PI*M_PI) )
          *( 1./( (1.-xi)*(1.-xi) ) ) 
          ; 


  double *Pi  = (double*)malloc(2*sizeof(double));

  Pi[0] = creal( res );
  Pi[1] = cimag( res );

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  /*printf("%.3f + i %.3f \n", creal(1./(o+k+r)), cimag(1./(o+k-r) - 1./(o-k+r)) );*/

  return Pi;
};

double *PI_qcd(double o, double q, pol X) {


  double                        *Pi  = (double*)malloc(2*sizeof(double));
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
  gsl_integration_qagp (&re_aux,  pts, 4, 0, tol, calls, WS1, &res, &err);        Pi[0] = res;
  gsl_integration_qagp (&im_aux,  pts, 4, 0, tol, calls, WS1, &res, &err);        Pi[1] = res;
  /*WS = gsl_integration_workspace_alloc(calls);*/
  /*gsl_integration_qags (&re_aux,  0, 1, 0, tol, calls, WS1, &res, &err);        Pi[0] = res;*/
  /*gsl_integration_qags (&im_aux,  0, 1, 0, tol, calls, WS2, &res, &err);        Pi[1] = res;*/
  gsl_integration_workspace_free (WS1);
  gsl_integration_workspace_free (WS2);
  /*free(pts);*/

  /*printf("%.3f  ::  %.3f + i %.3f \n", o, Pi[0], Pi[1]);*/
  return Pi;
}

double D_inv(double o, void *params) {
  struct Qpol * Q = (struct Qpol *)params;
  double q = Q->q;
  pol X = Q->X;

  switch (X) {
    case L:  return (  q*q + Pi_htl(o/q,X)[0]  );
    case T:  return (  o*o - q*q - Pi_htl(o/q,X)[0]  );
  }
}


double omega_g(double q, pol X) {

  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0, r_expected = sqrt (5.0);
  double o_lo = q+1e-6, o_hi = q + 1.;
  gsl_function F;
  /*struct quadratic_params params = {1.0, 0.0, -5.0};*/
  struct Qpol Q = {0., q, X};

  F.function = &D_inv;
  F.params = &Q;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, o_lo, o_hi);

  printf ("using %s method\n", 
          gsl_root_fsolver_name (s));

  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "root", 
          "err", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      o_lo = gsl_root_fsolver_x_lower (s);
      o_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (o_lo, o_hi,
                                       0, 0.001);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
              iter, o_lo, o_hi,
              r, r - r_expected, 
              o_hi - o_lo);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}
