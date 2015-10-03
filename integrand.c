#include "core.h"
#include <gsl/gsl_integration.h>
size_t calls=1e4; double tol=1e-2;

double fBose( float e ) {
  return 1. / (exp( e ) - 1.);
};

double fFermi( float e ) {
  return 1. / (exp( e ) + 1.);
};

double *help_qed(double q, double o, double k, double r, pol X) {
  double r2=r*r, r3=r*r2, r4=r*r3, q2=q*q, o2=o*o, k2=k*k, ko2=(2.*k+o)*(2.*k+o);
  double complex ll = clog(k+r+o*(1+I*1e-3)), D;

  switch(X) {
    case L:
      D =         -.5*r2 + r*(3.*k+o) - (ko2-q2)*(ll) ;          break;
    case T:
      D =         + r4/4. - r3*(k+o)/3. + r2*(o2-k2+2.*k*o)/2.
                  + r*(k*k2 + 4.*k*q2 - k2*o - 3.*k*o2 - o*o2)
                  - (q2-o2)*(ko2+q2)*(ll) ;                      break;
  }

  double *res     = (double *)malloc( 2*sizeof(double) );
          res[0]  = creal(D);
          res[1]  = cimag(D);

  return res;
}

double *Igd_PI_qed(double xi, void *params) {         // the integrand:

  /* members */
  struct Qpol * Q = (struct Qpol *)params;

  double o = Q->o;
  double q = Q->q;
  pol    X = Q->X;

  /*printf("%d\n",X);*/

  /* variable changes */
  double
    k   = xi/(1.-xi),
    /*k   = xi,*/
    k2  = k*k,
    q2  = q*q,
    fk  = fFermi(k);

  double rU = q+k, rL = fabs(q-k);
  double complex res;

  double *e_int  = (double*)malloc(2*sizeof(double));

  double sr, so;
  for (int j=0; j<4; j++) {     sr = (double) (2*(j%2)-1);  so = (double) (2*(j/2)-1);
    e_int = help_qed(q,so*o,k,sr*rU,X);
    res +=  1*(
            e_int[0] +   //    \__,{ upper
          I*e_int[1] ) ;   //    /
    e_int = help_qed(q,so*o,k,sr*rL,X);
    res -=  1*(
            e_int[0] +   //    \__
          I*e_int[1] ) ;   //    /  `{ lower
  }

  free(e_int);

  res *=   (-fk/q )
          *( 3./(4.*M_PI*M_PI) )
          *( 1./( (1.-xi)*(1.-xi) ) ) 
          ; 

  if (X==T) res /= -q2;
  /*
  switch(X) {
    case L:
    res *=   (-fk/q )
            *( 3./(4.*M_PI*M_PI) )
            *( 1./( (1.-xi)*(1.-xi) ) ) ; break;
    case T:
    res *=   (+fk/q )
            *( 3./(4.*M_PI*M_PI) )
            *( 1./( (1.-xi)*(1.-xi) ) ) ; break;
  }*/

  double *Pi  = (double*)malloc(2*sizeof(double));

  Pi[0] = creal( res );
  Pi[1] = cimag( res );

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  /*printf("%.3f + i %.3f \n", creal(1./(o+k+r)), cimag(1./(o+k-r) - 1./(o-k+r)) );*/

  return Pi;
};

double re_PI(double xi, void *params) { 
  double *res = Igd_PI_qed(xi,params); 
  double a = res[0];
  free(res);
  return a;
};

double im_PI(double xi, void *params) { 
  double *res = Igd_PI_qed(xi,params); 
  double a = res[1];
  free(res);
  return a;
};

/*double im_PI(double xi, void *params) { return Igd_PI_qed(xi,params)[1]; }*/
/*gsl_integration_workspace *WS;*/

double *PI_qed(double o, double q, pol X) {

  double                        *Pi  = (double*)malloc(2*sizeof(double));
  gsl_integration_workspace     *WS1 = gsl_integration_workspace_alloc(calls);
  gsl_integration_workspace     *WS2 = gsl_integration_workspace_alloc(calls);

  double res, err; struct Qpol Q = {o,q,X};

  /*printf("%.5f  \n", Q.o);*/

  gsl_function  re_aux; re_aux.function=&re_PI; re_aux.params=&Q;       // nicer to package, re_aux={&...,&...}???
  gsl_function  im_aux; im_aux.function=&im_PI; im_aux.params=&Q;

  double omq = fabs(o-q), opq = fabs(o+q);
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

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  return Pi;
}

