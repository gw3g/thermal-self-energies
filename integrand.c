#include "core.h"
#include <gsl/gsl_integration.h>

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
    case L: { D = -.5*r2 + r*(3.*k+o) - (ko2-q2)*(ll) ;          break; }
    case T: { D =   r4/4. - r3*(k+o)/3. + r2*(o2-k2+2.*k*o)/2.
                  + r*(k*k2 + 4.*k*q2 - k2*o - 3.*k*o2 - o*o2)
                  - (q2-o2)*(ko2+q2)*(ll) ;                      break; }
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

  /* variable changes */
  double
    k   = xi/(1.-xi),
    k2  = k*k,
    q2  = q*q,
    fk  = fFermi(k);

  double rU = q+k, rL = fabs(q-k);
  double complex res;

  double sr, so;
  for (int j=0; j<8; j++) {     sr = (double) (2*(j%2)-1);  so = (double) (2*(j/2)-1);
    res +=
            help_qed(q,so*o,k,sr*rU,X)[0] +   //    \__,{ upper
          I*help_qed(q,so*o,k,sr*rU,X)[1] -   //    /
            help_qed(q,so*o,k,sr*rL,X)[0] -   //    \__
          I*help_qed(q,so*o,k,sr*rL,X)[1] ;   //    /  `{ lower
  }

  switch(X) {
    case L:
    res *=   (-2.*fk/q )
            *( 1./(4.*M_PI*M_PI) )
            *( 1./( (1.-xi)*(1.-xi) ) ) ; break;
    case T:
    res *=   (+fk/pow(q,3) )
            *( 1./(2.*M_PI*M_PI) )
            *( 1./( (1.-xi)*(1.-xi) ) ) ; break;
  }

  double *Pi  = (double*)malloc(2*sizeof(double));

  Pi[0] = creal( res );
  Pi[1] = cimag( res );

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  /*printf("%.3f + i %.3f \n", creal(1./(o+k+r)), cimag(1./(o+k-r) - 1./(o-k+r)) );*/

  return Pi;
};

/*double re_PiL(double xi, void *params) {*/
  /*double *res = PiL(xi,params);*/
  /*double a = res[0];*/
  /*free(res);*/
  /*return a;*/
/*};*/
/*double im_PiL(double xi, void *params) {*/
  /*double *res = PiL(xi,params);*/
  /*double a = res[1];*/
  /*free(res);*/
  /*return a;*/
/*};*/


/*double *PiT(double xi, void *params) {         // the integrand:*/

  /*[> members <]*/
  /*struct pair * Q = (struct pair *)params;*/

  /*double complex o = Q->o;*/
  /*double complex q = Q->q;*/

  /*[> variable changes <]*/
  /*double*/
    /*k   = xi/(1.-xi),*/
    /*ko2=(2.*k+o)*(2.*k+o),*/
    /*k2  = k*k,*/
    /*q2  = q*q,*/
    /*fk  = fFermi(k);*/

  /*double rU = q+k, rL = fabs(q-k);*/

  /*double complex res =*/

    /*( ( i2(q,o,k,rU)[0] + i2(q,o,k,-rU)[0] + i2(q,-o,k,rU)[0] + i2(q,-o,k,-rU)[0] ) -*/
      /*( i2(q,o,k,rL)[0] + i2(q,o,k,-rL)[0] + i2(q,-o,k,rL)[0] + i2(q,-o,k,-rL)[0] ) ) + */

  /*I*(+( i2(q,o,k,rU)[1] + i2(q,o,k,-rU)[1] + i2(q,-o,k,rU)[1] + i2(q,-o,k,-rU)[1] ) -*/
      /*( i2(q,o,k,rL)[1] + i2(q,o,k,-rL)[1] + i2(q,-o,k,rL)[1] + i2(q,-o,k,-rL)[1] ) )   ;*/

  /*res *=   (+fk/pow(q,3) )*/
          /**( 1./(2.*M_PI*M_PI) )*/
          /**( 1./( (1.-xi)*(1.-xi) ) ) ;*/

  /*double *Pi  = (double*)malloc(2*sizeof(double));*/

  /*Pi[0] = creal( res );*/
  /*Pi[1] = cimag( res );*/

  /*[>printf("%.3f + i %.3f \n", Pi[0], Pi[1]);<]*/
  /*[>printf("%.3f + i %.3f \n", creal(1./(o+k+r)), cimag(1./(o+k-r) - 1./(o-k+r)) );<]*/

  /*return Pi;*/
/*};*/
size_t calls=1e3; double tol=1e-4;

double *PI_qed(double o, double q, pol X) {

  double                        *Pi  = (double*)malloc(2*sizeof(double));
  gsl_integration_workspace     *WS1 = gsl_integration_workspace_alloc(calls);
  gsl_integration_workspace     *WS2 = gsl_integration_workspace_alloc(calls);

  double res, err; struct Qpol Q = {o,q,X};

  double re_PI(double xi, void *params) { return Igd_PI_qed(xi,params)[0]; }
  double im_PI(double xi, void *params) { return Igd_PI_qed(xi,params)[1]; }

  gsl_function  re_aux; re_aux.function=&re_PI; re_aux.params=&Q;       // nicer to package, re_aux={&...,&...}???
  gsl_function  im_aux; im_aux.function=&im_PI; im_aux.params=&Q;

  gsl_integration_qags (&re_aux, 0, 1, 0, tol, calls, WS1, &res, &err);        Pi[0] = res;
  gsl_integration_qags (&im_aux, 0, 1, 0, tol, calls, WS2, &res, &err);        Pi[1] = res;
  gsl_integration_workspace_free (WS1);gsl_integration_workspace_free (WS2);

  /*printf("%.3f + i %.3f \n", Pi[0], Pi[1]);*/
  return Pi;
}

/*double re_PiT(double xi, void *params) {*/
  /*double *res = PiT(xi,params);*/
  /*double a = res[0];*/
  /*free(res);*/
  /*return a;*/
/*};*/
/*double im_PiT(double xi, void *params) {*/
  /*double *res = PiT(xi,params);*/
  /*double a = res[1];*/
  /*free(res);*/
  /*return a;*/
/*};*/

