#include "core.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

const gsl_root_fsolver_type *rst;
gsl_root_fsolver *rs;
int   HTL;

double D_inv(double q2, void *params) {
  struct Qpol * Q = (struct Qpol *)params;

  double o2  = Q->o;
  pol X     = Q->X;

  double complex o = csqrt(o2);
  double complex q = csqrt(fabs(q2));
  q*=csqrt(q2/fabs(q2));

  o*=g;q*=g;

  double g2 = g*g;   double re_Pi;

       if (HTL==1) { re_Pi = Pi_htl(o/q,X)[0]; }
  else if (HTL==2) { re_Pi = Pi_qed(o,q,X)[0]; }
  else if (HTL==3) { re_Pi = Pi_qcd(o,q,X)[0]; }

    /*printf(" (%.5f, %.5f) : rePi = %.9f \n", o2, q2, re_Pi);*/

  switch (X) {  case L:  return (  (g2*(q2     ) - re_Pi)  );
                case T:  return (  (g2*(o2 - q2) - re_Pi)  );   }
}


double disp(double o2, pol X) {

  int iter = 0, max_iter = 1000;

  double r = 0;

  double q2_lo, q2_hi;
  
  // root finding withing [a,b] ... need to be careful to avoid spurious regions...
  if ( (o2<1.)&&(HTL==1) ) {
  switch (X) {  case L:  q2_hi = -1e-3; q2_lo =  q2_hi - 3.2; 
                case T:  q2_hi = -1e-3; q2_lo =  q2_hi - 3.2;   }
  }
  /*else if (o2<1.1) {*/
  /*switch (X) {  case L:  q2_hi = 1.3; q2_lo = -.5; */
                /*case T:  q2_hi = 1.3; q2_lo = -.5;  return 0;}*/
  /*}*/
  else {
  switch (X) {  case L:  q2_hi = o2-1e-8; q2_lo =  1e-8; 
                case T:  q2_hi = o2-1e-8; q2_lo =  1e-8;   }
  }

  /*q2_lo = ( (o2>0) ? q2 : 0.) + 1e-6; q2_hi = q2_lo+1.;*/

  gsl_function F;
  /*struct quadratic_params params = {1.0, 0.0, -5.0};*/
  struct Qpol Q = {o2, 0., X};

  F.function = &D_inv;
  F.params = &Q;

  rst = gsl_root_fsolver_brent;
  rs = gsl_root_fsolver_alloc (rst);
  gsl_root_fsolver_set (rs, &F, q2_lo, q2_hi);

  do    {
            iter++;
            gsl_root_fsolver_iterate (rs);
            r = gsl_root_fsolver_root (rs);
            q2_lo = gsl_root_fsolver_x_lower (rs);
            q2_hi = gsl_root_fsolver_x_upper (rs);
            /*printf("%d ---q: %.4f, lo: %.4f, hi: %.4f\n", HTL, o2, q2_lo, q2_hi );*/
            /*if ((fabs(q2_hi)<1e-9)&&( fabs(q2_lo)<1e-9)) { q2_lo = 1e-4; q2_hi = o2 + 2.; iter = 0;}*/
        }
  while ( iter < max_iter);

  gsl_root_fsolver_free (rs);

  return r;
}
