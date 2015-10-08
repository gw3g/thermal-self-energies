#include "core.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

const gsl_root_fsolver_type *rst;
gsl_root_fsolver *rs;
int HTL;

double D_inv(double o, void *params) {
  struct Qpol * Q = (struct Qpol *)params;

  double q  = Q->q;
  pol X     = Q->X;

  double o2 = o*o, q2=q*q;   double re_Pi;

       if (HTL==1) { re_Pi = Pi_htl(o/q,X)[0]; }
  else if (HTL==2) { re_Pi = Pi_qed(o,q,X)[0]; }
  else if (HTL==3) { re_Pi = Pi_qcd(o,q,X)[0]; }

  switch (X) {  case L:  return (  q2      - re_Pi  );
                case T:  return (  o2 - q2 - re_Pi  );   }
}


double disp(double q, pol X) {

  int iter = 0, max_iter = 1000;

  double r = 0;

  double o_lo, o_hi;
  o_lo = q + 1e-8; o_hi = q+0.8;

  gsl_function F;
  /*struct quadratic_params params = {1.0, 0.0, -5.0};*/
  struct Qpol Q = {0., q, X};

  F.function = &D_inv;
  F.params = &Q;

  rst = gsl_root_fsolver_brent;
  rs = gsl_root_fsolver_alloc (rst);
  gsl_root_fsolver_set (rs, &F, o_lo, o_hi);

  do    {
            iter++;
            gsl_root_fsolver_iterate (rs);
            r = gsl_root_fsolver_root (rs);
            o_lo = gsl_root_fsolver_x_lower (rs);
            o_hi = gsl_root_fsolver_x_upper (rs);
        }
  while ( iter < max_iter);

  gsl_root_fsolver_free (rs);

  return r;
}
