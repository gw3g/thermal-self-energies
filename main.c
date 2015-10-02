#include "core.h"
/*#include <gsl/gsl_integration.h>*/

/*size_t calls=1e3; double tol=1e-4;*/

/*double reL(double o, double q) {*/

  /*gsl_integration_workspace * WS = gsl_integration_workspace_alloc(calls);*/

  /*double res, err; struct pair Q = {o,q};*/

  /*gsl_function  f_aux                     ;*/
                /*f_aux.function  = &re_PiL ;*/
                /*f_aux.params    = &Q      ;*/

  /*gsl_integration_qags (&f_aux, 0, 1, 0, tol, calls, WS, &res, &err);*/
  /*gsl_integration_workspace_free (WS);*/

  /*return res;*/
/*}*/

/*double imL(double o, double q) {*/

  /*gsl_integration_workspace * WS = gsl_integration_workspace_alloc(calls);*/

  /*double res, err; struct pair Q = {o,q};*/

  /*gsl_function  f_aux                     ;*/
                /*f_aux.function  = &im_PiL ;*/
                /*f_aux.params    = &Q      ;*/

  /*gsl_integration_qags (&f_aux, 0, 1, 0, tol, calls, WS, &res, &err);*/
  /*gsl_integration_workspace_free (WS);*/

  /*return res;*/
/*}*/

/*double reT(double o, double q) {*/

  /*gsl_integration_workspace * WS = gsl_integration_workspace_alloc(calls);*/

  /*double res, err; struct pair Q = {o,q};*/

  /*gsl_function  f_aux                     ;*/
                /*f_aux.function  = &re_PiT ;*/
                /*f_aux.params    = &Q      ;*/

  /*gsl_integration_qags (&f_aux, 0, 1, 0, tol, calls, WS, &res, &err);*/
  /*gsl_integration_workspace_free (WS);*/

  /*return res;*/
/*}*/

/*double imT(double o, double q) {*/

  /*gsl_integration_workspace * WS = gsl_integration_workspace_alloc(calls);*/

  /*double res, err; struct pair Q = {o,q};*/

  /*gsl_function  f_aux                     ;*/
                /*f_aux.function  = &im_PiT ;*/
                /*f_aux.params    = &Q      ;*/

  /*gsl_integration_qags (&f_aux, 0, 1, 0, tol, calls, WS, &res, &err);*/
  /*gsl_integration_workspace_free (WS);*/

  /*return res;*/
/*}*/

int main() {

  FILE *f;
  f = fopen("out/PI_htl, omega.csv","w+");
  fprintf(f,"# HTL photon self-energy\n");
  fprintf(f,"# o/T, Re(Pi_L),  Im(Pi_L),  Re(Pi_T),  Im(Pi_T)\n");

  int N = 100;
  double o;
  double q=1.;
  for(int i=0;i<N;i++) {
    o = 1.5*( (double) i )/( (double) N );

    fprintf(f,
          "%.5f, %.5f, %.5f, %.5f, %.5f\n",
          creal( o ),
          Pi_htl(o/q,L)[0], Pi_htl(o/q,L)[1], Pi_htl(o/q,T)[0], Pi_htl(o/q,T)[1] 
      );
  }
  fclose(f);

  f = fopen("out/TEST__PI, omega.csv","w+");
  fprintf(f,"# HTL photon self-energy\n");
  fprintf(f,"# o/T, Re(Pi_L),  Im(Pi_L),  Re(Pi_T),  Im(Pi_T)\n");

  for(int i=0;i<N;i++) {
    o = 1.5*( (double) i )/( (double) N );

    fprintf(f,
          "%.5f, %.5f, %.5f, %.5f, %.5f\n",
          creal( o ),
          1.5*PI_qed(o,q,L)[0],
          1.5*PI_qed(o,q,L)[1],
          1.5*PI_qed(o,q,T)[0],
          1.5*PI_qed(o,q,T)[1]
      );
  }
  fclose(f);

  return 0;
}
