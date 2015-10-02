#include "core.h"
/*#include <gsl/gsl_integration.h>*/
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
  double *piL;
  double *piT;
  for(int i=0;i<N;i++) {
    o = 1.5*( (double) i )/( (double) N );
    piL=Pi_htl(o/q,L);
    piT=Pi_htl(o/q,T);

    fprintf(f,
          "%.5f, %.5f, %.5f, %.5f, %.5f\n",
          creal( o ),
          piL[0], piL[1], piT[0], piT[1] 
      );
    free(piL);free(piT);
  }
  fclose(f);

  f = fopen("out/PI, omega.csv","w+");
  fprintf(f,"# HTL photon self-energy\n");
  fprintf(f,"# o/T, Re(Pi_L),  Im(Pi_L),  Re(Pi_T),  Im(Pi_T)\n");

  for(int i=0;i<N;i++) {
    o = 1.5*( (double) i )/( (double) N );
    piL=PI_qed(o,q,L);
    piT=PI_qed(o,q,T);

    fprintf(f,
          "%.5f, %.5f, %.5f, %.5f, %.5f\n",
          creal( o ),
          piL[0], piL[1], piT[0], piT[1] 
      );
    free(piL);free(piT);
  }
  fclose(f);

  return 0;
}
