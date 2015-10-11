/*
 * Author: greg jackson
 * Date: Oct 08 2015
 * thermal loops
 *
 */

#include "core.h"

/*-----------------------------------------------------------------------------------------------*/

double    tol = 1e-1;
size_t  calls = 1e7 ;
int       HTL = 1   ;  // =1 for HTL, =2 for QED, =3 for QCD
double      g = 1. ;

void   eval_PI(double,double); void   eval_disp(double,double); int   points;  // See after main() ...

/*-----------------------------------------------------------------------------------------------*/

int main() {

  points = 20;

  HTL = 1;    eval_PI(0.01, 2.);    eval_disp(0.001, 1.);
  HTL = 2;    eval_PI(0.01, 2.);    eval_disp(0.001, 1.);
  HTL = 3;    eval_PI(0.01, 2.5);    eval_disp(0.05, 1.);
  /*eval_disp(0.001, 2.);*/

  return 0;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

void eval_PI(double o_min, double o_max)
{

       if (HTL==1) sprintf(fname, "out/Pi, htl.csv"                                               );
  else if (HTL==2) sprintf(fname, "out/Pi, qed.csv"                                               );
  else if (HTL==3) sprintf(fname, "out/Pi, qcd.csv"                                               );

  file = fopen(fname,"w+");

       if (HTL==1) fprintf(file, "# HTL approx\n"                                                 );
  else if (HTL==2) fprintf(file, "# QED photon\n"                                                 );
  else if (HTL==3) fprintf(file, "# QCD gluon\n"                                                  );

  fprintf(file,   "#\n"                                                                           );
  fprintf(file,   "# o/T, Re(Pi_L),  Im(Pi_L),  Re(Pi_T),  Im(Pi_T)\n"                            );

  double o;
  double q=1.; 
  double *piL;   double *piT;

  for(int i=0; i<points; i++) {
    o = o_min + (o_max-o_min)*( (double) i )/( (double) points );

         if (HTL==1) {    piL=Pi_htl(o/q,L),    piT=Pi_htl(o/q,T);    }
    else if (HTL==2) {    piL=Pi_qed(o,q,L);    piT=Pi_qed(o,q,T);    }
    else if (HTL==3) {    piL=Pi_qcd(o,q,L);    piT=Pi_qcd(o,q,T);    }

    fprintf(file,
          "%.5f, %.5f, %.5f, %.5f, %.5f\n",
          creal( o ),
          piL[0], piL[1], piT[0], piT[1] 
      );
    free(piL);free(piT);
  }
  fclose(file);
}


void eval_disp(double q_min, double q_max)
{

       if (HTL==1) sprintf(fname, "out/omega(q), htl.csv"                                         );
  else if (HTL==2) sprintf(fname, "out/omega(q), qed.csv"                                         );
  else if (HTL==3) sprintf(fname, "out/omega(q), qcd.csv"                                         );

  file = fopen(fname,"w+");

       if (HTL==1) fprintf(file, "# HTL approx\n"                                                 );
  else if (HTL==2) fprintf(file, "# QED photon\n"                                                 );
  else if (HTL==3) fprintf(file, "# QCD gluon\n"                                                  );

  fprintf(file,   "# pole in propagator\n"                                                        );
  fprintf(file,   "#\n"                                                                           );
  /*fprintf(file,   "# o/m, q^2_L(o)/m, q^2_T(o)/m\n"                                           );*/
  fprintf(file,   "# q/m, omega_L(q)/m, omega_T(q)/m\n"                                           );

  double q; 

  for(int i=0; i<points; i++) {
    q = q_min + (q_max-q_min)*( (double) i )/( (double) points );

         /*if (HTL==1) {    piL=Pi_htl(o/q,L),    piT=Pi_htl(o/q,T);    }*/
    /*else if (HTL==2) {    piL=Pi_qed(o,q,L);    piT=Pi_qed(o,q,T);    }*/
    /*else if (HTL==3) {    piL=Pi_qcd(o,q,L);    piT=Pi_qcd(o,q,T);    }*/

    fprintf(file,
          "%.5f, %.5f, %.5f\n",
          q, disp(q,L), disp(q,T)
      );
  }
  fclose(file);
}

