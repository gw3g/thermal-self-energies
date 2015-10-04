#include "core.h"

double tol = 1e-2;
size_t calls = 1e5;

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
    o = 2.*( (double) i )/( (double) N );
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
  fprintf(f,"# HTL gluon self-energy\n");
  fprintf(f,"# o/T, Re(Pi_L),  Im(Pi_L),  Re(Pi_T),  Im(Pi_T)\n");

  for(int i=0;i<N;i++) {
    o = 2.*( (double) i+1 )/( (double) N );
    piL=PI_qcd(o,q,L);
    piT=PI_qcd(o,q,T);
    if (fabs(o-q)<1e-3) continue;

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
