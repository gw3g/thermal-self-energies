#include "core.h"

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
  return 0;
}
