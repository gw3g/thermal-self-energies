#include "core.h"
/* distribution functions */

double fBose( float e ) {
  return 1. / (exp( e ) - 1.);
};

double fFermi( float e ) {
  return 1. / (exp( e ) + 1.);
};

double f(float e, p_type X) {
  switch(X) {
    case F: return fFermi(e);
    case B: return fBose(e);
  }
};

double bf(float e, p_type X) {
  switch(X) {
    case F: return 1.-fFermi(e);
    case B: return 1.+fBose(e);
  }
};
