#include <occupation_p.h>

//==============================================================================================

double sico::occupation_p::PP ( const double AA, const double Amin, const double sigma_logA ) {

  double xx = ( std::log10( AA ) - std::log10( Amin ) ) / sigma_logA;
  double erf = gsl_sf_erf( xx );

  return 0.5 * ( 1. + erf );

}

//==============================================================================================

