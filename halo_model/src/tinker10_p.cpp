#include <tinker10_p.h>

//==============================================================================================

double scam::tinker10_p::Ncen ( const double Mhalo ) {
  
  double Nc = PP( Mhalo, A_min, sigma_logA );
  
  return ( Nc < 0 || std::isnan( Nc ) ) ? 0. : Nc;
  
}

//==============================================================================================

double scam::tinker10_p::Nsat ( const double Mhalo ) {

  const double Pc = PP( Mhalo, 2. * A_min, sigma_logA );
  const double Ns = Pc * std::pow( ( Mhalo / A_sat ), alpha_sat );
  
  return ( Ns < 0 || std::isnan( Ns ) ) ? 0. : Ns;
  
}

//==============================================================================================
