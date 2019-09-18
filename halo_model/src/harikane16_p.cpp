#include <harikane16_p.h>

//==============================================================================================

double sico::harikane16_p::Ncen ( const double Mhalo ) {
  
  double Nc = PP( Mhalo, M_min, sigma_logM );
  
  return ( Nc < 0 || std::isnan( Nc ) ) ? 0. : Nc;
  
}

//==============================================================================================

double sico::harikane16_p::Nsat ( const double Mhalo ) {

  const double Nc = sico::harikane16_p::Ncen( Mhalo );
  const double Ns = Nc * std::pow( ( ( Mhalo - M0 ) / M1 ), alpha );
  
  return ( Ns < 0 || std::isnan( Ns ) ) ? 0. : Ns;
  
}

//==============================================================================================
