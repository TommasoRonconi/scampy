#ifndef __TINKER10_P__
#define __TINKER10_P__

// internal includes
#include <occupation_p.h>

namespace sico {

  struct tinker10_p : public occupation_p {
    
    double A_min;
    double sigma_logA;
    double A_sat;
    double alpha_sat;

    tinker10_p () = default;

    tinker10_p ( const double & Amin,
		 const double & siglogA,
		 const double & Asat,
		 const double & alpsat ) :
      A_min { Amin }, sigma_logA { siglogA }, A_sat { Asat }, alpha_sat { alpsat } {}

    tinker10_p ( const tinker10_p & handler ) = default;

    ~tinker10_p () = default;

    void set_parameters ( const double A_min,
			  const double sigma_logA,
			  const double A_sat,
			  const double alpha_sat );

    double Ncen ( const double Mhalo ) override;
    
    double Nsat ( const double Mhalo ) override;
    
  }; //endstruct tinker10_p

} // endnamespace sico

#endif //__TINKER10_P__
