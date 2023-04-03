/**
 *  @file cosmology/include/power_spectrum.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __POWER_SPECTRUM__
#define __POWER_SPECTRUM__

#include <map>
#include <string>

#include <utilities.h>
#include <interpolation.h>

/// General namespace of the library 
namespace scam {

  struct power_spec {

    using interp_lin = class scam::utl::interpolator< scam::utl::gsl_lin_interp >;
    interp_lin P0 {};
    double ss8_P0 = -1;
    std::size_t thinness = 200;

    power_spec ( const std::vector< double > & kh0,
		 const std::vector< double > & Pk0,
		 const std::size_t thinness = 200 ) :
      P0 { interp_lin { kh0, Pk0 } }, thinness { thinness } {

      compute_ss8_P0();

    }

    power_spec () = default;
    ~power_spec () = default;

    void compute_ss8_P0 ();
    
    double Pk_comoving ( const double & kk, const double & zz = 1.e-7 );
    double Pk ( const double & kk, const double & zz = 1.e-7 );

    // radius in Mpc/h 
    double sigma2R_comoving ( const double & Radius, const double & zz = 1.e-7 );
    double sigma2R ( const double & Radius, const double & zz = 1.e-7 );

    // radius in Mpc/h 
    double dsigma2RdR_comoving ( const double & Radius, const double & zz = 1.e-7 );
    double dsigma2RdR ( const double & Radius, const double & zz = 1.e-7 );

    // mass in solar masses ( not M_sol/h )
    double sigma2M ( const double & mm, const double & zz = 1.e-7 );

    // mass in solar masses
    double dsigma2MdM ( const double & mm, const double & zz = 1.e-7 );

    /// @} End of power spectrum dependent functions

  }; // endclass power_spec

} //endnamespace scam

#endif //__POWER_SPECTRUM__
