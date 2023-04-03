#include <power_spectrum.h>
#include <iostream>

//==============================================================================================

void scam::power_spectrum::compute_ss8_P0 () {

  if ( ss8_P0 < 0. ) {

    auto f_integrand = [ & ] ( double kk ) {

      double WF = scam::utl::TopHat_WF( 8 * kk );
      return kk * kk * P0( kk ) * WF * WF;

    };

    scam::power_spectrum::interp_lin int_f ( f_integrand,
					     1.e-4,
					     100.,
					     thinness );

    ss8_P0 = 0.5 * scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

  }

  return;

}

//==============================================================================================

double scam::power_spectrum::Pk_comoving ( const double & kk, const double & zz ) {

  const double Dz = DD( zz );
  const double D0 = DD( 1.e-7 );
  
  return param->sigma8 * param->sigma8 * Dz * Dz * P0_f( kk ) / ( D0 * D0 * ss8_P0 );

}

double scam::power_spectrum::Pk ( const double & kk, const double & zz ) {

  const double Dz = DD( zz );
  const double D0 = DD( 1.e-7 );
  
  return h_3 * param->sigma8 * param->sigma8 * Dz * Dz * P0_f( kk * h_1 ) / ( D0 * D0 * ss8_P0 );

}

//==============================================================================================

double scam::power_spectrum::sigma2R_comoving ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    return  kk * kk * Pk_comoving( kk, zz ) * WF * WF;

  };

  scam::power_spectrum::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );
  
  return 0.5 * scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}


double scam::power_spectrum::sigma2R ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    return  kk * kk * Pk( kk, zz ) * WF * WF;

  };
  

  scam::power_spectrum::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );
  
  return 0.5 * scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}

//==============================================================================================

double scam::power_spectrum::dsigma2RdR_comoving ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    double dWF = scam::utl::TopHat_WF_D1( rr * kk );
    return kk * kk * kk * Pk_comoving( kk, zz ) * WF * dWF;

  };
  
  scam::power_spectrum::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );

  return scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}

double scam::power_spectrum::dsigma2RdR ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    double dWF = scam::utl::TopHat_WF_D1( rr * kk );
    return  kk * kk * kk * Pk( kk, zz ) * WF * dWF;

  };

  scam::power_spectrum::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );

  return scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}

//==============================================================================================

double scam::power_spectrum::sigma2M ( const double & mm, const double & zz ) {

  const double rho = scam::power_spectrum::OmegaM( 1.e-7 ) * rho_crit( 1.e-7 );
  const double rad = std::pow( 0.75 * scam::ip * mm / rho, 0.33333333333333333 );
  
  return scam::power_spectrum::sigma2R( rad, zz );

}

//==============================================================================================

double scam::power_spectrum::dsigma2MdM ( const double & mm, const double & zz ) {

  const double rho = scam::power_spectrum::OmegaM( zz ) * rho_crit( zz );
  const double rad = std::pow( 0.75 * scam::ip * mm / rho, 0.33333333333333333 );
  
  return 0.33333333333333333 * rad * scam::power_spectrum::dsigma2RdR( rad, zz ) / mm;

}

//==============================================================================================
