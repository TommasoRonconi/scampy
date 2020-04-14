#include <cosmology.h>
#include <interpolation.h>
#include <gsl_interpolation_interface.h>
#include "gsl/gsl_sf_expint.h"

//==============================================================================================

// model from Giocoli et al., 2012
double scam::cosmology::z_form ( const double & Mass, const double & ff,
				 const double & z_now, const double & z_max ) {
  
  double alpha_f = 0.815 * std::exp( -2. * ff * ff * ff ) * std::pow( ff, - 0.707 );
  double omega_f = std::sqrt( 2. * std::log( 1. + alpha_f ) );
  double sigma2_0 = sigma2M( Mass );
  double sigma2_f = sigma2M( Mass * ff );
  double Z0 = 1.e-7;
  double D0 = DD( Z0 );
  double delta = deltac( Z0 ) * D0 / DD( z_now ) + omega_f * std::sqrt( sigma2_f - sigma2_0 );

  auto func = [ & ] ( double zz ) {
    
    return deltac( Z0 ) * D0 / DD( zz ) - delta;
    
  };

  const unsigned int nthin = 100;
  std::vector< double > zspace = scam::utl::lin_vector( nthin, z_now, z_max );
  std::vector< double > fspace ( nthin );
  for ( unsigned int ii = 0; ii < nthin; ++ii ) fspace[ ii ] = func( zspace[ ii ] );
  scam::utl::interpolator< scam::utl::gsl_lin_interp > z_form_f { fspace, zspace };
  
  // find a zero (0.) of function (func) between lower and upper guess:
  return z_form_f( 0. );

}

//==============================================================================================

// Eq. 13 from Zhao et al., 2009
double scam::cosmology::concentration_zhao09 ( const double Mass, const double Redshift ) {
  
  double t_now = cosmic_time( Redshift );
  double t_04 = cosmic_time( z_form( Mass, 0.04, Redshift, Redshift + 10 ) );

  return 4. * pow( 1. + pow( t_now / ( 3.75 * t_04 ), 8.4 ), 0.125 );

}

// Eq. 4 from Shimizu et al., 2003
double scam::cosmology::concentration_shimizu03 ( const double Mass, const double Redshift ) {

  const double cnorm = 8.;
  return cnorm / ( 1 + Redshift ) * std::pow( 1.0204e-14 * Mass * param->hh, -0.13 );

}

//==============================================================================================

// double scam::cosmology::density_profile_FS ( const double kk, const double Mass,
// 					     const double Redshift ) {

//   double zz = Redshift;

//   const double conc = concentration( Mass, zz );

//   const double V_vir = Mass / ( Delta_c( zz ) * rho_crit( zz ) );

//   const double r_vir = pow( 0.75 * V_vir * scam::ip, 0.333333333333333333333 );
  
//   const double r_s = r_vir / conc;
  
//   const double mu = kk * r_s;

//   const double cosi_part = gsl_sf_Ci( mu + mu * conc ) - gsl_sf_Ci( mu );
//   const double sini_part = gsl_sf_Si( mu + mu * conc ) - gsl_sf_Si( mu );
//   const double sin_part = sin( mu * conc ) / ( mu + mu * conc );

//   const double tmp = cos( mu ) * cosi_part + sin( mu ) * sini_part - sin_part;
  
//   return tmp / ( log( 1. + conc ) - conc / ( 1. + conc ) );
  
// }

//==============================================================================================

double scam::cosmology::density_profile_FS ( const double kk, const double Mass,
					     const double Redshift ) {

  double zz = 1.e-7;

  const double conc = concentration( Mass, Redshift );

  const double V_vir = Mass / ( Delta_c( zz ) * rho_crit( zz ) );

  const double r_vir = pow( 0.75 * V_vir * scam::ip, 0.333333333333333333333 );
  
  const double r_s = r_vir / conc;
  
  const double mu = kk * r_s;

  const double cosi_part = gsl_sf_Ci( mu + mu * conc ) - gsl_sf_Ci( mu );
  const double sini_part = gsl_sf_Si( mu + mu * conc ) - gsl_sf_Si( mu );
  const double sin_part = sin( mu * conc ) / ( mu + mu * conc );

  const double tmp = cos( mu ) * cosi_part + sin( mu ) * sini_part - sin_part;
  
  return tmp / ( log( 1. + conc ) - conc / ( 1. + conc ) );
  
}

//==============================================================================================
