#include <cosmological_model.h>
#include <iostream>
//==============================================================================================

scam::cosmo_model::cosmo_model ( const double Om_M,
				 const double Om_b,
				 const double Om_L,
				 const double Om_n,
				 const double Om_r,
				 const double Om_K,
				 const double hh,
				 const double sigma8,
				 const double w_0,
				 const double w_a,
				 const double zmin,
				 const double zmax,
				 const size_t thin ) :
  z_min { zmin }, z_max { zmax },
  thinness { thin } 
{

  // build map
  param["Om_M"]   = Om_M;
  param["Om_b"]   = Om_b;
  param["Om_L"]   = Om_L;
  param["Om_n"]   = Om_n;
  param["Om_r"]   = Om_r;
  param["Om_K"]   = Om_K;
  param["hh"]     = hh;
  param["sigma8"] = sigma8;
  param["w_0"]    = w_0;
  param["w_a"]    = w_a;
  param["zmin"]   = zmin;
  param["zmax"]   = zmax;
  param["thin"]   = thin;
  
  // build associated quantities
  H0 *= param["hh"];
  scam::cosmo_model::set_internal();

}

//==============================================================================================

void scam::cosmo_model::set_internal () {

  // when z_min = 0. resets it to a small value greater than 0.
  // ( necessary for logarithmic interpolation )
  std::vector< double > zz;
  if ( z_min > 0. ) {
    zz = utl::log_vector( thinness, z_min, z_max );
  }
  else {
    zz = utl::log_vector( thinness, 1.e-3, z_max );
    zz[ 0 ] = 0.0;
  }

  // builds an interpolated function object for E( z )
  std::vector< double > Ez; Ez.reserve( thinness );
  for ( auto && _z : zz ) 
    Ez.emplace_back( std::sqrt( Ez2( _z ) ) );
  Ez_f = interp_lin { zz, Ez };

  // builds an interpolated function object for 1/E( z )
  // without re-computing E( z ).
  // ( useful for most of the cosmographic applications of E( z ) 
  std::vector< double > zE;
  zE.reserve( thinness );
  for ( auto && ez : Ez_f.get_fv() ) zE.emplace_back( 1 / ez );
  zE_f = interp_lin { zz, zE };

  // computes the Hubble time in yr
  t_H0 = 1.e+3 * utl::cnst::pc / ( utl::cnst::yr * H0 );
  
  // computes the Hubble-orizon distance in Mpc
  d_H0 = 1.e-3 * utl::cnst::cc / H0;

  return;

}

//==============================================================================================

double scam::cosmo_model::Ea2 ( const double zz ) {

  const double aa = 1 / ( 1 + zz );

  double sum =					\
    ( ( param["Om_L"] * aa			\
	* aa + param["Om_K"] )			\
      * aa + param["Om_M"] )			\
    * aa + ( param["Om_r"] + param["Om_n"] );
  
  return sum;

}

//==============================================================================================

double scam::cosmo_model::Ez2 ( const double zz ) {

  const double red = 1 + zz;
  double sum =				      \
    ( ( param["Om_r"] + param["Om_n"] )	      \
      * red + param["Om_M"] )		      \
    * red + param["Om_K"];

  return red * red * sum + param["Om_L"];

}

//==============================================================================================

double scam::cosmo_model::comoving_volume_unit ( const double zz ) {

  const double zE = zE_f.integrate( 0., zz );

  return d_H0 * d_H0 * d_H0 * zE_f( zz ) * zE * zE;

}

//==============================================================================================

double scam::cosmo_model::comoving_volume ( const double zz ) {

  const double DC = scam::cosmo_model::d_C( zz );
  
  return 1.33333333333 * utl::cnst::pi * DC * DC * DC;
  
}

//==============================================================================================

double scam::cosmo_model::cosmic_time ( const double zz ) {

  std::vector< double > zv = zE_f.get_xv();
  std::vector< double > fv = zE_f.get_fv();

  for ( size_t ii = 0; ii < zE_f.get_thinness(); ++ii ) fv[ ii ] /= ( 1 + zv[ ii ] );

  scam::cosmo_model::interp_lin func { zv, fv };

  return 1.e-9 * t_H0 * func.integrate( zz, z_max );

}

//==============================================================================================

double scam::cosmo_model::rho_crit_comoving ( const double zz ) {

  return 3.75e+18 * utl::cnst::ip * utl::cnst::nG * utl::cnst::solM * utl::cnst::pc * Ez2( zz );

}

//==============================================================================================

double scam::cosmo_model::deltac ( const double & zz ) {

  const double Om_Mz = OmegaM( zz );
  
  return 1.686 * ( 1 + 0.012299 * std::log10( Om_Mz ) );

}

//==============================================================================================

double scam::cosmo_model::Delta_c_BryanNorman98 ( const double & zz ) {

  const double xx = OmegaM( zz ) - 1;
    
  return ( -39 * xx + 82 ) * xx + 18 * utl::cnst::pi * utl::cnst::pi;
  
}

//==============================================================================================

double scam::cosmo_model::Delta_c_NakamuraSuto98 ( const double & zz ) {

  const double xx = std::pow( 1 - param["Om_M"], 0.33333333333333 ) / ( 1 + zz );

  return 18 * utl::cnst::pi * ( 1 + 0.4093 * std::pow( xx, 2.7152 ) ) * OmegaM( zz );
  
}

//==============================================================================================

double scam::cosmo_model::DD ( const double & zz ) noexcept {
  
  std::vector< double > zv = zE_f.get_xv();
  std::vector< double > zE = zE_f.get_fv();
  std::vector< double > fv ( zE_f.get_thinness() );

  for ( size_t ii = 0; ii < zE_f.get_thinness(); ++ii )
    fv[ ii ] = ( 1 + zv[ ii ] ) * zE[ ii ] * zE[ ii ] * zE[ ii ];
  
  scam::cosmo_model::interp_lin func { zv, fv };

  return 2.5 * param["Om_M"] * Ez_f( zz ) * func.integrate( zz, z_max );

}

//==============================================================================================
