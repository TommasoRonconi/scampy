#include <cosmology_interface.h>

//==============================================================================================

scam::cosmo_model::cosmo_model ( const std::unique_ptr< cosmo_p > & parameters,
				 const std::vector< double > & kh0,
				 const std::vector< double > & Pk0,
				 const double zmin, const double zmax,
				 const size_t thin ) :
  param { new cosmo_p { * parameters } },
  z_min { zmin }, z_max { zmax },
  thinness { thin }, 
  P0_f { interp_log { kh0, Pk0 } }
{

  H0 *= param->hh;
  h_1 /= param->hh;
  h_2 *= h_1 * h_1;
  h_3 *= h_2 * h_1;
  scam::cosmo_model::set_internal();
  scam::cosmo_model::compute_ss8_P0();

}

//==============================================================================================

void scam::cosmo_model::set_internal () {

  // when z_min = 0. resets it to a small value greater than 0.
  // ( necessary for logarithmic interpolation )
  if ( ! ( z_min > 0. ) ) z_min = 1.e-7;

  // builds an interpolated function object for E( z )
  auto f_Ez = [ & ] ( double zz ) { return std::sqrt( Ez2( zz ) ); };
  Ez_f = interp_log { f_Ez, z_min, z_max, thinness };

  // builds an interpolated function object for 1/E( z )
  // without re-computing E( z ).
  // ( useful for most of the cosmographic applications of E( z ) 
  std::vector< double > xv = Ez_f.get_xv(), fv;
  fv.reserve( thinness );
  for ( auto && ez : Ez_f.get_fv() ) fv.push_back( 1 / ez );
  zE_f = interp_log { xv, fv };

  // computes the Hubble time in yr/h
  t_H0 = 1.e+3 * scam::pc / ( scam::yr * H0 );
  
  // computes the Hubble-orizon distance in Mpc/h
  d_H0 = 1.e-3 * scam::cc / H0;

  return;

}

//==============================================================================================

void scam::cosmo_model::compute_ss8_P0 () {

  if ( ss8_P0 < 0. ) {

    auto f_integrand = [ & ] ( double kk ) {

      double WF = scam::utl::TopHat_WF( 8 * kk );
      return kk * kk * P0_f( kk ) * WF * WF;

    };

    // scam::cosmo_model::interp_log int_f ( f_integrand,
    // 					  P0_f.get_xmin(),
    // 					  P0_f.get_xmax(),
    // 					  thinness );
    scam::cosmo_model::interp_log int_f ( f_integrand,
    					  1.e-4,
    					  100.,
    					  thinness );

    ss8_P0 = 0.5 * scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

  }

  return;

}

//==============================================================================================

double scam::cosmo_model::Ea2 ( const double zz ) {

  const double aa = 1 / ( 1 + zz );

  double sum =					\
    ( ( param->Om_L * aa			\
	* aa + param->Om_K )			\
      * aa + param->Om_M )			\
    * aa + ( param->Om_r + param->Om_n );
  
  return sum;

}

//==============================================================================================

double scam::cosmo_model::Ez2 ( const double zz ) {

  const double red = 1 + zz;
  double sum =				  \
    ( ( param->Om_r + param->Om_n ) * red \
      * red + param->Om_M )		  \
    * red + param->Om_K;

  return red * red * sum + param->Om_L;

}

//==============================================================================================

double scam::cosmo_model::comoving_volume_unit ( const double zz ) {

  const double zE = zE_f.integrate( 1.e-7, zz );

  return d_H0 * d_H0 * d_H0 * zE_f( zz ) * zE * zE;

}

//==============================================================================================

double scam::cosmo_model::comoving_volume ( const double zz ) {

  const double DC = scam::cosmo_model::d_C( zz );
  
  return 1.33333333333 * scam::pi * DC * DC * DC;
  
}

//==============================================================================================

double scam::cosmo_model::cosmic_time ( const double zz ) {

  std::vector< double > zv = zE_f.get_xv();
  std::vector< double > fv = zE_f.get_fv();

  for ( size_t ii = 0; ii < zE_f.get_thinness(); ++ii ) fv[ ii ] /= ( 1 + zv[ ii ] );

  scam::cosmo_model::interp_log func { zv, fv };

  return 1.e-9 * t_H0 * func.integrate( zz, 1.e+7 );

}

//==============================================================================================

double scam::cosmo_model::rho_crit_comoving ( const double zz ) {

  return 3.75e+18 * scam::ip * scam::nG * scam::solM * scam::pc * Ez2( zz );

}

//==============================================================================================

double scam::cosmo_model::deltac ( const double & zz ) {

  const double Om_Mz = OmegaM( zz );
  
  return 1.686 * ( 1 + 0.012299 * std::log10( Om_Mz ) );

}

//==============================================================================================

double scam::cosmo_model::Delta_c_BryanNorman98 ( const double & zz ) {

  const double xx = OmegaM( zz ) - 1;
    
  return ( -39 * xx + 82 ) * xx + 18 * scam::pi * scam::pi;
  
}

//==============================================================================================

double scam::cosmo_model::Delta_c_NakamuraSuto98 ( const double & zz ) {

  const double xx = std::pow( 1 - param->Om_M, 0.33333333333333 ) / ( 1 + zz );

  return 18 * scam::pi * ( 1 + 0.4093 * std::pow( xx, 2.7152 ) ) * OmegaM( zz );
  
}

//==============================================================================================

double scam::cosmo_model::DD ( const double & zz ) noexcept {
  
  std::vector< double > zv = zE_f.get_xv();
  std::vector< double > fv ( zE_f.get_thinness() );

  for ( size_t ii = 0; ii < zE_f.get_thinness(); ++ii ) fv[ ii ] = 1 + zv[ ii ];
  
  scam::cosmo_model::interp_log func { zv, fv };
  func *= zE_f * zE_f * zE_f;

  return 2.5 * param->Om_M * Ez_f( zz ) * func.integrate( zz, 1.e+7 );

}

//==============================================================================================

double scam::cosmo_model::Pk_comoving ( const double & kk, const double & zz ) {

  const double Dz = DD( zz );
  const double D0 = DD( 1.e-7 );
  
  return param->sigma8 * param->sigma8 * Dz * Dz * P0_f( kk ) / ( D0 * D0 * ss8_P0 );

}

double scam::cosmo_model::Pk ( const double & kk, const double & zz ) {

  const double Dz = DD( zz );
  const double D0 = DD( 1.e-7 );
  
  return h_3 * param->sigma8 * param->sigma8 * Dz * Dz * P0_f( kk * h_1 ) / ( D0 * D0 * ss8_P0 );

}

//==============================================================================================

double scam::cosmo_model::sigma2R_comoving ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    return  kk * kk * Pk_comoving( kk, zz ) * WF * WF;

  };
  
  // scam::cosmo_model::interp_log int_f ( f_integrand,
  // 					P0_f.get_xmin(),
  // 					0.1 * P0_f.get_xmax(),
  // 					thinness );  
  scam::cosmo_model::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );
  
  // return int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );
  return 0.5 * scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}


double scam::cosmo_model::sigma2R ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    return  kk * kk * Pk( kk, zz ) * WF * WF;

  };
  
  // scam::cosmo_model::interp_log int_f ( f_integrand,
  // 					P0_f.get_xmin(),
  // 					0.1 * P0_f.get_xmax(),
  // 					thinness );  
  scam::cosmo_model::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );
  
  // return int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );
  return 0.5 * scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}

//==============================================================================================

double scam::cosmo_model::dsigma2RdR_comoving ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    double dWF = scam::utl::TopHat_WF_D1( rr * kk );
    // return 2 * kk * kk * kk * Pk_comoving( kk, zz ) * WF * dWF;
    return kk * kk * kk * Pk_comoving( kk, zz ) * WF * dWF;

  };
  
  // scam::cosmo_model::interp_log int_f ( f_integrand,
  // 					P0_f.get_xmin(),
  // 					0.1 * P0_f.get_xmax(),
  // 					thinness );  
  scam::cosmo_model::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );

  return scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}

double scam::cosmo_model::dsigma2RdR ( const double & rr, const double & zz ) {

  auto f_integrand = [ & ] ( double kk ) {

    double WF = scam::utl::TopHat_WF( rr * kk );
    double dWF = scam::utl::TopHat_WF_D1( rr * kk );
    return  kk * kk * kk * Pk( kk, zz ) * WF * dWF;

  };
  
  // scam::cosmo_model::interp_log int_f ( f_integrand,
  // 					P0_f.get_xmin(),
  // 					0.1 * P0_f.get_xmax(),
  // 					thinness );  
  scam::cosmo_model::interp_log int_f ( f_integrand,
					1.e-4,
				        100.,
					thinness );

  // return int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );
  return scam::ip * scam::ip * int_f.integrate( int_f.get_xmin(), int_f.get_xmax() );

}

//==============================================================================================

double scam::cosmo_model::sigma2M ( const double & mm, const double & zz ) {

  const double rho = scam::cosmo_model::OmegaM( 1.e-7 ) * rho_crit( 1.e-7 );
  const double rad = std::pow( 0.75 * scam::ip * mm / rho, 0.33333333333333333 );
  
  return scam::cosmo_model::sigma2R( rad, zz );

}

//==============================================================================================

double scam::cosmo_model::dsigma2MdM ( const double & mm, const double & zz ) {

  const double rho = scam::cosmo_model::OmegaM( zz ) * rho_crit( zz );
  const double rad = std::pow( 0.75 * scam::ip * mm / rho, 0.33333333333333333 );
  
  return 0.33333333333333333 * rad * scam::cosmo_model::dsigma2RdR( rad, zz ) / mm;

}

//==============================================================================================

double scam::cosmo_model::dndM_ST01 ( const double & mm, const double & zz ) {

  const double gfact = scam::cosmo_model::DD( 1.e-7 ) / scam::cosmo_model::DD( zz );
  const double ss = std::sqrt( scam::cosmo_model::sigma2M( mm, 1.e-7 ) );
  const double nu = gfact * scam::cosmo_model::deltac( zz ) / ss;

  const double fnu =					\
    0.84083292 * std::sqrt( 2 * scam::ip ) * 0.3222 *	\
    ( 1 + std::pow( 0.84083292 * nu, - 0.6 ) ) * nu *	\
    std::exp( - 0.5 * 0.707 * nu * nu );
  
  const double dls =				\
    - 0.5 * scam::cosmo_model::dsigma2MdM( mm ) \
    / scam::cosmo_model::sigma2M( mm );
  
  const double rho0 = param->Om_M * rho_crit( 1.e-7 );

  return rho0 * dls * fnu / mm * h_1;

}

double scam::cosmo_model::dndM_Tinker08 ( const double & mm, const double & zz ) {

  const double mass = mm * h_1;
  const double gfact = scam::cosmo_model::DD( zz ) / scam::cosmo_model::DD( 1.e-7 );
  const double rho0 = param->Om_M * scam::cosmo_model::rho_crit_comoving( 1.e-7 );
  const double rad = std::pow( 0.75 * scam::ip * mass / rho0, 0.33333333333333333 );
  const double s2 = scam::cosmo_model::sigma2R_comoving( rad, 1.e-7 );
  const double ss = std::sqrt( s2 );
  const double dls =						\
    - 0.5 * 0.33333333333333333 * rad *				\
    scam::cosmo_model::dsigma2RdR_comoving( rad, 1.e-7 )	\
    / ( s2 * mass );
  
  const double nu = gfact * ss;
  const double Delta = 200, A0 = 0.186, a0 = 1.47, b0 = 2.57, c0 = 1.19;
  const double alpha = 1.0676e-2; //<= std::pow( 10., -std::pow( 0.75/std::log10( Delta/75 ), 1.2 ) )
  const double Az = A0 * std::pow( 1 + zz, - 0.14 );
  const double az = a0 * std::pow( 1 + zz, - 0.06 );
  const double bz = b0 * std::pow( 1 + zz, - alpha );

  const double fnu =				\
    Az * ( std::pow( nu / bz, - az ) + 1 ) *	\
    std::exp( - c0 / ( nu * nu ) );

  return rho0 * dls * fnu / mass;

}

double scam::cosmo_model::dndM_Behroozi13 ( const double & mm, const double & zz ) {

  const double ascale = 1. / ( 1 + zz );
  const double logNorm = 0.144 / ( 1 + std::exp( 14.79 * ( ascale - 0.213 ) ) );
  const double high_z_correction = std::pow( mm * 3.162278e-12,
					     0.5 / ( 1 + std::exp( 6.5 * ascale ) ) );

  return std::pow( 10., logNorm * high_z_correction ) * dndM_Tinker08( mm, zz );

}

//==============================================================================================

double scam::cosmo_model::hbias_SMT01 ( const double & mm, const double & zz ) {
  
  const double dd = scam::cosmo_model::deltac( zz );
  const double gfact = scam::cosmo_model::DD( 1.e-7 ) / scam::cosmo_model::DD( zz );
  const double nu2 = dd * dd / scam::cosmo_model::sigma2M( mm, 1.e-7 ) * gfact * gfact; 
  const double an2 = 0.707 * nu2;
  const double an2p = std::pow( an2, 0.6 );
  
  return 1. + 1. / dd *					\
    ( an2 + 0.5 * an2 / an2p -				\
      an2p / ( 0.84084292 * ( an2p + 0.14 ) ) );
  
}

//==============================================================================================

double scam::cosmo_model::hbias_Tinker10 ( const double & mm, const double & zz ) {
  
  const double dd = scam::cosmo_model::deltac( zz );
  const double gfact = scam::cosmo_model::DD( 1.e-7 ) / scam::cosmo_model::DD( zz );
  const double nu = dd / std::sqrt( scam::cosmo_model::sigma2M( mm, 1.e-7 ) ) * gfact;
  const double yy = 2.301; //<= std::log10( Delta = 200 )
  const double exy = std::exp( - 256 / ( yy * yy * yy * yy ) );
  const double AA = 1 + 0.24 * yy * exy;
  const double aa = 0.44 * yy - 0.88;
  const double BB = 0.183;
  const double bb = 1.5;
  const double CC = 1.9e-2 + 0.107 * yy + 0.19 * exy;
  const double cc = 2.4;

  return 1 - AA *							\
    std::pow( nu, aa ) / ( std::pow( nu, aa ) + std::pow( dd, aa ) ) +	\
    BB * std::pow( nu, bb ) + CC * std::pow( nu, cc );
  
}

//==============================================================================================
  
// double scam::cosmo_model::EzDE2 ( const double zz ) {

//   const double red = 1 + zz;
//   double sum = ( ( Om_r + Om_n ) * red  * red + Om_M ) * red + Om_K;
//   const double denerg = Om_L * std::pow( red, 1 + 3 * ( w_0 + w_a * zz / red ) );

//   return red * red * ( sum + Om_L ); 

// }

//==============================================================================================
