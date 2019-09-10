#include <error_handling.h>
#include <halo_model.h>
#include <functional>

//==============================================================================================

sico::halo_model::halo_model ( sico::halo_model_handler handler ) : _handler{ handler } {

  std::cout << "Computing interpolated functions ...";
  auto f_dndM = [ & ] ( double Mh ) { return _handler.cosmo.dndM( Mh, _handler.redshift ); };
  dndM_f = sico::halo_model::interp_func { f_dndM,
					   _handler.mass_integ_lim.inf,
					   _handler.mass_integ_lim.sup,
					   _handler.thinness };
  auto f_hbias = [ & ] ( double Mh ) { return _handler.cosmo.hbias( Mh, _handler.redshift ); };
  hbias_f = sico::halo_model::interp_func { f_hbias,
					    _handler.mass_integ_lim.inf,
					    _handler.mass_integ_lim.sup,
					    _handler.thinness };
  auto f_Ncen = [ & ] ( double Mh ) { return sico::halo_model::Ncen( Mh ); };
  Ncen_f = sico::halo_model::interp_func { f_Ncen,
					   _handler.mass_integ_lim.inf,
					   _handler.mass_integ_lim.sup,
					   _handler.thinness };
  auto f_Nsat = [ & ] ( double Mh ) { return sico::halo_model::Nsat( Mh ); };
  Nsat_f = sico::halo_model::interp_func { f_Nsat,
					   _handler.mass_integ_lim.inf,
					   _handler.mass_integ_lim.sup,
					   _handler.thinness };

  _Mv = Ncen_f.get_xv();
  std::vector< double > f_DC ( _handler.thinness, _handler.DC );
  DC_f = sico::halo_model::interp_func { _Mv, f_DC };
  std::vector< double > f_const2 ( _handler.thinness, 2 );
  const2_f = sico::halo_model::interp_func { _Mv, f_const2 };

  _kv = sico::utl::log_vector( _handler.thinness,
  			       _handler.wavk_integ_lim.inf,
  			       _handler.wavk_integ_lim.sup );

  density_profile_FS = std::vector< utl::interpolator< utl::gsl_log_interp > > { _handler.thinness };

#pragma omp parallel for
  for ( size_t _k = 0; _k < _handler.thinness; ++_k ) {
    
    auto f_denFS = [ & ] ( double Mh ) {
      return _handler.cosmo.density_profile_FS( _kv[_k], Mh, 0. );
    };
    
    density_profile_FS[ _k ] =
      sico::halo_model::interp_func {
      f_denFS,
      _handler.mass_integ_lim.inf,
      _handler.mass_integ_lim.sup,
      _handler.thinness
    };    
    
  }
    
  std::cout << " done." << std::endl;
  
}

//==============================================================================================

void sico::halo_model::set_parameters ( const double DC,
					const double Mmin,
					const double sigma_logM,
					const double M0,
					const double M1,
					const double alpha ) {

  _handler.DC = DC;
  _handler.M_min = Mmin;
  _handler.sigma_logM = sigma_logM;
  _handler.M0 = M0;
  _handler.M1 = M1;
  _handler.alpha = alpha;
  auto f_Ncen = [ & ] ( double Mh ) { return sico::halo_model::Ncen( Mh ); };
  Ncen_f = sico::halo_model::interp_func { f_Ncen,
					   _handler.mass_integ_lim.inf,
					   _handler.mass_integ_lim.sup,
					   _handler.thinness };
  auto f_Nsat = [ & ] ( double Mh ) { return sico::halo_model::Nsat( Mh ); };
  Nsat_f = sico::halo_model::interp_func { f_Nsat,
					   _handler.mass_integ_lim.inf,
					   _handler.mass_integ_lim.sup,
					   _handler.thinness };
  std::vector< double > f_DC ( _handler.thinness, _handler.DC );
  DC_f = sico::halo_model::interp_func { _Mv, f_DC };

  return;
  
}

//==============================================================================================

double sico::halo_model::PP (const double AA, const double Amin, const double sigma_logA) {

  double xx = ( std::log10( AA ) - std::log10( Amin ) ) / sigma_logA;
  double erf = gsl_sf_erf( xx );

  return 0.5 * ( 1. + erf );

}

//==============================================================================================


double sico::halo_model::Ncen ( const double Mhalo ) {
  
  double Nc = PP( Mhalo, _handler.M_min, _handler.sigma_logM );
  
  return ( Nc < 0 || std::isnan( Nc ) ) ? 0. : Nc;
  
}


//==============================================================================================

double sico::halo_model::Nsat ( const double Mhalo ) {

  const double Nc = sico::halo_model::Ncen( Mhalo );
  const double Ns = Nc * std::pow( ( ( Mhalo - _handler.M0 ) / _handler.M1 ), _handler.alpha );
  
  return ( Ns < 0 || std::isnan( Ns ) ) ? 0. : Ns;
  
}

//==============================================================================================

double sico::halo_model::ng () {
  
  auto integrand = ( Ncen_f + Nsat_f ) * DC_f * dndM_f;
    
  // integration limits of ng are from 0. up to some Mmax
  // substitute 0 with m_Mh_min? ----------------------v ?!
  return integrand.integrate( _handler.mass_integ_lim.inf, _handler.mass_integ_lim.sup );

}

// ============================================================================================

double sico::halo_model::bias () {

  auto integrand = ( Ncen_f + Nsat_f ) * DC_f * dndM_f * hbias_f;

  return 1 / ng() * integrand.integrate( _handler.mass_integ_lim.inf, _handler.mass_integ_lim.sup );
  
}

// ============================================================================================

double sico::halo_model::Mhalo () {

  sico::halo_model::interp_func Mh_f { _Mv, _Mv };
  auto integrand = ( Ncen_f + Nsat_f ) * DC_f * dndM_f * Mh_f;

  return 1 / ng() * integrand.integrate( _handler.mass_integ_lim.inf, _handler.mass_integ_lim.sup );
  
}

// ============================================================================================

double sico::halo_model::dngdM ( const double Mhalo ) {

  return ( Ncen( Mhalo ) + Nsat( Mhalo ) ) * _handler.DC * _handler.cosmo.dndM( Mhalo );

}

// ============================================================================================

double sico::halo_model::Pk_cc_integrand ( const double Mh, const double kk ) {
  
  // mean number of central-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = Ncen( Mh );

  // halo mass function -> it depends on cosmology
  const double dndM = _handler.cosmo.dndM( Mh );
  
  // density profile -> it depends on cosmology
  const double uk = _handler.cosmo.density_profile_FS( kk, Mh, _handler.redshift );
  const double ukp = ( NN * NN > 1 ) ? uk * uk : uk;
  
  return NN * NN * dndM * ukp;

}

// ============================================================================================

double sico::halo_model::Pk_cs_integrand ( const double Mh, const double kk ) {
  
  // mean number of central-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = Ncen( Mh ) * Nsat( Mh );

  // halo mass function -> it depends on cosmology
  const double dndM = _handler.cosmo.dndM( Mh );
  
  // density profile -> it depends on cosmology
  const double uk = _handler.cosmo.density_profile_FS( kk, Mh, _handler.redshift );
  const double ukp = ( NN > 1 ) ? uk * uk : uk;
  
  return NN * dndM * ukp;

}

// ============================================================================================

double sico::halo_model::Pk_ss_integrand (const double Mh, const double kk)
{
  
  // mean number of satellite-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = Nsat( Mh );

  // halo mass function -> it depends on cosmology
  const double dndM = _handler.cosmo.dndM( Mh );
      
  // density profile -> it depends on cosmology
  const double uk = _handler.cosmo.density_profile_FS( kk, Mh, _handler.redshift );
  const double ukp = ( NN * NN > 1 ) ? uk * uk : uk;
      
  return NN * NN * dndM * ukp;
  
}

// ============================================================================================

double sico::halo_model::Pk_1halo ( const size_t ii, const double fact_ng2 ) {

  // // density profile -> it depends on cosmology
  // const double uk = _handler.cosmo.density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    DC_f * DC_f							\
    * ( const2_f * Ncen_f + Nsat_f )				\
    * Nsat_f * dndM_f						\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact_ng2 * integrand.integrate( _handler.mass_integ_lim.inf,
					 _handler.mass_integ_lim.sup );

}

double sico::halo_model::Pk_cs ( const size_t ii, const double fact_ng2 ) {

  // // density profile -> it depends on cosmology
  // const double uk = _handler.cosmo.density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    DC_f * DC_f * const2_f					\
    * Ncen_f * Nsat_f * dndM_f					\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact_ng2 * integrand.integrate( _handler.mass_integ_lim.inf,
					 _handler.mass_integ_lim.sup );

}

double sico::halo_model::Pk_ss ( const size_t ii, const double fact_ng2 ) {

  // // density profile -> it depends on cosmology
  // const double uk = _handler.cosmo.density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    DC_f * DC_f							\
    * Nsat_f * Nsat_f * dndM_f					\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact_ng2 * integrand.integrate( _handler.mass_integ_lim.inf,
					 _handler.mass_integ_lim.sup );

}
  
// ============================================================================================

double sico::halo_model::Pk_2halo ( const size_t ii, const double fact_ng2 ) {

  auto integrand = ( Ncen_f + Nsat_f ) * DC_f * dndM_f * hbias_f * density_profile_FS[ ii ];

  double integral = integrand.integrate( _handler.mass_integ_lim.inf,
					 _handler.mass_integ_lim.sup );
  
  return _handler.cosmo.Pk( _kv[ ii ], _handler.redshift ) * fact_ng2 * integral * integral;

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Pk_1halo () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _handler.thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _handler.thinness; ++ii )
      Pk_gxy[ ii ] = Pk_1halo( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

std::vector< double > sico::halo_model::model_Pk_cs () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _handler.thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _handler.thinness; ++ii )
      Pk_gxy[ ii ] = Pk_cs( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

std::vector< double > sico::halo_model::model_Pk_ss () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _handler.thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _handler.thinness; ++ii )
      Pk_gxy[ ii ] = Pk_ss( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > sico::halo_model::model_Pk_2halo () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _handler.thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _handler.thinness; ++ii )
      Pk_gxy[ ii ] = Pk_2halo( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > sico::halo_model::model_Pk () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _handler.thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _handler.thinness; ++ii )
      Pk_gxy[ ii ] = Pk_1halo( ii, fact_ng2 ) + Pk_2halo( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > sico::halo_model::model_Xi_1halo ( const std::vector< double > & rad ) {

  std::vector< double > Pk = sico::halo_model::model_Pk_1halo();
  sico::utl::fftlog_3Dspace fft { _kv, Pk };

  return fft.transform( rad );

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Xi_2halo ( const std::vector< double > & rad ) {

  std::vector< double > Pk = sico::halo_model::model_Pk_2halo();
  sico::utl::fftlog_3Dspace fft { _kv, Pk };

  return fft.transform( rad );

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Xi ( const std::vector< double > & rad ) {

  std::vector< double > Pk = sico::halo_model::model_Pk();
  double rM = std::exp( 0.5 * ( std::log( rad.back() ) + std::log( rad.front() ) ) );
  double kM = std::exp( 0.5 * ( std::log( _kv.back() ) + std::log( _kv.front() ) ) );
  sico::utl::fftlog_3Dspace fft { _kv, Pk, rM * kM };
  
  return fft.transform( rad );

}

//==============================================================================================
    
std::vector< double > sico::halo_model::model_Wr_1halo ( const std::vector< double > & radp ) {

  std::vector< double > Pk = sico::halo_model::model_Pk_1halo();
  sico::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}
    
std::vector< double > sico::halo_model::model_Wr_cs ( const std::vector< double > & radp ) {

  std::vector< double > Pk = sico::halo_model::model_Pk_cs();
  sico::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}
    
std::vector< double > sico::halo_model::model_Wr_ss ( const std::vector< double > & radp ) {

  std::vector< double > Pk = sico::halo_model::model_Pk_ss();
  sico::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}

//==============================================================================================
    
std::vector< double > sico::halo_model::model_Wr_2halo ( const std::vector< double > & radp ) {

  std::vector< double > Pk = sico::halo_model::model_Pk_2halo();
  sico::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}

//==============================================================================================
    
std::vector< double > sico::halo_model::model_Wr ( const std::vector< double > & radp ) {

  std::vector< double > Pk = sico::halo_model::model_Pk();
  sico::utl::fftlog_projected fft { _kv, Pk };
  return fft.transform( radp );

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Wt_1halo ( const std::vector< double > & theta ) {

  const double r_z = _handler.cosmo.d_C( _handler.redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return sico::halo_model::model_Wr_1halo( radp );

}

std::vector< double > sico::halo_model::model_Wt_cs ( const std::vector< double > & theta ) {

  const double r_z = _handler.cosmo.d_C( _handler.redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return sico::halo_model::model_Wr_cs( radp );

}

std::vector< double > sico::halo_model::model_Wt_ss ( const std::vector< double > & theta ) {

  const double r_z = _handler.cosmo.d_C( _handler.redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return sico::halo_model::model_Wr_ss( radp );

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Wt_2halo ( const std::vector< double > & theta ) {

  const double r_z = _handler.cosmo.d_C( _handler.redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return sico::halo_model::model_Wr_2halo( radp );

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Wt ( const std::vector< double > & theta ) {

  const double r_z = _handler.cosmo.d_C( _handler.redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return sico::halo_model::model_Wr( radp );

}

//==============================================================================================

std::vector< double > sico::halo_model::model_Pk_large_scale () {

  const double b_gal = bias();
  const double b2 = b_gal * b_gal;
  
  std::vector< double > Pk_gxy( _handler.thinness );  
#pragma omp parallel for
  for ( size_t ii = 0; ii < _handler.thinness; ++ii )
    Pk_gxy[ ii ] = b2 * _handler.cosmo.Pk( _kv[ ii ], _handler.redshift );
    
  return Pk_gxy;
  
}

//==============================================================================================
    
std::vector< double > sico::halo_model::model_Wt_large_scale ( const std::vector< double > & theta ) {

  const double r_z = _handler.cosmo.d_C( _handler.redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  std::vector< double > Pk = sico::halo_model::model_Pk_large_scale();
    
  sico::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}

//==============================================================================================
