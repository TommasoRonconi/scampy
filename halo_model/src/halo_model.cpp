#include <error_handling.h>
#include <halo_model.h>
#include <functional>

//==============================================================================================

scam::halo_model::halo_model ( const std::shared_ptr< scam::occupation_p > & ocp,
			       const std::shared_ptr< scam::cosmology > & cosmo,
			       const double redshift,
			       const size_t thinness ) : _handler{ ocp },
							 _cosmo{ cosmo },
							 _redshift{ redshift },
							 _thinness{ thinness } {
							     
  std::cout << "Computing interpolated functions ...";
  auto f_dndM = [ & ] ( double Mh ) { return _cosmo->dndM( Mh, _redshift ); };
  dndM_f = scam::halo_model::interp_func { f_dndM,
					   _mass_integ_lim.inf,
					   _mass_integ_lim.sup,
					   _thinness };
  auto f_hbias = [ & ] ( double Mh ) { return _cosmo->hbias( Mh, _redshift ); };
  hbias_f = scam::halo_model::interp_func { f_hbias,
					    _mass_integ_lim.inf,
					    _mass_integ_lim.sup,
					    _thinness };
  auto f_Ncen = [ & ] ( double Mh ) { return _handler->Ncen( Mh ); };
  Ncen_f = scam::halo_model::interp_func { f_Ncen,
					   _mass_integ_lim.inf,
					   _mass_integ_lim.sup,
					   _thinness };
  auto f_Nsat = [ & ] ( double Mh ) { return _handler->Nsat( Mh ); };
  Nsat_f = scam::halo_model::interp_func { f_Nsat,
					   _mass_integ_lim.inf,
					   _mass_integ_lim.sup,
					   _thinness };

  _Mv = Ncen_f.get_xv();

  std::vector< double > f_const2 ( _thinness, 2 );
  const2_f = scam::halo_model::interp_func { _Mv, f_const2 };

  _kv = scam::utl::log_vector( _thinness,
  			       _wavk_integ_lim.inf,
  			       _wavk_integ_lim.sup );

  density_profile_FS = std::vector< utl::interpolator< utl::gsl_log_interp > > { _thinness };

#pragma omp parallel for
  for ( size_t _k = 0; _k < _thinness; ++_k ) {
    
    auto f_denFS = [ & ] ( double Mh ) {
      return _cosmo->density_profile_FS( _kv[_k], Mh, 0. );
    };
    
    density_profile_FS[ _k ] =
      scam::halo_model::interp_func {
      f_denFS,
      _mass_integ_lim.inf,
      _mass_integ_lim.sup,
      _thinness
    };    
    
  }
    
  std::cout << " done." << std::endl;
  
}

//==============================================================================================

void scam::halo_model::set_parameters ( const std::shared_ptr< scam::occupation_p > & ocp ) {

  _handler = ocp;
  // _handler = std::unique_ptr< scam::occupation_p >{ ocp };
  
  auto f_Ncen = [ & ] ( double Mh ) { return _handler->Ncen( Mh ); };
  Ncen_f = scam::halo_model::interp_func { f_Ncen,
					   _mass_integ_lim.inf,
					   _mass_integ_lim.sup,
					   _thinness };
  auto f_Nsat = [ & ] ( double Mh ) { return _handler->Nsat( Mh ); };
  Nsat_f = scam::halo_model::interp_func { f_Nsat,
					   _mass_integ_lim.inf,
					   _mass_integ_lim.sup,
					   _thinness };

  return;
  
}

//==============================================================================================

double scam::halo_model::ng () {
  
  auto integrand = ( Ncen_f + Nsat_f ) * dndM_f;
    
  // integration limits of ng are from 0. up to some Mmax
  // substitute 0 with m_Mh_min? ----------------------v ?!
  return integrand.integrate( _mass_integ_lim.inf, _mass_integ_lim.sup );

}

// ============================================================================================

double scam::halo_model::bias () {

  auto integrand = ( Ncen_f + Nsat_f ) * dndM_f * hbias_f;

  return 1 / ng() * integrand.integrate( _mass_integ_lim.inf, _mass_integ_lim.sup );
  
}

// ============================================================================================

double scam::halo_model::Mhalo () {

  scam::halo_model::interp_func Mh_f { _Mv, _Mv };
  auto integrand = ( Ncen_f + Nsat_f ) * dndM_f * Mh_f;

  return 1 / ng() * integrand.integrate( _mass_integ_lim.inf, _mass_integ_lim.sup );
  
}

// ============================================================================================

double scam::halo_model::dngdM ( const double Mhalo ) {

  return ( Ncen_f( Mhalo ) + Nsat_f( Mhalo ) ) * _cosmo->dndM( Mhalo );

}

// ============================================================================================

double scam::halo_model::Pk_cc_integrand ( const double Mh, const double kk ) {
  
  // mean number of central-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = Ncen_f( Mh );

  // halo mass function -> it depends on cosmology
  const double dndM = _cosmo->dndM( Mh );
  
  // density profile -> it depends on cosmology
  const double uk = _cosmo->density_profile_FS( kk, Mh, _redshift );
  const double ukp = ( NN * NN > 1 ) ? uk * uk : uk;
  
  return NN * NN * dndM * ukp;

}

// ============================================================================================

double scam::halo_model::Pk_cs_integrand ( const double Mh, const double kk ) {
  
  // mean number of central-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = Ncen_f( Mh ) * Nsat_f( Mh );

  // halo mass function -> it depends on cosmology
  const double dndM = _cosmo->dndM( Mh );
  
  // density profile -> it depends on cosmology
  const double uk = _cosmo->density_profile_FS( kk, Mh, _redshift );
  const double ukp = ( NN > 1 ) ? uk * uk : uk;
  
  return NN * dndM * ukp;

}

// ============================================================================================

double scam::halo_model::Pk_ss_integrand (const double Mh, const double kk)
{
  
  // mean number of satellite-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = Nsat_f( Mh );

  // halo mass function -> it depends on cosmology
  const double dndM = _cosmo->dndM( Mh );
      
  // density profile -> it depends on cosmology
  const double uk = _cosmo->density_profile_FS( kk, Mh, _redshift );
  const double ukp = ( NN * NN > 1 ) ? uk * uk : uk;
      
  return NN * NN * dndM * ukp;
  
}

// ============================================================================================

double scam::halo_model::Pk_1halo ( const size_t ii, const double fact_ng2 ) {

  // // density profile -> it depends on cosmology
  // const double uk = _cosmo->density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    ( const2_f * Ncen_f + Nsat_f )				\
    * Nsat_f * dndM_f						\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact_ng2 * integrand.integrate( _mass_integ_lim.inf,
					 _mass_integ_lim.sup );

}

double scam::halo_model::Pk_cs ( const size_t ii, const double fact_ng2 ) {

  // // density profile -> it depends on cosmology
  // const double uk = _cosmo->density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    const2_f							\
    * Ncen_f * Nsat_f * dndM_f					\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact_ng2 * integrand.integrate( _mass_integ_lim.inf,
					 _mass_integ_lim.sup );

}

double scam::halo_model::Pk_ss ( const size_t ii, const double fact_ng2 ) {

  // // density profile -> it depends on cosmology
  // const double uk = _cosmo->density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    Nsat_f * Nsat_f * dndM_f					\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact_ng2 * integrand.integrate( _mass_integ_lim.inf,
					 _mass_integ_lim.sup );

}
  
// ============================================================================================

double scam::halo_model::Pk_2halo ( const size_t ii, const double fact_ng2 ) {

  auto integrand = ( Ncen_f + Nsat_f ) * dndM_f * hbias_f * density_profile_FS[ ii ];

  double integral = integrand.integrate( _mass_integ_lim.inf,
					 _mass_integ_lim.sup );
  
  return _cosmo->Pk( _kv[ ii ], _redshift ) * fact_ng2 * integral * integral;

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Pk_1halo () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_1halo( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

std::vector< double > scam::halo_model::model_Pk_cs () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_cs( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

std::vector< double > scam::halo_model::model_Pk_ss () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_ss( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > scam::halo_model::model_Pk_2halo () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_2halo( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > scam::halo_model::model_Pk () {

  // multiplicative factor
  const double fact_ng = 1 / ng();
  const double fact_ng2 = fact_ng * fact_ng;

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_1halo( ii, fact_ng2 ) + Pk_2halo( ii, fact_ng2 );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > scam::halo_model::model_Xi_1halo ( const std::vector< double > & rad ) {

  std::vector< double > Pk = scam::halo_model::model_Pk_1halo();
  scam::utl::fftlog_3Dspace fft { _kv, Pk };

  return fft.transform( rad );

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Xi_2halo ( const std::vector< double > & rad ) {

  std::vector< double > Pk = scam::halo_model::model_Pk_2halo();
  scam::utl::fftlog_3Dspace fft { _kv, Pk };

  return fft.transform( rad );

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Xi ( const std::vector< double > & rad ) {

  std::vector< double > Pk = scam::halo_model::model_Pk();
  double rM = std::exp( 0.5 * ( std::log( rad.back() ) + std::log( rad.front() ) ) );
  double kM = std::exp( 0.5 * ( std::log( _kv.back() ) + std::log( _kv.front() ) ) );
  scam::utl::fftlog_3Dspace fft { _kv, Pk, rM * kM };
  
  return fft.transform( rad );

}

//==============================================================================================
    
std::vector< double > scam::halo_model::model_Wr_1halo ( const std::vector< double > & radp ) {

  std::vector< double > Pk = scam::halo_model::model_Pk_1halo();
  scam::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}
    
std::vector< double > scam::halo_model::model_Wr_cs ( const std::vector< double > & radp ) {

  std::vector< double > Pk = scam::halo_model::model_Pk_cs();
  scam::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}
    
std::vector< double > scam::halo_model::model_Wr_ss ( const std::vector< double > & radp ) {

  std::vector< double > Pk = scam::halo_model::model_Pk_ss();
  scam::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}

//==============================================================================================
    
std::vector< double > scam::halo_model::model_Wr_2halo ( const std::vector< double > & radp ) {

  std::vector< double > Pk = scam::halo_model::model_Pk_2halo();
  scam::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}

//==============================================================================================
    
std::vector< double > scam::halo_model::model_Wr ( const std::vector< double > & radp ) {

  std::vector< double > Pk = scam::halo_model::model_Pk();
  scam::utl::fftlog_projected fft { _kv, Pk };
  return fft.transform( radp );

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Wt_1halo ( const std::vector< double > & theta ) {

  const double r_z = _cosmo->d_C( _redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return scam::halo_model::model_Wr_1halo( radp );

}

std::vector< double > scam::halo_model::model_Wt_cs ( const std::vector< double > & theta ) {

  const double r_z = _cosmo->d_C( _redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return scam::halo_model::model_Wr_cs( radp );

}

std::vector< double > scam::halo_model::model_Wt_ss ( const std::vector< double > & theta ) {

  const double r_z = _cosmo->d_C( _redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return scam::halo_model::model_Wr_ss( radp );

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Wt_2halo ( const std::vector< double > & theta ) {

  const double r_z = _cosmo->d_C( _redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return scam::halo_model::model_Wr_2halo( radp );

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Wt ( const std::vector< double > & theta ) {

  const double r_z = _cosmo->d_C( _redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  return scam::halo_model::model_Wr( radp );

}

//==============================================================================================

std::vector< double > scam::halo_model::model_Pk_large_scale () {

  const double b_gal = bias();
  const double b2 = b_gal * b_gal;
  
  std::vector< double > Pk_gxy( _thinness );  
#pragma omp parallel for
  for ( size_t ii = 0; ii < _thinness; ++ii )
    Pk_gxy[ ii ] = b2 * _cosmo->Pk( _kv[ ii ], _redshift );
    
  return Pk_gxy;
  
}

//==============================================================================================
    
std::vector< double > scam::halo_model::model_Wt_large_scale ( const std::vector< double > & theta ) {

  const double r_z = _cosmo->d_C( _redshift );
  std::vector< double > radp ( theta.size() );
  for ( size_t ii = 0; ii < radp.size(); ++ii ) radp[ ii ] = theta[ ii ] * r_z;

  std::vector< double > Pk = scam::halo_model::model_Pk_large_scale();
    
  scam::utl::fftlog_projected fft { _kv, Pk };

  return fft.transform( radp );

}

//==============================================================================================
