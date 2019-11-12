#include <error_handling.h>
#include <cross_halo_model.h>
#include <functional>

//==============================================================================================

sico::cross_halo_model::cross_halo_model ( const std::shared_ptr< sico::occupation_p > & ocp1,
					   const std::shared_ptr< sico::occupation_p > & ocp2,
					   const std::shared_ptr< sico::cosmology > & cosmo,
					   const double redshift,
					   const size_t thinness ) : _pop1{ ocp1 },
								     _pop2{ ocp2 },
								     _cosmo{ cosmo },
								     _redshift{ redshift },
								     _thinness{ thinness } {
								       
  std::cout << "Computing interpolated functions ...";

  // cosmological
  auto f_dndM = [ & ] ( double Mh ) { return _cosmo->dndM( Mh, _redshift ); };
  dndM_f = sico::cross_halo_model::interp_func { f_dndM,
						 _mass_integ_lim.inf,
						 _mass_integ_lim.sup,
						 _thinness };
  auto f_hbias = [ & ] ( double Mh ) { return _cosmo->hbias( Mh, _redshift ); };
  hbias_f = sico::cross_halo_model::interp_func { f_hbias,
						  _mass_integ_lim.inf,
						  _mass_integ_lim.sup,
						  _thinness };

  // Population 1
  auto f_Ng1 = [ & ] ( double Mh ) { return _pop1->Ncen( Mh ) + _pop1->Nsat( Mh ); };
  Ng1_f = sico::cross_halo_model::interp_func { f_Ng1,
						_mass_integ_lim.inf,
						_mass_integ_lim.sup,
						_thinness };
  // auto f_Ncen1 = [ & ] ( double Mh ) { return _pop1->Ncen( Mh ); };
  // Ncen1_f = sico::cross_halo_model::interp_func { f_Ncen1,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  // auto f_Nsat1 = [ & ] ( double Mh ) { return _pop1->Nsat( Mh ); };
  // Nsat1_f = sico::cross_halo_model::interp_func { f_Nsat1,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };

  // Population 2
  auto f_Ng2 = [ & ] ( double Mh ) { return _pop2->Ncen( Mh ) + _pop2->Nsat( Mh ); };
  Ng2_f = sico::cross_halo_model::interp_func { f_Ng2,
						_mass_integ_lim.inf,
						_mass_integ_lim.sup,
						_thinness };
  // auto f_Ncen2 = [ & ] ( double Mh ) { return _pop2->Ncen( Mh ); };
  // Ncen2_f = sico::cross_halo_model::interp_func { f_Ncen2,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  // auto f_Nsat2 = [ & ] ( double Mh ) { return _pop2->Nsat( Mh ); };
  // Nsat2_f = sico::cross_halo_model::interp_func { f_Nsat2,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  
  _Mv = Ng1_f.get_xv();

  // computing density profile vector
  _kv = sico::utl::log_vector( _thinness,
  			       _wavk_integ_lim.inf,
  			       _wavk_integ_lim.sup );

  density_profile_FS = std::vector< utl::interpolator< utl::gsl_log_interp > > { _thinness };

#pragma omp parallel for
  for ( size_t _k = 0; _k < _thinness; ++_k ) {
    
    auto f_denFS = [ & ] ( double Mh ) {
      return _cosmo->density_profile_FS( _kv[_k], Mh, 0. );
    };
    
    density_profile_FS[ _k ] =
      sico::cross_halo_model::interp_func {
      f_denFS,
      _mass_integ_lim.inf,
      _mass_integ_lim.sup,
      _thinness
    };    
    
  }
    
  std::cout << " done." << std::endl;
  
}

//==============================================================================================

void sico::cross_halo_model::set_parameters_pop1 ( const std::shared_ptr< sico::occupation_p > & ocp1 ) {

  _pop1 = ocp1;

  auto f_Ng1 = [ & ] ( double Mh ) { return _pop1->Ncen( Mh ) + _pop1->Nsat( Mh ); };
  Ng1_f = sico::cross_halo_model::interp_func { f_Ng1,
						_mass_integ_lim.inf,
						_mass_integ_lim.sup,
						_thinness };
  
  // auto f_Ncen = [ & ] ( double Mh ) { return _pop1->Ncen( Mh ); };
  // Ncen1_f = sico::cross_halo_model::interp_func { f_Ncen,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  // auto f_Nsat = [ & ] ( double Mh ) { return _pop1->Nsat( Mh ); };
  // Nsat1_f = sico::cross_halo_model::interp_func { f_Nsat,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  
  return;
  
}

//==============================================================================================

void sico::cross_halo_model::set_parameters_pop2 ( const std::shared_ptr< sico::occupation_p > & ocp2 ) {

  _pop2 = ocp2;
  
  auto f_Ng2 = [ & ] ( double Mh ) { return _pop2->Ncen( Mh ) + _pop2->Nsat( Mh ); };
  Ng2_f = sico::cross_halo_model::interp_func { f_Ng2,
						_mass_integ_lim.inf,
						_mass_integ_lim.sup,
						_thinness };
  
  // auto f_Ncen = [ & ] ( double Mh ) { return _pop2->Ncen( Mh ); };
  // Ncen2_f = sico::cross_halo_model::interp_func { f_Ncen,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  // auto f_Nsat = [ & ] ( double Mh ) { return _pop2->Nsat( Mh ); };
  // Nsat2_f = sico::cross_halo_model::interp_func { f_Nsat,
  // 						  _mass_integ_lim.inf,
  // 						  _mass_integ_lim.sup,
  // 						  _thinness };
  
  return;
  
}

//==============================================================================================

double sico::cross_halo_model::ng1 () {
  
  auto integrand = Ng1_f * dndM_f;
    
  // integration limits of ng are from 0. up to some Mmax
  // substitute 0 with m_Mh_min? ----------------------v ?!
  return integrand.integrate( _mass_integ_lim.inf, _mass_integ_lim.sup );

}

//==============================================================================================

double sico::cross_halo_model::ng2 () {
  
  auto integrand = Ng2_f * dndM_f;
    
  // integration limits of ng are from 0. up to some Mmax
  // substitute 0 with m_Mh_min? ----------------------v ?!
  return integrand.integrate( _mass_integ_lim.inf, _mass_integ_lim.sup );

}

// ============================================================================================

double sico::cross_halo_model::Pk_1halo ( const size_t ii, const double fact ) {

  // // density profile -> it depends on cosmology
  // const double uk = _cosmo->density_profile_FS( kk, Mh );
  // const double ukp_cs = ( Nc * Ns > 1 ) ? uk * uk : uk;
  // const double ukp_ss = ( Ns * Ns > 1 ) ? uk * uk : uk;

  auto integrand =						\
    Ng1_f * Ng2_f * dndM_f					\
    * density_profile_FS[ ii ]					\
    * density_profile_FS[ ii ];

  return fact * integrand.integrate( _mass_integ_lim.inf,
				     _mass_integ_lim.sup );

}
  
// ============================================================================================

double sico::cross_halo_model::Pk_2halo ( const size_t ii, const double fact ) {

  auto integrand1 = Ng1_f * dndM_f * hbias_f * density_profile_FS[ ii ];

  double integral1 = integrand1.integrate( _mass_integ_lim.inf,
					   _mass_integ_lim.sup );

  auto integrand2 = Ng2_f * dndM_f * hbias_f * density_profile_FS[ ii ];

  double integral2 = integrand2.integrate( _mass_integ_lim.inf,
					   _mass_integ_lim.sup );
  
  return _cosmo->Pk( _kv[ ii ], _redshift ) * fact * integral1 * integral2;

}

//==============================================================================================

std::vector< double > sico::cross_halo_model::model_Pk_1halo () {

  // multiplicative factor
  const double fact = 1 / ( ng1() * ng2() );

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_1halo( ii, fact );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > sico::cross_halo_model::model_Pk_2halo () {

  // multiplicative factor
  const double fact = 1 / ( ng1() * ng2() );

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_2halo( ii, fact );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > sico::cross_halo_model::model_Pk () {

  // multiplicative factor
  const double fact = 1 / ( ng1() * ng2() );

  std::vector< double > Pk_gxy( _thinness );

#pragma omp parallel for
    for ( size_t ii = 0; ii < _thinness; ++ii )
      Pk_gxy[ ii ] = Pk_1halo( ii, fact ) + Pk_2halo( ii, fact );
    
  return Pk_gxy;
  
}

//==============================================================================================

std::vector< double > sico::cross_halo_model::model_Xi_1halo ( const std::vector< double > & rad ) {

  std::vector< double > Pk = sico::cross_halo_model::model_Pk_1halo();
  sico::utl::fftlog_3Dspace fft { _kv, Pk };

  return fft.transform( rad );

}

//==============================================================================================

std::vector< double > sico::cross_halo_model::model_Xi_2halo ( const std::vector< double > & rad ) {

  std::vector< double > Pk = sico::cross_halo_model::model_Pk_2halo();
  sico::utl::fftlog_3Dspace fft { _kv, Pk };

  return fft.transform( rad );

}

//==============================================================================================

std::vector< double > sico::cross_halo_model::model_Xi ( const std::vector< double > & rad ) {

  std::vector< double > Pk = sico::cross_halo_model::model_Pk();
  double rM = std::exp( 0.5 * ( std::log( rad.back() ) + std::log( rad.front() ) ) );
  double kM = std::exp( 0.5 * ( std::log( _kv.back() ) + std::log( _kv.front() ) ) );
  sico::utl::fftlog_3Dspace fft { _kv, Pk, rM * kM };
  
  return fft.transform( rad );

}

// ============================================================================================
