#include <fftlog_handler.h>

// ============================================================================

scam::utl::fftlog_handler::fftlog_handler ( const std::vector< double> & xx,
					    const std::vector< double> & fx,
					    const double kr , const int dir,
					    const double qq, const double mu
					    ) : _dir { dir }, _kr { kr }, _q { qq }, _mu { mu } {

  // std::cout << "Parameters:"
  // 	    << "\n\t- dir =\t" << _dir
  // 	    << "\n\t- kr =\t" << _kr
  // 	    << "\n\t- q =\t" << _q
  // 	    << "\n\t- mu =\t" << _mu
  // 	    << std::endl;
    
  scam::utl::interpolator< scam::utl::gsl_log_interp > interp { xx, fx };
  _size = interp.size();
  _dlnx = std::log( xx.back() / xx.front() ) / ( _size + 1 );
  
  _ci = 0.5 * ( _size + 1 );
  _lnx_med = 0.5 * std::log( xx.back() * xx.front() );
  _lnk_med = std::log( _kr ) - _lnx_med;
  _xx.resize( _size );
  _kk.resize( _size );
  _fx.resize( _size );
  _fk.resize( _size );
  for ( int ii = 0; ii < _size; ++ii ) {
    _xx[ ii ] = std::exp( _lnx_med + ( ii + 1 - _ci ) * _dlnx );
    _fx[ ii ] = interp( _xx[ ii ] );
    _kk[ ii ] = std::exp( _lnk_med + ( ii + 1 - _ci ) * _dlnx );
  }
  
  const int NMAX = 4096;
  _wsave = new double[ 2 * NMAX + 3 * ( NMAX / 2 ) + 19 ];

  scam::fhti_ ( &_size, &_mu, &_q, &_dlnx, &_kr, &_kropt, _wsave, &_ok );
  
  _ap = new double[ _size ];
  
}

// ============================================================================

scam::utl::fftlog_3Dspace::fftlog_3Dspace (  const std::vector< double> & xx,
					     const std::vector< double> & fx,
					     const double kr )
  : fftlog_handler{ xx, fx, kr, 1, 0., 0.5 }
{

  for ( int ii = 0; ii < _size; ++ii )
    _ap[ ii ] = _fact * _fx[ ii ] * _xx[ ii ] * std::sqrt( _xx[ ii ] );

}

// ============================================================================

void scam::utl::fftlog_3Dspace::transform () {

  _transform();
  for ( size_t ii = 0; ii < _kk.size(); ++ii )
    _fk[ ii ] = _ap[ ii ] / ( std::sqrt( _kk[ ii ] ) * _kk[ ii ] );
  fk = scam::utl::interpolator< scam::utl::gsl_log_interp > { _kk, _fk };
  
  return;

}

// ============================================================================

std::vector< double > scam::utl::fftlog_3Dspace::transform ( const std::vector< double > & kk ) {

  transform();
  std::vector< double > ret_vec;
  ret_vec.reserve( kk.size() );
  for ( auto && _k : kk )
    ret_vec.push_back( fk( _k ) );

  return ret_vec;

}

// ============================================================================

scam::utl::fftlog_projected::fftlog_projected ( const std::vector< double> & xx,
						const std::vector< double> & fx,
						const double kr  )
  : fftlog_handler{ xx, fx, kr, 1, 0., 0. }
{

  for ( int ii = 0; ii < _size; ++ii )
    _ap[ ii ] = _fact * _fx[ ii ] * _xx[ ii ];  

}

// ============================================================================

void scam::utl::fftlog_projected::transform () {

  _transform();
  for ( size_t ii = 0; ii < _kk.size(); ++ii )
    _fk[ ii ] = _ap[ ii ] / _kk[ ii ];
  fk = scam::utl::interpolator< scam::utl::gsl_log_interp > { _kk, _fk };
  
  return;

}

// ============================================================================

std::vector< double > scam::utl::fftlog_projected::transform ( const std::vector< double > & kk ) {

  transform();
  std::vector< double > ret_vec;
  ret_vec.reserve( kk.size() );
  for ( auto && _k : kk )
    ret_vec.push_back( fk( _k ) );

  return ret_vec;

}

// ============================================================================
