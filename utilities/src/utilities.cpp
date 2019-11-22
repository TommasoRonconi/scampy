#include <utilities.h>

using namespace scam;

//==============================================================================================

size_t scam::utl::lines_in_file ( std::ifstream & fin ) {

  return std::count( std::istreambuf_iterator< char >( fin ),
		     std::istreambuf_iterator< char >(), '\n' );

}

//==============================================================================================

size_t scam::utl::lines_in_file ( const std::string & input_file ) {

  std::ifstream fin ( input_file.c_str() );
  size_t counter = utl::lines_in_file( fin );
  fin.clear(); fin.close();
  
  return counter;

}

//==============================================================================================

double scam::utl::gen_func ( const double xx, void * prm ) {
  
  utl::STR_generic_func_GSL * pp = ( utl::STR_generic_func_GSL * ) prm;
  return pp->f( xx );

}

//==============================================================================================

double scam::utl::integrate_qng ( std::function< double( double ) > func,
				  const double aa, const double bb, const double prec ) {

  // build gsl_function
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = &gen_func;
  Func.params = &params;

  // integrate:
  gsl_set_error_handler_off();

  double Int, error;
  size_t neval;
  gsl_integration_qng(&Func, aa, bb, 0., prec, &Int, &error, &neval);

  if ( error / fabs( Int ) > prec ) 
    std::cerr << std::setprecision(4)
	      << "Warning in integrate_qng: unable to reach the requested precision!\n\t"
	      << "Measured:\t" << std::setw(8) << error / fabs( Int )
	      << " (" << error << "/" << fabs( Int ) << ")\n\t"
	      << "Required:\t" << std::setw(8) << prec << "\n\t"
	      << "Evaluations:\t" << std::setw(8) << neval << std::endl;
  
  return Int;

}

//==============================================================================================

double scam::utl::integrate_qag ( std::function< double( double ) > func,
				  const double aa, const double bb,
				  const double prec, const int limit_size, const int rule ) {

  // build gsl_function
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = &gen_func;
  Func.params = &params;

  // integrate:
  gsl_set_error_handler_off();

  double Int, error;
  gsl_integration_workspace * ww = gsl_integration_workspace_alloc( limit_size );
  gsl_integration_qag( &Func, aa, bb, 0., prec, limit_size, rule, ww, &Int, &error ); 
  gsl_integration_workspace_free( ww );

  if ( error / fabs( Int ) > prec ) 
    std::cerr << std::setprecision(4)
	      << "Warning in integrate_qag: unable to reach the requested precision!\n\t"
	      << "Measured:\t" << std::setw(8) << error / fabs( Int )
	      << " (" << error << "/" << fabs( Int ) << ")\n\t"
	      << "Required:\t" << std::setw(8) << prec << std::endl;
  
  return Int;

}

//==============================================================================================

double scam::utl::integrate_trap ( std::function< double( double ) > func,
				   const double aa, const double bb, const size_t size ) {
 
  double integral = 0.;
  std::vector< double > xx = utl::log_vector( size, aa, bb );
  for ( size_t ii = 0; ii < size - 1; ++ii )
    integral += ( func( xx[ ii ] ) + func( xx[ ii + 1 ] ) ) * fabs( xx[ ii + 1 ] - xx[ ii ] );
  
  return integral;

}

//==============================================================================================

double scam::utl::gen_root ( const double xx, void * prm ) {
  
  utl::STR_generic_func_GSL * pp = ( utl::STR_generic_func_GSL * ) prm;
  return pp->f( xx ) - pp->xx0;

}

//==============================================================================================

double scam::utl::root_brent (  std::function< double( double ) > func, const double xx0,
				const double low_guess, const double up_guess,
				const double rel_err, const double abs_err ) {


  // build gsl_function
  STR_generic_func_GSL params;
  params.f = func;
  params.xx0 = xx0;

  gsl_function Func;  
  Func.function = &gen_root;
  Func.params = &params;

  // find root:
  gsl_set_error_handler_off();

  int status = 0, iter = 0, max_iter = 10000;
  const gsl_root_fsolver_type * T;
  double r;

  gsl_root_fsolver * s;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc( T );
  
  double x_lo = low_guess;
  double x_hi = up_guess;

  gsl_root_fsolver_set( s, &Func, x_lo, x_hi ); 

  do {
    ++ iter;
    status = gsl_root_fsolver_iterate( s );

    if ( (status != GSL_SUCCESS ) && ( status != GSL_CONTINUE ) )
      throw scam_err::gsl_fail {
	"Error in gsl routine gsl_root_fsolver_iterate() used in scam::utl::root_brent()."
	  };
    
    r = gsl_root_fsolver_root( s );
    x_lo = gsl_root_fsolver_x_lower( s );
    x_hi = gsl_root_fsolver_x_upper( s );
    
    status = gsl_root_test_interval( x_lo, x_hi, abs_err, rel_err );

    if ( ( status != GSL_SUCCESS ) && ( status != GSL_CONTINUE ) )
      throw scam_err::gsl_fail {
	"Error in gsl routine gsl_root_test_interval() used in scam::utl::root_brent()."
	  };
    
  } while ( status == GSL_CONTINUE && iter < max_iter );

  gsl_root_fsolver_free( s );

  if ( status != GSL_SUCCESS )
    throw scam_err::gsl_fail {
      "Error in scam::utl::root_brent(): routine was not able to find root."
	};
  
  return r;

}

//==============================================================================================

std::vector<double> scam::utl::lin_vector ( const size_t nbin, const double min, const double max ) {

  double bin = ( max - min ) / ( nbin - 1 );
  std::vector< double > vec;
  for (size_t ii = 0; ii < nbin; ++ii )
    vec.emplace_back( min + bin * ii );
  
  return vec;

}

//==============================================================================================

std::vector<double> scam::utl::log_vector ( const size_t nbin, const double min, const double max ) {

  double bin = ( log( max ) - log( min ) )/( nbin - 1 );
  std::vector< double > vec;
  for ( size_t ii = 0; ii < nbin; ++ii )
    vec.emplace_back( exp( log( min ) + bin * ii ) );
  
  return vec;

}

//==============================================================================================
