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

double scam::utl::integrate_trap ( std::function< double( double ) > func,
				   const double aa, const double bb, const size_t size ) {
 
  double integral = 0.;
  std::vector< double > xx = utl::log_vector( size, aa, bb );
  for ( size_t ii = 0; ii < size - 1; ++ii )
    integral += ( func( xx[ ii ] ) + func( xx[ ii + 1 ] ) ) * fabs( xx[ ii + 1 ] - xx[ ii ] );
  
  return integral;

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
