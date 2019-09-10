#ifndef __ABEL_TRANSFORM__
#define __ABEL_TRANSFORM__

#include <utilities.h>
#include <interpolation.h>

namespace sico {
  namespace utl {

    template < class T = gsl_log_interp >
    std::vector< double > abel_transform ( std::vector< double > rp_vec, interpolator< T > func ) {

      std::cout << "Calling Abel-Transform\n";
      size_t size = func.size();
      std::vector< double > xv = func.get_xv();
      std::vector< double > fv = func.get_fv();
      std::vector< double > transformed;
      transformed.reserve( size );
    
      for ( auto && rp : rp_vec ) {

	std::vector< double > integrand ( size );
	for ( size_t ii = 0; ii < size; ++ii ) {
	  if ( xv[ ii ] > rp ) 
	    integrand[ ii ] = 2 * fv[ ii ] * xv[ ii ] / std::sqrt( xv[ ii ] * xv[ ii ] - rp * rp );
	  else
	    integrand[ ii ] = 0.;
	}
	interpolator< T > integrand_f { xv, integrand };
	std::cout << "check integrand_f ( rp = " << rp << ")\n\t"
		  << integrand_f.get_xmin() << "\t" << integrand_f.get_xmax()
		  << "\n\t"
		  << rp << "\t" << 0.999999999 * integrand_f.get_xmax()
		  << std::endl;
	transformed.push_back( integrand_f.integrate( rp, 0.999999999 * integrand_f.get_xmax() ) );
	
      }

      // for ( auto && rp : rp_vec ) {

      // 	std::vector< double > x_integrand, y_integrand;
      // 	x_integrand.reserve( size );
      // 	y_integrand.reserve( size );
      // 	for ( size_t ii = 0; ii < size; ++ii )
      // 	  if ( xv[ ii ] > rp ) {
      // 	    x_integrand.push_back( xv[ ii ] );
      // 	    y_integrand.push_back( 2 * fv[ ii ] * xv[ ii ] / std::sqrt( xv[ ii ] * xv[ ii ] - rp * rp ) );
      // 	  }
      // 	if ( ( x_integrand.size() == y_integrand.size() ) && ( x_integrand.size() > 5 ) ) {
      // 	  interpolator< T > integrand_f { x_integrand, y_integrand };
      // 	  std::cout << "check integrand_f\n\t"
      // 		    << integrand_f.get_xmin() << "\t" << integrand_f.get_xmax()
      // 		    << "\n\t"
      // 		    << x_integrand.front() << "\t" << x_integrand.back()
      // 		    << std::endl;
      // 	  transformed.push_back( integrand_f.integrate( x_integrand.front(), x_integrand.back() ) );
      // 	}
      // 	else transformed.push_back( 0. );

      // }

      std::cout << "Exiting Abel-transform\n";
      return transformed;

    }

  } // end namespace utl
} // end namespace sico

#endif //__ABEL_TRANSFORM__
