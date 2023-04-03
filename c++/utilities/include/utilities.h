/**
 *  @file include/utilities.h
 *
 *  @brief Utilities header containing support functions 
 *
 *  This file declares a set of support functions useful for
 *  several purposes.
 *  It implements an interface for constants, conversion factors
 *  integration methods and random selections
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */


#ifndef __UTILITIES__
#define __UTILITIES__

// STL includes
#include <vector>
#include <cmath>
#include <string>

namespace utl::cnst {
    
  /**
   * @name some static constants
   *
   * @{
   */

  /// value of \f$\pi\f$
  static const double pi = 3.1415926535897932;

  /// value of \f$\pi^{-1}\f$
  static const double ip = 1/3.1415926535897932;

  /// value of \f$e^{ - 10 }\f$
  static const double pxe_10 = exp( -10 );

  /// value of \f$\ln( 10 )\f$
  static const double ln_10 = log( 10 );

  /// value of the speed of light \f$c\f$ in \f$[ m / s ]\f$
  static const double cc = 2.99792458e+8;

  /// value of 1 parsec in meters \f$[ m ]\f$
  static const double pc = 3.08567758130573e+16;

  /// value of 1 year in seconds \f$[ s ]\f$
  static const double yr = 3600 * 24 * 365.256363004;

  /// value of the newtonian gravitational constant \f$[ m^3 kg^{ -1 } s^{ -2 } ]\f$
  static const double Gn = 6.6738480e-11;

  /// inverse of the newtonian gravitational constant \f$[ m^{ -3 } kg s^2 ]\f$
  static const double nG = 1.49838594e+10;

  /// mass of the sun in grams
  static const double Msol = 1.98855e+33;

  /// inverse of the sun mass in grams
  static const double solM = 1./Msol;

  /// @} End of static constants

} // endnamespace utl::cnst
  
namespace utl_err {

  struct size_invalid {

    std::string message;
    size_invalid ( const std::string & s ) : message{ s } {}
      
  };

  struct out_of_bounds {

    std::string message;
    out_of_bounds ( const std::string & s ) : message{ s } {}
      
  };

  struct type_invalid {

    std::string message;
    type_invalid ( const std::string & s ) : message{ s } {}

  };
  
  struct not_sorted {
    
    std::string message;
    not_sorted( const std::string & s ) : message{ s } {}
    
  };
  
} //endnamespace utl_err

namespace utl {

  /// type-safe C++ implementation of the sign (a.k.a. signum) function
  /// from first comment:
  /// stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
  template < typename T >
  int sgn ( T val ) {
    return ( T( 0 ) < val ) - ( val < T( 0 ) );
  }

  template < typename T, typename V >
  T heaviside ( V x, V delay ) {
    return 0.5 * ( 1 + sgn( x - delay ) );
  }

  /**
   * @brief Function that returns a linearly spaced vector in the range [min,max]
   *
   * @param nbin number of sub-intervals in which to divide the interval
   * 
   * @param min minimum of the interval
   *
   * @param max maximum of the interval
   *
   * @return a vector with nbin linearly spaced numbers between min and max
   */
  template < typename T >
  std::vector< T > lin_vector ( const size_t nbin, const T min, const T max ) {
    
    T bin = ( max - min ) / ( nbin - 1 );
    std::vector< T > vec ( nbin );
    for ( std::size_t ii = 0; ii < nbin; ++ii )
      vec[ ii ] = min + bin * ii;
  
    return vec;

  }

  /**
   * @brief Function that returns a logarithmically spaced vector in the range [min, max]
   *	                                                                     
   * @param nbin number of sub-intervals in which to divide the interval
   *	                                                                     
   * @param min minimum of the interval
   *	                                                                     
   * @param max maximum of the interval
   *	                                                                     
   * @return a vector with nbin logarithmically spaced numbers between min and max<br>
   *         <b>Note:</b> the interval [ln(min), ln(max)] is divided in nbin linearly
   *         spaced numbers
   */
  template < typename T >
  std::vector< T > log_vector ( const size_t nbin, const T min, const T max ) {

    T bin = ( std::log( max ) - std::log( min ) )/( nbin - 1 );
    std::vector< T > vec ( nbin );
    for ( std::size_t ii = 0; ii < nbin; ++ii )
      vec[ ii ] = std::exp( std::log( min ) + bin * ii );
  
    return vec;

  }

  /**
   * @brief Naive (a.k.a. brute force) implementation of search within sorted vector of the 
   *        lower index for interpolation/extrapolation between two values
   *
   * @param val value to search for
   * @param vec sorted vector
   * @param start_indx starting point in vector (default = 0, first element of vector)
   * 
   * @result std::size_t index of vector vec.
   *         if val is smaller than the first element, returns the first index (0)
   *         if val is larger than the last element, returns the index before 
   *         the last one ( vec.size()-2 )
   * @warning the algorithm assumes:
   *          - vec.size >= 2,
   *          - vec is sorted in ascending order
   *          - start_indx <= vec.size() 
   *          no exception nor error is raised, 
   *          if the above is not respected, behaviour is undefined
   */  
  inline std::size_t
  find_low ( const double val,
	     const std::vector< double > & vec,
	     const std::size_t start_indx = 0 ) noexcept {
    
    std::size_t indx = start_indx;
    while ( indx < vec.size() - 1 ) {
      if ( val < vec[ indx + 1 ] )
	return indx;
      ++indx;
    }
    return indx - 1;

  }

  /**
   * @brief Naive (a.k.a. brute force) implementation of search within sorted C-style array of the 
   *        lower index for interpolation/extrapolation between two values
   *
   * @param val value to search for
   * @param arr sorted C-style array
   * @param asize size of arr array
   * @param start_indx starting point in array (default = 0, first element of arrat)
   * 
   * @result std::size_t index of array arr.
   *         if val is smaller than the first element, returns the first index (0)
   *         if val is larger than the last element, returns the index before 
   *         the last one ( asize-2 )
   * @warning the algorithm assumes:
   *          - asize >= 2,
   *          - arr is sorted in ascending order
   *          - start_indx <= asize
   *          no exception nor error is raised, 
   *          if the above is not respected, behaviour is undefined
   */  
  inline std::size_t
  find_low ( const double val,
	     const  double * const & arr,
	     const std::size_t asize,
	     const std::size_t start_indx = 0 ) noexcept {
    
    std::size_t indx = start_indx;
    while ( indx < asize - 1 ) {
      if ( val < arr[ indx + 1 ] )
	return indx;
      ++indx;
    }
    return indx - 1;

  }
  
  inline double
  line_from_2points ( const double xx,
		      const double x1, const double y1,
		      const double x2, const double y2 ) noexcept
  { return ( y2 - y1 ) / ( x2 - x1 ) * ( xx - x1 ) + y1; }

  /**
   * @brief Integrate on the whole x-domain with 
   *        the trapezoid rule
   *
   * @param xx x-domain grid
   * @param fx f(x) values
   * 
   * @result the integral of f(x) along the whole x-domain
   */
  inline double integrate_trap ( const std::vector< double > xx,
				 const std::vector< double > fx ) noexcept {

    double integral = 0.;
    for ( std::size_t ii = 1; ii < xx.size(); ++ii ) 
      integral += ( fx[ ii ] + fx[ ii-1 ] ) * ( xx[ ii ] - xx[ ii-1 ] );
    return 0.5 * integral;
    
  }

} //endnamespace utl

#endif //__UTILITIES__
