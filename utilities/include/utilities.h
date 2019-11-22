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

/// std library includes:
#include <functional> // std::function
#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cin, std::err, std::getline
#include <random>     // std::mt19937, std::generator
#include <fstream>    // std::ifstream, std::ofstream
#include <iomanip>    // std::setw, std::setprecision
#include <string>     // std::string
#include <algorithm>  // std::shuffle

/// gsl includes:
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"

/// Cuba (4.2) includes
// #include <cuba.h>

/// internal includes:
#include <error_handling.h>

namespace scam {
    
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

  /// @} End of public functions of the class

  namespace conv {

    /// converts arcseconds to radians
    inline double arcsec_to_rad ( const double angle_arcsec ) {
      return angle_arcsec * pi * 1.543e-6;
    }

    /// converts radians to arcseconds
    inline double rad_to_arcsec ( const double angle_rad ) {
      return angle_rad * ip * 6.48e+5;
    }

  } //endnamespace conv

  namespace utl {

    /**
     * @class Structure to handle functions to be integrated with gsl
     */
    struct STR_generic_func_GSL {
      
      std::function< double( double ) > f;
      double xx0;

      std::function< double( std::vector< double > ) > fmin;

      std::function< double( std::vector< double > & ) > fmin_return;
      std::vector< double > parameters_return;

    };

    /**
     * @brief This function returns the result at position <code>xx</code> of the function
     *        stored within an object of type <code>scam::utl::STR_generic_func_GSL</code>
     *        that is casted in input to a pointer to void and casted back to its original
     *        pointer type within the function itself.
     *
     * @param xx position at which to compute the value of the function
     *
     * @param params pointer to type STR_generic_func_GSL casted to pointer to void 
     *
     * @return the estimated value of the function stored within <code>params</code> struct
     */
    double gen_func ( const double xx, void * params );

    /**
     * @brief Function to handle integration with <code>gsl_integrate_qng</code> function.
     * 
     * @param func an object of std::function type taking a double as argument and
     *             returning a double as result
     * 
     * @param aa lower limit of integral
     * 	                                
     * @param bb upper limit of integral
     * 
     * @param prec required precision of the integral
     *
     * @return the double precision result of the integral
     */
    double integrate_qng ( std::function< double( double ) > func, const double aa, const double bb, const double prec = 1.e-2 );

    /**
     * @brief Function to handle integration with <code>gsl_integrate_qag</code> function.
     * 	                                             
     * @param func an object of std::function type taking a double as argument and
     * 	           returning a double as result
     *
     * @param aa lower limit of integral
     *
     * @param bb upper limit of integral
     *
     * @param prec required precision of the integral                                       
     *
     * @param limit_size number of maximum iterations
     * 
     * @param rule the integration rule
     *
     * @return the double precision result of the integral
     *
     * @warning should expand this section of doc with further infos from gsl doc page
     */
    double integrate_qag ( std::function< double( double ) > func, const double aa, const double bb, const double prec = 1.e-2, const int limit_size = 1000, const int rule = 6 );

    /** 
     * @brief Naive implementation of the trapezoidal integration method. 
     *
     * @param func an object of std::function type taking a double as argument and
     *             returning a double as result
     *
     * @param aa lower limit of integral
     *
     * @param bb upper limit of integral
     *
     * @param size thiness of the integral in the interval [ aa, bb ]
     *
     * @return the double precision result of the integral
     */
    double integrate_trap ( std::function< double( double ) > func, const double aa, const double bb, const size_t size = 100 );

    /**
     * @brief This function returns the result at position <code>xx</code> of the function
     *        from which 
     *        stored within an object of type <code>scam::utl::STR_generic_func_GSL</code>
     *        that is casted in input to a pointer to void and casted back to its original
     *        pointer type within the function itself.
     *
     * @param xx position at which to compute the value of the function
     *
     * @param params pointer to type STR_generic_func_GSL casted to pointer to void 
     *
     * @return the estimated value of the function stored within <code>params</code> struct
     *
     * @warning check that it actually does what is supposed to do, I guess it should return
     *         \f[ f(x) - f(x_0) \f]
     *         not
     *         \f[ f(x) - x_0 \f]
     */
    double gen_root ( const double xx, void * prm );
    double root_brent ( std::function< double( double ) > func, const double xx0,
			const double low_guess, const double up_guess,
			const double rel_err = 1.e-3, const double abs_err = 1.e-6 );
 
    /**
     * @brief Function to count lines in a file, it uses std::count from the STL
     *
     * @param fin reference to object of type input stream
     *
     * @return the number of lines in the file, an integer of type size_t
     *
     * @warning not templated, but it should be
     */
    size_t lines_in_file ( std::ifstream & fin );

    /**
     * @brief Function to count lines in a file, 
     *        it calls scam::utl::lines_in_file( std::ifstream & fin )
     *
     * @param input_file reference to string "/path/to/filename.ext" 
     *
     * @return the number of lines in the file, an integer of type size_t
     *
     * @warning not templated, but it should be
     */
    size_t lines_in_file ( const std::string & input_file );

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
     *
     * @warning not templated, but it should be
     */
    std::vector<double> lin_vector ( const size_t nbin, const double min, const double max );

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
     *
     * @warning not templated, but it should be
     */    
    std::vector<double> log_vector ( const size_t nbin, const double min, const double max );

    /** @brief Bitwise-based function isNegative
     *
     *  @param xx input integer
     * 
     *  @return a double precision floating point number:
     *          <ul>
     *          <li> 1 if xx \f$\leq\f$ 0;</li>
     *          <li> 0 if xx \f$>\f$ 0.</li>
     *          </ul>
     */
    inline double isNegative ( const int xx ) {
      return ( ( xx & 0x80000000 ) >> 31 | !xx );
    }

    /** @brief Bitwise-based function isPositive
     *
     *  @param xx input integer
     * 
     *  @return a double precision floating point number:
     *          <ul>
     *          <li> 0 if xx \f$\leq\f$ 0;</li>
     *          <li> 1 if xx \f$>\f$ 0.</li>
     *          </ul>
     */
    inline double isPositive ( const int xx ) {
      return !( ( xx & 0x80000000 ) >> 31 | !xx );
    }

    /**
     * @brief Templeted ( on class T ) and computationally efficient function 
     *        to compute the value of an N-order polynomial at some value of 
     *        the x-variable. This function is guaranteed to not throw exceptions
     *        ( unlike the 'poly_3' function ). It executes the least amount of
     *        operations possible.
     *
     * @param xx the value of the x-variable ( T type )
     *
     * @param coeff vector of the T-type coefficients of the polynomial, note that
     *              the order is important: the first element of the vector is the
     *              coefficient of the higher order term of the polynomial, the last
     *              element of the vector is the constant value of the polynomial
     * 
     * @result returns the value y, defined as:
     *         \f[ y = \sum_{i = 0}^N a_i * x^i \f]
     */
    template < class T >
    T poly_N ( const T & xx, const std::vector< T > & coeff ) noexcept {

      T ret_val = 0;
      for ( auto && cc : coeff )
	ret_val *= xx + cc;

      return ret_val;
      
    }

    /**
     * @brief
     *
     * @param
     * 
     * @result
     */
    template < class T >
    T poly_3 ( const T & xx, const std::vector< T > & coeff ) {

      if ( coeff.size() != 4 )
	throw scam_err::size_invalid {
	  "The coefficients vector must have size = 4. You gave me " +
	    std::to_string( coeff.size() ) + "! :("
	    };

      return coeff[ 0 ] + coeff[ 1 ] * xx + coeff[ 2 ] * xx * xx + coeff[ 3 ] * xx * xx * xx;
      
    }
    
    /**
     * @brief Templated (on typename T) function to extract a randomly selected number of object
     *        from an input vector: first it shuffles the original vector then returns the firs
     *        size elements.
     *        It uses the random number generator <code>std::mt19937</code> along with function  
     *        <code>std::shuffle</code> of the stdlib.
     *
     * @param vec the vector from which to extract a random subsample
     * 
     * @param size the size of the subsample
     * 
     * @param seed the seed used in the random engine, if a negative number is provided, it will
     *             use function <code>std::random_device</code> to generate a random seed
     *
     * @return An <code>std::vector</code> of same type of the input vector containing a randomly
     *          selected subsample of the original one. If provided size is <= of the input vector 
     *         sizeit returns the input vector itself without any action but printing a warning on
     *         <code>stderr</code>.
     */
    template < typename T >
      std::vector< T > random_selection ( const std::vector< T > vec, const size_t size, const int seed = 555666 )
      {

	if (vec.size()>size) {
	  // random engine std::generator
	  // (either you provide the seed or it will be chosen randomly)
	  std::mt19937 generator { ( seed < 0 ) ? std::random_device{}() : seed };
	  
	  // shuffle the vector:
	  std::shuffle( vec.begin(), vec.end(), generator );
	  
	  //return only first size elements
	  vec.resize( size );
	  return vec;
	}
	else {
	  std::cerr << "Warning in random_selection of utilities.h: provided vector smaller than required size!" << std::endl;
	  return vec;
	}
      }

    template < typename T > 
    T TopHat_WF ( const T kR ) {
      
      return 3. * ( std::sin( kR ) - kR * std::cos( kR ) ) / ( kR * kR * kR );
      
    }

    template < typename T > 
    T TopHat_WF_D1 ( const T kR ) 
    {
      return ( 3. * ( kR * kR - 3. ) * std::sin( kR ) + 9. * kR * std::cos( kR ) ) \
	/ ( kR * kR * kR * kR );
    }
    
  } //endnamespace utl
} //endnamespace scam

#endif //__UTILITIES__
