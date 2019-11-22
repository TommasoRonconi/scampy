#ifndef __INTERPOLATOR__
#define __INTERPOLATOR__

/// gsl includes
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

/// internal includes
// #include <error_handling.h>
#include <utilities.h>
#include <interpolation_interface.h>
#include <gsl_interpolation_interface.h>

namespace scam {

  namespace utl {

    template< class T = gsl_log_interp >
    class interpolator {

    private:

      T _interface;

    public:

      interpolator () = default;

      interpolator ( const T & interface )
	: _interface{ interface } {}

      interpolator ( std::function< double ( double ) > func,
		     const double x_min, const double x_max,
		     const size_t thinness ) noexcept
	: _interface{ func, x_min, x_max, thinness } {}

      interpolator ( const std::vector< double > & xv,
      		     const std::vector< double > & fv )
      	: _interface{ xv, fv } {}

      ~interpolator () = default;

      double operator() ( const double xx ) {
  
	return _interface.eval( xx );
  
      }

      double integrate ( const double aa, const double bb ) {

	return _interface.integrate( aa, bb );

      }

      size_t get_thinness () { return _interface.get_thinness(); }
      
      double get_xmin () { return _interface.get_xmin(); }
      
      double get_xmax () { return _interface.get_xmax(); }
      
      std::vector< double > get_xv () { return _interface.get_xv(); }

      std::vector< double > get_fv () { return _interface.get_fv(); }

      size_t size () { return _interface.size(); }

      interpolator & operator+= ( const interpolator & rhs ) {

	_interface += rhs._interface;
	
	return * this;

      }
      
      friend interpolator operator+ ( interpolator lhs,
				      const interpolator & rhs ) {

	lhs += rhs;
	return lhs;

      }

      interpolator & operator*= ( const interpolator & rhs ) {

	_interface *= rhs._interface;
	
	return * this;

      }
      
      friend interpolator operator* ( interpolator lhs,
				      const interpolator & rhs ) {

	lhs *= rhs;
	return lhs;

      }

    }; //endclass interpolator

  } //endnamespace utl

} //endnamespace scam

#endif //__INTERPOLATOR__
