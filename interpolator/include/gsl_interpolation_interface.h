#ifndef __GSL_INTERPOLATOR_INTERFACE__
#define __GSL_INTERPOLATOR_INTERFACE__

/// gsl includes
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

/// internal includes
#include <interpolation_interface.h>

namespace sico {

  namespace utl {

    class gsl_interpolator_interface : public interpolator_interface { 

      using uniquePtrAcc = std::unique_ptr< gsl_interp_accel, void( * )( gsl_interp_accel * ) >;
      
      using uniquePtrSpl = std::unique_ptr< gsl_spline, void( * )( gsl_spline * ) >;

    private:

      uniquePtrAcc _accel { gsl_interp_accel_alloc(), gsl_interp_accel_free };

      uniquePtrSpl _spline { nullptr, gsl_spline_free };

    protected:
      
      std::vector< double > _gv;
      
      virtual void _alloc ( const gsl_interp_type * interp_type ) noexcept {
	
	_spline.reset( gsl_spline_alloc( interp_type, _thinness ) );
	gsl_spline_init( _spline.get(), _xv.data(), _gv.data(), _thinness );

	return;

      }

    public:

      gsl_interpolator_interface () = default;

      gsl_interpolator_interface ( const double x_min, const double x_max,
				   const size_t thinness )
	: interpolator_interface{ x_min, x_max, thinness } {}

      /// move constructor
      gsl_interpolator_interface ( gsl_interpolator_interface && gsl_ii )
	: interpolator_interface{ gsl_ii },
	  _accel{ std::move( gsl_ii._accel ) }, _spline{ std::move( gsl_ii._spline ) },
	  _gv{ std::move( gsl_ii._gv) } {}

      /// copy constructor
      gsl_interpolator_interface ( const gsl_interpolator_interface & gsl_ii )
	: interpolator_interface{ gsl_ii } {

	_gv = gsl_ii._gv;
	_alloc( gsl_interp_cspline );

      }

      ~gsl_interpolator_interface () = default;

      /// move assignment
      gsl_interpolator_interface & operator= ( gsl_interpolator_interface && gsl_ii ) noexcept = default;
      
      /// copy assignment   
      gsl_interpolator_interface & operator= ( gsl_interpolator_interface gsl_ii ) {

	std::swap( _accel, gsl_ii._accel );
	std::swap( _spline, gsl_ii._spline );
	std::swap( _gv, gsl_ii._gv);

	gsl_ii.swap( *this );
	
      	return * this;

      }

      virtual double eval ( const double xx ) const noexcept override {

	return gsl_spline_eval( _spline.get(), xx, _accel.get() );
	
      }

      virtual double integrate ( const double aa, const double bb ) const noexcept override {

	return gsl_spline_eval_integ( _spline.get(), aa, bb , _accel.get() );
	
      }

      /// overload of operator += for same type add
      virtual gsl_interpolator_interface & operator+= ( const gsl_interpolator_interface & rhs ) {
	
      	if ( _thinness != rhs._thinness)
      	  throw sico_err::size_invalid {
      	    "Error in multiplication: right hand side has different size from left hand side!"
      	      };

      	for ( size_t ii = 0; ii < _thinness; ++ii ) {
	  _fv[ ii ] += rhs._fv[ ii ];
	  _gv[ ii ] += rhs._gv[ ii ];
	}

      	_spline.reset();
      	_alloc( gsl_interp_cspline );
	
      	return *this;
	
      }

      /// overload of operator *= for same type mult
      virtual gsl_interpolator_interface & operator*= ( const gsl_interpolator_interface & rhs ) {
	
      	if ( _thinness != rhs._thinness)
      	  throw sico_err::size_invalid {
      	    "Error in multiplication: right hand side has different size from left hand side!"
      	      };

      	for ( size_t ii = 0; ii < _thinness; ++ii ) {
	  _fv[ ii ] *= rhs._fv[ ii ];
	  _gv[ ii ] *= rhs._fv[ ii ];
	}

      	_spline.reset();
      	_alloc( gsl_interp_cspline );
	
      	return *this;
	
      }

    }; // endclass gsl_interpolator_interface

    struct gsl_lin_interp : public gsl_interpolator_interface {

      gsl_lin_interp () = default;
      
      gsl_lin_interp ( std::function< double ( double ) > func,
		       const double x_min, const double x_max,
		       const size_t thinness ) noexcept
	: gsl_interpolator_interface{ x_min, x_max, thinness } {
	
	_xv = lin_vector( thinness, x_min, x_max );
	
	for ( auto&& _x : _xv )
	  _fv.push_back( func( _x ) );
	_gv = _fv;
	
	_alloc( gsl_interp_cspline );
	
      }

      gsl_lin_interp ( const std::vector< double > & xv,
		       const std::vector< double > & fv )
	: gsl_interpolator_interface{ xv.front(), xv.back(), xv.size() } {

	_xv = xv;
	_fv = fv;
	_gv = fv;
	_alloc( gsl_interp_cspline );

      }

      ~gsl_lin_interp () = default;
      
      double eval ( const double xx ) const noexcept override {
	
	return gsl_interpolator_interface::eval( xx ); }
      
      double integrate ( const double aa, const double bb ) const noexcept override {
	
	return gsl_interpolator_interface::integrate ( aa, bb );
	
      }
      
    }; // endstruct gsl_lin_interp

    struct gsl_log_interp : public gsl_interpolator_interface {
      
      gsl_log_interp () = default;

      gsl_log_interp ( std::function< double ( double ) > func,
		       const double x_min, const double x_max,
		       const size_t thinness )
	: gsl_interpolator_interface{ x_min, x_max, thinness } {

	std::vector< double > xv = log_vector( thinness, x_min, x_max ); 
	_xv = lin_vector( thinness, std::log( x_min ), std::log( x_max ) );

	for ( auto&& _x : xv ) {
	  _fv.push_back( func( _x ) );
	  _gv.push_back( _x * _fv.back() );
	}
	
	_alloc( gsl_interp_cspline );
	
      }
      
      gsl_log_interp ( const std::vector< double > & xv,
		       const std::vector< double > & fv )
	: gsl_interpolator_interface{ xv.front(), xv.back(), xv.size() } {
	  
	  _xv = lin_vector( _thinness, std::log( xv.front() ), std::log( xv.back() ) );
	  _fv = fv;
	  
	  _gv.resize( _thinness );
	  for ( size_t ii = 0; ii < _thinness; ++ii ) _gv[ ii ] = xv[ ii ] * fv[ ii ];
	  
	  _alloc( gsl_interp_cspline );
	  
	}

      ~gsl_log_interp () = default;

      double eval ( const double xx ) const noexcept override {
  
        return gsl_interpolator_interface::eval( std::log( xx ) ) / xx;

      }

      double integrate ( const double aa, const double bb ) const noexcept override {

	return gsl_interpolator_interface::integrate( std::log( aa ), std::log( bb ) );

      }

      std::vector< double > get_xv () override {

	std::vector< double > exp_xv ( _thinness );
	for ( size_t ii = 0; ii < _thinness; ++ii )
	  exp_xv[ ii ] = std::exp( _xv[ ii ] );

	return exp_xv;

      }

    }; // endstruct gsl_log_interpolator

  } // endnamespace utl

} //endnamespace sico


#endif //__GSL_INTERPOLATOR_INTERFACE__
