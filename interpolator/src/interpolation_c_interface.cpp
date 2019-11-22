#include <interpolation_c_interface.h>
#include <interpolation.h>

extern "C" {

  // ========================================================================================
  // ========================================= LIN ==========================================
  // ========================================================================================
  
  lin_interpolator_t create_lin_interpolator ( const double * xv, const double * fv,
					       const size_t thinness ) {

    std::vector< double > xv_v { xv, xv + thinness };
    std::vector< double > fv_v { fv, fv + thinness };

    return new scam::utl::interpolator< scam::utl::gsl_lin_interp > { xv_v, fv_v };

  }

  // ========================================================================================
  
  void free_lin_interpolator ( lin_interpolator_t intrp ) {

    delete static_cast< scam::utl::interpolator< scam::utl::gsl_lin_interp > * >( intrp );
    
  }

  // ========================================================================================

  double lin_interpolator_eval ( const double xx, const lin_interpolator_t intrp ) {

    return
      ( * static_cast<
	scam::utl::interpolator< scam::utl::gsl_lin_interp > *
	>( intrp ) )( xx );

  }

  // ========================================================================================

  double lin_interpolator_integrate ( const double aa, const double bb,
				      const lin_interpolator_t intrp ) {

    return
      static_cast<
      scam::utl::interpolator< scam::utl::gsl_lin_interp > *
      >( intrp )->integrate( aa, bb );

  }

  // ========================================================================================
  // ========================================= LOG ==========================================
  // ========================================================================================

  log_interpolator_t create_log_interpolator ( const double * xv, const double * fv,
					       const size_t thinness  ) {

    std::vector< double > xv_v { xv, xv + thinness };
    std::vector< double > fv_v { fv, fv + thinness };

    return new scam::utl::interpolator< scam::utl::gsl_log_interp > { xv_v, fv_v };

  }

  // ========================================================================================
  
  void free_log_interpolator ( log_interpolator_t intrp ) {

    delete static_cast< scam::utl::interpolator< scam::utl::gsl_log_interp > * >( intrp );
    
  }

  // ========================================================================================

  double log_interpolator_eval ( const double xx, const lin_interpolator_t intrp ) {

    return
      ( * static_cast<
	scam::utl::interpolator< scam::utl::gsl_log_interp > *
	>( intrp ) )( xx );

  }

  // ========================================================================================

  double log_interpolator_integrate ( const double aa, const double bb,
				      const log_interpolator_t intrp ) {

    return
      static_cast<
      scam::utl::interpolator< scam::utl::gsl_log_interp > *
      >( intrp )->integrate( aa, bb );

  }
  
  // ========================================================================================

} // endextern "C"
