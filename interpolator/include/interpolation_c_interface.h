#ifndef __INTERPOLATOR_C_INTERFACE__
#define __INTERPOLATOR_C_INTERFACE__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

  // ========================================================================================
  // ========================================= LIN ==========================================
  // ========================================================================================

  typedef void * lin_interpolator_t;
  
  lin_interpolator_t create_lin_interpolator ( const double * xv, const double * fv,
					       const size_t thinness );
  
  void free_lin_interpolator ( lin_interpolator_t intrp );

  double lin_interpolator_eval ( const double xx, const lin_interpolator_t intrp );

  double lin_interpolator_integrate ( const double aa, const double bb,
				      const lin_interpolator_t intrp );

  // ========================================================================================
  // ========================================= LOG ==========================================
  // ========================================================================================

  typedef void * log_interpolator_t; 

  log_interpolator_t create_log_interpolator ( const double * xv, const double * fv,
					       const size_t thinness  );
  
  void free_log_interpolator ( log_interpolator_t intrp );

  double log_interpolator_eval ( const double xx, const log_interpolator_t intrp );

  double log_interpolator_integrate ( const double aa, const double bb,
				      const log_interpolator_t intrp );

  // ========================================================================================
  
#ifdef __cplusplus
} // endextern "C"
#endif  

#endif //__INTERPOLATOR_C_INTERFACE__
