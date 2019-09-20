#include <occupation_c_interface.h>
#include <harikane16_p.h>
#include <tinker10_p.h>

extern "C" {

  
  // ========================================================================================
  // ========================================================================================
  // ========================================================================================
  // Handler

  H16_occupation_t create_H16_occupation ( double DC,
					   double Mmin,
					   double sigma_logM,
					   double M0,
					   double M1,
					   double alpha ) {

    return new sico::harikane16_p { DC, Mmin, sigma_logM, M0, M1, alpha };

  }

  

  T10_occupation_t create_T10_occupation ( double Amin,
					   double siglogA,
					   double Asat,
					   double alpsat ) {

    return new sico::tinker10_p { Amin, siglogA, Asat, alpsat };

  }
  
  // hm_handler_t create_hm_handler ( double DC,
  // 				   double Mmin,
  // 				   double sigma_logM,
  // 				   double M0,
  // 				   double M1,
  // 				   double alpha,
  // 				   double redshift,
  // 				   cosmology_t cosmo ) {

  //   if ( cosmo )
  //     return new sico::halo_model_handler { DC, Mmin, sigma_logM, M0, M1, alpha, redshift,
  // 	  * static_cast< sico::cosmology * >( cosmo ) };
  //   else
  //     return new sico::halo_model_handler { DC, Mmin, sigma_logM, M0, M1, alpha, redshift };

  // }
  
  // ========================================================================================

  double Ncen_H16_ocp ( double Mh, H16_occupation_t ocp_h16 ) {

    return static_cast< sico::harikane16_p * >( ocp_h16 )->Ncen( Mh );

  }
  
  double Nsat_H16_ocp ( double Mh, H16_occupation_t ocp_h16 ) {

    return static_cast< sico::harikane16_p * >( ocp_h16 )->Nsat( Mh );

  }

  double Ncen_T10_ocp ( double Mh, T10_occupation_t ocp_t10 ) {

    return static_cast< sico::tinker10_p * >( ocp_t10 )->Ncen( Mh );

  }
  
  double Nsat_T10_ocp ( double Mh, T10_occupation_t ocp_t10 ) {

    return static_cast< sico::tinker10_p * >( ocp_t10 )->Nsat( Mh );

  }
  
  // ========================================================================================

  // unsigned int hm_thinness ( hm_handler_t hm_h ) {

  //   return static_cast< sico::halo_model_handler * >( hm_h )->thinness;

  // }
  
  // ========================================================================================
  
  void free_H16_occupation ( H16_occupation_t ocp_h16 ) {

    delete static_cast< sico::harikane16_p * >( ocp_h16 );

    return;

  }
  
  void free_T10_occupation ( T10_occupation_t ocp_t10 ) {

    delete static_cast< sico::tinker10_p * >( ocp_t10 );

    return;

  }

  // void free_hm_handler ( hm_handler_t hm_h ) {

  //   delete static_cast< sico::halo_model_handler * >( hm_h );

  //   return;

  // }

  
  // ========================================================================================

} // endextern "C"
