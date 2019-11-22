#include <cosmology.h>
#include <cross_halo_model.h>
#include <harikane16_p.h>
#include <tinker10_p.h>
#include <cross_halo_model_c_interface.h>

extern "C" {
  
  // ========================================================================================
  // ========================================================================================
  // ========================================================================================
  // Class Cross Halo Model

  cross_halo_model_t create_cross_halo_model_H16 ( H16_occupation_t ocp1_h16,
						   H16_occupation_t ocp2_h16,
						   cosmology_t cosmo,
						   const double redshift,
						   const size_t thinness ) {

    return new scam::cross_halo_model {
      std::make_shared< scam::harikane16_p >( ( *static_cast< scam::harikane16_p * >( ocp1_h16 ) ) ),
	std::make_shared< scam::harikane16_p >( ( *static_cast< scam::harikane16_p * >( ocp2_h16 ) ) ),
	std::make_shared< scam::cosmology >( ( *static_cast< scam::cosmology * >( cosmo ) ) ),
	redshift, thinness }; 

  }

  cross_halo_model_t create_cross_halo_model_T10 ( T10_occupation_t ocp1_t10,
						   T10_occupation_t ocp2_t10,
						   cosmology_t cosmo,
						   const double redshift,
						   const size_t thinness ) {

    return new scam::cross_halo_model {
      std::make_shared< scam::tinker10_p >( ( *static_cast< scam::tinker10_p * >( ocp1_t10 ) ) ),
	std::make_shared< scam::tinker10_p >( ( *static_cast< scam::tinker10_p * >( ocp2_t10 ) ) ),
	std::make_shared< scam::cosmology >( ( *static_cast< scam::cosmology * >( cosmo ) ) ),
	redshift, thinness };

  }
  
  // ========================================================================================

  void free_cross_halo_model ( cross_halo_model_t chm ) {

    delete static_cast< scam::cross_halo_model * >( chm );

    return;

  }
  
  // ========================================================================================

  void set_parameters_pop1_chm_H16 ( double DC,
				     double Mmin,
				     double sigma_logM,
				     double M0,
				     double M1,
				     double alpha,
				     cross_halo_model_t chm ) {

    static_cast< scam::cross_halo_model * >( chm )->
      set_parameters_pop1( std::make_shared< scam::harikane16_p >( DC, Mmin,
								   sigma_logM,
								   M0, M1,
								   alpha  ) );
    
  }

  void set_parameters_pop2_chm_H16 ( double DC,
				     double Mmin,
				     double sigma_logM,
				     double M0,
				     double M1,
				     double alpha,
				     cross_halo_model_t chm ) {

    static_cast< scam::cross_halo_model * >( chm )->
      set_parameters_pop2( std::make_shared< scam::harikane16_p >( DC, Mmin,
								   sigma_logM,
								   M0, M1,
								   alpha  ) );
    
  }

  void set_parameters_pop1_chm_T10 ( double Amin,
				     double siglogA,
				     double Asat,
				     double alpsat,
				     cross_halo_model_t chm ) {

    static_cast< scam::cross_halo_model * >( chm )->
      set_parameters_pop1( std::make_shared< scam::tinker10_p >( Amin, siglogA,
								 Asat, alpsat ) );
    
  }

  void set_parameters_pop2_chm_T10 ( double Amin,
				     double siglogA,
				     double Asat,
				     double alpsat,
				     cross_halo_model_t chm ) {

    static_cast< scam::cross_halo_model * >( chm )->
      set_parameters_pop2( std::make_shared< scam::tinker10_p >( Amin, siglogA,
								 Asat, alpsat ) );
    
  }
  
  // ========================================================================================

  size_t get_thinness_chm ( cross_halo_model_t chm ) {

    return static_cast< scam::cross_halo_model * >( chm )->get_thinness();

  }
  
  // ========================================================================================

  double * get_kv_chm ( cross_halo_model_t chm ) {

    std::vector< double > kv = static_cast< scam::cross_halo_model * >( chm )->get_kv();
    
    return kv.data();

  }
  
  // ========================================================================================

  double ng1_chm ( cross_halo_model_t chm ) {

    return static_cast< scam::cross_halo_model * >( chm )->ng1();

  }
  
  // ========================================================================================

  double ng2_chm ( cross_halo_model_t chm ) {

    return static_cast< scam::cross_halo_model * >( chm )->ng2();

  }
  
  // ========================================================================================

  void model_Pk_chm ( double * kv, double * Pk, cross_halo_model_t chm ) {

    std::vector< double > kv_v = static_cast< scam::cross_halo_model * >( chm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::cross_halo_model * >( chm )->model_Pk();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_1halo_chm ( double * kv, double * Pk, cross_halo_model_t chm ) {

    std::vector< double > kv_v = static_cast< scam::cross_halo_model * >( chm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::cross_halo_model * >( chm )->model_Pk_1halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_2halo_chm ( double * kv, double * Pk, cross_halo_model_t chm ) {

    std::vector< double > kv_v = static_cast< scam::cross_halo_model * >( chm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::cross_halo_model * >( chm )->model_Pk_2halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Xi_chm ( double * rr, double * Xi, unsigned int size, cross_halo_model_t chm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< scam::cross_halo_model * >( chm )->model_Xi( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_1halo_chm ( double * rr, double * Xi, unsigned int size, cross_halo_model_t chm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< scam::cross_halo_model * >( chm )->model_Xi_1halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_2halo_chm ( double * rr, double * Xi, unsigned int size, cross_halo_model_t chm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< scam::cross_halo_model * >( chm )->model_Xi_2halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================
  
} // end extern "C"
