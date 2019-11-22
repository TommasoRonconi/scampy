#include <cosmology.h>
#include <halo_model.h>
#include <cross_halo_model.h>
#include <harikane16_p.h>
#include <tinker10_p.h>
#include <halo_model_c_interface.h>

extern "C" {
  
  // ========================================================================================
  // ========================================================================================
  // ========================================================================================
  // Class Halo Model

  halo_model_t create_halo_model_H16 ( H16_occupation_t ocp_h16,
				       cosmology_t cosmo,
				       const double redshift,
				       const size_t thinness ) {

    return new scam::halo_model {
      std::make_shared< scam::harikane16_p >( ( *static_cast< scam::harikane16_p * >( ocp_h16 ) ) ),
	std::make_shared< scam::cosmology >( ( *static_cast< scam::cosmology * >( cosmo ) ) ),
	redshift, thinness }; 

  }

  halo_model_t create_halo_model_T10 ( T10_occupation_t ocp_t10,
				       cosmology_t cosmo,
				       const double redshift,
				       const size_t thinness ) {

    return new scam::halo_model {
      std::make_shared< scam::tinker10_p >( ( *static_cast< scam::tinker10_p * >( ocp_t10 ) ) ),
	std::make_shared< scam::cosmology >( ( *static_cast< scam::cosmology * >( cosmo ) ) ),
	redshift, thinness };

  }
  
  // ========================================================================================

  void free_halo_model ( halo_model_t hm ) {

    delete static_cast< scam::halo_model * >( hm );

    return;

  }
  
  // ========================================================================================

  void set_parameters_hm_H16 ( double DC,
			       double Mmin,
			       double sigma_logM,
			       double M0,
			       double M1,
			       double alpha,
			       halo_model_t hm ) {

    static_cast< scam::halo_model * >( hm )->
      set_parameters( std::make_shared< scam::harikane16_p >( DC, Mmin,
							      sigma_logM,
							      M0, M1,
							      alpha  ) );
    // static_cast< scam::halo_model * >( hm )->set_parameters( new scam::harikane16_p {
    // 	DC, Mmin, sigma_logM, M0, M1, alpha } );
    
  }

  void set_parameters_hm_T10 ( double Amin,
			       double siglogA,
			       double Asat,
			       double alpsat,
			       halo_model_t hm ) {

    static_cast< scam::halo_model * >( hm )->
      set_parameters( std::make_shared< scam::tinker10_p >( Amin, siglogA,
							    Asat, alpsat ) );
    // static_cast< scam::halo_model * >( hm )->set_parameters( new scam::tinker10_p {
    //     Amin, siglogA, Asat, alpsat } );
    
  }
  
  // ========================================================================================

  size_t get_thinness_hm ( halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->get_thinness();

  }
  
  // ========================================================================================

  double Ncen_hm ( double Mh, halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->Ncen( Mh );

  }
  
  // ========================================================================================

  double Nsat_hm ( double Mh, halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->Nsat( Mh );

  }
  
  // ========================================================================================

  double ng_hm ( halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->ng();

  }
  
  // ========================================================================================

  double bias_hm ( halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->bias();

  }
  
  // ========================================================================================

  double Mhalo_hm ( halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->Mhalo();

  }
  
  // ========================================================================================

  double dngdM_hm ( double Mh, halo_model_t hm ) {

    return static_cast< scam::halo_model * >( hm )->dngdM( Mh );

  }
  
  // ========================================================================================

  void model_Pk_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< scam::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::halo_model * >( hm )->model_Pk();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_1halo_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< scam::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::halo_model * >( hm )->model_Pk_1halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }

  void model_Pk_cs_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< scam::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::halo_model * >( hm )->model_Pk_cs();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }

  void model_Pk_ss_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< scam::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::halo_model * >( hm )->model_Pk_ss();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_2halo_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< scam::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< scam::halo_model * >( hm )->model_Pk_2halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Xi_hm ( double * rr, double * Xi, unsigned int size, halo_model_t hm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< scam::halo_model * >( hm )->model_Xi( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_1halo_hm ( double * rr, double * Xi, unsigned int size, halo_model_t hm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< scam::halo_model * >( hm )->model_Xi_1halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_2halo_hm ( double * rr, double * Xi, unsigned int size, halo_model_t hm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< scam::halo_model * >( hm )->model_Xi_2halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wr_hm ( double * rp, double * Wr, unsigned int size, halo_model_t hm ) {

    std::vector< double > rp_v { rp, rp + size };
    std::vector< double > Wr_v = static_cast< scam::halo_model * >( hm )->model_Wr( rp_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wr[ ii ] = Wr_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wr_1halo_hm ( double * rp, double * Wr, unsigned int size, halo_model_t hm ) {

    std::vector< double > rp_v { rp, rp + size };
    std::vector< double > Wr_v = static_cast< scam::halo_model * >( hm )->model_Wr_1halo( rp_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wr[ ii ] = Wr_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wr_2halo_hm ( double * rp, double * Wr, unsigned int size, halo_model_t hm ) {

    std::vector< double > rp_v { rp, rp + size };
    std::vector< double > Wr_v = static_cast< scam::halo_model * >( hm )->model_Wr_2halo( rp_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wr[ ii ] = Wr_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< scam::halo_model * >( hm )->model_Wt( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_1halo_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< scam::halo_model * >( hm )->model_Wt_1halo( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }

  void model_Wt_cs_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< scam::halo_model * >( hm )->model_Wt_cs( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }

  void model_Wt_ss_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< scam::halo_model * >( hm )->model_Wt_ss( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_2halo_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< scam::halo_model * >( hm )->model_Wt_2halo( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_large_scale_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< scam::halo_model * >( hm )->model_Wt_large_scale( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  double * get_kv_hm ( halo_model_t hm ) {

    std::vector< double > kv = static_cast< scam::halo_model * >( hm )->get_kv();
    
    return kv.data();

  }
  
  // ========================================================================================
  
  // ========================================================================================
  // ========================================================================================
  // ========================================================================================
  // Class Cross Halo Model

  cross_halo_model_t create_cross_halo_model_H16 ( H16_occupation_t ocp1_h16,
						   H16_occupation_t ocp2_h16,
						   cosmology_t cosmo,
						   const double redshift,
						   const size_t thinness ) {

    return new sico::cross_halo_model {
      std::make_shared< sico::harikane16_p >( ( *static_cast< sico::harikane16_p * >( ocp1_h16 ) ) ),
	std::make_shared< sico::harikane16_p >( ( *static_cast< sico::harikane16_p * >( ocp2_h16 ) ) ),
	std::make_shared< sico::cosmology >( ( *static_cast< sico::cosmology * >( cosmo ) ) ),
	redshift, thinness }; 

  }

  cross_halo_model_t create_cross_halo_model_T10 ( T10_occupation_t ocp1_t10,
						   T10_occupation_t ocp2_t10,
						   cosmology_t cosmo,
						   const double redshift,
						   const size_t thinness ) {

    return new sico::cross_halo_model {
      std::make_shared< sico::tinker10_p >( ( *static_cast< sico::tinker10_p * >( ocp1_t10 ) ) ),
	std::make_shared< sico::tinker10_p >( ( *static_cast< sico::tinker10_p * >( ocp2_t10 ) ) ),
	std::make_shared< sico::cosmology >( ( *static_cast< sico::cosmology * >( cosmo ) ) ),
	redshift, thinness };

  }
  
  // ========================================================================================

  void free_cross_halo_model ( cross_halo_model_t chm ) {

    delete static_cast< sico::cross_halo_model * >( chm );

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

    static_cast< sico::cross_halo_model * >( chm )->
      set_parameters_pop1( std::make_shared< sico::harikane16_p >( DC, Mmin,
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

    static_cast< sico::cross_halo_model * >( chm )->
      set_parameters_pop2( std::make_shared< sico::harikane16_p >( DC, Mmin,
								   sigma_logM,
								   M0, M1,
								   alpha  ) );
    
  }

  void set_parameters_pop1_chm_T10 ( double Amin,
				     double siglogA,
				     double Asat,
				     double alpsat,
				     cross_halo_model_t chm ) {

    static_cast< sico::cross_halo_model * >( chm )->
      set_parameters_pop1( std::make_shared< sico::tinker10_p >( Amin, siglogA,
								 Asat, alpsat ) );
    
  }

  void set_parameters_pop2_chm_T10 ( double Amin,
				     double siglogA,
				     double Asat,
				     double alpsat,
				     cross_halo_model_t chm ) {

    static_cast< sico::cross_halo_model * >( chm )->
      set_parameters_pop2( std::make_shared< sico::tinker10_p >( Amin, siglogA,
								 Asat, alpsat ) );
    
  }
  
  // ========================================================================================

  size_t get_thinness_chm ( cross_halo_model_t chm ) {

    return static_cast< sico::cross_halo_model * >( chm )->get_thinness();

  }
  
  // ========================================================================================

  double * get_kv_chm ( cross_halo_model_t chm ) {

    std::vector< double > kv = static_cast< sico::cross_halo_model * >( chm )->get_kv();
    
    return kv.data();

  }
  
  // ========================================================================================

  double ng1_chm ( cross_halo_model_t chm ) {

    return static_cast< sico::cross_halo_model * >( chm )->ng1();

  }
  
  // ========================================================================================

  double ng2_chm ( cross_halo_model_t chm ) {

    return static_cast< sico::cross_halo_model * >( chm )->ng2();

  }
  
  // ========================================================================================

  void model_Pk_chm ( double * kv, double * Pk, cross_halo_model_t chm ) {

    std::vector< double > kv_v = static_cast< sico::cross_halo_model * >( chm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::cross_halo_model * >( chm )->model_Pk();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_1halo_chm ( double * kv, double * Pk, cross_halo_model_t chm ) {

    std::vector< double > kv_v = static_cast< sico::cross_halo_model * >( chm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::cross_halo_model * >( chm )->model_Pk_1halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_2halo_chm ( double * kv, double * Pk, cross_halo_model_t chm ) {

    std::vector< double > kv_v = static_cast< sico::cross_halo_model * >( chm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::cross_halo_model * >( chm )->model_Pk_2halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Xi_chm ( double * rr, double * Xi, unsigned int size, cross_halo_model_t chm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< sico::cross_halo_model * >( chm )->model_Xi( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_1halo_chm ( double * rr, double * Xi, unsigned int size, cross_halo_model_t chm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< sico::cross_halo_model * >( chm )->model_Xi_1halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_2halo_chm ( double * rr, double * Xi, unsigned int size, cross_halo_model_t chm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< sico::cross_halo_model * >( chm )->model_Xi_2halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================
  
} // end extern "C"
