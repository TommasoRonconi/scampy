#include <cosmology.h>
#include <halo_model_c_interface.h>
#include <halo_model.h>

#ifdef USE_CBL
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;
#endif

extern "C" {

  
  // ========================================================================================
  // ========================================================================================
  // ========================================================================================
  // Handler
  
  hm_handler_t create_hm_handler ( double DC,
				   double Mmin,
				   double sigma_logM,
				   double M0,
				   double M1,
				   double alpha,
				   double redshift,
				   cosmology_t cosmo ) {

    if ( cosmo )
      return new sico::halo_model_handler { DC, Mmin, sigma_logM, M0, M1, alpha, redshift,
	  * static_cast< sico::cosmology * >( cosmo ) };
    else
      return new sico::halo_model_handler { DC, Mmin, sigma_logM, M0, M1, alpha, redshift };

  }
  
  // ========================================================================================

  unsigned int hm_thinness ( hm_handler_t hm_h ) {

    return static_cast< sico::halo_model_handler * >( hm_h )->thinness;

  }
  
  // ========================================================================================

  void free_hm_handler ( hm_handler_t hm_h ) {

    delete static_cast< sico::halo_model_handler * >( hm_h );

    return;

  }

  
  // ========================================================================================
  // ========================================================================================
  // ========================================================================================
  // Class Halo Model

  halo_model_t create_halo_model ( hm_handler_t hm_h ) {

    // std::cout << "Building halo_model..." << std::endl;
    return new sico::halo_model { * static_cast< sico::halo_model_handler * >( hm_h ) };

  }
  
  // ========================================================================================

  void free_halo_model ( halo_model_t hm ) {

    // std::cout << "Deleting halo_model..." << std::endl;
    delete static_cast< sico::halo_model * >( hm );

    return;

  }
  
  // ========================================================================================

  void set_parameters_hm ( double DC,
			   double Mmin,
			   double sigma_logM,
			   double M0,
			   double M1,
			   double alpha,
			   halo_model_t hm ) {
    
    static_cast< sico::halo_model * >( hm )->set_parameters( DC, Mmin, sigma_logM, M0, M1, alpha );
    return;
    
  }
  
  // ========================================================================================

  double Ncen_hm ( double Mh, halo_model_t hm ) {

    return static_cast< sico::halo_model * >( hm )->Ncen( Mh );

  }
  
  // ========================================================================================

  double Nsat_hm ( double Mh, halo_model_t hm ) {

    return static_cast< sico::halo_model * >( hm )->Nsat( Mh );

  }
  
  // ========================================================================================

  double ng_hm ( halo_model_t hm ) {

    return static_cast< sico::halo_model * >( hm )->ng();

  }
  
  // ========================================================================================

  double bias_hm ( halo_model_t hm ) {

    return static_cast< sico::halo_model * >( hm )->bias();

  }
  
  // ========================================================================================

  double Mhalo_hm ( halo_model_t hm ) {

    return static_cast< sico::halo_model * >( hm )->Mhalo();

  }
  
  // ========================================================================================

  double dngdM_hm ( double Mh, halo_model_t hm ) {

    return static_cast< sico::halo_model * >( hm )->dngdM( Mh );

  }
  
  // ========================================================================================

  void model_Pk_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< sico::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::halo_model * >( hm )->model_Pk();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_1halo_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< sico::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::halo_model * >( hm )->model_Pk_1halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }

  void model_Pk_cs_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< sico::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::halo_model * >( hm )->model_Pk_cs();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }

  void model_Pk_ss_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< sico::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::halo_model * >( hm )->model_Pk_ss();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Pk_2halo_hm ( double * kv, double * Pk, halo_model_t hm ) {

    std::vector< double > kv_v = static_cast< sico::halo_model * >( hm )->get_kv();
    std::vector< double > Pk_v = static_cast< sico::halo_model * >( hm )->model_Pk_2halo();

    for ( size_t ii = 0; ii < Pk_v.size(); ++ii ) {
      kv[ ii ] = kv_v[ ii ];
      Pk[ ii ] = Pk_v[ ii ];
    }

    return;

  }
  
  // ========================================================================================

  void model_Xi_hm ( double * rr, double * Xi, unsigned int size, halo_model_t hm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< sico::halo_model * >( hm )->model_Xi( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_1halo_hm ( double * rr, double * Xi, unsigned int size, halo_model_t hm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< sico::halo_model * >( hm )->model_Xi_1halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Xi_2halo_hm ( double * rr, double * Xi, unsigned int size, halo_model_t hm ) {

    std::vector< double > rr_v { rr, rr + size };
    std::vector< double > Xi_v = static_cast< sico::halo_model * >( hm )->model_Xi_2halo( rr_v );

    for ( size_t ii = 0; ii < size; ++ii ) Xi[ ii ] = Xi_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wr_hm ( double * rp, double * Wr, unsigned int size, halo_model_t hm ) {

    std::vector< double > rp_v { rp, rp + size };
    std::vector< double > Wr_v = static_cast< sico::halo_model * >( hm )->model_Wr( rp_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wr[ ii ] = Wr_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wr_1halo_hm ( double * rp, double * Wr, unsigned int size, halo_model_t hm ) {

    std::vector< double > rp_v { rp, rp + size };
    std::vector< double > Wr_v = static_cast< sico::halo_model * >( hm )->model_Wr_1halo( rp_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wr[ ii ] = Wr_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wr_2halo_hm ( double * rp, double * Wr, unsigned int size, halo_model_t hm ) {

    std::vector< double > rp_v { rp, rp + size };
    std::vector< double > Wr_v = static_cast< sico::halo_model * >( hm )->model_Wr_2halo( rp_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wr[ ii ] = Wr_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< sico::halo_model * >( hm )->model_Wt( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_1halo_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< sico::halo_model * >( hm )->model_Wt_1halo( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }

  void model_Wt_cs_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< sico::halo_model * >( hm )->model_Wt_cs( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }

  void model_Wt_ss_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< sico::halo_model * >( hm )->model_Wt_ss( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_2halo_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< sico::halo_model * >( hm )->model_Wt_2halo( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  void model_Wt_large_scale_hm ( double * tt, double * Wt, unsigned int size, halo_model_t hm ) {

    std::vector< double > tt_v { tt, tt + size };
    std::vector< double > Wt_v = static_cast< sico::halo_model * >( hm )->model_Wt_large_scale( tt_v );

    for ( size_t ii = 0; ii < size; ++ii ) Wt[ ii ] = Wt_v[ ii ];

    return;
    
  }
  
  // ========================================================================================

  double * get_kv_hm ( halo_model_t hm ) {

    std::vector< double > kv = static_cast< sico::halo_model * >( hm )->get_kv();
    
    return kv.data();

  }
  
  // ========================================================================================
  
} // end extern "C"
