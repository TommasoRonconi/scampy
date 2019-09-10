#ifndef __HALO_MODEL__
#define __HALO_MODEL__

/// gsl includes
#include <gsl/gsl_sf_erf.h>

/// internal includes
#include <utilities.h>
#include <cosmology.h>
#include <interpolation.h>
#include <fftlog_handler.h>
#include <abel_transform.h>

namespace sico {

  struct halo_model_handler {

    double DC;
    double M_min;
    double sigma_logM;
    double M0;
    double M1;
    double alpha;
    double redshift = 1.e-7;

    struct {
      double inf = 1.e+5;
      double sup = 1.e+17;
    } mass_integ_lim {};

    struct {
      double inf = 1.e-3;
      double sup = 1.e+3;
    } wavk_integ_lim {};

    size_t thinness = 200;

    cosmology cosmo {};

    halo_model_handler () = default;

    halo_model_handler ( const double & DC,
			 const double & Mmin,
			 const double & siglogM,
			 const double & M0,
			 const double & M1,
			 const double & alpha,
			 const double & zz = 1.e-7,
			 const cosmology & cosmo = {} ) :
      DC { DC }, M_min { Mmin }, sigma_logM { siglogM }, M0 { M0 }, M1 { M1 }, alpha { alpha },
      redshift { zz }, cosmo { cosmo } {}

    halo_model_handler ( const halo_model_handler & handler ) = default;

    ~halo_model_handler() = default;
    
  }; //endstruct halo_model_handler

  class halo_model {

  private:

    halo_model_handler _handler {};

    /**
     *  @name Interpolated Functions
     *  
     *  @{
     */

    using interp_func = class sico::utl::interpolator< sico::utl::gsl_log_interp >;

    interp_func dndM_f {}, hbias_f {}, Ncen_f {}, Nsat_f {}, DC_f {}, const2_f {};

    std::vector< interp_func > density_profile_FS;

    /// @} End of Interpolated Functions
    
    std::vector< double > _kv, _Mv;
    
    double PP ( const double AA, const double Amin, const double sigma_logA );

    double Pk_cc_integrand ( const double Mh, const double kk );

    double Pk_cs_integrand ( const double Mh, const double kk );

    double Pk_ss_integrand ( const double Mh, const double kk );

  public:

    double Ncen ( const double Mhalo );
    
    double Nsat ( const double Mhalo );

    halo_model () = default;

    halo_model ( halo_model_handler handler );

    ~halo_model () = default;

    std::vector< double > get_kv () { return _kv; }

    void set_parameters ( const double DC,
			  const double Mmin,
			  const double sigma_logM,
			  const double M0,
			  const double M1,
			  const double alpha );
    // void set_parameters ( const double A_min,
    // 			  const double sigma_logA,
    // 			  const double A_sat,
    // 			  const double alpha_sat );

    double ng ();
    
    double bias ();

    double Mhalo ();
    
    double dngdM ( const double Mhalo );

    double Pk_1halo ( const size_t ii, const double fact_ng2 );
    double Pk_cs ( const size_t ii, const double fact_ng2 );
    double Pk_ss ( const size_t ii, const double fact_ng2 );
    double Pk_2halo ( const size_t ii, const double fact_ng2 );

    std::vector< double > model_Pk_1halo ();
    std::vector< double > model_Pk_cs ();
    std::vector< double > model_Pk_ss ();
    std::vector< double > model_Pk_2halo ();
    std::vector< double > model_Pk ();
    
    std::vector< double > model_Xi_1halo ( const std::vector< double > & rad );
    std::vector< double > model_Xi_2halo ( const std::vector< double > & rad );
    std::vector< double > model_Xi ( const std::vector< double > & rad );
    
    std::vector< double > model_Wr_1halo ( const std::vector< double > & radp );
    std::vector< double > model_Wr_cs ( const std::vector< double > & radp );
    std::vector< double > model_Wr_ss ( const std::vector< double > & radp );
    std::vector< double > model_Wr_2halo ( const std::vector< double > & radp );
    std::vector< double > model_Wr ( const std::vector< double > & radp );
    
    std::vector< double > model_Wt_1halo ( const std::vector< double > & theta );
    std::vector< double > model_Wt_cs ( const std::vector< double > & theta );
    std::vector< double > model_Wt_ss ( const std::vector< double > & theta );
    std::vector< double > model_Wt_2halo ( const std::vector< double > & theta );
    std::vector< double > model_Wt ( const std::vector< double > & theta );

    
    std::vector< double > model_Pk_large_scale ();
    std::vector< double > model_Wt_large_scale ( const std::vector< double > & theta );    

  }; //endclass halo_model

} //endnamespace sico

#endif //__HALO_MODEL__
