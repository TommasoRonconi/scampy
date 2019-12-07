/**
 *  @file halo_model/include/halo_model.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __HALO_MODEL__
#define __HALO_MODEL__

/// gsl includes
#include <gsl/gsl_sf_erf.h>

#include <memory>

/// internal includes
#include <utilities.h>
#include <cosmology.h>
#include <interpolation.h>
#include <fftlog_handler.h>
#include <abel_transform.h>
#include <occupation_p.h>

namespace scam {

  class halo_model {

  private:

    std::shared_ptr< occupation_p > _handler {};
    std::shared_ptr< cosmology > _cosmo {};
    
    double _redshift = 1.e-7;

    struct {
      double inf = 1.e+5;
      double sup = 1.e+17;
    } _mass_integ_lim {};

    struct {
      double inf = 1.e-3;
      double sup = 1.e+3;
    } _wavk_integ_lim {};

    size_t _thinness = 200;
    
    /**
     *  @name Interpolated Functions
     *  
     *  @{
     */

    using interp_func = class scam::utl::interpolator< scam::utl::gsl_log_interp >;

    interp_func dndM_f {}, hbias_f {}, Ncen_f {}, Nsat_f {}, const2_f {};

    std::vector< interp_func > density_profile_FS;

    /// @} End of Interpolated Functions
    
    std::vector< double > _kv, _Mv;

    double Pk_cc_integrand ( const double Mh, const double kk );

    double Pk_cs_integrand ( const double Mh, const double kk );

    double Pk_ss_integrand ( const double Mh, const double kk );

  public:

    double Ncen ( const double Mhalo ) { return _handler->Ncen( Mhalo ); }
    
    double Nsat ( const double Mhalo ) { return _handler->Nsat( Mhalo ); }

    halo_model () = default;

    halo_model ( const std::shared_ptr< occupation_p > & ocp,
		 const std::shared_ptr< cosmology > & cosmo,
		 const double redshift = 1.e-7,
		 const size_t thinness = 50 );

    ~halo_model () = default;

    size_t get_thinness () { return _thinness; }

    std::vector< double > get_kv () { return _kv; }

    void set_parameters ( const std::shared_ptr< occupation_p > & ocp );

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

} //endnamespace scam

#endif //__HALO_MODEL__
