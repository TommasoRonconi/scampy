#ifndef __CROSS_HALO_MODEL__
#define __CROSS_HALO_MODEL__

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

  class cross_halo_model {

  private:

    std::shared_ptr< occupation_p > _pop1 {}, _pop2 {};
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

    interp_func dndM_f {}, hbias_f {};
    interp_func Ng1_f {}, Ng2_f{};

    std::vector< interp_func > density_profile_FS;

    /// @} End of Interpolated Functions
    
    std::vector< double > _kv, _Mv;

  public:

    cross_halo_model () = default;

    cross_halo_model ( const std::shared_ptr< occupation_p > & ocp1,
		       const std::shared_ptr< occupation_p > & ocp2,
		       const std::shared_ptr< cosmology > & cosmo,
		       const double redshift = 1.e-7,
		       const size_t thinness = 50 );

    ~cross_halo_model () = default;

    size_t get_thinness () { return _thinness; }

    std::vector< double > get_kv () { return _kv; }

    void set_parameters_pop1 ( const std::shared_ptr< occupation_p > & ocp1 );
    void set_parameters_pop2 ( const std::shared_ptr< occupation_p > & ocp2 );

    double ng1 ();
    double ng2 ();

    double Pk_1halo ( const size_t ii, const double fact );
    double Pk_2halo ( const size_t ii, const double fact );

    std::vector< double > model_Pk_1halo ();
    std::vector< double > model_Pk_2halo ();
    std::vector< double > model_Pk ();
    
    std::vector< double > model_Xi_1halo ( const std::vector< double > & rad );
    std::vector< double > model_Xi_2halo ( const std::vector< double > & rad );
    std::vector< double > model_Xi ( const std::vector< double > & rad );

  }; //endclass cross_halo_model

} //endnamespace scam

#endif //__CROSS_HALO_MODEL__
