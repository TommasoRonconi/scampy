/**
 *  @file cosmology/include/cosmology.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __COSMOLOGY_H__
#define __COSMOLOGY_H__

// internal includes
#include <utilities.h> // scam::ln_10 scam::utl::integrate_qng

#ifdef USE_CBL

#include <cbl_cosmology_interface.h>

/**
 *  @addtogroup scam
 *
 *  @{
 */

namespace scam { typedef struct cbl_cosmo_model cosmo_model; }

/** 
 * @} End of Doxygen Groups
 */

#else

#include <cosmology_interface.h>

/**
 *  @addtogroup scam
 *
 *  @{
 */

namespace scam { typedef struct cosmo_model cosmo_model; }

/** 
 * @} End of Doxygen Groups
 */

#endif //USE_CBL

/**
 *  @addtogroup scam
 *
 *  @{
 */

namespace scam {
  
  /**
   * @name cosmology class inheriting from scam::cosmo_model defined in file include/cosmology.h
   */
  class cosmology : public cosmo_model {

  private:

    //cosmo_p p_csmp;

    // struct phi_param_t {
    //   std::vector< double > phi_star;
    //   std::vector< double > L_star;
    //   std::vector< double > alpha;
    // };
    // static const phi_param_t Bouwens15 { { 0.47, -0.27 }, { -20.95, 0.01 }, { -1.87, -0.10 } };
    // static const phi_param_t Bouwens16 { { 0.45, -0.21 }, { -20.97, 0.17 }, { -1.91, -0.13 } };

  public:

    /**
     * @name Constructors/Destructors
     *
     * @{
     */

    /// inherit all constructors of base class
    using cosmo_model::cosmo_model;

    /**
     * @brief default constructor
     */
    cosmology () = default;

#ifndef USE_CBL
    /**
     * @brief default destructor
     */
    ~cosmology () override = default;
#endif

    /// @} End of Ctor/Dtor
    
    /**
     * @name Luminosity functions defined in src/luminosity_function.cpp
     *
     * @{
     */

    /**
     * @brief \f$\phi(L, z)\f$ luminosity function from Bouwens et al., 2015
     *
     * Parameterisation of the Schechter Luminosity function, 
     * valid in redshift range \f$4 \leq z \leq 8\f$:
     * \f[\phi(L, z) = \phi^*\; 2.5 \ln(10)\;  10^{- 0.4 ( L - L^* )( \alpha + 1 )} 
     *                 \exp\bigl(- 10^{ - 0.4 ( L - L^*) }\bigr)\f]
     * with
     * \f[ L_{UV}^* = (-20.95 \pm 0.10) + (0.01 \pm 0.06)(z - 6) \f]
     * \f[ \phi^* = ( 0.47_{-0.10}^{+0.11} ) 10^{ ( -0.27 \pm 0.05 )( z - 6 )} 10^{-3} Mpc^{-3} \f]
     * \f[ \alpha = ( -1.87 \pm 0.05 ) + ( -0.10 \pm 0.03 ) ( z - 6 ) \f]
     * 
     * @param LL UV luminosity in magnitudes
     * 
     * @param zz redshift
     *
     * @return
     */
    double phi_Bouwens15 ( const double LL, const double zz );

    /**
     * @brief \f$\phi(L, z)\f$ luminosity function from Bouwens et al., 2015
     *
     * Parameterisation of the Schechter Luminosity function, 
     * valid in redshift range \f$4 \leq z \leq 8\f$:
     * \f[\phi(L, z) = \phi^*\; 2.5 \ln(10)\;  10^{- 0.4 ( L - L^* )( \alpha + 1 )} 
     *                 \exp\bigl(- 10^{ - 0.4 ( L - L^*) }\bigr)\f]
     * with
     * \f[ L_{UV}^* = (-20.95 \pm 0.10) + (0.01 \pm 0.06)(z - 6) \f]
     * \f[ \phi^* = ( 0.47_{-0.10}^{+0.11} ) 10^{ ( -0.27 \pm 0.05 )( z - 6 )} 10^{-3} Mpc^{-3} \f]
     * \f[ \alpha = ( -1.87 \pm 0.05 ) + ( -0.10 \pm 0.03 ) ( z - 6 ) \f]
     * 
     * @param LL UV luminosity in magnitudes
     * 
     * @param zz redshift
     *
     * @return
     */
    double phi_Bouwens16 ( const double LL, const double zz );

    /**
     * @brief \f$\phi(L, z)\f$ luminosity function from Lapi et al., 2017
     *
     * Parameterisation of , 
     * valid in redshift range \f$4 \leq z \leq 8\f$:
     * with
     * 
     * @param LL UV luminosity in magnitudes
     * 
     * @param zz redshift
     *
     * @param param vector containing the parameters for computing the luminosity function
     *
     * @return
     */
    double phi_Lapi17 ( const double & LL, const double & zz,
			const std::vector< double > & param );
    // double phi_Lapi17 ( const double & LL, const double & zz,
    // 			   const std::vector< double > & philog,
    // 			   const std::vector< double > & sfrlog,
    // 			   const std::vector< double > & alphap );
    
    double phi_Lapi17_uv ( const double & LL, const double & zz );
    
    double phi_Lapi17_uvir ( const double & LL, const double & zz );

    double phi (const double & LL, const double & zz, const std::string & modelLF);

    /// @} End of luminosity functions
    
    /**
     * @name [ ... ] defined in src/[ ... ].cpp
     *
     * @{
     */

    double norm_z ( const double & zz,
		    const double & Lmin, const double & Lmax,
		    const double & zmin = 1.e-7, const double & zmax = 30. );

    double z_form ( const double & Mass, const double & ff,
		    const double & z_now, const double & z_max );

    double concentration_zhao09 ( const double Mass, const double Redshift );

    double concentration_shimizu03 ( const double Mass, const double Redshift );

    double concentration ( const double Mass, const double Redshift ) {

      return concentration_shimizu03 ( Mass, Redshift );

    }

    double density_profile_FS ( const double kk, const double Mass, const double Redshift = 1.e-7 );
    
    /// @} End of [ ... ]
  
  }; // endclass cosmology

} // endnamespace scam

/** 
 * @} End of Doxygen Groups
 */

#endif //__COSMOLOGY_H__
