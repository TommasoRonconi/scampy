/**
 *  @file cosmology/include/cosmology_interface.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __COSMOLOGY_INTERFACE__
#define __COSMOLOGY_INTERFACE__

#include <utilities.h>
#include <interpolation.h>

/**
 *  @addtogroup scam
 *
 *  @{
 */

/// General namespace of the library 
namespace scam {

  struct cosmo_p {

    /**
     * @name Density parameters
     *
     * @{ 
     */

    /// \f$\Omega_M\f$ matter density parameter
    double Om_M { 0.3 };
    
    /// \f$\Omega_b\f$ baryonic matter density parameter
    double Om_b { 0.045 };
    
    /// \f$\Omega_\Lambda\f$ dark-energy density parameter
    double Om_L { 0.7 };
    
    /// \f$\Omega_\nu\f$ neutrinos density parameter
    double Om_n { 0. };
    
    /// \f$\Omega_r\f$ radiation density parameter
    double Om_r { 0. };
    
    /// \f$\Omega_K\f$ curvature density parameter
    double Om_K { 0. };
    
    /// @} End of density parameters

    /// \f$h\f$ Hubble parameter, \f$H/(100\; [km \cdot s^{-1} \cdot Mpc^{-1}])\f$ (at \f$z = 0\f$)
    double hh { 0.7 };

    /// \f$\sigma_8\f$ normalization at \f$z = 0\f$
    double sigma8 { 0.8 };

    /**
     * @name DE equation of state 
     *
     * @brief CPL parameterisation: \f$ w_{DE} = w_0 + w_a \cdot z/(1+z) \f$
     *
     * @{ 
     */

    /// \f$w_0\f$ constant term
    double w_0 { -1. };

    /// \f$w_a\f$ slope
    double w_a { 0. };
    
    /// @} End of DE eq. of state parameters
    
  };

  // struct matter_power_spec {
  
  //   /**
  //    * @name Power-spectrum properties
  //    *
  //    * @{ 
  //    */

  //   /// \f$\mathcal{A}_s\f$ initial scalar amplitude
  //   double A_s { 2.1e-9 };

  //   /// \f$k_p\f$ scalar pivot in \f$Mpc^{-1}\f$
  //   double pivot { 0.05 };

  //   /// \f$n_s\f$ primordial spectral index
  //   double n_s { 0.96 };

  //   /// \f$\tau\f$ Thomson scattering optical depth due to reionization
  //   double tau { 0.09 };
    
  //   /// @} End of power-spectrum properties

  //   /**
  //    * @name Power spectrum at \f$z = 0\f$
  //    *
  //    * @brief stored values of k/h and P(k)
  //    *
  //    * @{ 
  //    */

  //   /// sigma_8( z = 0 ) of the power spectrum
  //   double sigma8_cur { 1. };

  //   /// wavenumber values in cosmological units \f$[ h / Mpc ]\f$
  //   std::vector< double > kh0v {};

  //   /// corresponding value of the linear matter power spectrum at \f$z = 0\f$
  //   std::vector< double > Pk0v {};
    
  //   /// @} End of power spectrum at \f$z = 0\f$


  // }; // end struct matter_power_spec

  struct cosmo_model {

    std::unique_ptr< cosmo_p > param {};
    double z_min { 1.e-7 }, z_max { 1.e+7 };
    size_t thinness { 200 };
    
    double H0 = 100., ss8_P0 = -1;
    double h_1 = 1., h_2 = 1., h_3 = 1.;
    double t_H0, d_H0;

    /**
     * @name Interpolated functions
     *
     * @{ 
     */

    using interp_lin = class scam::utl::interpolator< scam::utl::gsl_lin_interp >;
    using interp_log = class scam::utl::interpolator< scam::utl::gsl_log_interp >;

    interp_log Ez_f {}, zE_f {};
    interp_log P0_f {};

    /// @} End of Interpolated functions

    /**
     * @name Constructors/Destructors
     *
     * @{
     */

    cosmo_model () {

      param = std::unique_ptr< cosmo_p > { new cosmo_p };
      std::vector< double > kh0 = scam::utl::log_vector( 200, 1.e-4, 1.e+4 );
      P0_f = interp_log { kh0, kh0 };
      set_internal();

    }

    cosmo_model ( const std::unique_ptr< cosmo_p > & parameters,
		  const std::vector< double > & kh0,
		  const std::vector< double > & Pk0,
		  const double zmin = 1.e-7, const double zmax = 1.e+7,
		  const size_t thin = 200 );

    cosmo_model ( const cosmo_model & other ) :
      param { new cosmo_p { * other.param } },
      z_min { other.z_min }, z_max { other.z_max },
      thinness { other.thinness } {
	
	P0_f = other.P0_f;
	Ez_f = other.Ez_f;
	zE_f = other.zE_f;
	H0 = other.H0;
	h_1 = other.h_1;
	h_2 = other.h_2;
	h_3 = other.h_3;
	ss8_P0 = other.ss8_P0;
	t_H0 = other.t_H0;
	d_H0 = other.d_H0;

      }

    virtual ~cosmo_model () = default;

    /// @} End of Ctor/Dtor

    

    void set_internal ();
    void compute_ss8_P0 ();
    
    /**
     * @name Cosmographic functions
     *
     * @{
     */

    /**
     * @brief Function to compute \f$E^2( z )\f$
     * 
     * The squared time-derivative of the scale factor logarithm \f$\log a( t )\t$.
     * Namely it computes
     * \f[ E^2(z) = ( 1 + z )^2 \bigl[ \sum_i \Omega_i ( 1 + z )^{ 1 + 3 w_i } \bigr] \f]
     * where \f$ w_M = 0 \f$, \f$ w_r = w_\nu = 1/3 \f$, \f$ w_K = -1/3 \f$
     *
     * @param zz the redshift at which to compute \f$E(z)\f$
     *
     * @return \f$ E( z ) \f$ computed in the most effecient way possible. In the case of
     *         flat cosmology and cosmological constant this function requires 2 memory allocations,
     *         5 multiplications and 5 additions.
     */
    virtual double Ez2 ( const double zz );
    virtual double Ea2 ( const double zz );

    double H_z ( const double zz ) { return H0 * Ez_f( zz ); }
    double HH ( const double zz ) { return H_z( zz ); }

    /**
     * @brief Function that returns the Hubble-distance \f$H_0\ [ pc ]\f$ at given redshift \f$z\f$
     *
     */
    double d_H ( const double zz ) { return 1.e-3 * scam::cc / H_z( zz ); }
    
    double d_C ( const double zz ) { return d_H0 * zE_f.integrate( 1.e-7, zz ); }
    double dd_C ( const double zz ) { return d_H0 * zE_f( zz ); }

    double d_A ( const double zz ) { return d_C( zz ) / ( 1 + zz ); }

    /// \f$\biggl[ \dfrac{dV}{dzd\Omega} \biggr] = [ Mpc^3 h^{-3} \text{rad}^{-1} ]\f$
    double comoving_volume_unit ( const double zz );
    
    /// \f$\biggl[ \dfrac{dV}{dzd\Omega} \biggr] = [ Mpc^3 h^{-3} \text{rad}^{-1} ]\f$
    double dV_dZdOmega ( const double zz,
			 const bool place_hold __attribute__(( unused )) = true ) {

      return comoving_volume_unit( zz );
      
    }

    /// \f$[ V ( z ) ] = [ Mpc^3 h^{ -3 } ]\f$
    double comoving_volume ( const double zz );

    /// the age of the Universe at the time the photons were emitted in \f$[ Gyr ]\f$
    double cosmic_time ( const double zz );

    /// @} End of cosmographic functions
    
    /**
     * @name Matter density related functions
     *
     * @{
     */

    /// critical density \f$[ M_\odot Mpc^{ -3 } h^2 ]\f$
    double rho_crit_comoving ( const double zz );

    /// critical density \f$[ M_\odot Mpc^{ -3 } ]\f$
    double rho_crit ( const double zz ) { return rho_crit_comoving( zz ) * param->hh * param->hh; }

    double OmegaM ( const double zz ) {

      return  param->Om_M / ( Ea2( zz ) * ( 1 + zz ) );
      
    }

    double Omegab ( const double zz ) {

      return  param->Om_b / ( Ea2( zz ) * ( 1 + zz ) );
      
    }

    double deltac ( const double & zz );

    double Delta_c_BryanNorman98 ( const double & zz );
    double Delta_c_NakamuraSuto98 ( const double & zz );
    double Delta_c ( const double & zz,
		     const std::string & author __attribute__(( unused )) = "NakamuraSuto98" ) {

      return Delta_c_NakamuraSuto98( zz );

    }

    /**
     * @brief computes the amplitude of the growing mode \f$D(z)\f$ as a function of redshift
     *
     * Readapted from Hamilton2001 ( arXiv:astro-ph/0006089v3, note to Eq. 1 ) computes
     * \f[ D(z) \equiv \dfrac{5 \Omega_{M,0}}{2} E( z ) \int_z^\infty \dfrac{(1 + z) dz}{E(z)} \f]
     *
     * @param zz the redshift at which to compute \f$D(z)\f$
     *
     * @return not normalized value of \f$D(z)\f$
     *
     * @warning this function uses stored logarithmic interpolators for computing \f$E(z)\f$
     *          thus it breaks at \f$z = 0\f$. Compute at a low value of redshift instead 
     *          (values accepted have to be \f$10^{-7} \le z \le 10^7\f$).
     *          The function does not throw any exception.
     */
    double DD ( const double & zz ) noexcept;

    /// growth factor (see function scam::cosmo_model::DD)
    double gz ( const double & zz ) noexcept { return ( 1 + zz ) * DD( zz ); }

    /// @} End of matter density related functions
    
    /**
     * @name Power spectrum dependent functions
     *
     * @{
     */

    double Pk_comoving ( const double & kk, const double & zz = 1.e-7 );
    double Pk ( const double & kk, const double & zz = 1.e-7 );

    // radius in Mpc/h 
    double sigma2R_comoving ( const double & Radius, const double & zz = 1.e-7 );
    double sigma2R ( const double & Radius, const double & zz = 1.e-7 );

    // radius in Mpc/h 
    double dsigma2RdR_comoving ( const double & Radius, const double & zz = 1.e-7 );
    double dsigma2RdR ( const double & Radius, const double & zz = 1.e-7 );

    // mass in solar masses ( not M_sol/h )
    double sigma2M ( const double & mm, const double & zz = 1.e-7 );

    // mass in solar masses
    double dsigma2MdM ( const double & mm, const double & zz = 1.e-7 );

    // Sheth & Tormen mass function
    double dndM_ST01 ( const double & mm, const double & zz = 1.e-7 );
    double dndM_Tinker08 ( const double & mm, const double & zz );
    double dndM_Behroozi13 ( const double & mm, const double & zz );
    double dndM ( const double & mm, const double & zz = 1.e-7 ) {

      return dndM_Behroozi13( mm, zz );
      // return dndM_Tinker08( mm, zz );
      // return dndM_ST01( mm, zz );

    }

    // Sheth, Mo & Tormen 2001 halo-bias
    double hbias_SMT01 ( const double & mm, const double & zz = 1.e-7 );
    // Tinker et al., 2010 halo-bias
    double hbias_Tinker10 ( const double & mm, const double & zz = 1.e-7 );
    double hbias ( const double & mm, const double & zz = 1.e-7 ) {

      return hbias_Tinker10( mm, zz );
      // return hbias_SMT01( mm, zz );

    }

    /// @} End of power spectrum dependent functions
    

  }; //endstruct cosmo_model

} //endnamespace scam

/** 
 * @} End of Doxygen Groups
 */

#endif //__COSMOLOGY_INTERFACE__
