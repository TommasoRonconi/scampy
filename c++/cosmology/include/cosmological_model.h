/**
 *  @file cosmology/include/cosmological_model.h
 *
 *  @brief 
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */

#ifndef __COSMOLOGICAL_MODEL__
#define __COSMOLOGICAL_MODEL__

#include <map>
#include <string>

#include <utilities.h>
#include <interpolation.h>

/// General namespace of the library 
namespace scam {

  // struct cosmo_p {

  //   /**
  //    * @name Density parameters
  //    *
  //    * @{ 
  //    */

  //   /// \f$\Omega_M\f$ matter density parameter
  //   double Om_M { 0.3 };
    
  //   /// \f$\Omega_b\f$ baryonic matter density parameter
  //   double Om_b { 0.045 };
    
  //   /// \f$\Omega_\Lambda\f$ dark-energy density parameter
  //   double Om_L { 0.7 };
    
  //   /// \f$\Omega_\nu\f$ neutrinos density parameter
  //   double Om_n { 0. };
    
  //   /// \f$\Omega_r\f$ radiation density parameter
  //   double Om_r { 0. };
    
  //   /// \f$\Omega_K\f$ curvature density parameter
  //   double Om_K { 0. };
    
  //   /// @} End of density parameters

  //   /// \f$h\f$ Hubble parameter, \f$H/(100\; [km \cdot s^{-1} \cdot Mpc^{-1}])\f$ (at \f$z = 0\f$)
  //   double hh { 0.7 };

  //   /// \f$\sigma_8\f$ normalization at \f$z = 0\f$
  //   double sigma8 { 0.8 };

  //   /**
  //    * @name DE equation of state 
  //    *
  //    * @brief CPL parameterisation: \f$ w_{DE} = w_0 + w_a \cdot z/(1+z) \f$
  //    *
  //    * @{ 
  //    */

  //   /// \f$w_0\f$ constant term
  //   double w_0 { -1. };

  //   /// \f$w_a\f$ slope
  //   double w_a { 0. };
    
  //   /// @} End of DE eq. of state parameters
    
  // };

  struct cosmo_model {

    std::map< std::string, float > param;
    double z_min { 0. }, z_max { 1.e+3 };
    size_t thinness { 10000 };
    
    double H0 = 100.;
    double t_H0, d_H0;

    /**
     * @name Interpolated functions
     *
     * @{ 
     */

    using interp_lin = class utl::interpolator< utl::lin_interp >;

    interp_lin Ez_f {}, zE_f {};

    /// @} End of Interpolated functions

    /**
     * @name Constructors/Destructors
     *
     * @{
     */

    cosmo_model ( const double Om_M = 0.3,
		  const double Om_b = 0.045,
		  const double Om_L = 0.7,
		  const double Om_n = 0.,
		  const double Om_r = 0.,
		  const double Om_K = 0.,
		  const double hh = 0.7,
		  const double sigma8 = 0.8,
		  const double w_0 = -1.,
		  const double w_a = 0.,
		  const double zmin = 0.,
		  const double zmax = 1.e+2,
		  const size_t thin = 10000 );

    cosmo_model ( const cosmo_model & other ) :
      param { other.param },
      z_min { other.z_min }, z_max { other.z_max },
      thinness { other.thinness } {
        
	Ez_f = other.Ez_f;
	zE_f = other.zE_f;
	H0 = other.H0;
	t_H0 = other.t_H0;
	d_H0 = other.d_H0;

      }

    virtual ~cosmo_model () = default;

    /// @} End of Ctor/Dtor

    void set_internal ();
    
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
    double d_H ( const double zz ) { return 1.e-3 * utl::cnst::cc / H_z( zz ); }
    
    double d_C ( const double zz ) { return d_H0 * zE_f.integrate( 0.0, zz ); }
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
    double rho_crit ( const double zz ) { return rho_crit_comoving( zz ) * param["hh"] * param["hh"]; }

    double OmegaM ( const double zz ) {

      return  param["Om_M"] / ( Ea2( zz ) * ( 1 + zz ) );
      
    }

    double Omegab ( const double zz ) {

      return  param["Om_b"] / ( Ea2( zz ) * ( 1 + zz ) );
      
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

  }; //endstruct cosmo_model

} //endnamespace scam

/** 
 * @} End of Doxygen Groups
 */

#endif //__COSMOLOGICAL_MODEL__
