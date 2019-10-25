#ifndef __COSMO_C_INTERFACE__
#define __COSMO_C_INTERFACE__

#include <stddef.h>

struct cosmo_params {

  /**
   * @name Density parameters
   *
   * @{ 
   */

  /// \f$\Omega_M\f$ matter density parameter
  double Om_M;
    
  /// \f$\Omega_b\f$ baryonic matter density parameter
  double Om_b;
    
  /// \f$\Omega_\Lambda\f$ dark-energy density parameter
  double Om_L;
    
  /// \f$\Omega_\nu\f$ neutrinos density parameter
  double Om_n;
    
  /// \f$\Omega_r\f$ radiation density parameter
  double Om_r;
    
  /// \f$\Omega_K\f$ curvature density parameter
  double Om_K;
    
  /// @} End of density parameters

  /// \f$h\f$ Hubble parameter, \f$H/(100\; [km \cdot s^{-1} \cdot Mpc^{-1}])\f$ (at \f$z = 0\f$)
  double hh;

  /// \f$\sigma_8\f$ normalization at \f$z = 0\f$
  double sigma8;

} cosmo_params_default = { 0.3, 0.045, 0.7, 0., 0., 0., 0.7, 0.8 };

typedef struct cosmo_params cosmo_params_t;

typedef void * cosmology_t;

#ifdef __cplusplus
extern "C" {
#endif

  cosmology_t create_cosmology ( const cosmo_params_t * const csmp,
				 const double * kh0, const double * Pk0,
				 const size_t size_k,
				 const double zmin, const double zmax,
				 const size_t thin );
  
  void free_cosmology ( cosmology_t cosmo );

  double cosmo_Hz ( const double zz, const cosmology_t cosmo );

  double cosmo_dC ( const double zz, const cosmology_t cosmo );

  double cosmo_ddC ( const double zz, const cosmology_t cosmo );

  double cosmo_comoving_volume_unit ( const double zz, const cosmology_t cosmo );

  double cosmo_comoving_volume ( const double zz, const cosmology_t cosmo );

  double cosmo_cosmic_time ( const double zz, const cosmology_t cosmo );

  double cosmo_rho_crit ( const double zz, const cosmology_t cosmo );

  double cosmo_rho_crit_comoving ( const double zz, const cosmology_t cosmo );

  double cosmo_OmegaM ( const double zz, const cosmology_t cosmo );

  double cosmo_deltac ( const double zz, const cosmology_t cosmo );

  double cosmo_Deltac_BN98 ( const double zz, const cosmology_t cosmo );

  double cosmo_Deltac_NS98 ( const double zz, const cosmology_t cosmo );

  double cosmo_DD ( const double zz, const cosmology_t cosmo );

  double cosmo_gz ( const double zz, const cosmology_t cosmo );

  double cosmo_Pk ( const double kk, const double zz, const cosmology_t cosmo );

  double cosmo_sigma2M ( const double mm, const double zz, const cosmology_t cosmo );

  double cosmo_dndM ( const double mm, const double zz, const cosmology_t cosmo );

  double cosmo_hbias ( const double mm, const double zz, const cosmology_t cosmo );

  double cosmo_density_profile_FS ( const double kk,
				    const double mm,
				    const double zz,
				    const cosmology_t cosmo );

  double cosmo_dphidL ( const double ll, const double zz, const cosmology_t cosmo );

  double cosmo_dphidL_Bouwens15 ( const double ll, const double zz, const cosmology_t cosmo );

  double cosmo_dphidL_Bouwens16 ( const double ll, const double zz, const cosmology_t cosmo );

  double cosmo_dphidL_Lapi17 ( const double ll, const double zz,
			       const double * param, const cosmology_t cosmo );
  
#ifdef __cplusplus
}
#endif

#endif //__COSMO_C_INTERFACE__
