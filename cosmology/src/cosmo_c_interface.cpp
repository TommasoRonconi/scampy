#include <cosmo_c_interface.h>
#include <cosmology.h>

extern "C" {
  
  // ========================================================================================

  cosmology_t create_cosmology ( const cosmo_params_t * const csmp,
				 const double * kh0, const double * Pk0,
				 const size_t size_k,
				 const double zmin, const double zmax,
				 const size_t thin ) {

    std::vector< double > kh0_v { kh0, kh0 + size_k };
    std::vector< double > Pk0_v { Pk0, Pk0 + size_k };
    std::unique_ptr< sico::cosmo_p > params {
      new sico::cosmo_p {
	csmp->Om_M,
	  csmp->Om_b, csmp->Om_L,
	  csmp->Om_n, csmp->Om_r, csmp->Om_K,
	  csmp->hh, csmp->sigma8 } };

    return new sico::cosmology { params, kh0_v, Pk0_v, zmin, zmax, thin };
    
  }    
  
  // ========================================================================================

  void free_cosmology ( cosmology_t cosmo ) {

    delete static_cast< sico::cosmology * >( cosmo );

  }
  
  // ========================================================================================

  double cosmo_Hz  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->H_z( zz );

  }
  
  // ========================================================================================

  double cosmo_dC  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->d_C( zz );

  }
  
  // ========================================================================================

  double cosmo_ddC  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->dd_C( zz );

  }
  
  // ========================================================================================

  double cosmo_comoving_volume_unit  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->comoving_volume_unit( zz );

  }
  
  // ========================================================================================

  double cosmo_comoving_volume  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->comoving_volume( zz );

  }
  
  // ========================================================================================

  double cosmo_cosmic_time  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->cosmic_time( zz );

  }
  
  // ========================================================================================

  double cosmo_rho_crit ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->rho_crit( zz );

  }
  
  // ========================================================================================

  double cosmo_rho_crit_comoving ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->rho_crit_comoving( zz );

  }
  
  // ========================================================================================

  double cosmo_OmegaM  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->OmegaM( zz );

  }
  
  // ========================================================================================

  double cosmo_deltac  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->deltac( zz );

  }
  
  // ========================================================================================

  double cosmo_Deltac_BN98  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->Delta_c_BryanNorman98( zz );

  }
  
  // ========================================================================================

  double cosmo_Deltac_NS98  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->Delta_c_NakamuraSuto98( zz );

  }
  
  // ========================================================================================

  double cosmo_DD  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->DD( zz );

  }
  
  // ========================================================================================

  double cosmo_gz  ( const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->gz( zz );

  }
  
  // ========================================================================================

  double cosmo_Pk  ( const double kk, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->Pk( kk, zz );

  }
  
  // ========================================================================================

  double cosmo_sigma2M  ( const double mm, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->sigma2M( mm, zz );

  }
  
  // ========================================================================================

  double cosmo_dndM  ( const double mm, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->dndM( mm, zz );

  }
  
  // ========================================================================================

  double cosmo_hbias  ( const double mm, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->hbias( mm, zz );

  }
  
  // ========================================================================================

  double cosmo_density_profile_FS  ( const double kk,
				     const double mm,
				     const double zz,
				     const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->density_profile_FS( kk, mm, zz );

  }
  
  // ========================================================================================

  double cosmo_dphidL  ( const double ll, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->phi_Bouwens15( ll, zz );

  }
  
  // ========================================================================================

  double cosmo_dphidL_Bouwens15  ( const double ll, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->phi_Bouwens15( ll, zz );

  }
  
  // ========================================================================================

  double cosmo_dphidL_Bouwens16  ( const double ll, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->phi_Bouwens16( ll, zz );

  }
  
  // ========================================================================================

  double cosmo_dphidL_Lapi17_uv  ( const double ll, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->phi_Lapi17_uv( ll, zz );

  }
  
  // ========================================================================================

  double cosmo_dphidL_Lapi17_uvir  ( const double ll, const double zz, const cosmology_t cosmo ) {

    return static_cast< sico::cosmology * >( cosmo )->phi_Lapi17_uvir( ll, zz );

  }
  
  // ========================================================================================

  double cosmo_dphidL_Lapi17  ( const double ll, const double zz,
				const double * param, const cosmology_t cosmo ) {

    std::vector< double > param_v { param, param + 12 };
    return static_cast< sico::cosmology * >( cosmo )->phi_Lapi17( ll, zz, param_v );

  }
  
  // ========================================================================================


} // end extern "C"
