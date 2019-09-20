#ifndef _OCP_C_INTERFACE_
#define _OCP_C_INTERFACE_

// typedef void * hm_handler_t;
typedef void * H16_occupation_t;
typedef void * T10_occupation_t;

#ifdef __cplusplus
extern "C" {
#endif

  /**
   *  @name Halo Model Handler C-wrapping
   *  
   *  @{
   */

  H16_occupation_t create_H16_occupation ( double DC,
					   double Mmin,
					   double sigma_logM,
					   double M0,
					   double M1,
					   double alpha );

  double Ncen_H16_ocp ( double Mh, H16_occupation_t ocp_h16 );
  
  double Nsat_H16_ocp ( double Mh, H16_occupation_t ocp_h16 );
  
  void free_H16_occupation ( H16_occupation_t ocp_h16 );

  T10_occupation_t create_T10_occupation ( double Amin,
					   double siglogA,
					   double Asat,
					   double alpsat );

  double Ncen_T10_ocp ( double Mh, T10_occupation_t ocp_t10 );
  
  double Nsat_T10_ocp ( double Mh, T10_occupation_t ocp_t10 );
  
  void free_T10_occupation ( T10_occupation_t ocp_t10);
  
  // hm_handler_t create_hm_handler ( double DC,
  // 				   double Mmin,
  // 				   double sigma_logM,
  // 				   double M0,
  // 				   double M1,
  // 				   double alpha,
  // 				   double redshift,
  // 				   cosmology_t cosmo );

  // unsigned int hm_thinness ( hm_handler_t hm_h );
  
  // void free_hm_handler ( hm_handler_t hm_h );

  /// @} End of Halo Model Handler C-wrapping


#ifdef __cplusplus
}
#endif

#endif //_OCP_C_INTERFACE_
