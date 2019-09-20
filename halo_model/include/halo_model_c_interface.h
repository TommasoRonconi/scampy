#ifndef _HM_C_INTERFACE_
#define _HM_C_INTERFACE_

#include <cosmo_c_interface.h>
#include <occupation_c_interface.h>

// typedef void * hm_handler_t;
typedef void * H16_occupation_t;
typedef void * T10_occupation_t;
typedef void * halo_model_t;

#ifdef __cplusplus
extern "C" {
#endif

  // /**
  //  *  @name Halo Model Handler C-wrapping
  //  *  
  //  *  @{
  //  */

  // H16_occupation_t create_H16_occupation ( double DC,
  // 					   double Mmin,
  // 					   double sigma_logM,
  // 					   double M0,
  // 					   double M1,
  // 					   double alpha );

  // double Ncen_h16_ocp ( double Mh, H16_occupation_t ocp_h16 );
  
  // double Nsat_h16_ocp ( double Mh, H16_occupation_t ocp_h16 );
  
  // void free_H16_occupation_t ( H16_occupation_t ocp_h16 );

  // T10_occupation_t create_T10_occupation ( double Amin,
  // 					   double siglogA,
  // 					   double Asat,
  // 					   double alpsat );

  // double Ncen_t10_ocp ( double Mh, T10_occupation_t ocp_t10 );
  
  // double Nsat_t10_ocp ( double Mh, T10_occupation_t ocp_t10 );
  
  // void free_T10_occupation_t ( T10_occupation_t ocp_t10);
  
  // // hm_handler_t create_hm_handler ( double DC,
  // // 				   double Mmin,
  // // 				   double sigma_logM,
  // // 				   double M0,
  // // 				   double M1,
  // // 				   double alpha,
  // // 				   double redshift,
  // // 				   cosmology_t cosmo );

  // // unsigned int hm_thinness ( hm_handler_t hm_h );
  
  // // void free_hm_handler ( hm_handler_t hm_h );

  // /// @} End of Halo Model Handler C-wrapping

  /**
   *  @name Halo Model C-wrapping
   *  
   *  @{
   */

  halo_model_t create_halo_model_H16 ( H16_occupation_t ocp_h16,
				       cosmology_t cosmo,
				       const double redshift = 1.e-7,
				       const size_t thinness = 50 );

  halo_model_t create_halo_model_T10 ( T10_occupation_t ocp_t10,
				       cosmology_t cosmo,
				       const double redshift = 1.e-7,
				       const size_t thinness = 50 );
  
  // halo_model_t create_halo_model ( hm_handler_t hm_h );

  void free_halo_model ( halo_model_t hm );

  void set_parameters_hm_H16 ( double DC,
			       double Mmin,
			       double sigma_logM,
			       double M0,
			       double M1,
			       double alpha,
			       halo_model_t hm );

  void set_parameters_hm_T10 ( double Amin,
			       double siglogA,
			       double Asat,
			       double alpsat,
			       halo_model_t hm );

  // void set_parameters_hm ( double DC,
  // 			   double Mmin,
  // 			   double sigma_logM,
  // 			   double M0,
  // 			   double M1,
  // 			   double alpha,
  // 			   halo_model_t hm );

  size_t get_thinness_hm ( halo_model_t hm );

  double Ncen_hm ( double Mh, halo_model_t hm );
  
  double Nsat_hm ( double Mh, halo_model_t hm );

  double ng_hm ( halo_model_t hm );

  double bias_hm ( halo_model_t hm );

  double Mhalo_hm ( halo_model_t hm );

  double dngdM_hm ( double Mh, halo_model_t hm );

  void model_Pk_hm ( double * kv,
		     double * Pk,
		     halo_model_t hm );

  void model_Pk_1halo_hm ( double * kv,
			   double * Pk,
			   halo_model_t hm );

  void model_Pk_cs_hm ( double * kv,
			double * Pk,
			halo_model_t hm );
  
  void model_Pk_ss_hm ( double * kv,
			double * Pk,
			halo_model_t hm );
  
  void model_Pk_2halo_hm ( double * kv,
			   double * Pk,
			   halo_model_t hm );

  void model_Xi_hm ( double * rr,
		     double * Xi,
		     unsigned int size,
		     halo_model_t hm );

  void model_Xi_1halo_hm ( double * rr,
			   double * Xi,
			   unsigned int size,
			   halo_model_t hm );

  void model_Xi_2halo_hm ( double * rr,
			   double * Xi,
			   unsigned int size,
			   halo_model_t hm );

  void model_Wr_hm ( double * rp,
		     double * Wr,
		     unsigned int size,
		     halo_model_t hm );

  void model_Wr_1halo_hm ( double * rp,
			   double * Wr,
			   unsigned int size,
			   halo_model_t hm );

  void model_Wr_2halo_hm ( double * rp,
			   double * Wr,
			   unsigned int size,
			   halo_model_t hm );

  void model_Wt_hm ( double * tt,
		     double * Wt,
		     unsigned int size,
		     halo_model_t hm );

  void model_Wt_1halo_hm ( double * tt,
			   double * Wt,
			   unsigned int size,
			   halo_model_t hm );

  void model_Wt_cs_hm ( double * tt,
			double * Wt,
			unsigned int size,
			halo_model_t hm );
  
  void model_Wt_ss_hm ( double * tt,
			double * Wt,
			unsigned int size,
			halo_model_t hm );
  
  void model_Wt_2halo_hm ( double * tt,
			   double * Wt,
			   unsigned int size,
			   halo_model_t hm );

  void model_Wt_large_scale_hm ( double * tt,
				 double * Wt,
				 unsigned int size,
				 halo_model_t hm );

  double * get_kv_hm ( halo_model_t hm );

  /// @} End of Halo Model C-wrapping  

#ifdef __cplusplus
}
#endif

#endif //_HM_C_INTERFACE_
