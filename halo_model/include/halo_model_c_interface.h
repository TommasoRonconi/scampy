#ifndef _HM_C_INTERFACE_
#define _HM_C_INTERFACE_

#include <cosmo_c_interface.h>
#include <occupation_c_interface.h>

// typedef void * hm_handler_t;
typedef void * H16_occupation_t;
typedef void * T10_occupation_t;
typedef void * halo_model_t;
typedef void * cross_halo_model_t;

#ifdef __cplusplus
extern "C" {
#endif

  /**
   *  @name Halo Model C-wrapping
   *  
   *  @{
   */

  halo_model_t create_halo_model_H16 ( H16_occupation_t ocp_h16,
				       cosmology_t cosmo,
				       const double redshift,
				       const size_t thinness );

  halo_model_t create_halo_model_T10 ( T10_occupation_t ocp_t10,
				       cosmology_t cosmo,
				       const double redshift,
				       const size_t thinness );

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

  /**
   *  @name Cross Halo Model C-wrapping
   *  
   *  @{
   */

  cross_halo_model_t create_cross_halo_model_H16 ( H16_occupation_t ocp1_h16,
						   H16_occupation_t ocp2_h16,
						   cosmology_t cosmo,
						   const double redshift,
						   const size_t thinness );

  cross_halo_model_t create_cross_halo_model_T10 ( T10_occupation_t ocp1_t10,
						   T10_occupation_t ocp2_t10,
						   cosmology_t cosmo,
						   const double redshift,
						   const size_t thinness );

  void free_cross_halo_model ( cross_halo_model_t chm );

  void set_parameters_pop1_chm_H16 ( double DC,
				     double Mmin,
				     double sigma_logM,
				     double M0,
				     double M1,
				     double alpha,
				     cross_halo_model_t chm );

  void set_parameters_pop2_chm_H16 ( double DC,
				     double Mmin,
				     double sigma_logM,
				     double M0,
				     double M1,
				     double alpha,
				     cross_halo_model_t chm );
  
  void set_parameters_pop1_chm_T10 ( double Amin,
				     double siglogA,
				     double Asat,
				     double alpsat,
				     cross_halo_model_t chm );
  
  void set_parameters_pop2_chm_T10 ( double Amin,
				     double siglogA,
				     double Asat,
				     double alpsat,
				     cross_halo_model_t chm );

  size_t get_thinness_chm ( cross_halo_model_t chm );

  double * get_kv_chm ( cross_halo_model_t chm );

  double ng1_chm ( cross_halo_model_t chm );

  double ng2_chm ( cross_halo_model_t chm );

  void model_Pk_chm ( double * kv,
		      double * Pk,
		      cross_halo_model_t chm );

  void model_Pk_1halo_chm ( double * kv,
			    double * Pk,
			    cross_halo_model_t chm );
  
  void model_Pk_2halo_chm ( double * kv,
			    double * Pk,
			    cross_halo_model_t chm );

  void model_Xi_chm ( double * rr,
		      double * Xi,
		      unsigned int size,
		      cross_halo_model_t chm );

  void model_Xi_1halo_chm ( double * rr,
			    double * Xi,
			    unsigned int size,
			    cross_halo_model_t chm );

  void model_Xi_2halo_chm ( double * rr,
			    double * Xi,
			    unsigned int size,
			    cross_halo_model_t chm );

  /// @} End of Cross Halo Model C-wrapping  

#ifdef __cplusplus
}
#endif

#endif //_HM_C_INTERFACE_
