#ifndef _CHM_C_INTERFACE_
#define _CHM_C_INTERFACE_

// #include <cosmo_c_interface.h>
// #include <occupation_c_interface.h>
#include <halo_model_c_interface.h>

typedef void * cross_halo_model_t;

#ifdef __cplusplus
extern "C" {
#endif

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

#endif //_CHM_C_INTERFACE_
