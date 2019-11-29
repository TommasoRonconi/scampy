#include <halo_model_c_interface.h>
#include <occupation_c_interface.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main () {

  // Build object cosmology (with P_0(k) from file):
  cosmo_params_t * csmp = &cosmo_params_default;

  size_t size_k = 2048;
  double var1, var2;
  double * kh0 = ( double * ) malloc( sizeof( double ) * size_k );
  double * pk0 = ( double * ) malloc( sizeof( double ) * size_k );
  
  FILE * pkfile = fopen( "/home/tomi/phd/get_Pk/pk_harikane16_camb_extrapolated_z0.dat", "r" );
  
  size_k = 0;
  while ( fscanf( pkfile, "%lf\t%lf", &var1, &var2 ) == 2 ) {
    
    *( kh0 + size_k ) = var1;
    *( pk0 + size_k ) = var2;
    ++size_k;
    
  }

  cosmology_t cosmo = create_cosmology( csmp, kh0, pk0, size_k, 1.e-7, 1.e+7, 200 );

  free( kh0 );
  free( pk0 );

  // Build object halo_model:
  
  harikane16_t ocp_h16 = create_H16_occupation( 0.5, 1.e+11, 1., 1.e+10, 1.e+10, 1. );
  tinker10_t ocp_t10 = create_T10_occupation( 1.e+11, 1., 1.e+10, 1. );

  halo_model_t hm_h16 = create_halo_model_H16( ocp_h16, cosmo );
  halo_model_t hm_t10 = create_halo_model_T10( ocp_t10, cosmo );

  // compute galaxy density
  printf( "%f\t%f\n", ng_hm( hm_h16 ), ng_hm( hm_t10 ) );

  // compute galaxy power spectrum
  double * kk_h16 = ( double * ) malloc( sizeof( double ) * get_thinness_hm( hm_h16 ) );
  double * Pk_t10 = ( double * ) malloc( sizeof( double ) * get_thinness_hm( hm_h16 ) );
  double * kk_h16 = ( double * ) malloc( sizeof( double ) * get_thinness_hm( hm_h16 ) );
  double * Pk_t10 = ( double * ) malloc( sizeof( double ) * get_thinness_hm( hm_h16 ) );

  model_Pk_hm( kk_h16, Pk_h16, hm_h16 );
  model_Pk_hm( kk_t10, Pk_t10, hm_t10 );

  printf( "# k\tP_h16(k)\tP_t10(k)\n" );
  for ( unsigned int ii = 0; ii < get_thinness_hm( hm_h16 ); ii += 10 )
    printf( "%f\t%f\t%f\t%f\n",
	    *( kk_h16 + ii ), *( Pk_h16 + ii ),
	    *( kk_t10 + ii ), *( Pk_t10 + ii ) );

  /* // compute galaxy clustering */
  /* unsigned int nbin = 10; */
  /* double * rr = ( double * ) malloc( sizeof( double ) * nbin ); */
  /* double * Xi = ( double * ) malloc( sizeof( double ) * nbin ); */

  /* double min = 1.e-2, max = 1.e+2; */
  
  /* double bin = ( log( max ) - log( min ) )/( nbin - 1 ); */
  /* for ( unsigned int ii = 0; ii < nbin; ++ii ) *( rr + ii ) = exp( log( min ) + bin * ii ); */

  /* model_Xi_hm( rr, Xi, nbin, hm ); */

  /* printf( "# r\txi(r)\n" ); */
  /* for ( unsigned int ii = 0; ii < nbin; ++ii ) fprintf( stderr, "%f\t%f\n", *( rr + ii ), *( Xi + ii ) ); */

  /* // compute galaxy clustering */
  /* double * tt = ( double * ) malloc( sizeof( double ) * nbin ); */
  /* double * Wt = ( double * ) malloc( sizeof( double ) * nbin ); */

  /* double mint = 1.e-6, maxt = 1.e+0; */
  
  /* double bint = ( log( maxt ) - log( mint ) )/( nbin - 1 ); */
  /* for ( unsigned int ii = 0; ii < nbin; ++ii ) *( tt + ii ) = exp( log( mint ) + bint * ii ); */

  /* model_Wt_ss_hm( tt, Wt, nbin, hm ); */

  /* printf( "# t\tw(t)\n" ); */
  /* for ( unsigned int ii = 0; ii < nbin; ++ii ) fprintf( stderr, "%f\t%f\n", *( tt + ii ), *( Wt + ii ) ); */

  // free allocated stuff
  free( kk_h16 );
  free( Pk_h16 );
  free( kk_t10 );
  free( Pk_t10 );
  /* free( rr ); */
  /* free( Xi ); */
  /* free( tt ); */
  /* free( Wt ); */

  free_halo_model_H16( hm_h16 );
  free_halo_model_T10( hm_t10 );
  free_H16_occupation( ocp_h16 );
  free_T10_occupation( ocp_t10 );
  free_cosmology( cosmo );

  return 0;

}
