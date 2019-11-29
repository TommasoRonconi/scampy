#include <cosmo_c_interface.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main () {

  cosmo_params_t * csmp = &cosmo_params_default;

  printf( "Omega Matter =\t%f\n", csmp->Om_M );

  size_t size_k = 1000;
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

  fprintf( stdout, "H( z = 2 )\t=\t%e\n", cosmo_Hz( 2., cosmo ) );
  fprintf( stdout, "tt( z = 3 )\t=\t%e\n", cosmo_cosmic_time( 3., cosmo ) );
  fprintf( stdout, "dV( z = 2 )\t=\t%e\n", cosmo_comoving_volume_unit( 2., cosmo ) );
  fprintf( stdout, "V( z = 2 )\t=\t%e\n", cosmo_comoving_volume( 2., cosmo ) );

  free( kh0 );
  free( pk0 );
  free_cosmology( cosmo );

  return 0;

}
