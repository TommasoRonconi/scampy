#include <memory>
#include <cosmology.h>
#include <halo_model.h>
#include <harikane16_p.h>
#include <tinker10_p.h>

int main () {
  
  std::string filein = "/home/tomi/phd/get_Pk/pk_harikane16_camb_extrapolated_z0.dat";
  std::ifstream fin ( filein.c_str() );
  std::vector< double > kh0, Pk0;
  double var1, var2;
  while ( fin >> var1 >> var2 ) {
    
    kh0.push_back( var1 );
    Pk0.push_back( var2 );

  }
  fin.clear(); fin.close();
  std::unique_ptr< scam::cosmo_p > csmp { new scam::cosmo_p };
  scam::cosmology * pcosmo = new scam::cosmology { csmp, kh0, Pk0 };

  scam::tinker10_p * pocp_tinkr = new scam::tinker10_p { 1.e+11, 1., 1.e+11, 1. };
  scam::halo_model hm { std::make_shared< scam::tinker10_p >( ( *pocp_tinkr ) ),
      std::make_shared< scam::cosmology >( ( *pcosmo ) ), 1.e-7, 20 };

  hm.set_parameters( std::make_shared< scam::tinker10_p >( 1.e+10, 1., 1.e+11, 1. ) );

  std::cout << hm.ng() << "\t" << hm.Mhalo() << "\t" << hm.bias() << std::endl;

  delete pocp_tinkr;
  delete pcosmo;
  
  return 0;

}

