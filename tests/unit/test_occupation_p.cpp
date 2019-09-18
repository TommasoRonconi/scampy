#include <harikane16_p.h>
#include <tinker10_p.h>

int main () {

  sico::harikane16_p ocp { 0.5, 1.e+11, 1., 1.e+11, 1.e+11, 1. };

  std::unique_ptr< sico::occupation_p > Pocp_harik { new sico::harikane16_p{ ocp } };

  double Mh = 1.e+12;
  
  std::cout << Mh << "\t" << Pocp_harik->Ncen( Mh ) << "\t" << Pocp_harik->Nsat( Mh ) << std::endl;
  
  std::unique_ptr< sico::occupation_p > Pocp_tinkr { new sico::tinker10_p{ 1.e+11, 1., 1.e+11, 1. } };
  
  std::cout << Mh << "\t" << Pocp_tinkr->Ncen( Mh ) << "\t" << Pocp_tinkr->Nsat( Mh ) << std::endl;
  
  return 0;

}

