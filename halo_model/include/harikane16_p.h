#ifndef __HARIKANE16_P__
#define __HARIKANE16_P__

/// internal includes
#include <occupation_p.h>

namespace sico {

  struct harikane16_p : public occupation_p {

    double DC;
    double M_min;
    double sigma_logM;
    double M0;
    double M1;
    double alpha;
    
    harikane16_p () = default;

    harikane16_p ( const double & DC,
		   const double & Mmin,
		   const double & siglogM,
		   const double & M0,
		   const double & M1,
		   const double & alpha ) :
      DC { DC }, M_min { Mmin }, sigma_logM { siglogM }, M0 { M0 }, M1 { M1 }, alpha { alpha } {}

    harikane16_p ( const harikane16_p & handler ) = default;

    ~harikane16_p() = default;

    double Ncen ( const double Mhalo ) override;

    double Nsat ( const double Mhalo ) override;

  }; // endclass harikane16_p

} // endnamespace sico

#endif //__HARIKANE16_P__
