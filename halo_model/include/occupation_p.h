#ifndef __OCCUPATION_P__
#define __OCCUPATION_P__

/// gsl includes
#include <gsl/gsl_sf_erf.h>

/// internal includes
#include <utilities.h>

namespace sico {

  struct occupation_p {
    
    occupation_p () = default;

    occupation_p ( const occupation_p & ocp ) = default;

    virtual ~occupation_p () = default;

    virtual double PP ( const double AA, const double Amin, const double sigma_logA );

    virtual double Ncen ( const double AA ) = 0;

    virtual double Nsat ( const double AA ) = 0;

  }; // endstruct occupation_p
    
} // endnamespace sico

#endif //__OCCUPATION_P__
