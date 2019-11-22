#ifndef __ERROR_HANDLING__
#define __ERROR_HANDLING__

#include <string>

namespace scam_err {

  struct size_invalid {

    std::string message;
    size_invalid ( const std::string & s ) : message{ s } {}
      
  };

  struct out_of_bounds {

    std::string message;
    out_of_bounds ( const std::string & s ) : message{ s } {}
      
  };

  struct type_invalid {

    std::string message;
    type_invalid ( const std::string & s ) : message{ s } {}

  };

  struct gsl_fail {

    std::string message;
    gsl_fail ( const std::string & s ) : message{ s } {}
    
  };

  struct not_sorted {

    std::string message;
    not_sorted( const std::string & s ) : message{ s } {}
    
  };

} //endnamespace scam_err

#endif //__ERROR_HANDLING__
