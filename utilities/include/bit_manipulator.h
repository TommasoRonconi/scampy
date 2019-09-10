#ifndef __BITMANIP__
#define __BITMANIP__

#include <bitset>
#include <utilities.h>

namespace sico {

  namespace utl {

    static const size_t bitset_size = 64;

    /// prints the bits of KK
    void showbits ( const long int & KK );

    /// return true if pos bit of KK is 1
    bool bittest ( const long int & KK, const int pos );

    /// set pos bit of KK to 1
    long int bitset ( const long int & KK, const int pos );

    /// set pos bit of KK to 0
    long int bitclr ( const long int & KK, const int pos );

    /// from pos bit of KK extract lenght bits
    long int bits ( const long int & KK, const int pos, const int lenght );

    /// move the bits of KK to the right ( pos > 0 ) or to the left ( pos < 0 )
    long int bitshift ( const long int & KK, const int pos );

  } // end namespace utl

} // end namespace sico


#endif //__BITMANIP__
