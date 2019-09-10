#include <bit_manipulator.h>

void sico::utl::showbits ( const long int & KK ) {
  
  std::bitset< bitset_size > manip_this { size_t( KK ) };
  
  std::cout << "value = " << KK << "\tbits = " << manip_this << "\n";

  return;

}  

bool sico::utl::bittest ( const long int & KK, const int pos ) {

  std::bitset< bitset_size > manip_this { size_t( KK ) };
  
  return manip_this[ pos ];

}

long int sico::utl::bitset ( const long int & KK, const int pos ) {

  std::bitset< bitset_size > manip_this { size_t( KK ) };

  manip_this.set( pos );

  return manip_this.to_ulong();

}

long int sico::utl::bitclr ( const long int & KK, const int pos ) {

  std::bitset< bitset_size > manip_this { size_t( KK ) };

  manip_this.reset( pos );

  return manip_this.to_ulong();

}

long int sico::utl::bits ( const long int & KK, const int pos, const int lenght ) {

  std::bitset< bitset_size > manip_this { size_t( KK ) }, ret_set;

  for ( int ii = 0; ii < lenght; ++ii )
    ret_set[ ii ] = manip_this[ ( ii + pos ) ];
  
  return ret_set.to_ulong();

}

long int sico::utl::bitshift ( const long int & KK, const int pos ) {

  if ( pos > 0 )
    return ( KK << pos );
  else if ( pos < 0 )
    return ( KK >> abs( pos ) );
  else
    return KK;

}
