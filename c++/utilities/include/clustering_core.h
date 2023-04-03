#ifndef __CLUSTERING_CORE__
#define __CLUSTERING_CORE__

// STL includes
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>

namespace utl {

  // inline float dist ( const float d ) { return d; }

  // inline float dist_periodic ( const float d, const float box ) {
  //   return d < 0.5*box ? d : box-d;
  // }
  
  //==================================================================================

  // template < typename fdist, typename ... Args >
  // std::vector< std::size_t > d2D_DD_generic ( const std::vector< float > & XX,
  // 					      const std::vector< float > & YY,
  // 					      const std::vector< float > & rbin,
  // 					      Args ... args ) {
    
  //   std::size_t size = XX.size();
  
  //   std::vector< std::size_t > NDD ( rbin.size() );
  //   float rmin = rbin.front(), rmax = rbin.back();
  //   float delta = std::log10(rmax/rmin)/rbin.size();

  //   for ( std::size_t ii = 0; ii < size; ++ii ) {
  //     float dx, dy, rr;
  //     std::size_t ib;
  //     for ( std::size_t jj = ii+1; jj < size; ++jj ) {
  // 	dx = fdist( XX[ii]-XX[jj], args ... );
  // 	dy = fdist( YY[ii]-YY[jj], args ... );
  // 	rr = std::sqrt( dx*dx + dy*dy );

  // 	if ( rmin <= rr && rr <= rmax ) {
  // 	  ib = int( std::log10( rr / rmin ) / delta );
  // 	  NDD[ ib ] += 1;
  // 	}
      
  //     } // endfor jj
  //   } // endfor ii
  
  //   return NDD;
  
  // }
  
  //==================================================================================
  //======================================= 2D =======================================
  //==================================================================================

  std::vector< std::size_t > d2D_DD ( const std::vector< float > & XX,
				      const std::vector< float > & YY,
				      const std::vector< float > & rbin );

  std::vector< std::size_t > d2D_DD_omp ( const std::vector< float > & XX,
					  const std::vector< float > & YY,
					  const std::vector< float > & rbin );

  std::vector< std::size_t > d2D_DR ( const std::vector< float > & X1,
				      const std::vector< float > & Y1,
				      const std::vector< float > & X2,
				      const std::vector< float > & Y2,
				      const std::vector< float > & rbin );

  std::vector< std::size_t > d2D_DR_omp ( const std::vector< float > & X1,
					  const std::vector< float > & Y1,
					  const std::vector< float > & X2,
					  const std::vector< float > & Y2,
					  const std::vector< float > & rbin );

  //==================================================================================
  //======================================= 3D =======================================
  //==================================================================================

  std::vector< std::size_t > d3D_DD ( const std::vector< float > & XX,
				      const std::vector< float > & YY,
				      const std::vector< float > & ZZ,
				      const std::vector< float > & rbin );

  std::vector< std::size_t > d3D_DD_omp ( const std::vector< float > & XX,
					  const std::vector< float > & YY,
					  const std::vector< float > & ZZ,
					  const std::vector< float > & rbin );

  std::vector< std::size_t > d3D_DR ( const std::vector< float > & X1,
				      const std::vector< float > & Y1,
				      const std::vector< float > & Z1,
				      const std::vector< float > & X2,
				      const std::vector< float > & Y2,
				      const std::vector< float > & Z2,
				      const std::vector< float > & rbin );

  std::vector< std::size_t > d3D_DR_omp ( const std::vector< float > & X1,
					  const std::vector< float > & Y1,
					  const std::vector< float > & Z1,
					  const std::vector< float > & X2,
					  const std::vector< float > & Y2,
					  const std::vector< float > & Z2,
					  const std::vector< float > & rbin );

} //endnamespace utl

#endif //__CLUSTERING_CORE__
