#include <clustering_core.h>
#include <omp.h>

//==================================================================================
//======================================= 2D =======================================
//==================================================================================

std::vector< std::size_t > utl::d2D_DD ( const std::vector< float > & XX,
					 const std::vector< float > & YY,
					 const std::vector< float > & rbin ) {
  
  std::size_t size = XX.size();
  
  std::vector< std::size_t > NDD ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

  for ( std::size_t ii = 0; ii < size; ++ii ) {
    float dx, dy, rr;
    std::size_t ib;
    for ( std::size_t jj = ii+1; jj < size; ++jj ) {
      dx = XX[ii]-XX[jj];
      dy = YY[ii]-YY[jj];
      rr = std::sqrt( dx*dx + dy*dy );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
	NDD[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDD;
  
}

std::vector< std::size_t > utl::d2D_DD_omp ( const std::vector< float > & XX,
					     const std::vector< float > & YY,
					     const std::vector< float > & rbin ) {
  
  std::size_t size = XX.size();
  
  std::vector< std::size_t > NDD ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

#pragma omp parallel for shared(NDD)
  for ( std::size_t ii = 0; ii < size; ++ii ) {
    float dx, dy, rr;
    std::size_t ib;
    for ( std::size_t jj = ii+1; jj < size; ++jj ) {
      dx = XX[ii]-XX[jj];
      dy = YY[ii]-YY[jj];
      rr = std::sqrt( dx*dx + dy*dy );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
#pragma omp atomic
	NDD[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDD;
  
}

std::vector< std::size_t > utl::d2D_DR ( const std::vector< float > & X1,
					 const std::vector< float > & Y1,
					 const std::vector< float > & X2,
					 const std::vector< float > & Y2,
					 const std::vector< float > & rbin ) {
  
  std::size_t size1 = X1.size();
  std::size_t size2 = X2.size();
  
  std::vector< std::size_t > NDR ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

  for ( std::size_t ii = 0; ii < size1; ++ii ) {
    float dx, dy, rr;
    std::size_t ib;
    for ( std::size_t jj = 0; jj < size2; ++jj ) {
      dx = X1[ii]-X2[jj];
      dy = Y1[ii]-Y2[jj];
      rr = std::sqrt( dx*dx + dy*dy );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
	NDR[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDR;
  
}

std::vector< std::size_t > utl::d2D_DR_omp ( const std::vector< float > & X1,
					     const std::vector< float > & Y1,
					     const std::vector< float > & X2,
					     const std::vector< float > & Y2,
					     const std::vector< float > & rbin ) {
  
  std::size_t size1 = X1.size();
  std::size_t size2 = X2.size();
  
  std::vector< std::size_t > NDR ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

#pragma omp parallel for shared(NDR)
  for ( std::size_t ii = 0; ii < size1; ++ii ) {
    float dx, dy, rr;
    std::size_t ib;
    for ( std::size_t jj = 0; jj < size2; ++jj ) {
      dx = X1[ii]-X2[jj];
      dy = Y1[ii]-Y2[jj];
      rr = std::sqrt( dx*dx + dy*dy );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
#pragma omp atomic
	NDR[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDR;
  
}

//==================================================================================
//======================================= 3D =======================================
//==================================================================================

std::vector< std::size_t > utl::d3D_DD ( const std::vector< float > & XX,
					 const std::vector< float > & YY,
					 const std::vector< float > & ZZ,
					 const std::vector< float > & rbin ) {
  
  std::size_t size = XX.size();
  
  std::vector< std::size_t > NDD ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

  for ( std::size_t ii = 0; ii < size; ++ii ) {
    float dx, dy, dz, rr;
    std::size_t ib;
    for ( std::size_t jj = ii+1; jj < size; ++jj ) {
      dx = XX[ii]-XX[jj];
      dy = YY[ii]-YY[jj];
      dz = ZZ[ii]-ZZ[jj];
      rr = std::sqrt( dx*dx + dy*dy + dz*dz );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
	NDD[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDD;
  
}

std::vector< std::size_t > utl::d3D_DD_omp ( const std::vector< float > & XX,
					     const std::vector< float > & YY,
					     const std::vector< float > & ZZ,
					     const std::vector< float > & rbin ) {
  
  std::size_t size = XX.size();
  
  std::vector< std::size_t > NDD ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

#pragma omp parallel for shared(NDD)
  for ( std::size_t ii = 0; ii < size; ++ii ) {
    float dx, dy, dz, rr;
    std::size_t ib;
    for ( std::size_t jj = ii+1; jj < size; ++jj ) {
      dx = XX[ii]-XX[jj];
      dy = YY[ii]-YY[jj];
      dz = ZZ[ii]-ZZ[jj];
      rr = std::sqrt( dx*dx + dy*dy + dz*dz );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
#pragma omp atomic
	NDD[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDD;
  
}

std::vector< std::size_t > utl::d3D_DR ( const std::vector< float > & X1,
					 const std::vector< float > & Y1,
					 const std::vector< float > & Z1,
					 const std::vector< float > & X2,
					 const std::vector< float > & Y2,
					 const std::vector< float > & Z2,
					 const std::vector< float > & rbin ) {
  
  std::size_t size1 = X1.size();
  std::size_t size2 = X2.size();
  
  std::vector< std::size_t > NDR ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

  for ( std::size_t ii = 0; ii < size1; ++ii ) {
    float dx, dy, dz, rr;
    std::size_t ib;
    for ( std::size_t jj = 0; jj < size2; ++jj ) {
      dx = X1[ii]-X2[jj];
      dy = Y1[ii]-Y2[jj];
      dz = Z1[ii]-Z2[jj];
      rr = std::sqrt( dx*dx + dy*dy + dz*dz );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
	NDR[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDR;
  
}

std::vector< std::size_t > utl::d3D_DR_omp ( const std::vector< float > & X1,
					     const std::vector< float > & Y1,
					     const std::vector< float > & Z1,
					     const std::vector< float > & X2,
					     const std::vector< float > & Y2,
					     const std::vector< float > & Z2,
					     const std::vector< float > & rbin ) {
  
  std::size_t size1 = X1.size();
  std::size_t size2 = X2.size();
  
  std::vector< std::size_t > NDR ( rbin.size() );
  float rmin = rbin.front(), rmax = rbin.back();
  float delta = std::log10(rmax/rmin)/rbin.size();

#pragma omp parallel for shared(NDR)
  for ( std::size_t ii = 0; ii < size1; ++ii ) {
    float dx, dy, dz, rr;
    std::size_t ib;
    for ( std::size_t jj = 0; jj < size2; ++jj ) {
      dx = X1[ii]-X2[jj];
      dy = Y1[ii]-Y2[jj];
      dz = Z1[ii]-Z2[jj];
      rr = std::sqrt( dx*dx + dy*dy + dz*dz );

      if ( rmin <= rr && rr <= rmax ) {
	ib = int( std::log10( rr / rmin ) / delta );
#pragma omp atomic
	NDR[ ib ] += 1;
      }
      
    } // endfor jj
  } // endfor ii
  
  return NDR;
  
}

//==================================================================================
//==================================================================================
