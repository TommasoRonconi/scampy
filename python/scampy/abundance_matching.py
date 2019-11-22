import numpy
from scipy.integrate import quad

from mockpy.interpolator import lin_interpolator

def differential_counts ( var, binning, fact ) :
    
    counts, bins = numpy.histogram( var, bins = binning )

    size = bins.shape[ 0 ]
    M = 0.5 * ( bins[ :size - 1 ] + bins[ 1: ] )
    dM = ( bins[ 1: ] - bins[ :size - 1 ] )
    dndM = fact * counts / dM
    dndM_er = fact * numpy.sqrt( counts ) / dM
    
    return M, dndM, dndM_er

def cumulative_counts ( var, binning, fact ) :
    
    size = binning.shape[ 0 ]
    
    dN = numpy.empty( size )
    dN_er = numpy.empty( size )
    for ii in range( size ) :
        ww, = numpy.where( [ _b >= binning[ ii ] for _b in var ] )
        dN[ ii ] = ww.shape[ 0 ] * fact
        dN_er[ ii ] = numpy.sqrt( ww.shape[ 0 ] ) * fact
        
    return dN, dN_er

def cumulative_from_differential ( func, bins, fact ) :

    size = len( bins )
    
    dN = numpy.empty( size )
    for ii in range( size ) :
        dN[ ii ], _ = quad( func, -numpy.inf, bins[ ii ] )
        dN[ ii ] *= fact
        
    return dN

def abundance_matching ( gxy_array, lum_func, minL = -30, maxL = -15,
                         minP = 1.e-20, maxP = 1., nbinL = 200, factL = 1.,
                         minM = 1.e+9, maxM = 1.e+15, nbinM = 50, factM = 1. ) :

    # Get cumulative luminosity function
    LL = numpy.linspace( minL, maxL, nbinL )
    dPhi = cumulative_from_differential( lum_func, LL, factL )

    # Get limits for inversion
    wP = numpy.where( [ ( minP <= _p ) & ( _p <= maxP ) for _p in dPhi ] )

    # Get inverse dPhi 
    invPhi_f = lin_interpolator( dPhi[ wP ], LL[ wP ] )

    # Get cumulative sub-halo mass function
    lminM = numpy.log10( minM )
    lmaxM = numpy.log10( maxM )
    invbinM = ( nbinM - 1 ) / ( lmaxM - lminM )
    MM = numpy.logspace( lminM, lmaxM, nbinM )
    dN, _ = cumulative_counts( numpy.array( [ gxy.mass for gxy in gxy_array ] ), MM, factM )

    # Variational bin matching
    for obj in gxy_array :
        idx = int( ( numpy.log10( obj.mass ) - lminM ) * invbinM ) + 1
        Nl = dN[ idx ]
        Nh = dN[ idx - 1 ]
        if ( invPhi_f.get_xmin() < Nl ) & ( Nh < invPhi_f.get_xmax() ) :
            L_low = invPhi_f( Nh )
            L_high = invPhi_f( Nl )
            obj.luminosity = numpy.random.uniform( L_low, L_high )


    return gxy_array
