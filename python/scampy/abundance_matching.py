import numpy
from scipy.integrate import quad

from scampy.interpolator import lin_interpolator

def differential_counts ( var, binning, fact = 1. ) :
    """ Function to obtain the differential distribution of a variable :math:`X`.
    It computes :math:`\\frac{dN}{dX} \\cdot const`, where :math:`dN` is the number
    of :math:`X` values in the bin, :math:`dX` it the size of the bin and :math:`const`
    is a constant value.

    Parameters
    ----------
    var : array-like
       the list of :math:`X` values
    binning : array-like
       it defines the bin edges, including the rightmost edge, 
       allowing for non-uniform bin widths.
    fact : scalar
       the value of :math:`const`, constant factor by which to multiply the result

    Returns
    -------
    array-like
       bin center (shape = binning.shape[ 0 ] - 1)
    array-like
       sequence of values with the values of 
       :math:`\\frac{dN}{dX} \\cdot const` (shape = binning.shape[ 0 ] - 1)
    array-like
       Poisson errors on the differential counts, defined as
       :math:`\\frac{\\sqrt{dN}}{dX} \\cdot const` (shape = binning.shape[ 0 ] - 1)    
    """
    
    counts, bins = numpy.histogram( var, bins = binning )

    size = bins.shape[ 0 ]
    M = 0.5 * ( bins[ :size - 1 ] + bins[ 1: ] )
    dM = ( bins[ 1: ] - bins[ :size - 1 ] )
    dndM = fact * counts / dM
    dndM_er = fact * numpy.sqrt( counts ) / dM
    
    return M, dndM, dndM_er

def cumulative_counts ( var, binning, fact = 1. ) :
    """ Function to obtain the cumulative distribution of a variable :math:`X`.
    It computes :math:`N(X>x_i) \\cdot const`, where :math:`N(X>x_i)` is the number
    of :math:`X` values greater than :math:`x_i` and :math:`const`
    is a constant value.

    Parameters
    ----------
    var : array-like
       the list of :math:`X` values
    binning : array-like
       it defines the binned values :math:`x_i`, allowing for non-uniform bin widths.
    fact : scalar
       the value of :math:`const`, constant factor by which to multiply the result

    Returns
    -------
    array-like
       sequence of values with the result of 
       :math:`N(X>x_i) \\cdot const` (shape = binning.shape)
    array-like
       Poisson errors on the cumulative counts, defined as
       :math:`\\sqrt{N(X>x_i)} \\cdot const` (shape = binning.shape)    
    """
    
    size = binning.shape[ 0 ]
    
    dN = numpy.empty( size )
    dN_er = numpy.empty( size )
    for ii in range( size ) :
        ww, = numpy.where( [ _b >= binning[ ii ] for _b in var ] )
        dN[ ii ] = ww.shape[ 0 ] * fact
        dN_er[ ii ] = numpy.sqrt( ww.shape[ 0 ] ) * fact
        
    return dN, dN_er

def cumulative_from_differential ( func, bins, fact = 1. ) :
    """ Function to obtain the cumulative distribution from a differential distribution :math:`f(X)`
    It computes 

    .. math:: F(X > x) \\equiv const \\cdot \\int_{-\\infty}^x f(X) dX

    where :math:`F(X>x)` is the cumulative distribution and :math:`const` is a constant value.

    Parameters
    ----------
    func : function
       the differential function :math:`f(X)`, it should take only one argument (a.k.a. the 
       value of :math:`X`). It allows :code:`lambda` expressions.
    bins : array-like
       it defines the list :math:`x_i` values, allowing for non-uniform bin widths.
    fact : scalar
       the value of :math:`const`, constant factor by which to multiply the result

    Returns
    -------
    array-like
       sequence of values with the result of :math:`F(X>x_i)` (shape = bins.shape)
    """

    size = len( bins )
    
    dN = numpy.empty( size )
    for ii in range( size ) :
        dN[ ii ], _ = quad( func, -numpy.inf, bins[ ii ] )
        dN[ ii ] *= fact
        
    return dN

def abundance_matching ( gxy_array, lum_func, minL = -30, maxL = -15,
                         minP = 1.e-20, maxP = 1., nbinL = 200, factL = 1.,
                         minM = None, maxM = None, nbinM = 10, factM = 1. ) :
    """ This function implements the Sub-halo Abundance Matching algorithm to an array of 
    objects of type :code:`galaxy` that are supposed to have been obtained from the HOD
    prescription of the :code:`catalogue` class.

    Parameters
    ----------
    gxy_array : array of galaxy objects
    
    lum_func : function

    minL : scalar

    maxL : scalar

    minP : scalar
    
    maxP : scalar
    
    nbinL : int
    
    factL : scalar
    
    minM : scalar or None
    
    maxM : scalar or None
    
    nbinM : int
    
    factM : scalar

    Returns
    -------
    array of galaxy objects
    
    Note
    ----
    stability not guaranteed, results have some randomness
    """

    # some checks and set-ups:
    if nbinL < 3 :
        raise Exception( 'nbinL should not be smaller than 3, you gave {:d}'.format( nbinL ) )
    masses = numpy.array( [ gxy.mass for gxy in gxy_array ] )
    if minM is None :
        minM = numpy.min( masses )
        minM -= 0.01 * minM
    if maxM is None :
        maxM = numpy.max( masses )
        maxM -= 0.01 * maxM
    
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
    dN, _ = cumulative_counts( numpy.array( [ gxy.mass for gxy in gxy_array ] ), MM, int( 1 ) )
    dN *= factM

    if dN.min() < numpy.min( dPhi[ wP ] ) : 
        raise Exception( 'Minimum cumulative number of objects ({0:.3e}) is smaller than minimum probability ({1:.3e})'.format( numpy.min( dN ),
                                                                                                                                numpy.min( dPhi[ wP ] ) )  )
    if dN.max() > numpy.max( dPhi[ wP ] ) :
        raise Exception( 'Maximum cumulative number of objects ({0:.3e}) is greater than maximum probability ({1:.3e})'.format( numpy.max( dN ),
                                                                                                                                numpy.max( dPhi[ wP ] ) )  )    
    # Variational bin matching
    for obj in gxy_array :

        # get index of higher probability:
        idx = int( ( numpy.log10( obj.mass ) - lminM ) * invbinM )
        
        # this exception here might be useless:
        if ( idx < 0 ) | ( nbinM <= idx ) :
            raise Exception( 'Index {0:d} is out of bounds. This galaxy has mass {1:.3e} Msol/h.'.format( idx, obj.mass ) )

        # if last bin, set to minimum probability:
        if idx + 1 != nbinM :
            Nl = dN[ idx + 1 ]
        else :
            Nl = invPhi_f.get_xmin()

        # get higher probability bin:
        Nh = dN[ idx ]

        # get random value between limits of binned observables: 
        L_low = invPhi_f( Nh )
        L_high = invPhi_f( Nl )
        obj.luminosity = numpy.random.uniform( L_low, L_high )


    return gxy_array
