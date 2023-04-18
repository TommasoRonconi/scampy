import numpy

############################################################################################
# Utility functions

def dist_periodic ( coord1, coord2, side = None ) :
    """ Returns the distance in a periodic box
    
    Parameters
    ----------
    coord1, coord2 : array-like or scalar
        position of the 2 sets of coordinates
    side : scalar
        side lenght of the periodic box, if it
        is not None (default), the periodic distance
        is computed
    
    Returns
    -------
    : array-like or scalar
        distance between coord1 and coord2
    """
    
    coord1 = numpy.asarray( coord1 )
    coord2 = numpy.asarray( coord2 )
    diff = abs(coord1 - coord2)
    if side is not None :
        med = 0.5 * side
        diff = numpy.array( [ [ side-d if d > med else d for d in dd ] for dd in diff ] )
    return numpy.sqrt(numpy.sum( diff * diff, axis=-1 ))

def str_format_deg ( theta, rad = False ) :
    """String formatting of angles in (degrees-primes-seconds)
    
    Parameters
    ----------
    theta : scalar float
        angle to be formatted
    rad : bool
        (Optional, default = False) 
        whether theta has been provided in radians or not
    
    Returns
    -------
    : str
        formatted string
    """
    if rad :
        theta /= dtor
    deg = int( theta )
    pri = theta - deg
    sec = ( pri - int( pri ) ) * 60
    pri = int( pri * 60 )
    return f'{deg:d}° {pri:d}\' {sec:.1f}\"'

############################################################################################
# Matemathical functions

def FT_tophat ( kR ) :
    """Fourier space transform of a top-hat filter function
    """
    return 3 * ( numpy.sin(kR) - kR * numpy.cos(kR) ) / kR**3

def FT_tophat_D1 ( kR ) :
    """Fourier space transform of a top-hat filter function, first order derivative
    """
    return ( 3 * ( kR**2 - 3. ) * numpy.sin( kR ) + 9 * kR * numpy.cos( kR ) ) / kR**4 

############################################################################################

def trap_int ( xx, yy ) :
    """ Trapezoid integration

    Parameters
    ----------
    xx : array-like
       x-domain grid
    yy : array-like
       y-domain grid
    
    Returns
    -------
    : float 
       the integral along the whole x-domain
    """
    xx = numpy.asarray(xx)
    yy = numpy.asarray(yy)
    return numpy.sum( 0.5 * ( yy[:-1] + yy[1:] ) * ( xx[1:] - xx[:-1] ), axis = 0 )

############################################################################################

def rejection_sampling ( pdf, 
                         xmin = 0.0, xmax = 1.0, 
                         umin = 0.0, umax = 1.0,
                         logu = False,
                         size = 1, rng = None, 
                         verbose = False
                       ) :
    """Samples the ``X`` random variable from a custom target PDF using rejection sampling.
    Proposals for the random variable are drawn from a uniform distribution defined between 
    (xmin, xmax).
    A proposal ``x_j`` is accepted if the corresponding value of the proposal PDF ``u_j``
    is lower than the value of the target PDF computed at the proposal random variable value, 
    i.e. if ``u_j < PDF(x_j)``.
    Values of the proposal PDF are drawn from a uniform distribution defined between (umin, umax).
    The higher the value of umax with respect to the maximum of the target PDF, the lower
    the acceptance fraction (but a higher value of umax might produce more reliable samples).
    
    Parameters
    ----------
    pdf : callable
        the callable function defining the custom PDF (it must have a ``__call__`` method)
    xmin : float
        (Optional, default = 0.0) minimum value of the random variable X
    xmax : float
        (Optional, default = 1.0) maximum value of the random variable X
    umin : float
        (Optional, default = 0.0) minimum value of the proposal PDF 
    umax : float
        (Optional, default = 1.0) maximum value of the proposal PDF
    logu : bool
        (Optional, default = False) if True, ``umin`` and ``umax`` are interpreted 
        as log10(umin) and log10(umax)
    size : int
        (Optional, default = 1) size of the output sample
    rng : int or random generator
        (Optional, default = None) either an integer or a random number generator.
        If an integer is passed, a generator will be instantiated with ``seed=rng``.
        If None is passed, a random number generator will be instantiated 
        with randomly generated seed.
    verbose : bool
        (Optional, default = False) if True prints on stdout informations about the sampling
    
    Returns
    -------
    : 1d-array
        A 1D array with ``size`` samples of the random variable extracted from the target PDF
        
    Note
    ----
    One can run the Kolmogorov-Smirnov Test to verify that the sampled random variable has
    a PDF that approximates the target PDF.
    """
    
    if not hasattr( pdf, '__call__' ) :
        raise ValueError( 'pdf argument should be callable' )
    
    if rng is None :
        rng = numpy.random.default_rng()
    elif isinstance( rng, int ) :
        rng = numpy.random.default_rng( seed = rng )
    elif not isinstance( rng, numpy.random.Generator ) :
        raise ValueError( 'rng argument should be one among None, integer, random generator')
    
    sample = numpy.zeros( size )
    notassigned = numpy.ones_like( sample, dtype=bool )
    it = 0
    if verbose : acceptfrac = 0.0
    while notassigned.any() :
        x = rng.uniform( low = xmin, high = xmax, size = size )
        u = rng.uniform( low = umin, high = umax, size = size )
        if logu : u = 10**u
        assign = ( u < pdf(x) )
        sample[assign] = x[assign]
        notassigned &= ~assign
        it += 1
        if verbose :
            acceptfrac += assign.sum()/assign.size
            print(f'Assigned {100*(~notassigned).sum()/notassigned.size:.2f} %', end='\r')
    if verbose : 
        print(
            f'\nDone in {it} iterations'
            f'\nAverage percent acceptance: {100*acceptfrac/it:.2f}%'
         )
    
    return sample

############################################################################################
