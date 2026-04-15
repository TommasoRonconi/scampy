"""Mathematical and utility helper functions.

Provides coordinate distance, angle formatting, Fourier-space window
functions, trapezoidal integration, rejection sampling, and a
simple linear interpolator used throughout the halo-model pipeline.
"""

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
    """Fourier transform of a spherical top-hat window function.

    .. math::

        \\tilde{W}(kR) = \\frac{3\\,[\\sin(kR) - kR\\cos(kR)]}{(kR)^3}.

    Parameters
    ----------
    kR : scalar or array-like
        Dimensionless product of wavenumber and smoothing radius.

    Returns
    -------
    ndarray
        Window function values, same shape as ``kR``.
    """
    return 3 * ( numpy.sin(kR) - kR * numpy.cos(kR) ) / kR**3

def FT_tophat_D1 ( kR ) :
    """Derivative of the spherical top-hat window with respect to :math:`kR`.

    .. math::

        \\tilde{W}'(kR) = \\frac{3\\,[(kR)^2 - 3]\\sin(kR)
        + 9\\,kR\\cos(kR)}{(kR)^4}.

    Used by :func:`~scampy.power_spectrum.power_spectrum.dsigma2RdR`
    to evaluate :math:`\\mathrm{d}\\sigma^2/\\mathrm{d}R`.

    Parameters
    ----------
    kR : scalar or array-like
        Dimensionless product of wavenumber and smoothing radius.

    Returns
    -------
    ndarray
        Derivative values, same shape as ``kR``.
    """
    return ( 3 * ( kR**2 - 3. ) * numpy.sin( kR ) + 9 * kR * numpy.cos( kR ) ) / kR**4

def mod_erf ( x, factor=1.0 ) :
    """Normalised cumulative-like transform of an array via the error function.

    Centres ``x`` at its midpoint and rescales by ``factor * std(x)``,
    then applies the standard CDF of the normal distribution:

    .. math::

        \\mathrm{mod\\_erf}(x) =
        \\frac{1}{2}\\left[1 + \\mathrm{erf}\\!
        \\left(\\frac{x - x_M}{\\mathrm{factor}\\cdot\\sigma_x}
        \\right)\\right],

    where :math:`x_M = (x_\\min + x_\\max)/2` and
    :math:`\\sigma_x = \\mathrm{std}(x)`.

    Parameters
    ----------
    x : array-like
        Input array.
    factor : float, optional
        Scaling factor applied to the standard deviation (default: ``1.0``).

    Returns
    -------
    ndarray
        Values in :math:`[0, 1]`, same shape as ``x``.
    """
    from scipy.special import erf
    x = numpy.array( x )
    xM = 0.5*(x.min()+x.max())
    sigma = factor*x.std()
    t = (x-xM)/sigma
    return 0.5*(1+erf(t))

def truncated_gaussian(mean, std, lower, upper, size, rng = None, kw_rng = {'seed' : 555}):
    """Draw samples from a Gaussian distribution truncated to ``[lower, upper)``.

    Internally converts the bounds to the standard-normal parameterisation
    and delegates to :func:`scipy.stats.truncnorm`.

    Parameters
    ----------
    mean : float
        Mean of the untruncated Gaussian.
    std : float
        Standard deviation of the untruncated Gaussian.
    lower : float
        Lower bound of the truncation interval.
    upper : float
        Upper bound of the truncation interval (exclusive).
    size : int
        Number of samples to draw.
    rng : numpy.random.Generator or None, optional
        Random number generator.  If ``None`` (default) one is created
        via ``numpy.random.default_rng(**kw_rng)``.
    kw_rng : dict, optional
        Keyword arguments forwarded to ``numpy.random.default_rng``
        (default: ``{'seed': 555}``).

    Returns
    -------
    ndarray
        1-D array of ``size`` samples drawn from
        :math:`\\mathcal{N}(\\mu,\\sigma^2)` restricted to ``[lower, upper)``.
    """
    from scipy.stats import truncnorm
    if rng is None :
        rng = numpy.random.default_rng(**kw_rng)
    a, b = (lower - mean) / std, (upper - mean) / std  # Convert to standard normal space
    return truncnorm.rvs(a, b, loc=mean, scale=std, size=size, random_state=rng)

############################################################################################
# integration 

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

def linear_interpolation ( x, xgrid, ygrid ) :
    """Interpolates values of some input array on a user defined grid, where the
    x-domain must be evenly spaced.
    Specifically, computes:

    .. math::
    
        y = y_i + ( x - x_i ) \dfrac{y_{i+1}-y_i}{x_{i+1}-x_i}
    
    Parameters
    ----------
    x : 1d-array
        x-values
    xg : 1d-array
        x-axis of the interpolation grid
    yg : 1d-array
        y-axis of the interpolation grid
    
    Returns
    -------
    y : 1d-array
        interpolated y-values
    """

    x = numpy.asarray( x )
    xg = numpy.array( xgrid )
    yg = numpy.array( ygrid )
    if xg.ndim != 1 :
        raise AttributeError(
            "Currently working only on 1-dimensional arrays"
        )
    if not numpy.all( xg == sorted(xg) ) :
        raise RuntimeError(
            "The input grid must be sorted in ascending xg-order"
        )
    if xg.size != yg.size :
        raise AttributeError(
            "Arrays defining the interpolation grid should have same shape"
        )

    # Get extremes of interpolation interval
    # Note that this will work also for future N-dimensional broadcasting
    xmin, xmax = xg.min(axis=-1), xg.max(axis=-1)
    
    # Compute constant x-step
    dxg = numpy.abs( xmax - xmin ) / ( xg.size - 1 )
    
    # Get index of closest interval (accounting for extrapolation)
    idx = numpy.floor( ( x - xmin ) / dxg ).astype( int )
    idx = numpy.where( idx < xg.size - 1, idx, xg.size - 2 )
    idx = numpy.where( idx >= 0, idx, 0 )
    
    # Final sanity check
    if idx.size != x.size :
        raise RuntimeError( "Something went wrong." )
    
    # return interpolated values
    return (
        yg[idx] + ( x - xg[idx] ) * ( yg[ idx + 1 ] - yg[ idx ] ) /
        ( xg[ idx + 1 ] - xg[ idx ] )
    )

############################################################################################

def repeated_mask ( trues, falses ) :
    """Generate a boolean mask with repetitions
    (equivalent to numpy.repeat with alternating sequences of booleans)
    
    Parameters
    ----------
    trues : ndarray of int
        number of true values
    false : ndarray of int
        number of false values
    
    Returns
    -------
    mask : ndarray of bool
    
    Examples
    --------
    >>> from scampy.utilities.functions import repeated_mask
    >>> repeated_mask( [1,2], [2,3] )
    [ True False False True True False False False ]
    """

    trues  = numpy.array(trues)
    falses = numpy.array(falses)
    if trues.shape != falses.shape :
        raise RuntimeError( 'the two input arrays should have the same shape' )
    if numpy.any(trues<0) :
        raise RuntimeError( 'negative values in trues array' )
    if numpy.any(falses<0) :
        raise RuntimeError( 'negative values in falses array' )

    output = numpy.empty( numpy.sum(trues)+numpy.sum(falses), dtype = bool )
    index = 0
    for t, f in zip( trues, falses ) :
        output[index:index+t] = True
        index += t
        output[index:index+f] = False
        index += f
    return output

############################################################################################
