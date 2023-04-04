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
