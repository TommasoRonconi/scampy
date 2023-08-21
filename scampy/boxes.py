"""Functions and classes to manipulate N-Dimensional boxes
"""

##############################################################################

# External imports
import numpy

##############################################################################

def shift_periodic_box ( box, shift, Lbox ) :
    """Given the cartesian coordinates from a periodic box applies a shift

    Parameters
    ----------
    box : ndarray
        Input comoving cartesian coordinates of particles in the box.
        Expects an array with with shape (Nparticles, Ndim).
        It will assume the coordinates are expressed in units 
        of [cMpc/h] and belong to the interval (0,Lbox)
    shift : scalar or iterable
        shift to apply along all or each direction.
        if it is a scalar, the same shift will be applied to
        all directions. If it is an iterable it should have
        the same length of Ndim, specifying the shift to apply
        along each dimension.
    Lbox : scalar
        box side lenght in [Mpc/h]
    
    Returns
    -------
    : ndarray
        Copy of the input array with the coordinate shift applied
    """
    Nobj, Ndim = box.shape
    box = numpy.array(box)
    if hasattr(shift, '__len__') and len(shift) != Ndim : 
        raise AttributeError( 'attribute ``shift`` has the wrong size' )
    box += shift
    wlower = box < 0.0
    box[wlower] = Lbox+box[wlower]
    wupper = box > Lbox
    box[wupper] = box[wupper]-Lbox
    return box

def centre_box ( box, Lbox ) :
    """Changes the coordinates of the particles in a comoving cubic box 
    moving the origin to the box centre.
    
    Parameters
    ----------
    box : ndarray
        Input comoving coordinates of particles in the box.
        Expects an array with with shape (Nparticles, Ndim).
        It will assume the coordinates are expressed in units 
        of [Mpc/h] and belong to the interval (0,Lbox)
    Lbox : scalar
        box side lenght in [Mpc/h]
    
    Returns
    -------
    : ndarray
        Copy of the input array with the coordinate centred
    """
    box = numpy.array(box)
    return box - 0.5*Lbox

def cartesian_to_polar ( coords, centre = (0.0,0.0,0.0) ) :
    """Converts cartesian coordinates (x,y,z) to polar coordinates (rho, theta, phi) 
    with physics convention (https://en.wikipedia.org/wiki/Spherical_coordinate_system).
    
    Parameters
    ----------
    coords : iterable
        an array or list or tuple with shape (3,Npoints), where Npoints is the number
        of coordinates to convert.
    
    centre : iterable
        (Optional, default = (0,0,0)) 
        an array or list or tuple with shape (3,) containing the coordinates of the
        centre of the polar coordinates system in cartesian coordinates
        
    Returns
    -------
    : iterable
        spherical coordinates of the input positions
        
    See Also
    --------
    polar_to_cartesian : inverse transform
    """

    try :
        x, y, z = numpy.asarray(coords)
    except Exception :
        raise
    try :
        x0, y0, z0 = numpy.asarray(centre)
    except Exception :
        raise

    rho   = numpy.sqrt( ( x-x0 )**2 + ( y-y0 )**2 + ( z-z0 )**2 )
    phi   = numpy.arctan2( ( y-y0 ), ( x-x0 ) )
    theta = numpy.arccos( ( z-z0 ) / rho )

    return rho, theta, phi

def polar_to_cartesian ( coords, centre = (0.0,0.0,0.0) ) :
    """Converts polar coordinates (rho, theta, phi) to cartesian coordinates (x,y,z)
    with physics convention (https://en.wikipedia.org/wiki/Spherical_coordinate_system).
    
    Parameters
    ----------
    coords : iterable
        an array or list or tuple with shape (3,Npoints), where Npoints is the number
        of coordinates to convert.
    centre : iterable
        (Optional, default = (0,0,0)) 
        an array or list or tuple with shape (3,) containing the coordinates of the
        centre of the polar coordinates system in cartesian coordinates
        
    Returns
    -------
    : iterable
        cartesian coordinates of the input positions
        
    See Also
    --------
    cartesian_to_polar : inverse transform
    """
    
    try :
        rho, theta, phi = coords
    except Exception :
        raise
    try :
        x0, y0, z0 = centre
    except Exception :
        raise
    
    x = rho * numpy.sin( theta ) * numpy.cos( phi )
    y = rho * numpy.sin( theta ) * numpy.sin( phi )
    z = rho * numpy.cos( theta )
    
    return x+x0, y+y0, z+z0

def cartesian_to_equatorial ( coords, centre = (0.0,0.0,0.0) ) : 
    """Converts cartesian coordinates (x,y,z) to astronomical coordinates (rho, theta, phi) 
    with equatorial convention.
    
    Parameters
    ----------
    coords : iterable
        an array or list or tuple with shape (3,Npoints), where Npoints is the number
        of coordinates to convert.
    
    centre : iterable
        (Optional, default = (0,0,0)) 
        an array or list or tuple with shape (3,) containing the coordinates of the
        centre of the equatorial coordinates system (the observer's position)
        in cartesian coordinates
        
    Returns
    -------
    : iterable
        equatorial coordinates of the input positions
        
    See Also
    --------
    equatorial_to_cartesian : inverse transform
    """
    
    try :
        x, y, z    = coords
    except Exception :
        raise
    try :
        x0, y0, z0  = centre
    except Exception :
        raise

    dd  = numpy.sqrt( ( x-x0 )**2 + ( y-y0 )**2 + ( z-z0 )**2 )
    ra  = numpy.arctan2( ( y-y0 ), ( x-x0 ) )
    dec = numpy.arcsin( ( z-z0 ) / dd )

    return dd, ra, dec

def equatorial_to_cartesian ( coords, centre = (0.0,0.0,0.0) ) :
    """Converts astronomical coordinates (d, ra, dec) to cartesian coordinates (x,y,z)
    with equatorial convention.
    
    Parameters
    ----------
    coords : iterable
        an array or list or tuple with shape (3,Npoints), where Npoints is the number
        of coordinates to convert.
    centre : iterable
        (Optional, default = (0,0,0)) 
        an array or list or tuple with shape (3,) containing the coordinates of the
        centre of the equatorial coordinates system (the observer's position) 
        in cartesian coordinates
        
    Returns
    -------
    : iterable
        cartesian coordinates of the input positions
        
    See Also
    --------
    cartesian_to_equatorial : inverse transform
    """
    
    try :
        dd, ra, dec = coords
    except Exception :
        raise
    try :
        x0, y0, z0  = centre
    except Exception :
        raise
    
    x = dd * numpy.cos( dec ) * numpy.cos( ra )
    y = dd * numpy.cos( dec ) * numpy.sin( ra )
    z = dd * numpy.sin( dec )
    
    return x+x0, y+y0, z+z0

##############################################################################
