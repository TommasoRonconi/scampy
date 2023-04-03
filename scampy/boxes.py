"""Functions and classes to manipulate N-Dimensional boxes
"""

##############################################################################

# External imports
import numpy

##############################################################################

def shift_periodic_box ( box, shift, Lbox ) :
    """Given coordinates from a periodic box applies a shift

    Parameters
    ----------
    box : ndarray
        Input comoving coordinates of particles in the box.
        Expects an array with with shape (Nparticles, Ndim).
        It will assume the coordinates are expressed in units 
        of [Mpc/h] and belong to the interval (0,Lbox)
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

##############################################################################
