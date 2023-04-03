"""Functions and classes to generate and manipulate lightcones 
"""

##############################################################################

# External imports
import numpy

# Internal imports
import scampy.cosmology as c
from scampy.boxes import centre_box, shift_periodic_box
from scampy.utilities.interpolation import lin_interp as lint
from scampy.utilities.constants import degrees_to_radians as dtor

##############################################################################
# Functions

def sequential_lightcones_limits ( zmax, Lbox, d2z = None, centre = True ) :
    """Computes the central redshift of boxes of given size necessary
    to build a lightcone from the resulting snapshots up to some given 
    redshift.
    
    Parameters
    ----------
    zmax : float
        Maximum redshift of the lightcone
    Lbox : float
        comoving side size of the simulation box
    d2z : function or dictionary or None
        (Optional, default = ``'None'``) either a function 
        taking a comoving distance as argument and returning 
        the corresponding redshift or a dictionary with fields
        ``'cosmo'`` with an instance of object of class scampy.cosmology.model
        and ``'zgrid'`` the redshift grid on witch to interpolate
        the necessary function d2z (the thinner the grid, 
        the more accurate the interpolation) 
    centre : bool
        (Optional, default = ``True``) if true associate the values
        computed to the centre of the box, otherwise, they will
        be associated to the closest side.

    Returns
    -------
    zeta : 1d-array
        array with the sequential redshifts of boxes necessary
        to build a lightcone
    dbox : 1d-array
        array with the sequential comoving distances of the boxes 
        necessary to build a lightcone
    """
    
    # Preparatory steps, checking inputs and generating grids if necessary
    if d2z is None or not hasattr( d2z, "__call__" ) :
        cosmo = None
        zz = None
        if isinstance( d2z, dict ) :
            try :
                cosmo = d2z['cosmo']
                if not isinstance( cosmo, c.model ) :
                    raise TypeError( "Element 'cosmo' of attribute d2z "
                                     "should be an instance of class scampy.cosmology.model" )
            except KeyError :
                cosmo = c.model()
            try :
                zz = d2z['zgrid']
                if not hasattr( zz, '__len__' ) :
                    raise TypeError( "Element 'zgrid' of attribute d2z "
                                     "should be an iterable" )
            except KeyError :
                zz = numpy.arange( 0.0, 20.0, 1.e-2 )
        else :
            cosmo = c.model()
            zz = numpy.arange( 0.0, 20.0, 1.e-2 )
            
        dC = cosmo.dC(zz)
        d2z = lint( dC, zz )
        
    # Actual computation
    dbox = [int(centre) * 0.5 * Lbox]
    zeta = [d2z(dbox[-1])]
    while zeta[-1] < zmax :
        dbox += [dbox[-1]+Lbox]
        zeta += [d2z(dbox[-1])]
        
    return numpy.asarray( zeta ), numpy.asarray( dbox )

##############################################################################

def ang_size ( lC, z, cosmo, rad = True ) :
    """Converts a comoving distance, perpendicular to the LoS,
    at a given redshift to an angle.
    
    Parameters
    ----------
    lC : scalar or 1d-array
        Comoving distance in [Mpc/h], assumed perpendicular to the LoS.
        Note that if this is a 1D-array, then argument ```z``` must be 
        a scalar.
    z : scalar or 1d-array
        redshift along the LoS.
        Note that if this is a 1D-array, then argument ```lC``` must be 
        a scalar.
    cosmo : scampy.cosmology.model
        Instance of an object of type scampy.cosmology.model
    rad : bool
        (Optional, default = True) whether to output an angle in
        radians (rad = True), or degrees (rad = False)
    
    Returns
    -------
    : scalar or nd-array
        Angle corresponding to the input distances at the input redshifts
    
    Warning
    -------
    Does not support broadcasting on two dimensions. Only one among
    lC and z can be a 1d-array, not both.
    """
    theta = lC / cosmo.dC(z)
    if not rad :
        return theta / dtor
    return theta

##############################################################################

def size_ang ( theta, z, cosmo, rad = True ) :
    """Converts an angle at a given redshift to a
    comoving distance, perpendicular to the LoS.
    
    Parameters
    ----------
    theta : scalar or 1d-array
        Viewing angle. 
        Note that if this is a 1D-array, then argument ```z``` must be 
        a scalar.
    z : scalar or 1d-array
        redshift along the LoS. 
        Note that if this is a 1D-array, then argument ```theta``` must be 
        a scalar.
    cosmo : scampy.cosmology.model
        Instance of an object of type scampy.cosmology.model
    rad : bool
        (Optional, default = True) whether the input angle is
        given in radians (rad = True) or degrees (rad = False)
    
    Returns
    -------
    : scalar or nd-array
        Comoving distances in [Mpc/h], assumed perpendicular to the LoS, 
        computed at the input redshifts.
    
    Warning
    -------
    Does not support broadcasting on two dimensions. Only one among
    theta and z can be a 1d-array, not both.
    """
    if not rad :
        theta *= dtor
    return theta * cosmo.dC(z)

##############################################################################

def to_cone ( box, zbox, Lbox, cosmo,
              funcs = None, theta = None,
              zmin = None, zmax = None,
              return_indices = False ) :
    """Converts the coordinates of particles in a comoving cubic box 
    to coordinates in the lightcone.
    
    Paramters
    ---------
    box : 3d-array
        Input comoving coordinates of particles in the box.
        Expects an array with with shape (Nparticles, 3).
        It will assume the coordinates are expressed in units 
        of [Mpc/h] and belong to the interval (0,Lbox)
    zbox : scalar
        Redshift associated to the centre of the box
    Lbox : scalar
        box side lenght in [Mpc/h]
    cosmo : scampy.cosmology.model
        Instance of an object of type scampy.cosmology.model
    funcs : tuple
        (Optional, default = None) tuple (function,inverse) computing comoving distance
        from redshift and vice-versa. The functions should be vectorizable
        and returning scalar of numpy.ndarray. If the tuple is not provided
        attribute ```cosmo``` will be used to allocate a couple 
        of interpolators serving this purpose.
    theta : scalar
        (Optional, default = None) viewing angle in radians.
        If None is passed (default behaviour), the largest possible viewing angle
        associated with the provided Lbox and zbox is assumed.
    zmin : scalar or None
        (Optional, default = None) if passed, only selects coordinates with redshift
        larger than this value (exclusive)
    zmax : scalar or None
        (Optional, default = None) if passed, only selects coordinates with redshift
        smaller than this value (inclusive)
    return_indices : bool
        (Optional, default = False) when True, also return the mask that, applied to
        the input ``box`` argument, selects only coordinates within the light-cone
    
    Returns
    -------
    cone_coords : 3d-array
        Spherical coordinates of the generated lightcone with shape (Nview, 3),
        where Nview <= Nparticles is the number of particles within the field of view.
        The light-cone provides, for each particle within the field of view, 
        a (radians, radians, redshift) array.
    cone_indices : bool 1d-array
        the boolean mask that, applied to the input ``box`` argument, selects only
        coordinates within the light cone. Only provided if ``return_indices`` is True.

    Warning
    -------
    This will just compute the radial coordinates of objects WITHOUT the effect of
    gravitational lensing, no light rays will be computed. 
    """

    # Preparatory steps, checking inputs and generating grids if necessary
    if not isinstance( cosmo, c.model ) :
        raise TypeError( "Attribute cosmo should be an instance of class scampy.cosmology.model" )
    if funcs is None :
        zz = numpy.arange(0.0, 20.0, 1.e-2)
        dC = cosmo.dC(zz)
        z2d = lint( zz, dC )
        d2z = lint( dC, zz )
    else :
        try :
            z2d, d2z = funcs
        except :
            raise
        if not hasattr( z2d, "__call__" ) :
            raise TypeError( "First element of argument funcs should be callable" )
        if not hasattr( d2z, "__call__" ) :
            raise TypeError( "Second element of argument funcs should be callable" )

    # Find comoving distance of box centre
    dbox = z2d( zbox )
    # Find viewing angle (if not provided as input)
    if theta is None :
        zmin, zmax  = d2z(dbox+numpy.array([-0.5,0.5])*Lbox)
        theta = ang_size( Lbox, zmax, cosmo )
    
    # Centre the box to the observers LOS
    x, y, z = centre_box( box, Lbox ).T
    
    # Get redshift
    z = d2z( dbox + z )

    # Apply eventual redshift cut
    view = numpy.ones_like( z, dtype = bool)
    if zmin is not None :
        view &= zmin < z
    if zmax is not None :
        view &= z <= zmax
    
    # Get associated comoving distance
    dC = cosmo.dC( z[view] )
    
    # Convert to angles
    x[view] /= dC
    y[view] /= dC

    # Select objects within FOV
    view[view] = ( abs(x[view]) <= 0.5*theta ) & ( abs(y[view]) <= 0.5*theta )

    # Build cone
    cone = numpy.array( [ x[view], y[view], z[view] ] ).T

    # Return
    if return_indices :
        return cone, view
    return cone

##############################################################################

def collate_lightcones ( zlist, zmin, zmax, Lbox,
                         initshift = 0.0, cosmo = None ) :
    """Computes automathically the tuple of quantities necessary for
    building a lightcone by collating coordinates from different 
    snapshots of a cosmological N-body simulation.
    By shifting the position of the observer it is possible to 
    obtain different realizations.
    
    Parameters
    ----------
    zlist : iterable
        list of redshifts associated to each box
    zmin : scalar float
        minimum redshift of final lightcone
    zmax : scalar float
        maximum redshift of final lightcone
    Lbox : scalar float
        size of the periodic box
    initshift : scalar or iterable
        (Optional, default = 0.0) either a scalar float or an iterable
        that can be unpacked in 3 values.
        In case of a float, it will be replicated along the 3 axis,
        shifting the periodic box by initshift Mpc/h in each dimension.
        In case iterable, it should strictly have length = 3, the first
        two elements are the shifts associated with the perpendicular
        direction while the 3rd is the shift applied along the LoS direction
    cosmo : scampy.cosmology.model
        (Optional, default = None) Instance of an object of type scampy.cosmology.model.
        if none is passed, it will be built with default cosmological parameters
        (not an optimal solution)
    
    Returns
    -------
    zcen : 1d-array
        list of the central redshifts for each box (this differs from the input
        zlist iterable only by the first and last value, that are computed in order
        to have zmin < zcen[0] and zcen[-1] < zmax )
    dcen : 1d-array
        comoving distance of the observer from the centre of each box
    shifts : 2d-array
        shifts that are to be applied to each direction of each box
        in order to preserve continuity of the simulation (shape=(Nbox,3))
    zlim : 2d-array
        redshift limits of each box composing the lightcone (shape=(Nbox,2)). 
        For each box j, these should be zlim[j,0] < zcen[j] < zlim[j,1]
    """

    #####################################################################################
    # Preparatory steps, checking inputs and generating grids if necessary
    
    if hasattr( initshift, "__len__" ) :
        try :
            x1shift, x2shift, x3shift = numpy.array( initshift ) % Lbox
        except ValueError :
            raise ValueError(
                "initshift argument received wrong number of dimensions, "
                "requires either 0 or 3 dimensions"
            )
    else :
        try :
            x1shift, x2shift, x3shift = float( initshift ) % Lbox * numpy.ones( shape=(3,) )
        except ValueError :
            raise ValueError(
                "could not cast provided initshift to float"
            )
            
    # cosmology build or check
    if cosmo is None :
        cosmo = c.model()
    elif not isinstance( cosmo, c.model ) :
        raise TypeError( "Attribute cosmo should be an instance of class scampy.cosmology.model" )
            
    #####################################################################################
    # Here we build stuff

    # Compute central redshifts
    zcen = numpy.array( zlist )
    zcen[ 0 ]  = 0.5 * ( zcen[ 1 ] + zmin  )
    zcen[ -1 ] = 0.5 * ( zmax + zcen[ -2 ] )
    
    # Get number of boxes
    Nbox, = zcen.shape
    
    # Compute central distances
    dcen = cosmo.dC( zcen )

    # Compute box shifts
    # 1. allocate zero-padded array
    shift = numpy.zeros_like( zcen )
    # 2. first shift is user-defined
    shift[ 0 ] = x3shift
    # 3. all the others are given by the centres
    #    (comoving) cumulative distances from each other
    #    guaranteeing periodicity
    shift[ 1: ] = dcen[1:] - dcen[:-1]
    shift = numpy.cumsum(shift)
    shift %= Lbox
    # 4. build 2D shift matrix shape=(Nbox,3)
    shifts = numpy.array( [ x1shift * numpy.ones_like(zcen),
                            x2shift * numpy.ones_like(zcen),
                            -shift ] ).T

    # Compute limits of redshift slices
    # 1. allocate empty 2D array
    zlim = numpy.empty( shape = (Nbox, 2), dtype = float )
    # 2. set the minimum and maximum limits to the
    #    required zmin and zmax (user-defined)
    zlim[0,0]  = zmin
    zlim[-1,1] = zmax
    # 3. compute the upper limits == mid-distance
    #    between two box centres
    zlim[:-1,1] = 0.5 * ( zcen[:-1] + zcen[1:] )
    # 4. deduce lower limits
    zlim[1:,0] = zlim[:-1,1]
            
    #####################################################################################
    # Final check and return

    # Error if any of the boxes has wrongly defined redshift slice limits
    if numpy.any( [ z1 < z0 for (z0, z1) in zlim ] ) :
        raise RuntimeError(
            "something went wrong while computing the redshift limits" )

    return zcen, dcen, shifts, zlim

##############################################################################
