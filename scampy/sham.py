##########################################################################################

# External imports
import numpy
import warnings

# Internal imports
from scampy.utilities.functions import mod_erf, truncated_gaussian
from scampy.utilities.interpolation import lin_interp as lint
from scampy.measure.abundances import cumulative_counts, cumulative_from_differential

##########################################################################################

def scatter_array ( arr, method = 'normal', log = False,
                    rng = None, kw_rng = { 'seed' : 555 }, **kwargs  ) :
    """Given an input array adds scattering to the input values.
    Different methods are available:
    
    * 'normal' : Gaussian scatter around the input array

    Parameters
    ----------
    arr : iterable
    
    method : str
       (Optional, default = ``'normal'``)

    log : bool
       (Optional, default = ``False``)
    
    rng : random number generator
        (Optional, default = ``None``) if None is passed, it will build one using
        function ``numpy.random.default_rng(**kw_rng)``
    
    kw_rng : dict
        (Optional, default = { 'seed' : 555 }) keyword arguments to pass to
        the default random number generator. (if rng is not None, this input is ignored) 
    
    Keyword Arguments
    -----------------
    Depending on the method selected
        
    * method = ``'normal'``:
       scale : float or iterable
          (Optional, default = ``1.0``) scattering magnitude (in dex). The default value is ``0.0``
       loc : float or iterable
          if not zero, adds a bias to the relation. The default value is ``0.0``.

    Returns
    -------
    : ndarray
    a copy of the original array with scatter added
    """

    arr = numpy.array( arr )
    if log :
        arr = numpy.log10( arr )
    if rng is None :
        rng = numpy.random.default_rng(**kw_rng)
        
    # Add scatter
    if method == 'normal' :
        arr += rng.normal(size = arr.size, **kwargs)
        #arr *= 1.0+rng.normal(size = arr.size, **kwargs)
        
    if method == 'normal_decaying' :
        scale = (1-mod_erf( arr, kwargs.get( 'factor', 1.0 ) )) * kwargs.get( 'scale', 0.0 )
        arr *= 1.0+rng.normal(size = arr.size, scale = scale)
        
    if method == 'uniform' :
        scale = kwargs.get( 'scale', 0.0 )
        arr *= 1.0+scale*rng.uniform(low = -1., high = 1., size = arr.size)
        
    if method == 'uniform_decaying' :
        scale = (1-mod_erf( arr, kwargs.get( 'factor', 1.0 ) )) * kwargs.get( 'scale', 0.0 )
        arr *= 1.0+scale*rng.uniform(low = -1., high = 1., size = arr.size)
        
    if method == 'rank_ordering' :
        from warnings import warn
        warn('`rank_ordering` method not implemented yet, passing without adding any scatter')

    if log :
        return 10**arr
    return arr

##########################################################################################

def matching_sample ( InSize, OutSize,
                      InProb = None,
                      rng = None, kw_rng = { 'seed' : 555 } ) :
    """Returns a boolean 1D mask of given size, with given number of occurrencies unmasked.
    
    Parameters
    ----------
    InSize : int
        Required mask size
    
    OutSize : int
        Required number of un-masked occurencies in output mask
    
    InProb : iterable
        (Optional, default = ``None``) An iterable of length = ``InSize`` with 
        probability associated to each occurrency in the output mask.
        If not given, the sample assumes a uniform distribution over all
        entries in the output mask.
    
    rng : random number generator
        (Optional, default = ``None``) if None is passed, it will build one using
        function ``numpy.random.default_rng(**kw_rng)``
    
    kw_rng : dict
        (Optional, default = { 'seed' : 555 }) keyword arguments to pass to
        the default random number generator. (if rng is not None, this input is ignored) 
    
    Returns
    -------
    : 1d-array
    boolean mask with size = ``InSize`` and summing to ``OutSize``.
    """
    
    if rng is None :
        rng = numpy.random.default_rng(**kw_rng)

    mask = numpy.zeros( InSize, dtype = bool )
    try : 
        sample = rng.choice( range(InSize), size=OutSize,
                             replace=False, p=InProb, shuffle=False )
    except :
        raise
    mask[sample] = True
    
    return mask

##########################################################################################

def match_distribution ( var, func, minF = -20, maxF = -10,
                         minP = 1.e-20, maxP = 1., nbinF = 200, factF = 1.,
                         minV = None, maxV = None, nbinV = 100, factV = 1., 
                         scatter_type = None, kw_scatter = {},
                         rng = None, kw_rng = { 'seed' : 555 } ) :
    """Implements the Abundance Matching algorithm.
    A value in the distribution ``func`` is associated to each element of the array ``var``.
    Note that, the underline assumption is that all the elements of ``var`` can be linked
    to an element in the distribution ``func``.
    This will not perform an HOD nor it will consider the possibility of not having
    counter parts in the ``func`` distribution.

    Parameters
    ----------
    var : iterable
        Input variables. Their cumulative distribution is computed in 
        ``nbinV`` bins between ``minV`` and ``maxV``. This cumulative
        abundance is then matched to the cumulative abundance of ``func``.    
    func : function
        Differential distribution (i.e. density function) of the properties
        that have to be associated to variables in the input iterable ``var``.
        The corresponding cumulative distribution is computed in ``nbinF`` bins
        between ``minF`` and ``maxF`` using function ``cumulative_from_differential``
        of module ``scampy.measure.abundances``
    minF : scalar

    maxF : scalar

    minP : scalar
    
    maxP : scalar
    
    nbinF : int
    
    factF : scalar
    
    minV : scalar or None
    
    maxV : scalar or None
    
    nbinV : int
    
    factV : scalar
    
    scatter_type : str
        (Optional, default = ``None``) Type of scatter to add to the output array of properties.
        Currently implemented:

        * ``'normal'``: Gaussian scatter around the predicted value

        * ``None`` : none of the above, no scatter is added

    kw_scatter : dict
        (Optional, default = ``{}``) Dictionary with the keyword arguments to pass to the function
        computing the scatter. 
        
        * method = ``'normal'``:
          'scale' = scattering magnitude (in dex). The default value is ``1.0``
          'loc' = if not zero, adds a bias to the relation. The default value is ``0.0``
    
    rng : random number generator
        (Optional, default = ``None``) if None is passed, it will build one using
        function ``numpy.random.default_rng(**kw_rng)``
    
    kw_rng : dict
        (Optional, default = { 'seed' : 555 }) keyword arguments to pass to
        the default random number generator. (if rng is not None, this input is ignored) 

    Returns
    -------
    : ndarray
        array of properties sampled from func and matching in abundance the values in var 
    
    Note
    ----
    The output reproducibility (i.e. the possibility to get the same result running 
    the function twice on the same inputs) depends on whether a scattering has been 
    added and, in that case, is guaranteed only if the random number generator passed
    is in the same state. 
    """

    
    # some checks and set-ups:
    if nbinF <= 0 :
        raise ValueError( 'we will need at least one bin here, '
                          f'you gave {nbinF:d}' )
    var = numpy.array( var )
    if minV is None :
        minV = var.min()
        minV -= 0.01 * minV
    if maxV is None :
        maxV = var.max()
        maxV += 0.01 * maxV
    print(f'Vmin = {minV:.3e}, Vmax = {maxV:.3e}')
        
    if rng is None :
        rng = numpy.random.default_rng( **kw_rng )
    
    # Get cumulative function
    fbin = numpy.linspace( minF, maxF, nbinF )
    Pf = cumulative_from_differential( func, fbin, factF )

    # Get limits for inversion
    wPf = ( minP <= Pf ) & ( Pf <= maxP )
    print( f'There are {wPf.sum():d} valid points '
           'function in the probability grid' )

    # Get inverse dPhi 
    f_invPf = lint( Pf[ wPf ], fbin[ wPf ] )

    # Get cumulative sub-halo mass function
    lminV = numpy.log10( minV )
    lmaxV = numpy.log10( maxV )
    lvbin = numpy.logspace( lminV, lmaxV, nbinV )
    Pv, _ = cumulative_counts( var, lvbin, factV )
    
    if ( (Pv.min() < numpy.min(f_invPf.get_x())) | 
         (numpy.max(f_invPf.get_x()) < Pv.max()) ) :
        warnings.warn( 'Some of the values are outside the '
                       'modelled probability interval, '
                       'the model will be extrapolated.' )
    
    # Direct function inversion matching
    f_Pv = lint( lvbin, Pv )
    out = f_invPf( f_Pv( var ) )

    if scatter_type is not None :
        out = scatter_array( out, method = scatter_type, rng = rng, **kw_scatter )
    
    return out

##########################################################################################

def match_cross ( X, Y,
                  Xlim = None, Xnbin = 100, factX = 1.0, logX = False,
                  Ylim = None, Ynbin = 100, factY = 1.0, logY = False,
                  Xscatter = None, kw_Xscatter = {},
                  Yscatter = None, kw_Yscatter = {},
                  rng = None, kw_rng = { 'seed' : 555 } ) :
    """Find indices of the elements in the X array that match in abundance with elements
    of the Y array. Returns the two list of indices in bijective relation. 

    Parameters
    ----------
    X : 1darray

    Y : 1darray

    Xlim : tuple

    Xnbin : int

    factX : float

    logX : bool

    Ylim : tuple

    Ynbin : int

    factY : float

    logY : bool
    
    Xscatter : str
        (Optional, default = ``None``) Type of scatter to add to the output array of properties.
        Currently implemented:

        * ``'normal'``: Gaussian scatter around the predicted value

        * ``None`` : none of the above, no scatter is added to the X-array

    kw_Xscatter : dict
        (Optional, default = ``{}``) Dictionary with the keyword arguments to pass to the function
        computing the scatter. 
        
        * method = ``'normal'``:
          'scale' = scattering magnitude (in dex). The default value is ``1.0``
          'loc' = if not zero, adds a bias to the relation. The default value is ``0.0``
    
    Yscatter : str
        (Optional, default = ``None``) same as ``Xscatter`` for the ``Y`` array

    kw_Yscatter : dict
        (Optional, default = ``{}``) same as ``kw_Xscatter`` for the ``Y`` array

    rng : random number generator
        (Optional, default = ``None``) if None is passed, it will build one using
        function ``numpy.random.default_rng(**kw_rng)``
    
    kw_rng : dict
        (Optional, default = { 'seed' : 555 }) keyword arguments to pass to
        the default random number generator. (if rng is not None, this input is ignored) 

    Returns
    -------
    idx_sort : int ndarray
    idy_sort : int ndarray
    
    Note
    ----
    The output reproducibility (i.e. the possibility to get the same result running 
    the function twice on the same inputs) depends on whether a scattering has been 
    added and, in that case, is guaranteed only if the random number generator passed
    is in the same state on all runs. 
    """

    if rng is None :
        rng = numpy.random.default_rng(**kw_rng)

    # Copy the 2 arrays
    X = numpy.array( X )
    Y = numpy.array( Y )

    # manage input X limits
    if Xlim is None :
        Xmin = X.min(); Xmin -= 0.01 * Xmin
        Xmax = X.max(); Xmax += 0.01 * Xmax
    else :
        try :
            Xmin, Xmax = Xlim
        except :
            raise
    print(f'Xmin = {Xmin:.3e}, Xmax = {Xmax:.3e}')

    # Generate bins in X-space:
    if logX :
        xbin = numpy.logspace( numpy.log10(Xmin), numpy.log10(Xmax), Xnbin )
    else :
        xbin = numpy.linspace( Xmin, Xmax, Xnbin )

    # manage input Y limits
    if Ylim is None :
        Ymin = Y.min(); Ymin -= 0.01 * Ymin
        Ymax = Y.max(); Ymax += 0.01 * Ymax
    else :
        try :
            Ymin, Ymax = Ylim
        except :
            raise
    print(f'Ymin = {Ymin:.3e}, Ymax = {Ymax:.3e}')

    # Generate bins in Y-space:
    if logY :
        ybin = numpy.logspace( numpy.log10(Ymin), numpy.log10(Ymax), Ynbin )
    else :
        ybin = numpy.linspace( Ymin, Ymax, Ynbin )

    # Compute the cumulative abundances of the two input arrays:
    Nx, _ = cumulative_counts( X, xbin, fact = factX )
    Ny, _ = cumulative_counts( Y, ybin, fact = factY )

    # Allocate interpolators for the abundances above
    f_Px = lint(xbin, Nx)
    f_Py = lint(ybin, Ny)

    # Find limits for existence of counterparts
    Pmin = numpy.max([Nx.min(), Ny.min()])
    Pmax = numpy.min([Nx.max(), Ny.max()])

    # Compute probability object-by-object
    Px = f_Px(X)
    Py = f_Py(Y)

    # Add scatter to input arrays
    if Xscatter is not None :
        Px = scatter_array( Px, method = Xscatter, log = logX, rng = rng, **kw_Xscatter ) 
    if Yscatter is not None :
        Py = scatter_array( Py, method = Yscatter, log = logY, rng = rng, **kw_Yscatter )

    # find mask for counterparts
    wx = (Pmin <= Px)&(Px <=Pmax)
    wy = (Pmin <= Py)&(Py <=Pmax)

    # match the two samples
    if wx.sum() > wy.sum() :
        wx[wx] &= matching_sample( wx.sum(), wy.sum(), Px[wx]/Px[wx].sum(), rng = rng )
    elif wx.sum() < wy.sum() :
        wy[wy] &= matching_sample( wy.sum(), wx.sum(), Py[wy]/Py[wy].sum(), rng = rng )
    else :
        pass

    # Indices sorted w.r.t. probability
    idx_sort = numpy.argsort(Px)
    idy_sort = numpy.argsort(Py)

    # 1-to-1 association:
    return idx_sort[wx[idx_sort]], idy_sort[wy[idy_sort]]

##########################################################################################
