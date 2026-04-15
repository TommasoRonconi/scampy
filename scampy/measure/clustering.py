"""Two-point correlation function estimators.

Provides standard and Landy–Szalay estimators for the two-point
correlation function :math:`\\xi(r)` (Eq. 23 of Ronconi et al. 2020),
together with a bootstrap error estimator.  Pair counting is delegated
to the compiled C++ extension :mod:`scampy.measure.clustering_core`.
"""

##################################################################################
# External imports
import numpy

#Internal imports
import scampy.measure.clustering_core as cc
from scampy.boxes import equatorial_to_cartesian, angular_to_euclidean_dist

##################################################################################

def _kernel_DD (data, Nd, rbins, omp = True#, angular = False
                ) :
    """Count data–data pairs in each separation bin."""
    
    if omp :
        if Nd == 2 :
            # if angular :
            #     return numpy.array( cc.dA2D_DD_omp( *data, rbins ) )
            return numpy.array( cc.d2D_DD_omp( *data, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DD_omp( *data, rbins ) )
    else :
        if Nd == 2 :
            # if angular :
            #     return numpy.array( cc.dA2D_DD( *data, rbins ) )
            return numpy.array( cc.d2D_DD( *data, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DD( *data, rbins ) )
            
    return None

def _kernel_DR (data1, data2, Nd, rbins, omp = True#, angular = False
                ) :
    """Count data–random cross-pairs in each separation bin."""
    
    if omp :
        if Nd == 2 :
            # if angular :
            #     return numpy.array( cc.dA2D_DR_omp( *data1, *data2, rbins ) )
            return numpy.array( cc.d2D_DR_omp( *data1, *data2, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DR_omp( *data1, *data2, rbins ) )
    else :
        if Nd == 2 :
            # if angular :
            #     return numpy.array( cc.dA2D_DR( *data1, *data2, rbins ) )
            return numpy.array( cc.d2D_DR( *data1, *data2, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DR( *data1, *data2, rbins ) )
            
    return None

##################################################################################

def _kernel_standard ( DD, RR ) :
    """Apply the standard estimator: :math:`\\xi = DD/RR - 1`."""

    if DD.size != RR.size :
        raise ValueError( 'distance vectors have different sizes' )

    ww = RR > 0
    out = numpy.zeros_like( ww, dtype = float )
    out[ww] = DD[ww] / RR[ww] - 1

    return out

##################################################################################

def two_point_standard ( data, rand, rbins, omp = True, angular = False
                        ) :
    """Two-point correlation function with the standard estimator.

    Computes :math:`\\xi(r) = DD/RR - 1`, where :math:`DD` and
    :math:`RR` are the normalised data–data and random–random pair
    counts in each separation bin.

    Parameters
    ----------
    data : ndarray, shape ``(Ndim, Nobj)``
        Coordinates of the data catalogue.  ``Ndim`` must be 2 or 3.
    rand : ndarray, shape ``(Ndim, Nrand)``
        Coordinates of the random catalogue, same dimensionality as
        ``data``.
    rbins : array-like
        Bin edges for the separation :math:`r`.
    omp : bool, optional
        Use the OpenMP-parallel pair counter (default: ``True``).
    angular : bool, optional
        Reserved for future angular-space counting (currently unused).

    Returns
    -------
    ndarray
        Correlation function :math:`\\xi(r)` in each bin, shape
        ``(rbins.size - 1,)``.
    """

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    if not NobjD > 0 or not NobjR > 0 :
        raise ValueError(
            "Cannot compute clustering if one of the two catalogues does not have at least 2 elements"
        )

    normDD = 2.0 / ( NobjD * ( NobjD - 1 ) )
    normRR = 2.0 / ( NobjR * ( NobjR - 1 ) )

    # DD = _kernel_DD( data, NdimD, rbins, omp, angular ) * normDD
    DD = _kernel_DD( data, NdimD, rbins, omp ) * normDD
    # RR = _kernel_DD( rand, NdimD, rbins, omp, angular ) * normRR
    RR = _kernel_DD( rand, NdimD, rbins, omp ) * normRR

    return _kernel_standard( DD, RR )
    
##################################################################################

def _kernel_landy_szalay ( DD, RR, DR ) :
    """Apply the Landy–Szalay estimator: :math:`\\xi = (DD - 2DR)/RR + 1`."""

    if DD.size != RR.size or RR.size != DR.size :
        raise ValueError( 'distance vectors have different sizes' )

    ww = RR > 0
    out = numpy.zeros_like( ww, dtype = float )
    out[ww] = ( DD[ww] - 2.0 * DR[ww] ) / RR[ww] + 1

    return out
    
##################################################################################

def _kernel_error_landy_szalay () :
    """Analytical error on the Landy–Szalay estimator (not yet implemented)."""
    pass
    
##################################################################################

def two_point_landyszalay ( data, rand, rbins, omp = True, return_error = False, angular = False ) :
    """Two-point correlation function with the Landy–Szalay estimator.

    Implements Eq. 23 of Ronconi et al. (2020):

    .. math::

        \\xi(r) = \\frac{DD(r) - 2\\,DR(r) + RR(r)}{RR(r)},

    where :math:`DD`, :math:`DR`, and :math:`RR` are the normalised
    data–data, data–random, and random–random pair counts.  This is the
    preferred estimator used in the paper's validation tests.

    Parameters
    ----------
    data : ndarray, shape ``(Ndim, Nobj)``
        Coordinates of the data catalogue.  ``Ndim`` must be 2 or 3.
    rand : ndarray, shape ``(Ndim, Nrand)``
        Coordinates of the random catalogue.
    rbins : array-like
        Bin edges for the separation :math:`r`.
    omp : bool, optional
        Use the OpenMP-parallel pair counter (default: ``True``).
    return_error : bool, optional
        If ``True``, also return a per-bin error estimate derived from
        the Poisson fluctuations of the pair counts (default: ``False``).
    angular : bool, optional
        Reserved for future angular-space counting (currently unused).

    Returns
    -------
    xi : ndarray
        Correlation function :math:`\\xi(r)`, shape ``(rbins.size - 1,)``.
    exi : ndarray
        Per-bin error estimate, only returned when ``return_error=True``.
    """

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    if not NobjD > 0 or not NobjR > 0 :
        raise ValueError(
            "Cannot compute clustering if one of the two catalogues does not have at least 2 elements"
        )

    normDD = 2.0 / ( NobjD * ( NobjD - 1 ) )
    normRR = 2.0 / ( NobjR * ( NobjR - 1 ) )
    normDR = 1.0 / ( NobjD * NobjR )

    # DD = _kernel_DD( data, NdimD, rbins, omp, angular )
    DD = _kernel_DD( data, NdimD, rbins, omp )
    DDn = DD * normDD
    # RR = _kernel_DD( rand, NdimD, rbins, omp, angular )
    RR = _kernel_DD( rand, NdimD, rbins, omp )
    RRn = RR * normRR
    # DR = _kernel_DR( data, rand, NdimD, rbins, omp, angular )
    DR = _kernel_DR( data, rand, NdimD, rbins, omp )
    DRn = DR * normDR

    # compute baseline clustering
    xi = _kernel_landy_szalay( DDn, RRn, DRn )
    
    if return_error :
        fact = (
            normRR / normDD * RR * (1.0 + xi ) +
            4.0 / NobjD * ( normRR * RR / normDD * ( 1.0 + xi ) )**2
        )
        exi = numpy.zeros_like( xi )
        ww = RRn > 0
        exi[ww] = normDD / RRn[ww] * numpy.sqrt( fact[ww] ) * numpy.sqrt( 3 )
        return xi, exi
    return xi

##################################################################################

def bootstrap_two_point ( data, rand, rbins,
                          standard = True,
                          Nboots = 10, return_boots = False,
                          omp = True, verbose = True, angular = False,
                          rng = None, kw_rng = { 'seed' : 555 } ) :
    """Two-point correlation function with bootstrap error estimate.

    Computes a baseline :math:`\\xi(r)` and estimates its uncertainty by
    resampling the data catalogue with replacement ``Nboots`` times.  The
    random catalogue is fixed across all resamples (only :math:`RR` is
    computed once).  Either the standard or the Landy–Szalay estimator
    can be used for each resample.

    Parameters
    ----------
    data : ndarray, shape ``(Ndim, Nobj)``
        Coordinates of the data catalogue.  ``Ndim`` must be 2 or 3,
        or 2 when ``angular=True``.
    rand : ndarray, shape ``(Ndim, Nrand)``
        Coordinates of the random catalogue.
    rbins : array-like
        Bin edges for the separation :math:`r` (or angular separation
        in radians when ``angular=True``).
    standard : bool, optional
        If ``True`` (default) use the standard estimator; otherwise use
        the Landy–Szalay estimator.
    Nboots : int, optional
        Number of bootstrap resamples (default: ``10``).
    return_boots : bool, optional
        If ``True``, also return the full array of bootstrap realisations
        (default: ``False``).
    omp : bool, optional
        Use the OpenMP-parallel pair counter (default: ``True``).
    verbose : bool, optional
        Print progress messages (default: ``True``).
    angular : bool, optional
        If ``True``, interpret ``data`` and ``rand`` as equatorial
        coordinates ``(ra, dec)`` in radians, convert to 3D Cartesian
        unit vectors, and convert ``rbins`` from angular to chord
        distances (default: ``False``).
    rng : numpy.random.Generator or None, optional
        Random number generator.  If ``None`` (default) one is created
        via ``numpy.random.default_rng(**kw_rng)``.
    kw_rng : dict, optional
        Keyword arguments forwarded to ``numpy.random.default_rng``
        (default: ``{'seed': 555}``).

    Returns
    -------
    xi : ndarray
        Baseline correlation function, shape ``(rbins.size - 1,)``.
    sigma : ndarray
        Bootstrap standard deviation per bin, same shape as ``xi``.
    boots : ndarray, shape ``(Nboots, rbins.size - 1)``
        Individual bootstrap realisations; only returned when
        ``return_boots=True``.
    """

    if rng is None :
        rng = numpy.random.default_rng( **kw_rng )
    data = numpy.array( data )
    rand = numpy.array( rand )

    if angular :
        rbins = angular_to_euclidean_dist( rbins, r = 1.0 )
        if data.shape[0] != 2 or rand.shape[0] != 2 :
            raise RuntimeError(
                'Input spaces ``data`` and ``rand`` should be 2D when ``angular = True`` is chosen'
            )
        data = numpy.array(
            equatorial_to_cartesian(
                ( numpy.ones_like( data[ 0 ] ), data[ 0 ], data[ 1 ] )
            )
        )
        rand = numpy.array(
            equatorial_to_cartesian(
                ( numpy.ones_like( rand[ 0 ] ), rand[ 0 ], rand[ 1 ] )
            )
        )

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    if not NobjD > 0 or not NobjR > 0 :
        raise ValueError(
            "Cannot compute clustering if one of the two catalogues does not have at least 2 elements"
        )

    normDD = 2.0 / ( NobjD * ( NobjD - 1 ) )
    normRR = 2.0 / ( NobjR * ( NobjR - 1 ) )
    normDR = 1.0 / ( NobjD * NobjR )
    
    if verbose : print( 'Computing RR ...' )
    # RR = _kernel_DD( rand, NdimR, rbins, omp, angular ) * normRR
    RR = _kernel_DD( rand, NdimR, rbins, omp ) * normRR
    if verbose : print( '... done RR.' )

    # get the baseline estimate
    if verbose : print( 'Computing baseline ...' )
    # DD = _kernel_DD( data, NdimD, rbins, omp, angular ) * normDD
    DD = _kernel_DD( data, NdimD, rbins, omp ) * normDD
    if standard :
        tpt = _kernel_standard( DD, RR )
    else :
        # DR = _kernel_DR( data, rand, NdimD, rbins, omp, angular ) * normDR
        DR = _kernel_DR( data, rand, NdimD, rbins, omp ) * normDR
        tpt = _kernel_landy_szalay( DD, RR, DR )
    if verbose : print( '... done baseline.' )

    # get bootstraps
    boots = numpy.zeros( ( Nboots, tpt.size ) )
    for ii in range( Nboots ):
        if verbose : print( f'Computing bootstrap {ii+1:d} of {Nboots} ...', end='\r' )
        iD = rng.integers( NobjD, size = NobjD )
        # DD = _kernel_DD( data[:,iD], NdimD, rbins, omp, angular ) * normDD
        DD = _kernel_DD( data[:,iD], NdimD, rbins, omp ) * normDD
        if standard :
            boots[ii] = _kernel_standard( DD, RR )
        else :
            # DR = _kernel_DR( data[:,iD], rand, NdimD, rbins, omp, angular ) * normDR
            DR = _kernel_DR( data[:,iD], rand, NdimD, rbins, omp ) * normDR
            boots[ii] = _kernel_landy_szalay( DD, RR, DR )
    if verbose : print( '\n... done bootstraps.' )

    if return_boots:
        return tpt, boots.std(axis=0), boots
    else:
        return tpt, boots.std(axis=0)

