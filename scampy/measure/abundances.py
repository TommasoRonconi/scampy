"""Functions computing abundance statistics of different quantities that can be extracted 
from catalogues.
"""

import numpy
from scipy.stats import binned_statistic

############################################################################################

def get_abundances ( Ax, Nc, Ns, bins, statistic = 'mean' ) :
    """ Compute the mean number of central and satellite objects
    within bins along a given axis.
    
    Parameters:
    -----------
    Ax : array-like
        the axis along which we want to compute 
        the average values of Ncen and Nsat
    Nc, Ns : array-like
        the value of Ncen and Nsat for all the
        values in Ax
    bins : array-like
        edges of the bins
    statistic : string or callable
        (Optional, default = ``'mean'``) 
        argument statistic to be passed to 
        function ``scipy.stats.binned_statistic`` 
        used to compute the abundances
    
        
    Returns:
    --------
    : array-like
        bin-centre (computed linearly)
    <Nc>, <Ns> : array-like
        tuple with the chosen statistic of Ncen and Nsat 
        (respectively, default is mean), 
        within the provided bins.
        It has the same shape of bins
    """

    [ Nc_binned, Ns_binned ], binCntr, *_ = binned_statistic( Ax, values = ( Nc, Ns ),
                                                              statistic = statistic,
                                                              bins = bins )
    Nc_binned[ numpy.isnan( Nc_binned ) ] = 0.
    Ns_binned[ numpy.isnan( Ns_binned ) ] = 0.
    binCntr = 0.5 * ( binCntr[:-1] + binCntr[1:] )
    
    return binCntr, Nc_binned, Ns_binned

############################################################################################

def differential_counts ( var, bins, fact = 1. ) :
    """Obtain the differential distribution of a variable :math:`v`.
    Computes :math:`\\frac{dN}{dx} \\cdot const`, where :math:`dN` is the number
    of :math:`v` values in the bin, :math:`dv` it the size of the bin and :math:`const`
    is a constant value.

    Parameters
    ----------
    var : array-like
       the list of :math:`X` values
    bins : array-like
       it defines the bin edges, including the rightmost edge, 
       allowing for non-uniform bin widths.
    fact : scalar
       the value of :math:`const`, constant factor by which to multiply the result

    Returns
    -------
    array-like
       bin center (shape = len(bins) - 1)
    array-like
       sequence of values with the values of 
       :math:`\\frac{dN}{dv} \\cdot const` (shape = len(bins) - 1)
    array-like
       Poisson errors on the differential counts, defined as
       :math:`\\frac{\\sqrt{dN}}{dv} \\cdot const` (shape = len(bins) - 1)    
    """
    
    counts, bins = numpy.histogram( var, bins = bins )

    v = 0.5 * ( bins[ :-1 ] + bins[ 1: ] )
    dv = ( bins[ 1: ] - bins[ :-1 ] )
    dndv = fact * counts / dv
    dndv_er = fact * numpy.sqrt( counts ) / dv
    
    return v, dndv, dndv_er

############################################################################################

def cumulative_counts ( X, bins, fact = 1. ) :
    """Obtain the cumulative distribution of a variable :math:`X`.
    Computes :math:`N(X>x_i) \\cdot const`, where :math:`N(X>x_i)` is the number
    of :math:`X` values greater than :math:`x_i` and :math:`const`
    is a constant value.

    Parameters
    ----------
    var : array-like
       the list of :math:`X` values
    binning : array-like
       defines the binned values :math:`x_i`, allowing for non-uniform bin widths.
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
    
    Nx = numpy.zeros_like(bins)
    Nx[:-1], *_ = numpy.histogram( X, bins )
    Nx = numpy.flip(numpy.flip(Nx).cumsum())
    
    return Nx*fact, numpy.sqrt(Nx)*fact

############################################################################################

def cumulative_from_differential ( func, bins, fact = 1. ) :
    """Obtain the cumulative distribution from a differential distribution :math:`f(X)`.
    Computes 

    .. math:: F(X > x) \\equiv const \\cdot \\int_{-\\infty}^x f(X) dX

    where :math:`F(X>x)` is the cumulative distribution and :math:`const` is a constant value.

    Parameters
    ----------
    func : function
       the differential function :math:`f(X)`, it should take only one argument (a.k.a. the 
       value of :math:`X`). It allows :code:`lambda` expressions.
    bins : array-like
       defines the list :math:`x_i` values, allowing for non-uniform bin widths.
    fact : scalar
       the value of :math:`const`, constant factor by which to multiply the result

    Returns
    -------
    array-like
       sequence of values with the result of :math:`F(X>x_i)` (shape = bins.shape)
    """
    
    dx = (bins[1:] - bins[:-1])
    dN = func(bins)
    cN = numpy.zeros_like(bins)
    cN[1:] = ( 0.5 * (dN[1:]+dN[:-1]) * dx ).cumsum() * fact
        
    return cN

############################################################################################
