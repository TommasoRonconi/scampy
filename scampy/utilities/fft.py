"""Fast Hankel transforms on log-spaced grids.

Provides wrappers around :func:`scipy.fft.fht` for computing the
3D Fourier transform pair :math:`P(k) \\leftrightarrow \\xi(r)` and
the projected (Abel) transform pair
:math:`P(k) \\leftrightarrow w(r_p)`, both on log-spaced grids.
"""

#############################################################################################

# External imports
from scipy import fft
import numpy

#############################################################################################

def fstl ( xn, yn, lk0 = 0.0, bias = 0.0, mu = 0.5 ) :
    """Fast Spherical Transform on a log-spaced grid (direct).

    Evaluates the 3D two-point correlation function from the power
    spectrum via a Fast Hankel Transform (FHT) of order
    :math:`\\mu = 1/2`, implementing Eq. 6 of Ronconi et al. (2020):

    .. math::

        \\xi(r) = \\frac{1}{2\\pi^2}
        \\int_0^\\infty k^2\\,P(k)\\,\\frac{\\sin(kr)}{kr}\\,\\mathrm{d}k.

    The input grid ``xn`` must be log-spaced; the output :math:`r`-grid
    is the conjugate grid returned by :func:`scipy.fft.fht`.

    Parameters
    ----------
    xn : array-like
        Log-spaced wavenumber grid :math:`k` in
        :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
    yn : array-like
        Power spectrum :math:`P(k)` values, shape ``(xn.size,)`` or
        ``(xn.size, Nz)`` for multiple redshifts.
    lk0 : float, optional
        Initial log-offset passed to :func:`scipy.fft.fhtoffset`
        (default: ``0.0``).
    bias : float, optional
        Power-law bias for the FHT (default: ``0.0``).
    mu : float, optional
        Order of the Hankel transform (default: ``0.5``, corresponding
        to the spherical Bessel function :math:`j_0`).

    Returns
    -------
    kn : ndarray
        Output :math:`r`-grid in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
    xi : ndarray
        Correlation function :math:`\\xi(r)`, shape ``(xn.size,)`` or
        ``(Nz, xn.size)``.
    """

    # copy input and prep for transform
    lxn = numpy.log( xn )
    yn = numpy.array( yn ).T    
    # log and copy input
    lxn = numpy.log( xn )
    
    # Find log-spacing and check equal
    dlr = None
    try :
        # Note that the `numpy.around` function is less
        # accurate than the equivalent python built-in 
        # function `round` 
        # (rounds to the closest machine precision float)
        dlr = numpy.unique( 
            numpy.around( lxn[1:]-lxn[:-1], 4 ) 
        ).item()
    except :
        raise
        
    # Prepare input function for transform
    fct = 0.5 * numpy.sqrt( 0.5 * numpy.pi ) / numpy.pi**2
    an  = fct * yn * xn * numpy.sqrt( xn )
        
    # Find optimal centre for transform
    kr = fft.fhtoffset( dlr, initial = lk0, 
                        mu = mu, bias = bias )
    
    # Compute wavenumber grid
    kn = numpy.exp( kr - lxn[::-1] )
    
    # Compute sin transform of an
    bn = fft.fht( an, dlr, mu = mu, bias = bias, offset = kr )
    
    # Return a tuple with
    # - wavenumber grid
    # - converted sin-transform of input function
    return kn, ( bn / ( kn * numpy.sqrt( kn ) ) ).T

#############################################################################################

def fpstl ( xn, yn, lk0 = 0.0, bias = 0.0, mu = 0.0 ) :
    """Fast Projected Spherical Transform on a log-spaced grid (direct).

    Evaluates the projected two-point correlation function from the
    power spectrum via a zeroth-order Fast Hankel Transform (FHT),
    implementing the Abel transform of Eq. 15 of Ronconi et al. (2020):

    .. math::

        w(r_p) = \\frac{1}{2\\pi}
        \\int_0^\\infty k\\,P(k)\\,J_0(k\\,r_p)\\,\\mathrm{d}k,

    where :math:`J_0` is the zeroth-order Bessel function of the first
    kind.

    Parameters
    ----------
    xn : array-like
        Log-spaced wavenumber grid :math:`k` in
        :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
    yn : array-like
        Power spectrum :math:`P(k)` values, shape ``(xn.size,)`` or
        ``(xn.size, Nz)`` for multiple redshifts.
    lk0 : float, optional
        Initial log-offset passed to :func:`scipy.fft.fhtoffset`
        (default: ``0.0``).
    bias : float, optional
        Power-law bias for the FHT (default: ``0.0``).
    mu : float, optional
        Order of the Hankel transform (default: ``0.0``, corresponding
        to :math:`J_0`).

    Returns
    -------
    kn : ndarray
        Output :math:`r_p`-grid in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
    wp : ndarray
        Projected correlation function :math:`w(r_p)`, shape
        ``(xn.size,)`` or ``(Nz, xn.size)``.
    """
    
    # copy input and prep for transform
    lxn = numpy.log( xn )
    yn = numpy.array( yn ).T
    
    # Find log-spacing and check equal
    dlr = None
    try :
        # Note that the `numpy.around` function is less
        # accurate than the equivalent python built-in 
        # function `round` 
        # (rounds to the closest machine precision float)
        dlr = numpy.unique( 
            numpy.around( lxn[1:]-lxn[:-1], 4 ) 
        ).item()
    except :
        raise
        
    # Prepare input function for transform
    fct = 0.5  / numpy.pi
    an  = numpy.array( fct * yn * xn )
    
    # Find optimal centre for transform
    kr = fft.fhtoffset( dlr, initial = lk0, 
                        mu = mu, bias = bias )
    
    # Compute wavenumber grid
    kn = numpy.exp( kr - lxn[::-1] )
    
    # Compute sin transform of an
    bn = fft.fht( an, dlr, mu = mu, bias = bias, offset = kr )
    
    # Return a tuple with
    # - wavenumber grid
    # - converted sin-transform of input function
    return kn, ( bn / kn ).T

#############################################################################################
