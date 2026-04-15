"""Halo bias functions.

Provides fitting functions for the large-scale linear halo bias
:math:`b(M_h, z)`, defined such that the halo–matter power spectrum
equals :math:`b^2(M_h)\\,P_\\mathrm{lin}(k,z)` on large scales.
"""

#############################################################################################

# External imports
import numpy

# Internal imports
from scampy.power_spectrum import power_spectrum

#############################################################################################

def ShethMoTormen01 ( mm, zz, pk, comoving = True ) :
    """Sheth, Mo & Tormen (2001) halo bias.

    Implements the fitting formula

    .. math::

        b(M_h, z) = 1 + \\frac{1}{\\delta_c}
        \\left[
            a\\nu^2 + \\frac{a\\nu^2}{(a\\nu^2)^{0.6}} -
            \\frac{(a\\nu^2)^{0.6}}{(a\\nu^2)^{0.6} + 0.14}
        \\right],\\quad a = 0.707,

    where :math:`\\nu = \\delta_c(z)\\,/\\,[D(z)\\,\\sigma(M_h)]`
    is the peak-height parameter and :math:`\\delta_c(z)` is the linear
    collapse threshold.

    Parameters
    ----------
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s) at which to evaluate the bias.
        The output is broadcast over ``(mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object carrying the cosmological model and growth factor.
    comoving : bool, optional
        If ``True`` (default) masses are in comoving units (currently unused,
        reserved for future physical-unit support).

    Returns
    -------
    ndarray
        Array of shape ``(mm.size,)`` or ``(mm.size, zz.size)`` containing
        :math:`b(M_h, z)`.
    """

    mm = numpy.asarray(mm)
    if mm.ndim == 0 :
        mm = mm[None]
    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
    # temporaly not included
    # h_1 = 1.0
    # if not comoving :
    #     h_1 /= self.h0
  
    dd = pk.cosmo.deltac( zz )
    nu2 = ( dd / pk.D( zz ) )**2 / pk.sigma2M( mm, 0.0 )[:,numpy.newaxis]
    an2 = 0.707 * nu2
    an2p = an2**0.6
    
    return 1. + 1. / dd * numpy.squeeze( 
        an2 + 0.5 * an2 / an2p -
        an2p / ( 0.84084292 * ( an2p + 0.14 ) ) 
    )

#############################################################################################

def Tinker10 ( mm, zz, pk, Delta = 200, comoving = True ) :
    """Tinker et al. (2010) halo bias.

    Implements the fitting formula

    .. math::

        b(M_h, z) = 1
        - A\\,\\frac{\\nu^a}{\\nu^a + \\delta_c^a}
        + B\\,\\nu^b
        + C\\,\\nu^c,

    with coefficients depending on the overdensity threshold
    :math:`\\Delta` via :math:`y = \\log_{10}\\Delta`:

    .. math::

        A = 1 + 0.24\\,y\\,e^{-(4/y)^4},\\quad
        a = 0.44\\,y - 0.88,\\quad
        B = 0.183,\\quad b = 1.5,\\quad
        C = 0.019 + 0.107\\,y + 0.19\\,e^{-(4/y)^4},\\quad c = 2.4.

    Parameters
    ----------
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s) at which to evaluate the bias.
        The output is broadcast over ``(mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object carrying the cosmological model and growth factor.
    Delta : float, optional
        Overdensity threshold with respect to the mean matter density
        (default: 200).
    comoving : bool, optional
        If ``True`` (default) masses are in comoving units (currently unused,
        reserved for future physical-unit support).

    Returns
    -------
    ndarray
        Array of shape ``(mm.size,)`` or ``(mm.size, zz.size)`` containing
        :math:`b(M_h, z)`.
    """

    mm = numpy.asarray(mm)
    if mm.ndim == 0 :
        mm = mm[None]
    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
    # temporaly not included
    # h_1 = 1.0
    # if not comoving :
    #     h_1 /= self.h0
    
    dd = pk.cosmo.deltac( zz )
    nu = dd / pk.D( zz ) / numpy.sqrt( pk.sigma2M( mm, 0.0 ) )[:,numpy.newaxis]
    yy = numpy.log10( Delta )
    exy = numpy.exp( -( 4/yy )**4 )
    AA = 1 + 0.24 * yy * exy
    aa = 0.44 * yy - 0.88
    BB = 0.183
    bb = 1.5
    CC = 1.9e-2 + 0.107 * yy + 0.19 * exy
    cc = 2.4
    
    return numpy.squeeze( 1 - 
        AA * nu**aa / ( nu**aa + dd**aa ) +
        BB * nu**bb + 
        CC * nu**cc
    )

#############################################################################################
