"""Halo mass functions.

Provides fitting functions for the halo mass function :math:`n(M_h, z)`,
i.e. the comoving number density of dark matter haloes per unit mass interval.
All functions share the same call signature and return
:math:`\\mathrm{d}n/\\mathrm{d}M_h` in units of
:math:`[M_\\odot^{-1}\\,h^3\\,\\mathrm{Mpc}^{-3}]` (comoving) or
:math:`[M_\\odot^{-1}\\,\\mathrm{Mpc}^{-3}]` (physical).
"""

#############################################################################################

# External imports
import numpy

# Internal imports
from scampy.power_spectrum import power_spectrum

#############################################################################################

def ShethTormen01 ( mm, zz, pk, comoving = True ) :
    """Sheth & Tormen (2001) halo mass function.

    Implements the fitting formula

    .. math::

        \\frac{\\mathrm{d}n}{\\mathrm{d}M} =
        \\rho_0 \\left|\\frac{\\mathrm{d}\\ln\\sigma}{\\mathrm{d}M}\\right|
        f(\\nu)\\,/\\,M

    where

    .. math::

        f(\\nu) = A\\sqrt{\\frac{2}{\pi}}\\,0.3222
        \\left(1 + (a\\nu^2)^{-0.3}\\right)\\nu\\,
        e^{-a\\nu^2/2},
        \\quad a = 0.707,\\quad A \\approx 0.8409,

    and :math:`\\nu = \\delta_c(z)\\,/\\,[D(z)\\,\\sigma(M)]` is the
    peak-height parameter.

    Parameters
    ----------
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s) at which to evaluate the mass function.
        The output is broadcast over ``(mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object carrying the cosmological model and the
        normalised growth factor.
    comoving : bool, optional
        If ``True`` (default) masses and densities are in comoving units.

    Returns
    -------
    ndarray
        Array of shape ``(mm.size,)`` or ``(mm.size, zz.size)`` containing
        :math:`\\mathrm{d}n/\\mathrm{d}M_h`.
    """

    mm = numpy.asarray(mm)
    if mm.ndim == 0 :
        mm = mm[None]
    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
    h_1 = 1.0
    if not comoving :
        h_1 /= self.h0
    
    s20 = pk.sigma2M( mm, 0.0 )
    nu = pk.cosmo.deltac( zz ) / pk.D( zz ) / numpy.sqrt( s20 )[:,numpy.newaxis]
    
    fnu = (
        0.84083292 * numpy.sqrt( 2 / numpy.pi ) * 0.3222 *
        ( 1 + ( 0.84083292 * nu )**( -0.6 ) ) * nu *
        numpy.exp( - 0.5 * 0.707 * nu**2 )
    )
    
    dls = -0.5 * pk.dsigma2MdM( mm, 0.0 ) / s20
    
    rho0 = pk.cosmo.param['Om_M'] * pk.cosmo.critical_density( 0.0 )
    
    return numpy.squeeze(
        rho0 * dls[:,numpy.newaxis] * fnu / mm[:,numpy.newaxis] * h_1
    )

#############################################################################################

def Tinker08 ( mm, zz, pk, comoving = True ) :
    """Tinker et al. (2008) halo mass function.

    Implements the fitting formula

    .. math::

        \\frac{\\mathrm{d}n}{\\mathrm{d}M} =
        \\rho_0 \\left|\\frac{\\mathrm{d}\\ln\\sigma}{\\mathrm{d}M}\\right|
        f(\\sigma)\\,/\\,M,

    where

    .. math::

        f(\\sigma) = A_z
        \\left[\\left(\\frac{\\sigma}{b_z}\\right)^{-a_z} + 1\\right]
        e^{-c_0/\\sigma^2}

    with redshift-dependent coefficients

    .. math::

        A_z = A_0\\,(1+z)^{-0.14},\\quad
        a_z = a_0\\,(1+z)^{-0.06},\\quad
        b_z = b_0\\,(1+z)^{-\\alpha},

    calibrated for a :math:`\\Delta = 200` overdensity threshold
    (:math:`A_0=0.186,\\,a_0=1.47,\\,b_0=2.57,\\,c_0=1.19`).

    Parameters
    ----------
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s) at which to evaluate the mass function.
        The output is broadcast over ``(mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object carrying the cosmological model and growth factor.
    comoving : bool, optional
        If ``True`` (default) masses and densities are in comoving units.

    Returns
    -------
    ndarray
        Array of shape ``(mm.size,)`` or ``(mm.size, zz.size)`` containing
        :math:`\\mathrm{d}n/\\mathrm{d}M_h`.
    """

    mm = numpy.asarray(mm)
    if mm.ndim == 0 :
        mm = mm[None]
    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
    h_1 = 1.0
    if not comoving :
        h_1 /= self.h0
        
    mass  = mm * h_1
    gfact = pk.D( zz )
    rho0  = pk.cosmo.param['Om_M'] * pk.cosmo.critical_density_comoving( 0.0 )
    rr = (0.75 * mm / ( rho0 * numpy.pi ))**(1/3)
    s2 = pk.sigma2R( rr, 0.0, comoving = True )
    ss = numpy.sqrt( s2 )
    dls = (
        - 0.5 * rr *
        pk.dsigma2RdR( rr, 0.0, comoving = True )
        / ( 3 * s2 * mass )
    )

    nu = gfact * ss[:,numpy.newaxis]
    A0 = 0.186; a0 = 1.47; b0 = 2.57; c0 = 1.19
    alpha = 1.0676e-2 # = 10.^(- (0.75/log10( Delta/75 ))**1.2 ) ) with Delta=200
    Az = A0 * ( 1 + zz )**( -0.14 )
    az = a0 * ( 1 + zz )**( -0.06 )
    bz = b0 * ( 1 + zz )**( -alpha )

    fnu = (
        Az * ( ( nu/bz )**( -az ) + 1 ) *
        numpy.exp( -c0 / nu**2 )
    )

    return numpy.squeeze( rho0 * dls[:,numpy.newaxis] * fnu / mass[:,numpy.newaxis] )

#############################################################################################

def Behroozi13 ( mm, zz, pk, comoving = True ) :
    """Behroozi, Wechsler & Conroy (2013) halo mass function.

    Applies a high-redshift correction on top of :func:`Tinker08`:

    .. math::

        n_\\text{B13}(M,z) =
        10^{N(a)}\\cdot\\left(3.16\\times10^{-12}\\,M\\right)^{p(a)}
        \\cdot n_\\text{T08}(M,z),

    where :math:`a = 1/(1+z)` is the scale factor,

    .. math::

        N(a) = \\frac{0.144}{1 + e^{14.79\\,(a - 0.213)}},\\quad
        p(a) = \\frac{0.5}{1 + e^{6.5\\,a}}.

    The correction is negligible at low redshift and boosts the abundance
    of low-mass haloes at high redshift.

    Parameters
    ----------
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s) at which to evaluate the mass function.
        The output is broadcast over ``(mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object carrying the cosmological model and growth factor.
    comoving : bool, optional
        If ``True`` (default) masses and densities are in comoving units.

    Returns
    -------
    ndarray
        Array of shape ``(mm.size,)`` or ``(mm.size, zz.size)`` containing
        :math:`\\mathrm{d}n/\\mathrm{d}M_h`.
    """

    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
    
    ascale = 1. / ( 1 + zz )
    logNorm = 0.144 / ( 1 + numpy.exp( 14.79 * ( ascale - 0.213 ) ) )
    high_z_correction = ( 3.162278e-12 * mm[:,numpy.newaxis] )**(
        0.5 / ( 1 + numpy.exp( 6.5 * ascale ) )
    )

    return numpy.squeeze(
        10.**( logNorm * high_z_correction )
    ) * Tinker08( mm, zz, pk, comoving )

#############################################################################################
