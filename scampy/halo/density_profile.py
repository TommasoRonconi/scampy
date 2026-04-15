"""NFW halo density profile and halo concentration models.

Provides the Fourier transform of the NFW density profile and several
empirical fitting functions for the halo concentration–mass relation
:math:`c(M_h, z)`.
"""

#############################################################################################

# External imports
import numpy
from scipy.special import sici

# Internal imports
from scampy.power_spectrum import power_spectrum

#############################################################################################

# Eq. 7, 11, 13 of Giocoli+2012 
# def z_form ( mm, zz, pk, ff = 0.04, z_max = 100.0, thin = 100 ) :
        
#     alpha_f = 0.815 * numpy.exp( -2. * ff**3 ) * ff**( -0.707 )
#     omega_f = numpy.sqrt( 2.*numpy.log( 1. + alpha_f ) )
#     delta_sigma = numpy.array(
#         numpy.sqrt( 
#             pk.sigma2M( mm * ff, 0.0 ) - pk.sigma2M( mm, 0.0 ) 
#         )
#     )
#     if delta_sigma.ndim == 0 :
#         delta_sigma = delta_sigma[None]
#     dc0 = pk.cosmo.deltac( 0.0 )

#     zspace = numpy.linspace( zz, z_max, thin )
#     print(dc0 * ( 1 / pk.D( zspace ) - 1 ))
#     fspace = numpy.squeeze( 
#         dc0 * ( 1 / pk.D( zspace ) - 1 ) - 
#         omega_f * delta_sigma[:,numpy.newaxis]
#     )
#     return zspace, fspace.T

#############################################################################################

def concentration_giocoli12 ( mm, zz, pk, comoving = True ) :
    """Giocoli et al. (2012) halo concentration–mass relation.

    .. note::
        Not yet implemented; returns ``None``.
    """
    pass

#############################################################################################

def concentration_zhao09 ( mm, zz, pk, comoving = True ) :
    """Zhao et al. (2009) halo concentration–mass relation.

    .. note::
        Not yet implemented; returns ``None``.
    """
    pass

#############################################################################################

def concentration_shimizu03 ( mm, zz, pk, comoving = True ) :
    """Shimizu et al. (2003) halo concentration–mass relation.

    Implements Eq. 4 of Shimizu et al. (2003):

    .. math::

        c(M_h, z) = \\frac{8}{1+z}\\,
        \\left(1.0204\\times10^{-14}\\,M_h\\,h\\right)^{-0.13}.

    Parameters
    ----------
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s).  The output is broadcast over ``(mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object; used to retrieve :math:`h` when
        ``comoving=False``.
    comoving : bool, optional
        If ``True`` (default) masses are already in
        :math:`M_\\odot\\,h^{-1}` and no :math:`h`-conversion is applied.

    Returns
    -------
    ndarray
        Array of shape ``(mm.size,)`` or ``(mm.size, zz.size)`` containing
        the dimensionless concentration :math:`c = r_\\mathrm{vir}/r_s`.
    """
        
    mm = numpy.asarray(mm)
    if mm.ndim == 0 :
        mm = mm[None]
    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
    hh = 1.
    if not comoving :
        hh = pk.h0
        
    return numpy.squeeze(
        8.0 / ( 1 + zz ) *
        ( ( 1.0204e-14 * mm * hh )**(-0.13) )[:,numpy.newaxis]
    )

#############################################################################################

def density_profile_FT ( kk, mm, zz, pk, comoving = True ) :
    """Fourier transform of the NFW density profile.

    Returns the normalised Fourier transform :math:`\\tilde{u}(k,z|M_h)`
    of the NFW profile for a halo of mass :math:`M_h`, as given by
    Eqs. 8–9 of Ronconi et al. (2020):

    .. math::

        \\tilde{u}(k, z|M_h) =
        \\frac{
            \\cos(\\mu)\\,[\\mathrm{Ci}(\\mu(1+c)) - \\mathrm{Ci}(\\mu)]
            + \\sin(\\mu)\\,[\\mathrm{Si}(\\mu(1+c)) - \\mathrm{Si}(\\mu)]
            - \\dfrac{\\sin(c\\mu)}{\\mu(1+c)}
        }{\\ln(1+c) - c/(1+c)},

    where :math:`\\mu = k\\,r_s`, :math:`c = c(M_h, z)` is the halo
    concentration, :math:`r_s = r_\\mathrm{vir}/c` is the NFW scale radius,
    and :math:`\\mathrm{Si}`, :math:`\\mathrm{Ci}` are the sine and cosine
    integrals.  The concentration is computed via
    :func:`concentration_shimizu03`.

    Parameters
    ----------
    kk : scalar or array-like
        Wavenumbers in :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
    mm : scalar or array-like
        Halo masses :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.
    zz : scalar or array-like
        Redshift(s).  The output is broadcast over ``(kk, mm, zz)``.
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object; used to retrieve the cosmological model
        (critical density, :math:`\\Delta_c`).
    comoving : bool, optional
        If ``True`` (default) the virial volume uses the comoving critical
        density; otherwise the physical critical density is used.

    Returns
    -------
    ndarray
        Normalised profile transform :math:`\\tilde{u}(k,z|M_h)`,
        dimensionless and in the range :math:`[0, 1]`.
    """
    
    kk = numpy.asarray(kk)
    if kk.ndim == 0 :
        kk = kk[None]
    mm = numpy.asarray(mm)
    if mm.ndim == 0 :
        mm = mm[None]
    zz = numpy.asarray(zz)
    if zz.ndim == 0 :
        zz = zz[None]
        
    conc = concentration_shimizu03( mm, zz, pk, comoving )
    if comoving :
        Vvir = numpy.squeeze(
            mm[:,numpy.newaxis] /
            ( pk.cosmo.Delta_c( zz ) * pk.cosmo.critical_density_comoving( zz ) )
        )
    else :
        Vvir = numpy.squeeze(
            mm[:,numpy.newaxis] /
            ( pk.cosmo.Delta_c( zz ) * pk.cosmo.critical_density( zz ) )
        )
    rvir = ( 0.75 * Vvir / numpy.pi )**(1/3)
    rs = rvir / conc
    mu = numpy.squeeze( kk[:,numpy.newaxis,numpy.newaxis] * rs )
    Si1, Ci1 = sici( mu )
    Si2, Ci2 = sici( mu + mu * conc )
    return numpy.squeeze( 
        ( 
            numpy.cos( mu ) * ( Ci2 - Ci1 ) + 
            numpy.sin( mu ) * ( Si2 - Si1 ) -
            numpy.sin( mu * conc ) / ( mu + mu * conc )
        ) /
        ( numpy.log( 1 + conc ) - conc / ( 1 + conc ) )
    )

#############################################################################################
