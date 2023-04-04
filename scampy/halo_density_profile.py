""" 
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
    pass

#############################################################################################

def concentration_zhao09 ( mm, zz, pk, comoving = True ) :
    pass

#############################################################################################

# Eq. 4 from Shimizu et al. 2003
def concentration_shimizu03 ( mm, zz, pk, comoving = True ) :
    hh = 1.
    if not comoving :
        hh = pk.h0
        
    return 8.0 / ( 1 + zz ) * ( 1.0204e-14 * mm * hh )**(-0.13)

#############################################################################################

# NFW Fourier Transform
def density_profile_FT ( kk, mm, zz, pk, comoving = True ) :
    
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
        Vvir = mm / ( pk.cosmo.Delta_c( zz ) * pk.cosmo.critical_density_comoving( zz ) )
    else :
        Vvir = mm / ( pk.cosmo.Delta_c( zz ) * pk.cosmo.critical_density( zz ) )
    rvir = ( 0.75 * Vvir / numpy.pi )**(1/3)
    rs = rvir / conc
    mu = kk * rs[:,numpy.newaxis]
    Si1, Ci1 = sici( mu )
    Si2, Ci2 = sici( mu + mu * conc[:,numpy.newaxis] )
    return (
        ( 
            numpy.cos( mu ) * ( Ci2 - Ci1 ) + 
            numpy.sin( mu ) * ( Si2 - Si1 ) -
            numpy.sin( mu * conc[:,numpy.newaxis] ) / ( mu + mu * conc[:,numpy.newaxis] )
        ) /
        ( numpy.log( 1 + conc ) - conc / ( 1 + conc ) )[:,numpy.newaxis]
    )

#############################################################################################
