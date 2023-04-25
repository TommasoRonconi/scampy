""" 
"""

#############################################################################################

# External imports
import numpy

# Internal imports
from scampy.power_spectrum import power_spectrum

#############################################################################################

def ShethTormen01 ( mm, zz, pk, comoving = True ) :
        
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
