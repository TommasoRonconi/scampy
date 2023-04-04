""" 
"""

#############################################################################################

# External imports
import numpy

# Internal imports
from scampy.power_spectrum import power_spectrum

#############################################################################################

def ShethMoTormen01 ( mm, zz, pk, comoving = True ) :
  
    dd = pk.cosmo.deltac( zz )
    gfact = 1.
    nu2 = dd**2 / pk.sigma2M( mm, 0.0 ) / pk.D( zz )**2 
    an2 = 0.707 * nu2
    an2p = an2**0.6
    
    return 1. + 1. / dd * ( 
        an2 + 0.5 * an2 / an2p -
        an2p / ( 0.84084292 * ( an2p + 0.14 ) ) 
    )

#############################################################################################

def Tinker10 ( mm, zz, pk, Delta = 200, comoving = True ) :
    
    dd = pk.cosmo.deltac( zz )
    nu = dd / numpy.sqrt( pk.sigma2M( mm, 0.0 ) ) / pk.D( zz )
    yy = numpy.log10( Delta )
    exy = numpy.exp( -( 4/yy )**4 )
    AA = 1 + 0.24 * yy * exy
    aa = 0.44 * yy - 0.88
    BB = 0.183
    bb = 1.5
    CC = 1.9e-2 + 0.107 * yy + 0.19 * exy
    cc = 2.4
    
    return ( 1 - 
        AA * nu**aa / ( nu**aa + dd**aa ) +
        BB * nu**bb + 
        CC * nu**cc
    )

#############################################################################################
