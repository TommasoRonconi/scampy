""" 
"""

#############################################################################################

# External imports
import numpy

# Internal imports
from scampy.power_spectrum import power_spectrum
from scampy.halo.mass_function import Tinker08
from scampy.halo.bias import Tinker10
from scampy.halo.density_profile import density_profile_FT
from scampy.hod import HOD
from scampy.utilities.functions import trap_int
from scampy.utilities.fft import fstl
from scampy.utilities.interpolation import lin_interp

#############################################################################################

class halo_model () :
    
    def __init__ ( self, pk, comoving = True, Mlim = ( 5, 17 ), 
                   klim = ( -3, +3 ), thin = 256,
                   # these strings are just placeholders for now
                   hmf = 'Tinker08', hbias = 'Tinker10', 
                   density = 'NFW' 
                 ) :
        if not isinstance( pk, power_spectrum ) :
            raise ValueError( 
                'argument pk should be an instance of type scampy.power_spectrum' 
            )
        self.pk = pk
        self.hmf = Tinker08
        self.hbs = Tinker10
        self.hdp = density_profile_FT
        self.Mh_min, self.Mh_max = Mlim
        self.k_min, self.k_max = klim
        self.Mh_grid = numpy.logspace( *Mlim, thin )
        self.kh_grid = numpy.logspace( *klim, thin ) 
        
    def ng ( self, hod, zz = 0.0 ) :
        return trap_int( 
            self.Mh_grid, 
            ( hod.Pcen( self.Mh_grid ) + 
              hod.Psat( self.Mh_grid ) ) *
            self.hmf( self.Mh_grid, zz, self.pk )
        )
    
    def bias ( self, hod, zz = 0.0 ) :
        return trap_int(
            self.Mh_grid,
            ( hod.Pcen( self.Mh_grid ) + 
              hod.Psat( self.Mh_grid ) ) *
            self.hmf( self.Mh_grid, zz, self.pk ) * 
            self.hbs( self.Mh_grid, zz, self.pk )
        ) / self.ng( hod, zz )
    
    def Mhalo ( self, hod, zz = 0.0 ) :
        return trap_int(
            self.Mh_grid,
            ( hod.Pcen( self.Mh_grid ) + 
              hod.Psat( self.Mh_grid ) ) *
            self.hmf( self.Mh_grid, zz, self.pk ) * 
            self.Mh_grid
        ) / self.ng( hod, zz )
    
    def dngdM ( self, mm, hod, zz = 0.0 ) :
        return ( hod.Pcen( mm ) + hod.Psat( mm ) ) * self.hmf( mm, zz, self.pk )
    
    ####################################################################################
    # Power spectrum
    
    def Pk_1halo ( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        if fact > 0.0 :
            pass
        else :
            fact = 1. / self.ng( hod, zz )**2
        return fact * trap_int(
            self.Mh_grid[:,numpy.newaxis],
            ( ( 2 * hod.Pcen( self.Mh_grid ) + 
                hod.Psat( self.Mh_grid ) ) * hod.Psat( self.Mh_grid ) *
              self.hmf( self.Mh_grid, zz, self.pk )
            )[:,numpy.newaxis] *
            self.hdp( kk, self.Mh_grid, zz, self.pk )**2
        )
    
    def Pk_2halo( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        if fact > 0.0 :
            pass
        else :
            fact = 1. / self.ng( hod, zz )**2
        return fact * self.pk.Pz( kk, zz ) * trap_int(
            self.Mh_grid[:,numpy.newaxis],
            ( ( hod.Pcen( self.Mh_grid ) + 
                hod.Psat( self.Mh_grid ) ) *
              self.hmf( self.Mh_grid, zz, self.pk ) *
              self.hbs( self.Mh_grid, zz, self.pk )
            )[:,numpy.newaxis] *
            self.hdp( kk, self.Mh_grid, zz, self.pk )
        )**2
    
    def Pk ( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        if fact > 0.0 :
            pass
        else :
            fact = 1. / self.ng( hod, zz )**2
        return (
            self.Pk_1halo( kk, hod, zz, fact ) + 
            self.Pk_2halo( kk, hod, zz, fact )
        )
    
    ####################################################################################
    # 3D-space 2-point correlation function
    
    def Xi_1halo ( self, rr, hod, zz = 0.0, fact = -1.0 ) :
        rn, xi1h = fstl(
            self.kh_grid,
            self.Pk_1halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.5
        )
        return lin_interp( rn, xi1h )( rr )
    
    def Xi_2halo ( self, rr, hod, zz = 0.0, fact = -1.0 ) :
        rn, xi2h = fstl(
            self.kh_grid,
            self.Pk_2halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.5
        )
        return lin_interp( rn, xi2h )( rr )
    
    def Xi ( self, rr, hod, zz = 0.0, fact = -1.0 ) :
        rn, xi1h = fstl(
            self.kh_grid,
            self.Pk_1halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.5
        )
        _,  xi2h = fstl(
            self.kh_grid,
            self.Pk_2halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.5
        )
        return lin_interp( rn, xi1h + xi2h )( rr )

#############################################################################################
