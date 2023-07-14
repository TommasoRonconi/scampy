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
from scampy.utilities.functions import trap_int, linear_interpolation
from scampy.utilities.fft import fstl, fpstl
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
        zz = numpy.asarray( zz )
        mgrid = self.Mh_grid
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
        return trap_int(
            mgrid,
            ( hod.Pcen( mgrid ) + 
              hod.Psat( mgrid ) ) *
            self.hmf( self.Mh_grid, zz, self.pk )
        )
    
    def bias ( self, hod, zz = 0.0 ) :
        zz = numpy.asarray( zz )
        mgrid = self.Mh_grid
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
        return trap_int(
            mgrid,
            ( hod.Pcen( mgrid ) + 
              hod.Psat( mgrid ) ) *
            self.hmf( self.Mh_grid, zz, self.pk ) * 
            self.hbs( self.Mh_grid, zz, self.pk )
        ) / self.ng( hod, zz )
    
    def Mhalo ( self, hod, zz = 0.0 ) :
        zz = numpy.asarray( zz )
        mgrid = self.Mh_grid
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
        return trap_int(
            mgrid,
            ( hod.Pcen( mgrid ) + hod.Psat( mgrid ) ) *
            self.hmf( self.Mh_grid, zz, self.pk ) * mgrid
        ) / self.ng( hod, zz )
    
    def dngdM ( self, mm, hod, zz = 0.0 ) :
        zz = numpy.asarray( zz )
        mgrid = mm
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
        return ( hod.Pcen( mgrid ) + hod.Psat( mgrid ) ) * self.hmf( mm, zz, self.pk )
    
    ####################################################################################
    # Power spectrum
    
    def Pk_1halo ( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        zz = numpy.asarray( zz )
        mgrid = self.Mh_grid
        intshape = [1,0]
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
            intshape += [2]
        if numpy.all( fact > 0.0 ) :
            pass
        else :
            fact = 1. / self.ng( hod, zz )**2
        return fact * trap_int(
            mgrid[:,numpy.newaxis],
            numpy.transpose(
                ( ( 2 * hod.Pcen( mgrid ) + 
                    hod.Psat( mgrid ) ) * hod.Psat( mgrid ) *
                  self.hmf( self.Mh_grid, zz, self.pk )
                )[numpy.newaxis,:] *
                self.hdp( kk, self.Mh_grid, zz, self.pk )**2,
                axes=intshape
            )
        )
    
    def Pk_2halo( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        zz = numpy.asarray( zz )
        mgrid = self.Mh_grid
        intshape = [1,0]
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
            intshape += [2]
        if numpy.all( fact > 0.0 ) :
            pass
        else :
            fact = 1. / self.ng( hod, zz )**2
        return fact * self.pk.Pz( kk, zz ) * trap_int(
            mgrid[:,numpy.newaxis],
            numpy.transpose(
                ( hod.Pcen( mgrid ) + 
                  hod.Psat( mgrid ) ) *
                self.hmf( self.Mh_grid, zz, self.pk ) *
                self.hbs( self.Mh_grid, zz, self.pk ) *
                self.hdp( kk, self.Mh_grid, zz, self.pk ),
                axes=intshape
            )
        )**2
    
    def Pk ( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        if numpy.all( fact > 0.0 ) :
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
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rn, xi1h = fstl(
            self.kh_grid,
            self.Pk_1halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.5
        )
        return lin_interp( rn, xi1h )( rr )
    
    def Xi_2halo ( self, rr, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rn, xi2h = fstl(
            self.kh_grid,
            self.Pk_2halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.5
        )
        return lin_interp( rn, xi2h )( rr )
    
    def Xi ( self, rr, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
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
    
    ####################################################################################
    # Projected 2-point correlation function

    def Wr_1halo ( self, rp, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rn, wr1h = fpstl(
            self.kh_grid,
            self.Pk_1halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.0
        )
        return lin_interp( rn, wr1h )( rp )

    def Wr_2halo ( self, rp, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rn, wr2h = fpstl(
            self.kh_grid,
            self.Pk_2halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.0
        )
        return lin_interp( rn, wr2h )( rp )
    
    def Wr ( self, rp, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rn, wr1h = fpstl(
            self.kh_grid,
            self.Pk_1halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.0
        )
        _,  wr2h = fpstl(
            self.kh_grid,
            self.Pk_2halo( self.kh_grid, hod, zz, fact ),
            lk0 = 0.0, bias = 0.0, mu = 0.0
        )
        return lin_interp( rn, wr1h + wr2h )( rp )
    
    ####################################################################################
    # Angular 2-point correlation function

    def Wt_1halo ( self, th, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rp = th * self.pk.cosmo.dC( zz )
        return self.Wr_1halo( rp, hod, zz, fact )

    def Wt_2halo ( self, th, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rp = th * self.pk.cosmo.dC( zz )
        return self.Wr_2halo( rp, hod, zz, fact )

    def Wt ( self, th, hod, zz = 0.0, fact = -1.0 ) :
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rp = th * self.pk.cosmo.dC( zz )
        return self.Wr( rp, hod, zz, fact )

    def Wt_zdist ( self, th, hod, zbins, Nzdist, fact = -1, valid_rp = (1.e-2, 8.e+1) ) :
        """
        Parameters
        ----------
        theta : 1d-array
        hod : hod instance
        zbins : 1d-array
        Nzdist : 1d-array
            normalized redshift distribution of sources (# in bin/# tot).
            Nzdist.size == zbins.size - 1
        fact : float or 1d-array
        valid_rp : tuple
            interval of projected radii valid for interpolation of the
            angular 2pt at fixed redshift (change at own risk)
        
        Returns
        -------
        """
    
        # input stuff
        th = numpy.array( th )
        zmed = 0.5 * ( zbins[1:] + zbins[:-1] )
        fact = numpy.array( fact )
        if zbins.size - 1 != Nzdist.size :
            raise RuntimeError(
                'argument ``zbins`` should have one element more than argument ``Nzdist``'
            )
    
        # compute power spectrum
        pks = self.Pk(self.kh_grid, hod, zmed, fact = fact)
        
        # fourier-transform to projected
        rint, wtint = fpstl(
            self.kh_grid, pks,
            lk0 = 0.0, bias = 0.0, mu = 0.0
        )
        vr = (valid_rp[0] < rint) & (rint < valid_rp[1])
        
        # derive theta array for fourier-transformed x-space
        thint = rint[:,numpy.newaxis] / self.pk.cosmo.dC(zmed)
        
        # compute angular 2-point on input array
        wt = numpy.array( [ 
            numpy.exp( 
                linear_interpolation(
                    numpy.log( th ),
                    numpy.log( _t[vr] ),
                    numpy.log( _w[vr] ) ) 
            ) 
            for _t, _w in zip( thint.T, wtint.T ) ] 
        ).T
        
        # derivative of the comoving distance
        ddC = self.pk.cosmo.ddC( zbins )
        
        # compute approximated averaged 2pt angular
        return numpy.sum( ( 
            ( zbins[1:] - zbins[:-1] ) * 
            Nzdist**2 / 
            numpy.abs( ddC[1:] - ddC[:-1] )
        ) * wt, axis = 1 )
    
#############################################################################################
