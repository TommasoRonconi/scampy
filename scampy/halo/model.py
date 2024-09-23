""" 
"""

#############################################################################################

# External imports
import numpy
from scipy.special import erf

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
        self._lkhdamp = numpy.log10( self.pk.kh[ self.pk.pk0.argmax() ] )
        self._sigmalk = 0.25
        
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

    def _damp_1halo ( self, kk ) :
        return 0.5 * ( 1. + erf( ( numpy.log10(kk) - self._lkhdamp ) / self._sigmalk ) )
    
    def Pk_1halo ( self, kk, hod, zz = 0.0, fact = -1.0 ) :
        zz = numpy.asarray( zz )
        mgrid = self.Mh_grid
        intshape = [1,0]
        damp = self._damp_1halo(kk)
        if zz.size > 1 :
            mgrid = mgrid[:, numpy.newaxis]
            damp  = damp[:, numpy.newaxis]
            intshape += [2]
        if numpy.all( fact > 0.0 ) :
            pass
        else :
            fact = 1. / self.ng( hod, zz )**2
        return fact * damp * trap_int(
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

    def Nz ( self, hod, zbins, solid_angle = 4 * numpy.pi ) :
        """Normalised redshift distribution
        
        N(z) = Number of sources in z-bin / ( width of z-bin * total number of sources )
        
        Parameters
        ----------
        hod :
        zbins : array of floats
        solid_angle : float
        """

        dz = numpy.diff(zbins)
        zz = 0.5*(zbins[1:]+zbins[:-1])
        zhist = self.ng(hod, zz) * self.pk.cosmo.comoving_volume_unit( zz ) * solid_angle * dz
        
        return zhist / dz / zhist.sum()

    def Wt_zdist ( self, th, hod, zbins,
                   Nzdist = None, fact = -1,
                   valid_rp = (1.e-2, 8.e+1),
                   solid_angle = 4 * numpy.pi,
                   component = 'total' ) :
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
        solid_angle : float
            solid angle sub-tended the field of view 
        component : str
            (Optional, default = 'total') which component of the correlation function to combine
            Available options are: 
            'total' compute the sum of 1-halo and 2-halo term, 
            '1halo' compute the 1-halo term only, 
            '2halo' compute the 2-halo term only.
        
        Returns
        -------
        """
    
        # input stuff
        th = numpy.array( th )
        zmed = 0.5 * ( zbins[1:] + zbins[:-1] )
        fact = numpy.array( fact )
        if Nzdist is None :
            Nzdist = self.Nz( hod, zbins, solid_angle )
        if zbins.size - 1 != Nzdist.size :
            raise RuntimeError(
                'argument ``zbins`` should have one element more than argument ``Nzdist``'
            )
    
        # compute power spectrum
        if component == 'total' :
            pks = self.Pk(self.kh_grid, hod, zmed, fact = fact)
        elif component == '1halo' :
            pks = self.Pk_1halo(self.kh_grid, hod, zmed, fact = fact)
        elif component == '2halo' :
            pks = self.Pk_2halo(self.kh_grid, hod, zmed, fact = fact)
        else :
            raise RuntimeError( 'Chosen component invalid, valid options are "total", "1halo", "2halo"' )
        
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
                    numpy.log( _w[vr],
                               out = -30 * numpy.ones_like(_w[vr]),
                               where = (_w[vr]>0.0) )
                ) 
            ) 
            for _t, _w in zip( thint.T, wtint.T ) ] 
        ).T
        
        integrand = ( Nzdist**2 / self.pk.cosmo.ddC(zmed) ) * wt
        return trap_int( zmed[:,numpy.newaxis], integrand.T )
    
#############################################################################################
