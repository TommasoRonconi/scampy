""" 
"""

#############################################################################################

# External imports
import numpy

# Internal imports
import scampy.cosmology as c
from scampy.utilities.functions import trap_int
from scampy.utilities.interpolation import lin_interp
from scampy.utilities.functions import FT_tophat, FT_tophat_D1
from scampy.utilities.fft import fstl

#############################################################################################
# Classes

class power_spectrum () :
    
    def __init__ ( self,  kh, pk0,
                   kmin = 1.e-4, kmax = 100,
                   cosmo = None, thinness = 1000 ) :

        #####################################################################################
        # Checking inputs and generating objects if necessary
        
        self.kh = numpy.array( kh )
        self.kmin = kmin
        if kmin is None :
            self.kmin = self.kh.min()
        self.kmax = kmax
        if kmax is None :
            self.kmax = self.kh.max()
            
        # cosmology build or check
        if cosmo is None :
            self.cosmo = c.model()
        elif not isinstance( cosmo, c.model ) :
            raise TypeError(
                "Attribute cosmo should be an instance of class scampy.cosmology.model"
            )
        else :
            self.cosmo = cosmo

        #####################################################################################
        #
        
        self.pk0 = numpy.array( pk0 )
        self.P0 = lin_interp( self.kh, self.pk0 )
        self.kh = numpy.logspace( numpy.log10(self.kmin), numpy.log10(self.kmax), 1024 )
        self.pk0 = self.P0(self.kh)
        self.h0 = self.cosmo.param['hh']
        self.sigma8 = -1
        _ = self.compute_sigma8()
        self.zz = numpy.linspace( self.cosmo.zmin, self.cosmo.zmax, thinness )
        self.D = lin_interp( self.zz, self.cosmo.D(self.zz)/self.cosmo.D(0.0) )
        
    def compute_sigma8 ( self ) :
        if self.sigma8 < 0.0 :
            integrand = self.kh**2 * self.pk0 * FT_tophat( 8.0*self.kh )**2
            self.sigma8 = 0.5 * trap_int(self.kh, integrand ) / numpy.pi**2
        self.sigma8correction = self.cosmo.param['sigma8'] / self.sigma8
        return self.sigma8
    
    def Pz ( self, kk, zz = 0.0, comoving = True ) :
        """
        Returns
        -------
        : ndarray
            array with shape (kk.size, zz.size)
        """
        
        kk = numpy.asarray(kk)
        if kk.ndim == 0 :
            kk = kk[None]
        zz = numpy.asarray(zz)
        if zz.ndim == 0 :
            zz = zz[None]
        h_1 = 1.0
        h_3 = 1.0
        if not comoving :
            h_1 /= self.h0
            h_3 = h_1**3
            
        fact = h_3 * self.sigma8correction * self.cosmo.param['sigma8'] * self.D( zz )**2
        
        if numpy.asarray( fact ).size > 1 :
            return numpy.squeeze(fact * self.P0( kk * h_1 )[:,numpy.newaxis])
        return fact * self.P0( kk * h_1 )
    
    def sigma2R ( self, rr, zz, comoving = True ) :
        """
        Returns
        -------
        : ndarray
            array with shape (rr.size, zz.size)
        """
        
        rr = numpy.asarray(rr)
        if rr.ndim == 0 :
            rr = rr[None]
        zz = numpy.asarray(zz)
        if zz.ndim == 0 :
            zz = zz[None]
            
        integrand = (
            (self.kh**2) * 
            self.Pz( self.kh, zz, comoving = comoving ).T * 
            (FT_tophat( self.kh * rr[:,numpy.newaxis] )**2)[:,numpy.newaxis,:] 
        )
        
        return numpy.squeeze( 
            0.5 * trap_int( self.kh[:,numpy.newaxis, numpy.newaxis], 
                            numpy.transpose(integrand, axes=(2,0,1)) ) / 
            numpy.pi**2 
        )
    
    def dsigma2RdR ( self, rr, zz, comoving = True ) :
        """
        Returns
        -------
        : ndarray
            array with shape (rr.size, zz.size)
        """
        
        rr = numpy.asarray(rr)
        if rr.ndim == 0 :
            rr = rr[None]
        zz = numpy.asarray(zz)
        if zz.ndim == 0 :
            zz = zz[None]
            
        integrand = (
            (self.kh**3) * 
            self.Pz( self.kh, zz, comoving = comoving ).T * 
            ( FT_tophat( self.kh * rr[:,numpy.newaxis] ) *
              FT_tophat_D1( self.kh * rr[:,numpy.newaxis] ) )[:,numpy.newaxis,:] )
        
        return numpy.squeeze( 
            trap_int( self.kh[:,numpy.newaxis, numpy.newaxis], 
                      numpy.transpose(integrand, axes=(2,0,1)) ) / 
            numpy.pi**2 
        )        
    
    def sigma2M ( self, mm, zz, comoving = True ) :
        
        mm = numpy.asarray(mm)
        if mm.ndim == 0 :
            mm = mm[None]
        zz = numpy.asarray(zz)
        if zz.ndim == 0 :
            zz = zz[None]
        
        # at redshift z = 0.0 because 
        # the redshift dependence is in the sigma2R function:
        if comoving :
            rhoM = self.cosmo.critical_density_comoving(0.0) * self.cosmo.OmegaM(0.0)
        else :
            rhoM = self.cosmo.critical_density(0.0) * self.cosmo.OmegaM(0.0)
        
        # radius of the equivalent sphere
        rr = (0.75 * mm / ( rhoM * numpy.pi ))**(1/3)
        
        return self.sigma2R( rr, zz, comoving = comoving )
        
    def dsigma2MdM ( self, mm, zz, comoving = True ) :
        
        # at redshift z = 0.0 because 
        # the redshift dependence is in the sigma2R function:
        if comoving :
            rhoM = self.cosmo.critical_density_comoving(0.0) * self.cosmo.OmegaM(0.0)
        else :
            rhoM = self.cosmo.critical_density(0.0) * self.cosmo.OmegaM(0.0)
        
        # radius of the equivalent sphere
        rr = (0.75 * mm / ( rhoM * numpy.pi ))**(1/3)
        
        return numpy.squeeze(
            ( rr / ( 3. * mm ) ) * 
            self.dsigma2RdR( rr, zz, comoving = comoving ).T
        ).T

    def Xi ( self, rr, zz = 0.0, comoving = True ) :
        """Linear 2-point correlation function
        obtained by Fourier-Transforming the linear power spectrum
        """
        rn, xin = fstl( self.kh, self.Pz( self.kh, zz, comoving ) )
        return lin_interp( rn, xin )( rr )                        

#############################################################################################
