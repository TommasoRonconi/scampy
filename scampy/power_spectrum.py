"""Linear matter power spectrum and related statistics.

Provides the :class:`power_spectrum` class, which wraps a tabulated
linear :math:`P(k)` at :math:`z=0`, normalises it to the cosmological
:math:`\\sigma_8`, evolves it to arbitrary redshift via the linear growth
factor, and exposes derived quantities used by the halo model: the
variance :math:`\\sigma^2(R,z)`, its mass-keyed counterpart
:math:`\\sigma^2(M,z)`, and the real-space two-point correlation function
:math:`\\xi(r,z)`.
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
    """Container for the linear matter power spectrum and its derived statistics.

    On construction the input tabulated :math:`P(k)` is interpolated onto
    a uniform log-spaced grid, normalised so that
    :math:`\\sigma_8` matches the value stored in the cosmological model,
    and the linear growth factor :math:`D(z)/D(0)` is pre-computed on a
    redshift grid for fast evaluation.

    Parameters
    ----------
    kh : array-like
        Wavenumbers of the input power spectrum in
        :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
    pk0 : array-like
        Linear power spectrum at :math:`z=0` in
        :math:`[h^{-3}\\,\\mathrm{Mpc}^3]`, evaluated at the wavenumbers
        ``kh``.
    kmin : float, optional
        Minimum wavenumber of the internal interpolation grid
        (default: ``1e-4``).
    kmax : float, optional
        Maximum wavenumber of the internal interpolation grid
        (default: ``100``).
    cosmo : scampy.cosmology.model or None, optional
        Cosmological model instance.  If ``None`` a default model is
        constructed.
    thinness : int, optional
        Number of redshift points used to pre-compute the growth factor
        (default: ``1000``).
    """

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
        """Compute and store :math:`\\sigma_8` from the input power spectrum.

        Evaluates

        .. math::

            \\sigma_8^2 = \\frac{1}{2\\pi^2}
            \\int_0^\\infty k^2\\,P(k)\\,
            \\tilde{W}^2(8\\,h^{-1}\\,\\mathrm{Mpc}\\cdot k)\\,
            \\mathrm{d}k,

        where :math:`\\tilde{W}` is the Fourier transform of a top-hat
        window of radius :math:`8\\,h^{-1}\\,\\mathrm{Mpc}`.
        The result is stored as ``self.sigma8`` and a correction factor
        ``self.sigma8correction = cosmo.sigma8 / sigma8`` is set so that
        subsequent calls to :meth:`Pz` return correctly normalised spectra.

        Returns
        -------
        float
            The unnormalised :math:`\\sigma_8` computed from the raw input
            ``pk0``.
        """
        if self.sigma8 < 0.0 :
            integrand = self.kh**2 * self.pk0 * FT_tophat( 8.0*self.kh )**2
            self.sigma8 = 0.5 * trap_int(self.kh, integrand ) / numpy.pi**2
        self.sigma8correction = self.cosmo.param['sigma8'] / self.sigma8
        return self.sigma8
    
    def Pz ( self, kk, zz = 0.0, comoving = True ) :
        """Linear matter power spectrum evolved to redshift :math:`z`.

        .. math::

            P(k, z) = \\sigma_8^{\\mathrm{cosmo},2}\\,
            \\frac{D^2(z)}{D^2(0)}\\,
            \\frac{P_0(k)}{\\sigma_8^2},

        where :math:`P_0(k)` is the input spectrum at :math:`z=0` and
        :math:`D(z)` is the linear growth factor.

        Parameters
        ----------
        kk : scalar or array-like
            Wavenumbers in :math:`[h\\,\\mathrm{Mpc}^{-1}]` (comoving)
            or :math:`[\\mathrm{Mpc}^{-1}]` (physical).
        zz : scalar or array-like, optional
            Redshift(s) (default: ``0.0``).
        comoving : bool, optional
            If ``True`` (default) wavenumbers and the output are in comoving
            units; if ``False`` a factor of :math:`h` is applied.

        Returns
        -------
        ndarray
            Power spectrum values with shape ``(kk.size,)`` or
            ``(kk.size, zz.size)``.
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
        """Variance of the linear density field smoothed on scale :math:`R`.

        .. math::

            \\sigma^2(R, z) = \\frac{1}{2\\pi^2}
            \\int_0^\\infty k^2\\,P(k,z)\\,
            \\tilde{W}^2(kR)\\,\\mathrm{d}k,

        where :math:`\\tilde{W}(x) = 3(\\sin x - x\\cos x)/x^3` is the
        Fourier transform of a spherical top-hat of radius :math:`R`.

        Parameters
        ----------
        rr : scalar or array-like
            Smoothing radii in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
        zz : scalar or array-like
            Redshift(s).  The output is broadcast over ``(rr, zz)``.
        comoving : bool, optional
            Passed to :meth:`Pz` (default: ``True``).

        Returns
        -------
        ndarray
            Variance with shape ``(rr.size,)`` or ``(rr.size, zz.size)``.
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
        """Derivative of :math:`\\sigma^2(R,z)` with respect to :math:`R`.

        .. math::

            \\frac{\\mathrm{d}\\sigma^2}{\\mathrm{d}R}(R, z) =
            \\frac{1}{\\pi^2}
            \\int_0^\\infty k^3\\,P(k,z)\\,
            \\tilde{W}(kR)\\,\\tilde{W}'(kR)\\,\\mathrm{d}k,

        where :math:`\\tilde{W}'` denotes the derivative of the top-hat
        window with respect to its argument.
        Used by :meth:`sigma2M` and :meth:`dsigma2MdM` via the chain rule.

        Parameters
        ----------
        rr : scalar or array-like
            Smoothing radii in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
        zz : scalar or array-like
            Redshift(s).  The output is broadcast over ``(rr, zz)``.
        comoving : bool, optional
            Passed to :meth:`Pz` (default: ``True``).

        Returns
        -------
        ndarray
            Derivative with shape ``(rr.size,)`` or ``(rr.size, zz.size)``.
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
        """Variance of the density field on the mass scale :math:`M`.

        Converts a halo mass to the radius of an equivalent sphere of mean
        matter density,

        .. math::

            R(M) = \\left(\\frac{3M}{4\\pi\\rho_0}\\right)^{1/3},

        and returns :math:`\\sigma^2(R(M), z)` via :meth:`sigma2R`.

        Parameters
        ----------
        mm : scalar or array-like
            Halo masses in :math:`[M_\\odot\\,h^{-1}]`.
        zz : scalar or array-like
            Redshift(s).  The output is broadcast over ``(mm, zz)``.
        comoving : bool, optional
            If ``True`` (default) the mean density :math:`\\rho_0` is the
            comoving matter density; otherwise the physical value is used.

        Returns
        -------
        ndarray
            Variance with shape ``(mm.size,)`` or ``(mm.size, zz.size)``.
        """

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
        """Derivative of :math:`\\sigma^2(M,z)` with respect to :math:`M`.

        Applies the chain rule through the :math:`R(M)` relation:

        .. math::

            \\frac{\\mathrm{d}\\sigma^2}{\\mathrm{d}M} =
            \\frac{R}{3M}\\,\\frac{\\mathrm{d}\\sigma^2}{\\mathrm{d}R}.

        Used internally by halo mass function routines.

        Parameters
        ----------
        mm : scalar or array-like
            Halo masses in :math:`[M_\\odot\\,h^{-1}]`.
        zz : scalar or array-like
            Redshift(s).  The output is broadcast over ``(mm, zz)``.
        comoving : bool, optional
            If ``True`` (default) uses the comoving matter density.

        Returns
        -------
        ndarray
            Derivative with shape ``(mm.size,)`` or ``(mm.size, zz.size)``.
        """

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
        """Linear two-point correlation function.

        Obtained by Fourier-transforming the linear power spectrum
        (Eq. 6 of Ronconi et al. 2020):

        .. math::

            \\xi(r, z) = \\frac{1}{2\\pi^2}
            \\int_0^\\infty k^2\\,P(k,z)\\,
            \\frac{\\sin(kr)}{kr}\\,\\mathrm{d}k.

        The integral is evaluated via a Fast Spherical Transform.

        Parameters
        ----------
        rr : scalar or array-like
            Comoving separations in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
        zz : float, optional
            Redshift (default: ``0.0``).  Broadcasting over ``zz`` is not
            yet supported.
        comoving : bool, optional
            Passed to :meth:`Pz` (default: ``True``).

        Returns
        -------
        ndarray
            Correlation function :math:`\\xi(r,z)` evaluated at ``rr``.
        """
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rn, xin = fstl( self.kh, self.Pz( self.kh, zz, comoving ) )
        return lin_interp( rn, xin )( rr )                        

#############################################################################################
