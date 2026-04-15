"""Halo model for the two-point statistics of a galaxy population.

Implements the halo-model framework (Cooray & Sheth 2002) as described in
Sec. 2.1 of Ronconi et al. (2020).  Given a power spectrum and an HOD, the
:class:`halo_model` class provides the galaxy number density, effective bias,
and the 3D, projected, and angular two-point correlation functions decomposed
into their 1-halo and 2-halo contributions.
"""

#############################################################################################

# External imports
import numpy
from scipy.special import erf

# Internal imports
from scampy.power_spectrum import power_spectrum
from scampy.halo import mass_function
from scampy.halo import bias
from scampy.halo.density_profile import density_profile_FT
from scampy.hod import HOD
from scampy.utilities.functions import trap_int, linear_interpolation
from scampy.utilities.fft import fstl, fpstl
from scampy.utilities.interpolation import lin_interp

#############################################################################################

class halo_model () :
    """Halo model for galaxy two-point statistics.

    Encapsulates the halo-model integrals over the halo mass function,
    bias, and NFW density profile (see Ronconi et al. 2020, Sec. 2.1).
    All power-spectrum and correlation-function methods accept an HOD
    instance that defines the galaxy occupation probabilities
    :math:`\\langle N_\\mathrm{cen}\\rangle(M_h)` and
    :math:`\\langle N_\\mathrm{sat}\\rangle(M_h)`.

    Parameters
    ----------
    pk : scampy.power_spectrum.power_spectrum
        Power-spectrum object carrying the cosmological model and growth
        factor.
    comoving : bool, optional
        Whether to work in comoving coordinates (default: ``True``).
    Mlim : tuple of float, optional
        Log10 of the minimum and maximum halo mass for the integration
        grid, in :math:`M_\\odot\\,h^{-1}` (default: ``(5, 17)``).
    klim : tuple of float, optional
        Log10 of the minimum and maximum wavenumber for the internal
        :math:`k`-grid in :math:`h\\,\\mathrm{Mpc}^{-1}`
        (default: ``(-3, +3)``).
    thin : int, optional
        Number of points in both the mass and wavenumber grids
        (default: ``256``).
    hmf : str, optional
        Halo mass function label (placeholder, currently always
        :func:`~scampy.halo.mass_function.Tinker08`).
    hbias : str, optional
        Halo bias label (placeholder, currently always
        :func:`~scampy.halo.bias.Tinker10`).
    density : str, optional
        Density profile label (placeholder, currently always NFW via
        :func:`~scampy.halo.density_profile.density_profile_FT`).
    """

    def __init__ ( self, pk, comoving = True, Mlim = ( 5, 17 ),
                   klim = ( -3, +3 ), thin = 256,
                   hmf = 'Tinker08', hbias = 'Tinker10', 
                   # these strings are just placeholders for now
                   density = 'NFW' 
                 ) :
        if not isinstance( pk, power_spectrum ) :
            raise ValueError( 
                'argument pk should be an instance of type scampy.power_spectrum' 
            )
        self.pk = pk
        self.hmf = getattr( mass_function, hmf )
        self.hbs = getattr( bias, hbias )
        self.hdp = density_profile_FT
        self.Mh_min, self.Mh_max = Mlim
        self.k_min, self.k_max = klim
        self.Mh_grid = numpy.logspace( *Mlim, thin )
        self.kh_grid = numpy.logspace( *klim, thin )
        self._lkhdamp = numpy.log10( self.pk.kh[ self.pk.pk0.argmax() ] )
        self._sigmalk = 0.25
        
    def ng ( self, hod, zz = 0.0 ) :
        """Mean galaxy number density.

        Integrates the HOD-weighted halo mass function (Eq. 5 of
        Ronconi et al. 2020):

        .. math::

            n_g(z) = \\int_{M_\\mathrm{min}}^{M_\\mathrm{max}}
            \\langle N_g \\rangle(M_h)\\,n(M_h, z)\\,\\mathrm{d}M_h,

        where :math:`\\langle N_g\\rangle = \\langle N_\\mathrm{cen}\\rangle
        + \\langle N_\\mathrm{sat}\\rangle`.

        Parameters
        ----------
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).

        Returns
        -------
        float or ndarray
            Galaxy number density in
            :math:`[h^3\\,\\mathrm{Mpc}^{-3}]`.
        """
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
        """Effective large-scale halo bias of the galaxy population.

        .. math::

            b_g(z) = \\frac{1}{n_g(z)}
            \\int \\langle N_g\\rangle(M_h)\\,
            n(M_h, z)\\,b(M_h, z)\\,\\mathrm{d}M_h.

        Parameters
        ----------
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).

        Returns
        -------
        float or ndarray
            Dimensionless effective bias :math:`b_g(z)`.
        """
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
        """Mean host halo mass of the galaxy population.

        .. math::

            \\langle M_h \\rangle(z) = \\frac{1}{n_g(z)}
            \\int \\langle N_g\\rangle(M_h)\\,
            n(M_h, z)\\,M_h\\,\\mathrm{d}M_h.

        Parameters
        ----------
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).

        Returns
        -------
        float or ndarray
            Mean host halo mass in :math:`[M_\\odot\\,h^{-1}]`.
        """
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
        """Differential galaxy number density per unit halo mass.

        Returns the integrand of :meth:`ng` evaluated at specific mass
        values:

        .. math::

            \\frac{\\mathrm{d}n_g}{\\mathrm{d}M}(M_h, z) =
            \\langle N_g\\rangle(M_h)\\,n(M_h, z).

        Parameters
        ----------
        mm : scalar or array-like
            Halo masses in :math:`[M_\\odot\\,h^{-1}]`.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).

        Returns
        -------
        ndarray
            Differential density with shape ``(mm.size,)`` or
            ``(mm.size, zz.size)`` in
            :math:`[M_\\odot^{-1}\\,h^4\\,\\mathrm{Mpc}^{-3}]`.
        """
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
        """1-halo term of the galaxy power spectrum.

        Pairs within the same halo (central–satellite and
        satellite–satellite), Eq. 10 of Ronconi et al. (2020):

        .. math::

            P_\\mathrm{1h}(k,z) = \\frac{1}{n_g^2(z)}
            \\int \\left(2\\langle N_\\mathrm{cen}\\rangle +
            \\langle N_\\mathrm{sat}\\rangle\\right)
            \\langle N_\\mathrm{sat}\\rangle\\,
            n(M_h,z)\\,|\\tilde{u}(k,z|M_h)|^2\\,\\mathrm{d}M_h.

        Parameters
        ----------
        kk : scalar or array-like
            Wavenumbers in :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).
        fact : float or array-like, optional
            Pre-computed normalisation factor :math:`1/n_g^2`.  If
            negative (default) it is computed internally.

        Returns
        -------
        ndarray
            1-halo power spectrum with shape ``(kk.size,)`` or
            ``(kk.size, zz.size)`` in
            :math:`[h^{-3}\\,\\mathrm{Mpc}^3]`.
        """
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
        """2-halo term of the galaxy power spectrum.

        Pairs from different haloes, Eq. 11 of Ronconi et al. (2020):

        .. math::

            P_\\mathrm{2h}(k,z) = \\frac{P_\\mathrm{lin}(k,z)}{n_g^2(z)}
            \\left[\\int \\langle N_g\\rangle(M_h)\\,
            n(M_h,z)\\,b(M_h,z)\\,
            \\tilde{u}(k,z|M_h)\\,\\mathrm{d}M_h\\right]^2.

        Parameters
        ----------
        kk : scalar or array-like
            Wavenumbers in :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).
        fact : float or array-like, optional
            Pre-computed normalisation factor :math:`1/n_g^2`.  If
            negative (default) it is computed internally.

        Returns
        -------
        ndarray
            2-halo power spectrum with shape ``(kk.size,)`` or
            ``(kk.size, zz.size)`` in
            :math:`[h^{-3}\\,\\mathrm{Mpc}^3]`.
        """
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
        """Total halo-model galaxy power spectrum.

        Sum of the 1-halo and 2-halo terms (Eq. 7 of Ronconi et al. 2020):

        .. math::

            P(k,z) = P_\\mathrm{1h}(k,z) + P_\\mathrm{2h}(k,z).

        Parameters
        ----------
        kk : scalar or array-like
            Wavenumbers in :math:`[h\\,\\mathrm{Mpc}^{-1}]`.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float or array-like, optional
            Redshift(s) (default: ``0.0``).
        fact : float or array-like, optional
            Pre-computed normalisation factor :math:`1/n_g^2`.  If
            negative (default) it is computed internally.

        Returns
        -------
        ndarray
            Total power spectrum with shape ``(kk.size,)`` or
            ``(kk.size, zz.size)`` in
            :math:`[h^{-3}\\,\\mathrm{Mpc}^3]`.
        """
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
        """1-halo term of the 3D two-point correlation function.

        Obtained by Fourier-transforming :meth:`Pk_1halo` via a Fast
        Spherical Transform (Eq. 6 of Ronconi et al. 2020).

        Parameters
        ----------
        rr : scalar or array-like
            Comoving separations in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float, optional
            Redshift (default: ``0.0``).  Broadcasting over ``zz`` is
            not yet supported.
        fact : float, optional
            Pre-computed normalisation factor.  If negative (default)
            it is computed internally.

        Returns
        -------
        ndarray
            1-halo correlation function :math:`\\xi_\\mathrm{1h}(r)`
            evaluated at ``rr``.
        """
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
        """2-halo term of the 3D two-point correlation function.

        Fourier transform of :meth:`Pk_2halo`.  Same signature as
        :meth:`Xi_1halo`.
        """
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
        """Total 3D two-point correlation function.

        Sum of :meth:`Xi_1halo` and :meth:`Xi_2halo`.  Same signature as
        :meth:`Xi_1halo`.
        """
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
        """1-halo term of the projected two-point correlation function.

        Obtained by Abel-transforming :meth:`Pk_1halo` via a Fast
        Projected Spherical Transform (Eq. 15 of Ronconi et al. 2020,
        zeroth-order Hankel transform).

        Parameters
        ----------
        rp : scalar or array-like
            Projected comoving separations in :math:`[h^{-1}\\,\\mathrm{Mpc}]`.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float, optional
            Redshift (default: ``0.0``).  Broadcasting over ``zz`` is
            not yet supported.
        fact : float, optional
            Pre-computed normalisation factor.  If negative (default)
            it is computed internally.

        Returns
        -------
        ndarray
            1-halo projected correlation function
            :math:`w_\\mathrm{1h}(r_p)` evaluated at ``rp``.
        """
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
        """2-halo term of the projected two-point correlation function.

        Abel transform of :meth:`Pk_2halo`.  Same signature as
        :meth:`Wr_1halo`.
        """
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
        """Total projected two-point correlation function.

        Sum of :meth:`Wr_1halo` and :meth:`Wr_2halo`.  Same signature as
        :meth:`Wr_1halo`.
        """
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
        """1-halo term of the angular two-point correlation function.

        Converts the angular separation :math:`\\theta` to a projected
        comoving separation :math:`r_p = \\theta\\,d_C(z)` and returns
        :meth:`Wr_1halo` (flat-sky / Limber approximation,
        Eqs. 15–16 of Ronconi et al. 2020).

        Parameters
        ----------
        th : scalar or array-like
            Angular separations in radians.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zz : float, optional
            Redshift (default: ``0.0``).  Broadcasting over ``zz`` is
            not yet supported.
        fact : float, optional
            Pre-computed normalisation factor.  If negative (default)
            it is computed internally.

        Returns
        -------
        ndarray
            1-halo angular correlation function
            :math:`\\omega_\\mathrm{1h}(\\theta)` evaluated at ``th``.
        """
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rp = th * self.pk.cosmo.dC( zz )
        return self.Wr_1halo( rp, hod, zz, fact )

    def Wt_2halo ( self, th, hod, zz = 0.0, fact = -1.0 ) :
        """2-halo term of the angular two-point correlation function.

        Converts :math:`\\theta` to :math:`r_p` and returns
        :meth:`Wr_2halo`.  Same signature as :meth:`Wt_1halo`.
        """
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rp = th * self.pk.cosmo.dC( zz )
        return self.Wr_2halo( rp, hod, zz, fact )

    def Wt ( self, th, hod, zz = 0.0, fact = -1.0 ) :
        """Total angular two-point correlation function.

        Sum of :meth:`Wt_1halo` and :meth:`Wt_2halo`.  Same signature as
        :meth:`Wt_1halo`.
        """
        if hasattr( zz, '__len__' ) :
            raise RuntimeError(
                'function not supporting broadcasting to 2nd dimension yet' )
        rp = th * self.pk.cosmo.dC( zz )
        return self.Wr( rp, hod, zz, fact )

    def Nz ( self, hod, zbins, solid_angle = 4 * numpy.pi ) :
        """Normalised redshift distribution of the galaxy population.

        Computes the fraction of sources per unit redshift interval:

        .. math::

            N(z) = \\frac{n_g(z)\\,\\frac{\\mathrm{d}V_C}{\\mathrm{d}z}\\,
            \\Omega\\,\\Delta z}{\\Delta z\\,\\sum_i n_g(z_i)\\,
            \\frac{\\mathrm{d}V_C}{\\mathrm{d}z_i}\\,\\Omega\\,\\Delta z_i},

        where :math:`\\mathrm{d}V_C/\\mathrm{d}z` is the comoving volume
        element and :math:`\\Omega` is the solid angle.

        Parameters
        ----------
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zbins : array-like
            Bin edges in redshift.  The bin centres are taken as the
            midpoints.
        solid_angle : float, optional
            Solid angle of the survey in steradians
            (default: :math:`4\\pi`, full sky).

        Returns
        -------
        ndarray
            Normalised :math:`N(z)` with ``zbins.size - 1`` elements,
            such that :math:`\\sum_i N(z_i)\\,\\Delta z_i = 1`.
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
        """Angular two-point correlation function averaged over a redshift distribution.

        Implements Eq. 17 of Ronconi et al. (2020) in the Limber approximation:

        .. math::

            \\omega(\\theta) = \\int \\frac{\\mathrm{d}V}{\\mathrm{d}z}\\,
            N^2(z)\\,
            \\omega\\!\\left[r_p(\\theta, z),\\, z\\right]\\mathrm{d}z,

        where :math:`N(z)` is the normalised source redshift
        distribution, :math:`r_p(\\theta, z) = \\theta\\,d_C(z)` is the
        projected comoving separation at the comoving distance :math:`d_C(z)`,
        and :math:`\\mathrm{d}V/\\mathrm{d}z` is the comoving volume element.
        The integral is evaluated by computing :math:`\\omega(r_p, z)` via a
        Fast Projected Spherical Transform at each redshift bin centre and
        summing with trapezoidal weights.

        See Also
        --------
        Nz : Computes :math:`N(z)` from the HOD and cosmology; called
             internally when ``Nzdist`` is not provided.

        Parameters
        ----------
        th : array-like
            Angular separations in radians.
        hod : scampy.hod.HOD or compatible
            Object providing ``Pcen(M)`` and ``Psat(M)`` methods.
        zbins : array-like
            Redshift bin edges.  Bin centres are taken as midpoints.
        Nzdist : array-like or None, optional
            Normalised redshift distribution :math:`N(z_i)` as returned
            by :math:`Nz`,
            with ``Nzdist.size == zbins.size - 1``.  If ``None``
            (default) it is computed internally via :math:`Nz`.
        fact : float or array-like, optional
            Pre-computed normalisation factor :math:`1/n_g^2`.  If
            negative (default) it is computed internally.
        valid_rp : tuple of float, optional
            Range of projected separations :math:`(r_{p,\\min}, r_{p,\\max})`
            in :math:`[h^{-1}\\,\\mathrm{Mpc}]` over which the interpolation
            of the fixed-:math:`z` angular correlation is trusted
            (default: ``(1e-2, 8e1)``).
        solid_angle : float, optional
            Survey solid angle in steradians, used when ``Nzdist`` is
            computed internally (default: :math:`4\\pi`, full sky).
        component : str, optional
            Which halo-model term to include: ``'total'`` (default),
            ``'1halo'``, or ``'2halo'``.

        Returns
        -------
        ndarray
            Angular correlation function :math:`\\omega(\\theta)` with the
            same shape as ``th``.
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
