"""Flat or curved FLRW cosmological model.

Provides the :class:`model` class, which pre-tabulates cosmographic
quantities on a log-spaced redshift grid at construction time for fast
evaluation.

Available methods of :class:`model`
-------------------------------------

**Cosmographic distances**

.. list-table::
   :widths: 30 70

   * - :meth:`~model.Hz`
     - Hubble parameter H(z) [km/s/Mpc].
   * - :meth:`~model.dH`
     - Hubble distance c/H(z) [Mpc].
   * - :meth:`~model.dC`
     - Line-of-sight comoving distance [Mpc].
   * - :meth:`~model.ddC`
     - Derivative of comoving distance dd_C/dz [Mpc].
   * - :meth:`~model.dA`
     - Angular diameter distance [Mpc].

**Volume and time**

.. list-table::
   :widths: 30 70

   * - :meth:`~model.comoving_volume_unit`
     - Comoving volume element :math:`dV/(dz\\,d\\Omega)` :math:`[(\\mathrm{Mpc}/h)^3\\,\\mathrm{sr}^{-1}]`.
   * - :meth:`~model.comoving_volume`
     - Total comoving volume V(z) :math:`[(\\mathrm{Mpc}/h)^3]`.
   * - :meth:`~model.cosmic_time`
     - Age of the Universe at redshift z [Gyr].

**Density and matter content**

.. list-table::
   :widths: 30 70

   * - :meth:`~model.critical_density_comoving`
     - Critical density in comoving units :math:`[M_\\odot\\,\\mathrm{Mpc}^{-3}\\,h^2]`.
   * - :meth:`~model.critical_density`
     - Critical density in physical units :math:`[M_\\odot\\,\\mathrm{Mpc}^{-3}]`.
   * - :meth:`~model.OmegaM`
     - Total matter density parameter :math:`\\Omega_M(z)`.
   * - :meth:`~model.Omegab`
     - Baryonic matter density parameter :math:`\\Omega_b(z)`.

**Structure growth**

.. list-table::
   :widths: 30 70

   * - :meth:`~model.deltac`
     - Linear critical overdensity :math:`\\delta_c(z)`.
   * - :meth:`~model.D`
     - Linear growth factor D(z).
   * - :meth:`~model.gz`
     - Unnormalised growth factor (1+z) D(z).
   * - :meth:`~model.Delta_c`
     - Virial overdensity :math:`\\Delta_c(z)` (Nakamura & Suto 1998).
"""

import numpy
from scipy.integrate import cumulative_trapezoid

from scampy.utilities.interpolation import lin_interp

# Physical constants in SI units
_cc  = 2.99792458e8          # speed of light [m/s]
_pc  = 3.08567758130573e16   # parsec [m]
_yr  = 3600 * 24 * 365.256363004  # year [s]

# Critical density coefficient: rho_crit,0 for H0 = 100 km/s/Mpc
# = 3 H0^2 / (8 pi G)  in  [M_sun Mpc^-3 h^2]
# Verified analytically from the C++ constants (utilities.h):
#   3.75e18 * (1/pi) * (1/G_Newton) * (1/M_sun) * pc = 2.775e11
_RHO_CRIT0 = 2.7753586541352e11 #2.775e11


class model:
    """Flat or curved FLRW cosmological model.

    Pre-tabulates E(z), 1/E(z), the comoving-distance integrand, the
    linear growth factor, and the cosmic-time integrand on a log-spaced
    redshift grid at construction time.  All public methods are then fast
    interpolations or closed-form expressions.

    Parameters
    ----------
    Om_M : float, optional
        Total matter density parameter :math:`\\Omega_M` (default: 0.3).
    Om_b : float, optional
        Baryonic matter density parameter :math:`\\Omega_b` (default: 0.045).
    Om_L : float, optional
        Dark-energy density parameter :math:`\\Omega_\\Lambda` (default: 0.7).
    Om_n : float, optional
        Massive-neutrino density parameter (default: 0.0).
    Om_r : float, optional
        Radiation density parameter (default: 0.0).
    Om_K : float, optional
        Curvature density parameter (default: 0.0).
    hh : float, optional
        Dimensionless Hubble constant :math:`h = H_0/100` (default: 0.7).
    sigma8 : float, optional
        RMS matter-fluctuation amplitude at :math:`8\\,h^{-1}\\,\\mathrm{Mpc}`
        (default: 0.8).
    w_0 : float, optional
        Dark-energy equation-of-state constant term (default: -1.0).
    w_a : float, optional
        Dark-energy equation-of-state linear term
        :math:`w(a) = w_0 + w_a(1-a)` (default: 0.0).
    zmin : float, optional
        Minimum redshift of the internal tabulation grid (default: 0.0).
    zmax : float, optional
        Maximum redshift of the internal tabulation grid (default: 100.0).
    thin : int, optional
        Number of points in the internal tabulation grid (default: 1000).
    """

    def __init__(self,
                 Om_M=0.3, Om_b=0.045, Om_L=0.7,
                 Om_n=0.0, Om_r=0.0,  Om_K=0.0,
                 hh=0.7, sigma8=0.8,
                 w_0=-1.0, w_a=0.0,
                 zmin=0.0, zmax=100.0, thin=1000):

        self.param = {
            'Om_M': float(Om_M), 'Om_b': float(Om_b), 'Om_L': float(Om_L),
            'Om_n': float(Om_n), 'Om_r': float(Om_r), 'Om_K': float(Om_K),
            'hh':   float(hh),   'sigma8': float(sigma8),
            'w_0':  float(w_0),  'w_a':  float(w_a),
            'zmin': float(zmin), 'zmax': float(zmax), 'thin': float(thin),
        }

        self.zmin  = float(zmin)
        self.zmax  = float(zmax)

        self.H0   = 100.0 * hh                    # km/s/Mpc
        self.t_H0 = 1e3 * _pc / (_yr * self.H0)   # Hubble time [yr]
        self.d_H0 = 1e-3 * _cc / self.H0          # Hubble distance [Mpc]

        self._build_tables(int(thin))

    # ------------------------------------------------------------------
    # Internal setup

    def _build_tables(self, thin):
        # Log-spaced redshift grid (matches C++ logic exactly)
        if self.zmin > 0.0:
            zz = numpy.geomspace(self.zmin, self.zmax, thin)
        else:
            zz = numpy.geomspace(1e-3, self.zmax, thin)
            zz[0] = 0.0

        self._zz   = zz
        Ez2_tab    = self._Ez2(zz)
        Ez_tab     = numpy.sqrt(Ez2_tab)
        zE_tab     = 1.0 / Ez_tab

        self._Ez_f = lin_interp(zz, Ez_tab)   # E(z)
        self._zE_f = lin_interp(zz, zE_tab)   # 1/E(z)

        # Comoving-distance integral: int_0^z dz'/E(z')
        dC_int = cumulative_trapezoid(zE_tab, zz, initial=0.0)
        self._dC_int = lin_interp(zz, dC_int)

        # Growth-factor integral from right: int_z^zmax (1+z')/E(z')^3 dz'
        D_intgd   = (1.0 + zz) * zE_tab ** 3
        D_int     = (-1.0) * cumulative_trapezoid(D_intgd[::-1], zz[::-1], initial=0.0)[::-1]
        D_int[-1] = 0.0
        D_tab     = 2.5 * self.param['Om_M'] * Ez_tab * D_int
        self._D_f = lin_interp(zz, D_tab)

        # Cosmic-time integral from right: int_z^zmax dz'/[(1+z')*E(z')]
        ct_intgd  = zE_tab / (1.0 + zz)
        ct_int    = (-1.0) * cumulative_trapezoid(ct_intgd[::-1], zz[::-1], initial=0.0)[::-1]
        ct_int[-1] = 0.0
        self._ct_f = lin_interp(zz, ct_int)

    def _Ez2(self, zz):
        """E^2(z) = [H(z)/H0]^2 as a scalar or array."""
        p   = self.param
        red = 1.0 + numpy.asarray(zz, dtype=float)
        return ( ( p['Om_r'] + p['Om_n'] ) * red ** 4
                 + p['Om_M'] * red ** 3
                 + p['Om_K'] * red ** 2
                 + p['Om_L'] )

    def _Ea2(self, zz):
        """E^2(z) expressed in terms of the scale factor a = 1/(1+z)."""
        p  = self.param
        aa = 1.0 / (1.0 + numpy.asarray(zz, dtype=float))
        return ( p['Om_L'] * aa ** 4
                 + p['Om_K'] * aa ** 2
                 + p['Om_M'] * aa
                 + p['Om_r'] + p['Om_n'] )

    # ------------------------------------------------------------------
    # Cosmographic functions

    def Hz(self, zz):
        """Hubble parameter H(z) in km/s/Mpc.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            H(z) [km/s/Mpc].
        """
        return self.H0 * self._Ez_f(zz)

    def dH(self, zz):
        """Hubble distance c/H(z) in Mpc.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            c/H(z) [Mpc].
        """
        return 1e-3 * _cc / self.Hz(zz)

    def dC(self, zz):
        """Line-of-sight comoving distance in Mpc.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            d_C(z) [Mpc].
        """
        return self.d_H0 * self._dC_int(zz)

    def ddC(self, zz):
        """Derivative of the comoving distance with respect to redshift, dd_C/dz.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            dd_C/dz [Mpc].
        """
        return self.d_H0 * self._zE_f(zz)

    def dA(self, zz):
        """Angular diameter distance in Mpc.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            d_A(z) [Mpc].
        """
        return self.dC(zz) / (1.0 + numpy.asarray(zz, dtype=float))

    def comoving_volume_unit(self, zz):
        """Comoving volume element per unit redshift per unit steradian, dV/(dz dOmega).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            dV/(dz dOmega) [(Mpc/h)^3 / sr].
        """
        dC_int = self._dC_int(zz)
        return self.d_H0 ** 3 * self._zE_f(zz) * dC_int ** 2

    def comoving_volume(self, zz):
        """Total comoving volume enclosed within redshift z in (Mpc/h)^3.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            V_C(z) [(Mpc/h)^3].
        """
        dc = self.dC(zz)
        return (4.0 / 3.0) * numpy.pi * dc ** 3

    def cosmic_time(self, zz):
        """Age of the Universe at redshift z (lookback time from z to infinity) in Gyr.

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            t(z) [Gyr].
        """
        return 1e-9 * self.t_H0 * self._ct_f(zz)

    # ------------------------------------------------------------------
    # Density and structure-growth functions

    def critical_density_comoving(self, zz):
        """Critical density in comoving units [Msol Mpc^-3 h^2].

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            rho_crit(z) [Msol Mpc^-3 h^2].
        """
        return _RHO_CRIT0 * self._Ez2(zz)

    def critical_density(self, zz):
        """Critical density in physical units [Msol Mpc^-3].

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            rho_crit(z) [Msol Mpc^-3].
        """
        return self.critical_density_comoving(zz) * self.param['hh'] ** 2

    def OmegaM(self, zz):
        """Total matter density parameter Omega_M(z) = rho_M(z)/rho_crit(z).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            Omega_M(z).
        """
        return self.param['Om_M'] / (self._Ea2(zz) * (1.0 + numpy.asarray(zz, dtype=float)))

    def Omegab(self, zz):
        """Baryonic matter density parameter Omega_b(z).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            Omega_b(z).
        """
        return self.param['Om_b'] / (self._Ea2(zz) * (1.0 + numpy.asarray(zz, dtype=float)))

    def deltac(self, zz):
        """Linear critical overdensity for spherical collapse delta_c(z).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            delta_c(z).
        """
        return 1.686 * (1.0 + 0.012299 * numpy.log10(self.OmegaM(zz)))

    def D(self, zz):
        """Linear growth factor D(z), unnormalised (divide by D(0) for D(0)=1 convention).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            D(z).
        """
        return self._D_f(zz)

    def gz(self, zz):
        """Unnormalised growth factor g(z) = (1+z)*D(z).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            (1+z)*D(z).
        """
        return (1.0 + numpy.asarray(zz, dtype=float)) * self.D(zz)

    def Delta_c(self, zz):
        """Virial overdensity Delta_c(z) from Nakamura & Suto (1998).

        Parameters
        ----------
        zz : float or array-like
            Redshift.

        Returns
        -------
        float or ndarray
            Delta_c(z).
        """
        xx = (1.0 - self.param['Om_M']) ** (1.0 / 3.0) / (1.0 + numpy.asarray(zz, dtype=float))
        return 18.0 * numpy.pi * (1.0 + 0.4093 * xx ** 2.7152) * self.OmegaM(zz)
