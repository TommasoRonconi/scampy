import numpy
from .internal.cwrap import *

#Wrap class cosmology:
class cosmology () :
    """ Class to handle cosmology parameterisation and functions

    Parameters
    ----------
    kh0 : array-like object, same dimension of pk0
      binned scale parameter :math:`k/h` (logarithmically equispaced)
    pk0 : array-like object, same dimension of kh0
      input power spectrum computed in the positions of `kh0` at redshift z = 0
      ( the :math:`P(k)` will be authomatically re-normalized at the required value of :math:`\\sigma_8` )

    zmin : float
      minimum value of redshift for interpolation (note: must be :math:`> 0`)

    zmax : float
      maximum value of redshift for interpolation

    thinness : int
      number of bins of the interpolation in redshift-space
        
    Note
    ----
    If the `hh` parameter has been set to `1` (e.g. using function cosmology.set_parameters())
    then the output units of functions (if not conversely specified) are equivalent to their 
    comoving counterpart (e.g. :math:`[\\text{Mpc}]\\ \\rightarrow\\ [\\text{Mpc}/h]`
    
    Warning
    -------
    Internally functions use an interpolated method on a logarithmic bin in redshift, 
    except where conversely specified.
    (i.e. do not try to compute it with :math:`z = 0`, use a value :math:`\\ge 10^{-7}`)
    """

    params = { 'Om_M' : 0.3,
               'Om_b' : 0.045,
               'Om_L' : 0.7,
               'Om_n' : 0.00,
               'Om_r' : 0.00,
               'Om_K' : 0.00,
               'hh' : 0.7,
               'sigma8' : 0.8 }
    
    def __init__ ( self, kh0, pk0,
                   zmin = 1.e-7, zmax = 1.e+7,
                   thinness = 200 ) :
        """ Constructor of the class cosmology
        """
        if len( kh0 ) == len( pk0 ) :
            self._size_k = len( kh0 )
        else :
            raise ValueError( "list kh0 and list pk0 must have the same length" )
        
        self._c_par = c_cosmo_params_t( self.params[ 'Om_M' ],
                                        self.params[ 'Om_b' ],
                                        self.params[ 'Om_L' ],
                                        self.params[ 'Om_n' ],
                                        self.params[ 'Om_r' ],
                                        self.params[ 'Om_K' ],
                                        self.params[ 'hh' ],
                                        self.params[ 'sigma8' ] )
        
        self._kh0 = ( c_double * len( kh0 ) )( *[ _k for _k in kh0 ] )
        self._pk0 = ( c_double * len( pk0 ) )( *[ _p for _p in pk0 ] )
        self._zmin = zmin
        self._zmax = zmax
        self._thinness = thinness
        self.obj = lib_cosmo.create_cosmology( self._c_par,
                                               self._kh0, self._pk0,
                                               c_size_t( self._size_k ),
                                               c_double( self._zmin ), c_double( self._zmax ),
                                               c_size_t( self._thinness ) )


    def __del__ ( self ) :
        """ Destructor of the class cosmology
        It calls the corresponding destructor wrapped from C
        """

        # Python call to cosmology dtor:
        lib_cosmo.free_cosmology( self.obj )

    def set_parameters ( self,
                         Om_M = 0.3,
                         Om_b = 0.045,
                         Om_L = 0.7,
                         Om_n = 0.,
                         Om_r = 0.,
                         Om_K = 0.,
                         hh = 0.7,
                         sigma8 = 0.8 ) :
        """ Function to set cosmological parameters

        Parameters
        ----------
        Om_M : float
          :math:`\\Omega_M` matter density parameter
        Om_b : float
          :math:`\\Omega_b` baryonic matter density parameter
        Om_L : float
          :math:`\\Omega_\\Lambda` dark-energy density parameter
        Om_n : float
          :math:`\\Omega_\\nu` neutrinos density parameter
        Om_r : float
          :math:`\\Omega_r` radiation density parameter
        Om_K : float
          :math:`\\Omega_K` curvature density parameter
        hh : float
          :math:`h` Hubble parameter :math:`H/(100\; [km \cdot s^{-1} \cdot Mpc^{-1}])` at :math:`z = 0`
        sigma8 : float
          :math:`\\sigma_8` normalization at :math:`z = 0`

        Returns 
        -------
        None
        """
        self._c_par = c_cosmo_params_t( Om_M, Om_b, Om_L, Om_n, Om_r, Om_K, hh, sigma8 )
        self.params= { 'Om_M' : Om_M,
                       'Om_b' : Om_b,
                       'Om_L' : Om_L,
                       'Om_n' : Om_n,
                       'Om_r' : Om_r,
                       'Om_K' : Om_K,
                       'hh' : hh,
                       'sigma8' : sigma8 }
        self.obj = lib_cosmo.create_cosmology( self._c_par,
                                                     self._kh0, self._pk0,
                                                     c_size_t( self._size_k ),
                                                     c_double( self._zmin ), c_double( self._zmax ),
                                                     c_size_t( self._thinness ) )
        return

    def Hz ( self, zz ) :
        """ Computes the Hubble-parameter ad a given redshift.\n
        
        It computes
        
        .. math:: H(z) = H_0 \\cdot E(z)
        
        with
        
        .. math:: E^2(z) = ( 1 + z )^2 \\bigl[ \sum_i \Omega_i ( 1 + z )^{ 1 + 3 w_i } \\bigr]
        and with :math:`w_M = 0`, :math:`w_r = w_\\nu = 1/3`, :math:`w_K = -1/3`

        Parameters
        ----------
        zz : float
           redshift

        Returns
        -------
        float
           Value of the Hubble parameter at given redshift
        """

        return lib_cosmo.cosmo_Hz( c_double( zz ), self.obj )

    def dC ( self, zz ) :
        """ Comoving distance at given redshift.\n
        
        It computes the integral
        
        .. math:: d_C(z) = d_{H, 0} \\int_0^{z} \\dfrac{d z}{E(z)}
        
        where :math:`d_{H, 0} = c/H_0` is the Hubble distance at :math:`z = 0`

        Parameters
        ----------
        zz : float
           redshift

        Returns
        -------
        float
          the comoving distance at redshift `zz` in units of :math:`[\\text{Mpc}]`
        """

        return lib_cosmo.cosmo_dC( c_double( zz ), self.obj )

    def ddCdz ( self, zz ) :
        """ Derivative of the comoving distance at given redshift.\n
        
        It computes
        
        .. math:: d_C(z) = d_{H, 0} \\dfrac{d z}{E(z)}
        
        where :math:`d_{H, 0} = c/H_0` is the Hubble distance at :math:`z = 0`

        Parameters
        ----------
        zz : float
           redshift

        Returns
        -------
        float
          the comoving distance at redshift `zz` in units of :math:`[\\text{Mpc}]`
        """

        return lib_cosmo.cosmo_ddC( c_double( zz ), self.obj )

    def comoving_volume_unit ( self, zz ) :
        """ Comoving volume in unit of redshift and solid angle :math:`d\\Omega`.\n
        
        It computes
        
        .. math:: \\dfrac{dV(z)}{dz d\\Omega} = \\dfrac{d_C^2(z)}{d_{H, 0}} \\int_0^{z} \\dfrac{d z}{E(z)}
        
        where :math:`d_{H, 0} = c/H_0` and `d_C(z)` are the Hubble distance at :math:`z = 0` and
        the comoving distance at redshift :math:`z`, respectively.

        Parameters
        ----------
        zz : float
           redshift

        Returns
        -------
        float
          the comoving volume unit at redshift `zz` in units of 
          :math:`[ \\text{Mpc}^3 \\text{rad}^{-1} ]`
        """

        return lib_cosmo.cosmo_comoving_volume_unit( c_double( zz ), self.obj )

    def comoving_volume ( self, zz ) :
        """ Comoving volume at given redshift.\n
        
        It computes
        
        .. math:: V_C = \\dfrac{4 \\pi}{3} d_C^3(z)

        where :math:`d_C(z)` is the comoving distance computed as in cosmology.dC( z )
        
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          the all-sky, total comoving volume up to redshift `zz`
        """

        return lib_cosmo.cosmo_comoving_volume( c_double( zz ), self.obj )

    def cosmic_time ( self, zz ) :
        """ Cosmic time at given redshift.

        It gives the age of the Universe at a given redshift

        .. math:: t_C(z) = t_{H, 0} \\int_z^{\\infty} \\dfrac{dz'}{(1 + z') E(z')}

        where :math:`t_{H, 0}` is the Hubble time at :math:`z = 0`, which is computed as
        :math:`\\propto 1/ H_0`.


        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          Age of the Universe in :math:`[\\text{Gyr}]`
        """

        return lib_cosmo.cosmo_cosmic_time( c_double( zz ), self.obj )

    def rho_crit ( self, zz ) :
        """ Critical matter density at redshift :math:`z = 0`
        
        .. math:: \\rho_\\text{crit}(z) = \\dfrac{8 \\pi G}{3 H(z)}

        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          the critical overdensity :math:`\\rho_\\text{crit}` 
          in units of :math:`[M_\\odot \\text{Mpc}^{-3}]`
        """

        return lib_cosmo.cosmo_rho_crit( c_double( zz ), self.obj )

    def rho_crit_comoving ( self, zz ) :
        """ Comoving critical matter density at redshift :math:`z = 0`
        
        .. math:: \\rho_\\text{crit}(z) = \\dfrac{8 \\pi G}{3 H(z)}

        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          the comoving value of the critical matter overdensity :math:`\\rho_\\text{crit}` 
          in units of :math:`[M_\\odot \\text{Mpc}^{-3} h^2]`
          
        Note
        ----
        The returned value is independent to the chosen parameter `hh`
        """

        return lib_cosmo.cosmo_rho_crit_comoving( c_double( zz ), self.obj )

    def OmegaM ( self, zz ) :
        """ The dimensioneless matter density parameter at given redshift.

        It is computed as
        
        .. math:: \\Omega_M(z) = \\dfrac{\\Omega_{M, 0}}{E^2(a) (1+z) }

        where :math:`\\Omega_{M, 0}` is the matter density parameter at redshift :math:`z=0`
        and :math:`E^2(a)` is given by the sum over all the :math:`i` energy components:

        .. math:: E^2(a) = \\sum_i \\Omega_i \\biggl(\\dfrac{a}{a_0}\\biggr)^{w_i}

        with :math:`a/a_0` the scale factor at the target cosmic time normalized with respect
        to its value today 
        
        
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          :math:`\\Omega_M(z)`
        """

        return lib_cosmo.cosmo_OmegaM( c_double( zz ), self.obj )

    def Omegab ( self, zz ) :
        """ The dimensioneless baryonic matter density parameter at given redshift.

        It is computed as
        
        .. math:: \\Omega_b(z) = \\dfrac{\\Omega_{b, 0}}{E^2(a) (1+z) }

        where :math:`\\Omega_{b, 0}` is the baryon density parameter at redshift :math:`z=0`
        and :math:`E^2(a)` is given by the sum over all the :math:`i` energy components:

        .. math:: E^2(a) = \\sum_i \\Omega_i \\biggl(\\dfrac{a}{a_0}\\biggr)^{w_i}

        with :math:`a/a_0` the scale factor at the target cosmic time normalized with respect
        to its value today 
        
        
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          :math:`\\Omega_b(z)`
        """

        return lib_cosmo.cosmo_Omegab( c_double( zz ), self.obj )

    def deltac ( self, zz ) :
        """ Linear extrapolation of the critical density contrast at given redshift
        
        It is computed as

        .. math:: \\delta_c = 1.686 * ( 1 + 0.012299 * \\log( \\Omega_M(z) ) )
        
        where :math:`\\Omega_M(z)` is obtained with function OmegaM
        
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          The value of the spherical collapse density threshold, :math:`\\delta_c`, 
          by using the approximated equation (C.28) provided by Nakamura & Suto (1997)
        """

        return lib_cosmo.cosmo_deltac( c_double( zz ), self.obj )

    def Deltac_BN98 ( self, zz ) :
        """
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Deltac_BN98( c_double( zz ), self.obj )

    def Deltac_NS98 ( self, zz ) :
        """
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Deltac_NS98( c_double( zz ), self.obj )

    def DD ( self, zz ) :
        """ Computes the amplitude of the growing mode :math:`D(z)` as a function of redshift
        
        Readapted from Hamilton2001 ( arXiv:astro-ph/0006089v3, note to Eq. 1 ) computes

        .. math:: D(z) \equiv \dfrac{5 \Omega_{M,0}}{2} E( z ) \int_z^\infty \dfrac{(1 + z) dz}{E(z)}

        Parameters
        ----------
        zz : float
          redshift at which to compute :math:`D( z )`

        Returns
        -------
        float
          not normalized value of :math:`D(z)`
        """

        return lib_cosmo.cosmo_DD( c_double( zz ), self.obj )

    def gz ( self, zz ) :
        """ Growth factor at redshift `zz`
        
        See cosmology.DD for further information
        
        Parameters
        ----------
        zz : float
          redshift at which to compute :math:`g( z )`

        Returns
        -------
        float
          growth factor at redshift `zz`
        """

        return lib_cosmo.cosmo_gz( c_double( zz ), self.obj )

    def Pk ( self, kk, zz ) :
        """
        Parameters
        ----------
        kk : float
          scale 
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Pk( c_double( kk ), c_double( zz ), self.obj )

    def sigma2M ( self, mm, zz ) :
        """
        Parameters
        ----------
        mm : float
          mass
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_sigma2M( c_double( mm ), c_double( zz ), self.obj )

    def dndM ( self, mm, zz ) :
        """
        Parameters
        ----------
        mm : float
          mass
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dndM( c_double( mm ), c_double( zz ), self.obj )

    def hbias ( self, mm, zz ) :
        """
        Parameters
        ----------
        mm : float
          mass
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_hbias( c_double( mm ), c_double( zz ), self.obj )

    def density_profile_FS ( self, kk, mm, zz ) :
        """
        Parameters
        ----------
        kk : float
          scale
        mm : float
          mass
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_density_profile_FS( c_double( kk ), c_double( mm ), c_double( zz ), self.obj )

    def dphidL ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : float
          luminosity [ mag ]
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_B15 ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : float
          luminosity [ mag ]
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Bouwens15( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_B16 ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : float
          luminosity [ mag ]
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Bouwens16( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_L17_uv ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : float
          luminosity [ mag ]
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Lapi17_uv( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_L17_uvir ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : float
          luminosity [ mag ]
        zz : float
          redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Lapi17_uvir( c_double( ll ), c_double( zz ), self.obj )

    def logsfr_from_Muv ( self, Muv ) :
        """
        Parameters
        ----------
        Muv : float
          UV luminosity in [ mag ]

        Returns
        -------
        float
          :math:`\\log( S.F.R. )` in :math:`[ M_sol / yr ]`
        """

        return - 7.4 - 0.4 * Muv

    def ionizing_photons ( self, sfr, fesc, kion ) :
        """
        Parameters
        ----------
        sfr : float
          star formation rate (S.F.R.)
        fesc : float
          escape fraction ( float between 0 and 1 )
        kion : float
          average production of ionizing photons per unit S.F.R.

        Returns
        -------
        float
          number of ionizing photons produced by the source per unit time [ sec^-1 ]
        """

        return sfr * fesc * kion

    def hydrogen_density ( self, zz ) :
        """
        Parameters
        ----------
        zz : float
          redshift

        Returns
        -------
        float
          Hydrogen mean density at given redshift [ cm^-3 ]
        """      
        return 2.e-7 * self.Omegab( zz ) * 45.454545455


    def recombination_time ( self, zz, cHII = 3. ) :
        """
        Parameters
        ----------
        zz : float
          redshift
        cHII : float
          clumping factor of HII regions ( default = 3 )
        
        Returns
        -------
        float
          Recombination time in [ Gyr ]
        """

        red = 7. / ( 1. + zz )
        
        return 3.2 * red * red * red / cHII

    def stromgren_radius ( self, tt, Nion, nH, trec ) :
        """
        Parameters
        ----------
        tt : float
          time elapsed since the production of ionizing photons started ( in [ Gyr ] )
        Nion : float
          number of ionizing photons produced in unit time
        nH : float
          average density of neutral hydrogen around source
        trec : float
          recombination time in [ Gyr ]

        Returns
        -------
        float
          radius of the completelly ionized region around source (i.e. Stromgren sphere )
        """

        return ( 0.75 * ( numpy.pi * 1.e+7 * Nion ) * ( 1.e+9 * trec ) * ( 1. - numpy.exp( - tt / trec ) ) / ( nH * numpy.pi ) )**0.3333333333333333 * 3.2407789e-25
    

    
