import numpy
from .cwrap.cwrap import *

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
        
        self._c_par = c_cosmo_params_t( params[ 'Om_M' ],
                                        params[ 'Om_b' ],
                                        params[ 'Om_L' ],
                                        params[ 'Om_n' ],
                                        params[ 'Om_r' ],
                                        params[ 'Om_K' ],
                                        params[ 'hh' ],
                                        params[ 'sigma8' ] )
        
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
        Om_b : float
        Om_L : float
        Om_n : float
        Om_r : float
        Om_K : float
        hh : float
        sigma8 : float

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
        """ Computes the Hubble-parameter ad a given redshift.
        Namely it computes
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

        Warning
        -------
        Internally it uses an interpolated method on a logarithmic bin in redshift
        (i.e. do not try to compute it with :math:`z = 0`, use a value :math:`\\ge 10^{-7}`)
        """

        return lib_cosmo.cosmo_Hz( c_double( zz ), self.obj )

    def dC ( self, zz ) :
        """ Comoving distance at given redshift

        Parameters
        ----------
        zz : float
           redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dC( c_double( zz ), self.obj )

    def ddCdz ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_ddC( c_double( zz ), self.obj )

    def comoving_volume_unit ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_comoving_volume_unit( c_double( zz ), self.obj )

    def comoving_volume ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_comoving_volume( c_double( zz ), self.obj )

    def cosmic_time ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_cosmic_time( c_double( zz ), self.obj )

    def rho_crit ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_rho_crit( c_double( zz ), self.obj )

    def rho_crit_comoving ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_rho_crit_comoving( c_double( zz ), self.obj )

    def OmegaM ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_OmegaM( c_double( zz ), self.obj )

    def Omegab ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Omegab( c_double( zz ), self.obj )

    def deltac ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_deltac( c_double( zz ), self.obj )

    def Deltac_BN98 ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Deltac_BN98( c_double( zz ), self.obj )

    def Deltac_NS98 ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Deltac_NS98( c_double( zz ), self.obj )

    def DD ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_DD( c_double( zz ), self.obj )

    def gz ( self, zz ) :
        """
        Parameters
        ----------
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_gz( c_double( zz ), self.obj )

    def Pk ( self, kk, zz ) :
        """
        Parameters
        ----------
        kk : scale 
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_Pk( c_double( kk ), c_double( zz ), self.obj )

    def sigma2M ( self, mm, zz ) :
        """
        Parameters
        ----------
        mm : mass
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_sigma2M( c_double( mm ), c_double( zz ), self.obj )

    def dndM ( self, mm, zz ) :
        """
        Parameters
        ----------
        mm : mass
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dndM( c_double( mm ), c_double( zz ), self.obj )

    def hbias ( self, mm, zz ) :
        """
        Parameters
        ----------
        mm : mass
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_hbias( c_double( mm ), c_double( zz ), self.obj )

    def density_profile_FS ( self, kk, mm, zz ) :
        """
        Parameters
        ----------
        kk : scale
        mm : mass
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_density_profile_FS( c_double( kk ), c_double( mm ), c_double( zz ), self.obj )

    def dphidL ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : luminosity [ mag ]
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_B15 ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : luminosity [ mag ]
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Bouwens15( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_B16 ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : luminosity [ mag ]
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Bouwens16( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_L17_uv ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : luminosity [ mag ]
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Lapi17_uv( c_double( ll ), c_double( zz ), self.obj )

    def dphidL_L17_uvir ( self, ll, zz ) :
        """
        Parameters
        ----------
        ll : luminosity [ mag ]
        zz : redshift

        Returns
        -------
        float
        """

        return lib_cosmo.cosmo_dphidL_Lapi17_uvir( c_double( ll ), c_double( zz ), self.obj )

    def logsfr_from_Muv ( self, Muv ) :
        """
        Parameters
        ----------
        Muv : UV luminosity in [ mag ]

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
        sfr : star formation rate (S.F.R.)
        fesc : escape fraction ( float between 0 and 1 )
        kion : average production of ionizing photons per unit S.F.R.

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
        zz : redshift

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
        zz : redshift
        cHII : clumping factor of HII regions ( default = 3 )
        
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
        tt : time elapsed since the production of ionizing photons started ( in [ Gyr ] )
        Nion : number of ionizing photons produced in unit time
        nH : average density of neutral hydrogen around source
        trec : recombination time in [ Gyr ]

        Returns
        -------
        float
          radius of the completelly ionized region around source (i.e. Stromgren sphere )
        """

        return ( 0.75 * ( np.pi * 1.e+7 * Nion ) * ( 1.e+9 * trec ) * ( 1. - np.exp( - tt / trec ) ) / ( nH * np.pi ) )**0.3333333333333333 * 3.2407789e-25
    

        

if __name__ == '__main__' :

    input_dir = "../tests/integration_tests/input/"

    kh0, pk0 = numpy.genfromtxt( input_dir + "not-norm_pk_lcdm_camb.dat",
                                 unpack = True )
    cosmo = cosmology( kh0, pk0 )

    print( "H( z = 2 )\t=\t{:e}".format( cosmo.Hz( 2. ) ) )
    print( "tt( z = 3 )\t=\t{:e}".format( cosmo.cosmic_time( 3. ) ) )
    print( "dV( z = 2 )\t=\t{:e}".format( cosmo.comoving_volume_unit( 2. ) ) )
    print( "V( z = 2 )\t=\t{:e}".format( cosmo.comoving_volume( 2. ) ) )
    
