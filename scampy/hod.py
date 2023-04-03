import numpy
from scipy.special import erf

from scampy.catalogue import catalogue

def get_hosts ( cat, Pcen, Psat, 
                kw_pcen = {}, 
                kw_psat = {}, 
                rng = None, kw_rng = {} ) :
    """
    """

    # Check input validity
    if not isinstance( cat, catalogue ) :
        raise ValueError( 'Argument cat should be an instance of type scampy.catalogue' )
    
    # Build random number generator
    if rng is None :
        rng = np.random.default_rng(**kw_rng)
        
    # Central probability
    pcg = np.ones(self.haloes.size)
    if hasattr(Pcen, '__call__') :
        pcg = Pcen( self.haloes.Mhalo, **kw_pcen )
    elif hasattr( Pcen, '__len__' ) :
        if len(Pcen) != self.haloes.size :
            raise ValueError( 'len(Pcen) != haloes.size' )
        pcg = Pcen
    else :
        raise ValueError( 'wrong value passed to argument Pcen' )
    
    # Satellites probability
    psg = np.ones(self.haloes.size)
    if hasattr(Psat, '__call__') :
        Nsg = Psat( self.haloes.Mhalo, **kw_psat, **kw_pcen )
        Nsh = self.Nsat()
        wws = ( 0.0 < Nsh ) & ( Nsh > Nsg )
        psg[wws] = Nsg[wws] / Nsh[wws]
    elif hasattr( Psat, '__len__' ) :
        if len(Psat) != self.haloes.size :
            raise ValueError( 'len(Psat) != haloes.size' )
        psg = Psat
    else :
        raise ValueError( 'wrong value passed to argument Psat' )
        
    # Find centrals mask
    cen_idx = self.haloes.centrals()
    cen_gxy = np.zeros(self.subhaloes.size, dtype=bool)
    cen_gxy[cen_idx] = rng.binomial(1, pcg[self.subhaloes.Parent[cen_idx]])
        
    # Find satellites mask
    sat_idx = self.haloes.satellites()
    sat_gxy = np.zeros(self.subhaloes.size, dtype=bool)
    sat_gxy[sat_idx] = rng.binomial(1, psg[self.subhaloes.Parent[sat_idx]])
    
    return cen_gxy, sat_gxy


class HOD () :
    """Halo Occupation Distribution class, given a set of parameters
    computes the average number of central and satellite galaxies hosted by 
    a halo of a given mass.

    .. math::
      
      <N_\\text{cen}(M_h)> = 
      <N_\\text{sat}(M_h)> = 
    
    Parameters
    ----------
    Mmin : float
    sigma : float
    M0 : float
    M1 : float
    alpha : float
    """

    def __init__ ( self,
                   Mmin = 1.e+10,
                   sigma = 1.0,
                   M0 = 1.e+10,
                   M1 = 1.e+10,
                   alpha = 1.0 ) :

        self.Mmin = Mmin
        self.sigma = sigma
        self.M0 = M0
        self.M1 = M1
        self.alpha = alpha

    def Pcen ( self, Mh, **kwargs ) :
        """ Computes

        .. math::
      
          <N_\\text{cen}(M_h)> = 

        Parameters
        ----------
        Mh : scalar or array-like
          Mass of the halo
    
        Keyword Arguments
        -----------------
        Mmin : float
        sigma : float
        M0 : float
        M1 : float
        alpha : float
        
        Returns
        -------
        : scalar or array-like
        """
        self.__dict__.update( kwargs )
        ret = numpy.zeros_like( Mh )
        ww = Mh > 0.0
        ret[ww] = 0.5 * (
            1. + erf(
                (numpy.log10( Mh[ww] ) - numpy.log10( self.Mmin )) /
                self.sigma
            )
        )
        return ret

    def Psat ( self, Mh, **kwargs ) :
        """

        .. math::
        
          <N_\\text{sat}(M_h)> = 
        
        Parameters
        ----------
        Mh : scalar or array-like
          Mass of the halo
        
        Keyword Arguments
        -----------------
        Mmin : float
        sigma : float
        M0 : float
        M1 : float
        alpha : float
        
        Returns
        -------
        : scalar or array-like
        """
        self.__dict__.update( kwargs )
        psat = ( ( Mh - self.M0 ) / self.M1 )
        ret = numpy.zeros_like(Mh)
        ww = (psat>0.0)
        ret[ww] = self.Pcen( Mh[ww] ) * psat[ww]**self.alpha 
        return ret

    def get_hosts ( self, cat, rng = None, kw_rng = {}, **kwargs ) :
        """Applies the HOD parameterisation to find hosts among the
        subhaloes of a catalogues.
        
        Patameters
        ----------
        cat : scampy.catalogue.catalogue instance
        rng : 
        kw_rng : dict
        
        Keyword arguments
        -----------------
        Mmin : float
        sigma : float
        M0 : float
        M1 : float
        alpha : float
        
        Returns
        -------
        cen_gxy : int 1d-array
        sat_gxy : int 1d-array
        """

        # Check input validity
        if not isinstance( cat, catalogue ) :
            raise ValueError( 'Argument cat should be an instance '
                              'of type scampy.catalogue.catalogue' )

        # Update HOD parameters
        self.__dict__.update( **kwargs )

        # Build random number generator
        if rng is None :
            rng = numpy.random.default_rng(**kw_rng)
        
        # Central probability
        pcg = self.Pcen( cat.haloes.Mhalo )
    
        # Satellites probability
        psg = numpy.ones(cat.haloes.size)
        Nsg = self.Psat( cat.haloes.Mhalo )
        Nsh = cat.Nsat()
        wws = ( 0.0 < Nsh ) & ( Nsh >= Nsg )
        psg[wws] = Nsg[wws] / Nsh[wws]
        
        # Find centrals mask
        cen_idx = cat.haloes.centrals()
        cen_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        cen_gxy[cen_idx] = rng.binomial(1, pcg[cat.subhaloes.Parent[cen_idx]])
        
        # Find satellites mask
        sat_idx = cat.haloes.satellites()
        sat_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        sat_gxy[sat_idx] = rng.binomial(1, psg[cat.subhaloes.Parent[sat_idx]])
    
        return cen_gxy, sat_gxy
        

    
