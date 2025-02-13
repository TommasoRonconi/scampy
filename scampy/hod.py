import numpy
from scipy.special import erf

from scampy.catalogue import catalogue
from scampy.utilities.functions import repeated_mask

def get_hosts ( cat, Pcen, Psat,
                method = 'binomial',
                kw_pcen = {}, 
                kw_psat = {}, 
                rng = None, kw_rng = {} ) :
    """
    """
    from scampy.catalogue import catalogue

    # Check input validity
    if not isinstance( cat, catalogue ) :
        raise ValueError( 'Argument cat should be an instance of type scampy.catalogue' )
    
    # Build random number generator
    if rng is None :
        rng = numpy.random.default_rng(**kw_rng)
        
    # Central probability
    pcg = numpy.ones(cat.haloes.size)
    if hasattr(Pcen, '__call__') :
        pcg = Pcen( cat.haloes.Mhalo, **kw_pcen )
    elif hasattr( Pcen, '__len__' ) :
        if len(Pcen) != cat.haloes.size :
            raise ValueError( 'len(Pcen) != haloes.size' )
        pcg = Pcen
    else :
        raise ValueError( 'wrong value passed to argument Pcen' )
        
    # Find centrals mask
    cen_idx = cat.haloes.centrals()
    cen_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
    cen_gxy[cen_idx] = rng.binomial(1, pcg[cat.subhaloes.Parent[cen_idx]])
    
    # Satellites probability
    psg = numpy.ones(cat.haloes.size)
    Nsh = cat.Nsat()
    Nsg = numpy.array(Nsh) # copy the array
    if hasattr(Psat, '__call__') :
        Nsg = Psat( cat.haloes.Mhalo, **kw_psat, **kw_pcen )
        if method == 'rankorder' :
            Nsg = rng.poisson(Nsg)
        wws = ( 0.0 < Nsh ) & ( Nsh > Nsg )
        if method == 'binomial' :
            psg[wws] = Nsg[wws] / Nsh[wws]
    elif hasattr( Psat, '__len__' ) :
        if len(Psat) != cat.haloes.size :
            raise ValueError( 'len(Psat) != haloes.size' )
        if method == 'binomial' :
            psg = Psat
        if method == 'rankorder' :
            Nsg = rng.poisson(Psat*Nsh)
    else :
        raise ValueError( 'wrong value passed to argument Psat' )
        
    # Find satellites mask
    sat_idx = cat.haloes.satellites()
    sat_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
    if method == 'binomial' :
        sat_gxy[sat_idx] = rng.binomial(1, psg[cat.subhaloes.Parent[sat_idx]])
    if method == 'rankorder' :
        sat_gxy[sat_idx] = repeated_mask(Nsg.astype(int), (Nsh-Nsg).astype(int))
    
    return cen_gxy, sat_gxy

# def get_hosts ( cat, Pcen, Psat, 
#                 kw_pcen = {}, 
#                 kw_psat = {}, 
#                 rng = None, kw_rng = {} ) :
#     """
#     """

#     # Check input validity
#     if not isinstance( cat, catalogue ) :
#         raise ValueError( 'Argument cat should be an instance of type scampy.catalogue' )
    
#     # Build random number generator
#     if rng is None :
#         rng = np.random.default_rng(**kw_rng)
        
#     # Central probability
#     pcg = np.ones(self.haloes.size)
#     if hasattr(Pcen, '__call__') :
#         pcg = Pcen( self.haloes.Mhalo, **kw_pcen )
#     elif hasattr( Pcen, '__len__' ) :
#         if len(Pcen) != self.haloes.size :
#             raise ValueError( 'len(Pcen) != haloes.size' )
#         pcg = Pcen
#     else :
#         raise ValueError( 'wrong value passed to argument Pcen' )
    
#     # Satellites probability
#     psg = np.ones(self.haloes.size)
#     if hasattr(Psat, '__call__') :
#         Nsg = Psat( self.haloes.Mhalo, **kw_psat, **kw_pcen )
#         Nsh = self.Nsat()
#         wws = ( 0.0 < Nsh ) & ( Nsh > Nsg )
#         psg[wws] = Nsg[wws] / Nsh[wws]
#     elif hasattr( Psat, '__len__' ) :
#         if len(Psat) != self.haloes.size :
#             raise ValueError( 'len(Psat) != haloes.size' )
#         psg = Psat
#     else :
#         raise ValueError( 'wrong value passed to argument Psat' )
        
#     # Find centrals mask
#     cen_idx = self.haloes.centrals()
#     cen_gxy = np.zeros(self.subhaloes.size, dtype=bool)
#     cen_gxy[cen_idx] = rng.binomial(1, pcg[self.subhaloes.Parent[cen_idx]])
        
#     # Find satellites mask
#     sat_idx = self.haloes.satellites()
#     sat_gxy = np.zeros(self.subhaloes.size, dtype=bool)
#     sat_gxy[sat_idx] = rng.binomial(1, psg[self.subhaloes.Parent[sat_idx]])
    
#     return cen_gxy, sat_gxy


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

    def get_hosts ( self, cat, smask = None, method = 'binomial', rng = None, kw_rng = {}, **kwargs ) :
        """Applies the HOD parameterisation to find hosts among the
        subhaloes of a catalogues.

        Patameters
        ----------
        cat : scampy.catalogue.catalogue instance
        smask :
        method : str
        method for galaxy assignment:
        'binomial' = centrals are extracted from a binomial 
                    distribution and satellites from a Poisson distribution
        'rankorder' = centrals are extracted from a binomial distribution
                     while satellites are obtained by keeping the first Nsat
                     objects in the list, it assumes the satellites array is
                     ordered based on some property
                     (assumes the input array of satellites is mass-ordered)
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
        from scampy.catalogue import catalogue
        from scampy.utilities.functions import repeated_mask

        # Check input validity
        if not isinstance( cat, catalogue ) :
            raise ValueError( 'Argument cat should be an instance '
                              'of type scampy.catalogue.catalogue' )

        # Update HOD parameters
        self.__dict__.update( **kwargs )

        # Build random number generator
        if rng is None :
            rng = numpy.random.default_rng(**kw_rng)

        # Sub-halo mask
        if smask is not None and smask.size != cat.subhaloes.size :
            raise RuntimeError('subhalo mask should have the same size '
                               'of the subhalo catalogue')
        if smask is None :
            smask = numpy.ones(cat.subhaloes.size, dtype = bool)
        
        # divide mask between centrals and satellites 
        cen_mask = numpy.zeros_like(smask)
        cen_mask[cat.centrals()] = True
        sat_mask = (~cen_mask)
        
        # Central probability
        pcg = self.Pcen( cat.haloes.Mhalo )
        
        # Find centrals mask
        cen_idx = cat.centrals(smask=smask)
        cen_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        cen_gxy[cen_idx] = rng.binomial(1, pcg[cat.subhaloes.Parent[cen_idx]])

        # Satellites probability
        psg = numpy.ones(cat.haloes.size)
        Nsh = cat.Nsat()
        Nsg = self.Psat( cat.haloes.Mhalo )
        if method == 'binomial' :
            wws = ( 0.0 < Nsh ) & ( Nsh >= Nsg )
            psg[wws] = Nsg[wws] / Nsh[wws]
        if method == 'rankorder' :
            Nsg = rng.poisson(Nsg)

        # Find satellites mask
        sat_idx = cat.satellites()
        sat_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        if method == 'binomial' :
            sat_gxy[sat_idx] = rng.binomial(1, psg[cat.subhaloes.Parent[sat_idx]])
        if method == 'rankorder' :
            Nsg_mask = numpy.zeros_like(Nsg)
            if any((~smask)&sat_mask) :
                idmask, countsmask = numpy.unique(cat.subhaloes.Parent[(~smask)&sat_mask], return_counts=True)
                Nsg_mask[idmask] += countsmask
            sat_gxy[sat_idx] = repeated_mask((Nsg+Nsg_mask).astype(int), (Nsh-Nsg-Nsg_mask).astype(int))
            sat_gxy &= smask

        return cen_gxy, sat_gxy
        
class HOD_unconditioned_sat () :
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
        ret[ww] = psat[ww]**self.alpha 
        return ret

    def get_hosts ( self, cat, method = 'binomial', rng = None, kw_rng = {}, **kwargs ) :
        """Applies the HOD parameterisation to find hosts among the
        subhaloes of a catalogues.
        
        Patameters
        ----------
        cat : scampy.catalogue.catalogue instance
        method : str
           method for galaxy assignment:
           'binomial' = centrals are extracted from a binomial 
                        distribution and satellites from a Poisson distribution
           'rankorder' = centrals are extracted from a binomial distribution
                         while satellites are obtained by keeping the first Nsat
                         objects in the list, it assumes the satellites array is
                         ordered based on some property
                         (assumes the input array of satellites is mass-ordered)
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
        
        # Find centrals mask
        cen_idx = cat.haloes.centrals()
        cen_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        cen_gxy[cen_idx] = rng.binomial(1, pcg[cat.subhaloes.Parent[cen_idx]])

        # Satellites probability
        psg = numpy.ones(cat.haloes.size)
        Nsh = cat.Nsat()
        Nsg = self.Psat( cat.haloes.Mhalo )
        if method == 'binomial' :
            wws = ( 0.0 < Nsh ) & ( Nsh >= Nsg )
            psg[wws] = Nsg[wws] / Nsh[wws]
        if method == 'rankorder' :
            Nsg = rng.poisson(Nsg)

        # Find satellites mask
        sat_idx = cat.haloes.satellites()
        sat_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        if method == 'binomial' :
            sat_gxy[sat_idx] = rng.binomial(1, psg[cat.subhaloes.Parent[sat_idx]])
        if method == 'rankorder' :
            sat_gxy[sat_idx] = repeated_mask(Nsg.astype(int), (Nsh-Nsg).astype(int))

        return cen_gxy, sat_gxy
    
        # # Satellites probability
        # psg = numpy.ones(cat.haloes.size)
        # Nsh = cat.Nsat()
        # Nsg = self.Psat( cat.haloes.Mhalo )
        # wws = ( 0.0 < Nsh ) & ( Nsh >= Nsg )
        # psg[wws] = Nsg[wws] / Nsh[wws]
        
        # # Find satellites mask
        # sat_idx = cat.haloes.satellites()
        # sat_gxy = numpy.zeros(cat.subhaloes.size, dtype=bool)
        # sat_gxy[sat_idx] = rng.binomial(1, psg[cat.subhaloes.Parent[sat_idx]])
    
        # return cen_gxy, sat_gxy


class HOD_zdep () :
    
    zdep_parameters = ['Mmin', 'sigma', 'M0', 'M1', 'alpha']
    free_parameters = set([ 
        'lMmin', 'lMmin_b', 
        'lsigma', 'lsigma_b', 
        'lM0', 'lM0_b', 
        'lM1', 'lM1_b', 
        'lalpha', 'lalpha_b' 
    ])
    
    def __init__ ( self, 
                   redshift,
                   lMmin = +10,
                   lMmin_b = 0.0,
                   lsigma = 0.0,
                   lsigma_b = 0.0,
                   lM0 = +10,
                   lM0_b = 0.0,
                   lM1 = +10,
                   lM1_b = 0.0,
                   lalpha = 0.0,
                   lalpha_b = 0.0 ) :

        self.scalar = False
        self._set_z( redshift )
        
        self.lMmin = lMmin
        self.lMmin_b = lMmin_b
        self.lsigma = lsigma
        self.lsigma_b = lsigma_b
        self.lM0 = lM0
        self.lM0_b = lM0_b
        self.lM1 = lM1
        self.lM1_b = lM1_b
        self.lalpha = lalpha
        self.lalpha_b = lalpha_b
        self.set()
        
    def _set_z ( self, redshift ) :
        
        self.z = numpy.array(redshift)
        if self.z.ndim == 0 :
            self.z = self.z[None]
            self.scalar = True
        return;
        
    def set ( self, redshift = None, **kwargs ) :
        
        if redshift is not None :
            self._set_z( redshift )
                
        self.__dict__.update( {
            k : v for k, v in kwargs.items() if k in type(self).free_parameters
        } )
        self.params = dict(zip(
            type(self).zdep_parameters,
            10.**numpy.array([
                self.lMmin + self.lMmin_b * self.z,
                self.lsigma + self.lsigma_b * self.z,
                self.lM0 + self.lM0_b * self.z,
                self.lM1 + self.lM1_b * self.z,
                self.lalpha + self.lalpha_b * self.z,
            ])
        ))
        
        return;
    
    def _Pcen ( self, Mh, redshift = None, **kwargs ) :
        
        self.set( redshift, **kwargs )
        Mh = numpy.asarray( numpy.squeeze(Mh) )
        ret = numpy.zeros( shape = (Mh.size, self.z.size), dtype=float )
        ww = Mh > 0.0
        ret[ww] = 0.5 * (
            1. + erf(
                (numpy.log10( Mh[ww,numpy.newaxis] ) -
                 numpy.log10( self.params['Mmin'][numpy.newaxis,:] )) /
                self.params['sigma']
            )
        )
        return ret
    
    def Pcen ( self, Mh, redshift = None, **kwargs ) :
        return numpy.squeeze(self._Pcen( Mh, redshift, **kwargs ))
    
    def Psat ( self, Mh, redshift = None, **kwargs ) :
        
        self.set( redshift, **kwargs )
        Mh = numpy.asarray(numpy.squeeze(Mh))
        psat = ( 
            ( Mh[:,numpy.newaxis] - self.params['M0'][numpy.newaxis,:] ) / 
            self.params['M1'][numpy.newaxis,:] 
        )
        psat[~(psat>0.0)] = 0.0
        return numpy.squeeze( 
            self._Pcen( Mh ) * 
            psat**self.params['alpha'][numpy.newaxis,:] 
        )

