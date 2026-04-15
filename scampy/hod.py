"""Halo Occupation Distribution (HOD) models.

Provides classes for the 5-parameter HOD prescription of
Zheng, Coil & Zehavi (2007) and Zheng et al. (2009), as described in
Sec. 2.1 and Eqs. 19–20 of Ronconi et al. (2020), together with a
redshift-dependent extension.  All classes expose ``Pcen(M)`` and
``Psat(M)`` methods compatible with :class:`~scampy.halo.model.halo_model`.
"""

import numpy
from scipy.special import erf

from scampy.catalogue import catalogue
from scampy.utilities.functions import repeated_mask

def get_hosts ( cat, Pcen, Psat,
                method = 'binomial',
                kw_pcen = {}, 
                kw_psat = {}, 
                rng = None, kw_rng = {} ) :
    """Select host sub-haloes from a catalogue using arbitrary occupation functions.

    Low-level helper that applies user-supplied central and satellite
    occupation probabilities to a catalogue without requiring a full HOD
    class instance.  For the standard HOD workflow prefer the
    :meth:`~scampy.hod.HOD.get_hosts` method of :class:`HOD`.

    Parameters
    ----------
    cat : scampy.catalogue.catalogue
        Halo/sub-halo catalogue to populate.
    Pcen : callable or array-like
        Central occupation probability.  If callable it is called as
        ``Pcen(Mhalo, **kw_pcen)`` and must return an array of
        probabilities in ``[0, 1]``.  If array-like it must have length
        ``cat.haloes.size``.
    Psat : callable or array-like
        Desired mean number of satellites per halo.  If callable it is
        called as ``Psat(Mhalo, **kw_psat, **kw_pcen)`` and must return
        an array of non-negative values.  If array-like it must have
        length ``cat.haloes.size``.
    kw_pcen : dict, optional
        Keyword arguments forwarded to ``Pcen`` when callable.
    kw_psat : dict, optional
        Keyword arguments forwarded to ``Psat`` when callable.
    rng : numpy.random.Generator or None, optional
        Random number generator.  If ``None`` (default) a new generator
        is created via ``numpy.random.default_rng(**kw_rng)``.
    kw_rng : dict, optional
        Keyword arguments forwarded to ``numpy.random.default_rng`` when
        ``rng`` is ``None``.

    Returns
    -------
    cen_gxy : ndarray of bool
        Boolean mask over ``cat.subhaloes`` selecting central galaxies.
    sat_gxy : ndarray of bool
        Boolean mask over ``cat.subhaloes`` selecting satellite galaxies.
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
    """5-parameter Halo Occupation Distribution (Zheng et al. 2007).

    Computes the average number of central and satellite galaxies hosted
    by a halo of mass :math:`M_h` (Eqs. 19–20 of Ronconi et al. 2020):

    .. math::

        \\langle N_\\mathrm{cen}(M_h)\\rangle =
        \\frac{1}{2}\\left[1 + \\mathrm{erf}\\!
        \\left(\\frac{\\log M_h - \\log M_\\mathrm{min}}{\\sigma}
        \\right)\\right],

    .. math::

        \\langle N_\\mathrm{sat}(M_h)\\rangle =
        \\langle N_\\mathrm{cen}(M_h)\\rangle
        \\left(\\frac{M_h - M_0}{M_1}\\right)^\\alpha.

    The satellite term is conditioned on the central (i.e. a halo must
    host a central to host satellites).  For the unconditioned variant
    see :class:`HOD_unconditioned_sat`.

    Parameters
    ----------
    Mmin : float, optional
        Characteristic minimum halo mass :math:`M_\\mathrm{min}` for
        hosting a central galaxy (default: ``1e10``).
    sigma : float, optional
        Width :math:`\\sigma` of the central occupation transition in
        :math:`\\log_{10} M` (default: ``1.0``).
    M0 : float, optional
        Satellite mass offset :math:`M_0` (default: ``1e10``).
    M1 : float, optional
        Satellite mass normalisation :math:`M_1` (default: ``1e10``).
    alpha : float, optional
        Satellite power-law slope :math:`\\alpha` (default: ``1.0``).
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

    def __str__ ( self ) :
        sstr = ''
        for k, v in self.__dict__.items() :
            sstr += f'\t{k:s} = {v:.6e},\n'
        return f'{type(self).__name__}(\n{sstr})'

    def __repr__ ( self ) :
        sstr = ''
        for k, v in self.__dict__.items() :
            sstr += f'\t{k:s} = {v},\n'
        return f'{type(self).__name__}(\n{sstr})'

    def Pcen ( self, Mh, **kwargs ) :
        """Average number of central galaxies per halo of mass :math:`M_h`.

        .. math::

            \\langle N_\\mathrm{cen}(M_h)\\rangle =
            \\frac{1}{2}\\left[1 + \\mathrm{erf}\\!
            \\left(\\frac{\\log M_h - \\log M_\\mathrm{min}}{\\sigma}
            \\right)\\right].

        Parameters
        ----------
        Mh : scalar or array-like
            Halo mass :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.

        Keyword Arguments
        -----------------
        Mmin, sigma, M0, M1, alpha : float
            Override the instance parameters for this call.

        Returns
        -------
        ndarray
            Values in :math:`[0, 1]`, same shape as ``Mh``.
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
        """Average number of satellite galaxies per halo of mass :math:`M_h`.

        .. math::

            \\langle N_\\mathrm{sat}(M_h)\\rangle =
            \\langle N_\\mathrm{cen}(M_h)\\rangle
            \\left(\\frac{M_h - M_0}{M_1}\\right)^\\alpha.

        Returns zero when :math:`M_h \\leq M_0`.

        Parameters
        ----------
        Mh : scalar or array-like
            Halo mass :math:`M_h` in :math:`[M_\\odot\\,h^{-1}]`.

        Keyword Arguments
        -----------------
        Mmin, sigma, M0, M1, alpha : float
            Override the instance parameters for this call.

        Returns
        -------
        ndarray
            Non-negative values, same shape as ``Mh``.
        """
        self.__dict__.update( kwargs )
        psat = ( ( Mh - self.M0 ) / self.M1 )
        ret = numpy.zeros_like(Mh)
        ww = (psat>0.0)
        ret[ww] = self.Pcen( Mh[ww] ) * psat[ww]**self.alpha
        return ret

    def get_hosts (
            self, cat,
            smask = None, method = 'binomial',
            rng = None, kw_rng = {}, **kwargs
    ) :
        """Apply the HOD to select host sub-haloes from a catalogue.

        Parameters
        ----------
        cat : scampy.catalogue.catalogue
            Halo/sub-halo catalogue to populate.
        smask : ndarray of bool or None, optional
            Boolean mask of length ``cat.subhaloes.size`` pre-selecting
            sub-halo candidates.  ``True`` = include.  If ``None``
            (default) all sub-haloes are considered.
        method : {'binomial', 'rankorder'}, optional
            Galaxy assignment method (default: ``'binomial'``):

            * ``'binomial'`` — centrals drawn from a Binomial(1, Pcen)
              distribution; satellites drawn from a Binomial(1, Psat/Nsh)
              distribution.
            * ``'rankorder'`` — centrals as above; satellites chosen by
              keeping the first :math:`N_\\mathrm{sat}` sub-haloes in the
              input order (assumes sub-haloes are sorted by some property,
              e.g. mass).
        rng : numpy.random.Generator or None, optional
            Random number generator.  If ``None`` (default) one is created
            via ``numpy.random.default_rng(**kw_rng)``.
        kw_rng : dict, optional
            Keyword arguments forwarded to ``numpy.random.default_rng``
            (default: ``{}``).

        Keyword Arguments
        -----------------
        Mmin, sigma, M0, M1, alpha : float
            Override the instance HOD parameters for this call.

        Returns
        -------
        cen_gxy : ndarray of bool
            Boolean mask of length ``cat.subhaloes.size``; ``True`` where
            a central galaxy is assigned.
        sat_gxy : ndarray of bool
            Boolean mask of length ``cat.subhaloes.size``; ``True`` where
            a satellite galaxy is assigned.
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
        Nsh = cat.Nsat().astype(int)
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
                idmask, countsmask = numpy.unique(
                    cat.subhaloes.Parent[(~smask)&sat_mask], return_counts=True
                )
                Nsg_mask[idmask] += countsmask
            # Guarantee not more than the total available number of satellites is chosen:
            Nsg_tot = numpy.min([Nsg+Nsg_mask, Nsh], axis=0).astype(int)
            # Generate new mask accounting for availability of satellite hosts
            sat_gxy[sat_idx] = repeated_mask( Nsg_tot, Nsh-Nsg_tot )
            sat_gxy &= smask

        return cen_gxy, sat_gxy
        
class HOD_unconditioned_sat () :
    """5-parameter HOD with unconditioned satellite occupation.

    Identical to :class:`HOD` except that the satellite term is *not*
    conditioned on the central:

    .. math::

        \\langle N_\\mathrm{cen}(M_h)\\rangle =
        \\frac{1}{2}\\left[1 + \\mathrm{erf}\\!
        \\left(\\frac{\\log M_h - \\log M_\\mathrm{min}}{\\sigma}
        \\right)\\right],

    .. math::

        \\langle N_\\mathrm{sat}(M_h)\\rangle =
        \\left(\\frac{M_h - M_0}{M_1}\\right)^\\alpha.

    Parameters
    ----------
    Mmin : float, optional
        Characteristic minimum halo mass :math:`M_\\mathrm{min}`
        (default: ``1e10``).
    sigma : float, optional
        Width of the central occupation transition in
        :math:`\\log_{10} M` (default: ``1.0``).
    M0 : float, optional
        Satellite mass offset :math:`M_0` (default: ``1e10``).
    M1 : float, optional
        Satellite mass normalisation :math:`M_1` (default: ``1e10``).
    alpha : float, optional
        Satellite power-law slope :math:`\\alpha` (default: ``1.0``).
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
        """Average number of central galaxies per halo of mass :math:`M_h`.

        Same formula as :meth:`HOD.Pcen`.

        Parameters
        ----------
        Mh : scalar or array-like
            Halo mass in :math:`[M_\\odot\\,h^{-1}]`.

        Keyword Arguments
        -----------------
        Mmin, sigma, M0, M1, alpha : float
            Override the instance parameters for this call.

        Returns
        -------
        ndarray
            Values in :math:`[0, 1]`, same shape as ``Mh``.
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
        """Average number of satellite galaxies per halo of mass :math:`M_h`.

        .. math::

            \\langle N_\\mathrm{sat}(M_h)\\rangle =
            \\left(\\frac{M_h - M_0}{M_1}\\right)^\\alpha.

        Unlike :meth:`HOD.Psat`, this is *not* conditioned on the central
        occupation.  Returns zero when :math:`M_h \\leq M_0`.

        Parameters
        ----------
        Mh : scalar or array-like
            Halo mass in :math:`[M_\\odot\\,h^{-1}]`.

        Keyword Arguments
        -----------------
        Mmin, sigma, M0, M1, alpha : float
            Override the instance parameters for this call.

        Returns
        -------
        ndarray
            Non-negative values, same shape as ``Mh``.
        """
        self.__dict__.update( kwargs )
        psat = ( ( Mh - self.M0 ) / self.M1 )
        ret = numpy.zeros_like(Mh)
        ww = (psat>0.0)
        ret[ww] = psat[ww]**self.alpha
        return ret

    def get_hosts ( self, cat, method = 'binomial', rng = None, kw_rng = {}, **kwargs ) :
        """Apply the HOD to select host sub-haloes from a catalogue.

        Parameters
        ----------
        cat : scampy.catalogue.catalogue
            Halo/sub-halo catalogue to populate.
        method : {'binomial', 'rankorder'}, optional
            Galaxy assignment method (default: ``'binomial'``):

            * ``'binomial'`` — centrals drawn from a Binomial(1, Pcen)
              distribution; satellites drawn from a Binomial(1, Psat/Nsh)
              distribution.
            * ``'rankorder'`` — centrals as above; satellites chosen by
              keeping the first :math:`N_\\mathrm{sat}` sub-haloes in the
              input order (assumes sub-haloes are sorted by some property,
              e.g. mass).
        rng : numpy.random.Generator or None, optional
            Random number generator.  If ``None`` (default) one is created
            via ``numpy.random.default_rng(**kw_rng)``.
        kw_rng : dict, optional
            Keyword arguments forwarded to ``numpy.random.default_rng``
            (default: ``{}``).

        Keyword Arguments
        -----------------
        Mmin, sigma, M0, M1, alpha : float
            Override the instance HOD parameters for this call.

        Returns
        -------
        cen_gxy : ndarray of bool
            Boolean mask of length ``cat.subhaloes.size``; ``True`` where
            a central galaxy is assigned.
        sat_gxy : ndarray of bool
            Boolean mask of length ``cat.subhaloes.size``; ``True`` where
            a satellite galaxy is assigned.
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
    """Redshift-dependent 5-parameter HOD.

    Extends :class:`HOD` by allowing each parameter to evolve linearly
    with redshift in log-space:

    .. math::

        \\log_{10} X(z) = \\ell_X + \\ell_{X,b}\\,z,

    so that :math:`X(z) = 10^{\\ell_X + \\ell_{X,b}\\,z}`.
    Setting all :math:`\\ell_{X,b} = 0` recovers a redshift-independent
    :class:`HOD`.

    The occupation functions follow the same :class:`HOD` formulae with
    the redshift-evaluated parameters.

    Parameters
    ----------
    redshift : float or array-like
        Redshift(s) at which to evaluate the HOD parameters.
    lMmin : float, optional
        :math:`\\log_{10} M_\\mathrm{min}` at :math:`z=0`
        (default: ``10``).
    lMmin_b : float, optional
        Redshift slope of :math:`\\log_{10} M_\\mathrm{min}`
        (default: ``0.0``).
    lsigma : float, optional
        :math:`\\log_{10} \\sigma` at :math:`z=0` (default: ``0.0``).
    lsigma_b : float, optional
        Redshift slope of :math:`\\log_{10} \\sigma` (default: ``0.0``).
    lM0 : float, optional
        :math:`\\log_{10} M_0` at :math:`z=0` (default: ``10``).
    lM0_b : float, optional
        Redshift slope of :math:`\\log_{10} M_0` (default: ``0.0``).
    lM1 : float, optional
        :math:`\\log_{10} M_1` at :math:`z=0` (default: ``10``).
    lM1_b : float, optional
        Redshift slope of :math:`\\log_{10} M_1` (default: ``0.0``).
    lalpha : float, optional
        :math:`\\log_{10} \\alpha` at :math:`z=0` (default: ``0.0``).
    lalpha_b : float, optional
        Redshift slope of :math:`\\log_{10} \\alpha` (default: ``0.0``).
    """

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
        """Update redshift and/or free parameters, then recompute ``self.params``.

        Parameters
        ----------
        redshift : float or array-like or None, optional
            New redshift(s).  If ``None`` the current redshift is kept.
        **kwargs
            Any subset of the free parameters
            (``lMmin``, ``lMmin_b``, ``lsigma``, ``lsigma_b``,
            ``lM0``, ``lM0_b``, ``lM1``, ``lM1_b``,
            ``lalpha``, ``lalpha_b``).
            Unknown keys are silently ignored.
        """
        
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

