"""Handling catalogues of haloes and subhaloes
"""

# External includes
import numpy

# Internal includes
from .utilities.base_classes import fixedSizeDict

############################################################################################
# Kernel Catalogue class, inheriting from maskedDict

class kernelCat ( fixedSizeDict ) :
    """Base catalogue container enforcing mandatory coordinate fields.

    Extends :class:`~scampy.utilities.base_classes.fixedSizeDict` by
    requiring that the fields listed in the class attribute ``fields``
    (``'X'``, ``'Y'``, ``'Z'``) are always present.  Subclasses extend
    this list with additional mandatory fields.

    Parameters
    ----------
    *args, **kwargs
        Forwarded to :class:`~scampy.utilities.base_classes.fixedSizeDict`.
        Must include at least the mandatory fields in ``fields``.
    """

    fields = [ 'X', 'Y', 'Z' ]
    _internal = fixedSizeDict._internal + [ 'fields' ]

    def __init__ ( self, *args, **kwargs ) :
        super().__init__( *args, **kwargs )
        for key in type(self).fields :
            if key not in self.keys() :
                raise AttributeError( f"Mandatory field '{key}' not provided." )

    def coord ( self ) :
        """Return the 3D Cartesian coordinates as a ``(Nobj, 3)`` array."""
        return numpy.array( [ self.X, self.Y, self.Z ] ).T

    def return_sample ( self, X_fields = [], Y_fields = [] ) :
        """Return two arrays built from selected fields.

        Parameters
        ----------
        X_fields : list of str, optional
            Field names to stack into the first output array.
        Y_fields : list of str, optional
            Field names to stack into the second output array.

        Returns
        -------
        X : ndarray, shape ``(Nobj, len(X_fields))``
        Y : ndarray, shape ``(Nobj, len(Y_fields))``
        """

        return ( numpy.array( [ self[ f ] for f in X_fields ] ).T,
                 numpy.array( [ self[ f ] for f in Y_fields ] ).T )

    def __getstate__ ( self ) :
        return self.__dict__.copy()
    
############################################################################################
# Halo Catalogue class, inheriting from kernelCat

class haloCat ( kernelCat ) :
    """Catalogue of host haloes.

    Extends :class:`kernelCat` with the mandatory halo-specific fields:

    * ``Mhalo`` — halo mass :math:`[M_\\odot\\,h^{-1}]`
    * ``Rhalo`` — halo radius :math:`[h^{-1}\\,\\mathrm{Mpc}]`
    * ``firstSub`` — index of the first sub-halo in the sub-halo table
    * ``numSubs`` — total number of sub-haloes in the halo
    """

    fields = kernelCat.fields + [ 'Mhalo', 'Rhalo', 'firstSub', 'numSubs' ]

    def centrals ( self ) :
        """ Return list of indices in the sub-halo table locating the central 
        sub-halo of each halo. (takes no arguments)

        See Also
        ---------
        satellites : same for satellite sub-haloes
        """
        
        return self['firstSub'][self['numSubs'] > 0]

    
    def satellites ( self ) :
        """ Return list of indices in the sub-halo table locating the satellite
        sub-haloes of each halo. (takes no arguments)

        Warning
        -------
        Computes the slices only for un-masked haloes

        See Also
        ---------
        centrals : same for central sub-haloes
        """

        hassat = ( self['numSubs'] > 1 )
        sattot = self['numSubs'][hassat].sum() - hassat.sum()
        sat = -1 * numpy.ones(int(sattot), dtype = int)
        offset = 0
        for first, nsub in zip( self['firstSub'][hassat],
                                self['numSubs'][hassat] ) :
            append = nsub-1
            sat[offset:offset+append] = numpy.arange(first+1,
                                                     first+nsub,
                                                     dtype=int)
            offset += append

        if numpy.any(sat < 0) :
            raise RuntimeError( "Something went wrong in returning satellites" )
        
        return sat
        
############################################################################################
# Sub-Halo Catalogue class, inheriting from kernelCat
        
class subhaloCat ( kernelCat ) :
    """Catalogue of sub-haloes.

    Extends :class:`kernelCat` with the mandatory sub-halo-specific fields:

    * ``Msubh`` — sub-halo mass :math:`[M_\\odot\\,h^{-1}]`
    * ``Parent`` — index of the parent halo in the :class:`haloCat` table
    """

    fields = kernelCat.fields + [ 'Msubh', 'Parent' ]

    def Nsub ( self, mask = None ) :
        """ Function returning the number of un-masked sub-haloes in not-empty parent halo
        A halo is considered 'not-empty' when it has at least one sub-halo which is un-masked
        
        Parameters
        ----------
        mask : bool or iterable or None
          It accepts an iterable with lenght equal to the size of the sub-halo catalogue
          or 'None' to use all the objects in the sub-halo catalogue.

        Returns
        -------
        Nsub : numpy array
          array with the number of sub-haloes in each not-empty halo
        hidx : numpy array
          indices in the halo catalogue of the not-empty haloes 
        
        See Also
        --------
        catalogue.Ncen : function returning the valid central galaxy number
        catalogue.Nsat : function returning the valid satellite galaxies number
        catalogue.NcenNsat : function returning a tuple with the number of valid galaxies
        """

        if mask is not None :
            if len( mask ) != self.size :
                raise AttributeError(
                    "the provided mask should have the same size of the sub-halo catalogue: "
                    f"mask.size = {len(mask)} != {self.size} = catalogue.size")
        else :
            mask = numpy.ones( self.size, dtype = bool )
            
        hidx, Nobj = numpy.unique( self['Parent'][mask],
                                   return_counts = True )
        
        return Nobj, hidx
    
############################################################################################
# Catalogue class

class catalogue () :
    """ class catalogue for managing halo/subhalo hierarchies
    
    Parameters
    ----------
    haloes : haloCat object

    subhaloes : subhaloCat object
    
    **kwargs : 
      any argument passed as keyword argument is considered as
      a metadatum of the catalogue 
      (as long as it does not already exist in the internal
       dictionary of the class).
    """

    def __init__ ( self,
                   haloes,
                   subhaloes,
                   Lbox, 
                   **kwargs ) :

        self.haloes = haloes
        self.subhaloes = subhaloes
        self.Lbox = Lbox
        self.volume = self.Lbox**3
        self.meta = { 'Lbox' : self.Lbox, 'volume' : self.volume }
        for k,v in kwargs.items() :
            if k not in self.__dict__.keys() :
                self.meta[ k ] = v
            else :
                raise AttributeError( f"duplicate attribute {k}." )
    
    def save ( self, outPath, hard = False ) :
        """ Stores current state of catalogue in an hdf5 file
        
        Parameters
        ----------
        outPath : string
          '/path/to/name/of/file' the function will append the '.hdf5'
          extension to this path then check whether the file already exists
        hard : bool
          if True eventually overwrites an already existing file
          with the same name (default = False)

        See Also
        --------
        load : function for building catalogue from hdf5 file
        """
        import os, h5py

        outfile = f"{outPath}.hdf5"
        if not hard and os.path.isfile( outfile ) :
            raise AttributeError( f"file\n{outfile}\nalready exist, "
                                  "to overwrite it set argument hard=True." )
        
        with h5py.File( outfile, "w" ) as f :
            
            # store meta-data
            for k,v in self.meta.items() :
                f.attrs[k] = v
            
            # store halo catalogue
            halo_group = f.create_group('haloes')
            for kh in self.haloes.keys() :
                try :
                    halo_group.create_dataset( kh, data = self.haloes[kh] )
                except Exception as err :
                    print( f"- Skipping dataset '{kh}' of group 'haloes' as it raises "
                           f"exception {type(err).__name__} with message:\n  '{err}'" )
            
            # store subhalo catalogue
            subh_group = f.create_group('subhaloes')
            for ks in self.subhaloes.keys() :
                try :
                    subh_group.create_dataset( ks, data = self.subhaloes[ks] )
                except Exception as err :
                    print( f"- Skipping dataset '{ks}' of group 'subhaloes' as it raises "
                           f"exception {type(err).__name__} with message:\n  '{err}'" )

        print( f"Wrote on file {outfile}" )
        return;
            
    def Nsub ( self, mask_subhaloes = None ) :
        """ Function returning the number of un-masked sub-haloes in not-empty parent halo
        A halo is considered 'not-empty' when it has at least one sub-halo which is un-masked
        
        Parameters
        ----------
        mask_subhaloes : bool or iterable or None
          It accepts an iterable with lenght equal to the size of the sub-halo catalogue
          or 'None' to use all the objects in the sub-halo catalogue.

        Returns
        -------
        Nsub : numpy array
          array with the number of sub-haloes in each not-empty halo
        hidx : numpy array
          indices in the halo catalogue of the not-empty haloes 
        
        See Also
        --------
        Ncen : function returning the valid central galaxy number
        Nsat : function returning the valid satellite galaxies number
        NcenNsat : function returning a tuple with the number of valid galaxies
        """

        if mask_subhaloes is not None :
            if len( mask_subhaloes ) != self.subhaloes.size :
                raise AttributeError(
                    "the provided mask should have the same size of the sub-halo catalogue: "
                    f"mask_subhaloes.size = {len(mask_subhaloes)} != {self.subhaloes.size} "
                    "= catalogue.subhaloes.size")
        else :
            mask = numpy.ones( self.subhaloes.size, dtype = bool )
            
        hidx, Nobj = numpy.unique( self.subhaloes['Parent'][mask_subhaloes],
                                   return_counts = True )
        
        return Nobj, hidx

    def Ncen ( self, mask = None ) :
        """ Function returning the number of un-masked central galaxies
        
        Parameters
        ----------
        mask : bool or iterable or None
          Accepts an iterable with lenght equal to the lenght of the sub-halo catalogue mask
          or 'None' to use all the objects in the sub-halo catalogue.

        Returns
        -------
        : numpy array
          array of shape (haloes.size(),) with the number of centrals for each halo
        
        See Also
        --------
        Nsub : function returning the number of un-masked sub-haloes 
               in not-empty haloes
        Nsat : function returning the valid satellite galaxies number
        NcenNsat : function returning a tuple with the number of valid galaxies
        """

        if mask is not None :
            if len( mask ) != self.subhaloes.size :
                raise AttributeError(
                    "the provided mask should have the same size of the sub-halo catalogue: "
                    f"mask.size = {len(mask)} != {self.subhaloes.size} = catalogue.subhaloes.size"
                )
        else :
            mask = numpy.ones( self.subhaloes.size, dtype = bool )
        Nobj, hidx = self.Nsub( mask )
        
        Nc = numpy.zeros( (self.haloes.size,), dtype = float )
        numpy.put( Nc, hidx, (Nobj > 0).astype(float) )
        
        return Nc
        
    def Nsat ( self, mask = None ) :
        """ Function returning the number of un-masked satellite galaxies
        
        Parameters
        ----------
        add_mask : bool or iterable or None
          It accepts either a boolean (True = mask all the entries, 
          False = unmask all the entries),
          an iterable with lenght equal to the lenght of the sub-halo catalogue mask
          or 'None' to only use the default mask of the sub-halo catalogue.

        Returns
        -------
        : numpy array
          array of shape (haloes.size(),) with the number of satellites for each halo
        
        See Also
        --------
        Nsub : function returning the number of un-masked sub-haloes 
               in not-empty haloes
        Ncen : function returning the valid central galaxy number
        NcenNsat : function returning a tuple with the number of valid galaxies
        """

        if mask is not None :
            if len( mask ) != self.subhaloes.size :
                raise AttributeError(
                    "the provided mask should have the same size of the sub-halo catalogue: "
                    f"mask.size = {len(mask)} != {self.subhaloes.size} = catalogue.subhaloes.size")
        else :
            mask = numpy.ones( self.subhaloes.size, dtype = bool )            
        Nobj, hidx = self.Nsub( mask )
        
        Ns = numpy.zeros( (self.haloes.size, ), dtype = float )
        numpy.put( Ns, hidx, Nobj - 1 )
        Ns[ Ns < 0 ] = 0
        
        return Ns

    def NcenNsat ( self, mask = None ) :
        """ Function returning a tuple with the number of un-masked galaxies 
        (centrals and satellites)
        
        Parameters
        ----------
        add_mask : bool or iterable or None
          It accepts either a boolean (True = mask all the entries, 
          False = unmask all the entries),
          an iterable with lenght equal to the lenght of the sub-halo catalogue mask
          or 'None' to only use the default mask of the sub-halo catalogue.

        Returns
        -------
        Nc : numpy array
          array of shape (haloes.size(),) with the number of centrals for each halo
        Ns : numpy array
          array of shape (haloes.size(),) with the number of satellites for each halo
        
        See Also
        --------
        Nsub : function returning the number of un-masked sub-haloes 
               in not-empty haloes
        Ncen : function returning the valid central galaxy number
        Nsat : function returning the valid satellite galaxies number
        """

        if mask is not None :
            if len( mask ) != self.subhaloes.size :
                raise AttributeError(
                    "the provided mask should have the same size of the sub-halo catalogue: "
                    f"mask.size = {len(mask)} != {self.subhaloes.size} = catalogue.subhaloes.size")
        else :
            mask = numpy.ones( self.subhaloes.size, dtype = bool )
        Nobj, hidx = self.Nsub( mask )

        Nc, Ns = numpy.zeros( ( 2, self.haloes.size ), dtype = float )

        numpy.put( Nc, hidx, (Nobj > 0).astype(float) )
        numpy.put( Ns, hidx, Nobj - 1 ); Ns[ Ns < 0 ] = 0
        
        return Nc, Ns

    def centrals ( self, hmask = None, smask = None ) :
        """Return indices in the sub-halo table locating the central sub-halo
        of each halo.

        Parameters
        ----------
        hmask : ndarray of bool or None, optional
            Boolean mask of length ``haloes.size``; only centrals whose
            parent halo passes this mask are returned.  If ``None``
            (default) all haloes are included.
        smask : ndarray of bool or None, optional
            Boolean mask of length ``subhaloes.size``; only sub-haloes
            passing this mask are considered as centrals.  If ``None``
            (default) all sub-haloes are included.

        Returns
        -------
        cen : ndarray of int
            Indices into the sub-halo catalogue of the selected central
            sub-haloes.

        See Also
        --------
        satellites : same for satellite sub-haloes
        """
        
        # Halo mask
        if hmask is not None and hmask.size != self.haloes.size :
            raise RuntimeError('halo mask should have the same size '
                               'of the halo catalogue')
        if hmask is None :
            hmask = numpy.ones(self.haloes.size, dtype = bool)

        # Sub-halo mask
        if smask is not None and smask.size != self.subhaloes.size :
            raise RuntimeError('subhalo mask should have the same size '
                               'of the subhalo catalogue')
        if smask is None :
            smask = numpy.ones(self.subhaloes.size, dtype = bool)

        # Get all the centrals
        cen = self.haloes.centrals()

        # Retain only un-masked subhaloes
        cen = cen[smask[cen]]

        # Get all parents
        par = self.subhaloes.Parent[cen]

        # Retain only subhaloes in un-masked haloes
        cen = cen[hmask[par]]

        return cen    

    def satellites ( self, hmask = None, smask = None ) :
        """Return indices in the sub-halo table locating the satellite
        sub-haloes of each halo.

        Parameters
        ----------
        hmask : ndarray of bool or None, optional
            Boolean mask of length ``haloes.size``; only satellites whose
            parent halo passes this mask are returned.  If ``None``
            (default) all haloes are included.
        smask : ndarray of bool or None, optional
            Boolean mask of length ``subhaloes.size``; only sub-haloes
            passing this mask are considered as satellites.  If ``None``
            (default) all sub-haloes are included.

        Returns
        -------
        sat : ndarray of int
            Indices into the sub-halo catalogue of the selected satellite
            sub-haloes.

        Warning
        -------
        Computes the slices only for un-masked haloes.

        See Also
        --------
        centrals : same for central sub-haloes
        """
        
        # Halo mask
        if hmask is not None and hmask.size != self.haloes.size :
            raise RuntimeError('halo mask should have the same size '
                               'of the halo catalogue')
        if hmask is None :
            hmask = numpy.ones(self.haloes.size, dtype = bool)

        # Sub-halo mask
        if smask is not None and smask.size != self.subhaloes.size :
            raise RuntimeError('subhalo mask should have the same size '
                               'of the subhalo catalogue')
        if smask is None :
            smask = numpy.ones(self.subhaloes.size, dtype = bool)
        
        # Get all the satellites
        sat = self.haloes.satellites()

        # Retain only un-masked subhaloes
        sat = sat[smask[sat]]

        # Get all parents
        par = self.subhaloes.Parent[sat]

        # Retain only subhaloes in un-masked haloes
        sat = sat[hmask[par]]
                
        return sat

############################################################################################
