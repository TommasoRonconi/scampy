# external:
import numpy
import copy

# internal:
from scampy import gadget_file
from scampy.objects import *

##########################################################################################
# catalogue class

class catalogue () :
    """
    Class to handle catalogues of objects of type halo, host_halo, galaxy
    """
    
    def __init__ ( self, X = None, scale_lenght = 1.e-3, scale_mass = 1.e+10, boxsize = None ) :
        if X is None : 
            self.content = numpy.array([])
        else :
            self.content = X 
        self.scale_lenght = scale_lenght
        self.scale_mass = scale_mass
        self.boxsize = boxsize
            
        self.gadget = None

    def Nhost ( self, mask = None ) :
        """ Return the total number of host haloes (central + satellites)

        Parameters
        ----------
        mask : array-like
          Array mask for filtering the original catalogue

        Returns
        -------
        : int
        """

        if mask is None :
            return numpy.array( [ obj.Ncen + obj.Nsat for obj in self.content ] ).sum()
        else :
            return numpy.array( [ obj.Ncen + obj.Nsat for obj in self.content[ mask ] ] ).sum()            
    
    def set_content ( self, X ) :
        """ Add element(s) to the catalogue
        
        Parameters
        ----------
        X : array-like

        Returns
        -------
        None
        """
        self.content = numpy.append( self.content, X )
        return None

    def sub_sample ( self, nsample ) :
        """ Return a sub-sampled catalogue

        Parameters
        ----------
        nsample : int
          number of objects in the sub-sampled catalogue, must be in the interval [0.,Nhost)
          where Nhost is the number of sub-haloes of the current catalogue

        Returns
        -------
        : catalogue
          sub-sampled copy of the original catalogue        
        """

        # allocate memory for new catalogue content
        content = numpy.empty( self.content.shape, dtype = host_halo )

        # create sub-sampling mask
        mask = numpy.full( self.Nhost(), False, dtype = bool )
        mask[ :nsample ] = True
        numpy.random.shuffle( mask )

        start = 0
        for ii in range( len( self.content ) ) :
            Nhalo = self.content[ ii ].Ncen + self.content[ ii ].Nsat
            content[ ii ] = self.content[ ii ].mask( mask[ start:start+Nhalo ] )
            start += Nhalo

        ww = numpy.where( [ obj is not None for obj in content ] )

        return catalogue( X = content[ ww ], boxsize = self.boxsize )
            

    def read_hierarchy_from_gadget ( self, filebase, boxsize = None ) :
        """ Reads the halo/sub-halo hierarchy from a Subgroup gadget output
        
        Parameters
        ----------
        filebase :
        
        Returns
        -------
        None        
        """

        self.boxsize = boxsize
        
        if self.gadget is None :
            self.gadget = gadget_file.gadget_file( filebase )
            for ii in range( self.gadget.glob[ 'task' ] ) :
                self.gadget.read_file( ii, scale_mass = self.scale_mass,
                                       scale_lenght = self.scale_lenght,
                                       add_to_internal = True )
        
        self.content = numpy.empty( self.gadget.glob[ 'tot_groups' ], 
                                    dtype = host_halo )
        sub_offset = 0 
        sub_var = [ 'coords', 'mass', 'spin', 'veldisp' ]
        tot_sat = 0
        tot_cen = 0
        for ii in range( self.gadget.glob[ 'tot_groups' ] ) :
            
            if self.gadget.nsubs[ ii ] > 0 :
                central = numpy.array( [ halo( coords  = self.gadget.sub_coord[ sub_offset ],
                                               mass    = self.gadget.sub_mass[ sub_offset ],
                                               spin    = self.gadget.sub_spin[ sub_offset ],
                                               veldisp = self.gadget.sub_veldisp[ sub_offset ] ) ] )
                
                satellites = numpy.array( [ halo( coords  = self.gadget.sub_coord[ jj ], 
                                                  mass    = self.gadget.sub_mass[ jj ], 
                                                  spin    = self.gadget.sub_spin[ jj ], 
                                                  veldisp = self.gadget.sub_veldisp[ jj ] ) 
                                            for jj in range( sub_offset + len( central ), 
                                                             sub_offset + self.gadget.nsubs[ ii ] ) ] )
            
            else :
                central = numpy.array([])
                satellites = numpy.array([])

            sub_offset += self.gadget.nsubs[ ii ]
            
            self.content[ ii ] = host_halo( coords     = self.gadget.coord[ ii ], 
                                            mass       = self.gadget.mass[ ii ], 
                                            central    = central,
                                            satellites = satellites )
            
    def get_coord_cen ( self, store = False ) :
        """ Get the coordinates of all the central objects in the catalogue
        
        Parameters
        ----------
        store : whether to store the returned array into an internal variable ( default = False )

        Returns
        -------
        Array of central objects coordinates ( shape = ( n_central, 3 ) )
        """
        coord_cen = numpy.array( [ cen.centre for obj in self.content for cen in obj.central ] )
        if store : 
            self.coord_cen = coord_cen
        return coord_cen
        
    def get_coord_sat ( self, store = False ) :
        """ Get the coordinates of all the satellite objects in the catalogue
        
        Parameters
        ----------
        store : whether to store the returned array into an internal variable ( default = False )

        Returns
        -------
        Array of satellite objects coordinates ( shape = ( n_satellites, 3 ) )
        """
        coord_sat = numpy.array( [ sat.centre for obj in self.content for sat in obj.satellites ] )
        if store : 
            self.coord_sat = coord_sat
        return coord_sat
        
    def get_mass_cen ( self, store = False ) :
        """ Get the masses of all the central objects in the catalogue
        
        Parameters
        ----------
        store : whether to store the returned array into an internal variable ( default = False )

        Returns
        -------
        Array of central objects mass ( shape = ( n_central, ) )
        """
        mass_cen = numpy.array( [ cen.mass for obj in self.content for cen in obj.central ] )
        if store : 
            self.mass_cen = mass_cen
        return mass_cen
        
    def get_mass_sat ( self, store = False ) :
        """ Get the masses of all the satellite objects in the catalogue
        
        Parameters
        ----------
        store : whether to store the returned array into an internal variable ( default = False )

        Returns
        -------
        Array of satellite objects mass ( shape = ( n_satellites, ) )
        """
        mass_sat = numpy.array( [ sat.mass for obj in self.content for sat in obj.satellites ] )
        if store : 
            self.mass_sat = mass_sat
        return mass_sat

    def get_mass_halo ( self, store = False ) :
        """ Get the masses of 
        
        Parameters
        ----------
        store : whether to store the returned array into an internal variable ( default = False )

        Returns
        -------
        Array of ( shape = ( n_central + n_satellites, ) )
        """
        mass_halo_cen = numpy.array( [ obj.mass for obj in self.content for _ in obj.central ] )
        mass_halo_sat = numpy.array( [ obj.mass for obj in self.content for _ in obj.satellites ] )
        mass_halo = numpy.append( mass_halo_cen, mass_halo_sat )
        if store : 
            self.mass_halo = mass_halo
        return mass_halo
    
    def populate ( self, model, extract = False, rank_order = False ) :
        
        ll = numpy.array( sorted( copy.deepcopy( self.content ),
                                  key = lambda x : x.mass,
                                  reverse = True ) )

        masses = numpy.array( [ obj.mass for obj in ll ] )
        action = numpy.array( [ [ model.Ncen( _m ), model.Nsat( _m ) ] for _m in masses ] )

        gxy_cat = catalogue( X = ll,
                             scale_mass = self.scale_mass,
                             scale_lenght = self.scale_lenght,
                             boxsize = self.boxsize )

        if not extract :
            return populate( gxy_cat, action, rank_order, probs = False )
        else :
            return select_galaxies( gxy_cat, action, rank_order, probs = False )


##########################################################################################
# catalogue manipulation functions

def _set_probs( cata, actions, rank_order ) :
    
    # centrals
    all_cen = [ obj.central for obj in cata.content ]
    Nc = numpy.random.binomial( numpy.ones( len( actions.T[0] ), dtype=int ), actions.T[ 0 ] )
    
    # satellites
    if rank_order :
        all_sat, Nsat = np.array( [ [ numpy.array( sorted( obj.satellites,
                                                           key = lambda x : x.mass,
                                                           reverse = True ) ),
                                      obj.Nsat ]
                                    for obj in cata.content ], dtype = object ).T
    else :      
        all_sat, Nsat = np.array( [ [ obj.satellites, obj.Nsat ] 
                                    for obj in catalogue.content ], dtype = object ).T  
        list( map( numpy.random.shuffle, all_sat ) )
    Ns = numpy.random.poisson( actions.T[1] * Nsat.astype( int ) )

    return all_cen, Nc, all_sat, Ns 

def _set_numbs( cata, actions, rank_order ) :
    
    # centrals
    all_cen = [ obj.central for obj in cata.content ]
    Nc = numpy.random.binomial( numpy.ones( len( actions.T[0] ), dtype=int ), actions.T[ 0 ] )
    
    # satellites
    if rank_order :
        all_sat = [ numpy.array( sorted( obj.satellites,
                                         key = lambda x : x.mass,
                                         reverse = True ) )
                    for obj in cata.content ]
    else :      
        all_sat = [ obj.satellites for obj in cata.content ]  
        list( map( numpy.random.shuffle, all_sat ) )
    Ns = numpy.random.poisson( actions.T[1] )

    return all_cen, Nc, all_sat, Ns 
    
def select_galaxies ( cata, actions, rank_order = False, probs = True ) :

    if probs :
        _set = _set_probs
    else :
        _set = _set_numbs
    
    return numpy.block( [ numpy.append( Ac[:nc], As[:ns] ) 
                          for Ac, nc, As, ns 
                          in zip( *_set( cata, actions, rank_order ) ) ] )

def _set_halo( obj, cens, sats ) :
    
    return host_halo( coords = obj.coords, 
                      mass = obj.mass,
                      central = cens, 
                      satellites = sats )

def populate ( cata, actions, rank_order = False, probs = True ) :

    if probs :
        _set = _set_probs
    else :
        _set = _set_numbs

    # content with same shape as halo catalogue
    content = numpy.array( [ _set_halo( obj, Ac[:nc], As[:ns] ) 
                             for obj, Ac, nc, As, ns 
                             in zip( cata.content, *_set( cata, actions, rank_order ) ) ] )
    
    return catalogue( X = content,
                      scale_mass = cata.scale_mass,
                      scale_lenght = cata.scale_lenght,
                      boxsize = cata.boxsize )

def get_abundances ( cata, bins ) :
    """
    """
        
    # these extract several properties stored in the internal container of catalogue
    Mh = numpy.array( [ obj.mass for obj in cata.content ] ) # mass of the halo
    Nc = numpy.array( [ obj.Ncen for obj in cata.content ] ) # number of centrals
    Ns = numpy.array( [ obj.Nsat for obj in cata.content ] ) # number of satellites
    
    # we obtain the index of the bin each mass belongs
    dig = numpy.digitize( Mh, bins )

    # mean number of centrals in each mass-bin
    Nc_binned = numpy.array( [ Nc[ numpy.where( dig == idx ) ].mean() 
                               if numpy.isin( idx, dig ) else 0. 
                               for idx in range( len( bins ) ) ] ) 
    
    # mean number of satellites in each mass-bin
    Ns_binned = numpy.array( [ Ns[ numpy.where( dig == idx ) ].mean() 
                               if numpy.isin( idx, dig ) else 0. 
                               for idx in range( len( bins ) ) ] )
    
    return Nc_binned, Ns_binned

def extract_galaxies ( hhaloes, ngxy ) :
    """
    Extracts objects of type galaxy from an host halo catalogue
    Parameters
    ----------
    hhaloes : the input host-halo catalogue
    ngxy : the number of galaxies hosted by all the host haloes
           in the input catalogue
    
    Return
    ------
    numpy array containing ngxy galaxies
    """
    
    galaxies = numpy.empty( ngxy, dtype = galaxy )
    idx = 0
    for obj in hhaloes :

        for cen in obj.central :
            galaxies[ idx ] = galaxy.from_halo( cen )
            idx += 1

        for sat in obj.satellites :
            galaxies[ idx ] = galaxy.from_halo( sat )
            idx += 1
            
    return galaxies


