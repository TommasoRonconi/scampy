# external:
import numpy
import copy
from sklearn.neighbors import BallTree, KDTree

# internal:
from scampy import gadget_file
from scampy.objects import *

class catalogue () :
    
    def __init__ ( self, X = None, scale_lenght = 1.e-3, scale_mass = 1.e+10 ) :
        if X is None : 
            self.content = numpy.array([])
        else :
            self.content = X 
        self.scale_lenght = scale_lenght
        self.scale_mass = scale_mass
            
        self.gadget = None
        self.tree = None
    
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
    
    def read_balltree_from_gadget ( self, filebase ) :
        """ Creates instance of sklearn.BallTree from a Subgroup gadget output
        
        Parameters
        ----------
        filebase :
        
        Returns
        -------
        None        
        """
        
        if self.gadget is None :
            self.gadget = gadget_file.gadget_file( filebase )
            for ii in range( self.gadget.glob[ 'task' ] ) :
                self.gadget.read_file( ii, scale_mass = self.scale_mass,
                                       scale_lenght = self.scale_lenght,
                                       add_to_internal = True )
        
        self.tree = BallTree( self.gadget.sub_coord, leaf_size = 10 )
    
    def read_kdtree_from_gadget ( self, filebase ) :
        """ Creates instance of sklearn.BallTree from a Subgroup gadget output
        
        Parameters
        ----------
        filebase :
        
        Returns
        -------
        None        
        """
        
        if self.gadget is None :
            self.gadget = gadget_file.gadget_file( filebase )
            for ii in range( self.gadget.glob[ 'task' ] ) :
                self.gadget.read_file( ii, scale_mass = self.scale_mass,
                                       scale_lenght = self.scale_lenght,
                                       add_to_internal = True )
        
        self.tree = KDTree( self.gadget.sub_coord, leaf_size = 10 )        

    def read_hierarchy_from_gadget ( self, filebase ) :
        """ Reads the halo/sub-halo hierarchy from a Subgroup gadget output
        
        Parameters
        ----------
        filebase :
        
        Returns
        -------
        None        
        """
        
        
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

    # def get_num_cen ( self, store = False ) :

    #     num_cen += len( obj.central ) for obj in self.content
    #     if store :
    #         self.num_cen = num_cen

    #     return num_cen
            
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
        Array of central objects mass ( shape = ( n_central, 3 ) )
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
        Array of satellite objects mass ( shape = ( n_satellites, 3 ) )
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
        Array of ( shape = ( , 3 ) )
        """
        mass_halo_cen = numpy.array( [ obj.mass for obj in self.content for _ in obj.central ] )
        mass_halo_sat = numpy.array( [ obj.mass for obj in self.content for _ in obj.satellites ] )
        mass_halo = numpy.append( mass_halo_cen, mass_halo_sat )
        if store : 
            self.mass_halo = mass_halo
        return mass_halo
    
    def populate ( self, model, extract = False ) :
        
        ll = numpy.array( sorted( copy.deepcopy( self.content ),
                                  key = lambda x : x.mass,
                                  reverse = True ) )

        Ngxy = 0
        for obj in ll :
            
            Nc = model.Ncen( obj.mass )
            select = numpy.random.choice( [ True, False ],
                                          size = 1,
                                          replace = False, 
                                          p = [ Nc, 1 - Nc ] )
            if select :
                central = obj.central
            else :
                central = numpy.array([])
            obj.set_central( central = central )
            
            Ns = numpy.random.poisson( model.Nsat( obj.mass ) )
            if Ns < len( obj.satellites ) :
                satellites = numpy.random.choice( obj.satellites, 
                                                  size = int( round( Ns ) ), 
                                                  replace = False )
            else :
                satellites = obj.satellites
            obj.set_satellites( satellites = satellites )

            Ngxy += ( obj.Ncen + obj.Nsat )

        if extract :
            return extract_galaxies( ll, Ngxy )
        else :
            return catalogue( ll ), Ngxy

def extract_galaxies ( hhaloes, ngxy ) :
    
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

if __name__ == '__main__' :

    cat = catalogue()
    base = "/home/tomi/phd/gadget_catalogues/lcdm_mis0.125/output/groups_011/subhalo_tab_011"
    cat.read_hierarchy_from_gadget( base )

    coord_cen = cat.get_coord_cen()
    coord_sat = cat.get_coord_sat()

    print( coord_cen.shape, coord_sat.shape )

    mass_cen = cat.get_mass_cen()
    mass_sat = cat.get_mass_sat()

    print( mass_cen.shape, mass_sat.shape )


