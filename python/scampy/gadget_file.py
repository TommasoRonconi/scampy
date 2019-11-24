import numpy as np

class gadget_file () :
    """
    A class for reading GaDGET SUBFIND subgroup tables
    """
    
    local_elements = ( 'groups', 'ids', 'subs' )
    global_elements =  ( 'tot_groups', 'tot_ids', 'task', 'tot_subs' )
    mass        = np.empty( (0,),  dtype = 'float32' )
    coord       = np.empty( (0,3), dtype = 'float32' )
    nsubs       = np.empty( (0,),  dtype = 'uint32'  )
    sub_mass    = np.empty( (0,),  dtype = 'float32' )
    sub_coord   = np.empty( (0,3), dtype = 'float32' )
    sub_spin    = np.empty( (0,3), dtype = 'float32' )
    sub_veldisp = np.empty( (0,),  dtype = 'float32' )
    sub_vmax    = np.empty( (0,),  dtype = 'float32' )
    sub_vmaxrad = np.empty( (0,),  dtype = 'float32' )
    sub_group   = np.empty( (0,),  dtype = 'uint32'  )
    
    def __init__ ( self, filebase, byteorder = 'little', ids_lenght = 8, masstab = True ) :
        
        if not isinstance( filebase, str ) :
            raise TypeError( "filebase must be set to a string" )
        self.filebase = filebase
        
        self.byteorder = byteorder
        
        self.ids_lenght = ids_lenght
        
        self.masstab = masstab
        
        _, self.glob = self.read_header()
        self.loc = {}
        
    def _read_header ( self, stream ) :
        
        Ngroups = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        totNgroups = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        Nids = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        totNids = int.from_bytes( stream.read( 8 ), byteorder = self.byteorder )
        Ntask = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        Nsubs = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        totNsubs = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )

        return ( Ngroups, Nids, Nsubs) , ( totNgroups, totNids, Ntask, totNsubs )
        
    def read_header ( self, num = 0 ) :
        """
        Reads only the Header of the `num`^th file in `base`
        Parameters
        ----------
        num : the file number to read
        
        Returns
        -------
        loc : dictionary with meta-data local to file `num`
        glob : dictionary with meta-data global to all files in `base`
        """
        
        current = self.filebase + ".{:d}".format( num )
        with open( current, "rb" ) as f:
            _local, _global = self._read_header( f )
            loc = dict( ( name, elem ) for name, elem in zip( self.local_elements, _local ) )
            glob = dict( ( name, elem ) for name, elem in zip( self.global_elements, _global ) )
            
        return loc, glob
    
    def read_file ( self, num, scale_mass = 1.e+10, scale_lenght = 1.e-3, add_to_internal = False, verbose = False ) :
        
        current = self.filebase + ".{:d}".format( num )
        with open( current, "rb" ) as f :
            
            # Header:
            _local, _ = self._read_header( f )
            loc = dict( ( name, elem ) for name, elem in zip( self.local_elements, _local ) )
            
            ##############################
            ########### GROUPS ###########
            ##############################
            
            # jump lenght and offset blocks:
            f.seek( 2 * 4 * loc[ 'groups' ], 1 ) 
            
            # Halo mass:
            mass_block = scale_mass * np.fromfile( f, dtype = 'float32', count = loc[ 'groups' ] )
            
            # Halo coords:
            coord_block = scale_lenght * np.fromfile( f, dtype = 'float32', count = 3 * loc[ 'groups' ] ).reshape( loc[ 'groups' ], 3 )
            
            # jump estimates block (mass, radius, veldisp)*(m200, c200, t200)
            f.seek( 9 * 4 * loc[ 'groups' ], 1 ) 
            
            # jump 2 blocks I don't know what they are
            f.seek( 2 * 4 * loc[ 'groups' ], 1 )
            
            # Sub-Halos number:
            nsubs_block = np.fromfile( f, dtype = 'uint32', count = loc[ 'groups' ] )
            
            # jump last groups block
            f.seek( 4 * loc[ 'groups' ], 1 )
            
            ##############################
            ######### SUB-GROUPS #########
            ##############################
            
            # jump lenght, offset, parent block
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # Sub-Halo mass:
            sub_mass_block = scale_mass * np.fromfile( f, dtype = 'float32', count = loc[ 'subs' ] )
            
            # Sub-Halo coords:
            sub_coord_block = scale_lenght * np.fromfile( f, dtype = 'float32', count = 3 * loc[ 'subs' ] ).reshape( loc[ 'subs' ], 3 )
            
            # jump vel block:
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # jump CM block:
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # Sub-Halo spin:
            sub_spin_block = np.fromfile( f, dtype = 'float32', count = 3 * loc[ 'subs' ] ).reshape( loc[ 'subs' ], 3 )
            
            # Sub-Halo veldisp:
            sub_veldisp_block = np.fromfile( f, dtype = 'float32', count = loc[ 'subs' ] )
            
            # Sub-Halo Vmax:
            sub_vmax_block = np.fromfile( f, dtype = 'float32', count = loc[ 'subs' ] )
            
            # Sub-Halo Vmax-rad:
            sub_vmaxrad_block = np.fromfile( f, dtype = 'float32', count = loc[ 'subs' ] )
            
            # jump half-mass radius block:
            f.seek( 4 * loc[ 'subs' ], 1 )
            
            # jump ID most bound block:
            f.seek( self.ids_lenght * loc[ 'subs' ], 1 )
            
            # Group number:
            sub_group_block = np.fromfile( f, dtype = 'uint32', count = loc[ 'subs' ] )
            
            # jump mass tab (if present):
            if self.masstab :
                f.seek( 6 * 4 * loc[ 'subs' ], 1 )
            
            # check end-of-file is reached:
            left = f.read()
            if len( left ) > 0 and verbose :
                print( "Did not reach EOF:\t {:.8f} bytes left.".format( len( left ) ) )
                
            if add_to_internal :
                self.mass = np.append( self.mass, mass_block )
                self.coord = np.append( self.coord, coord_block, axis = 0 )
                self.nsubs = np.append( self.nsubs, nsubs_block )
                self.sub_mass = np.append( self.sub_mass, sub_mass_block )
                self.sub_coord = np.append( self.sub_coord, sub_coord_block, axis = 0 )
                self.sub_spin = np.append( self.sub_spin, sub_spin_block, axis = 0 )
                self.sub_veldisp = np.append( self.sub_veldisp, sub_veldisp_block )
                self.sub_vmax = np.append( self.sub_vmax, sub_vmax_block )
                self.sub_vmaxrad = np.append( self.sub_vmaxrad, sub_vmaxrad_block )
                self.sub_group = np.append( self.sub_group, sub_group_block )

if __name__ == '__main__' :

    gadg = gadget_file( "/home/tomi/phd/gadget_catalogues/lcdm_mis0.125/output/groups_011/subhalo_tab_011" )
    gadg.read_header()
    
    for ii in range( gadg.glob['task'] ) :
        gadg.read_file( ii, add_to_internal = True, verbose = True )

    print( "{0:d} = {0:d} ?".format( gadg.glob[ 'tot_groups' ], gadg.mass.shape[ 0 ] ) )

    print( "{0:d} = {0:d} ?".format( gadg.glob[ 'tot_subs' ], gadg.sub_mass.shape[ 0 ] ) )
