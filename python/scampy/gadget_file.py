import numpy

class gadget_file () :
    """
    A class for reading GaDGET SUBFIND subgroup tables
    
    Parameters
    ----------
    filebase : string
      common name of the files to be read
    byteorder : string
      'little' or 'big' for little or big endian, respectively
    ids_lenght : int
      number of bytes of the integer storing the IDs (8 or 16)
    masstab : bool
      whether the file contains mass tables
    """
    
    local_elements = ( 'groups', 'ids', 'subs' )
    global_elements =  ( 'tot_groups', 'tot_ids', 'task', 'tot_subs' )
    mass        = numpy.empty( (0,),  dtype = 'float32' )
    coord       = numpy.empty( (0,3), dtype = 'float32' )
    nsubs       = numpy.empty( (0,),  dtype = 'uint32'  )
    sub_mass    = numpy.empty( (0,),  dtype = 'float32' )
    sub_coord   = numpy.empty( (0,3), dtype = 'float32' )
    sub_spin    = numpy.empty( (0,3), dtype = 'float32' )
    sub_veldisp = numpy.empty( (0,),  dtype = 'float32' )
    sub_vmax    = numpy.empty( (0,),  dtype = 'float32' )
    sub_vmaxrad = numpy.empty( (0,),  dtype = 'float32' )
    sub_group   = numpy.empty( (0,),  dtype = 'uint32'  )
    
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
    
    def read_file ( self, num,
                    scale_mass = 1.e+10,
                    scale_lenght = 1.e-3,
                    add_to_internal = False,
                    verbose = False ) :
        """ Reads the `num`^th file in `base`
        
        Parameters
        ----------
        num : int
          file to read
        scale_mass : float
          mass unit (in terms of solar masses, default = 1.e+10)
        scale_lenght : float
          lenght unit (in terms of Mpc/h, default = 1.e-3)
        add_to_internal : bool
          whether to add the data to the internal storage array of the class (default = False)
        verbose : bool
          whether to print on screen additional info (default = False)
        
        Returns
        -------
        None
        """
        
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
            mass_block = scale_mass * numpy.fromfile( f, dtype = 'float32',
                                                      count = loc[ 'groups' ] )
            
            # Halo coords:
            coord_block = scale_lenght * numpy.fromfile( f, dtype = 'float32',
                                                         count = 3 * loc[ 'groups' ] ).reshape( loc[ 'groups' ], 3 )
            
            # jump estimates block (mass, radius, veldisp)*(m200, c200, t200)
            f.seek( 9 * 4 * loc[ 'groups' ], 1 ) 
            
            # jump 2 blocks I don't know what they are
            f.seek( 2 * 4 * loc[ 'groups' ], 1 )
            
            # Sub-Halos number:
            nsubs_block = numpy.fromfile( f, dtype = 'uint32',
                                          count = loc[ 'groups' ] )
            
            # jump last groups block
            f.seek( 4 * loc[ 'groups' ], 1 )
            
            ##############################
            ######### SUB-GROUPS #########
            ##############################
            
            # jump lenght, offset, parent block
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # Sub-Halo mass:
            sub_mass_block = scale_mass * numpy.fromfile( f, dtype = 'float32',
                                                          count = loc[ 'subs' ] )
            
            # Sub-Halo coords:
            sub_coord_block = scale_lenght * numpy.fromfile( f, dtype = 'float32',
                                                             count = 3 * loc[ 'subs' ] ).reshape( loc[ 'subs' ], 3 )
            
            # jump vel block:
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # jump CM block:
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # Sub-Halo spin:
            sub_spin_block = numpy.fromfile( f, dtype = 'float32',
                                             count = 3 * loc[ 'subs' ] ).reshape( loc[ 'subs' ], 3 )
            
            # Sub-Halo veldisp:
            sub_veldisp_block = numpy.fromfile( f, dtype = 'float32',
                                                count = loc[ 'subs' ] )
            
            # Sub-Halo Vmax:
            sub_vmax_block = numpy.fromfile( f, dtype = 'float32',
                                             count = loc[ 'subs' ] )
            
            # Sub-Halo Vmax-rad:
            sub_vmaxrad_block = numpy.fromfile( f, dtype = 'float32',
                                                count = loc[ 'subs' ] )
            
            # jump half-mass radius block:
            f.seek( 4 * loc[ 'subs' ], 1 )
            
            # jump ID most bound block:
            f.seek( self.ids_lenght * loc[ 'subs' ], 1 )
            
            # Group number:
            sub_group_block = numpy.fromfile( f, dtype = 'uint32', count = loc[ 'subs' ] )
            
            # jump mass tab (if present):
            if self.masstab :
                f.seek( 6 * 4 * loc[ 'subs' ], 1 )
            
            # check end-of-file is reached:
            left = f.read()
            if len( left ) > 0 and verbose :
                print( "Did not reach EOF:\t {:.8f} bytes left.".format( len( left ) ) )
                
            if add_to_internal :
                self.mass = numpy.append( self.mass, mass_block )
                self.coord = numpy.append( self.coord, coord_block, axis = 0 )
                self.nsubs = numpy.append( self.nsubs, nsubs_block )
                self.sub_mass = numpy.append( self.sub_mass, sub_mass_block )
                self.sub_coord = numpy.append( self.sub_coord, sub_coord_block, axis = 0 )
                self.sub_spin = numpy.append( self.sub_spin, sub_spin_block, axis = 0 )
                self.sub_veldisp = numpy.append( self.sub_veldisp, sub_veldisp_block )
                self.sub_vmax = numpy.append( self.sub_vmax, sub_vmax_block )
                self.sub_vmaxrad = numpy.append( self.sub_vmaxrad, sub_vmaxrad_block )
                self.sub_group = numpy.append( self.sub_group, sub_group_block )

            return;

