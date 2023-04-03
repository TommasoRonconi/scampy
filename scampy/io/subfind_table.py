import numpy

class subfind_table () :
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
    global_elements =  ( 'tot_groups', 'tot_ids', 'Nfiles', 'tot_subs' )
    fields = ( 'GroupMass', 'GroupPos', 'GroupNsubs', 'GroupFirstSub',
               'GroupMassM200', 'GroupRadiusM200',
               'GroupMassC200', 'GroupRadiusC200',
               'GroupMassT200', 'GroupRadiusT200',
               'GroupVelDispM200', 'GroupVelDispC200', 'GroupVelDispT200', 
               'SubhaloMass', 'SubhaloPos', 'SubhaloSpin',
               'SubhaloVelDisp', 'SubhaloVmax', 'SubhaloVmaxRad',
               'SubhaloGrNr' )
    values = (
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,3), 'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'uint32'  },
        { 'shape' : (0,),  'dtype' : 'uint32'  },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,3), 'dtype' : 'float32' },
        { 'shape' : (0,3), 'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'float32' },
        { 'shape' : (0,),  'dtype' : 'uint32'  },
    )
    
    def __init__ ( self, filebase, byteorder = 'little', ids_lenght = 8, masstab = True ) :
        
        if not isinstance( filebase, str ) :
            raise TypeError( "filebase must be set to a string" )
        self.filebase = filebase
        
        self.byteorder = byteorder
        
        self.ids_lenght = ids_lenght
        
        self.masstab = masstab
        
        _, self.glob = self.read_header()
        self.loc = {}
        
        self.content = self.generate_content_dict()

        print(f"{self.glob['Nfiles']}\tfiles\n"
              f"{self.glob['tot_groups']}\thaloes\n"
              f"{self.glob['tot_subs']}\tsubhaloes\n")

    def generate_content_dict ( self ) :
        return dict( zip( self.fields, ( numpy.empty(**v) for v in self.values ) ) )
        
    def _read_header ( self, stream ) :
        
        Ngroups = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        totNgroups = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        Nids = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        totNids = int.from_bytes( stream.read( 8 ), byteorder = self.byteorder )
        Nfiles = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        Nsubs = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )
        totNsubs = int.from_bytes( stream.read( 4 ), byteorder = self.byteorder )

        return ( Ngroups, Nids, Nsubs) , ( totNgroups, totNids, Nfiles, totNsubs )
        
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
                    add_to_content = False,
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
        add_to_content : bool
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
            coord_block = scale_lenght * numpy.fromfile(
                f, dtype = 'float32',
                count = 3 * loc[ 'groups' ]
            ).reshape( loc[ 'groups' ], 3 )
            
            # Estimates blocks (mass, radius, veldisp)*(m200, c200, t200):
            # m200 : density > 200*Mean density of the Universe
            # c200 : density > 200*Critical density of the Universe
            # t200 : density in top-hat sphere with mean density Delta_c (from Bryan+1998)
            #        (the subscript '200' here can be ignored)
            mass_m200_block = scale_mass * numpy.fromfile( f, dtype = 'float32',
                                                           count = loc[ 'groups' ] )
            radius_m200_block = scale_lenght * numpy.fromfile( f, dtype = 'float32',
                                                               count = loc[ 'groups' ] )
            mass_c200_block = scale_mass * numpy.fromfile( f, dtype = 'float32',
                                                           count = loc[ 'groups' ] )
            radius_c200_block = scale_lenght * numpy.fromfile( f, dtype = 'float32',
                                                               count = loc[ 'groups' ] )
            mass_t200_block = scale_mass * numpy.fromfile( f, dtype = 'float32',
                                                           count = loc[ 'groups' ] )
            radius_t200_block = scale_lenght * numpy.fromfile( f, dtype = 'float32',
                                                               count = loc[ 'groups' ] )
            veldisp_m200_block = numpy.fromfile( f, dtype = 'float32',
                                                 count = loc[ 'groups' ] )
            veldisp_c200_block = numpy.fromfile( f, dtype = 'float32',
                                                 count = loc[ 'groups' ] )
            veldisp_t200_block = numpy.fromfile( f, dtype = 'float32',
                                                 count = loc[ 'groups' ] )

            # jump 2 empty blocks
            f.seek( 2 * 4 * loc[ 'groups' ], 1 )
            
            # Sub-Halos number:
            nsubs_block = numpy.fromfile( f, dtype = 'uint32',
                                          count = loc[ 'groups' ] )
            
            # jump last groups block
            first_sub_block = numpy.fromfile( f, dtype = 'uint32',
                                              count = loc[ 'groups' ] )
            
            ##############################
            ######### SUB-GROUPS #########
            ##############################
            
            # jump lenght, offset, parent block
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # Sub-Halo mass:
            sub_mass_block = scale_mass * numpy.fromfile( f, dtype = 'float32',
                                                          count = loc[ 'subs' ] )
            
            # Sub-Halo coords:
            sub_coord_block = scale_lenght * numpy.fromfile(
                f, dtype = 'float32',
                count = 3 * loc[ 'subs' ]
            ).reshape( loc[ 'subs' ], 3 )
            
            # jump vel block:
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # jump CM block:
            f.seek( 3 * 4 * loc[ 'subs' ], 1 )
            
            # Sub-Halo spin:
            sub_spin_block = numpy.fromfile(
                f, dtype = 'float32',
                count = 3 * loc[ 'subs' ]
            ).reshape( loc[ 'subs' ], 3 )
            
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
    
            if add_to_content :
                self.content[ 'GroupMass' ] = numpy.append(
                    self.content[ 'GroupMass' ], mass_block )
                self.content[ 'GroupPos' ] = numpy.append(
                    self.content[ 'GroupPos' ], coord_block, axis = 0 ) 
                self.content[ 'GroupNsubs' ] = numpy.append(
                    self.content[ 'GroupNsubs' ], nsubs_block )
                self.content[ 'GroupFirstSub' ] = numpy.append(
                    self.content[ 'GroupFirstSub' ], first_sub_block )
                self.content[ 'GroupMassM200' ] = numpy.append(
                    self.content[ 'GroupMassM200' ], mass_m200_block )
                self.content[ 'GroupRadiusM200' ] = numpy.append(
                    self.content[ 'GroupRadiusM200' ], radius_m200_block )
                self.content[ 'GroupMassC200' ] = numpy.append(
                    self.content[ 'GroupMassC200' ], mass_c200_block )
                self.content[ 'GroupRadiusC200' ] = numpy.append(
                    self.content[ 'GroupRadiusC200' ], radius_c200_block )
                self.content[ 'GroupMassT200' ] = numpy.append(
                    self.content[ 'GroupMassT200' ], mass_t200_block )
                self.content[ 'GroupRadiusT200' ] = numpy.append(
                    self.content[ 'GroupRadiusT200' ],  radius_t200_block )
                self.content[ 'GroupVelDispM200' ] = numpy.append(
                    self.content[ 'GroupVelDispM200' ], veldisp_m200_block )
                self.content[ 'GroupVelDispC200' ] = numpy.append(
                    self.content[ 'GroupVelDispC200' ], veldisp_c200_block )
                self.content[ 'GroupVelDispT200' ] = numpy.append(
                    self.content[ 'GroupVelDispT200' ], veldisp_t200_block )
                self.content[ 'SubhaloMass' ] = numpy.append(
                    self.content[ 'SubhaloMass' ], sub_mass_block )
                self.content[ 'SubhaloPos' ] = numpy.append(
                    self.content[ 'SubhaloPos' ], sub_coord_block, axis = 0 )
                self.content[ 'SubhaloSpin' ] = numpy.append(
                    self.content[ 'SubhaloSpin' ], sub_spin_block, axis = 0 )
                self.content[ 'SubhaloVelDisp' ] = numpy.append(
                    self.content[ 'SubhaloVelDisp' ], sub_veldisp_block )
                self.content[ 'SubhaloVmax' ] = numpy.append(
                    self.content[ 'SubhaloVmax' ], sub_vmax_block )
                self.content[ 'SubhaloVmaxRad' ] = numpy.append(
                    self.content[ 'SubhaloVmaxRad' ], sub_vmaxrad_block )
                self.content[ 'SubhaloGrNr' ] = numpy.append(
                    self.content[ 'SubhaloGrNr' ], sub_group_block )
                
            # End of file context manager
        return;
        
    def read_all_files ( self, reset = True, **kwargs ) :

        self.content = self.generate_content_dict()
        kwargs.update( add_to_content = True )
        for i in range(self.glob['Nfiles']) :
            self.read_file( i, **kwargs )
        return;

    # @classmethod
    # def load ( cls, infile ) :
    #     """ classmethod for loading from hdf5 file
        
    #     Parameters
    #     ----------
    #     infile : string
    #       the '/path/to/name/of/file.hdf5' to load
        
    #     Returns
    #     -------
    #     : catalogue
    #     """
    #     import os, h5py

    #     if not os.path.isfile( infile ) :
    #         raise IOError( f"file {infile} does not exist." )
        
    #     with h5py.File( infile, "r" ) as f :

    #         # load meta-data
    #         metadata = dict(f.attrs)
            
    #         # load halo catalogue
    #         haloes = haloCat(len(f['haloes']['mask']), 
    #                                 f['haloes'] )
    #         haloes.add_mask( f['haloes']['mask'] )
            
    #         # load subhalo catalogue
    #         subhaloes = subhaloCat( len(f['subhaloes']['mask']), f['subhaloes'] )
    #         subhaloes.add_mask( f['subhaloes']['mask'] )
            
    #     return cls( haloes, subhaloes, **metadata )
