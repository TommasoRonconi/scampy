from .cwrap.cwrap import *
from scampy.occupation_p import *

# Wrap class halo_model:
class halo_model () :
    """
    Class to handle the halo-model functions, 
    it depends on the classes occupation_p and cosmology
    
    Parameters
    ----------
    occupation : occupation_p
      object that defines the occupation probability functions
    cosmology : cosmology
      defines methods and functions dependent on the chosen cosmology
    redshift : float
      the redshift at which to compute the halo-model functions (should be greater than zero)
    thinness : int
      refinement of the interpolation grid (should be greater than 10)
    """

    def __init__ ( self,
                   occupation = None,
                   cosmology = None,
                   redshift = 1.e-7,
                   thinness = 50 ) :

        self.cosmology = cosmology
        self.thinness = c_size_t( thinness )
        self.redshift = c_double( redshift )
        self.handler = occupation
        if isinstance( self.handler, harikane16_p ) :
            self.obj = lib_hm.create_halo_model_H16( self.handler.obj,
                                                  self.cosmology.obj,
                                                  self.redshift,
                                                  self.thinness )
            self.set_parameters = self._set_parameters_H16
            
        elif isinstance( self.handler, tinker10_p ) :
            self.obj = lib_hm.create_halo_model_T10( self.handler.obj,
                                                  self.cosmology.obj,
                                                  self.redshift,
                                                  self.thinness )
            self.set_parameters = self._set_parameters_T10
            
        else :
            raise ValueError( 'Wrong value of occupation argument' ) 
        
    def __del__ ( self ) :
        """ Calls the destructor of the halo_model object
        """
        
        # Python call to halo_model dtor:
        lib_hm.free_halo_model( self.obj )

    def _set_parameters_H16 ( self,
                              DC = 0.5,
                              M_min = 1.e+11,
                              sigma_logM = 1.,
                              M0 = 1.e+10,
                              M1 = 1.e+10,
                              alpha = 1. ) :

        lib_hm.set_parameters_hm_H16( DC, M_min, sigma_logM, M0, M1, alpha, self.obj )

        return;

    def _set_parameters_T10 ( self,
                              Amin = 1.e+11,
                              siglogA = 1.,
                              Asat = 1.e+10,
                              alpsat = 1. ) :

        lib_hm.set_parameters_hm_T10( Amin, siglogA, Asat, alpsat, self.obj )

        return;

    def ng ( self ) :
        """ Average number of galaxies per unit volume.
        This is the halo-model estimate of the number density of galaxies at the redshift the model is computed:
          .. math:: n_g(z) = \int_{M_{min}}^{M_{max}} N_g(M_h) n(M_h) dM_h
        where :math:`N_g(M_h) = N_{cen}(M_h) + N_{sat}(M_h)`.

        Returns
        -------
        float
          :math:`n_g(z)`

        Note
        ----
        The values of :math:`z`, :math:`M_{min}` and :math:`M_{max}` are set by the constructor of the class halo_model
        """

        return lib_hm.ng_hm( self.obj )

    def bias ( self ) :

        return lib_hm.bias_hm( self.obj )

    def Mg ( self ) :

        return lib_hm.Mhalo_hm( self.obj )

    def dngdM ( self, Mh ) :

        return lib_hm.dngdM_hm( Mh, self.obj )

    def Pk ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk = ( c_double * self.thinness )()
        lib_hm.model_Pk_hm( kv, Pk, self.obj )
        
        return list( kv ), list( Pk )

    def Pk_1halo ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk_1h = ( c_double * self.thinness )()
        lib_hm.model_Pk_1halo_hm( kv, Pk_1h, self.obj )
        
        return list( kv ), list( Pk_1h )

    def Pk_cs ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk_cs = ( c_double * self.thinness )()
        lib_hm.model_Pk_cs_hm( kv, Pk_cs, self.obj )
        
        return list( kv ), list( Pk_cs )

    def Pk_ss ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk_ss = ( c_double * self.thinness )()
        lib_hm.model_Pk_ss_hm( kv, Pk_ss, self.obj )
        
        return list( kv ), list( Pk_ss )

    def Pk_2halo ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk_2h = ( c_double * self.thinness )()
        lib_hm.model_Pk_2halo_hm( kv, Pk_2h, self.obj )
        
        return list( kv ), list( Pk_2h )

    def Xi ( self, rv ) :

        import numpy as np
        
        rv = ( c_double * len( rv ) )( *[ np.float64( _r ) for _r in rv ] )
        Xi = ( c_double * len( rv ) )()
        size = c_uint( len( rv ) )
        lib_hm.model_Xi_hm( rv, Xi, size, self.obj )
        
        return list( Xi )

    def Xi_1halo ( self, rv ) :

        import numpy as np
        
        rv = ( c_double * len( rv ) )( *[ np.float64( _r ) for _r in rv ] )
        Xi_1h = ( c_double * len( rv ) )()
        size = c_uint( len( rv ) )
        lib_hm.model_Xi_1halo_hm( rv, Xi_1h, size, self.obj )
        
        return list( Xi_1h )

    def Xi_2halo ( self, rv ) :

        import numpy as np
        
        rv = ( c_double * len( rv ) )( *[ np.float64( _r ) for _r in rv ] )
        Xi_2h = ( c_double * len( rv ) )()
        size = c_uint( len( rv ) )
        lib_hm.model_Xi_2halo_hm( rv, Xi_2h, size, self.obj )
        
        return list( Xi_2h )

    def Wr ( self, rp ) :

        import numpy as np
        
        rp = ( c_double * len( rp ) )( *[ np.float64( _r ) for _r in rp ] )
        Wr = ( c_double * len( rp ) )()
        size = c_uint( len( rp ) )
        lib_hm.model_Wr_hm( rp, Wr, size, self.obj )
        
        return list( Wr )

    def Wr_1halo ( self, rp ) :

        import numpy as np
        
        rp = ( c_double * len( rp ) )( *[ np.float64( _r ) for _r in rp ] )
        Wr_1h = ( c_double * len( rp ) )()
        size = c_uint( len( rp ) )
        lib_hm.model_Wr_1halo_hm( rp, Wr_1h, size, self.obj )
        
        return list( Wr_1h )

    def Wr_2halo ( self, rp ) :

        import numpy as np
        
        rp = ( c_double * len( rp ) )( *[ np.float64( _r ) for _r in rp ] )
        Wr_2h = ( c_double * len( rp ) )()
        size = c_uint( len( rp ) )
        lib_hm.model_Wr_2halo_hm( rp, Wr_2h, size, self.obj )
        
        return list( Wr_2h )

    def Wt ( self, tt ) :

        import numpy as np
        
        tt = ( c_double * len( tt ) )( *[ np.float64( _t ) for _t in tt ] )
        Wt = ( c_double * len( tt ) )()
        size = c_uint( len( tt ) )
        lib_hm.model_Wt_hm( tt, Wt, size, self.obj )
        
        return list( Wt )

    def Wt_1halo ( self, tt ) :

        import numpy as np
        
        tt = ( c_double * len( tt ) )( *[ np.float64( _t ) for _t in tt ] )
        Wt_1h = ( c_double * len( tt ) )()
        size = c_uint( len( tt ) )
        lib_hm.model_Wt_1halo_hm( tt, Wt_1h, size, self.obj )
        
        return list( Wt_1h )

    def Wt_cs ( self, tt ) :

        import numpy as np
        
        tt = ( c_double * len( tt ) )( *[ np.float64( _t ) for _t in tt ] )
        Wt_cs = ( c_double * len( tt ) )()
        size = c_uint( len( tt ) )
        lib_hm.model_Wt_cs_hm( tt, Wt_cs, size, self.obj )
        
        return list( Wt_cs )

    def Wt_ss ( self, tt ) :

        import numpy as np
        
        tt = ( c_double * len( tt ) )( *[ np.float64( _t ) for _t in tt ] )
        Wt_ss = ( c_double * len( tt ) )()
        size = c_uint( len( tt ) )
        lib_hm.model_Wt_ss_hm( tt, Wt_ss, size, self.obj )
        
        return list( Wt_ss )

    def Wt_2halo ( self, tt ) :

        import numpy as np
        
        tt = ( c_double * len( tt ) )( *[ np.float64( _t ) for _t in tt ] )
        Wt_2h = ( c_double * len( tt ) )()
        size = c_uint( len( tt ) )
        lib_hm.model_Wt_2halo_hm( tt, Wt_2h, size, self.obj )
        
        return list( Wt_2h )

    def Wt_large_scale ( self, tt ) :

        import numpy as np
        
        tt = ( c_double * len( tt ) )( *[ np.float64( _t ) for _t in tt ] )
        Wt = ( c_double * len( tt ) )()
        size = c_uint( len( tt ) )
        lib_hm.model_Wt_large_scale_hm( tt, Wt, size, self.obj )
        
        return list( Wt )

#################################################################################
#################################################################################
#################################################################################
#################################################################################

# Wrap class cross_halo_model:
class cross_halo_model () :

    def __init__ ( self,
                   occupation1 = None,
                   occupation2 = None,
                   cosmology = None,
                   redshift = 1.e-7,
                   thinness = 50 ) :

        self.cosmology = cosmology
        self.thinness = c_size_t( thinness )
        self.redshift = c_double( redshift )
        self.ocp1 = occupation1
        self.ocp2 = occupation2
        if isinstance( self.ocp1, harikane16_p ) and isinstance( self.ocp2, harikane16_p ) :
            self.obj = lib_hm.create_cross_halo_model_H16( self.ocp1.obj,
                                                        self.ocp2.obj,
                                                        self.cosmology.obj,
                                                        self.redshift,
                                                        self.thinness )
            self.set_parameters_pop1 = self._set_parameters_pop1_H16
            self.set_parameters_pop2 = self._set_parameters_pop2_H16
            
        elif isinstance( self.ocp1, tinker10_p ) and isinstance( self.ocp2, tinker10_p ) :
            self.obj = lib_hm.create_cross_halo_model_T10( self.ocp1.obj,
                                                        self.ocp2.obj,
                                                        self.cosmology.obj,
                                                        self.redshift,
                                                        self.thinness )
            self.set_parameters_pop1 = self._set_parameters_pop1_T10
            self.set_parameters_pop2 = self._set_parameters_pop2_T10
            
        else :
            raise ValueError( 'Wrong value of occupation arguments' ) 
        
    def __del__ ( self ) :
        
        # Python call to cross_halo_model dtor:
        lib_hm.free_cross_halo_model( self.obj )

    def _set_parameters_pop1_H16 ( self,
                                   DC = 0.5,
                                   M_min = 1.e+11,
                                   sigma_logM = 1.,
                                   M0 = 1.e+10,
                                   M1 = 1.e+10,
                                   alpha = 1. ) :

        lib_hm.set_parameters_pop1_chm_H16( DC, M_min, sigma_logM, M0, M1, alpha, self.obj )

        return;

    def _set_parameters_pop2_H16 ( self,
                                   DC = 0.5,
                                   M_min = 1.e+11,
                                   sigma_logM = 1.,
                                   M0 = 1.e+10,
                                   M1 = 1.e+10,
                                   alpha = 1. ) :

        lib_hm.set_parameters_pop2_chm_H16( DC, M_min, sigma_logM, M0, M1, alpha, self.obj )

        return;

    def _set_parameters_pop1_T10 ( self,
                                   Amin = 1.e+11,
                                   siglogA = 1.,
                                   Asat = 1.e+10,
                                   alpsat = 1. ) :

        lib_hm.set_parameters_pop1_chm_T10( Amin, siglogA, Asat, alpsat, self.obj )

        return;

    def _set_parameters_pop2_T10 ( self,
                                   Amin = 1.e+11,
                                   siglogA = 1.,
                                   Asat = 1.e+10,
                                   alpsat = 1. ) :

        lib_hm.set_parameters_pop2_chm_T10( Amin, siglogA, Asat, alpsat, self.obj )

        return;

    def ng1 ( self ) :

        return lib_hm.ng1_chm( self.obj )

    def ng2 ( self ) :

        return lib_hm.ng2_chm( self.obj )

    def Pk ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk = ( c_double * self.thinness )()
        lib_hm.model_Pk_chm( kv, Pk, self.obj )
        
        return list( kv ), list( Pk )

    def Pk_1halo ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk_1h = ( c_double * self.thinness )()
        lib_hm.model_Pk_1halo_chm( kv, Pk_1h, self.obj )
        
        return list( kv ), list( Pk_1h )

    def Pk_2halo ( self ) :
        
        kv = ( c_double * self.thinness )()
        Pk_2h = ( c_double * self.thinness )()
        lib_hm.model_Pk_2halo_chm( kv, Pk_2h, self.obj )
        
        return list( kv ), list( Pk_2h )

    def Xi ( self, rv ) :

        import numpy as np
        
        rv = ( c_double * len( rv ) )( *[ np.float64( _r ) for _r in rv ] )
        Xi = ( c_double * len( rv ) )()
        size = c_uint( len( rv ) )
        lib_hm.model_Xi_chm( rv, Xi, size, self.obj )
        
        return list( Xi )

    def Xi_1halo ( self, rv ) :

        import numpy as np
        
        rv = ( c_double * len( rv ) )( *[ np.float64( _r ) for _r in rv ] )
        Xi_1h = ( c_double * len( rv ) )()
        size = c_uint( len( rv ) )
        lib_hm.model_Xi_1halo_chm( rv, Xi_1h, size, self.obj )
        
        return list( Xi_1h )

    def Xi_2halo ( self, rv ) :

        import numpy as np
        
        rv = ( c_double * len( rv ) )( *[ np.float64( _r ) for _r in rv ] )
        Xi_2h = ( c_double * len( rv ) )()
        size = c_uint( len( rv ) )
        lib_hm.model_Xi_2halo_chm( rv, Xi_2h, size, self.obj )
        
        return list( Xi_2h )

#################################################################################
#################################################################################
#################################################################################
#################################################################################
    

