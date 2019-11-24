import numpy
from .cwrap.cwrap import *

class occupation_p :

    # _param = None
    @classmethod
    def load () :
        pass

    def __init__ ( self ) :
        pass

    def Ncen ( self, *args, **kwargs ) :
        pass

    def Nsat ( self, *args, **kwargs ) :
        pass

    def save ( self, out_file ) :
        pass

class harikane16_p ( occupation_p ) :
    
    @classmethod
    def load ( cls, in_file ) :

        param = numpy.load( in_file )
        if len( param ) != 6 :
            raise RuntimeError( "The input parameter file passed ( " +
                                in_file +
                                " ) has the wrong number of parameters ( {:d}, should be 6 )".format( len( param ) ) )
        
        return cls( *param )

    def __init__ ( self,
                   DC = 0.5,
                   M_min = 1.e+11,
                   sigma_logM = 1.,
                   M0 = 1.e+10,
                   M1 = 1.e+10,
                   alpha = 1. ) :

        self._DC = c_double( DC )
        self._M_min = c_double( M_min )
        self._sigma_logM = c_double( sigma_logM )
        self._M0 = c_double( M0 )
        self._M1 = c_double( M1 )
        self._alpha = c_double( alpha )

        # Python call to occupation ctor:
        self.obj = lib_ocp.create_H16_occupation( self._DC, self._M_min,
                                                  self._sigma_logM, self._M0,
                                                  self._M1, self._alpha )
        
    def __del__ ( self ) :

        # Python call to occupation dtor:
        lib_ocp.free_H16_occupation( self.obj )

    # def __call__ ( self, idx ) :

    #     return _param[ idx ]

    def set_parameters ( self,
                         DC = None,
                         M_min = None,
                         sigma_logM = None,
                         M0 = None,
                         M1 = None,
                         alpha = None ) :

        self._DC = c_double( DC ) if DC is not None else self._DC
        self._M_min = c_double( M_min ) if M_min is not None else self._M_min
        self._sigma_logM = c_double( sigma_logM ) if sigma_logM is not None else self._sigma_logM
        self._M0 = c_double( M0 ) if M0 is not None else self._M0
        self._M1 = c_double( M1 ) if M1 is not None else self._M1
        self._alpha = c_double( alpha ) if alpha is not None else self._alpha
        
        self.obj = lib_ocp.create_H16_occupation( self._DC, self._M_min,
                                                  self._sigma_logM, self._M0,
                                                  self._M1, self._alpha )
        return
    
    def Ncen ( self, Mh ) :

        return lib_ocp.Ncen_H16_ocp( c_double( Mh ), self.obj )

    def Nsat ( self, Mh ) :

        return lib_ocp.Nsat_H16_ocp( c_double( Mh ), self.obj )

    def save ( self, out_file ) :

        numpy.save( out_file, numpy.array( [ self._DC.value,
                                             self._M_min.value,
                                             self._sigma_logM.value,
                                             self._M0.value,
                                             self._M1.value,
                                             self._alpha.value ] ) )

        return;

class tinker10_p ( occupation_p ) :

    # _param = { 'Amin' : _Amin,
    #            'siglogA' : _siglogA,
    #            'Asat' : _Asat,
    #            'alpsat' : _alpsat }
    
    @classmethod
    def load ( cls, in_file ) :

        param = numpy.load( in_file )
        if len( param ) != 4 :
            raise RuntimeError( "The input parameter file passed ( " +
                                in_file +
                                " ) has the wrong number of parameters ( {:d}, should be 4 )".format( len( param ) ) )
        
        return cls( *param )
    
    def __init__ ( self,
                   Amin = 1.e+11,
                   siglogA = 1.,
                   Asat = 1.e+10,
                   alpsat = 1. ) :

        self._Amin = c_double( Amin )
        self._siglogA = c_double( siglogA )
        self._Asat = c_double( Asat )
        self._alpsat = c_double( alpsat )

        # Python call to occupation ctor:
        self.obj = lib_ocp.create_T10_occupation( self._Amin, self._siglogA,
                                                  self._Asat, self._alpsat )
        
    def __del__ ( self ) :

        # Python call to occupation dtor:
        lib_ocp.free_T10_occupation( self.obj )

    # def __call__ ( self, idx ) :

    #     return _param[ idx ]

    def set_parameters ( self,
                         Amin = None,
                         siglogA = None,
                         Asat = None,
                         alpsat = None ) :

        self._Amin = c_double( Amin ) if Amin is not None else self._Amin
        self._siglogA = c_double( siglogA ) if siglogA is not None else self._siglogA
        self._Asat = c_double( Asat ) if Asat is not None else self._Asat
        self._alpsat = c_double( alpsat ) if alpsat is not None else self._alpsat
        
        self.obj = lib_ocp.create_T10_occupation( self._Amin, self._siglogA,
                                                  self._Asat, self._alpsat )
        return
    
    def Ncen ( self, Mh ) :

        return lib_ocp.Ncen_T10_ocp( c_double( Mh ), self.obj )

    def Nsat ( self, Mh ) :

        return lib_ocp.Nsat_T10_ocp( c_double( Mh ), self.obj )

    def save ( self, out_file ) :

        numpy.save( out_file, numpy.array( [ self._Amin.value,
                                             self._siglogA.value,
                                             self._Asat.value,
                                             self._alpsat.value ] ) )

        return;
