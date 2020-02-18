import numpy
from .cwrap.cwrap import *

class occupation_p :
    """ Base occupation_p class, it defines the conceptually-"virtual" methods that have
    to be overloaded in defining a derived occupation_p-type class.
    """

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
    """ c++ -wrapped 6-parameters occupation probability distribution (as in Harikane et al., 2016).
    It adds a linear :code:`DC` parameter to the classical 5-parameters HOD model 
    :math:`\\Rightarrow` by setting this additional parameter to :math:`1` it reduces to the
    standard 5-parameters model.
    The shapes of the occupation probabilities are:
    
    - average number of central galaxies per halo of mass :math:`M_h`:

    .. math:: \\langle N_\\text{cen}\\rangle ( M_h ) \\equiv \\text{DC} \\dfrac{1}{2} \\biggl[
              1 + \\text{erf}\\biggl(\\dfrac{\\log M_h - \\log M_\\text{min}}{\\sigma_{\\log M}}
              \\biggr)\\biggr]

    - average number of satellite galaxies per halo of mass :math:`M_h`:

    .. math:: \\langle N_\\text{sat}\\rangle ( M_h ) \\equiv \\text{DC} \\biggl(
              \\dfrac{M_h - M_0}{M_1}
              \\biggr)^\\alpha    

    Parameters
    ----------
    DC : scalar
       *duty cycle* parameter (:math:`DC`).
    M_min : scalar
       *cut-off mass* parameter (:math:`M_\\text{min}`)
    sigma_logM : scalar
       *transition phase* parameter (:math:`\\sigma_{\\log M}`)
    M0 : scalar
       *satellite mass bias* parameter (:math:`M_0`)
    M1 : scalar
       *satellite mass normalization* parameter (:math:`M_1`)
    alpha : scalar
       *spectral index* parameter (:math:`alpha`)
    """
    
    @classmethod
    def load ( cls, in_file ) :
        """ class-method for building the object from a :code:`.npy` file containing 
        an ordered array with the parameters of the model.

        Parameters
        ----------
        in_file : string
           string with :code:`/path/to/input_file.npy`
        
        Returns
        -------
        object of harikane16_p type
        
        Warning
        -------
        Raises a runtime error if the number of parameters in the file is :math:`\\neq 6`
        """

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
        """ Objects' denstructor. It calls the C-mangled destructor of the class.
        (it frees the heap memory slot in which the object has been stored)
        """

        # Python call to occupation dtor:
        lib_ocp.free_H16_occupation( self.obj )

    def set_parameters ( self,
                         DC = None,
                         M_min = None,
                         sigma_logM = None,
                         M0 = None,
                         M1 = None,
                         alpha = None ) :
        """ Convenience method for modifying the model parameters.
        It calls the C-mangled constructor of the object.    

        Parameters
        ----------
        DC : scalar or None
           *duty cycle* parameter (:math:`DC`). 
           If :code:`None` the old value is used.
        M_min : scalar or None
           *cut-off mass* parameter (:math:`M_\\text{min}`).
           If :code:`None` the old value is used.
        sigma_logM : scalar or None
           *transition phase* parameter (:math:`\\sigma_{\\log M}`).
           If :code:`None` the old value is used.
        M0 : scalar or None
           *satellite mass bias* parameter (:math:`M_0`).
           If :code:`None` the old value is used.
        M1 : scalar or None
           *satellite mass normalization* parameter (:math:`M_1`).
           If :code:`None` the old value is used.
        alpha : scalar or None
           *spectral index* parameter (:math:`alpha`).
           If :code:`None` the old value is used.
        """

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
        """ Function that returns the average number of central galaxies 
        per halo of mass :math:`M_h`:

        .. math:: \\langle N_\\text{cen}\\rangle ( M_h ) \\equiv \\text{DC} \\dfrac{1}{2} \\biggl[
                  1 + \\text{erf}\\biggl(\\dfrac{\\log M_h - \\log M_\\text{min}}{\\sigma_{\\log M}}
                  \\biggr)\\biggr]

        (It calls the C-mangled corresponding function)

        Parameters
        ----------
        Mh : scalar
           the mass :math:`M_h` of the host halo
        
        Returns
        -------
        scalar
           value of :math:`\\langle N_\\text{cen}\\rangle(M_h)`
        """

        return lib_ocp.Ncen_H16_ocp( c_double( Mh ), self.obj )

    def Nsat ( self, Mh ) :
        """ Function that returns the average number of satellite galaxies 
        per halo of mass :math:`M_h`:

        .. math:: \\langle N_\\text{sat}\\rangle ( M_h ) \\equiv \\text{DC} \\biggl(
                  \\dfrac{M_h - M_0}{M_1}
                  \\biggr)^\\alpha    

        (It calls the C-mangled corresponding function)

        Parameters
        ----------
        Mh : scalar
           the mass :math:`M_h` of the host halo
        
        Returns
        -------
        scalar
           value of :math:`\\langle N_\\text{sat}\\rangle(M_h)`
        """

        return lib_ocp.Nsat_H16_ocp( c_double( Mh ), self.obj )

    def save ( self, out_file ) :
        """ Convenience-method for saving the current parameterisation in :code:`.npy` format.
        It stores the current values of the defining parameters into a file :code:`out_file.npy`

        Paramters
        ---------
        out_file : string
           :code:`/path/to/output_file`, the :code:`numpy.save()` method is called

        Returns
        -------
        None
        """

        numpy.save( out_file, numpy.array( [ self._DC.value,
                                             self._M_min.value,
                                             self._sigma_logM.value,
                                             self._M0.value,
                                             self._M1.value,
                                             self._alpha.value ] ) )

        return;

class tinker10_p ( occupation_p ) :
    """ c++ -wrapped 4-parameters occupation probability distribution (as in Tinker et al., 2010).
    The shapes of the occupation probabilities are:
    
    - average number of central galaxies per halo of mass :math:`M_h`:

    .. math:: \\langle N_\\text{cen}\\rangle ( M_h ) \\equiv \\dfrac{1}{2} \\biggl[
              1 + \\text{erf}\\biggl(\\dfrac{\\log M_h - \\log M_\\text{min}}{\\sigma_{\\log M}}
              \\biggr)\\biggr]

    - average number of satellite galaxies per halo of mass :math:`M_h`:

    .. math:: \\langle N_\\text{sat}\\rangle ( M_h ) \\equiv \\biggl(
              \\dfrac{M_h}{M_\\text{sat}}
              \\biggr)^\\alpha    

    Parameters
    ----------
    Amin : scalar
       *cut-off mass* parameter (:math:`M_\\text{min}`)
    siglogA : scalar
       *transition phase* parameter (:math:`\\sigma_{\\log M}`)
    Asat : scalar
       *satellite mass normalization* parameter (:math:`M_\\text{sat}`)
    alpsat : scalar
       *spectral index* parameter (:math:`alpha`)
    """
    
    @classmethod
    def load ( cls, in_file ) :
        """ class-method for building the object from a :code:`.npy` file containing 
        an ordered array with the parameters of the model.

        Parameters
        ----------
        in_file : string
           string with :code:`/path/to/input_file.npy`
        
        Returns
        -------
        object of tinker10_p type
        
        Warning
        -------
        Raises a runtime error if the number of parameters in the file is :math:`\\neq 4`
        """

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

    def set_parameters ( self,
                         Amin = None,
                         siglogA = None,
                         Asat = None,
                         alpsat = None ) :
        """ Convenience method for modifying the model parameters.
        It calls the C-mangled constructor of the object.    

        Parameters
        ----------
        Amin : scalar or None
           *cut-off mass* parameter (:math:`M_\\text{min}`).
           If :code:`None` the old value is used.
        siglogA : scalar or None
           *transition phase* parameter (:math:`\\sigma_{\\log M}`).
           If :code:`None` the old value is used.
        Asat : scalar or None
           *satellite mass normalization* parameter (:math:`M_\\text{sat}`).
           If :code:`None` the old value is used.
        alpsat : scalar or None
           *spectral index* parameter (:math:`alpha`).
           If :code:`None` the old value is used.
        """

        self._Amin = c_double( Amin ) if Amin is not None else self._Amin
        self._siglogA = c_double( siglogA ) if siglogA is not None else self._siglogA
        self._Asat = c_double( Asat ) if Asat is not None else self._Asat
        self._alpsat = c_double( alpsat ) if alpsat is not None else self._alpsat
        
        self.obj = lib_ocp.create_T10_occupation( self._Amin, self._siglogA,
                                                  self._Asat, self._alpsat )
        return
    
    def Ncen ( self, Mh ) :
        """ Function that returns the average number of central galaxies 
        per halo of mass :math:`M_h`:

        .. math:: \\langle N_\\text{cen}\\rangle ( M_h ) \\equiv \\dfrac{1}{2} \\biggl[
                  1 + \\text{erf}\\biggl(\\dfrac{\\log M_h - \\log M_\\text{min}}{\\sigma_{\\log M}}
                  \\biggr)\\biggr]

        (It calls the C-mangled corresponding function)

        Parameters
        ----------
        Mh : scalar
           the mass :math:`M_h` of the host halo
        
        Returns
        -------
        scalar
           value of :math:`\\langle N_\\text{cen}\\rangle(M_h)`
        """

        return lib_ocp.Ncen_T10_ocp( c_double( Mh ), self.obj )

    def Nsat ( self, Mh ) :
        """ Function that returns the average number of satellite galaxies 
        per halo of mass :math:`M_h`:

        .. math:: \\langle N_\\text{sat}\\rangle ( M_h ) \\equiv \\biggl(
                  \\dfrac{M_h}{M_\\text{sat}}
                  \\biggr)^\\alpha    

        (It calls the C-mangled corresponding function)

        Parameters
        ----------
        Mh : scalar
           the mass :math:`M_h` of the host halo
        
        Returns
        -------
        scalar
           value of :math:`\\langle N_\\text{sat}\\rangle(M_h)`
        """

        return lib_ocp.Nsat_T10_ocp( c_double( Mh ), self.obj )

    def save ( self, out_file ) :
        """ Convenience-method for saving the current parameterisation in :code:`.npy` format.
        It stores the current values of the defining parameters into a file :code:`out_file.npy`

        Paramters
        ---------
        out_file : string
           :code:`/path/to/output_file`, the :code:`numpy.save()` method is called

        Returns
        -------
        None
        """

        numpy.save( out_file, numpy.array( [ self._Amin.value,
                                             self._siglogA.value,
                                             self._Asat.value,
                                             self._alpsat.value ] ) )

        return;
