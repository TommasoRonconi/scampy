import numpy
from .internal.cwrap import *

# ========================================================================================
# ========================================= LIN ==========================================
# ========================================================================================

# Wrap class interpolator< gsl_lin_interp >:
class lin_interpolator () :
    """ Class to handle linear interpolator objects.
    
    It defines a callable object to perform cubic spline interpolation
    on a linearly equispaced grid.
    
    Parameters
    ----------
    xv : array-like object
      should have the same dimension of fv.
      Defines the x-axis grid. Must be sorted and monothonic.
    fv : array-like object 
      Should have the same dimension of xv.
      Defines the y-axis values of the function to be interpolated.
    """

    @classmethod
    def from_function ( cls, func, xmin, xmax, thinness = 100 ) :
        """ class method to build object from function pointer
        
        Parameters
        ----------
        func : lambda
          function to interpolate
          Should take one argument (float, x-axis value) and return one argument
          (float, y-axis value)
        xmin : float
          lower limit of the interpolation grid
        xmax : float
          upper limit of the interpolation grid
          ( should be :math:`>x_\\text{min}` )
        thinness : int
          thinness of the interpolation grid (default = 100)
        
        Returns
        -------
        lin_interpolator
          Object of type lin_interpolator
        
        Warning
        -------
        thinness should always be greater than 10
        """

        xv = numpy.linspace( xmin, xmax, thinness )
        fv = numpy.array( [ func( _x ) for _x in xv ] )
        
        return cls( xv, fv )
    
    def __init__ ( self, xv, fv ) :
        """ Constructor of the class lin_interpolator
        """
        
        self._xmin = numpy.min( xv )
        self._xmax = numpy.max( xv )
        if len( xv ) == len( fv ) :
            self._thinness = len( xv )
        else :
            raise ValueError( "List xv and list fv must have the same length" )
        
        self._xv = ( c_double * len( xv ) )( *[ _x for _x in xv ] )
        self._fv = ( c_double * len( fv ) )( *[ _f for _f in fv ] )

        self.obj = lib_intrp.create_lin_interpolator( self._xv, self._fv,
                                                      c_size_t( self._thinness ) )

    def __del__ ( self ) :
        """ Destructor of the class lin_interpolator
        """

        # Python call to lin_interpolator dtor:
        lib_intrp.free_lin_interpolator( self.obj )

    def __call__ ( self, xx ) :
        """ Overload of the function-call operator
        
        Parameters
        ----------
        xx : float
          x-axis value in which to compute interpolated value of function
        
        Returns
        -------
        float
          :math:`f(x)` value of the interpolated function computed at position xx 
        
        Warning
        -------
        xx must be inside the interpolation domain
        """

        return lib_intrp.lin_interpolator_eval( c_double( xx ), self.obj )

    def get_xmin ( self ) :
        """ Get internal variable :math:`x_\\text{min}`
        
        Returns
        -------
        float
          lower limit of the interpolation domain
        """
        return self._xmin

    def get_xmax ( self ) :
        """ Get internal variable :math:`x_\\text{max}`
        
        Returns
        -------
        float
          upper limit of the interpolation domain
        """
        return self._xmax

    def integrate ( self, aa, bb ) :
        """ Integral of the interpolated function within interval
        
        It computes the integral
        
        .. math:: F(b) - F(a) = \\int_a^{b} f(x) d x

        where :math:`F(x)` is the integral of :math:`f(x)`.
        
        Parameters
        ----------
        aa : float
          lower limit of the integral
        bb : float
          upper limit of the integral
        
        Returns
        -------
        float
          integral of :math:`f(x)` within the defined integral 
        
        Note
        ----
        No numerical method for integration is used, the integral is obtained analytically 
        from the cubic polynomial that defines the cubic spline interpolation.
        
        Warning
        -------
        aa and bb must be inside the interpolation domain
        """

        return lib_intrp.lin_interpolator_integrate( c_double( aa ),
                                                     c_double( bb ),
                                                     self.obj )

# ========================================================================================
# ========================================= LOG ==========================================
# ========================================================================================

# Wrap class interpolator< gsl_log_interp >:
class log_interpolator () :
    """ Class to handle logarithm interpolator objects.
    
    It defines a callable object to perform cubic spline interpolation
    on a logarithmically equispaced grid.
    
    Parameters
    ----------
    xv : array-like object
      should have the same dimension of fv.
      Defines the x-axis grid. Must be sorted and monothonic.
      The system will assume the bin lenght is regular in the interval
      :math:`[\\log(x_0), \\log(x_N)]` where N is the lenght of the xv array
    fv : array-like object 
      Should have the same dimension of xv.
      Defines the y-axis values of the function to be interpolated.
    """

    @classmethod
    def from_function ( cls, func,
                        xmin = 1.e-7, xmax = 1.e+7,
                        thinness = 200 ) :
        """ class method to build object from function pointer
        
        Parameters
        ----------
        func : lambda
          function to interpolate
          Should take one argument (float, x-axis value) and return one argument
          (float, y-axis value)
        xmin : float
          lower limit of the interpolation grid 
          ( should be :math:`>0`, default = :math:`10^{-7}` )
        xmax : float
          upper limit of the interpolation grid
          ( should be :math:`>x_\\text{min}`, default = :math:`10^7` )
        thinness : int
          thinness of the interpolation grid (default = 200)
        
        Returns
        -------
        log_interpolator
          Object of type log_interpolator
        
        Warning
        -------
        thinness should always be greater than 10
        """

        xv = numpy.logspace( numpy.log10( xmin ), numpy.log10( xmax ), thinness )
        fv = numpy.array( [ func( _x ) for _x in xv ] )
        
        return cls( xv, fv )
    
    def __init__ ( self, xv, fv ) :
        """ Constructor of the class log_interpolator
        """
        
        self._xmin = numpy.min( xv )
        self._xmax = numpy.max( xv )
        if len( xv ) == len( fv ) :
            self._thinness = len( xv )
        else :
            raise ValueError( "List xv and list fv must have the same length" )
        
        self._xv = ( c_double * len( xv ) )( *[ _x for _x in xv ] )
        self._fv = ( c_double * len( fv ) )( *[ _f for _f in fv ] )

        self.obj = lib_intrp.create_log_interpolator( self._xv, self._fv,
                                                      c_size_t( self._thinness ) )

    def __del__ ( self ) :
        """ Destructor of the class log_interpolator
        """

        # Python call to log_interpolator dtor:
        lib_intrp.free_log_interpolator( self.obj )

    def __call__ ( self, xx ) :
        """ Overload of the function-call operator
        
        Parameters
        ----------
        xx : float
          x-axis value in which to compute interpolated value of function
        
        Returns
        -------
        float
          :math:`f(x)` value of the interpolated function computed at position xx 
        
        Warning
        -------
        xx must be inside the interpolation domain
        """

        return lib_intrp.log_interpolator_eval( c_double( xx ), self.obj )

    def get_xmin ( self ) :
        """ Get internal variable :math:`x_\\text{min}`
        
        Returns
        -------
        float
          lower limit of the interpolation domain
        """
        return self._xmin

    def get_xmax ( self ) :
        """ Get internal variable :math:`x_\\text{max}`
        
        Returns
        -------
        float
          upper limit of the interpolation domain
        """
        return self._xmax

    def integrate ( self, aa, bb ) :
        """ Integral of the interpolated function within interval
        
        It computes the integral
        
        .. math:: F(b) - F(a) = \\int_a^{b} f(x) d x

        where :math:`F(x)` is the integral of :math:`f(x)`.
        
        Parameters
        ----------
        aa : float
          lower limit of the integral
        bb : float
          upper limit of the integral
        
        Returns
        -------
        float
          integral of :math:`f(x)` within the defined integral 
        
        Note
        ----
        No numerical method for integration is used, the integral is obtained analytically 
        from the cubic polynomial that defines the cubic spline interpolation.
        
        Warning
        -------
        aa and bb must be inside the interpolation domain
        """

        return lib_intrp.log_interpolator_integrate( c_double( aa ),
                                                     c_double( bb ),
                                                     self.obj )


if __name__ == '__main__' :

    def lin_func ( xx ) :
        return xx

    lin_func_f = lin_interpolator.from_function( lin_func, 1., 10., 10 )
    print( "Lin-eval: {:f}".format( lin_func_f( 5. ) ), "should be ~ 5" )

    xv = numpy.logspace( numpy.log10( 1. ), numpy.log10( 10. ), 10 )
    fv = numpy.array( [ lin_func( _x ) for _x in xv ] )    
    log_func_f = log_interpolator( xv, fv )
    print( "Log-eval: {:f}".format( log_func_f( 5. ) ), "should be ~ 5" )

    print( "Lin-integrate: {:f}".format( lin_func_f.integrate( 2., 4. ) ) )
    print( "Log-integrate: {:f}".format( log_func_f.integrate( 2., 4. ) ) )

    from scipy.integrate import quad
    integral, error = quad( lin_func, 2., 4. )
    print( "Scipy.integrate quad: {:f}".format( integral ) )
