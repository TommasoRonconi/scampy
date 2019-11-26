import numpy
from .cwrap.cwrap import *

# ========================================================================================
# ========================================= LIN ==========================================
# ========================================================================================

# Wrap class interpolator< gsl_lin_interp >:
class lin_interpolator () :

    @classmethod
    def from_function ( cls, func, xmin, xmax, thinness = 100 ) :

        xv = numpy.linspace( xmin, xmax, thinness )
        fv = numpy.array( [ func( _x ) for _x in xv ] )
        
        return cls( xv, fv )
    
    def __init__ ( self, xv, fv ) :
        
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

        # Python call to lin_interpolator dtor:
        lib_intrp.free_lin_interpolator( self.obj )

    def __call__ ( self, xx ) :

        return lib_intrp.lin_interpolator_eval( c_double( xx ), self.obj )

    def get_xmin ( self ) :
        return self._xmin

    def get_xmax ( self ) :
        return self._xmax

    def integrate ( self, aa, bb ) :

        return lib_intrp.lin_interpolator_integrate( c_double( aa ),
                                                     c_double( bb ),
                                                     self.obj )

# ========================================================================================
# ========================================= LOG ==========================================
# ========================================================================================

# Wrap class interpolator< gsl_log_interp >:
class log_interpolator () :

    @classmethod
    def from_function ( cls, func,
                        xmin = 1.e-7, xmax = 1.e+7,
                        thinness = 200 ) :

        xv = numpy.logspace( numpy.log10( xmin ), numpy.log10( xmax ), thinness )
        fv = numpy.array( [ func( _x ) for _x in xv ] )
        
        return cls( xv, fv )
    
    def __init__ ( self, xv, fv ) :
        
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

        # Python call to log_interpolator dtor:
        lib_intrp.free_log_interpolator( self.obj )

    def __call__ ( self, xx ) :

        return lib_intrp.log_interpolator_eval( c_double( xx ), self.obj )

    def get_xmin ( self ) :
        return self._xmin

    def get_xmax ( self ) :
        return self._xmax

    def integrate ( self, aa, bb ) :

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