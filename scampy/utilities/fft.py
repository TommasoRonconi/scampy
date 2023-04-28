""" 
"""

#############################################################################################

# External imports
from scipy import fft
import numpy

#############################################################################################

# Fast sin-transform log (direct)
def fstl ( xn, yn, lk0 = 0.0, bias = 0.0, mu = 0.5 ) :

    # copy input and prep for transform
    lxn = numpy.log( xn )
    yn = numpy.array( yn ).T    
    # log and copy input
    lxn = numpy.log( xn )
    
    # Find log-spacing and check equal
    dlr = None
    try :
        # Note that the `numpy.around` function is less
        # accurate than the equivalent python built-in 
        # function `round` 
        # (rounds to the closest machine precision float)
        dlr = numpy.unique( 
            numpy.around( lxn[1:]-lxn[:-1], 4 ) 
        ).item()
    except :
        raise
        
    # Prepare input function for transform
    fct = 0.5 * numpy.sqrt( 0.5 * numpy.pi ) / numpy.pi**2
    an  = fct * yn * xn * numpy.sqrt( xn )
        
    # Find optimal centre for transform
    kr = fft.fhtoffset( dlr, initial = lk0, 
                        mu = mu, bias = bias )
    
    # Compute wavenumber grid
    kn = numpy.exp( kr - lxn[::-1] )
    
    # Compute sin transform of an
    bn = fft.fht( an, dlr, mu = mu, bias = bias, offset = kr )
    
    # Return a tuple with
    # - wavenumber grid
    # - converted sin-transform of input function
    return kn, ( bn / ( kn * numpy.sqrt( kn ) ) ).T

#############################################################################################

# Fast projected sin-transform log (direct)
def fpstl ( xn, yn, lk0 = 0.0, bias = 0.0, mu = 0.0 ) :
    
    # copy input and prep for transform
    lxn = numpy.log( xn )
    yn = numpy.array( yn ).T
    
    # Find log-spacing and check equal
    dlr = None
    try :
        # Note that the `numpy.around` function is less
        # accurate than the equivalent python built-in 
        # function `round` 
        # (rounds to the closest machine precision float)
        dlr = numpy.unique( 
            numpy.around( lxn[1:]-lxn[:-1], 4 ) 
        ).item()
    except :
        raise
        
    # Prepare input function for transform
    fct = 0.5  / numpy.pi
    an  = numpy.array( fct * yn * xn )
    
    # Find optimal centre for transform
    kr = fft.fhtoffset( dlr, initial = lk0, 
                        mu = mu, bias = bias )
    
    # Compute wavenumber grid
    kn = numpy.exp( kr - lxn[::-1] )
    
    # Compute sin transform of an
    bn = fft.fht( an, dlr, mu = mu, bias = bias, offset = kr )
    
    # Return a tuple with
    # - wavenumber grid
    # - converted sin-transform of input function
    return kn, ( bn / kn ).T

#############################################################################################
