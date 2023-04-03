##################################################################################
# External imports
import numpy

#Internal imports
import scampy.measure.clustering_core as cc

##################################################################################

def _kernel_DD (data, Nd, rbins, omp = True) :
    
    if omp :
        if Nd == 2 :
            return numpy.array( cc.d2D_DD_omp( *data, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DD_omp( *data, rbins ) )
    else :
        if Nd == 2 :
            return numpy.array( cc.d2D_DD( *data, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DD( *data, rbins ) )
            
    return None

def _kernel_DR (data1, data2, Nd, rbins, omp = True) :
    
    if omp :
        if Nd == 2 :
            return numpy.array( cc.d2D_DR_omp( *data1, *data2, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DR_omp( *data1, *data2, rbins ) )
    else :
        if Nd == 2 :
            return numpy.array( cc.d2D_DR( *data1, *data2, rbins ) )
        if Nd == 3 :
            return numpy.array( cc.d3D_DR( *data1, *data2, rbins ) )
            
    return None

##################################################################################

def two_point_standard ( data, rand, rbins, dfact = 1., omp = True ) :
    """
    """

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    
    fact = dfact * NobjR / NobjD
    DD = _kernel_DD( data, NdimD, rbins, omp ) / ( NdimD * ( NdimD - 1 ) )
    RR = _kernel_DD( rand, NdimD, rbins, omp ) / ( NdimR * ( NdimR - 1 ) )

    return DD / RR - 1
    #return fact * fact * DD / RR - 1

##################################################################################

def two_point_landyszalay ( data, rand, rbins, dfact = 1., omp = True ) :
    """
    """

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    
    fact = dfact * NobjR / NobjD
    DD = _kernel_DD( data, NdimD, rbins, omp ) / ( NdimD * ( NdimD - 1 ) )
    RR = _kernel_DD( rand, NdimD, rbins, omp ) / ( NdimR * ( NdimR - 1 ) )
    DR = _kernel_DR( data, rand, NdimD, rbins, omp ) / ( NdimD * NdimR )

    return ( DD - 2.0 * DR ) / RR + 1
    # return fact * ( fact * DD - 2.0 * DR ) / RR + 1

##################################################################################

def bootstrap_two_point ( data, rand, rbins,
                          dfact = 1.,
                          standard = True,
                          Nboots = 10, return_boots = False,
                          omp = True, verbose = True,
                          rng = None, kw_rng = { 'seed' : 555 } ) :
    """
    """

    if rng is None :
        rng = numpy.random.default_rng( **kw_rng )

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    
    fact = dfact * NobjR / NobjD

    if verbose : print( 'Computing RR ...' )
    RR = _kernel_DD( rand, NdimR, rbins, omp ) / ( NdimR * ( NdimR - 1 ) )
    if verbose : print( '... done RR.' )

    # get the baseline estimate
    if verbose : print( 'Computing baseline ...' )
    DD = _kernel_DD( data, NdimD, rbins, omp ) / ( NdimD * ( NdimD - 1 ) )
    if standard :
        tpt = DD / RR - 1
        # tpt = fact * fact * DD / RR - 1
    else :
        DR = _kernel_DR( data, rand, NdimD, rbins, omp ) / ( NdimD * NdimR )
        tpt = ( DD - 2.0 * DR ) / RR + 1
        # tpt = fact * ( fact * DD - 2.0 * DR ) / RR + 1
    if verbose : print( '... done baseline.' )

    # get bootstraps
    if verbose : print( 'Computing bootstraps ...' )
    boots = numpy.zeros( ( Nboots, tpt.size ) )
    for ii in range( Nboots ):
        iD = rng.integers( NobjD, size = NobjD )
        DD = _kernel_DD( data[:,iD], NdimD, rbins, omp ) / ( NdimD * ( NdimD - 1 ) )
        if standard :
            boots[ii] = DD / RR - 1
            # boots[ii] = fact * fact * DD / RR - 1
        else :
            DR = _kernel_DR( data[:,iD], rand, NdimD, rbins, omp ) / ( NdimD * NdimR )
            boots[ii] = ( DD - 2.0 * DR ) / RR + 1
            # boots[ii] = fact * ( fact * DD - 2.0 * DR ) / RR + 1
    if verbose : print( '... done bootstraps.' )

    if return_boots:
        return tpt, boots.std(axis=0), boots
    else:
        return tpt, boots.std(axis=0)
