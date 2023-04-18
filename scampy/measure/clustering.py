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

def _kernel_standard ( DD, RR ) :

    if DD.size != RR.size :
        raise ValueError( 'distance vectors have different sizes' )

    ww = RR > 0
    out = numpy.zeros_like( ww, dtype = float )
    out[ww] = DD[ww] / RR[ww] - 1

    return out

##################################################################################

def two_point_standard ( data, rand, rbins, omp = True ) :
    """
    """

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    if not NobjD > 0 or not NobjR > 0 :
        raise ValueError(
            "Cannot compute clustering if one of the two catalogues does not have at least 2 elements"
        )

    normDD = 2.0 / ( NobjD * ( NobjD - 1 ) )
    normRR = 2.0 / ( NobjR * ( NobjR - 1 ) )

    DD = _kernel_DD( data, NdimD, rbins, omp ) * normDD
    RR = _kernel_DD( rand, NdimD, rbins, omp ) * normRR

    return _kernel_standard( DD, RR )
    
##################################################################################

def _kernel_landy_szalay ( DD, RR, DR ) :

    if DD.size != RR.size or RR.size != DR.size :
        raise ValueError( 'distance vectors have different sizes' )

    ww = RR > 0
    out = numpy.zeros_like( ww, dtype = float )
    out[ww] = ( DD[ww] - 2.0 * DR[ww] ) / RR[ww] + 1

    return out
    
##################################################################################

def _kernel_error_landy_szalay () :
    pass
    
##################################################################################

def two_point_landyszalay ( data, rand, rbins, omp = True, return_error = False ) :
    """
    """

    NdimD, NobjD = data.shape
    NdimR, NobjR = rand.shape
    if NdimD != NdimR :
        raise ValueError( "Input spaces ``data`` and ``rand`` should have the same dimensions" )
    if not NobjD > 0 or not NobjR > 0 :
        raise ValueError(
            "Cannot compute clustering if one of the two catalogues does not have at least 2 elements"
        )

    normDD = 2.0 / ( NobjD * ( NobjD - 1 ) )
    normRR = 2.0 / ( NobjR * ( NobjR - 1 ) )
    normDR = 1.0 / ( NobjD * NobjR )

    DD = _kernel_DD( data, NdimD, rbins, omp )
    DDn = DD * normDD
    RR = _kernel_DD( rand, NdimD, rbins, omp )
    RRn = RR * normRR
    DR = _kernel_DR( data, rand, NdimD, rbins, omp )
    DRn = DR * normDR

    # compute baseline clustering
    xi = _kernel_landy_szalay( DDn, RRn, DRn )
    
    if return_error :
        fact = (
            normRR / normDD * RR * (1.0 + xi ) +
            4.0 / NobjD * ( normRR * RR / normDD * ( 1.0 + xi ) )**2
        )
        exi = numpy.zeros_like( xi )
        ww = RRn > 0
        exi[ww] = normDD / RRn[ww] * numpy.sqrt( fact[ww] ) * numpy.sqrt( 3 )
        return xi, exi
    return xi

##################################################################################

def bootstrap_two_point ( data, rand, rbins,
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
    if not NobjD > 0 or not NobjR > 0 :
        raise ValueError(
            "Cannot compute clustering if one of the two catalogues does not have at least 2 elements"
        )

    normDD = 2.0 / ( NobjD * ( NobjD - 1 ) )
    normRR = 2.0 / ( NobjR * ( NobjR - 1 ) )
    normDR = 1.0 / ( NobjD * NobjR )
    
    if verbose : print( 'Computing RR ...' )
    RR = _kernel_DD( rand, NdimR, rbins, omp ) * normRR
    if verbose : print( '... done RR.' )

    # get the baseline estimate
    if verbose : print( 'Computing baseline ...' )
    DD = _kernel_DD( data, NdimD, rbins, omp ) * normDD
    if standard :
        tpt = _kernel_standard( DD, RR )
    else :
        DR = _kernel_DR( data, rand, NdimD, rbins, omp ) * normDR
        tpt = _kernel_landy_szalay( DD, RR, DR )
    if verbose : print( '... done baseline.' )

    # get bootstraps
    boots = numpy.zeros( ( Nboots, tpt.size ) )
    for ii in range( Nboots ):
        if verbose : print( f'Computing bootstrap {ii+1:d} of {Nboots} ...', end='\r' )
        iD = rng.integers( NobjD, size = NobjD )
        DD = _kernel_DD( data[:,iD], NdimD, rbins, omp ) * normDD
        if standard :
            boots[ii] = _kernel_standard( DD, RR )
        else :
            DR = _kernel_DR( data[:,iD], rand, NdimD, rbins, omp ) * normDR
            boots[ii] = _kernel_landy_szalay( DD, RR, DR )
    if verbose : print( '\n... done bootstraps.' )

    if return_boots:
        return tpt, boots.std(axis=0), boots
    else:
        return tpt, boots.std(axis=0)

