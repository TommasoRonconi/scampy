import numpy

class halo ( object ) :
    
    def __init__ ( self, coords, mass, spin = None, veldisp = None, radius = None ) :
         
        self.centre = numpy.array( coords )
        self.mass = mass
        self.spin = numpy.array( spin )
        self.veldisp = veldisp
        self.radius = radius
        
    def __repr__ ( self ) :
        return 'Halo(centre = {0:s}, mass = {1:f})'.format( self.centre.__repr__(), self.mass )

# REMEMBER THAT: a galaxy IS NOT a halo, but a halo HAS a galaxy
#                ==> inheritance here is conceptually wrong!!
class galaxy ( halo ) :

    @classmethod
    def from_halo ( cls, halo, luminosity = None ) :
        return cls( coords = halo.centre,
                    mass = halo.mass,
                    spin = halo.spin,
                    veldisp = halo.veldisp,
                    radius = halo.radius,
                    luminosity = luminosity )
        
    
    def __init__ ( self,
                   coords, mass,
                   spin = None,
                   veldisp = None,
                   radius = None,
                   luminosity = None ) :

        self.luminosity = luminosity
        
        halo.__init__( self,
                       coords, mass,
                       spin = None,
                       veldisp = None,
                       radius = None )

    # def __repr__ ( self ) :
    #     return 'Galaxy({0:s}, luminosity = {1:f})'.format( super( galaxy, self ).__repr__() )
    #, self.luminosity )
        

class host_halo ( object ) :
    
    def __init__ ( self, coords, mass, central = numpy.array([]), satellites = numpy.array([]) ) :
        
        self.coords = numpy.array( coords )
        self.mass = mass
        self.central = central
        self.Ncen = len( self.central )
        self.satellites = satellites
        self.Nsat = len( self.satellites )
        
    def __repr__ ( self ) :
        out = 'Host-Halo(coords = {0:s}, mass = {1:f}, Ncen = {2:d}, Nsat = {3:d})'
        return out.format( repr( self.coords ), self.mass, self.Ncen, self.Nsat )
    
    def set_central ( self, central = numpy.array([]) ) :
        self.central = central
        self.Ncen = len( central )
        return
    
    def set_satellites ( self, satellites = numpy.array([]) ) :
        self.satellites = satellites
        self.Nsat = len(self.satellites)
        return

    def mask ( self, mask ) :
        """ Returns masked copy of self
        
        Parameters
        ----------
        mask : array-like of booleans
          The mask to be applied. The first element of the array masks the central 
          sub-halo, the remaining mask the satellite sub-haloes.
          The shape of the array should be ( Ncen + Nsat, )
        
        Returns
        -------
        : host_halo object or None
          If any of the elements of mask is True -> masked host_halo
          (with same coordinates and mass of the original, and masked sub-haloes).
          If all the elements of the mask are False -> None
        """
        if numpy.any( mask ) :
            return host_halo( coords = self.coords,
                              mass = self.mass,
                              central = self.central[
                                  mask[ 0:len( self.central ) ] ],
                              satellites = self.satellites[
                                  mask[ self.Ncen:self.Ncen+self.Nsat ] ] )
        else :
            return None
                          


if __name__ == '__main__' :

    hh = halo( [ 1., 2., 3. ],
               10.,
               spin = [ 1., 1., 1. ],
               radius = 7.,
               veldisp = 0.5 )
    print( hh )

    gxy = galaxy.from_halo( hh, luminosity = -20. )
    print( gxy )

    host = host_halo( coords = [1.,2.,3.],
                      mass=15. )
    print( host )

    host.set_central( numpy.array( [ hh ] ) )
    print( host )
