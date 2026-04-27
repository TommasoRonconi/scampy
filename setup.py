########################################################################################

import os
from setuptools import setup, find_packages

########################################################################################

def main () :

    setup( name = "scampy",
           version = "2.0.0",
           package_dir = {
               'scampy'              : 'scampy',
               'scampy.io'          : os.path.join( 'scampy', 'io' ),
               'scampy.halo'        : os.path.join( 'scampy', 'halo' ),
               'scampy.measure'     : os.path.join( 'scampy', 'measure' ),
               'scampy.plot'        : os.path.join( 'scampy', 'plot' ),
               'scampy.utilities'   : os.path.join( 'scampy', 'utilities' )
           },
           packages = [ 'scampy',
                        'scampy.io',
                        'scampy.halo',
                        'scampy.measure',
                        'scampy.plot',
                        'scampy.utilities'
           ],
           setup_requires = [ 'setuptools' ],
           install_requires = [ 'numpy', 'scipy', 'matplotlib', 'h5py' ]
    )

########################################################################################

if __name__ == "__main__" : main()
