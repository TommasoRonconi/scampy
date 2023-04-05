########################################################################################

import os
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

########################################################################################

def main () :

    ####################################################################################
    # Interpolation extension, currently providing the linear interpolator
    intrp_ext = Pybind11Extension(
        "scampy.utilities.interpolation",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_interpolation.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'utilities', 'include' ),
                                 os.path.join( 'c++', 'interpolator', 'include' ) ] ),
    )

    ####################################################################################
    # Cosmology extension, provides mainly cosmographic functions and functions
    # depending on the parameters and Friedmann equations
    cosmo_ext = Pybind11Extension(
        "scampy.cosmology",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_cosmology.cpp' ),
              os.path.join( 'c++', 'cosmology', 'src', 'cosmological_model.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'utilities', 'include' ),
                                 os.path.join( 'c++', 'interpolator', 'include' ),
                                 os.path.join( 'c++', 'cosmology', 'include' ) ] ),
    )

    ####################################################################################
    # Clustering extension providing compiled and OMP-parallel
    # functions to compute distances
    clust_ext = Pybind11Extension(
        "scampy.measure.clustering_core",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_clustering_core.cpp' ),
              os.path.join( 'c++', 'utilities', 'src', 'clustering_core.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'utilities', 'include' ) ] ),
        libraries = [ "m", "gomp" ],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-lgomp'],
    )

    ####################################################################################
    # Run setup
    setup( name = "scampy",
           package_dir = {
               'scampy' : 'scampy',
               'scampy.halo' : os.path.join( 'scampy', 'halo' ),
               'scampy.io' : os.path.join( 'scampy', 'io' ),
               'scampy.utilities' : os.path.join( 'scampy', 'utilities' ),
               'scampy.measure' : os.path.join( 'scampy', 'measure' )
           },
           packages = [ 'scampy',
                        'scampy.halo',
                        'scampy.io',
                        'scampy.utilities',
                        'scampy.measure'
           ],
           ext_modules = [ intrp_ext, cosmo_ext, clust_ext ],
           install_requires = [ 'numpy', 'scipy' ]
    )

########################################################################################

if __name__ == "__main__" : main()



    
