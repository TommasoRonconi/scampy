########################################################################################

import os, sys
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

########################################################################################

def main () :

    extra_OMP_compile_args = []
    extra_OMP_link_args = []
    if sys.platform == 'darwin' :
        print( '------------------> Running on MacOS' )
        extra_OMP_compile_args += [  '-I/usr/local/opt/libomp/include', '-Xpreprocessor', '-fopenmp' ]
        extra_OMP_link_args += [ '-L/usr/local/opt/libomp/lib' ]

    ####################################################################################
    # Interpolation extension, currently providing the linear interpolator
    intrp_ext = Pybind11Extension(
        "scampy.utilities.interpolation",
        sorted(
            [ os.path.join( 'pybind11', 'pyb11_interpolation.cpp' ) ]
        ),
        include_dirs = sorted( [ os.path.join( 'c++', 'utilities', 'include' ),
                                 os.path.join( 'c++', 'interpolator', 'include' ) ] ),
        extra_compile_args=[ '-std=c++17' ]
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
        extra_compile_args=[ '-std=c++17' ]
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
        libraries = [ "m", "omp" ],
        extra_compile_args=['-std=c++17'] + extra_OMP_compile_args, # , '-fopenmp'],
        extra_link_args=extra_OMP_link_args+['-lomp']
    )

    ####################################################################################
    # Run setup
    setup( name = "scampy",
           version='0.0.1',
           package_dir = {
               'scampy' : 'scampy',
               'scampy.io' : os.path.join( 'scampy', 'io' ),
               'scampy.halo' : os.path.join( 'scampy', 'halo' ),
               'scampy.measure' : os.path.join( 'scampy', 'measure' ),
               'scampy.plot' : os.path.join( 'scampy', 'plot' ),
               'scampy.utilities' : os.path.join( 'scampy', 'utilities' )
           },
           packages = [ 'scampy',
                        'scampy.io',
                        'scampy.halo',
                        'scampy.measure',
                        'scampy.plot',
                        'scampy.utilities'
           ],
           cmdclass={"build_ext": build_ext},
           ext_modules = [ intrp_ext, cosmo_ext, clust_ext ],
           setup_requires = [ 'setuptools', 'pybind11' ],
           install_requires = [ 'numpy', 'scipy', 'matplotlib', 'h5py' ]
    )

########################################################################################

if __name__ == "__main__" : main()



    
