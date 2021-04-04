#######################################################################
# Global project info

project( 'scampy',
         'cpp',
         default_options : [ 'cpp_std=c++14',
                             'buildtype=release',
                             'warning_level=0',
                             'libdir=lib' ],
         version : '0.2.0',
         license : 'GPL-3.0' )

root = meson.source_root()

cxx = meson.get_compiler( 'cpp' )
cxx_args = [ '-fopenmp' ]
foreach cxx_arg : cxx_args
if cxx.has_argument( cxx_arg ) 
add_global_arguments( cxx_arg,
                      language : 'cpp' )
endif
endforeach

#######################################################################
# Search for python installation

pymod = import( 'python' )
py = pymod.find_installation()
message( 'Installing python sector in ' + py.get_install_dir() )

#######################################################################
# Set library directory path, depends on the type of build chosen
# (i.e. in case of full-build=true, installs the whole package
#  in developer mode, otherwise only installs the python package
#  with the c++ sector compiled static libraries and the c-wrap
#  shared objects all contained in the sub-package scampy.internal)

if get_option( 'full-build' )
libdir = '/'.join( [ get_option( 'prefix' ), get_option( 'libdir' ) ] )
else
libdir = '/'.join( [ py.get_install_dir(), 'scampy', 'internal' ] )
endif

#######################################################################
# External dependencies

# OMP
omp_dep = dependency( 'openmp',
                      required : true )

# GSL
gsl_dep = dependency( 'gsl',
                      required : true,
                      version : '>=2' )

# FFTLog
fftlog_dep = dependency( 'fftlog',
                         fallback : [ 'fftlog', 'fftlog_dep' ],
                         required : true )

#######################################################################
# Empty list to append headers (from both C++ and C sectors)

if get_option( 'full-build' )
lib_headers = []
endif

#######################################################################
# Utilities Lib

subdir( 'utilities' )

#######################################################################
# Interpolation Lib

subdir( 'interpolator' )

#######################################################################
# FFTLog Wrap Lib

subdir( 'fftlog_wrap' )

#######################################################################
# Cosmology Lib

subdir( 'cosmology' )

#######################################################################
# Halo-Model Lib

subdir( 'halo_model' )

#######################################################################
# Further installation commands

if get_option( 'full-build' )
install_headers( lib_headers )
endif

#######################################################################
# Python module

subdir( 'python' )

#######################################################################
# Testing

if get_option( 'enable-test' )
    subdir( 'tests' )
endif

#######################################################################
# Documentation

if get_option( 'enable-doc' )
    subdir( 'doc' )
endif
