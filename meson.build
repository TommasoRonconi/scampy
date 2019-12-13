#######################################################################
# Global project info

project( 'master',
         'cpp',
         default_options : [ 'cpp_std=c++14',
                             'buildtype=release',
                             'prefix=' + meson.source_root() + '/install',
                             'libdir=lib' ],
         version : '0.1.0',
         license : 'SISSA' )

root = meson.source_root()

cxx = meson.get_compiler( 'cpp' )

cxx_args = [ '-Wall',
             '-Wextra',
             '-fopenmp' ]

foreach cxx_arg : cxx_args
if cxx.has_argument( cxx_arg ) 
add_global_arguments( cxx_arg,
                      language : 'cpp' )
endif
endforeach

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
# Internal includes

lib_headers = []

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

install_headers( lib_headers )

#######################################################################
# Python module

subdir( 'python') 

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
