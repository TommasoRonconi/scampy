#######################################################################
# Global project info

project( 'fftlog',
         'fortran',
         default_options : [ 'buildtype=release',
                             'warning_level=0' ] )

root = meson.source_root()

fc = meson.get_compiler( 'fortran' )

#######################################################################
# Make library

fftlog_lib = static_library( 'fftlog',
                             [ 'cdgamma.f',
                               'drfftb.f',
                               'drfftf.f',
                               'drffti.f',
                               'fftlog.f' ],
                             install : false )

fftlog_dep = declare_dependency( link_with : fftlog_lib )

#######################################################################

