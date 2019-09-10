#######################################################################
# Internal includes

inc_dir = 'include'

inc_fftlwrap = include_directories( inc_dir )

#######################################################################
# FFTLog Wrap Lib

fftlwrap_lib = library( 'fftlog',
                        [ 'src/fftlog_handler.cpp' ],
       	 	        include_directories : [ inc_fftlwrap,
                                                inc_intrp ],
		        dependencies : [ gsl_dep,
                                         fftlog_dep,
                                         utl_dep,
                                         intrp_dep ],
                        install : true )

fftlwrap_dep = declare_dependency( link_with : fftlwrap_lib,
                                   include_directories : [ inc_fftlwrap,
                                                           inc_intrp ] )

lib_headers += [ root + '/fftlog_wrap/' + inc_dir + '/fftlog_handler.h',
                 root + '/fftlog_wrap/' + inc_dir + '/abel_transform.h' ]

#######################################################################
