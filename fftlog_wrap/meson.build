#######################################################################
# Internal includes

inc_dir = 'include'

inc_fftlwrap = include_directories( inc_dir )

#######################################################################
# FFTLog Wrap Lib

if get_option( 'full-build' )
fftlwrap_lib = library( 'fftlog_wrap',
                        [ 'src/fftlog_handler.cpp' ],
       	 	        include_directories : [ inc_fftlwrap,
                                                inc_intrp ],
		        dependencies : [ gsl_dep,
                                         fftlog_dep,
                                         utl_dep,
                                         intrp_dep ],
                        install : true,
                        install_dir : libdir )
lib_headers += [ root + '/fftlog_wrap/' + inc_dir + '/fftlog_handler.h',
                 root + '/fftlog_wrap/' + inc_dir + '/abel_transform.h' ]
else
fftlwrap_lib = static_library( 'fftlog_wrap',
                               [ 'src/fftlog_handler.cpp' ],
       	 	               include_directories : [ inc_fftlwrap,
                                                       inc_intrp ],
		               dependencies : [ gsl_dep,
                                                fftlog_dep,
                                                utl_dep,
                                                intrp_dep ],
                               install : true,
                               install_dir : libdir )
endif


fftlwrap_dep = declare_dependency( link_with : fftlwrap_lib,
                                   include_directories : [ inc_fftlwrap,
                                                           inc_intrp ] )

#######################################################################
