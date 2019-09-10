#######################################################################
# Internal includes

inc_dir = 'include'

inc_intrp = include_directories( inc_dir )

#######################################################################
# Interpolation Lib

intrp_lib = library( 'interpolation',
                     [ 'src/interpolation.cpp' ],
       	 	     include_directories : inc_intrp,
		     dependencies : [ gsl_dep,
                                      utl_dep ],
                     install : true )

intrp_dep = declare_dependency( link_with : intrp_lib,
                                include_directories : inc_intrp )

lib_headers += [ root + '/interpolator/' + inc_dir + '/interpolation_interface.h',
                 root + '/interpolator/' + inc_dir + '/gsl_interpolation_interface.h',
                 root + '/interpolator/' + inc_dir + '/interpolation.h' ]

#######################################################################
