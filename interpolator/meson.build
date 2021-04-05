#######################################################################
# Internal includes

inc_dir = 'include'

inc_intrp = include_directories( inc_dir )

#######################################################################
# Interpolation Lib

if get_option( 'full-build' )
intrp_lib = library( 'interpolation',
                     [ 'src/interpolation.cpp' ],
       	 	     include_directories : inc_intrp,
		     dependencies : [ gsl_dep,
                                      utl_dep ],
                     install : true,
                     install_dir : libdir )
else
intrp_lib = static_library( 'interpolation',
                            [ 'src/interpolation.cpp' ],
       	 	            include_directories : inc_intrp,
		            dependencies : [ gsl_dep,
                                             utl_dep ],
                            install : true,
                            install_dir : libdir )
endif

intrp_dep = declare_dependency( link_with : intrp_lib,
                                include_directories : inc_intrp )

#######################################################################
# Interpolation Lib C interface

intrp_c_lib = library( 'interpolation_c_wrap',
                       [ 'src/interpolation_c_interface.cpp' ],
                       include_directories : inc_intrp,
                       dependencies : [ intrp_dep,
                                        gsl_dep,
                                        utl_dep ],
                       install : true,
                       install_dir : libdir )

intrp_c_dep = declare_dependency( link_with : intrp_c_lib,
                                  include_directories : [ inc_intrp ] )

intrp_c_lib_name = intrp_c_lib.full_path().split('/')[-1]
intrp_c_path = '/'.join( [ get_option( 'prefix' ),
                           get_option( 'libdir' ),
                           intrp_c_lib_name ] )

#######################################################################
# Add headers to list

if get_option( 'full-build' )
lib_headers += [ root + '/interpolator/' + inc_dir + '/interpolation_interface.h',
                 root + '/interpolator/' + inc_dir + '/gsl_interpolation_interface.h',
                 root + '/interpolator/' + inc_dir + '/interpolation.h',
                 root + '/interpolator/' + inc_dir + '/interpolation_c_interface.h' ]
endif

#######################################################################
